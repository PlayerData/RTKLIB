/*------------------------------------------------------------------------------
* accel_prn.c : per-axis IMU-accel sidecar — integrator over GPS epochs.
*
* Loads (time, a_E, a_N, a_U) rows; integrates Σ a²(τ) dτ per axis over a
* GPS epoch interval using a piecewise-constant model. udpos combines that
* with the configured baseline (rtk->opt.prn_imu_acch/_accv) — see
* rtkpos.c for the assembled Q.
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "rtklib.h"
#include "accel_prn.h"

typedef struct {
    gtime_t t;
    double  ae;   /* m/s² nav-frame east */
    double  an;   /* m/s² nav-frame north */
    double  au;   /* m/s² nav-frame up */
} accel_sample_t;

static accel_sample_t *g_samples = NULL;
static int             g_n        = 0;

static void die(const char *fmt, ...) {
    va_list ap;
    fprintf(stderr, "accel_prn: FATAL: ");
    va_start(ap, fmt); vfprintf(stderr, fmt, ap); va_end(ap);
    fputc('\n', stderr);
    exit(EXIT_FAILURE);
}

/* parse "YYYY-MM-DDTHH:MM:SS[.fff][Z]" into gtime_t (treated as GPST) */
static int parse_iso8601_gpst(const char *s, gtime_t *t) {
    double ep[6] = {0};
    int n = sscanf(s, "%lf-%lf-%lfT%lf:%lf:%lf",
                   &ep[0], &ep[1], &ep[2], &ep[3], &ep[4], &ep[5]);
    if (n != 6) return -1;
    *t = epoch2time(ep);
    return 0;
}

/* largest i such that g_samples[i].t <= t, or -1 if t precedes first sample */
static int find_idx(gtime_t t) {
    if (g_n == 0) return -1;
    if (timediff(t, g_samples[0].t) < 0.0) return -1;
    int lo = 0, hi = g_n - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (timediff(g_samples[mid].t, t) <= 0.0) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

int accel_prn_load(const char *path) {
    if (!path || !*path) return 0;

    FILE *fp = fopen(path, "r");
    if (!fp) die("cannot open %s: %s", path, strerror(errno));

    char line[512];
    int  lineno = 0, cap = 1024;
    g_samples = (accel_sample_t*)malloc((size_t)cap * sizeof(*g_samples));
    if (!g_samples) die("oom");
    g_n = 0;

    if (!fgets(line, sizeof(line), fp)) die("%s: empty file", path);
    lineno++;

    while (fgets(line, sizeof(line), fp)) {
        lineno++;
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '#') continue;

        char tbuf[64];
        double ae = 0.0, an = 0.0, au = 0.0;
        /* Positional parse: first 4 columns are (time, a_E, a_N, a_U) — column
           names not enforced. Trailing columns (if any) are ignored. */
        if (sscanf(line, "%63[^,],%lf,%lf,%lf", tbuf, &ae, &an, &au) != 4)
            die("%s:%d malformed row (need time,a_E,a_N,a_U)", path, lineno);

        gtime_t t;
        if (parse_iso8601_gpst(tbuf, &t) < 0)
            die("%s:%d unparseable time '%s'", path, lineno, tbuf);

        if (g_n > 0 && timediff(t, g_samples[g_n-1].t) <= 0.0)
            die("%s:%d timestamps not strictly increasing", path, lineno);

        if (g_n == cap) {
            cap *= 2;
            g_samples = (accel_sample_t*)realloc(
                g_samples, (size_t)cap * sizeof(*g_samples));
            if (!g_samples) die("oom");
        }
        g_samples[g_n].t  = t;
        g_samples[g_n].ae = ae;
        g_samples[g_n].an = an;
        g_samples[g_n].au = au;
        g_n++;
    }
    fclose(fp);

    if (g_n == 0) die("%s: no data rows", path);

    fprintf(stderr, "accel_prn: loaded %d IMU samples from %s\n", g_n, path);
    return g_n;
}

int accel_prn_loaded(void) { return g_n > 0; }

int accel_prn_integrate(gtime_t t_lo_in, gtime_t t_hi_in,
                        double *qe2, double *qn2, double *qu2) {
    if (g_n == 0) return 0;
    if (t_lo_in.time == 0 || t_hi_in.time == 0) return 0;

    gtime_t t_lo = t_lo_in, t_hi = t_hi_in;
    if (timediff(t_hi, t_lo) < 0.0) {
        gtime_t tmp = t_lo; t_lo = t_hi; t_hi = tmp;
    }

    /* Strict coverage. See header for rationale. */
    if (timediff(t_lo, g_samples[0].t)     < 0.0) return 0;
    if (timediff(t_hi, g_samples[g_n-1].t) > 0.0) return 0;

    /* Zero-width interval: integral is zero exactly. Treat as covered so
     * the caller still adds its own config²·dt term (which is also zero
     * when dt is 0). */
    if (timediff(t_hi, t_lo) == 0.0) {
        *qe2 = *qn2 = *qu2 = 0.0;
        return 1;
    }

    int i_start = find_idx(t_lo);
    int i_end   = find_idx(t_hi);
    if (i_start < 0) return 0;

    double sum_e = 0.0, sum_n = 0.0, sum_u = 0.0;
    for (int i = i_start; i <= i_end; i++) {
        gtime_t left  = (timediff(g_samples[i].t,   t_lo) > 0.0)
                          ? g_samples[i].t   : t_lo;
        gtime_t right = (timediff(g_samples[i+1].t, t_hi) < 0.0)
                          ? g_samples[i+1].t : t_hi;
        double dt = timediff(right, left);
        if (dt <= 0.0) continue;
        sum_e += g_samples[i].ae * g_samples[i].ae * dt;
        sum_n += g_samples[i].an * g_samples[i].an * dt;
        sum_u += g_samples[i].au * g_samples[i].au * dt;
    }
    *qe2 = sum_e;
    *qn2 = sum_n;
    *qu2 = sum_u;
    return 1;
}

void accel_prn_free(void) {
    free(g_samples);
    g_samples = NULL;
    g_n = 0;
}
