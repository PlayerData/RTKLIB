/*------------------------------------------------------------------------------
* accel_prn.c : per-epoch override of EKF acceleration process-noise stds
*
* Reads a sidecar CSV with columns: time, prn_e, prn_n, prn_u.
* Hooked from udpos() (rtkpos.c) to set the per-axis local-ENU process
* noise diagonal Q for the current epoch. When no CSV is loaded (or the
* epoch precedes the first sample), udpos falls back to the existing
* config-defined isotropic-horizontal opt.prn[3]/opt.prn[4] behaviour.
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
    double  prn_e;
    double  prn_n;
    double  prn_u;
} accel_prn_sample_t;

static accel_prn_sample_t *g_samples = NULL;
static int g_n = 0;

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
    g_samples = (accel_prn_sample_t*)malloc((size_t)cap * sizeof(*g_samples));
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
        double pe = 0.0, pn = 0.0, pu = 0.0;
        if (sscanf(line, "%63[^,],%lf,%lf,%lf", tbuf, &pe, &pn, &pu) != 4)
            die("%s:%d malformed row (need time,prn_e,prn_n,prn_u)", path, lineno);

        gtime_t t;
        if (parse_iso8601_gpst(tbuf, &t) < 0)
            die("%s:%d unparseable time '%s'", path, lineno, tbuf);

        if (pe < 0.0 || pn < 0.0 || pu < 0.0)
            die("%s:%d negative process-noise std", path, lineno);

        if (g_n > 0 && timediff(t, g_samples[g_n-1].t) <= 0.0)
            die("%s:%d timestamps not strictly increasing", path, lineno);

        if (g_n == cap) {
            cap *= 2;
            g_samples = (accel_prn_sample_t*)realloc(
                g_samples, (size_t)cap * sizeof(*g_samples));
            if (!g_samples) die("oom");
        }
        g_samples[g_n].t     = t;
        g_samples[g_n].prn_e = pe;
        g_samples[g_n].prn_n = pn;
        g_samples[g_n].prn_u = pu;
        g_n++;
    }
    fclose(fp);

    if (g_n == 0) die("%s: no data rows", path);

    fprintf(stderr, "accel_prn: loaded %d samples from %s\n", g_n, path);
    return g_n;
}

int accel_prn_loaded(void) { return g_n > 0; }

int accel_prn_integrate(gtime_t t_lo_in, gtime_t t_hi_in,
                        double *qe, double *qn, double *qu) {
    if (g_n == 0) return 0;
    if (t_lo_in.time == 0 || t_hi_in.time == 0) return 0;

    /* Normalise so t_lo <= t_hi regardless of processing direction. */
    gtime_t t_lo = t_lo_in, t_hi = t_hi_in;
    if (timediff(t_hi, t_lo) < 0.0) {
        gtime_t tmp = t_lo; t_lo = t_hi; t_hi = tmp;
    }

    /* Strict coverage: require [t_lo, t_hi] entirely inside the loaded
     * sample range. We could partial-integrate and fall back for the
     * uncovered part of the interval, but that mixes two semantics and
     * makes the resulting Q hard to reason about. Cleaner to fall back
     * entirely when coverage is incomplete — caller picks one path. */
    if (timediff(t_lo, g_samples[0].t)        < 0.0) return 0;
    if (timediff(t_hi, g_samples[g_n-1].t)    > 0.0) return 0;

    /* Zero-width interval: integrate to zero exactly. Treat as covered
     * (return 1 with zeros) — at the very first epoch rtk->tt is 0 and
     * the caller should add a zero process-noise contribution, matching
     * what the stock formula prn^2 * fabs(0) = 0 also produces. */
    if (timediff(t_hi, t_lo) == 0.0) {
        *qe = *qn = *qu = 0.0;
        return 1;
    }

    /* Sample i covers [t_i, t_{i+1}). i_start is the sample whose extent
     * contains t_lo; i_end is the one whose extent contains t_hi. We
     * already required t_hi <= t_{n-1} so i_end <= n-2 and we can always
     * read t_{i_end+1}. */
    int i_start = find_idx(t_lo);
    int i_end   = find_idx(t_hi);
    if (i_start < 0) return 0;  /* defensive — covered by checks above */

    double sum_e = 0.0, sum_n = 0.0, sum_u = 0.0;
    for (int i = i_start; i <= i_end; i++) {
        gtime_t left  = (timediff(g_samples[i].t,   t_lo) > 0.0)
                          ? g_samples[i].t   : t_lo;
        gtime_t right = (timediff(g_samples[i+1].t, t_hi) < 0.0)
                          ? g_samples[i+1].t : t_hi;
        double dt = timediff(right, left);
        if (dt <= 0.0) continue;
        sum_e += g_samples[i].prn_e * g_samples[i].prn_e * dt;
        sum_n += g_samples[i].prn_n * g_samples[i].prn_n * dt;
        sum_u += g_samples[i].prn_u * g_samples[i].prn_u * dt;
    }
    *qe = sum_e;
    *qn = sum_n;
    *qu = sum_u;
    return 1;
}

void accel_prn_free(void) {
    free(g_samples);
    g_samples = NULL;
    g_n = 0;
}
