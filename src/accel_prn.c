/*------------------------------------------------------------------------------
* accel_prn.c : IMU-coverage signal for EKF accel-state process noise.
*
* Tracks just the timestamps of loaded IMU samples; the prn-std values
* applied during/outside coverage live in prcopt_t.prn_imu_acch/_accv
* (RTKLIB conf options stats-prnaccelh-imu / stats-prnaccelv-imu) and
* prcopt_t.prn[3..4] (stats-prnaccelh / stats-prnaccelv) respectively.
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "rtklib.h"
#include "accel_prn.h"

static gtime_t *g_times = NULL;
static int      g_n     = 0;

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

int accel_prn_load(const char *path) {
    if (!path || !*path) return 0;

    FILE *fp = fopen(path, "r");
    if (!fp) die("cannot open %s: %s", path, strerror(errno));

    char line[512];
    int  lineno = 0, cap = 1024;
    g_times = (gtime_t*)malloc((size_t)cap * sizeof(*g_times));
    if (!g_times) die("oom");
    g_n = 0;

    if (!fgets(line, sizeof(line), fp)) die("%s: empty file", path);
    lineno++;

    while (fgets(line, sizeof(line), fp)) {
        lineno++;
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '#') continue;

        char tbuf[64];
        /* Parse only the first column. Trailing columns (if any) are ignored —
           this file's role is coverage-only; the prn values applied during
           coverage come from the RTKLIB conf, not the CSV. */
        if (sscanf(line, "%63[^,\n]", tbuf) != 1)
            die("%s:%d malformed row (need at least a time column)", path, lineno);

        gtime_t t;
        if (parse_iso8601_gpst(tbuf, &t) < 0)
            die("%s:%d unparseable time '%s'", path, lineno, tbuf);

        if (g_n > 0 && timediff(t, g_times[g_n-1]) <= 0.0)
            die("%s:%d timestamps not strictly increasing", path, lineno);

        if (g_n == cap) {
            cap *= 2;
            g_times = (gtime_t*)realloc(g_times, (size_t)cap * sizeof(*g_times));
            if (!g_times) die("oom");
        }
        g_times[g_n++] = t;
    }
    fclose(fp);

    if (g_n == 0) die("%s: no data rows", path);

    fprintf(stderr, "accel_prn: loaded %d coverage timestamps from %s\n",
            g_n, path);
    return g_n;
}

int accel_prn_loaded(void) { return g_n > 0; }

int accel_prn_covered(gtime_t t_lo_in, gtime_t t_hi_in) {
    if (g_n == 0) return 0;
    if (t_lo_in.time == 0 || t_hi_in.time == 0) return 0;

    gtime_t t_lo = t_lo_in, t_hi = t_hi_in;
    if (timediff(t_hi, t_lo) < 0.0) {
        gtime_t tmp = t_lo; t_lo = t_hi; t_hi = tmp;
    }

    if (timediff(t_lo, g_times[0])     < 0.0) return 0;
    if (timediff(t_hi, g_times[g_n-1]) > 0.0) return 0;
    return 1;
}

void accel_prn_free(void) {
    free(g_times);
    g_times = NULL;
    g_n = 0;
}
