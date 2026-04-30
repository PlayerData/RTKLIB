/*------------------------------------------------------------------------------
* t_accel_prn.c : low-hanging unit tests for accel_prn.{c,h}
*
* Build (standalone, no cmake / lapack required):
*   cc -std=c99 -Wall -O2 -I../../src -DTRACE \
*      t_accel_prn.c \
*      ../../src/accel_prn.c ../../src/rtkcmn.c ../../src/trace.c \
*      -lm -o t_accel_prn
*
* `accel_prn_mean` returns the time-weighted mean of (a_E, a_N, a_U) over
* [t_lo, t_hi] (m/s²). udimu() in rtkpos.c uses this as a 3-vector
* measurement of the EKF accel state in ENU.
*
* Tests covered:
*   M1   single-sample window    → returns the sample's value
*   M2   multi-sample window     → time-weighted mean across samples
*   M3   single sample, exact-time partial window → that sample's value
*   M4   straddles first sample  → no coverage
*   M5   straddles last sample   → no coverage
*   M6   reversed argument order → same result
*   M7   property test           → matches naive time-weighted mean
*   M8   trailing-column toleration
*   B1   frame identity at lat=0, lon=0
*   B3   trace + symmetry invariance under ENU->ECEF rotation
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include "rtklib.h"
#include "accel_prn.h"

static int g_failed = 0;

#define CHECK(cond, msg) do {                                                  \
    if (!(cond)) {                                                             \
        fprintf(stderr, "  FAIL [%s:%d] %s\n", __func__, __LINE__, msg);       \
        g_failed = 1;                                                          \
    }                                                                          \
} while (0)

#define CHECK_NEAR(a, b, tol, msg) do {                                        \
    double _da = (a), _db = (b);                                               \
    if (fabs(_da - _db) > (tol)) {                                             \
        fprintf(stderr, "  FAIL [%s:%d] %s: got %.12g expected %.12g (|d|=%.3e > %.3e)\n", \
                __func__, __LINE__, msg, _da, _db, fabs(_da-_db), (tol));      \
        g_failed = 1;                                                          \
    }                                                                          \
} while (0)

static gtime_t base_time(void) {
    double ep0[6] = {2024,1,15,12,0,0};
    return epoch2time(ep0);
}

static char *write_csv(const double *offsets, int n,
                       const double *e, const double *nn, const double *u) {
    static char path[64];
    snprintf(path, sizeof(path), "/tmp/t_accel_prn_%d.csv", (int)getpid());
    FILE *fp = fopen(path, "w");
    assert(fp);
    fputs("time,ae_mps2,an_mps2,au_mps2\n", fp);
    gtime_t t0 = base_time();
    for (int i = 0; i < n; i++) {
        gtime_t t = timeadd(t0, offsets[i]);
        char ts[64]; time2str(t, ts, 9);
        ts[4]='-'; ts[7]='-'; ts[10]='T';
        fprintf(fp, "%sZ,%.10f,%.10f,%.10f\n", ts, e[i], nn[i], u[i]);
    }
    fclose(fp);
    return path;
}

/* ---------- M1: single-sample window ----------------------------------- */
static void test_single_sample_window(void) {
    fprintf(stderr, "test_single_sample_window\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {-3.0, 4.0, 1.0};
    double u[] = {0.5, -1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double ze, zn, zu;
    /* [0.3, 0.7] entirely inside sample 0's [0,1) extent → mean = sample 0 */
    int ok = accel_prn_mean(timeadd(t0, 0.3), timeadd(t0, 0.7), &ze, &zn, &zu);
    CHECK(ok == 1, "M1: should be covered");
    CHECK_NEAR(ze, e[0],  1e-9, "M1: mean ze = e0");
    CHECK_NEAR(zn, nn[0], 1e-9, "M1: signed values pass through");
    CHECK_NEAR(zu, u[0],  1e-9, "M1: mean zu = u0");
    accel_prn_free(); unlink(p);
}

/* ---------- M2: multi-sample window ------------------------------------ */
static void test_multi_sample_window(void) {
    fprintf(stderr, "test_multi_sample_window\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double ze, zn, zu;
    /* [0.5, 1.8] → sample 0 covers [0.5,1)=0.5s, sample 1 covers [1,1.8)=0.8s */
    int ok = accel_prn_mean(timeadd(t0, 0.5), timeadd(t0, 1.8), &ze, &zn, &zu);
    CHECK(ok == 1, "M2: should be covered");
    double dt = 1.3;
    double exp_e = (e[0]*0.5 + e[1]*0.8) / dt;
    double exp_n = (nn[0]*0.5 + nn[1]*0.8) / dt;
    double exp_u = (u[0]*0.5 + u[1]*0.8) / dt;
    CHECK_NEAR(ze, exp_e, 1e-9, "M2: time-weighted mean E");
    CHECK_NEAR(zn, exp_n, 1e-9, "M2: time-weighted mean N");
    CHECK_NEAR(zu, exp_u, 1e-9, "M2: time-weighted mean U");
    accel_prn_free(); unlink(p);
}

/* ---------- M3: zero-width interval is rejected (no measurement) ------- */
static void test_zero_width(void) {
    fprintf(stderr, "test_zero_width\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double ze = 999, zn = 999, zu = 999;
    gtime_t t = timeadd(t0, 0.5);
    int ok = accel_prn_mean(t, t, &ze, &zn, &zu);
    CHECK(ok == 0, "M3: zero-width interval is not a valid measurement");
    CHECK(ze == 999 && zn == 999 && zu == 999,
          "M3: outputs untouched on zero-width");
    accel_prn_free(); unlink(p);
}

/* ---------- M4 / M5: straddle edges → no coverage ---------------------- */
static void test_straddle_edges(void) {
    fprintf(stderr, "test_straddle_edges\n");
    double offsets[] = {1.0, 2.0, 3.0};
    double e[] = {2,2,2}, nn[] = {1,1,1}, u[] = {1,1,1};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double ze = 999, zn = 999, zu = 999;
    int ok = accel_prn_mean(timeadd(t0, -1.0), timeadd(t0, 2.0), &ze, &zn, &zu);
    CHECK(ok == 0, "M4: straddles first sample → no coverage");
    CHECK(ze == 999, "M4: outputs untouched");
    ok = accel_prn_mean(timeadd(t0, 2.5), timeadd(t0, 4.0), &ze, &zn, &zu);
    CHECK(ok == 0, "M5: straddles last sample → no coverage");
    accel_prn_free(); unlink(p);
}

/* ---------- M6: reversed args ------------------------------------------ */
static void test_reversed_args(void) {
    fprintf(stderr, "test_reversed_args\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3, 4, 1}, u[] = {0.5, 1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double ze_f, zn_f, zu_f, ze_r, zn_r, zu_r;
    accel_prn_mean(timeadd(t0, 0.5), timeadd(t0, 1.8), &ze_f, &zn_f, &zu_f);
    accel_prn_mean(timeadd(t0, 1.8), timeadd(t0, 0.5), &ze_r, &zn_r, &zu_r);
    CHECK_NEAR(ze_f, ze_r, 1e-12, "M6: arg-order invariance E");
    CHECK_NEAR(zn_f, zn_r, 1e-12, "M6: arg-order invariance N");
    CHECK_NEAR(zu_f, zu_r, 1e-12, "M6: arg-order invariance U");
    accel_prn_free(); unlink(p);
}

/* ---------- M7: property test against naive time-weighted mean --------- */
static double naive_mean_axis(const double *offsets, const double *vals,
                              int n, double t_lo, double t_hi) {
    double num = 0.0, dt_total = t_hi - t_lo;
    for (int i = 0; i < n - 1; i++) {
        double left  = offsets[i]   > t_lo ? offsets[i]   : t_lo;
        double right = offsets[i+1] < t_hi ? offsets[i+1] : t_hi;
        if (right > left) num += vals[i] * (right - left);
    }
    return num / dt_total;
}

static void test_mean_property(void) {
    fprintf(stderr, "test_mean_property\n");
    const int N = 500;
    const int Q = 200;
    double *offsets = malloc(sizeof(double)*N);
    double *e = malloc(sizeof(double)*N);
    double *nn = malloc(sizeof(double)*N);
    double *u = malloc(sizeof(double)*N);
    srand(123);
    double cur = 0.0;
    for (int i = 0; i < N; i++) {
        cur += 0.05 + (rand()/(double)RAND_MAX) * 0.5;
        offsets[i] = cur;
        e[i]  = -10.0 + (rand()/(double)RAND_MAX) * 20.0;
        nn[i] = -10.0 + (rand()/(double)RAND_MAX) * 20.0;
        u[i]  = -10.0 + (rand()/(double)RAND_MAX) * 20.0;
    }
    char *p = write_csv(offsets, N, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();

    int mismatches = 0;
    double t_first = offsets[0], t_last = offsets[N-1];
    for (int q = 0; q < Q; q++) {
        double a = t_first + (rand()/(double)RAND_MAX) * (t_last - t_first - 0.001);
        double b = a + 0.001 + (rand()/(double)RAND_MAX) * (t_last - a - 0.001);
        if (b > t_last) b = t_last;
        double ze, zn, zu;
        int ok = accel_prn_mean(timeadd(t0, a), timeadd(t0, b), &ze, &zn, &zu);
        if (!ok) { mismatches++; continue; }
        double exp_e = naive_mean_axis(offsets, e,  N, a, b);
        double exp_n = naive_mean_axis(offsets, nn, N, a, b);
        double exp_u = naive_mean_axis(offsets, u,  N, a, b);
        /* Mean involves Σa·Δt / dt; the dt cancellation is not exact
           after the gtime_t round-trip, so allow a generous absolute
           floor on top of the relative tolerance. */
        double rtol = 1e-7, atol = 1e-7;
        if (fabs(ze - exp_e) > rtol*fabs(exp_e) + atol ||
            fabs(zn - exp_n) > rtol*fabs(exp_n) + atol ||
            fabs(zu - exp_u) > rtol*fabs(exp_u) + atol) {
            mismatches++;
            fprintf(stderr, "  q=%d a=%.6f b=%.6f ze=%.6f exp=%.6f\n",
                    q, a, b, ze, exp_e);
        }
    }
    CHECK(mismatches == 0, "M7: time-weighted mean must match naive");
    free(offsets); free(e); free(nn); free(u);
    accel_prn_free(); unlink(p);
}

/* ---------- M8: trailing column toleration ---------------------------- */
static void test_trailing_columns(void) {
    fprintf(stderr, "test_trailing_columns\n");
    char tmp[64];
    snprintf(tmp, sizeof(tmp), "/tmp/t_accel_prn_5col_%d.csv", (int)getpid());
    FILE *fp = fopen(tmp, "w");
    assert(fp);
    fputs("time,ae,an,au,extra\n", fp);
    gtime_t t0 = base_time();
    for (int i = 0; i < 3; i++) {
        gtime_t t = timeadd(t0, (double)i);
        char ts[64]; time2str(t, ts, 9);
        ts[4]='-'; ts[7]='-'; ts[10]='T';
        fprintf(fp, "%sZ,%.6f,%.6f,%.6f,99.9\n", ts, 1.0, 2.0, 3.0);
    }
    fclose(fp);
    int n = accel_prn_load(tmp);
    CHECK(n == 3, "M8: 5-column CSV loads (extra column ignored)");
    double ze, zn, zu;
    int ok = accel_prn_mean(timeadd(t0, 0.5), timeadd(t0, 1.5), &ze, &zn, &zu);
    CHECK(ok == 1 && fabs(ze - 1.0) < 1e-9 && fabs(zn - 2.0) < 1e-9
                  && fabs(zu - 3.0) < 1e-9,
          "M8: mean of constant samples = sample value");
    accel_prn_free(); unlink(tmp);
}

/* ---------- B1: frame identity at lat=0, lon=0 ------------------------- */
static void test_frame_at_origin(void) {
    fprintf(stderr, "test_frame_at_origin\n");
    double pos[3] = {0.0, 0.0, 0.0};
    double prn_e = 1.5, prn_n = 2.7, prn_u = 0.4;
    double Q[9] = {0};
    double Qv[9];
    Q[0] = prn_e * prn_e;
    Q[4] = prn_n * prn_n;
    Q[8] = prn_u * prn_u;
    covecef(pos, Q, Qv);
    CHECK_NEAR(Qv[0], prn_u*prn_u, 1e-12, "B1: Qv[X,X] = Qu");
    CHECK_NEAR(Qv[4], prn_e*prn_e, 1e-12, "B1: Qv[Y,Y] = Qe");
    CHECK_NEAR(Qv[8], prn_n*prn_n, 1e-12, "B1: Qv[Z,Z] = Qn");
    for (int i = 0; i < 9; i++) {
        if (i == 0 || i == 4 || i == 8) continue;
        CHECK_NEAR(Qv[i], 0.0, 1e-12, "B1: off-diagonals zero");
    }
}

/* ---------- B3: trace + symmetry invariance ---------------------------- */
static void test_trace_invariance(void) {
    fprintf(stderr, "test_trace_invariance\n");
    srand(7);
    for (int trial = 0; trial < 50; trial++) {
        double lat = (-89.0 + 178.0*(rand()/(double)RAND_MAX)) * (3.14159265358979/180.0);
        double lon = (-179.0 + 358.0*(rand()/(double)RAND_MAX)) * (3.14159265358979/180.0);
        double pos[3] = {lat, lon, 0.0};
        double pe = rand()/(double)RAND_MAX;
        double pn = rand()/(double)RAND_MAX;
        double pu = rand()/(double)RAND_MAX;
        double Q[9] = {0}, Qv[9];
        Q[0] = pe*pe; Q[4] = pn*pn; Q[8] = pu*pu;
        covecef(pos, Q, Qv);
        double tr_q = Q[0] + Q[4] + Q[8];
        double tr_qv = Qv[0] + Qv[4] + Qv[8];
        CHECK_NEAR(tr_qv, tr_q, 1e-10, "B3: trace preserved");
        CHECK_NEAR(Qv[1], Qv[3], 1e-10, "B3: symmetric (1,3)");
        CHECK_NEAR(Qv[2], Qv[6], 1e-10, "B3: symmetric (2,6)");
        CHECK_NEAR(Qv[5], Qv[7], 1e-10, "B3: symmetric (5,7)");
    }
}

int main(void) {
    traceopen("/dev/null");
    tracelevel(0);

    test_single_sample_window();
    test_multi_sample_window();
    test_zero_width();
    test_straddle_edges();
    test_reversed_args();
    test_mean_property();
    test_trailing_columns();
    test_frame_at_origin();
    test_trace_invariance();

    traceclose();
    if (g_failed) { fprintf(stderr, "FAILED\n"); return EXIT_FAILURE; }
    fprintf(stderr, "OK\n");
    return EXIT_SUCCESS;
}
