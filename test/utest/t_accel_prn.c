/*------------------------------------------------------------------------------
* t_accel_prn.c : low-hanging unit tests for accel_prn.{c,h}
*
* Build (standalone, no cmake / lapack required):
*   cc -std=c99 -Wall -O2 -I../../src -DTRACE \
*      t_accel_prn.c \
*      ../../src/accel_prn.c ../../src/rtkcmn.c ../../src/trace.c \
*      -lm -o t_accel_prn
*
* Run:
*   ./t_accel_prn   # exits 0 on pass, nonzero on first failure
*
* `accel_prn_integrate` returns Σᵢ aᵢ² · Δtᵢ per axis (m²/s³). It does
* NOT include the configured baseline prn_imu² · dt; that part lives in
* udpos. Tests therefore verify only the integral over CSV samples.
*
* Tests covered:
*   I1   single-sample window    → e²·dt etc.
*   I2   multi-sample window     → sum of rectangles
*   I3   zero-width interval     → integral = 0, returns covered
*   I4   straddles first sample  → no coverage
*   I5   straddles last sample   → no coverage
*   I6   reversed argument order → same result
*   I7   property test           → matches naive Riemann sum
*   I8   trailing-column toleration
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

/* ---------- I1: single-sample window ----------------------------------- */
static void test_single_sample_window(void) {
    fprintf(stderr, "test_single_sample_window\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {-3.0, 4.0, 1.0};
    double u[] = {0.5, -1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe2, qn2, qu2;
    int ok = accel_prn_integrate(timeadd(t0, 0.3), timeadd(t0, 0.7),
                                 &qe2, &qn2, &qu2);
    CHECK(ok == 1, "I1: should be covered");
    /* Sample 0 covers [0, 1); we integrate [0.3, 0.7] inside it. */
    CHECK_NEAR(qe2, e[0]*e[0]*0.4,  1e-9, "I1: qe² = e0²·0.4");
    CHECK_NEAR(qn2, nn[0]*nn[0]*0.4, 1e-9, "I1: signed values squared OK");
    CHECK_NEAR(qu2, u[0]*u[0]*0.4,  1e-9, "I1: qu² = u0²·0.4");
    accel_prn_free(); unlink(p);
}

/* ---------- I2: window spans multiple samples -------------------------- */
static void test_multi_sample_window(void) {
    fprintf(stderr, "test_multi_sample_window\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe2, qn2, qu2;
    int ok = accel_prn_integrate(timeadd(t0, 0.5), timeadd(t0, 1.8),
                                 &qe2, &qn2, &qu2);
    CHECK(ok == 1, "I2: should be covered");
    /* sample 0: [0,1)→[0.5,1)=0.5; sample 1: [1,2)→[1,1.8)=0.8 */
    double exp_e = e[0]*e[0]*0.5 + e[1]*e[1]*0.8;
    double exp_n = nn[0]*nn[0]*0.5 + nn[1]*nn[1]*0.8;
    double exp_u = u[0]*u[0]*0.5 + u[1]*u[1]*0.8;
    CHECK_NEAR(qe2, exp_e, 1e-9, "I2: qe² across samples");
    CHECK_NEAR(qn2, exp_n, 1e-9, "I2: qn² across samples");
    CHECK_NEAR(qu2, exp_u, 1e-9, "I2: qu² across samples");
    accel_prn_free(); unlink(p);
}

/* ---------- I3: zero-width interval ------------------------------------ */
static void test_zero_width(void) {
    fprintf(stderr, "test_zero_width\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe2 = 999, qn2 = 999, qu2 = 999;
    gtime_t t = timeadd(t0, 0.5);
    int ok = accel_prn_integrate(t, t, &qe2, &qn2, &qu2);
    CHECK(ok == 1, "I3: zero-width covered");
    CHECK_NEAR(qe2, 0.0, 0.0, "I3: qe² = 0 exactly");
    CHECK_NEAR(qn2, 0.0, 0.0, "I3: qn² = 0 exactly");
    CHECK_NEAR(qu2, 0.0, 0.0, "I3: qu² = 0 exactly");
    accel_prn_free(); unlink(p);
}

/* ---------- I4 / I5: straddle edges → no coverage ---------------------- */
static void test_straddle_edges(void) {
    fprintf(stderr, "test_straddle_edges\n");
    double offsets[] = {1.0, 2.0, 3.0};
    double e[] = {2,2,2}, nn[] = {1,1,1}, u[] = {1,1,1};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe2 = 999, qn2 = 999, qu2 = 999;
    int ok = accel_prn_integrate(timeadd(t0, -1.0), timeadd(t0, 2.0),
                                 &qe2, &qn2, &qu2);
    CHECK(ok == 0, "I4: straddles first sample → not covered");
    CHECK(qe2 == 999 && qn2 == 999 && qu2 == 999,
          "I4: outputs untouched on no-coverage");
    ok = accel_prn_integrate(timeadd(t0, 2.5), timeadd(t0, 4.0),
                             &qe2, &qn2, &qu2);
    CHECK(ok == 0, "I5: straddles last sample → not covered");
    accel_prn_free(); unlink(p);
}

/* ---------- I6: reversed args ------------------------------------------ */
static void test_reversed_args(void) {
    fprintf(stderr, "test_reversed_args\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3, 4, 1}, u[] = {0.5, 1.5, 2.5};
    char *p = write_csv(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe_f, qn_f, qu_f, qe_r, qn_r, qu_r;
    accel_prn_integrate(timeadd(t0, 0.5), timeadd(t0, 1.8),
                        &qe_f, &qn_f, &qu_f);
    accel_prn_integrate(timeadd(t0, 1.8), timeadd(t0, 0.5),
                        &qe_r, &qn_r, &qu_r);
    CHECK_NEAR(qe_f, qe_r, 1e-12, "I6: arg-order invariance");
    CHECK_NEAR(qn_f, qn_r, 1e-12, "I6: arg-order invariance");
    CHECK_NEAR(qu_f, qu_r, 1e-12, "I6: arg-order invariance");
    accel_prn_free(); unlink(p);
}

/* ---------- I7: property test against naive Riemann sum ---------------- */
static double naive_integrate_axis(const double *offsets, const double *vals,
                                   int n, double t_lo, double t_hi) {
    double sum = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double left  = offsets[i]   > t_lo ? offsets[i]   : t_lo;
        double right = offsets[i+1] < t_hi ? offsets[i+1] : t_hi;
        if (right > left) sum += vals[i] * vals[i] * (right - left);
    }
    return sum;
}

static void test_integrate_property(void) {
    fprintf(stderr, "test_integrate_property\n");
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
        e[i]  = -10.0 + (rand()/(double)RAND_MAX) * 20.0;  /* signed */
        nn[i] = -10.0 + (rand()/(double)RAND_MAX) * 20.0;
        u[i]  = -10.0 + (rand()/(double)RAND_MAX) * 20.0;
    }
    char *p = write_csv(offsets, N, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();

    int mismatches = 0;
    double t_first = offsets[0];
    double t_last  = offsets[N-1];
    for (int q = 0; q < Q; q++) {
        double a = t_first + (rand()/(double)RAND_MAX) * (t_last - t_first - 0.001);
        double b = a + 0.001 + (rand()/(double)RAND_MAX) * (t_last - a - 0.001);
        if (b > t_last) b = t_last;
        double qe2, qn2, qu2;
        int ok = accel_prn_integrate(timeadd(t0, a), timeadd(t0, b),
                                     &qe2, &qn2, &qu2);
        if (!ok) { mismatches++; continue; }
        double ee = naive_integrate_axis(offsets, e,  N, a, b);
        double nx = naive_integrate_axis(offsets, nn, N, a, b);
        double uu = naive_integrate_axis(offsets, u,  N, a, b);
        double rtol = 1e-7;
        if (fabs(qe2 - ee) > rtol*fabs(ee) + 1e-9 ||
            fabs(qn2 - nx) > rtol*fabs(nx) + 1e-9 ||
            fabs(qu2 - uu) > rtol*fabs(uu) + 1e-9) {
            mismatches++;
            fprintf(stderr, "  q=%d a=%.6f b=%.6f qe²=%.6f exp=%.6f rel=%.2e\n",
                    q, a, b, qe2, ee, fabs(qe2-ee)/fmax(fabs(ee),1e-12));
        }
    }
    CHECK(mismatches == 0, "I7: integrator must match naive Riemann sum");
    free(offsets); free(e); free(nn); free(u);
    accel_prn_free(); unlink(p);
}

/* ---------- I8: trailing column toleration ---------------------------- */
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
    CHECK(n == 3, "I8: 5-column CSV loads (extra column ignored)");
    double qe2, qn2, qu2;
    int ok = accel_prn_integrate(timeadd(t0, 0.5), timeadd(t0, 1.5),
                                 &qe2, &qn2, &qu2);
    CHECK(ok == 1 && fabs(qe2 - 1.0) < 1e-9 && fabs(qn2 - 4.0) < 1e-9
                  && fabs(qu2 - 9.0) < 1e-9,
          "I8: integration uses parsed values; trailing col ignored");
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
    test_integrate_property();
    test_trailing_columns();
    test_frame_at_origin();
    test_trace_invariance();

    traceclose();
    if (g_failed) { fprintf(stderr, "FAILED\n"); return EXIT_FAILURE; }
    fprintf(stderr, "OK\n");
    return EXIT_SUCCESS;
}
