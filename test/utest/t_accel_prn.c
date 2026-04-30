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
* The integrator semantics: piecewise-constant left-aligned, sample i
* covers [t_i, t_{i+1}). The integral over [t_lo, t_hi] is
* Σ_i σ_i² · clip([t_i, t_{i+1}), [t_lo, t_hi]).
*
* Tests covered:
*   I1   covered interval, single sample window     (rectangle area)
*   I2   covered interval spanning multiple samples (sum of rectangles)
*   I3   zero-width interval                        (integral = 0)
*   I4   interval before first sample               (no coverage → 0)
*   I5   interval past last sample                  (no coverage → 0)
*   I6   reversed argument order                    (same result, normalised)
*   I7   property test: random CSVs and intervals vs naive integrator
*   B1   frame identity at lat=0, lon=0
*   B3   trace invariance under ENU->ECEF rotation
*   B5   active-branch parity with stock isotropic-h/v formula
*        (for constant prn over the interval, ∫prn²dt == prn²·(t_hi - t_lo))
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

static char *write_temp_csv(const char *body) {
    static char path[64];
    snprintf(path, sizeof(path), "/tmp/t_accel_prn_%d.csv", (int)getpid());
    FILE *fp = fopen(path, "w");
    assert(fp);
    fputs(body, fp);
    fclose(fp);
    return path;
}

/* Build a CSV from arrays of (offset_s, e, n, u) */
static char *write_csv_from(const double *offsets, int n,
                            const double *e, const double *nn, const double *u) {
    static char path[64];
    snprintf(path, sizeof(path), "/tmp/t_accel_prn_arr_%d.csv", (int)getpid());
    FILE *fp = fopen(path, "w");
    assert(fp);
    fputs("time,prn_e,prn_n,prn_u\n", fp);
    double ep0[6] = {2024,1,15,12,0,0};
    gtime_t t0 = epoch2time(ep0);
    for (int i = 0; i < n; i++) {
        gtime_t t = timeadd(t0, offsets[i]);
        char ts[64]; time2str(t, ts, 9);   /* nanosecond precision */
        ts[4]='-'; ts[7]='-'; ts[10]='T';
        fprintf(fp, "%sZ,%.10f,%.10f,%.10f\n", ts, e[i], nn[i], u[i]);
    }
    fclose(fp);
    return path;
}

static gtime_t base_time(void) {
    double ep0[6] = {2024,1,15,12,0,0};
    return epoch2time(ep0);
}

/* ---------- I1: single covered interval inside one sample's extent ----- */
static void test_single_sample_window(void) {
    fprintf(stderr, "test_single_sample_window\n");
    /* Three samples at 0,1,2 s. Interval [0.3, 0.7] entirely inside sample 0,
       which covers [0, 1). Expected integral = (e^2)*0.4 etc. */
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv_from(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe, qn, qu;
    int ok = accel_prn_integrate(timeadd(t0, 0.3), timeadd(t0, 0.7), &qe, &qn, &qu);
    CHECK(ok == 1, "I1: should be covered");
    CHECK_NEAR(qe, e[0]*e[0]*0.4,  1e-9, "I1: qe = e0^2 * 0.4");
    CHECK_NEAR(qn, nn[0]*nn[0]*0.4, 1e-9, "I1: qn = n0^2 * 0.4");
    CHECK_NEAR(qu, u[0]*u[0]*0.4,  1e-9, "I1: qu = u0^2 * 0.4");
    accel_prn_free(); unlink(p);
}

/* ---------- I2: interval spans multiple samples ------------------------ */
static void test_multi_sample_window(void) {
    fprintf(stderr, "test_multi_sample_window\n");
    /* Samples at 0,1,2 s. Interval [0.5, 1.8] crosses sample 0 and 1.
       Sample 0 covers [0,1), so contribution from it is e0^2 * (1.0 - 0.5) = e0^2 * 0.5
       Sample 1 covers [1,2), contribution is e1^2 * (1.8 - 1.0) = e1^2 * 0.8 */
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv_from(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe, qn, qu;
    int ok = accel_prn_integrate(timeadd(t0, 0.5), timeadd(t0, 1.8), &qe, &qn, &qu);
    CHECK(ok == 1, "I2: should be covered");
    double expect_e = e[0]*e[0]*0.5 + e[1]*e[1]*0.8;
    double expect_n = nn[0]*nn[0]*0.5 + nn[1]*nn[1]*0.8;
    double expect_u = u[0]*u[0]*0.5 + u[1]*u[1]*0.8;
    CHECK_NEAR(qe, expect_e, 1e-9, "I2: qe across two samples");
    CHECK_NEAR(qn, expect_n, 1e-9, "I2: qn across two samples");
    CHECK_NEAR(qu, expect_u, 1e-9, "I2: qu across two samples");
    accel_prn_free(); unlink(p);
}

/* ---------- I3: zero-width interval ------------------------------------ */
static void test_zero_width(void) {
    fprintf(stderr, "test_zero_width\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0};
    double nn[] = {3.0, 4.0, 1.0};
    double u[] = {0.5, 1.5, 2.5};
    char *p = write_csv_from(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe = 999, qn = 999, qu = 999;
    gtime_t t = timeadd(t0, 0.5);
    int ok = accel_prn_integrate(t, t, &qe, &qn, &qu);
    CHECK(ok == 1, "I3: zero-width interval is still covered");
    CHECK_NEAR(qe, 0.0, 0.0, "I3: zero-width gives qe=0 exactly");
    CHECK_NEAR(qn, 0.0, 0.0, "I3: zero-width gives qn=0 exactly");
    CHECK_NEAR(qu, 0.0, 0.0, "I3: zero-width gives qu=0 exactly");
    accel_prn_free(); unlink(p);
}

/* ---------- I4: interval before first sample → no coverage ------------- */
static void test_before_first(void) {
    fprintf(stderr, "test_before_first\n");
    double offsets[] = {1.0, 2.0, 3.0};
    double e[] = {2.0, 5.0, 7.0}, nn[] = {1, 1, 1}, u[] = {1, 1, 1};
    char *p = write_csv_from(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe = 999, qn = 999, qu = 999;
    int ok = accel_prn_integrate(timeadd(t0, -1.0), timeadd(t0, 0.5), &qe, &qn, &qu);
    CHECK(ok == 0, "I4: interval straddling first sample → no coverage");
    CHECK(qe == 999 && qn == 999 && qu == 999,
          "I4: outputs untouched on no-coverage");
    accel_prn_free(); unlink(p);
}

/* ---------- I5: interval past last sample → no coverage ---------------- */
static void test_past_last(void) {
    fprintf(stderr, "test_past_last\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0}, nn[] = {1, 1, 1}, u[] = {1, 1, 1};
    char *p = write_csv_from(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe = 999;
    /* Interval [1.5, 3.0] straddles t_{n-1} = 2.0 → no coverage. */
    int ok = accel_prn_integrate(timeadd(t0, 1.5), timeadd(t0, 3.0), &qe, &qe, &qe);
    CHECK(ok == 0, "I5: interval past last sample → no coverage");
    accel_prn_free(); unlink(p);
}

/* ---------- I6: reversed argument order ------------------------------- */
static void test_reversed_args(void) {
    fprintf(stderr, "test_reversed_args\n");
    double offsets[] = {0.0, 1.0, 2.0};
    double e[] = {2.0, 5.0, 7.0}, nn[] = {3, 4, 1}, u[] = {0.5, 1.5, 2.5};
    char *p = write_csv_from(offsets, 3, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    double qe_fwd, qn_fwd, qu_fwd, qe_rev, qn_rev, qu_rev;
    accel_prn_integrate(timeadd(t0, 0.5), timeadd(t0, 1.8), &qe_fwd, &qn_fwd, &qu_fwd);
    accel_prn_integrate(timeadd(t0, 1.8), timeadd(t0, 0.5), &qe_rev, &qn_rev, &qu_rev);
    CHECK_NEAR(qe_fwd, qe_rev, 1e-12, "I6: backward argument order matches forward");
    CHECK_NEAR(qn_fwd, qn_rev, 1e-12, "I6: backward argument order matches forward");
    CHECK_NEAR(qu_fwd, qu_rev, 1e-12, "I6: backward argument order matches forward");
    accel_prn_free(); unlink(p);
}

/* ---------- I7: property test against naive integrator ------------------ */
/* Naive Riemann sum, same piecewise-constant convention. */
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
        e[i]  = (rand()/(double)RAND_MAX) * 10.0;
        nn[i] = (rand()/(double)RAND_MAX) * 10.0;
        u[i]  = (rand()/(double)RAND_MAX) * 10.0;
    }
    char *p = write_csv_from(offsets, N, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();

    int mismatches = 0;
    /* Sample query intervals strictly inside the sample range (avoid edge
       fall-back so we can compare to the naive integral). */
    double t_first = offsets[0];
    double t_last  = offsets[N-1];
    for (int q = 0; q < Q; q++) {
        double a = t_first + (rand()/(double)RAND_MAX) * (t_last - t_first - 0.001);
        double b = a + 0.001 + (rand()/(double)RAND_MAX) * (t_last - a - 0.001);
        if (b > t_last) b = t_last;
        double qe, qn, qu;
        int ok = accel_prn_integrate(timeadd(t0, a), timeadd(t0, b),
                                     &qe, &qn, &qu);
        if (!ok) {
            mismatches++;
            fprintf(stderr, "  unexpected no-coverage: a=%.6f b=%.6f\n", a, b);
            continue;
        }
        double ee = naive_integrate_axis(offsets, e,  N, a, b);
        double nn_e = naive_integrate_axis(offsets, nn, N, a, b);
        double uu = naive_integrate_axis(offsets, u,  N, a, b);
        /* Relative tolerance: ~1e-7 is the floor we can reasonably hit
           given the gtime_t double-precision sec field and CSV string
           round-trip. Anything beyond that is irreducible numerical noise. */
        double rtol = 1e-7;
        if (fabs(qe - ee) > rtol*fabs(ee) + 1e-9 ||
            fabs(qn - nn_e) > rtol*fabs(nn_e) + 1e-9 ||
            fabs(qu - uu) > rtol*fabs(uu) + 1e-9) {
            mismatches++;
            fprintf(stderr, "  q=%d a=%.6f b=%.6f  qe=%.6f exp=%.6f rel=%.2e\n",
                    q, a, b, qe, ee, fabs(qe-ee)/fmax(fabs(ee),1e-12));
        }
    }
    CHECK(mismatches == 0, "I7: integrator must match naive Riemann sum");
    free(offsets); free(e); free(nn); free(u);
    accel_prn_free(); unlink(p);
}

/* ---------- B1: frame identity at lat=0, lon=0 (unchanged) ------------- */
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
    CHECK_NEAR(Qv[0], prn_u*prn_u, 1e-12, "B1: Qv[X,X] must equal Qu");
    CHECK_NEAR(Qv[4], prn_e*prn_e, 1e-12, "B1: Qv[Y,Y] must equal Qe");
    CHECK_NEAR(Qv[8], prn_n*prn_n, 1e-12, "B1: Qv[Z,Z] must equal Qn");
    for (int i = 0; i < 9; i++) {
        if (i == 0 || i == 4 || i == 8) continue;
        CHECK_NEAR(Qv[i], 0.0, 1e-12, "B1: off-diagonals zero");
    }
}

/* ---------- B3: trace + symmetry invariance (unchanged) ---------------- */
static void test_trace_invariance(void) {
    fprintf(stderr, "test_trace_invariance\n");
    srand(7);
    for (int trial = 0; trial < 50; trial++) {
        double lat = (-89.0 + 178.0 * (rand()/(double)RAND_MAX)) * D2R;
        double lon = (-179.0 + 358.0 * (rand()/(double)RAND_MAX)) * D2R;
        double pos[3] = {lat, lon, 0.0};
        double pe = rand()/(double)RAND_MAX;
        double pn = rand()/(double)RAND_MAX;
        double pu = rand()/(double)RAND_MAX;
        double Q[9] = {0}, Qv[9];
        Q[0] = pe*pe; Q[4] = pn*pn; Q[8] = pu*pu;
        covecef(pos, Q, Qv);
        double tr_q  = Q[0]  + Q[4]  + Q[8];
        double tr_qv = Qv[0] + Qv[4] + Qv[8];
        CHECK_NEAR(tr_qv, tr_q, 1e-10, "B3: rotation preserves trace");
        CHECK_NEAR(Qv[1], Qv[3], 1e-10, "B3: Qv symmetric (1,3)");
        CHECK_NEAR(Qv[2], Qv[6], 1e-10, "B3: Qv symmetric (2,6)");
        CHECK_NEAR(Qv[5], Qv[7], 1e-10, "B3: Qv symmetric (5,7)");
    }
}

/* ---------- B5: active-branch parity with stock formula ----------------- */
/* When the CSV contains constant prn over the integration interval, the
 * integrated variance equals the stock formula prn² · dt for that axis.
 * This proves the new integrator reduces to stock when σ is constant. */
static void test_active_branch_parity(void) {
    fprintf(stderr, "test_active_branch_parity\n");
    /* 10 samples at 0.0, 0.1, ..., 0.9 s with constant prn values. */
    int N = 10;
    double offsets[10], e[10], nn[10], u[10];
    double prnh = 1.7, prnv = 0.6;
    for (int i = 0; i < N; i++) {
        offsets[i] = i * 0.1;
        e[i] = prnh; nn[i] = prnh; u[i] = prnv;
    }
    char *p = write_csv_from(offsets, N, e, nn, u);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    /* Integrate over [0.05, 0.85] — a non-aligned interval inside the range. */
    double t_lo_s = 0.05, t_hi_s = 0.85;
    double qe, qn, qu;
    int ok = accel_prn_integrate(timeadd(t0, t_lo_s), timeadd(t0, t_hi_s),
                                 &qe, &qn, &qu);
    CHECK(ok == 1, "B5: should be covered");
    double dt = t_hi_s - t_lo_s;
    CHECK_NEAR(qe, prnh*prnh*dt, 1e-10, "B5: integrated qe == prnh^2*dt");
    CHECK_NEAR(qn, prnh*prnh*dt, 1e-10, "B5: integrated qn == prnh^2*dt");
    CHECK_NEAR(qu, prnv*prnv*dt, 1e-10, "B5: integrated qu == prnv^2*dt");
    accel_prn_free(); unlink(p);
}

int main(void) {
    traceopen("/dev/null");
    tracelevel(0);

    test_single_sample_window();
    test_multi_sample_window();
    test_zero_width();
    test_before_first();
    test_past_last();
    test_reversed_args();
    test_integrate_property();
    test_frame_at_origin();
    test_trace_invariance();
    test_active_branch_parity();

    traceclose();
    if (g_failed) { fprintf(stderr, "FAILED\n"); return EXIT_FAILURE; }
    fprintf(stderr, "OK\n");
    return EXIT_SUCCESS;
}
