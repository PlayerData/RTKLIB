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
* The C-side accel_prn is now coverage-only: the loaded CSV provides
* timestamps marking IMU sample points, and accel_prn_covered(t_lo, t_hi)
* returns 1 iff the interval is fully bracketed by [first.t, last.t].
* The prn-std values applied during/outside coverage live in the RTKLIB
* conf (prn_imu_acch/_accv vs prn[3]/[4]).
*
* Tests covered:
*   C1   covered: interval strictly inside sample range
*   C2   covered: interval at exact sample boundaries
*   C3   not covered: interval before first sample
*   C4   not covered: interval past last sample
*   C5   not covered: interval straddling first sample
*   C6   not covered: interval straddling last sample
*   C7   reversed argument order produces same result
*   C8   trailing-column tolerance: CSV with extra columns loads cleanly
*   B1   frame identity at lat=0, lon=0 (covecef axis assignment)
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

/* Write a 1-column-only CSV (just timestamps) at the given offsets. */
static char *write_times_csv(const double *offsets, int n) {
    static char path[64];
    snprintf(path, sizeof(path), "/tmp/t_accel_prn_t_%d.csv", (int)getpid());
    FILE *fp = fopen(path, "w");
    assert(fp);
    fputs("time\n", fp);
    gtime_t t0 = base_time();
    for (int i = 0; i < n; i++) {
        gtime_t t = timeadd(t0, offsets[i]);
        char ts[64]; time2str(t, ts, 9);
        ts[4]='-'; ts[7]='-'; ts[10]='T';
        fprintf(fp, "%sZ\n", ts);
    }
    fclose(fp);
    return path;
}

/* Write a 4-column CSV: time + dummy values that should be ignored. */
static char *write_4col_csv(const double *offsets, int n) {
    static char path[64];
    snprintf(path, sizeof(path), "/tmp/t_accel_prn_4_%d.csv", (int)getpid());
    FILE *fp = fopen(path, "w");
    assert(fp);
    fputs("time,ae_mps2,an_mps2,au_mps2\n", fp);
    gtime_t t0 = base_time();
    for (int i = 0; i < n; i++) {
        gtime_t t = timeadd(t0, offsets[i]);
        char ts[64]; time2str(t, ts, 9);
        ts[4]='-'; ts[7]='-'; ts[10]='T';
        fprintf(fp, "%sZ,1.5,2.5,3.5\n", ts);
    }
    fclose(fp);
    return path;
}

/* ---------- C1: interval strictly inside sample range ------------------ */
static void test_covered_inside(void) {
    fprintf(stderr, "test_covered_inside\n");
    double offsets[] = {0.0, 1.0, 2.0, 3.0};
    char *p = write_times_csv(offsets, 4);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int cov = accel_prn_covered(timeadd(t0, 0.5), timeadd(t0, 2.5));
    CHECK(cov == 1, "C1: strict-inside interval is covered");
    accel_prn_free(); unlink(p);
}

/* ---------- C2: interval at exact sample boundaries -------------------- */
static void test_covered_at_boundaries(void) {
    fprintf(stderr, "test_covered_at_boundaries\n");
    double offsets[] = {0.0, 1.0, 2.0, 3.0};
    char *p = write_times_csv(offsets, 4);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int cov = accel_prn_covered(t0, timeadd(t0, 3.0));
    CHECK(cov == 1, "C2: interval matching first..last sample is covered");
    accel_prn_free(); unlink(p);
}

/* ---------- C3: interval entirely before first sample ------------------ */
static void test_before_first(void) {
    fprintf(stderr, "test_before_first\n");
    double offsets[] = {1.0, 2.0, 3.0};
    char *p = write_times_csv(offsets, 3);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int cov = accel_prn_covered(timeadd(t0, -2.0), timeadd(t0, -1.0));
    CHECK(cov == 0, "C3: interval before first sample is not covered");
    accel_prn_free(); unlink(p);
}

/* ---------- C4: interval entirely past last sample --------------------- */
static void test_past_last(void) {
    fprintf(stderr, "test_past_last\n");
    double offsets[] = {0.0, 1.0, 2.0};
    char *p = write_times_csv(offsets, 3);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int cov = accel_prn_covered(timeadd(t0, 3.0), timeadd(t0, 4.0));
    CHECK(cov == 0, "C4: interval past last sample is not covered");
    accel_prn_free(); unlink(p);
}

/* ---------- C5: interval straddling first sample ----------------------- */
static void test_straddle_first(void) {
    fprintf(stderr, "test_straddle_first\n");
    double offsets[] = {1.0, 2.0, 3.0};
    char *p = write_times_csv(offsets, 3);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int cov = accel_prn_covered(timeadd(t0, 0.5), timeadd(t0, 2.0));
    CHECK(cov == 0, "C5: interval starting before first sample is not covered");
    accel_prn_free(); unlink(p);
}

/* ---------- C6: interval straddling last sample ------------------------ */
static void test_straddle_last(void) {
    fprintf(stderr, "test_straddle_last\n");
    double offsets[] = {0.0, 1.0, 2.0};
    char *p = write_times_csv(offsets, 3);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int cov = accel_prn_covered(timeadd(t0, 1.0), timeadd(t0, 3.0));
    CHECK(cov == 0, "C6: interval ending after last sample is not covered");
    accel_prn_free(); unlink(p);
}

/* ---------- C7: reversed argument order -------------------------------- */
static void test_reversed_args(void) {
    fprintf(stderr, "test_reversed_args\n");
    double offsets[] = {0.0, 1.0, 2.0, 3.0};
    char *p = write_times_csv(offsets, 4);
    accel_prn_load(p);
    gtime_t t0 = base_time();
    int fwd = accel_prn_covered(timeadd(t0, 0.5), timeadd(t0, 2.5));
    int rev = accel_prn_covered(timeadd(t0, 2.5), timeadd(t0, 0.5));
    CHECK(fwd == 1 && rev == 1,
          "C7: reversed-arg coverage matches forward");
    accel_prn_free(); unlink(p);
}

/* ---------- C8: 4-column CSV (extra cols ignored) ---------------------- */
static void test_trailing_columns(void) {
    fprintf(stderr, "test_trailing_columns\n");
    double offsets[] = {0.0, 1.0, 2.0};
    char *p = write_4col_csv(offsets, 3);
    int n = accel_prn_load(p);
    CHECK(n == 3, "C8: 4-column CSV loads (extra columns ignored)");
    gtime_t t0 = base_time();
    CHECK(accel_prn_covered(timeadd(t0, 0.5), timeadd(t0, 1.5)) == 1,
          "C8: coverage works on 4-column CSV");
    accel_prn_free(); unlink(p);
}

/* ---------- B1: frame identity at lat=0, lon=0 ------------------------- */
/* At (lat=0, lon=0): East=+Y_ecef, North=+Z_ecef, Up=+X_ecef.
   Q in ENU diagonal → covecef rotates to ECEF; the rotation is the load-
   bearing axis-correctness check that East ends up on Y, North on Z, Up
   on X regardless of any future udpos changes. */
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

/* ---------- B3: trace + symmetry invariance ---------------------------- */
static void test_trace_invariance(void) {
    fprintf(stderr, "test_trace_invariance\n");
    srand(7);
    for (int trial = 0; trial < 50; trial++) {
        double lat = (-89.0 + 178.0 * (rand()/(double)RAND_MAX)) * (3.14159265358979/180.0);
        double lon = (-179.0 + 358.0 * (rand()/(double)RAND_MAX)) * (3.14159265358979/180.0);
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

int main(void) {
    traceopen("/dev/null");
    tracelevel(0);

    test_covered_inside();
    test_covered_at_boundaries();
    test_before_first();
    test_past_last();
    test_straddle_first();
    test_straddle_last();
    test_reversed_args();
    test_trailing_columns();
    test_frame_at_origin();
    test_trace_invariance();

    traceclose();
    if (g_failed) { fprintf(stderr, "FAILED\n"); return EXIT_FAILURE; }
    fprintf(stderr, "OK\n");
    return EXIT_SUCCESS;
}
