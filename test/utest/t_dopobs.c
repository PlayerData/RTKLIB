/*------------------------------------------------------------------------------
 * rtklib unit test driver : doppler-observation DD measurement update
 *
 * Tests the ddres_dopobs() helper (in rtkpos.c) in isolation, by hand-
 * constructing minimal rtk_t/obsd_t/prcopt_t structs and exercising
 * the doppler-DD math directly. Mirrors the printf+assert style of
 * t_ppp.c — no test framework.
 *
 * Run after building: ./t_dopobs
 *-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "../../src/rtklib.h"

extern int test_ddres_dopobs(rtk_t *rtk, const obsd_t *obs, const double *x,
                             const int *sat, const double *e, const double *azel,
                             const double *freq, const int *iu, const int *ir,
                             int ns, const double *rs, const double *dts,
                             double *v, double *H, double *Ri, double *Rj,
                             int *vflg, int *nb, int *b, int nv);

/* Minimal rtk_t initialiser. We only populate the fields ddres_dopobs
 * actually reads: opt, rb (base ECEF), nx, ssat[].sys, errbuf. */
static void init_rtk(rtk_t *rtk, double base_x, double base_y, double base_z) {
    memset(rtk, 0, sizeof(rtk_t));
    rtk->nx = 9;  /* dynamics=on → 3 pos + 3 vel + 3 accel */
    rtk->rb[0] = base_x;
    rtk->rb[1] = base_y;
    rtk->rb[2] = base_z;
    rtk->opt.mode = PMODE_KINEMA;
    rtk->opt.dopobs = 1;
    rtk->opt.dynamics = 1;
    rtk->opt.nf = 1;       /* L1-only */
    rtk->opt.ionoopt = IONOOPT_OFF;
    rtk->opt.err[4] = 1.0;     /* doppler error 1 Hz */
    rtk->opt.maxinno[0] = 1e9; /* effectively no gating in tests */
    rtk->opt.maxinno[1] = 1e9;
    /* enable a few satellites; sys field per sat */
    for (int i = 0; i < MAXSAT; i++) {
        rtk->ssat[i].sys = satsys(i+1, NULL);
    }
}

/* Build a single sat's [pos(0..2), vel(3..5)] in rs, plus dts (clk, drift). */
static void put_sat(double *rs, double *dts, int slot,
                    double sx, double sy, double sz,
                    double vx, double vy, double vz,
                    double clk, double drift) {
    rs[slot*6+0] = sx; rs[slot*6+1] = sy; rs[slot*6+2] = sz;
    rs[slot*6+3] = vx; rs[slot*6+4] = vy; rs[slot*6+5] = vz;
    dts[slot*2+0] = clk;
    dts[slot*2+1] = drift;
}

/* Compute and write LOS unit vector e from rcv at rr to sat pos at rs[slot]. */
static void put_los(double *e, int slot, const double *rs, const double *rr) {
    double dx = rs[slot*6+0] - rr[0];
    double dy = rs[slot*6+1] - rr[1];
    double dz = rs[slot*6+2] - rr[2];
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    /* RTKLIB convention (geodist): e points from receiver TO satellite,
     * i.e. normalized (rs - rr). */
    e[slot*3+0] = dx/r;
    e[slot*3+1] = dy/r;
    e[slot*3+2] = dz/r;
}

/* Two GPS sats overhead, rover at origin, base 1km away. Sats with small
 * but nonzero Doppler so the function doesn't skip them (D==0 is the
 * "no doppler available" guard). No motion anywhere → DD residual = 0
 * if Doppler is fed back through the predicted-rate inverse correctly. */
void utest1_static_zero_residual(void) {
    rtk_t rtk;
    init_rtk(&rtk, 0, 1000, 0);  /* base offset 1km in Y */

    int ns = 2;
    int sat[2] = {1, 2};      /* GPS sats G01, G02 */
    int iu[2] = {0, 1};       /* rover obs slots */
    int ir[2] = {2, 3};       /* base obs slots */

    obsd_t obs[4];
    memset(obs, 0, sizeof(obs));
    /* Use 1e-9 Hz: nonzero so the D==0 skip doesn't fire, but small enough
     * the residual contribution is negligible. */
    obs[0].sat = 1; obs[0].D[0] = 1e-9f;
    obs[1].sat = 2; obs[1].D[0] = 1e-9f;
    obs[2].sat = 1; obs[2].D[0] = 1e-9f;
    obs[3].sat = 2; obs[3].D[0] = 1e-9f;

    double rs[4*6] = {0};
    double dts[4*2] = {0};
    /* sat with nonzero velocity so the norm() check passes. Tiny velocity
     * so it doesn't materially affect residual. */
    put_sat(rs, dts, 0,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 1, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 2,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 3, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);

    double rover[3] = {0, 0, 0};
    double base[3]  = {0, 1000, 0};
    double e[4*3] = {0};
    put_los(e, 0, rs, rover);
    put_los(e, 1, rs, rover);
    put_los(e, 2, rs, base);
    put_los(e, 3, rs, base);

    /* azel selection: ref sat = highest elevation. Force sat 0 to be ref
     * by giving it slightly higher elevation than sat 1. */
    double azel[4*2] = {0, 1.5, 0, 1.4, 0, 1.5, 0, 1.4};
    double freq[4*1] = {FREQL1, FREQL1, FREQL1, FREQL1};
    double x[9] = {0, 0, 0,  0, 0, 0,  0, 0, 0};

    double v[8] = {0}, H[9*8] = {0}, Ri[8] = {0}, Rj[8] = {0};
    int vflg[8] = {0}, nb[10] = {0}, b = 0;

    int nv = test_ddres_dopobs(&rtk, obs, x, sat, e, azel, freq, iu, ir,
                               ns, rs, dts, v, H, Ri, Rj, vflg, nb, &b, 0);

    printf("utest1: nv=%d, b=%d, v[0]=%.6f\n", nv, b, v[0]);
    assert(nv == 1);  /* one DD per non-ref sat in one system */
    assert(b == 1);
    assert(fabs(v[0]) < 1e-3);

    /* iref = sat 0, non-ref = sat 1.
     * H row: position partials zero; velocity partials = -e_ref + e_other.
     * Same sign convention as RTKLIB's phase/code DD: H[k] is the partial
     * of the *predicted* measurement wrt state (not the residual). The
     * Kalman update is xp = x + K·v with v = measured - predicted, so K·v
     * pushes the state in the right direction for either sign convention
     * provided H matches the convention.
     * H accel partials zero. */
    printf("  H[0..8]=");
    for (int i = 0; i < 9; i++) printf("%.4f ", H[i]);
    printf("\n");
    assert(fabs(H[0]) < 1e-9);
    assert(fabs(H[1]) < 1e-9);
    assert(fabs(H[2]) < 1e-9);
    assert(fabs(H[3] - (-e[0*3+0] + e[1*3+0])) < 1e-9);
    assert(fabs(H[4] - (-e[0*3+1] + e[1*3+1])) < 1e-9);
    assert(fabs(H[5] - (-e[0*3+2] + e[1*3+2])) < 1e-9);
    assert(fabs(H[6]) < 1e-9);
    assert(fabs(H[7]) < 1e-9);
    assert(fabs(H[8]) < 1e-9);

    printf("utest1 OK\n\n");
}

/* Sat moving directly toward (vertically down). Rover and base both static
 * at ground level. Doppler should be the same on rover and base (short
 * baseline, high sat). The DD measured = 0 (both sats produce same SD).
 * Predicted should also be 0 if math is right. Residual ≈ 0. */
void utest2_radial_sat_velocity_no_residual(void) {
    rtk_t rtk;
    init_rtk(&rtk, 0, 1000, 0);

    int ns = 2;
    int sat[2] = {1, 2};
    int iu[2] = {0, 1};
    int ir[2] = {2, 3};

    /* sat 1 falling vertically at 1000 m/s. sat 2 same. Doppler will reflect
     * +1000 m/s closing rate ≈ + (1000 / λ_L1) Hz Doppler (positive=approaching).
     * Since meas = -D·c/f, rrate = -1000·c·c/(c·f·f) … wait no:
     * rrate (m/s) = -D · c/f, with +D = approaching = +Doppler-shift = increasing freq.
     * Sat approaching at 1000 m/s → D = +1000·f/c, so rrate = -1000 m/s.
     * That matches: positive rate = range increasing; sat approaching = range decreasing. */
    double doppler_hz = 1000.0 * FREQL1 / CLIGHT;  /* approaching */

    obsd_t obs[4];
    memset(obs, 0, sizeof(obs));
    obs[0].sat = 1; obs[0].D[0] = (float)doppler_hz;
    obs[1].sat = 2; obs[1].D[0] = (float)doppler_hz;
    obs[2].sat = 1; obs[2].D[0] = (float)doppler_hz;
    obs[3].sat = 2; obs[3].D[0] = (float)doppler_hz;

    /* sats falling straight down at 1000 m/s. v_z = -1000. */
    double rs[4*6] = {0}, dts[4*2] = {0};
    put_sat(rs, dts, 0,  3000000, 0, 19000000, 0, 0, -1000, 0, 0);
    put_sat(rs, dts, 1, -3000000, 0, 19000000, 0, 0, -1000, 0, 0);
    put_sat(rs, dts, 2,  3000000, 0, 19000000, 0, 0, -1000, 0, 0);
    put_sat(rs, dts, 3, -3000000, 0, 19000000, 0, 0, -1000, 0, 0);

    double rover[3] = {0, 0, 0};
    double base[3]  = {0, 1000, 0};
    double e[4*3] = {0};
    put_los(e, 0, rs, rover);
    put_los(e, 1, rs, rover);
    put_los(e, 2, rs, base);
    put_los(e, 3, rs, base);

    double azel[4*2] = {0, 1.5, 0, 1.4, 0, 1.5, 0, 1.4};
    double freq[4*1] = {FREQL1, FREQL1, FREQL1, FREQL1};

    double x[9] = {0, 0, 0,  0, 0, 0,  0, 0, 0};

    double v[8] = {0}, H[9*8] = {0}, Ri[8] = {0}, Rj[8] = {0};
    int vflg[8] = {0}, nb[10] = {0}, b = 0;

    int nv = test_ddres_dopobs(&rtk, obs, x, sat, e, azel, freq, iu, ir,
                               ns, rs, dts, v, H, Ri, Rj, vflg, nb, &b, 0);

    printf("utest2: nv=%d v[0]=%.6f (expect small — sat radial motion cancels in DD)\n",
           nv, v[0]);
    assert(nv == 1);
    /* with both sats radially moving identically, DD should be near zero.
     * Looser tolerance because LOS isn't perfectly radial. */
    assert(fabs(v[0]) < 0.5);

    printf("utest2 OK\n\n");
}

/* Rover moving at known velocity. Doppler measured matches the geometry.
 * Predicted velocity in state matches actual rover velocity → residual = 0. */
void utest3_rover_velocity_correct_prediction(void) {
    rtk_t rtk;
    init_rtk(&rtk, 0, 1000, 0);

    int ns = 2;
    int sat[2] = {1, 2};
    int iu[2] = {0, 1};
    int ir[2] = {2, 3};

    /* Two sats. Rover moving at +10 m/s in X. */
    double v_rover[3] = {10, 0, 0};

    obsd_t obs[4];
    memset(obs, 0, sizeof(obs));

    /* sat 1 at (3M, 0, 19M); sat 2 at (-3M, 0, 19M). Both static.
     * Rover at origin moving +10 X. LOS vectors:
     *   to sat 1: (3M, 0, 19M)/r → mostly +X with some +Z
     *   to sat 2: (-3M, 0, 19M)/r → mostly -X with some +Z
     * Range rate = e·(v_sat - v_rcv). With v_sat=0, range rate = -e·v_rover.
     *   sat 1: rate ≈ -(0.156)·10 ≈ -1.56 m/s (closing)
     *   sat 2: rate ≈ -(-0.156)·10 ≈ +1.56 m/s (opening)
     * Doppler (Hz, RTKLIB convention) = -rate · f/c
     *   sat 1: D = -(-1.56)·f/c = +1.56·f/c
     *   sat 2: D = -(+1.56)·f/c = -1.56·f/c
     * Sat velocity must be nonzero (norm() guard); use tiny radial value
     * so it doesn't materially affect the residual. */
    double rs[4*6] = {0}, dts[4*2] = {0};
    put_sat(rs, dts, 0,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 1, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 2,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 3, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);

    double rover[3] = {0, 0, 0};
    double base[3]  = {0, 1000, 0};
    double e[4*3] = {0};
    put_los(e, 0, rs, rover);
    put_los(e, 1, rs, rover);
    put_los(e, 2, rs, base);
    put_los(e, 3, rs, base);

    /* compute expected range rates per sat */
    double rate0 = -(e[0*3+0]*v_rover[0] + e[0*3+1]*v_rover[1] + e[0*3+2]*v_rover[2]);
    double rate1 = -(e[1*3+0]*v_rover[0] + e[1*3+1]*v_rover[1] + e[1*3+2]*v_rover[2]);
    /* rrate_meas = -D·c/f, so D_rover = -rrate · f/c */
    double D0 = -rate0 * FREQL1/CLIGHT;
    double D1 = -rate1 * FREQL1/CLIGHT;
    /* base is static, sats static → base sees zero range rate.
     * The function uses D==0 as a "no obs" sentinel and skips, so we use
     * a tiny stand-in (1e-9 Hz). */
    obs[0].sat = 1; obs[0].D[0] = (float)D0;
    obs[1].sat = 2; obs[1].D[0] = (float)D1;
    obs[2].sat = 1; obs[2].D[0] = 1e-9f;
    obs[3].sat = 2; obs[3].D[0] = 1e-9f;

    double azel[4*2] = {0, 1.5, 0, 1.4, 0, 1.5, 0, 1.4};
    double freq[4*1] = {FREQL1, FREQL1, FREQL1, FREQL1};

    /* state has rover velocity correctly populated */
    double x[9] = {0, 0, 0,  v_rover[0], v_rover[1], v_rover[2],  0, 0, 0};

    double v[8] = {0}, H[9*8] = {0}, Ri[8] = {0}, Rj[8] = {0};
    int vflg[8] = {0}, nb[10] = {0}, b = 0;

    int nv = test_ddres_dopobs(&rtk, obs, x, sat, e, azel, freq, iu, ir,
                               ns, rs, dts, v, H, Ri, Rj, vflg, nb, &b, 0);

    printf("utest3: nv=%d v[0]=%.6f (expect small — prediction matches measurement)\n",
           nv, v[0]);
    assert(nv == 1);
    /* Tolerance allows for the ~1mm/s earth rotation correction we don't
     * fully cancel between rover/base for non-perfectly-radial sats. */
    assert(fabs(v[0]) < 0.05);

    printf("utest3 OK\n\n");
}

/* Same as utest3 but with state holding WRONG rover velocity. The DD
 * residual should reflect the magnitude of the velocity error. */
void utest4_rover_velocity_wrong_prediction(void) {
    rtk_t rtk;
    init_rtk(&rtk, 0, 1000, 0);

    int ns = 2;
    int sat[2] = {1, 2};
    int iu[2] = {0, 1};
    int ir[2] = {2, 3};

    double v_actual[3] = {10, 0, 0};
    double v_state[3]  = {0, 0, 0};  /* state thinks rover is stationary */

    double rs[4*6] = {0}, dts[4*2] = {0};
    /* tiny nonzero sat velocity to pass the norm() guard */
    put_sat(rs, dts, 0,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 1, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 2,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 3, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);

    double rover[3] = {0, 0, 0};
    double base[3]  = {0, 1000, 0};
    double e[4*3] = {0};
    put_los(e, 0, rs, rover);
    put_los(e, 1, rs, rover);
    put_los(e, 2, rs, base);
    put_los(e, 3, rs, base);

    double rate0 = -(e[0*3+0]*v_actual[0] + e[0*3+1]*v_actual[1] + e[0*3+2]*v_actual[2]);
    double rate1 = -(e[1*3+0]*v_actual[0] + e[1*3+1]*v_actual[1] + e[1*3+2]*v_actual[2]);
    double D0 = -rate0 * FREQL1/CLIGHT;
    double D1 = -rate1 * FREQL1/CLIGHT;

    obsd_t obs[4];
    memset(obs, 0, sizeof(obs));
    obs[0].sat = 1; obs[0].D[0] = (float)D0;
    obs[1].sat = 2; obs[1].D[0] = (float)D1;
    obs[2].sat = 1; obs[2].D[0] = 1e-9f;
    obs[3].sat = 2; obs[3].D[0] = 1e-9f;

    double azel[4*2] = {0, 1.5, 0, 1.4, 0, 1.5, 0, 1.4};
    double freq[4*1] = {FREQL1, FREQL1, FREQL1, FREQL1};

    /* state thinks rover is stationary — but measurements come from moving rover */
    double x[9] = {0, 0, 0,  v_state[0], v_state[1], v_state[2],  0, 0, 0};

    double v[8] = {0}, H[9*8] = {0}, Ri[8] = {0}, Rj[8] = {0};
    int vflg[8] = {0}, nb[10] = {0}, b = 0;

    int nv = test_ddres_dopobs(&rtk, obs, x, sat, e, azel, freq, iu, ir,
                               ns, rs, dts, v, H, Ri, Rj, vflg, nb, &b, 0);

    /* expected DD residual: DD_meas - DD_pred
     *   DD_meas = (rate0 - 0) - (rate1 - 0) = rate0 - rate1
     *   DD_pred (v_state=0) = 0 - 0 - (0 - 0) = 0
     * so v[0] = (rate0 - rate1) - 0 = rate0 - rate1 */
    double expected = rate0 - rate1;
    printf("utest3: nv=%d v[0]=%.6f (expect %.6f = rate0-rate1)\n",
           nv, v[0], expected);
    assert(nv == 1);
    assert(fabs(v[0] - expected) < 0.05);

    /* H[3..5] should be -e_ref + e_other (RTKLIB H is ∂h/∂x, not ∂v/∂x). */
    assert(fabs(H[3] - (-e[0*3+0] + e[1*3+0])) < 1e-9);
    assert(fabs(H[4] - (-e[0*3+1] + e[1*3+1])) < 1e-9);
    assert(fabs(H[5] - (-e[0*3+2] + e[1*3+2])) < 1e-9);

    printf("utest4 OK\n\n");
}

/* Critical test: the H matrix sign must drive the EKF state in the
 * direction that REDUCES the residual, not increases it. If we have a
 * residual v from velocity error δv, then applying a small velocity
 * correction in the direction H'·v should produce a smaller |v|.
 *
 * This catches the most likely class of bug — a sign flip in the H
 * partials — that would otherwise pass numeric-magnitude assertions.
 *
 * Concretely: sets up rover with TRUE velocity (10, 0, 0). Filter state
 * starts at v=(0,0,0). Computes residual v_old, applies a small step
 * α·H[3..5]·v_old to the state velocity, then re-runs the function. The
 * residual should DECREASE if H has the right sign. */
void utest5_h_sign_drives_state_toward_truth(void) {
    rtk_t rtk;
    init_rtk(&rtk, 0, 1000, 0);

    int ns = 2;
    int sat[2] = {1, 2};
    int iu[2] = {0, 1};
    int ir[2] = {2, 3};

    double v_actual[3] = {10, 0, 0};

    double rs[4*6] = {0}, dts[4*2] = {0};
    put_sat(rs, dts, 0,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 1, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 2,  3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);
    put_sat(rs, dts, 3, -3000000, 0, 19000000, 1e-9, 0, 0, 0, 0);

    double rover[3] = {0, 0, 0};
    double base[3]  = {0, 1000, 0};
    double e[4*3] = {0};
    put_los(e, 0, rs, rover);
    put_los(e, 1, rs, rover);
    put_los(e, 2, rs, base);
    put_los(e, 3, rs, base);

    double rate0 = -(e[0*3+0]*v_actual[0] + e[0*3+1]*v_actual[1] + e[0*3+2]*v_actual[2]);
    double rate1 = -(e[1*3+0]*v_actual[0] + e[1*3+1]*v_actual[1] + e[1*3+2]*v_actual[2]);
    double D0 = -rate0 * FREQL1/CLIGHT;
    double D1 = -rate1 * FREQL1/CLIGHT;

    obsd_t obs[4];
    memset(obs, 0, sizeof(obs));
    obs[0].sat = 1; obs[0].D[0] = (float)D0;
    obs[1].sat = 2; obs[1].D[0] = (float)D1;
    obs[2].sat = 1; obs[2].D[0] = 1e-9f;
    obs[3].sat = 2; obs[3].D[0] = 1e-9f;

    double azel[4*2] = {0, 1.5, 0, 1.4, 0, 1.5, 0, 1.4};
    double freq[4*1] = {FREQL1, FREQL1, FREQL1, FREQL1};

    /* state starts with WRONG velocity */
    double x[9] = {0, 0, 0,  0, 0, 0,  0, 0, 0};

    double v[8] = {0}, H[9*8] = {0}, Ri[8] = {0}, Rj[8] = {0};
    int vflg[8] = {0}, nb[10] = {0}, b = 0;

    int nv = test_ddres_dopobs(&rtk, obs, x, sat, e, azel, freq, iu, ir,
                               ns, rs, dts, v, H, Ri, Rj, vflg, nb, &b, 0);
    assert(nv == 1);
    double v_old = v[0];
    double Hx_old = H[3], Hy_old = H[4], Hz_old = H[5];

    /* Apply a small step in the direction the EKF would update.
     *
     * RTKLIB's filter() in rtkcmn.c stores H as ∂(predicted_meas)/∂x
     * (NOT ∂(residual)/∂x), and applies xp = x + K·v with
     * K = P·H·(H'PH + R)^{-1} and v = measured - predicted.
     *
     * For our 1-observation case, K ∝ P·H, so δx ∝ +H·v. With the RTKLIB
     * H sign convention (h-partial), the residual v decreases when we
     * step in the +H·v direction.
     *
     * Note: a naive textbook "δx = -∂h/∂x · residual" approach would have
     * the opposite sign. The convention difference matters; this test
     * pins the RTKLIB convention. */
    double scale = 0.01;  /* small step to stay near linearisation */
    double dvx = scale * Hx_old * v_old;
    double dvy = scale * Hy_old * v_old;
    double dvz = scale * Hz_old * v_old;
    x[3] += dvx;
    x[4] += dvy;
    x[5] += dvz;

    /* re-run with updated state */
    double v2[8] = {0}, H2[9*8] = {0}, Ri2[8] = {0}, Rj2[8] = {0};
    int vflg2[8] = {0}, nb2[10] = {0}, b2 = 0;
    nv = test_ddres_dopobs(&rtk, obs, x, sat, e, azel, freq, iu, ir,
                           ns, rs, dts, v2, H2, Ri2, Rj2, vflg2, nb2, &b2, 0);

    printf("utest5: v_old=%.6f, v_new=%.6f (after dv=(%.4f,%.4f,%.4f))\n",
           v_old, v2[0], dvx, dvy, dvz);
    printf("        |v_new| should be < |v_old| if H sign is correct\n");
    assert(fabs(v2[0]) < fabs(v_old));

    printf("utest5 OK\n\n");
}

int main(void) {
    utest1_static_zero_residual();
    utest2_radial_sat_velocity_no_residual();
    utest3_rover_velocity_correct_prediction();
    utest4_rover_velocity_wrong_prediction();
    utest5_h_sign_drives_state_toward_truth();
    printf("all tests passed\n");
    return 0;
}
