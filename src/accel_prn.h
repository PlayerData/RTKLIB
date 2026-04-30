/*------------------------------------------------------------------------------
* accel_prn.h : per-axis IMU-accel sidecar for the EKF accel-state process
*               noise.
*
* Loads a CSV of per-IMU-sample nav-frame ENU acceleration values and lets
* udpos integrate Σᵢ aᵢ² · Δtᵢ over a GPS epoch interval. The udpos site
* combines that with the EKF config's process_noise_accel_{h,v} baseline:
*
*   covered:      Q[i,i] = process_noise_accel_axis² · dt
*                          + Σ a_axis_i² · Δt_i
*   not covered:  Q[i,i] = random_walk_accel_axis² · dt
*
* "Covered" means [t_lo, t_hi] is fully bracketed by [first.t, last.t] of
* the loaded CSV. Partial coverage falls back fully — split-mode would mix
* two semantics in one Q and is hard to reason about.
*-----------------------------------------------------------------------------*/
#ifndef ACCEL_PRN_H
#define ACCEL_PRN_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Load a time-sorted CSV with at least 4 columns (header required, names
 * not enforced):  time, a_E (m/s²), a_N (m/s²), a_U (m/s²).
 * Time is ISO8601 (YYYY-MM-DDTHH:MM:SS[.fff][Z]) treated as GPST.
 * Aborts (exit 1) on parse error or non-strictly-increasing timestamps.
 * Returns row count. NULL/empty path is a no-op (returns 0). */
int accel_prn_load(const char *path);

int accel_prn_loaded(void);

/* Integrate Σ a²(τ) dτ over [t_lo, t_hi] per axis using a piecewise-constant
 * model — sample i is held over [t_i, t_{i+1}). Argument order doesn't
 * matter; the function normalises.
 *
 * Returns 1 with *qe2, *qn2, *qu2 set to the per-axis integrals (units
 * m²/s³) when [t_lo, t_hi] is fully bracketed by [first.t, last.t].
 * Returns 0 when no CSV is loaded, when t.time == 0 on either bound, or
 * when the interval straddles either edge of the sample range; outputs
 * are left untouched and the caller falls back to random_walk_accel. */
int accel_prn_integrate(gtime_t t_lo, gtime_t t_hi,
                        double *qe2, double *qn2, double *qu2);

void accel_prn_free(void);

#ifdef __cplusplus
}
#endif

#endif
