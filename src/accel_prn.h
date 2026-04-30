/*------------------------------------------------------------------------------
* accel_prn.h : per-epoch override of EKF acceleration process-noise stds
*               from a measured-acceleration sidecar CSV.
*
* The CSV is loaded at full IMU rate. At each EKF epoch we integrate the
* per-axis variance contribution over the actual epoch interval [t_lo, t_hi]
* using a piecewise-constant model — sample i is held over [t_i, t_{i+1}).
* This is equivalent to a left-Riemann discretisation of the random-walk-
* accel process noise integral and removes any window-size or rate
* assumption: the math is correct regardless of GNSS / IMU rate.
*-----------------------------------------------------------------------------*/
#ifndef ACCEL_PRN_H
#define ACCEL_PRN_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Load a time-sorted CSV of per-axis measured process-noise stds.
 * Header row required: "time,prn_e,prn_n,prn_u" (case-insensitive).
 * Time column: ISO8601 (YYYY-MM-DDTHH:MM:SS[.fff][Z]) parsed as GPST,
 * matching our PPK pipeline convention (GPS time tagged Z).
 * Values: per-axis accel process-noise std in m/s^2 (local ENU frame).
 * Aborts (exit 1) on parse error or non-strictly-increasing timestamps.
 * Returns row count. NULL/empty path is a no-op (returns 0). */
int accel_prn_load(const char *path);

int accel_prn_loaded(void);

/* Integrate the per-axis variance contribution over [t_lo, t_hi] using the
 * piecewise-constant model where sample i is held over [t_i, t_{i+1}).
 *
 * Argument order doesn't matter — the function normalises to the smaller/
 * larger pair so it works for both forward and backward processing.
 *
 * Returns 1 with *qe, *qn, *qu set to the integrated variances (m^2/s^3
 * — i.e. already includes the dt factor; caller adds directly to Q[i,i]
 * with no further multiplication) when the interval is fully covered by
 * loaded samples. Returns 0 when no CSV is loaded, when the interval
 * straddles either edge of the sample range, or when t.time == 0 on
 * either bound; outputs are left untouched and the caller should fall
 * back to the existing config-defined opt.prn[3..4] formula. */
int accel_prn_integrate(gtime_t t_lo, gtime_t t_hi,
                        double *qe, double *qn, double *qu);

void accel_prn_free(void);

#ifdef __cplusplus
}
#endif

#endif
