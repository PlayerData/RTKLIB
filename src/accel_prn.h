/*------------------------------------------------------------------------------
* accel_prn.h : per-axis IMU-accel sidecar exposing a time-weighted mean
*               for loose-coupled IMU/GNSS measurement updates.
*
* The mean over [t_lo, t_hi] is fed to rtkpos.c as a 3-vector measurement of
* the EKF acceleration state in ENU. process_noise_accel_{h,v} (RTKLIB conf
* keys stats-prnaccelh-imu / stats-prnaccelv-imu) are the per-axis std of
* that measurement; random_walk_accel_{h,v} (stock stats-prnaccelh / -v)
* remain the random-walk-accel process noise applied unconditionally.
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

/* Time-weighted mean of (a_E, a_N, a_U) over [t_lo, t_hi] using a
 * piecewise-constant model: sample i is held over [t_i, t_{i+1}). Argument
 * order doesn't matter; the function normalises.
 *
 * Returns 1 with *ze, *zn, *zu set to the per-axis mean acceleration (m/s²)
 * when [t_lo, t_hi] is strictly covered by [first.t, last.t] and has
 * non-zero width. Returns 0 (outputs untouched) when no CSV is loaded,
 * t.time == 0 on either bound, the interval straddles either edge, or
 * the interval has zero width. */
int accel_prn_mean(gtime_t t_lo, gtime_t t_hi,
                   double *ze, double *zn, double *zu);

void accel_prn_free(void);

#ifdef __cplusplus
}
#endif

#endif
