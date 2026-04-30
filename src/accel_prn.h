/*------------------------------------------------------------------------------
* accel_prn.h : IMU-coverage signal for the EKF accel-state process noise.
*
* Loads a CSV listing the times at which IMU samples exist; the per-row
* values (if any) are ignored — only the time column is parsed. Used by
* udpos() to decide which pair of process-noise constants to apply at each
* epoch:
*   covered ([t_lo, t_hi] inside [t_first, t_last])  → opt.prn_imu_acch/v
*   not covered                                       → opt.prn[3] / prn[4]
*-----------------------------------------------------------------------------*/
#ifndef ACCEL_PRN_H
#define ACCEL_PRN_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Load a time-sorted CSV. Header row required; the first column is parsed
 * as ISO8601 (YYYY-MM-DDTHH:MM:SS[.fff][Z]) treated as GPST. Any further
 * columns are ignored — the file's role here is purely as a coverage
 * signal.
 *
 * Aborts (exit 1) on parse error or non-strictly-increasing timestamps.
 * Returns row count. NULL/empty path is a no-op (returns 0). */
int accel_prn_load(const char *path);

int accel_prn_loaded(void);

/* Return 1 if [t_lo, t_hi] (in either order) is fully bracketed by the
 * loaded sample timestamps — i.e. both endpoints fall in [first.t, last.t].
 * Return 0 when no CSV is loaded, when t.time == 0 on either bound, or
 * when the interval straddles either edge of the sample range. */
int accel_prn_covered(gtime_t t_lo, gtime_t t_hi);

void accel_prn_free(void);

#ifdef __cplusplus
}
#endif

#endif
