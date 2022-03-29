//
// Created by yuyttenh on 29/03/22.
//

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "space_getsid.h"
#include "timers.h"

#ifndef SWIFTSIM_RUNNER_DOIACT_GRID_HYDRO_H
#define SWIFTSIM_RUNNER_DOIACT_GRID_HYDRO_H

__attribute__((always_inline)) INLINE static void
runner_dopair_grid_flux_exchange(struct runner *restrict r,
                                 struct cell *restrict ci,
                                 struct cell *restrict cj) {
  /* TODO */
}

__attribute__((always_inline)) INLINE static void
runner_doself_grid_flux_exchange(struct runner *restrict r,
                                 struct cell *restrict ci) {
  /* TODO */
}

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_HYDRO_H
