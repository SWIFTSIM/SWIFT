/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_TIMESTEP_LIMITER_STRUCT_H
#define SWIFT_TIMESTEP_LIMITER_STRUCT_H

/**
 * @file src/chemistry_struct.h
 * @brief Branches between the different chemistry functions.
 */

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "timeline.h"

/**
 * @brief #part-carried quantities for the time-step limiter
 */
struct timestep_limiter_data {

  /* Need waking-up ? */
  timebin_t wakeup;

  /*! Minimal time-bin across all neighbours */
  timebin_t min_ngb_time_bin;

  /* Do we want this particle to be synched back on the time-line? */
  char to_be_synchronized;

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /* Weighted number of neighbours in the limiter loop */
  float n_limiter;

  /* Exact weighted number of neighbours in the limiter loop */
  float n_limiter_exact;

  /* Integer number of neighbours in the limiter loop */
  int N_limiter;

  /* Exact integer number of neighbours in the limiter loop */
  int N_limiter_exact;
#endif
};

#endif /* SWIFT_TIMESTEP_LIMITER_STRUCT_H */
