/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STRUCT_H
#define SWIFT_RT_STRUCT_H

/**
 * @file src/rt_struct.h
 * @brief Branches between the different radiative transfer structs.
 */

/* Config parameters. */
#include "timeline.h"

#include <config.h>

/* Import the right RT struct definition */
#if defined(RT_NONE)
#include "./rt/none/rt_struct.h"
#elif defined(RT_DEBUG)
#include "./rt/debug/rt_struct.h"
#elif defined(RT_GEAR)
#include "./rt/GEAR/rt_struct.h"
#elif defined(RT_SPHM1RT)
#include "./rt/SPHM1RT/rt_struct.h"
#else
#error "Invalid choice of radiation scheme"
#endif

/* Define a struct to contain all RT sub-cycling related
 * timestepping variables here. These variables need to be
 * identical for every scheme and users should never touch
 * them anyway, so hide them here. */

#if defined(RT_NONE)

/*! empty placeholder for RT timestepping data. */
struct rt_timestepping_data {
  union {
    /*! Time-bin this particle uses for RT interactions */
    timebin_t time_bin;

    /*! Minimal time-bin across all neighbours */
    timebin_t min_ngb_time_bin;
  };
};

#else

/*! data relevant to the sub-cycle timestepping of parts. */
struct rt_timestepping_data {

  /*! Time-bin this particle uses for RT interactions */
  timebin_t time_bin;

  /*! Minimal time-bin across all neighbours */
  timebin_t min_ngb_time_bin;
};
#endif

#endif /* SWIFT_RT_STRUCT_H */
