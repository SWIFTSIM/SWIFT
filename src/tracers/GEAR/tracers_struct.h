/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_TRACERS_STRUCT_NONE_H
#define SWIFT_TRACERS_STRUCT_NONE_H

/* Local includes */
#include "tracers_triggers.h"

/**
 * @brief Properties of the tracers stored in the extended particle data.
 */
struct tracers_xpart_data {

  /*! Radiation struct */
  struct {

    /*! Tag to ionize the part */
    char is_ionized;

    /*! Id of the star that ionized this particle */
    long long star_id;

    /*! Time when this particle is not ionized anymore */
    float end_time;

  } HII_region;
};

/* /\** */
/*  * @brief Properties of the tracers stored in the star particle data. */
/*  * */
/*  * Note: In this model, they are identical to the xpart data. */
/*  *\/ */
/* #define tracers_spart_data tracers_xpart_data */
struct tracers_spart_data {};

/**
 * @brief Properties of the tracers stored in the black hole particle data.
 */
struct tracers_bpart_data {};

/**
 * @brief Properties of the tracers stored in the sink particle data.
 */
struct tracers_sink_data {

  /*! Averaged SFR over N different time slices */
  float averaged_SFR[num_snapshot_triggers_sink];

  /*! Averaged accretion rate over N different time slices */
  float averaged_accretion_rate[num_snapshot_triggers_sink];
};

#endif /* SWIFT_TRACERS_STRUCT_NONE_H */
