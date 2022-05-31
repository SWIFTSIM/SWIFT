/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_SINK_PART_H
#define SWIFT_DEFAULT_SINK_PART_H

#include "timeline.h"

#define sink_need_unique_id 1

/**
 * @brief Particle fields for the sink particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct sink {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v[3];

  /*! Cut off radius. */
  float r_cut;

  /*! Sink particle mass */
  float mass;

  /*! Particle time bin */
  timebin_t time_bin;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

  /*! sink merger information (e.g. merging ID) */
  struct sink_sink_data merger_data;

#ifdef DEBUG_INTERACTIONS_SINKS
  /*! Number of interactions in merger SELF and PAIR */
  int num_ngb_merger;

  /*! List of interacting particles in merger SELF and PAIR */
  long long ids_ngbs_merger[MAX_NUM_OF_NEIGHBOURS_SINKS];

  /*! Number of interactions in compute formation SELF and PAIR */
  int num_ngb_formation;

  /*! List of interacting particles in compute formation SELF and PAIR */
  long long ids_ngbs_formation[MAX_NUM_OF_NEIGHBOURS_SINKS];

  /*! Number of interactions in compute formation SELF and PAIR */
  int num_ngb_accretion;

  /*! List of interacting particles in compute formation SELF and PAIR */
  long long ids_ngbs_accretion[MAX_NUM_OF_NEIGHBOURS_SINKS];
#endif
} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DEFAULT_SINK_PART_H */
