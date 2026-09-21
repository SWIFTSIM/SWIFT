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

  /*! Feedback received over this particle's whole lifetime, physical internal
   * units. */
  struct {

    /*! Cumulative |delta_p| per event (scalar sum, not vector: isotropic
        kicks would else cancel), physical. */
    float momentum_supernovae;
    float momentum_winds;

    /*! Cumulative specific internal energy received. */
    float energy_supernovae;
    float energy_winds;

    /*! Largest single-event kick velocity received (outflow diagnostic). */
    float max_kick_velocity_supernovae;
    float max_kick_velocity_winds;

  } feedback_cumulative;
};

/**
 * @brief Per-channel record of a star's own SN events over its lifetime.
 *
 * One event for a discrete star; for a population particle, "event" means
 * an active SN step (possibly fractional), not one discrete explosion.
 */
struct tracers_sn_event_data {

  /*! Number of events so far (fractional for a continuously-sampled
      population particle, matching feedback_data.number_snii/snia's own
      type) */
  float n_events;

  /*! Density at the most recent event, physical internal units. Whether
      this was the star's only event since the last snapshot is readable
      from n_events itself: compare it against the previous snapshot's
      value for the same star. */
  float density_at_last_event;

  /*! Scale-factor (cosmological runs) or time (non-cosmological), of the
      most recent event */
  union {
    float last_event_scale_factor;
    float last_event_time;
  };
};

/**
 * @brief Properties of the tracers stored in the star particle data.
 *
 */
struct tracers_spart_data {

  /*! SN event tracers, one per channel */
  struct tracers_sn_event_data snii_events;
  struct tracers_sn_event_data snia_events;
};

/**
 * @brief Properties of the tracers stored in the black hole particle data.
 */
struct tracers_bpart_data {

  /*! Averaged accretion rate over two different time slices */
  float averaged_accretion_rate[num_snapshot_triggers_bpart];
};

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
