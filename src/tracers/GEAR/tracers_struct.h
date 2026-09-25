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
 * @brief Record of a star's own stellar-wind budget over its lifetime.
 *
 * Winds are continuous, so the unit of record is an injection step: a
 * timestep whose wind budget the star hands to its gas neighbours. All
 * quantities are the budget the star ejected, summed once per star per
 * step, not a sum over the receiving gas. They exist to form the star-side
 * versus gas-side energy and momentum ratios.
 */
struct tracers_winds_data {

  /*! Cumulative mass ejected by winds (internal units). */
  double mass_ejected;

  /*! Cumulative wind energy ejected, after the winds efficiency factor
      (physical internal units). */
  double energy_ejected;

  /*! Cumulative wind momentum magnitude, sqrt(2 m_ej E_ej) per step, in the
      star's rest frame (physical internal units). */
  double momentum_ejected;
};

/**
 * @brief Properties of the tracers stored in the star particle data.
 *
 */
struct tracers_spart_data {
  /*! SN event tracers, one per channel */
  struct tracers_sn_event_data snii_events;
  struct tracers_sn_event_data snia_events;

  /*! Stellar-wind ejecta budget */
  struct tracers_winds_data winds;

  /* None of the three radiation channels is tracked here: HII's final
     extent lives in feedback_spart_data.radiation (star-side), radiation
     pressure's cumulative momentum/kick velocity live in
     feedback_xpart_data.radiation (gas-side, feedback_struct.h), and ISRF
     (photoelectric heating/LW dissociation) has no tracer at all. */
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
