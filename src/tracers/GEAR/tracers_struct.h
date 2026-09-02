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

  /*! Radiation struct. The tag core (is_ionized/star_id/end_time) lives on
      struct part's feedback_data instead (see
      src/feedback/GEAR_thermal/feedback_struct.h) for automatic
      MPI/restart coverage; only the owner-computed payload below stays
      here. */
  struct {

    /*! Mean photon energy above the 13.6 eV HI ionization threshold for
        the tagging star, frozen at tag time (only set when
        GEARFeedback:HII_couple_ionization_rate is on; 0 otherwise). Stored
        in cgs (erg), not internal units, since the absolute per-particle
        value underflows float precision in this project's internal unit
        system. */
    float excess_photon_energy_HI;

    /*! Photoionization rate coefficient Gamma_HI from the tagging star at
        this particle's location, frozen at tag time (internal 1/time;
        only set when GEARFeedback:HII_couple_ionization_rate is on, 0
        otherwise). */
    float photoionization_rate_HI;

  } HII_region;

  /*! Feedback received over this particle's whole lifetime, physical internal
   * units. */
  struct {

    /*! Cumulative momentum magnitude received from SN, stellar winds and
        radiation pressure. momentum_SN is comoving-frame (the SN branch in
        feedback_iact.h doesn't convert to physical velocities the way the
        winds branch does); exact for non-cosmological runs, revisit before
        trusting in a cosmological one. */
    float momentum_SN;
    float momentum_winds;
    float momentum_radiation;

    /*! Cumulative specific internal energy received from SN and stellar
        winds (radiation pressure carries no separate thermal channel) */
    float energy_SN;
    float energy_winds;

    /*! Largest single-event kick velocity received from each channel
        (outflow diagnostic; the peak coupling velocity near the source,
        before deceleration) */
    float max_kick_velocity_SN;
    float max_kick_velocity_winds;
    float max_kick_velocity_radiation;

  } feedback_cumulative;
};

/**
 * @brief Per-channel record of a star's own SN events over its lifetime.
 *
 * One event for a discrete (single_star) particle, possibly many for a
 * population particle (star_population/star_population_continuous_IMF)
 * sampled over its life. For a population particle, "event" means an
 * active step with nonzero fractional SN count, not a single discrete
 * explosion: the model spreads its SN injection continuously over an
 * extended window, so density_at_last_event there reads more like "density
 * at the star's last SN-active step" than "density at one specific blast".
 * Density is the star's own kernel-averaged local gas density
 * (feedback_data.enrichment_weight, comoving, converted to physical here),
 * not any one neighbour's.
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
  /*! Radius of the HII region before the star died or was not HII eligible
      for the rest of its lifetime */
  float final_HII_radius;

  /*! Ionized gas mass of that same final HII region */
  float final_HII_mass;

  /*! SN event tracers, one per channel */
  struct tracers_sn_event_data snii_events;
  struct tracers_sn_event_data snia_events;
};

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
