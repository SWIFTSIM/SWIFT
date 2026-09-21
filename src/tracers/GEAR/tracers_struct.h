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

    /*! Cumulative |delta_p| per event (scalar sum, not vector: isotropic
        kicks would else cancel), physical. */
    float momentum_supernovae;
    float momentum_winds;
    float momentum_radiation;

    /*! Cumulative specific internal energy received (radiation pressure has
        no separate thermal channel). */
    float energy_supernovae;
    float energy_winds;

    /*! Largest single-event kick velocity received (outflow diagnostic). */
    float max_kick_velocity_supernovae;
    float max_kick_velocity_winds;
    float max_kick_velocity_radiation;

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
 * cumulative quantities are the budget the star ejected, counted once per
 * star per step, not a sum over the receiving gas.
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

  /*! Number of injection steps so far */
  int n_injection_steps;

  /*! Density at the most recent injection step, physical internal units */
  float density_at_last_injection;

  /*! Scale-factor (cosmological runs) or time (non-cosmological), of the
      most recent injection step */
  union {
    float last_injection_scale_factor;
    float last_injection_time;
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

  /*! Stellar-wind injection tracers */
  /* In the end, I don't think we need all these */
  struct tracers_winds_data winds;

  /* Two of the three radiation channels already have a tracer elsewhere, so
     only ISRF (photoelectric heating/LW dissociation) is untracked: HII has
     final_HII_radius/final_HII_mass above (star-side); radiation pressure
     has feedback_cumulative.momentum_radiation/max_kick_velocity_radiation
     in tracers_xpart_data above (gas-side, via
     tracers_after_radiation_pressure_feedback_part()). */
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
