/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2021 Edo Altamura (edoardo.altamura@manchester.ac.uk)
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
#ifndef SWIFT_TRACERS_STRUCT_EAGLE_H
#define SWIFT_TRACERS_STRUCT_EAGLE_H

/* Local includes */
#include "tracers_triggers.h"

/*! The possible accretion modes every black hole can take. */
enum BH_accretion_modes {
  BH_thick_disc = 0,       /* At low Eddington ratios */
  BH_thin_disc,            /* At moderate Eddington ratios */
  BH_slim_disc,            /* Super-Eddington accretion */
  BH_accretion_modes_count /* Number of possible accretion modes */
};

/**
 * @brief Properties of the tracers stored in the extended particle data.
 */
struct tracers_xpart_data {

  /*! Maximum temperature achieved by this particle */
  float maximum_temperature;

  /*! Anonymous union for the cosmological non-cosmological runs distinction */
  union {

    /*! Scale-factor at which the maximal temperature was reached */
    float maximum_temperature_scale_factor;

    /*! Time at which the maximal temperature was reached */
    float maximum_temperature_time;
  };

  union {

    /*! Scale-factor at which the particle last received energy from AGN */
    float last_AGN_injection_scale_factor;

    /*! Time at which the particle last received energy from AGN */
    float last_AGN_injection_time;
  };

  union {

    /* Last scale factor this particle was kicked as part of jet feedback */
    float last_AGN_jet_feedback_scale_factor;

    /* Last time this particle was kicked as part of jet feedback */
    float last_AGN_jet_feedback_time;
  };

  /*! Averaged SFR over two different time slices */
  float averaged_SFR[num_snapshot_triggers_part];

  /*! Density of the gas before the last AGN feedback event
   * (physical internal units) */
  float density_before_last_AGN_feedback_event;

  /*! Entropy of the gas before the last AGN feedback event
   * (physical internal units) */
  float entropy_before_last_AGN_feedback_event;

  /*! Density of the gas at the last AGN feedback event
   * (physical internal units) */
  float density_at_last_AGN_feedback_event;

  /*! Entropy of the gas at the last AGN feedback event
   * (physical internal units) */
  float entropy_at_last_AGN_feedback_event;

  /*! Total amount of AGN feedback energy received by this particle
   * (physical units) */
  float AGN_feedback_energy;

  /*! Total jet feedback energy received by this particle */
  float jet_feedback_energy;

  /*! Counts how many times this particle has been kicked as part of
     jet feedback */
  char hit_by_jet_feedback;

  /*! Has this particle been hit by SNII feedback? */
  char hit_by_SNII_feedback;

  /*! Has this particle been hit by AGN feedback? */
  char hit_by_AGN_feedback;

  /*! Kick velocity at last AGN jet event */
  float last_jet_kick_velocity;

  /*! The accretion/feedback mode of the BH when this particle was last
   * kicked */
  enum BH_accretion_modes last_jet_kick_accretion_mode;

  /*! The ID of the BH that did the last kick */
  long long last_jet_kick_BH_id;
};

/**
 * @brief Properties of the tracers stored in the star particle data.
 *
 * Note: In this model, they are identical to the xpart data.
 */
#define tracers_spart_data tracers_xpart_data

/**
 * @brief Properties of the tracers stored in the black hole particle data.
 */
struct tracers_bpart_data {

  /*! Averaged accretion rate over two different time slices */
  float averaged_accretion_rate[num_snapshot_triggers_bpart];
};

#endif /* SWIFT_TRACERS_STRUCT_EAGLE_H */
