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

  /*! Has this particle been hit by SNII feedback? */
  char hit_by_SNII_feedback;

  /*! Has this particle been hit by AGN feedback? */
  char hit_by_AGN_feedback;
};

#endif /* SWIFT_TRACERS_STRUCT_EAGLE_H */
