/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_STRUCT_GEAR_H
#define SWIFT_FEEDBACK_STRUCT_GEAR_H

#include "chemistry_struct.h"

/**
 * @brief Feedback fields carried by each hydro particles
 */
struct feedback_part_data {};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! specific energy received from supernovae */
  float delta_u;

  /*! Momemtum received from a supernovae */
  float delta_p[3];

  /*! Indicator if the particle receives energy from SN specifically */
  char hit_by_SN;

  /*! Indicator if the particle receives energy from SW specifically */
  char hit_by_preSN;
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Is the star dead? */
  int is_dead;

  /*! Inverse of normalisation factor used for the enrichment. */
  float enrichment_weight;

  /*! Number of Ia supernovae */
  float number_snia;

  /*! Number of II supernovae */
  float number_snii;

  /*! Energy injected in the surrounding particles */
  float energy_ejected;

  /*! Total mass ejected by the supernovae */
  float mass_ejected;

  /*! Chemical composition of the mass ejected */
  double metal_mass_ejected[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Does the particle needs the feedback loop? */
  char will_do_feedback;

  /*! The relative velocity between the surrounding gas and the feedback emitter at the place of the emitter */
  float relative_velocity_gas[3];

  /*! The sound speed in the gas medium at the position of the feedback emitter */
  float sound_speed_gas;

  /*! The minimal smoothing length of the gas in the vicinity of the star */
  float minimal_h_gas;

  /*! The internal specific energy at the star localisation multiplied by the total mass in the kernel */
  float total_internal_energy_gas;

  /*! The kinetic energy of the gas at the star location */
  float total_kinetic_energy_gas;

  /*! Total mass gas in kernel */
  float total_gas_mass;

  /*! The total gas mass in the kernel */
  //float local_gas_mass;

  /*! Pre-SN data struct */
  struct {

    /*! Energy injected in the surrounding particles, needs to be double as the
     * energy is momentally passed as energy per unit time in erg/yr units and
     * is of order 10^40*/ /* TODO:change into float for memory but /!\ change the stellar_wind.c functions in accordance */
    double energy_ejected;

    /*! Energy injection rate, needed for timestep criterion */
    double energy_dot;
    /*! Energy injected in the surrounding particles */
    float energy_ejected;

    /*! Mass injected in the surrounding particles, needs to be double as the
     mass is currently in Msol units and can be of orders 10^-40 */ /* TODO: change into float for memory*/
    double mass_ejected;

    /*! Mass Loss rate, needed for timestep criterion*/
    double mass_dot;

  } preSN;
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
