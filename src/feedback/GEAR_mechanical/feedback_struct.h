/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H
#define SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H

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

  /*! Kinetic energy (not specific!) received from supernovae */
  float delta_E_kin;

  /*! Momemtum received from a supernovae */
  float delta_p[3];

  /* Number of supernovae affecting this particle */
  int number_SN;

};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

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

  /* Parameters to be accumulated in the feedback loops. Used to compute the
     vector weights */
  double f_plus_num[3];
  double f_plus_denom[3];
  double f_minus_num[3];
  double f_minus_denom[3];

  /* Accumulated value for the total energy available in the SN, taking into
     account gas-star motion */
  double E_total_accumulator;

  /* Parameters to determine the coupled energy, momentum and internal energy
     of the SN */
  double beta_1_accumulator; /* accumulated value for beta_1 */
  double beta_2_accumulator; /* accumulated value for beta_2 */

  /* Sum of the gas properties used to compute the mean of the properties */
  double sum_gas_density;
  double sum_gas_metallicity;

  /* Number of neighbours contributing to feedback */
  double density_wcount;

  /* Checks that the SN weighting is working (put debugging checks) */
  double delta_m_check;
  double delta_p_norm_check;
  double delta_p_check[3];
  double delta_p_tot[3];
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H */
