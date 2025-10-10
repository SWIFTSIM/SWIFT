/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
struct feedback_part_data {
  /* Save quantities computed in the #hydro density loop for feedback loop */
  struct {
    /*! Neighbour number count. */
    float wcount;
  } density;
};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! Mass received from supernovae */
  float delta_mass;

  /*! Physical specific energy received from supernovae */
  float delta_u;

  /*! Kinetic energy (not specific!) received from supernovae */
  float delta_E_kin;

  /*! Comoving momemtum received from a supernovae */
  float delta_p[3];

  /*! Number of supernovae affecting this particle */
  int number_SN;
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Is the star dead? */
  int is_dead;

  /*! Normalisation factor used for the enrichment. Corresponds to the
     denominator in eq (9) in https://arxiv.org/abs/1707.07010  */
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

  /*! Parameters to be accumulated in the feedback loops. Used to compute the
     vector weights (isotropic distribution) */
  double f_sum_plus_term[3];
  double f_sum_minus_term[3];

  struct {
    /*! Accumulated value for the total energy available in the SN, taking into
       account gas-star motion. This is eq (A4) (lower formula) sum terms in
       https://arxiv.org/abs/2404.16987, without the 0.5*m_ej. */
    double E_total;

    /*! Parameters to determine the coupled energy, momentum and internal energy
       of the SN */
    double beta_1; /* Accumulated value for beta_1 */
    double beta_2; /* Accumulated value for beta_2 */
  } accumulator;

  /*! Sum of the weighted gas properties used to compute terminal momentum */
  float weighted_gas_density;
  double weighted_gas_metallicity;

  /*! Checks that the sum of the fluxes is 0. These ensures the weights are
     properly constructed. */
#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  struct {
    double delta_m;
    double delta_p_norm;
    double delta_p[3];
  } fluxes_conservation_check;
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_MECHANICAL_H */
