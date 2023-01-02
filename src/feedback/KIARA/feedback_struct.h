/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_FEEDBACK_STRUCT_KIARA_H
#define SWIFT_FEEDBACK_STRUCT_KIARA_H

#include "chemistry_struct.h"

/**
 * @brief Feedback fields carried by each hydro particles
 */
struct feedback_part_data {
  /*! remaining time left for decoupling */
  float decoupling_delay_time;

  /*! Number of times decoupled */
  int number_of_times_decoupled;

  /*! The time to shut off cooling for this particle */
  float cooling_shutoff_delay_time;
};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Inverse of normalisation factor used for the enrichment */
  float enrichment_weight_inv;

  /*! Total mass (unweighted) of neighbouring gas particles */
  float ngb_mass;

  /*! Integer number of neighbouring gas particles */
  int num_ngbs;

  /*! SPH-weighted density of the neighbouring gas particles (internal
    * comoving units) */
  float ngb_rho;

  /*! SPH-weighted metallicity of the neighbouring gas particles
    * (dimensionless) */
  float ngb_Z;

  /*! Total (unweighted) number gas neighbours in the stellar kernel */
  int ngb_N;

  /*! Normalisation factor used for the enrichment */
  float enrichment_weight;

  /*! Mass released */
  float mass;

  /*! Total metal mass released */
  float total_metal_mass;

  /*! Total mass released by each element */
  float metal_mass[chemistry_element_count];

  /*! Energy change due to thermal and kinetic energy of ejecta */
  float energy;

  /*! Number of dark matter neighbours in the (gas) neighbourhood */
  int dm_ngb_N;

  /*! DM velocity dispersion in each direction */
  float dm_vel_diff2[3];

  /*! DM 1D vel. disp. from Vogelsberger et al (2013) equation 14. */
  float dm_vel_disp_1d;
};

#endif /* SWIFT_FEEDBACK_STRUCT_KIARA_H */
