/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.reavz@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_STRUCT_AGORA_H
#define SWIFT_FEEDBACK_STRUCT_AGORA_H

#include "chemistry_struct.h"

/**
 * @brief The stellar feedback type for each star type. Now, star particles can
 * represent a single star ("single_star"), a stellar population without SNII
 * feedback ("star_population_no_SNII") or a stellar population with SNII
 * feedback ("stellar population").
 */
enum star_feedback_type {
  single_star,             /* particle representing a single star */
  star_population_no_SNII, /* particle representing a population without SNII */
  star_population          /* particle representing a population (with SNII) */
};

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
};

/**
 * @brief Feedback fields carried by each star particles
 *
 * Nothing here since this is a no-feedback model.
 */
struct feedback_spart_data {

  /*! Inverse of normalisation factor used for the enrichment. */
  float enrichment_weight;

  /*! Energy injected in the surrounding particles */
  float energy_ejected;

  /*! Total mass ejected by the supernovae */
  float mass_ejected;

  /*! Metals ejected by the supernovae */
  float metal_mass_ejected[AGORA_CHEMISTRY_ELEMENT_COUNT];

  /*! Does the particle needs the feedback loop? */
  char will_do_feedback;

  /*! Set the stellar particle to idle regarding the feedback */
  char idle;

  /* Feedback type in function of the star particle type */
  enum star_feedback_type star_type;
};

#endif /* SWIFT_FEEDBACK_STRUCT_AGORA_H */
