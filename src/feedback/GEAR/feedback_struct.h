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

struct feedback_part_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! specific energy received from supernovae */
  float delta_u;

  /*! Momemtum received from a supernovae */
  float delta_p[3];
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Inverse of normalisation factor used for the enrichment. */
  float enrichment_weight;

  /* Assign two different names for clarification */
  union {
    /*! Number of supernovae */
    float number_sn;

    /*! Energy injected in the surrounding particles */
    float energy_ejected;
  };

  /*! Total mass ejected by the supernovae */
  float mass_ejected;

  /*! Chemical composition of the mass ejected */
  float metal_mass_ejected[GEAR_CHEMISTRY_ELEMENT_COUNT];
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
