/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_GEAR_MF_DIFFUSION_H
#define SWIFT_CHEMISTRY_STRUCT_GEAR_MF_DIFFUSION_H

#include "chemistry_properties.h"

/* Import the right chemistry_part_data definition */
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
#include "hyperbolic/chemistry_struct.h"
#else
#include "parabolic/chemistry_struct.h"
#endif

/**
 * @brief Properties of the chemistry function for #spart.
 */
struct chemistry_spart_data {

  /*! Fraction of the particle mass in a given element */
  double metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];
};

/**
 * @brief Chemical abundances traced by the #bpart in the GEAR model.
 */
struct chemistry_bpart_data {};

/**
 * @brief Chemical abundances traced by the #sink in the GEAR model.
 */
struct chemistry_sink_data {

  /*! Total mass of element in a particle. */
  double metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];
};

#endif /* SWIFT_CHEMISTRY_STRUCT_GEAR_MF_DIFFUSION_H */
