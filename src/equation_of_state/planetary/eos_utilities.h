/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_EOS_UTILITIES_H
#define SWIFT_PLANETARY_EOS_UTILITIES_H

/**
 * @file equation_of_state/planetary/eos_utilities.h
 *
 * Utitilies for the planetary equations of state.
 */

/* Local headers. */
#include "eos_setup.h"
#include "hm80.h"
#include "ideal_gas.h"
#include "linear.h"
#include "sesame.h"
#include "tillotson.h"

/**
 * @brief The parameters of the equation of state.
 *
 * Parameter structs for each EoS type and unit.
 */
struct eos_parameters {
  struct idg_params idg[eos_count_idg];
  struct Til_params Til[eos_count_Til];
  struct Til_params Til_custom[eos_count_Til_custom];
  struct HM80_params HM80[eos_count_HM80];
  struct SESAME_params SESAME[eos_count_SESAME];
  struct SESAME_params ANEOS[eos_count_ANEOS];
  struct linear_params linear[eos_count_linear];
  struct SESAME_params custom[eos_count_custom];

  struct mat_params mat_params[eos_count_total];
};

/*! Primary EoS parameter struct */
extern struct eos_parameters eos;

#endif /* SWIFT_PLANETARY_EOS_UTILITIES_H */
