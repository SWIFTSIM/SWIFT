/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2026 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_DAMAGE_DATA_H
#define SWIFT_DAMAGE_DATA_H

/**
 * @file strength/damage/damage_data.h
 * @brief The damage of a particle, stored on the particle and in the strength
 * state.
 */

/**
 * @brief The damage of a particle.
 */
struct strength_damage {

  /*! Total damage */
  float total;

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
  /*! Damage accumulated due to tension */
  float tensile;
#endif

#ifdef STRENGTH_DAMAGE_SHEAR_COLLINS
  /*! Damage accumulated due to shear */
  float shear;
#endif
};

#endif /* SWIFT_DAMAGE_DATA_H */
