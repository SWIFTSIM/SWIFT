/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_NONE_COOLING_SUBGRID_H
#define SWIFT_NONE_COOLING_SUBGRID_H


/**
 * @file src/cooling/grackle/cooling_none_subgrid.h
 * @brief Subgrid model for none cooling, independent from grackle.
 */

/**
 * @brief Update the gas properties with subgrid physics.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 */
INLINE static void cooling_update_part_subgrid(const struct phys_const *phys_const,
				 const struct unit_system *us,
				 const struct cosmology *cosmo,
				 const struct hydro_props *hydro_props,
				 const struct cooling_function_data *cooling,
				 struct part *p, struct xpart *xp, double dt,
				 double dt_therm) {}

#endif /* SWIFT_NONE_COOLING_SUBGRID_H */
