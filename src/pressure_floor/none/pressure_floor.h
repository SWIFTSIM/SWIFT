/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PRESSURE_FLOOR_NONE_H
#define SWIFT_PRESSURE_FLOOR_NONE_H

#include "adiabatic_index.h"
#include "cosmology.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "parser.h"
#include "part.h"
#include "units.h"

/**
 * @file src/pressure_floor/none/pressure_floor.h
 * @brief Model without pressure floor
 */

/**
 * @brief Properties of the pressure floor in the NONE model.
 */
struct pressure_floor_properties {};

/**
 * @brief Compute the physical pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param physical_pressure The physical pressure without any pressure floor.
 * @param cosmo The #cosmology model.
 *
 * @return The physical pressure with the floor.
 */
static INLINE float pressure_floor_get_physical_pressure(
    const struct part* p, const float physical_pressure, const struct cosmology *cosmo) {
  return physical_pressure;
}

/**
 * @brief Compute the comoving pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param comoving_pressure The comoving pressure without any pressure floor.
 * @param cosmo The #cosmology model.
 *
 * @return The comoving pressure with the floor.
 */
static INLINE float pressure_floor_get_comoving_pressure(
    const struct part* p, const float comoving_pressure, const struct cosmology *cosmo) {
  return comoving_pressure;
}

/**
 * @brief Initialise the pressure floor by reading the parameters and converting
 * to internal units.
 *
 * The input temperatures and number densities are converted to pressure and
 * density assuming a neutral gas of primoridal abundance.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_const The physical constants.
 * @param hydro_props The propoerties of the hydro scheme.
 * @param props The pressure floor properties to fill.
 */
static INLINE void pressure_floor_init(struct pressure_floor_properties* props,
                                       const struct phys_const* phys_const,
                                       const struct unit_system* us,
                                       const struct hydro_props* hydro_props,
                                       struct swift_params* params) {}

/**
 * @brief Print the properties of the pressure floor to stdout.
 *
 * @param props The pressure floor properties.
 */
static INLINE void pressure_floor_print(
    const struct pressure_floor_properties* props) {}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of pressure floor to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void pressure_floor_print_snapshot(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Pressure floor", "none");
}
#endif

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void pressure_floor_end_density(
    struct part* restrict p, const struct cosmology* cosmo) {}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
pressure_floor_part_has_no_neighbours(struct part* restrict p,
                                      struct xpart* restrict xp,
                                      const struct cosmology* cosmo) {}

/**
 * @brief Sets the pressure_floor properties of the (x-)particles to a valid
 * start state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void pressure_floor_init_part(
    struct part* restrict p, struct xpart* restrict xp) {}

/**
 * @brief Sets the pressure_floor properties of the (x-)particles to a valid
 * start state.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
pressure_floor_first_init_part(const struct phys_const* restrict phys_const,
                               const struct unit_system* restrict us,
                               const struct cosmology* restrict cosmo,
                               struct part* restrict p,
                               struct xpart* restrict xp) {}

#endif /* SWIFT_PRESSURE_FLOOR_NONE_H */
