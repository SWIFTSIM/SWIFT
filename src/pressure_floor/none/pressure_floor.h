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

/* Pre-declarations */
struct cosmology;
struct hydro_props;
struct phys_const;
struct part;
struct xpart;
struct swift_params;
struct unit_system;

/**
 * @file src/pressure_floor/none/pressure_floor.h
 * @brief Model without pressure floor
 */

/**
 * @brief Properties of the pressure floor in the NONE model.
 */
struct pressure_floor_props {};

/**
 * @brief Compute the comoving pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * Since this is the 'none' model for floor, there is no floor and
 * we just return the comoving pressure that was received.
 *
 * @param p The #part.
 * @param pfloor The pressure floor properties.
 * @param comoving_pressure The comoving pressure without any pressure floor.
 * @param cosmo The #cosmology model.
 *
 * @return The comoving pressure with the floor.
 */
static INLINE float pressure_floor_get_comoving_pressure(
    const struct part* p, const struct pressure_floor_props* pfloor,
    const float comoving_pressure, const struct cosmology* cosmo) {
  return comoving_pressure;
}

/**
 * @brief Initialise the pressure floor by reading the parameters and converting
 * to internal units.
 *
 * Nothing to do here.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_const The physical constants.
 * @param hydro_props The propoerties of the hydro scheme.
 * @param props The pressure floor properties to fill.
 */
static INLINE void pressure_floor_init(struct pressure_floor_props* props,
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
    const struct pressure_floor_props* props) {

  message("Pressure floor is 'none'");
}

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
 * @brief Write a pressure_floor struct to the given FILE as a stream of bytes.
 *
 * @param pressure_floor the struct
 * @param stream the file stream
 */
static INLINE void pressure_floor_struct_dump(
    const struct pressure_floor_props* pressure_floor, FILE* stream) {}

/**
 * @brief Restore a pressure_floor struct from the given FILE as a stream of
 * bytes.
 *
 * @param pressure_floor the struct
 * @param stream the file stream
 */
static INLINE void pressure_floor_struct_restore(
    struct pressure_floor_props* pressure_floor, FILE* stream) {}

#endif /* SWIFT_PRESSURE_FLOOR_NONE_H */
