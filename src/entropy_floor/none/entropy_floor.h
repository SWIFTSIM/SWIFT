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
#ifndef SWIFT_ENTROPY_FLOOR_NONE_H
#define SWIFT_ENTROPY_FLOOR_NONE_H

/**
 * @file src/entropy_floor/none/entropy_floor.h
 * @brief Empty functions used for simulations without entropy
 * floors.
 */

struct cosmology;
struct hydro_props;
struct part;

/**
 * @brief Properties of the entropy floor.
 *
 * Nothing here.
 */
struct entropy_floor_properties {};

/**
 * @brief Compute the entropy floor of a given #part.
 *
 * Simply return 0 (no floor).
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
static INLINE float entropy_floor(
    const struct part *p, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props) {

  return 0.f;
}

/**
 * @brief Compute the temperature from the entropy floor for a given #part
 *
 * Simply return 0 (no floor).
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
static INLINE float entropy_floor_temperature(
    const struct part *p, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props) {
  return 0.f;
}

/**
 * @brief Initialise the entropy floor by reading the parameters and converting
 * to internal units.
 *
 * Nothing to do here.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_cont The physical constants.
 * @param props The entropy floor properties to fill.
 */
static INLINE void entropy_floor_init(struct entropy_floor_properties *props,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct hydro_props *hydro_props,
                                      struct swift_params *params) {}

/**
 * @brief Print the properties of the entropy floor to stdout.
 *
 * @param props The entropy floor properties.
 */
static INLINE void entropy_floor_print(
    const struct entropy_floor_properties *props) {

  message("Entropy floor is 'no entropy floor'.");
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of entropy floor to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void entropy_floor_write_flavour(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Entropy floor", "None");
}
#endif

/**
 * @brief Write an entropy floor struct to the given FILE as a stream of bytes.
 *
 * Nothing to do here.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void entropy_floor_struct_dump(
    const struct entropy_floor_properties *props, FILE *stream) {}

/**
 * @brief Restore a entropy floor struct from the given FILE as a stream of
 * bytes.
 *
 * Nothing to do here.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void entropy_floor_struct_restore(
    struct entropy_floor_properties *props, FILE *stream) {}

#endif /* SWIFT_ENTROPY_FLOOR_NONE_H */
