/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_PRESSURE_FLOOR_GEAR_H
#define SWIFT_PRESSURE_FLOOR_GEAR_H

/* Forward declaration */
struct cosmology;
__attribute__((always_inline)) static INLINE float
pressure_floor_get_comoving_pressure(const struct part* p, const float pressure,
                                     const struct cosmology* cosmo);

#include "adiabatic_index.h"
#include "cosmology.h"
#include "dimension.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "parser.h"
#include "part.h"
#include "units.h"

/**
 * @file src/pressure_floor/GEAR/pressure_floor.h
 * @brief Pressure floor used in the GEAR model
 */

/**
 * @brief Properties of the pressure floor in the GEAR model.
 */
struct pressure_floor_props {

  /*! Jeans factor */
  float n_jeans;

  /*! The constants in internal units (4 G N_jeans^(2/3) / PI) */
  float constants;
};

/**
 * @brief Compute the comoving pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param pressure_comoving The comoving pressure without any pressure floor.
 * @param cosmo The #cosmology model.
 *
 * @return The comoving pressure with the floor.
 */
__attribute__((always_inline)) static INLINE float
pressure_floor_get_comoving_pressure(const struct part* p,
                                     const float pressure_comoving,
                                     const struct cosmology* cosmo) {

  const float a_coef = pow_three_gamma_minus_one(cosmo->a);
  const float rho = hydro_get_comoving_density(p);

  /* Compute the pressure floor */
  float floor = kernel_gamma * kernel_gamma * p->h * p->h * rho *
                pressure_floor_props.constants * cosmo->a_inv;
  floor *= a_coef * rho * hydro_one_over_gamma;

  return fmaxf(pressure_comoving, floor);
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
__attribute__((always_inline)) static INLINE void pressure_floor_init(
    struct pressure_floor_props* props, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    struct swift_params* params) {

  /* Read the Jeans factor */
  props->n_jeans =
      parser_get_param_float(params, "GEARPressureFloor:jeans_factor");

  /* Compute the constants */
  props->constants =
      4.0 * M_1_PI * phys_const->const_newton_G * pow(props->n_jeans, 2. / 3.);
}

/**
 * @brief Print the properties of the pressure floor to stdout.
 *
 * @param props The pressure floor properties.
 */
__attribute__((always_inline)) static INLINE void pressure_floor_print(
    const struct pressure_floor_props* props) {

  message("Pressure floor is 'GEAR' with:");
  message("Jeans factor  : %g", props->n_jeans);
  message("      constant: %f", props->constants);
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of pressure floor to the file
 * @param h_grp The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void pressure_floor_print_snapshot(
    hid_t h_grp) {

  io_write_attribute_s(h_grp, "Pressure floor", "GEAR");
}

#endif

/**
 * @brief Finishes the density calculation for the pressure floor properties.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void pressure_floor_end_density(
    struct part* restrict p, const struct cosmology* cosmo) {}

/**
 * @brief Sets all the pressure floor fields to sensible values when the #part
 * has 0 ngbs.
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
                               struct xpart* restrict xp) {

  pressure_floor_init_part(p, xp);
}

/**
 * @brief Write a pressure_floor struct to the given FILE as a stream of bytes.
 *
 * @param pressure_floor the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void pressure_floor_struct_dump(
    const struct pressure_floor_props* pressure_floor, FILE* stream) {

  restart_write_blocks((void*)pressure_floor,
                       sizeof(struct pressure_floor_props), 1, stream,
                       "pressure_floor", "pressure_floor");

  message("dumping pressure_floor...");
  pressure_floor_print(pressure_floor);
}

/**
 * @brief Restore a pressure_floor struct from the given FILE as a stream of
 * bytes.
 *
 * @param pressure_floor the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void pressure_floor_struct_restore(
    struct pressure_floor_props* pressure_floor, FILE* stream) {

  restart_read_blocks((void*)pressure_floor,
                      sizeof(struct pressure_floor_props), 1, stream, NULL,
                      "pressure_floor");

  message("restoring pressure_floor...");
  pressure_floor_print(pressure_floor);

  /* copy to the global variable */
  pressure_floor_props.n_jeans = pressure_floor->n_jeans;
  pressure_floor_props.constants = pressure_floor->constants;
}

#endif /* SWIFT_PRESSURE_FLOOR_GEAR_H */
