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
__attribute__((always_inline)) static INLINE float pressure_floor_get_pressure(
    const struct part *p, const float rho, const float pressure);

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
struct pressure_floor_properties {

  /*! Jeans factor */
  float n_jeans;

  /*! The constants in internal units (4 G N_jeans^(2/3) / PI) */
  float constants;
};

/**
 * @brief Compute the physical pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * The inputs can be either in physical or comoving coordinates (the output is
 * in the same coordinates).
 *
 * @param p The #part.
 * @param rho The physical or comoving density.
 * @param pressure The physical or comoving pressure without any pressure floor.
 *
 * @return The physical or comoving pressure with the floor.
 */
__attribute__((always_inline)) static INLINE float pressure_floor_get_pressure(
    const struct part *p, const float rho, const float pressure) {

  /* Compute pressure floor */
  float floor = p->h * p->h * rho * pressure_floor_props.constants
    - p->pressure_floor_data.sigma2;
  floor *= rho * hydro_one_over_gamma;

  return fmax(pressure, floor);
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
    struct pressure_floor_properties *props,
    const struct phys_const *phys_const,
    const struct unit_system *us,
    const struct hydro_props *hydro_props,
    struct swift_params *params) {

  /* Read the Jeans factor */
  props->n_jeans =
      parser_get_param_float(params, "GEARPressureFloor:Jeans_factor");

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
    const struct pressure_floor_properties *props) {

  message("Pressure floor is 'GEAR' with:");
  message("Jeans factor: %g", props->n_jeans);
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of pressure floor to the file
 * @param h_grp The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void pressure_floor_print_snapshot(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Pressure floor", "GEAR");
}

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void pressure_floor_end_density(
    struct part* restrict p,
    const struct cosmology* cosmo) {

  /* To finish the turbulence estimation we devide by the density */
  p->pressure_floor_data.sigma2 /=
      pow_dimension(p->h) * hydro_get_physical_density(p, cosmo);

  /* Add the cosmological factor */
  p->pressure_floor_data.sigma2 *= cosmo->a * cosmo->a;
}

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
                                      const struct cosmology* cosmo) {

  /* If part has 0 neighbours, the estimation of turbulence is 0 */
  p->pressure_floor_data.sigma2 = 0.f;
}

/**
 * @brief Sets the pressure_floor properties of the (x-)particles to a valid
 * start state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void pressure_floor_init_part(
    struct part* restrict p, struct xpart* restrict xp) {
  p->pressure_floor_data.sigma2 = 0.f;
}


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
__attribute__((always_inline)) INLINE static void pressure_floor_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    struct part* restrict p,
    struct xpart* restrict xp) {

  pressure_floor_init_part(p, xp);
}

#endif
#endif /* SWIFT_PRESSURE_FLOOR_GEAR_H */
