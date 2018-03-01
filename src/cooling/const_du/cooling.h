/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Richard Bower (r.g.bower@durham.ac.uk)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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
#ifndef SWIFT_COOLING_CONST_DU_H
#define SWIFT_COOLING_CONST_DU_H

/**
 * @file src/cooling/const_du/cooling.h
 * @brief Routines related to the "constant cooling" cooling function.
 *
 * This is the simplest possible cooling function. A constant cooling rate with
 * a minimal energy floor is applied. Should be used as a template for more
 * realistic functions.
 */

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "const.h"
#include "cooling_struct.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Cooling Model", "Constant du/dt");
}
#endif

/**
 * @brief Apply the cooling function to a particle.
 *
 * In this simple example we just apply the const cooling rate
 * and check that we don't go below the given floor.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {

  /* Internal energy floor */
  const float u_floor = cooling->min_energy;

  /* Get current internal energy */
  const float u_old = hydro_get_physical_internal_energy(p, cosmo);

  /* Current du_dt */
  const float hydro_du_dt = hydro_get_internal_energy_dt(p);

  /* Get cooling function properties */
  float cooling_du_dt = -cooling->cooling_rate;

  /* Integrate cooling equation to enforce energy floor */
  if (u_old + hydro_du_dt * dt < u_floor) {
    cooling_du_dt = 0.f;
  } else if (u_old + (hydro_du_dt + cooling_du_dt) * dt < u_floor) {
    cooling_du_dt = (u_old + dt * hydro_du_dt - u_floor) / dt;
  }

  /* Update the internal energy time derivative */
  hydro_set_internal_energy_dt(p, hydro_du_dt + cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy += -hydro_get_mass(p) * cooling_du_dt * dt;
}

/**
 * @brief Computes the cooling time-step.
 *
 * In this simple example, we return \f$ \alpha \frac{u}{du/dt} \f$.
 * This is used to compute the time-step of the particle. Cooling functions
 * that are solved implicitly can simply return FLT_MAX here.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us, const struct part* restrict p) {

  const float cooling_rate = cooling->cooling_rate;
  const float internal_energy = hydro_get_physical_internal_energy(p, cosmo);
  return cooling->cooling_tstep_mult * internal_energy / fabsf(cooling_rate);
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * In this case, we set the total radiated energy to 0. Note that the particle
 * structure is just passed in for cases where information needs to be read
 * from there.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct cooling_function_data* cooling) {

  xp->cooling_data.radiated_energy = 0.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * In this simple example we jsut return the quantity accumulated in the
 * #cooling_xpart_data of this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Initialises the cooling function properties from the parameter file
 *
 * In this example, we just read in the values from the YAML file without
 * doing any conversions or multiplying any constants in.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cooling The cooling properties to initialize
 */
static INLINE void cooling_init_backend(
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {

  cooling->cooling_rate =
      parser_get_param_double(parameter_file, "ConstCooling:cooling_rate");
  cooling->min_energy =
      parser_get_param_double(parameter_file, "ConstCooling:min_energy");
  cooling->cooling_tstep_mult = parser_get_param_double(
      parameter_file, "ConstCooling:cooling_tstep_mult");
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Constant cooling' with rate %f and floor %f.",
          cooling->cooling_rate, cooling->min_energy);
}

#endif /* SWIFT_COOLING_CONST_DU_H */
