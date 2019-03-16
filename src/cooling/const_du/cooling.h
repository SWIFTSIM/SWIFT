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
 * This is the simplest possible cooling function. A constant cooling rate
 * (du/dt) with a minimal energy floor is applied. Should be used as a template
 * for more realistic functions.
 *
 * This cooling model does NOT include cosmological terms and will hence yield
 * an incorrect answer when running in co-moving coordinates.
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "entropy_floor.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
INLINE static void cooling_update(const struct cosmology* cosmo,
                                  struct cooling_function_data* cooling,
                                  struct space* s) {
  // Add content if required.
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * In this simple example we just apply the const cooling rate
 * and check that we don't go below the given floor.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, const float dt,
    const float dt_therm) {

  /* Internal energy floor */
  const float u_floor = cooling->min_energy;

  /* Get current internal energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Current du_dt */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Get cooling function properties */
  float cooling_du_dt = -cooling->cooling_rate;

  /* Integrate cooling equation to enforce energy floor */
  if (u_old + hydro_du_dt * dt < u_floor) {
    cooling_du_dt = 0.f;
  } else if (u_old + (hydro_du_dt + cooling_du_dt) * dt < u_floor) {
    cooling_du_dt = (u_old + dt * hydro_du_dt - u_floor) / dt;
  }

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, hydro_du_dt + cooling_du_dt);

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
 * @param hydro_props The properties of the hydro scheme.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extedended particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props, const struct part* restrict p,
    const struct xpart* xp) {

  const float cooling_rate = cooling->cooling_rate;
  const float internal_energy =
      hydro_get_physical_internal_energy(p, xp, cosmo);
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
 * @param phys_const The physical constants in internal units.
 * @param cooling The properties of the cooling function.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
}

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
INLINE static float cooling_get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  /* Physical constants */
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  /* Gas properties */
  const double T_transition = hydro_props->hydrogen_ionization_temperature;
  const double mu_neutral = hydro_props->mu_neutral;
  const double mu_ionised = hydro_props->mu_ionised;

  /* Particle temperature */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Temperature over mean molecular weight */
  const double T_over_mu = hydro_gamma_minus_one * u * m_H / k_B;

  /* Are we above or below the HII -> HI transition? */
  if (T_over_mu > (T_transition + 1.) / mu_ionised)
    return T_over_mu * mu_ionised;
  else if (T_over_mu < (T_transition - 1.) / mu_neutral)
    return T_over_mu * mu_neutral;
  else
    return T_transition;
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
static INLINE void cooling_init_backend(struct swift_params* parameter_file,
                                        const struct unit_system* us,
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
 * @brief Restore cooling tables (if applicable) after
 * restart
 *
 * Nothing to do here
 *
 * @param cooling the cooling_function_data structure
 * @param cosmo cosmology structure
 */
static INLINE void cooling_restore_tables(struct cooling_function_data* cooling,
                                          const struct cosmology* cosmo) {}

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

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
static INLINE void cooling_clean(struct cooling_function_data* cooling) {}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * Nothing to do beyond writing the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
static INLINE void cooling_struct_dump(
    const struct cooling_function_data* cooling, FILE* stream) {
  restart_write_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
                       stream, "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Nothing to do beyond reading the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
static INLINE void cooling_struct_restore(struct cooling_function_data* cooling,
                                          FILE* stream,
                                          const struct cosmology* cosmo) {
  restart_read_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");
}

#endif /* SWIFT_COOLING_CONST_DU_H */
