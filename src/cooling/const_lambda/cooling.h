/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_CONST_LAMBDA_H
#define SWIFT_COOLING_CONST_LAMBDA_H

/**
 * @file src/cooling/const_lambda/cooling.h
 * @brief Routines related to the "constant lambda" cooling function.
 *
 * This model assumes a constant cooling rate Lambda irrespective of redshift
 * or density.
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "cooling_properties.h"
#include "cosmology.h"
#include "entropy_floor.h"
#include "error.h"
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
 * @param pressure_floor The properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
INLINE static void cooling_update(
    const struct cosmology* cosmo,
    const struct pressure_floor_props* pressure_floor,
    struct cooling_function_data* cooling, struct space* s) {
  // Add content if required.
}

/**
 * @brief Calculates du/dt in CGS units for a particle.
 *
 * The cooling rate is \f$\frac{du}{dt} = -\frac{\Lambda}{n_H^2}
 * \frac{n_H^2}{\rho} \f$, where \f$ \frac{\Lambda}{n_H^2} \f$ is a constant in
 * this model (lambda_nH2_cgs in #cooling_function_data).
 * The returned value is in physical [erg * g^-1 * s^-1].
 *
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @return The change in energy per unit mass due to cooling for this particle
 * in cgs units [erg * g^-1 * s^-1].
 */
__attribute__((always_inline)) INLINE static double cooling_rate_cgs(
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct cooling_function_data* cooling, const struct part* p) {

  /* Get particle density [g * cm^-3] */
  const double rho = hydro_get_physical_density(p, cosmo);
  const double rho_cgs = rho * cooling->conv_factor_density_to_cgs;

  /* Get Hydrogen mass fraction */
  const double X_H = hydro_props->hydrogen_mass_fraction;

  /* Hydrogen number density (X_H * rho / m_p) [cm^-3] */
  const double n_H_cgs = X_H * rho_cgs * cooling->proton_mass_cgs_inv;

  /* Calculate du_dt ((Lambda / n_H^2) * n_H^2 / rho) */
  const double du_dt_cgs =
      -cooling->lambda_nH2_cgs * n_H_cgs * n_H_cgs / rho_cgs;

  return du_dt_cgs;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props Properties of the entropy floor.
 * @param pressure_floor The properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time Time since Big Bang (or start of the simulation) in internal
 * units.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* phys_const, const struct unit_system* us,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct pressure_floor_props* pressure_floor,
    const struct cooling_function_data* cooling, struct part* p,
    struct xpart* xp, const float dt, const float dt_therm, const double time) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Current energy (in internal units) */
  const float u_old_com = hydro_get_comoving_internal_energy(p, xp);

  /* Y' = RHS of the comoving equation for du/dt that will be integrated
     forward in time using dt_therm */
  const float hydro_du_dt_com = hydro_get_comoving_internal_energy_dt(p);

  /* Calculate cooling du_dt (in cgs units) */
  const double cooling_du_dt_cgs =
      cooling_rate_cgs(cosmo, hydro_props, cooling, p);

  /* Convert to internal units */
  const float cooling_du_dt_physical =
      cooling_du_dt_cgs * cooling->conv_factor_energy_rate_from_cgs;

  /* Add cosmological term to get Y_cooling' */
  const float cooling_du_dt = cooling_du_dt_physical * cosmo->a * cosmo->a /
                              cosmo->a_factor_internal_energy;

  /* Y_total' */
  float total_du_dt = hydro_du_dt_com + cooling_du_dt;

  /* We now need to check that we are not going to go below any of the limits */

  /* Limit imposed by the entropy floor (comoving)
   * (Recall entropy is the same in physical and comoving frames) */
  const float A_floor_com = entropy_floor(p, cosmo, floor_props);
  const float rho_com = hydro_get_comoving_density(p);
  const float u_floor_com =
      gas_internal_energy_from_entropy(rho_com, A_floor_com);

  /* Absolute minimum */
  const float u_minimal_com =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Largest of both limits */
  const float u_limit_com = max(u_minimal_com, u_floor_com);

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_t. */
  if (u_old_com + total_du_dt * 1.5 * dt_therm < u_limit_com) {
    total_du_dt = (u_limit_com - u_old_com) / (1.5f * dt_therm);
  }

  /* Second, check whether the energy used in the prediction could get negative.
   * We need to check for the 1/2 dt kick followed by a full time-step drift
   * that could potentially be for a time-step twice as big. We hence check
   * for 2.5 delta_t but this time against 0 energy not the minimum */
  if (u_old_com + total_du_dt * 2.5 * dt_therm < 0.) {
    total_du_dt = -u_old_com / ((2.5f + 0.0001f) * dt_therm);
  }

  if (cooling->rapid_cooling) {
    const float u_new_com = u_old_com + total_du_dt * dt_therm;
    const float u_new_phys = u_new_com * cosmo->a_factor_internal_energy;
    hydro_set_physical_internal_energy(p, xp, cosmo, u_new_phys);
    hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor,
                                               u_new_phys);
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.);
  } else {
    /* Update the internal energy time derivative */
    hydro_set_comoving_internal_energy_dt(p, total_du_dt);
  }

  const float actual_cooling_du_dt = total_du_dt - hydro_du_dt_com;
  const float actual_cooling_du_dt_physical = actual_cooling_du_dt / cosmo->a /
                                              cosmo->a *
                                              cosmo->a_factor_internal_energy;
  /* Store the radiated energy (assuming dt will not change) */
  xp->cooling_data.radiated_energy +=
      -hydro_get_mass(p) * actual_cooling_du_dt_physical * dt;
}

/**
 * @brief Computes the time-step due to cooling for this particle.
 *
 * We compute a time-step \f$ \alpha \frac{u}{du/dt} \f$ in physical
 * coordinates. \f$\alpha\f$ is a parameter of the cooling function.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended data of the particle.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props, const struct part* restrict p,
    const struct xpart* restrict xp) {

  /* Start with the case where there is no limit */
  if (cooling->cooling_tstep_mult == FLT_MAX) return FLT_MAX;

  /* Get current internal energy and cooling rate */
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const double cooling_du_dt_cgs =
      cooling_rate_cgs(cosmo, hydro_props, cooling, p);

  /* Convert to internal units */
  const float cooling_du_dt = cooling_du_dt_cgs *
                              cooling->conv_factor_energy_rate_from_cgs /
                              cosmo->a_factor_internal_energy;

  /* If we are close to (or below) the limit, we ignore the condition */
  if (u < 1.01f * hydro_props->minimal_internal_energy || cooling_du_dt == 0.f)
    return FLT_MAX;
  else
    return cooling->cooling_tstep_mult * u / fabsf(cooling_du_dt);
}

/**
 * @brief Compute the electron pressure of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
__attribute__((always_inline)) INLINE static double
cooling_get_electron_pressure(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp) {
  return 0;
}

/**
 * @brief Compute the y-Compton contribution of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
__attribute__((always_inline)) INLINE static double cooling_get_ycompton(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp) {
  error("This cooling model does not compute Compton Y!");
  return 0.;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here. Just set the radiated energy counter to 0.
 *
 * @param phys_const The physical constants in internal units.
 * @param cooling The properties of the cooling function.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
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
  const double u = hydro_get_drifted_physical_internal_energy(p, cosmo);

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
 * @brief Returns the subgrid temperature of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_temperature(const struct part* p,
                                                    const struct xpart* xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

/**
 * @brief Returns the subgrid density of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_density(const struct part* p,
                                                const struct xpart* xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
static INLINE void cooling_split_part(struct part* p, struct xpart* xp,
                                      double n) {

  xp->cooling_data.radiated_energy /= n;
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The cooling properties to initialize
 */
static INLINE void cooling_init_backend(struct swift_params* parameter_file,
                                        const struct unit_system* us,
                                        const struct phys_const* phys_const,
                                        const struct hydro_props* hydro_props,
                                        struct cooling_function_data* cooling) {

  /* Read in the cooling parameters */
  cooling->lambda_nH2_cgs =
      parser_get_param_double(parameter_file, "LambdaCooling:lambda_nH2_cgs");
  cooling->cooling_tstep_mult = parser_get_opt_param_float(
      parameter_file, "LambdaCooling:cooling_tstep_mult", FLT_MAX);
  cooling->rapid_cooling = parser_get_opt_param_int(
      parameter_file, "LambdaCooling:rapid_cooling", 0);

  /* Some useful conversion values */
  cooling->conv_factor_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->conv_factor_energy_rate_from_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Useful constants */
  cooling->proton_mass_cgs_inv =
      1. / (phys_const->const_proton_mass *
            units_cgs_conversion_factor(us, UNIT_CONV_MASS));
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

  message(
      "Cooling function is 'Constant lambda' with Lambda/n_H^2=%g [erg * s^-1 "
      "* "
      "cm^3]",
      cooling->lambda_nH2_cgs);

  message("Lambda/n_H^2=%g [internal units]",
          cooling->lambda_nH2_cgs * cooling->conv_factor_energy_rate_from_cgs);

  if (cooling->cooling_tstep_mult == FLT_MAX)
    message("Cooling function time-step size is unlimited");
  else
    message("Cooling function time-step size limited to %f of u/(du/dt)",
            cooling->cooling_tstep_mult);

  if (cooling->rapid_cooling) {
    message("Using rapid cooling");
  } else {
    message("Using normal cooling");
  }
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

#endif /* SWIFT_COOLING_CONST_LAMBDA_H */
