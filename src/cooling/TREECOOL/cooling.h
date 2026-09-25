/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_TREECOOL_H
#define SWIFT_COOLING_TREECOOL_H

/**
 * @file src/cooling/TREECOOL/cooling.h
 * @brief Routines related to the TREECOOL cooling function, i.e. the
 * primordial cooling function of Katz, Weinberg & Hernquist (1996) combined
 * with a UV background tabulated in a TREECOOL file.
 *
 * The gas is assumed to be primordial (Hydrogen and Helium only) and in
 * ionization equilibrium with both the collisional processes and an optically
 * thin, spatially uniform UV background read from a TREECOOL table. The rate
 * coefficients and the list of cooling and heating processes are those of
 * Katz, Weinberg & Hernquist (1996), ApJS, 105, 19, Tables 1 and 2.
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "cooling_properties.h"
#include "cooling_rates.h"
#include "cooling_tables.h"
#include "cosmology.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "restart.h"
#include "units.h"

/*! Maximal number of iterations of the bisection scheme */
#define treecool_bisection_max_iterations 150

/*! Relative tolerance of the bisection scheme */
#define treecool_bisection_tolerance 1.e-6f

/*! Factor by which the bracket is widened at each attempt */
#define treecool_bracket_factor 1.5

/*! Relative change in energy below which the explicit solution is used */
#define treecool_explicit_tolerance 0.05f

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * We use this to interpolate the UV background to the current redshift, which
 * spares us from doing it once per particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param pressure_floor The properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 * @param time The current system time.
 */
INLINE static void cooling_update(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor,
    struct cooling_function_data *cooling, struct space *s, const double time) {

  /* Note that for non-cosmological runs the cosmology module reports z = 0, so
   * such runs use the z = 0 entry of the TREECOOL table. */
  const double redshift = cosmo->z;

  treecool_set_UV_background(cooling, redshift);

  /* Properties of the CMB the gas is Compton-cooling against */
  cooling->T_CMB_cgs =
      phys_const->const_T_CMB_0 * cooling->temperature_to_cgs * (1. + redshift);
  cooling->one_plus_z_to_the_4 =
      (1. + redshift) * (1. + redshift) * (1. + redshift) * (1. + redshift);
}

/**
 * @brief Solve the implicit cooling problem by bisection.
 *
 * We are looking for the energy u such that
 *
 *   f(u) = u - u_ini - ratefact * Lambda(u) * dt = 0,
 *
 * where Lambda is the net cooling rate per n_H^2. Since f(+inf) > 0 and
 * f(0) < 0 whenever the gas cools, the root always exists and bisection is
 * guaranteed to find it.
 *
 * @param u_ini_cgs The internal energy at the start of the step [erg * g^-1].
 * @param n_H_cgs The Hydrogen number density [cm^-3].
 * @param ratefact_cgs The factor n_H^2 / rho [cm^-3 * g^-1 * cm^3].
 * @param dt_cgs The time-step [s].
 * @param cooling The #cooling_function_data used in the run.
 * @param gas The #treecool_gas_state carrying the electron density guess.
 * @param id The ID of the particle (used for error messages).
 *
 * @return The internal energy at the end of the step [erg * g^-1].
 */
__attribute__((always_inline)) INLINE static double treecool_bisection_iter(
    const double u_ini_cgs, const double n_H_cgs, const double ratefact_cgs,
    const double dt_cgs, const struct cooling_function_data *cooling,
    struct treecool_gas_state *gas, const long long id) {

  double u_lower_cgs = max(u_ini_cgs, cooling->u_min_cgs);
  double u_upper_cgs = u_lower_cgs;

  double LambdaNet_cgs =
      treecool_cooling_rate_from_u(cooling, u_ini_cgs, n_H_cgs, gas);

  int i = 0;

  if (LambdaNet_cgs < 0.) {

    /* We are cooling: bracket the solution from below */
    u_lower_cgs =
        max(u_lower_cgs / treecool_bracket_factor, cooling->u_min_cgs);
    u_upper_cgs =
        max(u_upper_cgs * treecool_bracket_factor, cooling->u_min_cgs);

    LambdaNet_cgs =
        treecool_cooling_rate_from_u(cooling, u_lower_cgs, n_H_cgs, gas);

    while (u_lower_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs >
               0. &&
           i < treecool_bisection_max_iterations) {

      /* The implicit solution lies below the energy floor: stop there. (If
       * it lies between the floor and u_upper, the condition above fails once
       * u_lower reaches the floor and we proceed to the bisection.) */
      if (u_lower_cgs <= cooling->u_min_cgs) return cooling->u_min_cgs;

      u_lower_cgs =
          max(u_lower_cgs / treecool_bracket_factor, cooling->u_min_cgs);
      u_upper_cgs =
          max(u_upper_cgs / treecool_bracket_factor, cooling->u_min_cgs);

      LambdaNet_cgs =
          treecool_cooling_rate_from_u(cooling, u_lower_cgs, n_H_cgs, gas);

      ++i;
    }

  } else {

    /* We are heating: bracket the solution from above */
    u_lower_cgs /= treecool_bracket_factor;
    u_upper_cgs *= treecool_bracket_factor;

    LambdaNet_cgs =
        treecool_cooling_rate_from_u(cooling, u_upper_cgs, n_H_cgs, gas);

    while (u_upper_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs <
               0. &&
           i < treecool_bisection_max_iterations) {

      u_lower_cgs *= treecool_bracket_factor;
      u_upper_cgs *= treecool_bracket_factor;

      LambdaNet_cgs =
          treecool_cooling_rate_from_u(cooling, u_upper_cgs, n_H_cgs, gas);

      ++i;
    }
  }

  if (i >= treecool_bisection_max_iterations)
    error(
        "Particle %lld exceeded max iterations searching for the cooling "
        "bracket: n_H=%.4e cm^-3, u_ini=%.4e erg/g, dt=%.4e s",
        id, n_H_cgs, u_ini_cgs, dt_cgs);

  /* We now have a bracket: shrink it down to the requested tolerance */
  double u_next_cgs;
  i = 0;

  do {

    u_next_cgs = 0.5 * (u_lower_cgs + u_upper_cgs);

    LambdaNet_cgs =
        treecool_cooling_rate_from_u(cooling, u_next_cgs, n_H_cgs, gas);

    if (u_next_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs > 0.)
      u_upper_cgs = u_next_cgs;
    else
      u_lower_cgs = u_next_cgs;

    ++i;

  } while (fabs(u_upper_cgs - u_lower_cgs) / u_next_cgs >
               treecool_bisection_tolerance &&
           i < treecool_bisection_max_iterations);

  if (i >= treecool_bisection_max_iterations)
    error("Particle %lld failed to converge in the bisection scheme", id);

  return u_upper_cgs;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * We compute u_new such that u_new = u_old + dt * du/dt(u_new) at fixed
 * density and redshift. If the change in energy over the step is small enough
 * we simply use the explicit solution; otherwise we solve the implicit problem
 * by bisection.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props Properties of the entropy floor.
 * @param pressure_floor The properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle's extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time Time since Big Bang (or start of the simulation) in internal
 * units.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, const float dt, const float dt_therm, const double time) {

  /* Nothing to do over a zero time-step */
  if (dt == 0.) return;

  /* Internal energy at the last kick step */
  const float u_start = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Change in internal energy due to the hydro forces */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Internal energy at the end of the next kick step (assuming dt does not
   * change) */
  double u_0 = u_start + hydro_du_dt * dt_therm;
  u_0 = max(u_0, hydro_props->minimal_internal_energy);

  /* Convert to physical cgs units */
  const double u_0_cgs = u_0 * cooling->internal_energy_to_cgs;
  const double dt_cgs = dt * cooling->time_to_cgs;

  /* Hydrogen number density in physical cgs units */
  const double rho_cgs =
      hydro_get_physical_density(p, cosmo) * cooling->density_to_cgs;
  const double n_H_cgs = rho_cgs * cooling->X_H * cooling->inv_proton_mass_cgs;

  /* n_H^2 / rho, written so as to avoid a round-off prone division */
  const double ratefact_cgs =
      n_H_cgs * cooling->X_H * cooling->inv_proton_mass_cgs;

  /* Start the solvers from the electron density this particle had last time */
  struct treecool_gas_state gas;
  gas.n_e = xp->cooling_data.electron_fraction;

  /* Net cooling rate at the start of the step */
  const double LambdaNet_cgs =
      treecool_cooling_rate_from_u(cooling, u_0_cgs, n_H_cgs, &gas);

  double u_final_cgs;

  if (fabs(ratefact_cgs * LambdaNet_cgs * dt_cgs) <
      treecool_explicit_tolerance * u_0_cgs) {

    /* The change is small: the explicit solution is accurate enough */
    u_final_cgs = u_0_cgs + ratefact_cgs * LambdaNet_cgs * dt_cgs;

  } else {

    u_final_cgs = treecool_bisection_iter(u_0_cgs, n_H_cgs, ratefact_cgs,
                                          dt_cgs, cooling, &gas, p->id);
  }

  /* Remember the electron density for the next step */
  xp->cooling_data.electron_fraction = gas.n_e;

  /* Back to internal units */
  double u_final = u_final_cgs * cooling->internal_energy_from_cgs;

  /* We now need to check that we are not going below any of the limits */

  /* Absolute minimum */
  u_final = max(u_final, (double)hydro_props->minimal_internal_energy);

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho_physical = hydro_get_physical_density(p, cosmo);
  const double u_floor =
      gas_internal_energy_from_entropy(rho_physical, A_floor);
  u_final = max(u_final, u_floor);

  /* Expected change in energy over the next kick step (assuming dt does not
   * change) */
  const double delta_u = u_final - max((double)u_start, u_floor);

  /* Determine whether we are in the slow- or rapid-cooling regime by comparing
   * dt / t_cool to the threshold. Note that dt / t_cool = |delta_u| / u. */
  const double dt_over_t_cool = fabs(delta_u) / max((double)u_start, u_floor);

  if (cooling->rapid_cooling_threshold >= 0. &&
      dt_over_t_cool >= cooling->rapid_cooling_threshold) {

    /* Rapid-cooling regime: apply the change straight to the energy */
    hydro_set_physical_internal_energy(p, xp, cosmo, u_final);
    hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor,
                                               u_final);
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.);

  } else {

    /* Slow-cooling regime: update du/dt so that we can drift the energy */
    hydro_set_physical_internal_energy_dt(p, cosmo, delta_u / dt_therm);
  }

  /* Store the radiated energy. This is the change due to the cooling alone,
   * i.e. excluding the contribution of the hydro forces already folded into
   * u_0. */
  xp->cooling_data.radiated_energy -=
      hydro_get_mass(p) * (u_final - max(u_0, u_floor));
}

/**
 * @brief Computes the time-step due to cooling for this particle.
 *
 * We do not impose any limit on the time-step: the implicit solver is stable
 * for any time-step size.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended data of the particle.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props, const struct part *restrict p,
    const struct xpart *restrict xp) {

  return FLT_MAX;
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
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp) {

  const double u_cgs = hydro_get_drifted_physical_internal_energy(p, cosmo) *
                       cooling->internal_energy_to_cgs;

  const double rho_cgs =
      hydro_get_physical_density(p, cosmo) * cooling->density_to_cgs;
  const double n_H_cgs = rho_cgs * cooling->X_H * cooling->inv_proton_mass_cgs;

  struct treecool_gas_state gas;
  gas.n_e = xp->cooling_data.electron_fraction;

  return treecool_temperature_from_u(cooling, u_cgs, n_H_cgs, &gas);
}

/**
 * @brief Compute the electron number density of a #part.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 *
 * @return The physical electron number density in internal units
 * [U_L^-3].
 */
INLINE static double cooling_get_electron_density(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  /* Physical Hydrogen number density in internal units */
  const double n_H = hydro_get_physical_density(p, cosmo) * cooling->X_H /
                     phys_const->const_proton_mass;

  return xp->cooling_data.electron_fraction * n_H;
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
cooling_get_electron_pressure(const struct phys_const *phys_const,
                              const struct hydro_props *hydro_props,
                              const struct unit_system *us,
                              const struct cosmology *cosmo,
                              const struct cooling_function_data *cooling,
                              const struct part *p, const struct xpart *xp) {
  return 0.;
}

/**
 * @brief Compute the y-Compton contribution of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return an error.
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
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  error("This cooling model does not compute Compton Y!");
  return 0.;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * The electron fraction is set to its equilibrium value later on, in
 * cooling_post_init_part(), once the densities are known.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.electron_fraction = 0.f;
}

/**
 * @brief Perform additional init on the cooling properties of the
 * (x-)particles that requires the density to be known.
 *
 * We compute the equilibrium electron fraction so that the initial snapshot
 * and the first cooling step start from a meaningful value. Note that
 * cooling_update() has already been called at this point, so the UV
 * background is the one of the starting redshift.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_post_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *cooling, const struct part *restrict p,
    struct xpart *restrict xp) {

  const double u_cgs = hydro_get_physical_internal_energy(p, xp, cosmo) *
                       cooling->internal_energy_to_cgs;

  const double rho_cgs =
      hydro_get_physical_density(p, cosmo) * cooling->density_to_cgs;
  const double n_H_cgs = rho_cgs * cooling->X_H * cooling->inv_proton_mass_cgs;

  struct treecool_gas_state gas;
  gas.n_e = 1.;
  treecool_temperature_from_u(cooling, u_cgs, n_H_cgs, &gas);

  xp->cooling_data.electron_fraction = gas.n_e;
}

/**
 * @brief Returns the subgrid temperature of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_temperature(const struct part *p,
                                                    const struct xpart *xp) {
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
INLINE static float cooling_get_subgrid_density(const struct part *p,
                                                const struct xpart *xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Split the cooling content of a particle into n pieces.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
INLINE static void cooling_split_part(struct part *p, struct xpart *xp,
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
 * @param cooling The cooling properties to initialize.
 */
INLINE static void cooling_init_backend(struct swift_params *parameter_file,
                                        const struct unit_system *us,
                                        const struct phys_const *phys_const,
                                        const struct hydro_props *hydro_props,
                                        struct cooling_function_data *cooling) {

  /* Read the parameters */
  parser_get_param_string(parameter_file, "TREECOOLCooling:TREECOOL_file",
                          cooling->TREECOOL_file);

  cooling->rapid_cooling_threshold = parser_get_opt_param_float(
      parameter_file, "TREECOOLCooling:rapid_cooling_threshold", 0.333333f);

  cooling->UV_background_start_redshift = parser_get_opt_param_float(
      parameter_file, "TREECOOLCooling:UV_background_start_redshift", FLT_MAX);

  cooling->with_Compton_cooling = parser_get_opt_param_int(
      parameter_file, "TREECOOLCooling:with_Compton_cooling", 1);

  cooling->log10_T_min_cgs = parser_get_opt_param_double(
      parameter_file, "TREECOOLCooling:log10_T_min", 1.);

  cooling->log10_T_max_cgs = parser_get_opt_param_double(
      parameter_file, "TREECOOLCooling:log10_T_max", 9.);

  if (cooling->log10_T_max_cgs <= cooling->log10_T_min_cgs)
    error("TREECOOLCooling:log10_T_max must be larger than log10_T_min");

  /* Composition of the primordial gas */
  cooling->Y_He = phys_const->const_primordial_He_fraction;
  cooling->X_H = 1. - cooling->Y_He;
  cooling->y_He = 0.25 * cooling->Y_He / cooling->X_H;

  /* Conversion factors to and from cgs */
  cooling->internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  cooling->internal_energy_from_cgs = 1. / cooling->internal_energy_to_cgs;
  cooling->density_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->temperature_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Useful constants in cgs units */
  cooling->proton_mass_cgs = phys_const->const_proton_mass *
                             units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->inv_proton_mass_cgs = 1. / cooling->proton_mass_cgs;

  const float dimension_k[5] = {1, 2, -2, 0, -1};
  cooling->boltzmann_k_cgs =
      phys_const->const_boltzmann_k *
      units_general_cgs_conversion_factor(us, dimension_k);

  cooling->u_min_cgs =
      hydro_props->minimal_internal_energy * cooling->internal_energy_to_cgs;

  /* Build the tables of rate coefficients and read the UV background */
  treecool_make_rate_table(cooling);
  treecool_read_table(cooling);

  /* Provide sensible z = 0 values until cooling_update() is called for the
   * first time at the start of the run. */
  treecool_set_UV_background(cooling, /*redshift=*/0.);
  cooling->T_CMB_cgs = phys_const->const_T_CMB_0 * cooling->temperature_to_cgs;
  cooling->one_plus_z_to_the_4 = 1.;
}

/**
 * @brief Restore the cooling tables after a restart.
 *
 * The tables are stored in the cooling structure itself, so they are restored
 * alongside it and there is nothing to do here.
 *
 * @param cooling The #cooling_function_data.
 * @param cosmo The #cosmology structure.
 */
INLINE static void cooling_restore_tables(struct cooling_function_data *cooling,
                                          const struct cosmology *cosmo) {}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
INLINE static void cooling_print_backend(
    const struct cooling_function_data *cooling) {

  message(
      "Cooling function is 'TREECOOL' (primordial H/He in ionization "
      "equilibrium, following Katz, Weinberg & Hernquist 1996)");

  message(
      "UV background read from '%s' (%d redshift entries, z < %g)",
      cooling->TREECOOL_file, cooling->N_redshifts,
      exp10(cooling->TREECOOL_log10_1_plus_z[cooling->N_redshifts - 1]) - 1.);

  message("Rate coefficients tabulated for %g < log10(T/K) < %g in %d bins",
          cooling->log10_T_min_cgs, cooling->log10_T_max_cgs,
          treecool_cooling_N_temperature);

  message("Gas composition: X_H = %g, Y_He = %g (n_He / n_H = %g)",
          cooling->X_H, cooling->Y_He, cooling->y_He);

  if (cooling->UV_background_start_redshift == FLT_MAX)
    message("UV background switched on over the whole range of the table");
  else
    message("UV background switched on at z = %g",
            cooling->UV_background_start_redshift);

  if (cooling->with_Compton_cooling)
    message("Inverse Compton cooling off the CMB is included");
  else
    message("Inverse Compton cooling off the CMB is *not* included");

  if (cooling->rapid_cooling_threshold < 0.)
    message("Always using the slow-cooling regime");
  else
    message("Switching to rapid cooling for dt / t_cool > %g",
            cooling->rapid_cooling_threshold);
}

/**
 * @brief Clean-up the memory allocated for the cooling routines.
 *
 * Nothing to do here: the tables live inside the cooling structure.
 *
 * @param cooling The cooling data structure.
 */
INLINE static void cooling_clean(struct cooling_function_data *cooling) {}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * Nothing to do beyond writing the structure to the stream.
 *
 * @param cooling The struct.
 * @param stream The file stream.
 */
INLINE static void cooling_struct_dump(
    const struct cooling_function_data *cooling, FILE *stream) {

  restart_write_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                       stream, "cooling", "cooling function");
}

/**
 * @brief Restore a cooling struct from the given FILE as a stream of bytes.
 *
 * Nothing to do beyond reading the structure from the stream.
 *
 * @param cooling The struct.
 * @param stream The file stream.
 * @param cosmo The #cosmology structure.
 */
INLINE static void cooling_struct_restore(struct cooling_function_data *cooling,
                                          FILE *stream,
                                          const struct cosmology *cosmo) {

  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");
}

#endif /* SWIFT_COOLING_TREECOOL_H */
