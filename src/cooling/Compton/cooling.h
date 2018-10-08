/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_COMPTON_H
#define SWIFT_COOLING_COMPTON_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "const.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Compute the mean molecular weight as a function of temperature for
 * primordial gas.
 *
 * @param T The temperature of the gas [K].
 * @param H_mass_fraction The hydrogen mass fraction of the gas.
 * @param T_trans The temperature of the transition from HII to HI [K].
 */
__attribute__((always_inline, const)) INLINE static double
mean_molecular_weight(const double T, const double H_mass_fraction,
                      const double T_transition) {

  if (T > T_transition)
    return 4. / (8. - 5. * (1. - H_mass_fraction));
  else
    return 4. / (1. + 3. * H_mass_fraction);
}

/**
 * @brief Compute the temperature for a given internal energy per unit mass
 * assuming primordial gas.
 *
 * @param u_cgs The internal energy per unit mass of the gas [erg * g^-1].
 * @param H_mass_fraction The hydrogen mass fraction of the gas.
 * @param T_trans The temperature of the transition from HII to HI [K].
 * @param m_H_cgs The mass of the Hydorgen atom [g].
 * @param k_B_cgs The Boltzmann constant in cgs units [erg * K^-1]
 * @return The temperature of the gas [K]
 */
__attribute__((always_inline, const)) INLINE static double
temperature_from_internal_energy(const double u_cgs,
                                 const double H_mass_fraction,
                                 const double T_transition,
                                 const double m_H_cgs, const double k_B_cgs) {

  const double T_over_mu = hydro_gamma_minus_one * u_cgs * m_H_cgs / k_B_cgs;

  const double mu_high =
      mean_molecular_weight(T_transition + 1., H_mass_fraction, T_transition);
  const double mu_low =
      mean_molecular_weight(T_transition - 1., H_mass_fraction, T_transition);

  if (T_over_mu > (T_transition + 1.) / mu_high)
    return T_over_mu * mu_high;
  else if (T_over_mu < (T_transition - 1.) / mu_low)
    return T_over_mu * mu_low;
  else
    return T_transition;
}

/**
 * @brief Calculates du/dt in CGS units for a particle.
 *
 *
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The #cooling_function_data used in the run.
 * @param z The current redshift.
 * @param u The current internal energy in internal units.
 * @param p Pointer to the particle data.
 * @return The change in energy per unit mass due to cooling for this particle
 * in cgs units [erg * g^-1 * s^-1].
 */
__attribute__((always_inline)) INLINE static double Compton_cooling_rate_cgs(
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct cooling_function_data* cooling, const double z, const double u,
    const struct part* p) {

  /* Get particle density */
  const double rho = hydro_get_physical_density(p, cosmo);
  const double rho_cgs = rho * cooling->conv_factor_density_to_cgs;

  /* Powers of (1 + z) */
  const double zp1 = z + 1.;
  const double zp1p2 = zp1 * zp1;
  const double zp1p4 = zp1p2 * zp1p2; /* (1 + z)^4 */

  /* CMB temperature at this redshift */
  const double T_CMB = cooling->const_T_CMB_0 * zp1;

  /* Gas properties */
  const double H_mass_fraction = hydro_props->hydrogen_mass_fraction;
  const double T_transition = hydro_props->hydrogen_ionization_temperature;

  /* Particle temperature */
  const double u_cgs = u * cooling->conv_factor_energy_to_cgs;
  const double T = temperature_from_internal_energy(
      u_cgs, H_mass_fraction, T_transition, 1., 1.);  // MATTHIEU

  /* Temperature difference with the CMB */
  const double delta_T = T - T_CMB;

  /* Electron density */
  const double electron_density_cgs =
      rho_cgs * cooling->electron_abundance * cooling->proton_mass_cgs_inv;

  /* Compton formula */
  return cooling->const_Compton_rate_cgs * delta_T * zp1p4 *
         electron_density_cgs / rho_cgs;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, const float dt,
    const float dt_therm) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Internal energy floor */
  const float u_floor = hydro_props->minimal_internal_energy;

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Current du_dt in physical coordinates (internal units) */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Calculate cooling du_dt (in cgs units) */
  const double cooling_du_dt_cgs =
      Compton_cooling_rate_cgs(cosmo, hydro_props, cooling, cosmo->z, u_old, p);

  /* Convert to internal units */
  float cooling_du_dt =
      cooling_du_dt_cgs * cooling->conv_factor_energy_rate_from_cgs;

  /* Add cosmological term */
  cooling_du_dt *= cosmo->a * cosmo->a;

  float total_du_dt = hydro_du_dt + cooling_du_dt;

  /* We now need to check that we are not going to go below any of the limits */

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_t. */
  if (u_old + total_du_dt * 1.5 * dt_therm < u_floor) {
    total_du_dt = (u_floor - u_old) / (1.5f * dt_therm);
  }

  /* Second, check whether the energy used in the prediction could get negative.
   * We need to check for the 1/2 dt kick followed by a full time-step drift
   * that could potentially be for a time-step twice as big. We hence check
   * for 2.5 delta_t but this time against 0 energy not the minimum */
  if (u_old + total_du_dt * 2.5 * dt_therm < 0.) {
    total_du_dt = -u_old / ((2.5f + 0.0001f) * dt_therm);
  }

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, total_du_dt);

  /* Store the radiated energy (assuming dt will not change) */
  xp->cooling_data.radiated_energy +=
      -hydro_get_mass(p) * (total_du_dt - hydro_du_dt) * dt_therm;
}

/**
 * @brief Computes the time-step due to cooling for this particle.
 *
 * We impose no time-step limit.
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

  return FLT_MAX;
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
 * @param cosmo The current cosmological model.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
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
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Initialises the cooling properties.
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

  /* Some useful conversion values */
  cooling->conv_factor_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->conv_factor_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  cooling->conv_factor_energy_rate_from_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Useful constants */
  cooling->proton_mass_cgs_inv =
      1. / (phys_const->const_proton_mass *
            units_cgs_conversion_factor(us, UNIT_CONV_MASS));

  /* Temperature of the CMB in CGS */
  const double T_CMB_0 = phys_const->const_T_CMB_0 *
                         units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  cooling->const_T_CMB_0 = T_CMB_0; /* [K] */

  /* Compute the coefficient at the front of the Compton cooling expression */
  const double radiation_constant =
      4. * phys_const->const_stefan_boltzmann / phys_const->const_speed_light_c;
  const double compton_coefficient =
      4. * radiation_constant * phys_const->const_thomson_cross_section *
      phys_const->const_boltzmann_k /
      (phys_const->const_electron_mass * phys_const->const_speed_light_c);
  const float dimension_coefficient[5] = {1, 2, -3, 0, -5};

  /* This should be ~1.0178085e-37 [g cm^2 s^-3 K^-5] */
  const double compton_coefficient_cgs =
      compton_coefficient *
      units_general_cgs_conversion_factor(us, dimension_coefficient);

  /* And now the Compton rate [g cm^2 s^-3 K^-1] == [erg s^-1 K^-1]*/
  cooling->const_Compton_rate_cgs =
      compton_coefficient_cgs * T_CMB_0 * T_CMB_0 * T_CMB_0 * T_CMB_0;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Compton cooling'.");
}

#endif /* SWIFT_COOLING_COMPTON_H */
