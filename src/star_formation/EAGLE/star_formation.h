/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_EAGLE_STAR_FORMATION_H
#define SWIFT_EAGLE_STAR_FORMATION_H

/* Local includes */
#include "adiabatic_index.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "random.h"
#include "stars.h"
#include "units.h"

/**
 * @file src/star_formation/EAGLE/star_formation.h
 * @brief Star formation model used in the EAGLE model
 */

/**
 * @brief Properties of the EAGLE star formation model.
 */
struct star_formation {

  /*! Normalization of the KS star formation law (internal units) */
  double KS_normalization;

  /*! Normalization of the KS star formation law (Msun / kpc^2 / yr) */
  double KS_normalization_MSUNpYRpKPC2;

  /*! Slope of the KS law */
  double KS_power_law;

  /*! Slope of the high density KS law */
  double KS_high_den_power_law;

  /*! KS law High density threshold (internal units) */
  double KS_high_den_thresh;

  /*! KS high density normalization (internal units) */
  double KS_high_den_normalization;

  /*! KS high density normalization (H atoms per cm^3)  */
  double KS_high_den_thresh_HpCM3;

  /*! Critical overdensity */
  double min_over_den;

  /*! Dalla Vecchia & Schaye temperature criteria */
  double temperature_margin_threshold_dex;

  /*! 10^Tdex of Dalla Vecchia & SChaye temperature criteria */
  double ten_to_temperature_margin_threshold_dex;

  /*! gas fraction */
  double fgas;

  /*! Star formation law slope */
  double SF_power_law;

  /*! star formation normalization (internal units) */
  double SF_normalization;

  /*! star formation high density slope */
  double SF_high_den_power_law;

  /*! Star formation high density normalization (internal units) */
  double SF_high_den_normalization;

  /*! Density threshold to form stars (internal units) */
  double density_threshold;

  /*! Density threshold to form stars in user units */
  double density_threshold_HpCM3;

  /*! Maximum density threshold to form stars (internal units) */
  double density_threshold_max;

  /*! Maximum density threshold to form stars (H atoms per cm^3) */
  double density_threshold_max_HpCM3;

  /*! Reference metallicity for metal-dependant threshold */
  double Z0;

  /*! Inverse of reference metallicity */
  double Z0_inv;

  /*! critical density Metallicity power law (internal units) */
  double n_Z0;

  /*! Polytropic index */
  double EOS_polytropic_index;

  /*! EOS density norm (H atoms per cm^3) */
  double EOS_density_norm_HpCM3;

  /*! EOS Temperature norm (Kelvin)  */
  double EOS_temperature_norm_K;

  /*! EOS pressure norm, eq. 13 of Schaye & Dalla Vecchia 2008 (internal units)
   */
  double EOS_pressure_c;

  /*! EOS Temperarure norm, eq. 13 of Schaye & Dalla Vecchia 2008 (internal
   * units) */
  double EOS_temperature_c;

  /*! EOS density norm, eq. 13 of Schaye & Dalla Vecchia 2008 (internal units)
   */
  double EOS_density_c;

  /*! Inverse of EOS density norm (internal units) */
  double EOS_density_c_inv;

  /*! Max physical density (H atoms per cm^3)*/
  double max_gas_density_HpCM3;

  /*! Max physical density (internal units) */
  double max_gas_density;
};

/**
 * @brief Computes the density threshold for star-formation fo a given total
 * metallicity.
 *
 * Follows Schaye (2004) eq. 19 and 24 (see also Schaye et al. 2015, eq. 2).
 *
 * @param Z The metallicity (metal mass fraction).
 * @param starform The properties of the star formation model.
 * @param phys_const The physical constants.
 * @return The physical density threshold for star formation in internal units.
 */
INLINE static double star_formation_threshold(
    const double Z, const struct star_formation* starform,
    const struct phys_const* phys_const) {

  double density_threshold;

  /* Schaye (2004), eq. 19 and 24 */
  if (Z > 0.) {
    density_threshold = starform->density_threshold *
                        powf(Z * starform->Z0_inv, starform->n_Z0);
    density_threshold = min(density_threshold, starform->density_threshold_max);
  } else {
    density_threshold = starform->density_threshold_max;
  }

  /* Convert to mass density */
  return density_threshold * phys_const->const_proton_mass;
}

/**
 * @brief Compute the pressure on the polytropic equation of state for a given
 * Hydrogen number density.
 *
 * Schaye & Dalla Vecchia 2008, eq. 13.
 *
 * @param n_H The Hydrogen number density in internal units.
 * @param starform The properties of the star formation model.
 * @return The pressure on the equation of state in internal units.
 */
INLINE static double EOS_pressure(const double n_H,
                                  const struct star_formation* starform) {

  return starform->EOS_pressure_c *
         pow(n_H * starform->EOS_density_c_inv, starform->EOS_polytropic_index);
}

/**
 * @brief Compute the temperarue on the polytropic equation of state for a given
 * Hydrogen number density.
 *
 * Schaye & Dalla Vecchia 2008, eq. 13 rewritten for temperature
 *
 * @param n_H The Hydrogen number density in internal units.
 * @param starform The properties of the star formation model.
 * @return The temperature on the equation of state in internal units.
 */
INLINE static double EOS_temperature(const double n_H,
                                     const struct star_formation* starform) {

  return starform->EOS_temperature_c *
         pow(n_H, starform->EOS_polytropic_index - 1.);
}

/**
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 *
 */
INLINE static int star_formation_is_star_forming(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* Minimal density (converted from critical density) for star formation */
  const double rho_crit_times_min_over_den =
      cosmo->critical_density * starform->min_over_den;

  /* Physical density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Deside whether we should form stars or not,
   * first we deterime if we have the correct over density
   * if that is true we calculate if either the maximum density
   * threshold is reached or if the metallicity dependent
   * threshold is reached, after this we calculate if the
   * temperature is appropriate */
  if (physical_density < rho_crit_times_min_over_den) return 0;

  /* In this case there are actually multiple possibilities
   * because we also need to check if the physical density exceeded
   * the appropriate limit */

  const double Z = p->chemistry_data.smoothed_metal_mass_fraction_total;
  const double X_H = p->chemistry_data.smoothed_metal_mass_fraction[0];
  const double n_H = physical_density * X_H;

  /* Get the density threshold */
  const double density_threshold =
      star_formation_threshold(Z, starform, phys_const);

  /* Check if it exceeded the minimum density */
  if (n_H < density_threshold) return 0;

  /* Calculate the temperature */
  const double temperature = cooling_get_temperature(phys_const, hydro_props,
                                                     us, cosmo, cooling, p, xp);

  /* Temperature on the equation of state */
  const double temperature_eos = EOS_temperature(n_H, starform);

  /* Check the Scahye & Dalla Vecchia 2012 EOS-based temperature critrion */
  return (temperature <
          temperature_eos * starform->ten_to_temperature_margin_threshold_dex);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const double dt_star) {

  /* Abort early if time-step size is 0 */
  if (dt_star == 0.) {

    xp->sf_data.SFR = 0.f;
    return;
  }

  /* Hydrogen number density of this particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);
  const double X_H = p->chemistry_data.smoothed_metal_mass_fraction[0];
  const double n_H = physical_density * X_H / phys_const->const_proton_mass;

  /* Are we above the threshold for automatic star formation? */
  if (physical_density >
      starform->max_gas_density * phys_const->const_proton_mass) {

    xp->sf_data.SFR = hydro_get_mass(p) / dt_star;
    return;
  }

  /* Pressure on the effective EOS for this particle */
  const double pressure = EOS_pressure(n_H, starform);

  /* Calculate the specific star formation rate */
  double SFRpergasmass;
  if (hydro_get_physical_density(p, cosmo) <
      starform->KS_high_den_thresh * phys_const->const_proton_mass) {

    SFRpergasmass =
        starform->SF_normalization * pow(pressure, starform->SF_power_law);

  } else {

    SFRpergasmass = starform->SF_high_den_normalization *
                    pow(pressure, starform->SF_high_den_power_law);
  }

  /* Store the SFR */
  xp->sf_data.SFR = SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.
 *
 * Equation 21 of Schaye & Dalla Vecchia 2008.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct engine* e,
    const double dt_star) {

  /* Calculate the propability of forming a star */
  const double prob = xp->sf_data.SFR * dt_star / hydro_get_mass(p);

  /* Get a unique random number between 0 and 1 for star formation */
  const double random_number =
      random_unit_interval(p->id, e->ti_current, random_number_star_formation);

  /* Have we been lucky and need to form a star? */
  return (prob > random_number);
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct star_formation* starform, const int with_cosmology) {

  /* Check if it is the first time steps after star formation */
  if (xp->sf_data.SFR > 0.f) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      xp->sf_data.SFR = -e->cosmology->a;
    } else {
      xp->sf_data.SFR = -e->time;
    }
  }
}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 * @param starform the star formation law properties to use.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
INLINE static void star_formation_copy_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology) {

  /* Store the current mass */
  sp->mass = hydro_get_mass(p);

  /* Store the current mass as the initial mass */
  sp->mass_init = hydro_get_mass(p);

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Store the chemistry struct in the star particle */
  sp->chemistry_data = p->chemistry_data;

  /* Store the tracers data */
  sp->tracers_data = xp->tracers_data;

  /* Store the birth density in the star particle */
  sp->birth_density = hydro_get_physical_density(p, cosmo);
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units.
 * @param hydro_props The propertis of the hydro model.
 * @param starform the star formation law properties to initialize
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    struct star_formation* starform) {

  /* Get the Gravitational constant */
  const double G_newton = phys_const->const_newton_G;

  /* Initial Hydrogen abundance (mass fraction) */
  const double X_H = hydro_props->hydrogen_mass_fraction;

  /* Mean molecular weight assuming neutral gas */
  const double mean_molecular_weight = hydro_props->mu_neutral;

  /* Get the surface density unit Msun / pc^2 in internal units */
  const double Msun_per_pc2 =
      phys_const->const_solar_mass /
      (phys_const->const_parsec * phys_const->const_parsec);

  /* Get the SF surface density unit Msun / kpc^2 / yr in internal units */
  const double kpc = 1000. * phys_const->const_parsec;
  const double Msun_per_kpc2_per_year =
      phys_const->const_solar_mass / (kpc * kpc) / phys_const->const_year;

  /* Conversion of number density from cgs */
  const double number_density_from_cgs =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Quantities that have to do with the Normal Kennicutt-
   * Schmidt law will be read in this part of the code*/

  /* Load the equation of state for this model */
  starform->EOS_polytropic_index = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:EOS_gamma_effective");
  starform->EOS_temperature_norm_K = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:EOS_temperature_norm_K");
  starform->EOS_density_norm_HpCM3 = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:EOS_density_norm_H_p_cm3");
  starform->EOS_density_c =
      starform->EOS_density_norm_HpCM3 * number_density_from_cgs;
  starform->EOS_density_c_inv = 1. / starform->EOS_density_c;

  /* Calculate the EOS pressure normalization */
  starform->EOS_pressure_c =
      starform->EOS_density_c * starform->EOS_temperature_norm_K *
      phys_const->const_boltzmann_k / mean_molecular_weight / X_H;

  /* Normalisation of the temperature in the EOS calculatio */
  starform->EOS_temperature_c =
      starform->EOS_pressure_c / phys_const->const_boltzmann_k;
  starform->EOS_temperature_c *=
      pow(starform->EOS_density_c, starform->EOS_polytropic_index);

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:KS_min_over_density");

  /* Read the gas fraction from the file */
  starform->fgas = parser_get_opt_param_double(
      parameter_file, "EAGLEStarFormation:gas_fraction", 1.);

  /* Read the Kennicutt-Schmidt power law exponent */
  starform->KS_power_law =
      parser_get_param_double(parameter_file, "EAGLEStarFormation:KS_exponent");

  /* Calculate the power law of the corresponding star formation Schmidt law */
  starform->SF_power_law = (starform->KS_power_law - 1.) / 2.;

  /* Read the normalization of the KS law in KS law units */
  starform->KS_normalization_MSUNpYRpKPC2 = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:KS_normalisation");

  /* Convert to internal units */
  starform->KS_normalization =
      starform->KS_normalization_MSUNpYRpKPC2 * Msun_per_kpc2_per_year;

  /* Calculate the starformation pre-factor (eq. 12 of Schaye & Dalla Vecchia
   * 2008) */
  starform->SF_normalization =
      starform->KS_normalization * pow(Msun_per_pc2, -starform->KS_power_law) *
      pow(hydro_gamma * starform->fgas / G_newton, starform->SF_power_law);

  /* Read the high density Kennicutt-Schmidt power law exponent */
  starform->KS_high_den_power_law = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:KS_high_density_exponent");

  /* Calculate the SF high density power law */
  starform->SF_high_den_power_law = (starform->KS_high_den_power_law - 1.) / 2.;

  /* Read the high density criteria for the KS law in number density per cm^3 */
  starform->KS_high_den_thresh_HpCM3 = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:KS_high_density_threshold_H_p_cm3");

  /* Transform the KS high density criteria to simulation units */
  starform->KS_high_den_thresh =
      starform->KS_high_den_thresh_HpCM3 * number_density_from_cgs;

  /* Pressure at the high-density threshold */
  const double EOS_high_den_pressure =
      EOS_pressure(starform->KS_high_den_thresh, starform);

  /* Calculate the KS high density normalization
   * We want the SF law to be continous so the normalisation of the second
   * power-law is the value of the first power-law at the high-density threshold
   */
  starform->KS_high_den_normalization =
      starform->KS_normalization *
      pow(Msun_per_pc2,
          starform->KS_high_den_power_law - starform->KS_power_law) *
      pow(hydro_gamma * EOS_high_den_pressure * starform->fgas / G_newton,
          (starform->KS_power_law - starform->KS_high_den_power_law) * 0.5f);

  /* Calculate the SF high density normalization */
  starform->SF_high_den_normalization =
      starform->KS_high_den_normalization *
      pow(Msun_per_pc2, -starform->KS_high_den_power_law) *
      pow(hydro_gamma * starform->fgas / G_newton,
          starform->SF_high_den_power_law);

  /* Get the maximum physical density for SF */
  starform->max_gas_density_HpCM3 = parser_get_opt_param_double(
      parameter_file, "EAGLEStarFormation:KS_max_density_threshold_H_p_cm3",
      FLT_MAX);

  /* Convert the maximum physical density to internal units */
  starform->max_gas_density =
      starform->max_gas_density_HpCM3 * number_density_from_cgs;

  starform->temperature_margin_threshold_dex = parser_get_opt_param_double(
      parameter_file, "EAGLEStarFormation:KS_temperature_margin_dex", FLT_MAX);

  starform->ten_to_temperature_margin_threshold_dex =
      exp10(starform->temperature_margin_threshold_dex);

  /* Read the normalization of the metallicity dependent critical
   * density*/
  starform->density_threshold_HpCM3 = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:threshold_norm_H_p_cm3");

  /* Convert to internal units */
  starform->density_threshold =
      starform->density_threshold_HpCM3 * number_density_from_cgs;

  /* Read the scale metallicity Z0 */
  starform->Z0 = parser_get_param_double(parameter_file,
                                         "EAGLEStarFormation:threshold_Z0");
  starform->Z0_inv = 1. / starform->Z0;

  /* Read the power law of the critical density scaling */
  starform->n_Z0 = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:threshold_slope");

  /* Read the maximum allowed density for star formation */
  starform->density_threshold_max_HpCM3 = parser_get_param_double(
      parameter_file, "EAGLEStarFormation:threshold_max_density_H_p_cm3");

  /* Convert to internal units */
  starform->density_threshold_max =
      starform->density_threshold_max_HpCM3 * number_density_from_cgs;
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  message("Star formation law is EAGLE (Schaye & Dalla Vecchia 2008)");
  message(
      "With properties: normalization = %e Msun/kpc^2/yr, slope of the"
      "Kennicutt-Schmidt law = %e and gas fraction = %e ",
      starform->KS_normalization_MSUNpYRpKPC2, starform->KS_power_law,
      starform->fgas);
  message("At densities of %e H/cm^3 the slope changes to %e.",
          starform->KS_high_den_thresh_HpCM3, starform->KS_high_den_power_law);
  message(
      "The effective equation of state is given by: polytropic "
      "index = %e , normalization density = %e #/cm^3 and normalization "
      "temperature = %e K",
      starform->EOS_polytropic_index, starform->EOS_density_norm_HpCM3,
      starform->EOS_temperature_norm_K);
  message("Density threshold follows Schaye (2004)");
  message(
      "the normalization of the density threshold is given by"
      " %e #/cm^3, with metallicity slope of %e, and metallicity normalization"
      " of %e, the maximum density threshold is given by %e #/cm^3",
      starform->density_threshold_HpCM3, starform->n_Z0, starform->Z0,
      starform->density_threshold_max_HpCM3);
  message("Temperature threshold is given by Dalla Vecchia and Schaye (2012)");
  message("The temperature threshold offset from the EOS is given by: %e dex",
          starform->temperature_margin_threshold_dex);
  message("Running with a maximum gas density given by: %e #/cm^3",
          starform->max_gas_density_HpCM3);
}

#endif /* SWIFT_EAGLE_STAR_FORMATION_H */
