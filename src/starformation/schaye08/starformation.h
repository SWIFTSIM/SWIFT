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
#ifndef SWIFT_SCHAYE_STARFORMATION_H
#define SWIFT_SCHAYE_STARFORMATION_H

/* Some standard headers */
#include <stdlib.h>

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
#include "stars.h"
#include "units.h"

/* Starformation struct */
struct star_formation {

  /*! Normalization of the KS star formation law (internal units) */
  double KS_normalization;

  /*! Normalization of the KS star formation law in user units */
  double KS_normalization_MSUNpYRpKPC2;

  /*! Slope of the KS law */
  double KS_power_law;

  /*! Slope of the high density KS law */
  double KS_high_den_power_law;

  /*! KS law High density threshold (internal units) */
  double KS_high_den_thresh;

  /*! KS high density normalization (internal units) */
  double KS_high_den_normalization;

  /*! KS high density normalization (HpCM3)  */
  double KS_high_den_thresh_HpCM3;

  /*! Critical overdensity */
  double min_over_den;

  /*! Solar mass per square parsec */
  double Msunpsquaredpc;

  /*! Temperature threshold */
  double Temperature_threshold;

  /*! gas fraction */
  float fgas;

  /*! Star formation law slope */
  double SF_power_law;

  /*! star formation normalization of schaye+08 (internal units) */
  double SF_normalization;

  /*! star formation high density slope */
  double SF_high_den_power_law;

  /*! Star formation high density normalization (internal units) */
  double SF_high_den_normalization;

  /*! Inverse of RAND_MAX */
  double inv_RAND_MAX;

  /*! Density threshold to form stars (internal units) */
  double density_threshold;

  /*! Density threshold to form stars in user units */
  float density_threshold_HpCM3;

  /*! Maximum density threshold to form stars (internal units) */
  double density_threshold_max;

  /*! Maximum density threshold to form stars in user units */
  float density_threshold_max_HpCM3;

  /*! Scaling metallicity */
  double Z0;

  /*! one over the scaling metallicity */
  double Z0_inv;

  /*! critical density Metallicity power law (internal units) */
  double n_Z0;

  /*! Polytropic index */
  double polytropic_index;

  /*! EOS pressure norm */
  double EOS_pressure_norm;

  /*! EOS Temperature norm */
  double EOS_temperature_norm;

  /*! EOS density norm (internal units) */
  double EOS_density_norm;

  /*! EOS density norm in user units */
  float EOS_density_norm_HpCM3;
};

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
INLINE static int star_formation_potential_to_become_star(
    const struct star_formation* starform, const struct part* restrict p,
    const struct xpart* restrict xp, const struct phys_const* const phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* Read the critical overdensity factor and the critical density of
   * the universe to determine the critical density to form stars*/
  const double rho_crit_times_min_over_den =
      cosmo->critical_density * starform->min_over_den;
  const double particle_density = hydro_get_physical_density(p, cosmo);

  /* Deside whether we should form stars or not,
   * first we deterime if we have the correct over density
   * if that is true we calculate if either the maximum density
   * threshold is reached or if the metallicity dependent
   * threshold is reached, after this we calculate if the
   * temperature is appropriate */
  if (particle_density < rho_crit_times_min_over_den) return 0;

  /* In this case there are actually multiple possibilities
   * because we also need to check if the physical density exceeded
   * the appropriate limit */

  const double Z = p->chemistry_data.smoothed_metal_mass_fraction_total;
  double density_threshold_metal_dep;
  if (Z > 0) {
    density_threshold_metal_dep =
        starform->density_threshold * pow(Z * starform->Z0_inv, starform->n_Z0);
  } else {
    density_threshold_metal_dep = starform->density_threshold_max;
  }

  /* Calculate the maximum between both and convert to mass density instead of
   * number density*/
  const double density_threshold_current =
      min(density_threshold_metal_dep, starform->density_threshold_max) *
      phys_const->const_proton_mass;

  /* Check if it exceeded the maximum density */
  if (particle_density * p->chemistry_data.smoothed_metal_mass_fraction[0] <
      density_threshold_current)
    return 0;

  /* Calculate the temperature */
  const double temperature = cooling_get_temperature(phys_const, hydro_props,
                                                     us, cosmo, cooling, p, xp);

  /* Check the last criteria, if the temperature is satisfied */
  return (temperature < starform->Temperature_threshold);
}

/**
 * @brief Calculates if the gas particle gets converted
 *
 * @param the #engine
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 */
INLINE static int star_formation_convert_to_star(
    const struct engine* e, const struct star_formation* starform,
    const struct part* restrict p, struct xpart* restrict xp,
    const struct phys_const* const phys_const, const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling, const double dt_star,
    const int with_cosmology) {

  if (dt_star == 0.f) return 0;

  if (star_formation_potential_to_become_star(
          starform, p, xp, phys_const, cosmo, hydro_props, us, cooling)) {
    /* Get the pressure */
    
    const double pressure =
        starform->EOS_pressure_norm *
        pow(hydro_get_physical_density(p, cosmo) *
                p->chemistry_data.smoothed_metal_mass_fraction[0] /
                starform->EOS_density_norm / phys_const->const_proton_mass,
            starform->polytropic_index);
    

    double SFRpergasmass;
    if (hydro_get_physical_density(p, cosmo) <
        starform->KS_high_den_thresh * phys_const->const_proton_mass) {
      /* Calculate the star formation rate */
      SFRpergasmass = starform->SF_normalization *
                                   pow(pressure, starform->SF_power_law);
      //SFRpergasmass = starform->KS_normalization * pow(starform->Msunpsquaredpc,-starform->KS_power_law)*
      //    pow(hydro_gamma * starform->fgas / phys_const->const_newton_G * pressure, (starform->KS_power_law-1.f)/2.f );
    } else {
      SFRpergasmass = starform->SF_high_den_normalization *
                                   pow(pressure, starform->SF_high_den_power_law);
    }

    /* Store the SFR */
    xp->SFR = SFRpergasmass * p->mass;
    xp->SFRrate = SFRpergasmass;
    xp->sftracers_data.SFR = SFRpergasmass * p->mass;
    xp->sftracers_data.sSFR = SFRpergasmass;
    /* Calculate the propability of forming a star */
    const double prop = SFRpergasmass * dt_star;

    /* Calculate the seed */
    unsigned int seed = (p->id + e->ti_current) % 8191;

    /* Generate a random number between 0 and 1. */
    const double randomnumber = rand_r(&seed) * starform->inv_RAND_MAX;

    // message("Passed whole boundary thing! random number = %e, prop = %e
    // dt_star %e", randomnumber, prop,dt_star);

    /* Calculate if we form a star */
    return (prop > randomnumber);
  }

  /* Check if it is the first time steps after star formation */
  if (xp->SFR > 0.f) {
    if (with_cosmology) {
      xp->SFR = - cosmo->a;
      xp->SFRrate = 0.f;
      xp->sftracers_data.SFR = - cosmo->a;
      xp->sftracers_data.sSFR = 0.f;
    } else {
      xp->SFR = - e->time; 
      xp->sftracers_data.SFR = - e->time;
      xp->SFRrate = 0.f;
      xp->sftracers_data.sSFR = 0.f;
    }
  }

  return 0;
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
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
INLINE static void star_formation_copy_properties(
    const struct engine* e, const struct part* p, const struct xpart* xp,
    struct spart* sp, const struct star_formation* starform,
    const struct phys_const* const phys_const, const struct cosmology* cosmo,
    const int with_cosmology) {

  /* Store the current mass */
  sp->mass = p->mass;

  /* Store the current mass as the initial mass */
  sp->mass_init = p->mass;

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

  sp->new_star_flag = 1;

  message("A star has been formed!");
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param starform the star formation law properties to initialize
 *
 * */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, struct star_formation* starform) {

  /* Get the appropriate constant to calculate the
   * star formation constant */
  const double KS_const =
      phys_const->const_solar_mass /
      (1e6 * phys_const->const_parsec * phys_const->const_parsec) /
      phys_const->const_year;

  /* Get the Gravitational constant */
  const double G_newton = phys_const->const_newton_G;

  /* Get the surface density unit M_\odot / pc^2 */
  const double M_per_pc2 =
      phys_const->const_solar_mass /
      (phys_const->const_parsec * phys_const->const_parsec);

  starform->Msunpsquaredpc = M_per_pc2;
  /* Calculate inverse of RAND_MAX for the random numbers */
  starform->inv_RAND_MAX = 1.f / RAND_MAX;

  /* Conversion of number density from cgs */
  static const float dimension_numb_den[5] = {0, -3, 0, 0, 0};
  const double conversion_numb_density =
      1 / units_general_cgs_conversion_factor(us, dimension_numb_den);

  /* Quantities that have to do with the Normal Kennicutt-
   * Schmidt law will be read in this part of the code*/

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den =
      parser_get_param_double(parameter_file, "SchayeSF:thresh_MinOverDens");

  /* Read the critical temperature from the parameter file */
  starform->Temperature_threshold =
      parser_get_param_double(parameter_file, "SchayeSF:thresh_temp");

  /* Read the gas fraction from the file */
  starform->fgas = parser_get_param_double(parameter_file, "SchayeSF:fg");

  /* Read the normalization of the KS law in KS law units */
  starform->KS_normalization_MSUNpYRpKPC2 = parser_get_param_double(
      parameter_file, "SchayeSF:SchmidtLawCoeff_MSUNpYRpKPC2");

  /* Read the Kennicutt-Schmidt power law exponent */
  starform->KS_power_law =
      parser_get_param_double(parameter_file, "SchayeSF:SchmidtLawExponent");

  /* Calculate the power law of the star formation */
  starform->SF_power_law = (starform->KS_power_law - 1.f) / 2.f;

  /* Give the Kennicutt-Schmidt law the same units as internal units */
  starform->KS_normalization =
      starform->KS_normalization_MSUNpYRpKPC2 * KS_const;

  /* Calculate the starformation prefactor with most terms */
  starform->SF_normalization =
      starform->KS_normalization * pow(M_per_pc2, -starform->KS_power_law) *
      pow(hydro_gamma * starform->fgas / G_newton, starform->SF_power_law);

  /* Read the high density Kennicutt-Schmidt power law exponent */
  starform->KS_high_den_power_law = parser_get_param_double(
      parameter_file, "SchayeSF:SchmidtLawHighDensExponent");

  /* Read the high density criteria for the KS law in number density per cm^3 */
  starform->KS_high_den_thresh_HpCM3 = parser_get_param_double(
      parameter_file, "SchayeSF:SchmidtLawHighDens_thresh_HpCM3");

  /* Transform the KS high density criteria to simulation units */
  starform->KS_high_den_thresh =
      starform->KS_high_den_thresh_HpCM3 * conversion_numb_density;

  /* Calculate the SF high density power law */
  starform->SF_high_den_power_law =
      (starform->KS_high_den_power_law - 1.f) / 2.f;

  /* Load the equation of state for this model */
  starform->polytropic_index = parser_get_param_double(
      parameter_file, "SchayeSF:EOS_Jeans_GammaEffective");
  starform->EOS_temperature_norm = parser_get_param_double(
      parameter_file, "SchayeSF:EOS_Jeans_TemperatureNorm_K");
  starform->EOS_density_norm_HpCM3 = parser_get_param_double(
      parameter_file, "SchayeSF:EOS_JEANS_DensityNorm_HpCM3");
  starform->EOS_density_norm =
      starform->EOS_density_norm_HpCM3 * conversion_numb_density;

  /* Calculate the EOS pressure normalization */
  starform->EOS_pressure_norm = starform->EOS_density_norm *
                                starform->EOS_temperature_norm *
                                phys_const->const_boltzmann_k;

  const double EOS_high_den_pressure =
      starform->EOS_pressure_norm *
      pow(starform->KS_high_den_thresh / starform->EOS_density_norm,
          starform->polytropic_index);

  /* Calculate the KS high density normalization */
  starform->KS_high_den_normalization =
      starform->KS_normalization *
      pow(M_per_pc2, starform->KS_high_den_power_law - starform->KS_power_law) *
      pow(hydro_gamma * starform->fgas / G_newton * EOS_high_den_pressure,
          (starform->KS_power_law - starform->KS_high_den_power_law) / 2.f);

  /* Calculate the SF high density normalization */
  starform->SF_high_den_normalization =
      starform->KS_high_den_normalization *
      pow(M_per_pc2, -starform->KS_high_den_power_law) *
      pow(hydro_gamma * starform->fgas / G_newton,
          starform->SF_high_den_power_law);

  /* Use the Schaye (2004) metallicity dependent critical density
   * to form stars. */
  /* Read the normalization of the metallicity dependent critical
   * density*/
  starform->density_threshold_HpCM3 =
      parser_get_param_double(parameter_file, "SchayeSF:thresh_norm_HpCM3");

  starform->density_threshold =
      starform->density_threshold_HpCM3 * conversion_numb_density;

  /* Read the scale metallicity Z0 */
  starform->Z0 = parser_get_param_double(parameter_file, "SchayeSF:MetDep_Z0");

  /* Read the power law of the critical density scaling */
  starform->n_Z0 =
      parser_get_param_double(parameter_file, "SchayeSF:MetDep_SFthresh_Slope");

  /* Read the maximum allowed density for star formation */
  starform->density_threshold_max_HpCM3 =
      parser_get_param_double(parameter_file, "SchayeSF:thresh_max_norm_HpCM3");

  starform->density_threshold_max =
      starform->density_threshold_max_HpCM3 * conversion_numb_density;
  /* Claculate 1 over the metallicity */
  starform->Z0_inv = 1 / starform->Z0;
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  message("Star formation law is Schaye and Dalla Vecchia (2008)");
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
      "temperature = "
      "%e K",
      starform->polytropic_index, starform->EOS_density_norm_HpCM3,
      starform->EOS_temperature_norm);
  message("Density threshold is given by Schaye (2004)");
  message(
      "the normalization of the density threshold is given by"
      " %e #/cm^3, with metallicity slope of %e, and metallicity normalization"
      " of %e, the maximum density threshold is given by %e #/cm^3",
      starform->density_threshold_HpCM3, starform->n_Z0, starform->Z0,
      starform->density_threshold_max_HpCM3);
  message("Temperature threshold is given by Dalla Vecchia and Schaye (2012)");
  message("The temperature threshold is given by: %e K",
          starform->Temperature_threshold);
  message("DEBUG: PRINT COMPLETE STRUCT");
  message("KS_normalization = %e", starform->KS_normalization);
  // message("KS_normalization_MSUNpYRpKPC2 = %e",
  // starform->KS_normalization_MSUNpYRpKPC2);  message("KS_power_law = %e",
  // starform->KS_power_law);  message("KS_high_den_power_law = %e",
  // starform->KS_high_den_power_law);
  message("KS_high_den_thresh = %e", starform->KS_high_den_thresh);
  message("KS_high_den_normalization = %e",
          starform->KS_high_den_normalization);
  // message("min_over_den = %e", starform->min_over_den);
  // message("Temperature_threshold = %e", starform->Temperature_threshold);
  // message("fgas = %e", starform->fgas);
  // message("SF_power_law = %e", starform->SF_power_law);
  message("SF_normalization = %e", starform->SF_normalization);
  // message("SF_high_den_power_law = %e", starform->SF_high_den_power_law);
  message("SF_high_den_normalization = %e",
          starform->SF_high_den_normalization);
  // message("inv_RAND_MAX = %e", starform->inv_RAND_MAX);
  message("density_threshold = %e", starform->density_threshold);
  // message("density_threshold_HpCM3 = %e", starform->density_threshold_HpCM3);
  message("density_threshold_max = %e", starform->density_threshold_max);
  // message("density_threshold_max_HpCM3 = %e",
  // starform->density_threshold_max_HpCM3);  message("Z0 = %e", starform->Z0);
  // message("Z0_inv = %e", starform->Z0_inv);
  // message("n_Z0 = %e", starform->n_Z0);
  // message("polytropic_index = %e", starform->polytropic_index);
  message("EOS_pressure_norm = %e", starform->EOS_pressure_norm);
  // message("EOS_temperature_norm = %e", starform->EOS_temperature_norm);
  message("EOS_density_norm = %e", starform->EOS_density_norm);
  // message("EOS_density_norm_HpCM3 = %e", starform->EOS_density_norm_HpCM3);
  message("SFR normalization = %e",
          pow(starform->EOS_pressure_norm, starform->polytropic_index / 5.f) *
              starform->SF_normalization);
}

#endif /* SWIFT_SCHAYE_STARFORMATION_H */
