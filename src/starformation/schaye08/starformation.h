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
#include "cosmology.h"
#include "physical_constants.h"
#include "units.h"
#include "parser.h"
#include "equation_of_state.h"
#include "part.h"
#include "hydro.h"
#include "cooling.h"
#include "adiabatic_index.h"
#include "cell.h" 

/* Starformation struct */
struct star_formation {
  
  /*! Normalization of the KS star formation law */
  double KS_normalization;

  /*! Slope of the KS law */
  double KS_power_law;

  /*! Slope of the high density KS law */
  double KS_high_den_power_law;

  /*! KS law High density threshold */
  double KS_high_den_thresh;

  /*! KS high density normalization */
  double KS_high_den_normalization;

  /*! Critical overdensity */
  double min_over_den;

  /*! Critical temperature */
  double T_crit;

  /*! gas fraction */
  double fgas;

  /*! Star formation law slope */
  double SF_power_law;

  /*! star formation normalization of schaye+08 */
  double SF_normalization;

  /*! star formation high density slope */
  double SF_high_den_power_law;

  /*! Star formation high density normalization */
  double SF_high_den_normalization;

  /*! Inverse of RAND_MAX */
  double inv_RAND_MAX;

  /*! Critical density to form stars */
  double den_crit;

  /*! Maximum critical density to form stars */
  double den_crit_max;

  /*! Scaling metallicity */
  double Z0;

  /*! one over the scaling metallicity */
  double Z0_inv;

  /*! critical density Metallicity power law */
  double n_Z0;

  /*! Normalization of critical SF density of Schaye (2004) */
  double den_crit_star;

  /*! Metallicity dependent critical density from Schaye (2004) */
  int schaye04;
  
  /*! Polytropic index */
  double polytropic_index;

  /*! EOS pressure norm */
  double EOS_pressure_norm;

  /*! EOS density norm */
  double EOS_den0;
};

/*
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles
 * @param xp the additional properties of the gas particles
 * @param phys_const the physical constants in internal units
 * @param cosmo the cosmological parameters and properties
 *
 * */
INLINE static int star_formation_potential_to_become_star(
    const struct star_formation* starform, const struct part* restrict p,
    const struct xpart* restrict xp, const struct phys_const* const phys_const,
    const struct cosmology* cosmo){

  /* Read the critical overdensity factor and the critical density of 
   * the universe to determine the critical density to form stars*/
  const double rho_crit = cosmo->critical_density*starform->min_over_den; 
  
  /* double tempp = cooling_get_temperature() */
  double tempp;

  
  /* Deside whether we should form stars or not,
   * first we deterime if we have the correct over density
   * if that is true we calculate if either the maximum density
   * threshold is reached or if the metallicity dependent 
   * threshold is reached, after this we calculate if the 
   * temperature is appropriate */
  if (p->rho > rho_crit ) {
    /* In this case there are actually multiple possibilities
     * because we also need to check if the physical density exceeded
     * the appropriate limit */

    /* Check if it exceeded the maximum density */
    if (p->rho > starform->den_crit_max) {
      /* double tempp = cooling_get_temperature() */
      tempp = 5e3;
      /* Check the last criteria, if the temperature is satisfied */
      if (tempp < starform->T_crit) {
        return 1;
      } else {
        return 0;
      }
    /* If the maximum density is not exceeded check if the redshift dependent
     * one is exceeded */
    } else {
      /* NEED TO USE A PROPER WAY OF FINDING Z */
      double Z = 0.002;
      double den_crit_current = starform->den_crit * pow(Z*starform->Z0_inv, starform->n_Z0);
      if (p->rho > den_crit_current) {
        /* double tempp = cooling_get_temperature() */
        tempp = 5e3;
        /* Check the last criteria, if the temperature is satisfied */
        if (tempp < starform->T_crit) {
          return 1;
        } else {
          return 0;
        }
      } else {
        return 0;
      }
    }
  } else {
    return 0;
  }
}

/*
 * @brief Calculates if the gas particle gets converted
 *
 */
INLINE static int star_formation_convert_to_star(
    const struct star_formation* starform, const struct part* restrict p,
    const struct xpart* restrict xp, const struct phys_const* const phys_const,
    const struct cosmology* cosmo) {

  if (star_formation_potential_to_become_star(starform, p, xp, phys_const, cosmo)){
    /* Get the pressure */
    const double pressure = starform->EOS_pressure_norm 
    * pow(p->rho/starform->EOS_den0, starform->polytropic_index);
  
    /* Calculate the propability of forming a star */ 
    const double prop = starform->SF_normalization * pow(pressure,
    starform->SF_power_law) * p->time_bin;

    /* Calculate the seed */
    unsigned int seed = p->id;
    
    /* Generate a random number between 0 and 1. */
    const double randomnumber = rand_r(&seed)*starform->inv_RAND_MAX; 

    /* Calculate if we form a star */
    if (prop > randomnumber) {
      return 1;
    } else {
      return 0;
    }

  } else {
    return 0;
  }
}

INLINE static void star_formation_copy_properties(
    struct engine *e, struct cell *c, struct part* p,
    struct xpart* xp, const struct star_formation* starform, 
    const struct phys_const* const phys_const, const struct cosmology* cosmo) {
  
  struct spart *sp = cell_convert_part_to_spart(e, c, p, xp);
  sp->mass = p->mass;
  sp->mass_init = p->mass;
  message("Copy Properties");

}

/*
int banana(...) {

  if(star_formation_potential_to_become_ater(....) {

  //draw ranfom number
  //  return 1 or 0

  }else
  return 0;
  }
 */

/* 
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
  const double KS_const = phys_const->const_kennicutt_schmidt_units;

  /* Get the Gravitational constant */
  const double G_newton = phys_const->const_newton_G;

  /* Get the surface density unit M_\odot / pc^2 */
  const double M_per_pc2 = phys_const->const_solar_mass_per_parsec2;
  
  /* Calculate inverse of RAND_MAX for the random numbers */
  starform->inv_RAND_MAX = 1.f / RAND_MAX;

  /* Quantities that have to do with the Normal Kennicutt-
   * Schmidt law will be read in this part of the code*/

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den = parser_get_param_double(parameter_file, 
  "SchayeSF:thresh_MinOverDens");

  /* Read the critical temperature from the parameter file */
  starform->T_crit = parser_get_param_double(parameter_file,
  "SchayeSF:thresh_temp");

  /* Read the gas fraction from the file */
  starform->fgas = parser_get_param_double(parameter_file,
  "SchayeSF:fg");

  /* Read the normalization of the KS law in KS law units */
  const double normalization_MSUNpYRpKPC2 = parser_get_param_double(
  parameter_file, "SchayeSF:SchmidtLawCoeff_MSUNpYRpKPC2");

  /* Read the Kennicutt-Schmidt power law exponent */
  starform->KS_power_law = parser_get_param_double(
  parameter_file, "SchayeSF:SchmidtLawExponent");

  /* Calculate the power law of the star formation */
  starform->SF_power_law = (starform->KS_power_law - 1.f)/2.f;
  
  /* Give the Kennicutt-Schmidt law the same units as internal units */
  starform->KS_normalization = normalization_MSUNpYRpKPC2 * KS_const;
  
  /* Calculate the starformation prefactor with most terms */
  starform->SF_normalization = starform->KS_normalization * pow(M_per_pc2, -starform->KS_power_law) * 
  pow( hydro_gamma * starform->fgas / G_newton, starform->SF_power_law);

  /* Read the high density Kennicutt-Schmidt power law exponent */
  starform->KS_high_den_power_law = parser_get_param_double(
  parameter_file, "SchayeSF:SchmidtLawHighDensExponent");
  
  /* Read the high density criteria for the KS law in number density per cm^3 */
  const double KS_high_den_thresh_HpCM3 = parser_get_param_double(
  parameter_file, "SchayeSF:SchmidtLawHighDens_thresh_HpCM3");

  /* Transform the KS high density criteria to simulation units */
  starform->KS_high_den_thresh = KS_high_den_thresh_HpCM3 * UNIT_CONV_NUMBER_DENSITY;

  /* Calculate the SF high density power law */
  starform->SF_high_den_power_law = (starform->KS_high_den_power_law - 1.f)/2.f;

  /* Calculate the KS high density normalization */
  starform->KS_high_den_normalization = starform->KS_normalization * pow( M_per_pc2,
  starform->KS_high_den_power_law - starform->KS_power_law) * pow( hydro_gamma * 
  starform->fgas / G_newton * 1337.f, (starform->KS_power_law 
  - starform->KS_high_den_power_law)/2.f);

  /* Calculate the SF high density normalization */
  starform->SF_high_den_normalization = starform->KS_high_den_normalization 
  * pow(M_per_pc2, -starform->KS_high_den_power_law) * pow( hydro_gamma  
  * starform->fgas / G_newton, starform->SF_high_den_power_law);

  /* Read what kind of critical density we need to use
   * Schaye (2004) is metallicity dependent critical SF density*/
  starform->schaye04 = parser_get_param_double(
  parameter_file, "SchayeSF:Schaye2004");

  if (!starform->schaye04) {
    /* In the case that we do not use the Schaye (2004) critical
     * density to form stars but a constant value */
    starform->den_crit = parser_get_param_double(
    parameter_file, "SchayeSF:thresh_norm_HpCM3");
    starform->Z0 = 0.002;
    starform->n_Z0 = 0.0;
  } else {
    /* Use the Schaye (2004) metallicity dependent critical density
     * to form stars. */
    /* Read the normalization of the metallicity dependent critical 
     * density*/
    starform->den_crit = parser_get_param_double( 
    parameter_file, "SchayeSF:thresh_norm_HpCM3") *
    UNIT_CONV_NUMBER_DENSITY;

    /* Read the scale metallicity Z0 */
    starform->Z0 = parser_get_param_double(
    parameter_file, "SchayeSF:MetDep_Z0");

    /* Read the power law of the critical density scaling */
    starform->n_Z0 = parser_get_param_double(
    parameter_file, "SchayeSF:MetDep_SFthresh_Slope");

    /* Read the maximum allowed density for star formation */
    starform->den_crit_max = parser_get_param_double(
    parameter_file, "SchayeSF:thresh_max_norm_HpCM3");

  }
  /* Conversion of number density from cgs */
  static const float dimension_numb_den[5] = {0, -3, 0, 0, 0};
  const double conversion_numb_density = 1.f/
  units_general_cgs_conversion_factor(us, dimension_numb_den);

  /* Claculate 1 over the metallicity for speed up */
  starform->Z0_inv = 1/starform->Z0;

  /* Calculate the prefactor that is always common */
  /* !!!DONT FORGET TO DO THE CORRECT UNIT CONVERSION!!!*/
  starform->den_crit_star = starform->den_crit / pow(starform->Z0,
  starform->n_Z0) * conversion_numb_density;

  /* Load the equation of state for this model */
  starform->polytropic_index = parser_get_param_double(
  parameter_file, "SchayeSF:EOS_Jeans_GammaEffective");
  starform->EOS_pressure_norm = parser_get_param_double(
  parameter_file, "SchayeSF:EOS_Jeans_PressureNorm");
  starform->EOS_den0 = parser_get_param_double(
  parameter_file, "SchayeSF:EOS_JEANS_den0");

  
}

/* @brief Prints the used parameters of the star formation law 
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform){ 

  message("Star formation law is Schaye and Dalla Vecchia (2008)"
  " with properties, normalization = %e, slope of the Kennicutt"
  "-Schmidt law = %e, gas fraction = %e, critical "
  "density = %e and critical temperature = %e", starform->KS_normalization, 
  starform->KS_power_law, starform->fgas, starform->den_crit,
  starform->T_crit);
  if (!starform->schaye04) {
    message("Density threshold to form stars is constant and given"
    "by a density of %e", starform->den_crit_star);
  } else {
    message("Density threshold to form stars is given by Schaye "
    "(2004), the normalization of the star formation law is given by"
    " %e, with metallicity slope of %e, and metallicity normalization"
    "of %e", starform->den_crit_star, starform->n_Z0, starform->Z0);
  }

}


#endif /* SWIFT_SCHAYE_STARFORMATION_H */
