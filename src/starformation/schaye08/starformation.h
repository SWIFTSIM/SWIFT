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

/* Local includes */
#include "cosmology.h"
#include "physical_constants.h"
#include "units.h"
#include "parser.h"
#include "equation_of_state.h"

/* Starformation struct */
struct star_formation {
  
  /*! Normalization of the KS star formation law */
  double A;

  /*! Slope of the KS law */
  double nks;

  /*! Critical overdensity */
  double Delta_crit;

  /*! Critical temperature */
  double T_crit;

  /*! Ratio of the specific heats */
  double gamma;

  /*! gas fraction */
  double fg;

  /*! Star formation law slope */
  double nstar;

  /*! star formation normalization of schaye+08 */
  double Astar;
  
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
static int starformation_potential_to_become_star(
    const struct star_formation* starform, const struct parts* p,
    const struct xparts* xp, const struct phys_const* const phys_const,
    const struct cosmology* cosmo){

  /* Read the critical overdensity factor and the critical density of 
   * the universe to determine the critical density to form stars*/
  const double rho_crit = cosmo->critical_density*starform->Delta_crit; 
  
  /* Calculate the internal energy using the density and entropy */
  /* Ask Matthieu about p->entropy vs xp->entropy_full */
  const double internal_energy = gas_internal_energy_from_entropy(
  p->rho, p->entropy)

  /* Calculate the temperature over mu of the gas */
  const double T_over_mu = (starform->gamma - 1)*phys_const->const_proton_mass
  /phys_const->const_boltzmann_k;


  
  /* Deside whether we should form stars or not */
  if ((p->rho > rho_crit ) && (tempp < starform->T_crit)) {
    return 1;
  } else {
    return 0;
  }
}

/*
 * @brief Calculate if the gas particle is converted 
 *
 * */
static void starformation_convert_to_gas( 
    ){

}

/* 
 * @brief initialization of the star formation law 
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param starform the star formation law properties to initialize
 *
 * */
static void starformation_init_backend(
  struct swift_params* parameter_file, const struct phys_const* phys_const,
  const struct unit_system* us, const struct space* s,
  const struct star_formation* starform) {
  
  /* Default values for the normalization and the power law */
  static const double normalization_default = 2.5e-4;
  static const double KS_power_law_default = 1.4;

  /* Default value for the heat capacity ratio gamma */
  static const double gamma_default = 5.f/3.f;

  /* Read the critical density contrast from the parameter file*/
  starform->Delta_crit = parser_get_param_double(parameter_file, 
  "SchayeSF:Delta_crit");

  /* Read the critical temperature from the parameter file */
  starform->T_crit = parser_get_param_double(parameter_file,
  "SchayeSF:T_crit");

  /* Read the gas fraction from the file */
  starform->fg = parser_get_param_double(parameter_file,
  "SchayeSF:fg");

  /* Read the normalization */
  const double normalization = parser_get_opt_param_double(
  parameter_file, "SchayeSF:A", normalization_default);

  /* Read the Kennicutt-Schmidt power law exponent */
  starform->nks = parser_get_opt_param_double(
  parameter_file, "SchayeSF:nks", KS_power_law_default);

  /* Read the heat capacity ratio gamma */
  starform->gamma = parser_get_opt_param_double(
  parameter_file, "SchayeSF:gamma", gamma_default); 

  /* Calculate the power law of the star formation */
  starform->nstar = (starform->nks - 1.f)/2.f;

  /* Get the appropriate constant to calculate the 
   * star formation constant */ 
  const double G_newton = phys_const->const_newton_G;

  
}

/* @brief Prints the used parameters of the star formation law 
 *
 * @param starform the star formation law properties.
 * */
static void starformation_print_backend(
    const struct star_formation* starform){ 

  message("Star formation law is Schaye and Dalla Vecchia (2008)"
  " with properties, normalization = %e, slope of the Kennicutt"
  "-Schmidt law = %e, gamma = %e, gas fraction = %e, critical "
  "density = %e and critical temperature = %e", starform->A, 
  starform->nks, starform->gamma, starform->fg, starform->rho_crit,
  starform->T_crit);

}


#endif /* SWIFT_SCHAYE_STARFORMATION_H */


