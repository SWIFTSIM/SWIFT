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

/* Starformation struct */
struct star_formation {
  
  /* Normalization of the KS star formation law */
  double A;

  /* Slope of the KS law */
  double nks;

  /* Critical density */
  double rho_crit;

  /* Critical temperature */
  double T_crit;

  /* Ratio of the specific heats */
  double gamma;

  /* gas fraction */
  double fg;

  /* Star formation law slope */
  double nstar;

  /* star formation normalization of schaye+08 */
  double Astar;
  
};


/* Brief initialization of the star formation law 
 * 
 * */
static void starformation_init_backend(
    
    ) {
  /* Default values for the normalization and the power law */
  static const double normalization_default = 2.5e-4;
  static const double KS_power_law = 1.4;


  
}

/* @brief Prints the used parameters of the star formation law 
 *
 * */
static void starformation_print_backend(

    ){ 
  message("Star formation law is Schaye and Dalla Vecchia (2008)"
  " with properties, normalization = %e, slope of the Kennicutt"
  "-Schmidt law = %e, gamma = %e, gas fraction = %e, critical "
  "density = %e and critical temperature = %e", starform->A, 
  starform->nks, starform->gamma, starform->fg, starform->rho_crit,
  starform->T_crit);

}


#endif /* SWIFT_SCHAYE_STARFORMATION_H */


