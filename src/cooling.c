/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "potentials.h"

/**
 * @brief Initialises the cooling properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param us The current internal system of units
 * @param cooling  The cooling  properties to initialize
 */
void cooling_init(const struct swift_params* parameter_file,
                    struct UnitSystem* us,
                    struct cooling_data* cooling) {

#ifdef CONST_COOLING
  cooling->const_cooling.lambda = parser_get_param_double(parameter_file, "Cooling:lambda");
  cooling->const_cooling.min_energy = parser_get_param_double(parameter_file, "Cooling:min_energy");
  cooling->const_cooling.cooling_tstep_mult = parser_get_param_double(parameter_file, "Cooling:cooling_tstep_mult");
#endif /* CONST_COOLING */

}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param  cooling The cooling properties.
 */
void cooling_print(const struct cooling_data* cooling) {

#ifdef CONST_COOLING
  message(
      "Cooling properties are (lambda, min_energy, tstep multiplier) %g %g %g ",
      cooling->const_cooling.lambda,
      cooling->const_cooling.min_energy
      cooling->const_cooling.cooling_tstep_mult);
#endif /* CONST_COOLING */
}
int update_entropy(const struct cooling_data* cooling,
		   const struct phys_const* const phys_const, struct part* p, float dt){

  /*updates the entropy of a particle after integrating the cooling equation*/
  int status == 0;
  float u_old;
  float u_new;
  float new_entropy;
  float old_entropy = p->entropy;
  float rho = p->rho;

  u_old = old_entropy/(GAMMA_MINUS1) * pow(rho,GAMMA_MINUS1);
  status = calculate_new_thermal_energy(u_old,&u_new,dt,cooling):
  
  if (status == 0){
    new_entropy = u_new/pow(rho,GAMMA_MINUS1) * GAMMA_MINUS1;
    p->entropy = new_entropy
  }
  else
    message("Error with cooling, particle's entropy has not been updated");
 
  return status;
}

#ifdef CONST_COOLING
int calculate_new_thermal_energy(float u_old, float* u_new, float dt, const struct cooling_data* cooling){

  //This function integrates the cooling equation, given the initial thermal energy and the timestep dt.
  //Returns 0 if successful and 1 if not
  int status = 0;
  float du_dt = cooling->const_cooling.lambda;
  float u_floor = cooling->const_cooling.min_energy;
  if (u_old - du_dt*dt > min_energy):
    *u_new = u_old - du_dt*dt;
  else:
    *u_new = min_energy

  return status}

#endif /*CONST_COOLING
  

  
