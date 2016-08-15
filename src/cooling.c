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
#include "cooling.h"
#include "hydro.h"
#include "adiabatic_index.h"

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

#ifdef CREASEY_COOLING
  cooling->creasey_cooling.lambda = parser_get_param_double(parameter_file, "CreaseyCooling:lambda");
  cooling->creasey_cooling.min_energy = parser_get_param_double(parameter_file, "CreaseyCooling:min_energy");
  cooling->creasey_cooling.cooling_tstep_mult = parser_get_param_double(parameter_file, "CreaseyCooling:cooling_tstep_mult");
#endif /* CREASEY_COOLING */


}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param  cooling The cooling properties.
 */
void cooling_print(const struct cooling_data* cooling) {

#ifdef CONST_COOLING */
  message(
      "Cooling properties are (lambda, min_energy, tstep multiplier) %g %g %g ",
      cooling->const_cooling.lambda,
      cooling->const_cooling.min_energy,
      cooling->const_cooling.cooling_tstep_mult);
#endif /* CONST_COOLING */

#ifdef CREASEY_COOLING
  message(
      "Cooling properties for Creasey cooling are (lambda, min_energy, tstep multiplier) %g %g %g ",
      cooling->creasey_cooling.lambda,
      cooling->creasey_cooling.min_energy,
      cooling->creasey_cooling.cooling_tstep_mult);
#endif /* CREASEY_COOLING */
}

void update_entropy(const struct cooling_data* cooling,
		   const struct phys_const* const phys_const, struct part* p, 
		   double dt){

  /*updates the entropy of a particle after integrating the cooling equation*/
  float u_old;
  float u_new;
  float new_entropy;
  float old_entropy = p->entropy;
  float rho = p->rho;

  //  u_old = old_entropy/(GAMMA_MINUS1) * pow(rho,GAMMA_MINUS1);
  u_old = hydro_get_internal_energy(p,0); // dt = 0 because using current entropy
  u_new = calculate_new_thermal_energy(u_old,rho,dt,phys_const,cooling);
  new_entropy = u_new*pow_minus_gamma_minus_one(rho) * hydro_gamma_minus_one;
  p->entropy = new_entropy;
}

/*This function integrates the cooling equation, given the initial
  thermal energy, density and the timestep dt. Returns the final internal energy*/

float calculate_new_thermal_energy(float u_old, float rho, float dt, 
				   const struct phys_const* const phys_const, 
				   const struct cooling_data* cooling){
#ifdef CONST_COOLING
  /*du/dt = -lambda, independent of density*/
  float du_dt = -cooling->const_cooling.lambda;
  float u_floor = cooling->const_cooling.min_energy;
  float u_new;
  if (u_old - du_dt*dt > u_floor){
    u_new = u_old + du_dt*dt;
  }
  else{
    u_new = u_floor;
  }
#endif /*CONST_COOLING*/

#ifdef CREASEY_COOLING
  /* rho*du/dt = -lambda*n_H^2, where rho = m_p * n_H*/
  float m_p = phys_const->const_proton_mass;
  float n_H = rho / m_p;
  float lambda = cooling->creasey_cooling.lambda;
  float du_dt = -lambda * n_H * n_H / rho;
  float u_floor = cooling->creasey_cooling.min_energy;
  float u_new;
  if (u_old - du_dt*dt > u_floor){
    u_new = u_old + du_dt*dt;
  }
  else{
    u_new = u_floor;
  }
#endif /*CREASEY_COOLING*/
  return u_new;
}


  

  
