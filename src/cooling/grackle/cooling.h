/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_GRACKLE_H
#define SWIFT_COOLING_GRACKLE_H

/**
 * @file src/cooling/none/cooling.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* include the grackle wrapper */
#include "grackle_wrapper.h"





/*! This function computes the new entropy due to the cooling,
 *  between step t0 and t1.
 */

static INLINE double do_cooling_grackle(double energy, double density, double dtime, double *ne, double Z, double a_now)
{

  

  /*********************************************************************
   call to the main chemistry solver
   *********************************************************************/
  
  if (wrap_do_cooling(density, &energy, dtime,Z, a_now) == 0) {
    error("Error in do_cooling.\n");
    return 0;
  }
  
    
  return energy;


}


/**
 * @brief Apply the cooling function to a particle.
 *
 * We do nothing.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {
    



  /* Get current internal energy (dt=0) */
  const float u_old = hydro_get_internal_energy(p);
  /* Get current density */
  const float rho = hydro_get_density(p);
  /* Actual scaling fractor */
  const float a_now =  1. / (1. + cooling->GrackleRedshift); ;			/*  must be chaged !!! */

  double ne,Z;

  Z =  0.02041;					/* 0.02041 (= 1 Zsun in Grackle v2.0, but = 1.5761 Zsun in Grackle v2.1) */
  ne=  0.0;					/* mass fraction of eletron */          	/* useless for GRACKLE_CHEMISTRY = 0 */
 
 
  float u_new;
  float delta_u;
  
  u_new = do_cooling_grackle(u_old, rho, dt, &ne, Z, a_now);
  //u_new = u_old * 0.99;





  
  //if (u_new < 0)
  //if (p->id==50356)
  //  printf("WARNING !!! ID=%llu  u_old=%g  u_new=%g rho=%g dt=%g ne=%g Z=%g a_now=%g\n",p->id,u_old,u_new,rho,dt,ne,Z,a_now);

  
  delta_u = u_new - u_old;

  /* record energy lost */
  xp->cooling_data.radiated_energy += -delta_u * hydro_get_mass(p);  
  
 
  
  /* Update the internal energy */
  hydro_set_internal_energy_dt(p, delta_u / dt);  
    
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const struct part* restrict p) {

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_init_part(
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
static INLINE void cooling_init_backend(
    const struct swift_params* parameter_file,
    const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {
    


  char cloudytable[200];
  double units_density, units_length, units_time;
  int grackle_chemistry;
  int UVbackground;


  parser_get_param_string(parameter_file,"GrackleCooling:GrackleCloudyTable",cooling->GrackleCloudyTable);
  cooling->UVbackground = parser_get_param_int(parameter_file, "GrackleCooling:UVbackground");
  cooling->GrackleRedshift = parser_get_param_double(parameter_file, "GrackleCooling:GrackleRedshift");
  cooling->GrackleHSShieldingDensityThreshold = parser_get_param_double(parameter_file, "GrackleCooling:GrackleHSShieldingDensityThreshold");

  // FIXME : Why a strcpy ?
  strcpy(cloudytable,cooling->GrackleCloudyTable);
  
  
  UVbackground  =  cooling->UVbackground;
  grackle_chemistry = 0;			/* forced to be zero : read table */ 
  
  units_density = us->UnitMass_in_cgs/pow(us->UnitLength_in_cgs,3);
  units_length  = us->UnitLength_in_cgs;
  units_time    = us->UnitTime_in_cgs;
    
    

  printf("          ***************************************\n");
  printf("          initializing grackle cooling function\n");
  printf("          \n");
  
  printf("          CloudyTable                        = %s\n",cloudytable);
  printf("          UVbackground                       = %d\n",UVbackground);
  printf("          GrackleRedshift                    = %g\n",cooling->GrackleRedshift);
  printf("          GrackleHSShieldingDensityThreshold = %g\n",cooling->GrackleHSShieldingDensityThreshold);
  

  if(wrap_init_cooling(cloudytable,UVbackground,units_density, units_length, units_time, grackle_chemistry) != 1)
    {
      fprintf(stderr,"Error in initialize_chemistry_data.");
      exit(-1);
    }
  
  printf("          \n");
  printf("          ***************************************\n");



    
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle'.");
}

#endif /* SWIFT_COOLING_GRACKLE_H */
