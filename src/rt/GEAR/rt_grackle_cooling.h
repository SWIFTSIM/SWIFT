/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GEAR_GRACKLE_COOLING_H
#define SWIFT_RT_GEAR_GRACKLE_COOLING_H

#include "chemistry.h"

/**
 * @file src/rt/GEAR/rt_grackle_cooling.h
 * @brief header file for the GEAR M1 closure radiative transfer scheme
 * kiara grackle cooling related functions.
 */

/**
 * @brief Returns the value of G0 for given particle p
 *
 * @param p Pointer to the particle data.
 * @param cooling The properties of the cooling function.
 *
 */
__attribute__((always_inline)) INLINE float cooling_compute_G0(
                const struct part *restrict p,
                const struct rt_props* rt_props) {

      float G0 = 0.f;
      /* Determine ISRF in Habing units based on chosen method */
      if (rt_props->G0_computation_method==0) {
          G0 = 0.f;
      }
      else if (rt_props->G0_computation_method==1) {
          G0 = p->chemistry_data.local_sfr_density * rt_props->G0_factor1;
      }
      else if (rt_props->G0_computation_method==2) {
          G0 = p->group_data.ssfr * rt_props->G0_factor2;
      }
      else if (rt_props->G0_computation_method==3) {
          if (p->group_data.ssfr > 0.) {
              G0 = p->group_data.ssfr * rt_props->G0_factor2;
          }
          else {
              G0 = p->chemistry_data.local_sfr_density * rt_props->G0_factor1;
          }
      }
      return G0;
}

/**
 * @brief get the densities of all species and electrons.
 *
 * @param p particle to use
 * @param rho particle physical density
 * @param species_densities array to write densities in
 *
 **/
__attribute__((always_inline)) INLINE static void
rt_tchem_get_species_densities(const struct part* restrict p, gr_float rho,
                               gr_float species_densities[RT_N_SPECIES],
                               gr_float species_extra[rt_species_extra_count],
                               const struct rt_props* rt_props) {
    
    /* nHII = rho_HII / m_p
     * nHeII = rho_HeII / 4 m_p
     * nHeIII = rho_HeIII / 4 m_p
     * ne = nHII + nHeII + 2 * nHeIII
     * But: it is grackle convention to use rho_e = n_e * m_p */
  
  for (enum rt_species species = 0; species < rt_species_count; species++){
        species_densities[species] = p->rt_data.tchem.mass_fraction[species] * rho;
  }

#ifdef SWIFT_RT_GRACKLE_DUST
  /* Load dust and metal info */ 
  species_extra[rt_species_dust] = p->rt_data.cooling.dust_mass / p->mass * rho;
  species_extra[rt_species_SNe_ThisTimeStep] = p->feedback_data.SNe_ThisTimeStep;
  //if( chemistry_get_total_metal_mass_fraction_for_cooling(p)>0.f) message("Zsm= %g Zp= %g Z= %g Zd= %g",chemistry_get_total_metal_mass_fraction_for_cooling(p), p->chemistry_data.metal_mass_fraction_total, species_densities[19], species_densities[20]);
  
  /* Determine ISRF in Habing units based on chosen method */
  species_extra[rt_species_isrf_habing] = cooling_compute_G0(p, rt_props);
  
  /* Load gas metallicities NEED TO CHANGE TO RIGHT COUNT  */
  int i = 0;
  for (enum rt_species_extra species = rt_species_He_gas; species < rt_species_He_dust; species++){
        species_extra[species] = p->chemistry_data.metal_mass_fraction[i] * rho;
        i++;
  }
  
  int i = 0;
  for (enum rt_species_extra species = rt_species_He_dust; species < rt_species_extra_count; species++){
        species_extra[species] = p->rt_data.cooling.dust_mass_fraction[i] * rho;
        i++;
  }


#endif

}

/**
 * @brief copy the grackle data to particle data for grackle chemistry mode 1.
 *
 * @param p particle to use
 * @param rho particle physical density
 * @param particle_grackle_data The grackle_field_data structure from grackle.
 **/
__attribute__((always_inline)) INLINE static void
rt_tchem_copy_from_grackle1(struct part* p, gr_float rho,
                               grackle_field_data* data) {
  
  const gr_float one_over_rho = 1. / rho;

  p->rt_data.tchem.mass_fraction[rt_species_HI] =
      *data->HI_density * one_over_rho;

  p->rt_data.tchem.mass_fraction[rt_species_HII] =
      *data->HII_density * one_over_rho;

  p->rt_data.tchem.mass_fraction[rt_species_HeI] =
      *data->HeI_density * one_over_rho;
  
  p->rt_data.tchem.mass_fraction[rt_species_HeII] =
      *data->HeII_density * one_over_rho;
  
  p->rt_data.tchem.mass_fraction[rt_species_HeIII] =
      *data->HeIII_density * one_over_rho;
  
  p->rt_data.tchem.mass_fraction[rt_species_e] =
      *data->e_density * one_over_rho;

}


/**
 * @brief copy the grackle data to particle data for grackle chemistry mode 2.
 *
 * @param p particle to use
 * @param rho particle physical density
 * @param particle_grackle_data The grackle_field_data structure from grackle.
 **/
__attribute__((always_inline)) INLINE static void
rt_tchem_copy_from_grackle2(struct part* p, gr_float rho,
                               grackle_field_data* data) {
 
  const gr_float one_over_rho = 1. / rho;
  double rhodust;

  p->rt_data.tchem.mass_fraction[rt_species_HM] =
      *data->HM_density * one_over_rho;

  p->rt_data.tchem.mass_fraction[rt_species_H2I] =
      *data->H2I_density * one_over_rho;

  p->rt_data.tchem.mass_fraction[rt_species_H2II] =
      *data->H2II_density * one_over_rho;

#ifdef SWIFT_RT_GRACKLE_DUST
  /* Load gas metallicities */
      p->chemistry_data.metal_mass_fraction[1] = *data->He_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[2] = *data->C_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[3] = *data->N_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[4] = *data->O_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[5] = *data->Ne_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[6] = *data->Mg_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[7] = *data->Si_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[8] = *data->S_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[9] = *data->Ca_gas_metalDensity * one_over_rho;
      p->chemistry_data.metal_mass_fraction[10] = *data->Fe_gas_metalDensity * one_over_rho;
      /* Load dust metallicities */
      p->rt_data.cooling.dust_mass = *data->dust_density * p->mass * one_over_rho;
      p->rt_data.cooling.dust_mass_fraction[0] = 0.f;
      for (int i=1; i<chemistry_element_count; i++) {
	  if (i==1) rhodust = *data->He_dust_metalDensity;
	  if (i==2) rhodust = *data->C_dust_metalDensity;
	  if (i==3) rhodust = *data->N_dust_metalDensity;
	  if (i==4) rhodust = *data->O_dust_metalDensity;
	  if (i==5) rhodust = *data->Ne_dust_metalDensity;
	  if (i==6) rhodust = *data->Mg_dust_metalDensity;
	  if (i==7) rhodust = *data->Si_dust_metalDensity;
	  if (i==8) rhodust = *data->S_dust_metalDensity;
	  if (i==9) rhodust = *data->Ca_dust_metalDensity;
	  if (i==10) rhodust = *data->Fe_dust_metalDensity;
          p->rt_data.cooling.dust_mass_fraction[i] = rhodust * one_over_rho;
          p->rt_data.cooling.dust_mass_fraction[0] += p->rt_data.cooling.dust_mass_fraction[i];

#endif /* For dust model */
}

}

/**
 * @brief copy the grackle data to particle data.
 *
 * @param p particle to use
 * @param rho particle physical density
 * @param particle_grackle_data The grackle_field_data structure from grackle.
 **/
__attribute__((always_inline)) INLINE static void
rt_tchem_copy_from_grackle(struct part* p, gr_float rho,
                               grackle_field_data* data) {
  
  /* for primordial_chemistry >= 1 */
#if GEARRT_GRACKLE_MODE >= 1
  rt_tchem_copy_from_grackle1(p, rho, data);
#endif 

  /* for primordial_chemistry >= 2 */
#if GEARRT_GRACKLE_MODE >= 2
  rt_tchem_copy_from_grackle2(p, rho, data);
#endif /* For mode 2. */

}

#endif /* SWIFT_RT_GEAR_GRACKLE_COOLING_H */
