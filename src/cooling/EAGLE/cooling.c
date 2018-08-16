/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
/**
 * @file src/cooling/EAGLE/cooling.c
 * @brief EAGLE cooling functions
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "cooling.h"
#include "interpolate.h"
#include "chemistry.h"
#include "cooling_struct.h"
#include "eagle_cool_tables.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

const float explicit_tolerance = 0.05;
const float newton_tolerance = 1.0e-2; // lower than bisection_tol because using log(u)
const float bisection_tolerance = 1.0e-6;
const double bracket_factor = 1.0488088481701; // sqrt(1.1) to match EAGLE

/*
 * @brief calculates heating due to helium reionization
 *
 * @param z redshift
 * @param dz redshift offset
 */
__attribute__((always_inline)) INLINE double eagle_helium_reionization_extraheat(double z, double dz, const struct cooling_function_data *restrict cooling) {

  double extra_heating = 0.0;

/* dz is the change in redshift (start to finish) and hence *should* be < 0 */
#ifdef SWIFT_DEBUG_CHECKS
  if (dz > 0) {
    error(" formulation of helium reionization expects dz<0, whereas you have "
        "dz=%e\n", dz);
  }
#endif

  /* Helium reionization */
  if (cooling->he_reion == 1) {
    double he_reion_erg_pG = cooling->he_reion_ev_pH * eagle_ev_to_erg / eagle_proton_mass_cgs;
    extra_heating += he_reion_erg_pG *
                     (erf((z - dz - cooling->he_reion_z_center) /
                          (pow(2., 0.5) * cooling->he_reion_z_sigma)) -
                      erf((z - cooling->he_reion_z_center) /
                          (pow(2., 0.5) * cooling->he_reion_z_sigma))) /
                     2.0;
  } else
  extra_heating = 0.0;

  return extra_heating;
}

/*
 * @brief Calculates cooling rate by interpolating 4d tables
 * of cooling rates which depend on redshift, temperature (log base 10), 
 * hydrogen number density, helium fraction and metal abundances. 
 * Only the temperature changes when cooling a given particle, so
 * redshift, hydrogen number density and helium fraction indices and
 * offsets passed in.
 *
 * @param log_10_u Log base 10 of internal energy
 * @param dlambda_du Pointer to value to be set to rate
 * of change of cooling rate with respect to internal energy
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param element_lambda Pointer to array for printing contribution 
 * to cooling rate from each of the metals.
 * @param solar_ratio Array of ratios of particle metal abundances
 * to solar metal abundances
 */
__attribute__((always_inline)) INLINE double eagle_metal_cooling_rate(
    double log_10_u, double *dlambda_du, int z_index, float dz, int n_h_i, float d_n_h,
    int He_i, float d_He, const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const,
    double *element_lambda,
    float *solar_ratio) {
  double T_gam, solar_electron_abundance, solar_electron_abundance1, solar_electron_abundance2, elem_cool1, elem_cool2;
  double n_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                            internal_const) *
               cooling->number_density_scale; 
  double z = cosmo->z;
  double cooling_rate = 0.0, temp_lambda, temp_lambda1, temp_lambda2;
  float du;
  double h_plus_he_electron_abundance;
  double h_plus_he_electron_abundance1;
  double h_plus_he_electron_abundance2;

  int i, temp_i;
  double temp;
  float d_temp;

  *dlambda_du = 0.0;


  // interpolate to get temperature of particles, find where we are in
  // the temperature table.

  temp = eagle_convert_u_to_temp(log_10_u, &du, z_index, n_h_i, He_i,
                                 dz, d_n_h, d_He, 
				 cooling, cosmo);
  get_index_1d(cooling->Temp, cooling->N_Temp, temp, &temp_i, &d_temp);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  // contribution to cooling and electron abundance from H, He.
  if (z > cooling->reionisation_redshift) {
    temp_lambda = interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_heating,
                              n_h_i, He_i, temp_i, d_n_h, d_He,
                              d_temp, cooling->N_nH,
                              cooling->N_He, cooling->N_Temp,&temp_lambda2,&temp_lambda1);
    h_plus_he_electron_abundance = interpol_3d(
        cooling->table.photodissociation_cooling.H_plus_He_electron_abundance,
        n_h_i, He_i, temp_i, d_n_h, d_He, d_temp,
        cooling->N_nH, cooling->N_He, cooling->N_Temp,&h_plus_he_electron_abundance2,&h_plus_he_electron_abundance1);
  } else if (z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    temp_lambda = interpol_3d(cooling->table.no_compton_cooling.H_plus_He_heating,
                              n_h_i, He_i, temp_i, d_n_h, d_He,
                              d_temp, cooling->N_nH,
                              cooling->N_He, cooling->N_Temp,&temp_lambda2,&temp_lambda1);
    h_plus_he_electron_abundance = interpol_3d(
        cooling->table.no_compton_cooling.H_plus_He_electron_abundance,
        n_h_i, He_i, temp_i, d_n_h, d_He, d_temp,
        cooling->N_nH, cooling->N_He, cooling->N_Temp,&h_plus_he_electron_abundance2,&h_plus_he_electron_abundance1);
  } else {
    temp_lambda = interpol_4d(cooling->table.element_cooling.H_plus_He_heating,
                              z_index, n_h_i, He_i, temp_i, dz, d_n_h, d_He,
                              d_temp, cooling->N_Redshifts, cooling->N_nH,
                              cooling->N_He, cooling->N_Temp,&temp_lambda2,&temp_lambda1);
    h_plus_he_electron_abundance = interpol_4d(
        cooling->table.element_cooling.H_plus_He_electron_abundance, z_index,
        n_h_i, He_i, temp_i, dz, d_n_h, d_He, d_temp, cooling->N_Redshifts,
        cooling->N_nH, cooling->N_He, cooling->N_Temp,&h_plus_he_electron_abundance2,&h_plus_he_electron_abundance1);
  }
  cooling_rate += temp_lambda;
  *dlambda_du += (temp_lambda2 - temp_lambda1)/du;
  if (element_lambda != NULL) element_lambda[0] = temp_lambda;

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  // inverse Compton cooling is not in collisional table
  // before reionisation so add now 
  
  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {

    T_gam = eagle_cmb_temperature * (1 + z);

    temp_lambda = -eagle_compton_rate *
                  (temp - eagle_cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                  h_plus_he_electron_abundance / n_h;
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  // for each element the cooling rate is multiplied by the ratio of H, He electron
  // abundance to solar electron abundance then by the ratio of the particle metal
  // abundance to solar metal abundance.
  
  if (z > cooling->reionisation_redshift) {
    solar_electron_abundance =
        interpol_2d(cooling->table.photodissociation_cooling.electron_abundance,
                    n_h_i, temp_i, d_n_h, d_temp,
                    cooling->N_nH, cooling->N_Temp,&solar_electron_abundance2,&solar_electron_abundance1);

    for (i = 0; i < cooling->N_Elements; i++) {
      temp_lambda =
          interpol_3d(cooling->table.photodissociation_cooling.metal_heating,
                      n_h_i, temp_i, i, d_n_h, d_temp, 0.0,
                      cooling->N_nH, cooling->N_Temp, cooling->N_Elements,&elem_cool2,&elem_cool1) *
          (h_plus_he_electron_abundance / solar_electron_abundance)*
          solar_ratio[i+2];
      cooling_rate += temp_lambda;
      *dlambda_du +=
          (elem_cool2 * h_plus_he_electron_abundance2 /
               solar_electron_abundance2 -
           elem_cool1 * h_plus_he_electron_abundance1 /
               solar_electron_abundance1) / du * solar_ratio[i+2];
      if (element_lambda != NULL) element_lambda[i + 2] = temp_lambda;
    }
  } else if (z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    solar_electron_abundance =
        interpol_2d(cooling->table.no_compton_cooling.electron_abundance,
                    n_h_i, temp_i, d_n_h, d_temp,
                    cooling->N_nH, cooling->N_Temp,&solar_electron_abundance2,&solar_electron_abundance1);

    for (i = 0; i < cooling->N_Elements; i++) {
      temp_lambda =
          interpol_3d(cooling->table.no_compton_cooling.metal_heating,
                      n_h_i, temp_i, i, d_n_h, d_temp, 0.0,
                      cooling->N_nH, cooling->N_Temp, cooling->N_Elements,&elem_cool2,&elem_cool1) *
          (h_plus_he_electron_abundance / solar_electron_abundance)*
          solar_ratio[i+2];
      cooling_rate += temp_lambda;
      *dlambda_du +=
          (elem_cool2 * h_plus_he_electron_abundance2 /
               solar_electron_abundance2 -
           elem_cool1 * h_plus_he_electron_abundance1 /
               solar_electron_abundance1) / du * solar_ratio[i+2];
      if (element_lambda != NULL) element_lambda[i + 2] = temp_lambda;
    }
  } else {
    solar_electron_abundance =
        interpol_3d(cooling->table.element_cooling.electron_abundance, z_index,
                    n_h_i, temp_i, dz, d_n_h, d_temp, cooling->N_Redshifts,
                    cooling->N_nH, cooling->N_Temp,&solar_electron_abundance2,&solar_electron_abundance1);

    for (i = 0; i < cooling->N_Elements; i++) {
      temp_lambda =
          interpol_4d(cooling->table.element_cooling.metal_heating, z_index,
                      n_h_i, temp_i, i, dz, d_n_h, d_temp, 0.0, cooling->N_Redshifts,
                      cooling->N_nH, cooling->N_Temp, cooling->N_Elements,&elem_cool2,&elem_cool1) *
          (h_plus_he_electron_abundance / solar_electron_abundance)*
          solar_ratio[i+2];
      cooling_rate += temp_lambda;
      *dlambda_du +=
          (elem_cool2 * h_plus_he_electron_abundance2 /
               solar_electron_abundance2 -
           elem_cool1 * h_plus_he_electron_abundance1 /
               solar_electron_abundance1) / du * solar_ratio[i+2];
      if (element_lambda != NULL) element_lambda[i + 2] = temp_lambda;
    }
  }

  return cooling_rate;
}

/*
 * @brief Calculates cooling rate by interpolating precomputed 1d tables
 * of cooling rates which depend only on (log base 10 of) temperature.
 * May be used in bisection scheme when need to perform many table interpolations.
 *
 * @param log_10_u Log base 10 of internal energy
 * @param dlambda_du Pointer to value to be set to rate
 * @param H_plus_He_heat_table 1D table of heating rates due to H, He
 * @param H_plus_He_electron_abundance_table 1D table of electron abundances due
 * to H,He
 * @param element_cooling_table 1D table of heating rates due to metals
 * @param element_electron_abundance_table 1D table of electron abundances due
 * to metals
 * @param temp_table 1D table of temperatures
 * of change of cooling rate with respect to internal energy
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param element_lambda Pointer to array for printing contribution 
 * to cooling rate from each of the metals.
 * @param solar_ratio Array of ratios of particle metal abundances
 * to solar metal abundances
 */
__attribute__((always_inline)) INLINE double
eagle_metal_cooling_rate_1d_table(
    double log_10_u, double *dlambda_du, double *H_plus_He_heat_table,
    double *H_plus_He_electron_abundance_table, double *element_cooling_table,
    double *element_electron_abundance_table, double *temp_table, 
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, 
    double *element_lambda) {
  double T_gam, solar_electron_abundance;
  double n_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                            internal_const) *
               cooling->number_density_scale;  // chemistry_data
  double z = cosmo->z;
  double cooling_rate = 0.0, temp_lambda;
  float du;
  float h_plus_he_electron_abundance;

  int i, temp_i;
  double temp;
  float d_temp;

  *dlambda_du = 0.0;

  // interpolate to get temperature of particles, find where we are in
  // the temperature table.
  temp = eagle_convert_u_to_temp_1d_table(log_10_u, &du, temp_table, cooling);
  get_index_1d(cooling->Temp, cooling->N_Temp, temp, &temp_i, &d_temp);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  // contribution to cooling and electron abundance from H, He.
  temp_lambda = interpol_1d_dbl(H_plus_He_heat_table, temp_i, d_temp);
  h_plus_he_electron_abundance =
      interpol_1d_dbl(H_plus_He_electron_abundance_table, temp_i, d_temp);
  cooling_rate += temp_lambda;
  *dlambda_du += (((double)H_plus_He_heat_table[temp_i + 1]) -
                  ((double)H_plus_He_heat_table[temp_i])) /
                 ((double)du);
  if (element_lambda != NULL) element_lambda[0] = temp_lambda;

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  // inverse Compton cooling is not in collisional table
  // before reionisation so add now 
   
  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {
    /* inverse Compton cooling is not in collisional table
       before reionisation so add now */

    T_gam = eagle_cmb_temperature * (1 + z);

    temp_lambda = -eagle_compton_rate *
                  (temp - eagle_cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                  h_plus_he_electron_abundance / n_h;
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  // for each element the cooling rate is multiplied by the ratio of H, He electron
  // abundance to solar electron abundance. Note: already multiplied by ratio of 
  // particle to solar metal abundances when creating 1d tables.

  solar_electron_abundance =
      interpol_1d_dbl(element_electron_abundance_table, temp_i, d_temp);

    for (i = 0; i < cooling->N_Elements; i++) {
      temp_lambda = interpol_1d_dbl(element_cooling_table + i*cooling->N_Temp, temp_i, d_temp) *
                    (h_plus_he_electron_abundance / solar_electron_abundance); 
      cooling_rate += temp_lambda;
      *dlambda_du +=
          (((double)element_cooling_table[i*cooling->N_Temp + temp_i + 1]) *
               ((double)H_plus_He_electron_abundance_table[temp_i + 1]) /
               ((double)element_electron_abundance_table[temp_i + 1]) -
           ((double)element_cooling_table[i*cooling->N_Temp + temp_i]) *
               ((double)H_plus_He_electron_abundance_table[temp_i]) /
               ((double)element_electron_abundance_table[temp_i])) /
          ((double)du);
      if (element_lambda != NULL) element_lambda[i + 2] = temp_lambda;
    }

  return cooling_rate;
}

/*
 * @brief Wrapper function used to calculate cooling rate and dLambda_du. 
 * Table indices and offsets for redshift, hydrogen number density and 
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param logu Natural log of internal energy
 * @param dlambda_du Pointer to gradient of cooling rate (set in eagle_metal_cooling_rate)
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param abundance_ratio Ratio of element abundance to solar
 */
__attribute__((always_inline)) INLINE double eagle_cooling_rate(
    double logu, double *dLambdaNet_du, int z_index,
    float dz, int n_h_i, float d_n_h, int He_i, float d_He,
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const,
    float *abundance_ratio) {
  
  // set to NULL so will not print file of 
  // contributions to cooling from each element
  double *element_lambda = NULL;   
  double lambda_net = 0.0;

  lambda_net = eagle_metal_cooling_rate(
      logu/eagle_log_10, dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
      He_i, d_He, p, cooling, cosmo, internal_const, element_lambda, abundance_ratio);
  if (isnan(lambda_net)) printf("Eagle cooling.c lambda_net is nan id logu z_index dz %llu %.5e %d %.5e \n", p->id, logu, z_index, dz);
    
  return lambda_net;
}

/*
 * @brief Wrapper function used to calculate cooling rate and dLambda_du. 
 * Table indices and offsets for redshift, hydrogen number density and 
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param logu Natural log of internal energy
 * @param dlambda_du Pointer to gradient of cooling rate (set in eagle_metal_cooling_rate)
 * @param H_plus_He_heat_table 1D table of heating rates due to H, He
 * @param H_plus_He_electron_abundance_table 1D table of electron abundances due
 * to H,He
 * @param element_cooling_table 1D table of heating rates due to metals
 * @param element_electron_abundance_table 1D table of electron abundances due
 * to metals
 * @param temp_table 1D table of temperatures
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param abundance_ratio Ratio of element abundance to solar
 */
__attribute__((always_inline)) INLINE double eagle_cooling_rate_1d_table(
    double logu, double *dLambdaNet_du, double *H_plus_He_heat_table,
    double *H_plus_He_electron_abundance_table, double *element_cooling_table,
    double *element_electron_abundance_table, double *temp_table, int z_index,
    float dz, int n_h_i, float d_n_h, int He_i, float d_He,
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const,
    float *abundance_ratio) {
  double *element_lambda = NULL, lambda_net = 0.0;

  lambda_net = eagle_metal_cooling_rate_1d_table(
      logu/eagle_log_10, dLambdaNet_du, H_plus_He_heat_table,
      H_plus_He_electron_abundance_table, element_cooling_table,
      element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);

  return lambda_net;
}

/*
 * @brief Wrapper function used to calculate cooling rate and dLambda_du. 
 * Writes to file contribution from each element to cooling rate for testing purposes.
 * Table indices and offsets for redshift, hydrogen number density and 
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param abundance_ratio Ratio of element abundance to solar
 */
__attribute__((always_inline)) INLINE double
eagle_print_metal_cooling_rate(
    int z_index, float dz, int n_h_i, float d_n_h, int He_i,
    float d_He, const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const,
    float *abundance_ratio) {
  double *element_lambda, lambda_net = 0.0, dLambdaNet_du;
  element_lambda = malloc((cooling->N_Elements + 2) * sizeof(double));
  double u = hydro_get_physical_internal_energy(p, cosmo) *
             cooling->internal_energy_scale;

  char output_filename[25];
  FILE **output_file = malloc((cooling->N_Elements + 2) * sizeof(FILE *));
  for (int element = 0; element < cooling->N_Elements + 2; element++) {
    sprintf(output_filename, "%s%d%s", "cooling_element_", element, ".dat");
    output_file[element] = fopen(output_filename, "a");
    if (output_file == NULL) {
      printf("Error opening file!\n");
      exit(1);
    }
  }

  for (int j = 0; j < cooling->N_Elements + 2; j++) element_lambda[j] = 0.0;
  lambda_net = eagle_metal_cooling_rate(
      log10(u), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
      He_i, d_He, p, cooling, cosmo, internal_const, element_lambda, abundance_ratio);
  for (int j = 0; j < cooling->N_Elements + 2; j++) {
    fprintf(output_file[j], "%.5e\n", element_lambda[j]);
  }

  for (int i = 0; i < cooling->N_Elements + 2; i++) fclose(output_file[i]);

  return lambda_net;
}

/*
 * @brief Wrapper function used to calculate cooling rate and dLambda_du. 
 * Prints contribution from each element to cooling rate for testing purposes.
 * Table indices and offsets for redshift, hydrogen number density and 
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param H_plus_He_heat_table 1D table of heating rates due to H, He
 * @param H_plus_He_electron_abundance_table 1D table of electron abundances due
 * to H,He
 * @param element_cooling_table 1D table of heating rates due to metals
 * @param element_electron_abundance_table 1D table of electron abundances due
 * to metals
 * @param temp_table 1D table of temperatures
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param abundance_ratio Ratio of element abundance to solar
 */
__attribute__((always_inline)) INLINE double
eagle_print_metal_cooling_rate_1d_table(
    double *H_plus_He_heat_table, double *H_plus_He_electron_abundance_table,
    double *element_print_cooling_table, double *element_electron_abundance_table,
    double *temp_table, const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const) {
  double *element_lambda, lambda_net = 0.0;
  element_lambda = malloc((cooling->N_Elements + 2) * sizeof(double));
  double u = hydro_get_physical_internal_energy(p, cosmo) *
             cooling->internal_energy_scale;
  double dLambdaNet_du;

  char output_filename[25];
  FILE **output_file = malloc((cooling->N_Elements + 2) * sizeof(FILE *));
  for (int element = 0; element < cooling->N_Elements + 2; element++) {
    sprintf(output_filename, "%s%d%s", "cooling_element_", element, ".dat");
    output_file[element] = fopen(output_filename, "a");
    if (output_file == NULL) {
      printf("Error opening file!\n");
      exit(1);
    }
  }

  for (int j = 0; j < cooling->N_Elements + 2; j++) element_lambda[j] = 0.0;
  lambda_net = eagle_metal_cooling_rate_1d_table(
      log10(u), &dLambdaNet_du, H_plus_He_heat_table,
      H_plus_He_electron_abundance_table, element_print_cooling_table,
      element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);
  for (int j = 0; j < cooling->N_Elements + 2; j++) {
    fprintf(output_file[j], "%.5e\n", element_lambda[j]);
  }

  for (int i = 0; i < cooling->N_Elements + 2; i++) fclose(output_file[i]);

  return lambda_net;
}

/*
 * @brief Bisection integration scheme used when Newton Raphson method fails to
 * converge
 *
 * @param logu_init Initial guess for log(internal energy)
 * @param u_ini Internal energy at beginning of hydro step
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cosmo Cosmology structure
 * @param cooling Cooling data structure
 * @param phys_const Physical constants structure
 * @param dt Hydro timestep
 */
__attribute__((always_inline)) INLINE float bisection_iter(
    float logu_init, double u_ini, int z_index,
    float dz, int n_h_i, float d_n_h, int He_i, float d_He, float He_reion_heat,
    struct part *restrict p, const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    float *abundance_ratio, float dt) {
  double u_upper, u_lower, u_next, LambdaNet, dLambdaNet_du;
  int i = 0;
  double u_init = exp(logu_init);

  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction in Hydrogen number density */
  double inn_h =
      chemistry_get_number_density(p, cosmo, chemistry_element_H, phys_const) *
      cooling->number_density_scale;

  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  double ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  
  // Bracketing
  u_lower = u_init;
  u_upper = u_init;
  LambdaNet = (He_reion_heat / (dt * ratefact)) + eagle_cooling_rate(
      log(u_init), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
      He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);

  if (LambdaNet < 0){
    u_lower /= bracket_factor;
    u_upper *= bracket_factor;
    
    LambdaNet = (He_reion_heat / (dt * ratefact)) + eagle_cooling_rate(
        log(u_lower), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
        He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);
    while(u_lower - u_ini - LambdaNet*ratefact*dt > 0) {
      u_lower /= bracket_factor;
      u_upper /= bracket_factor;
      LambdaNet = (He_reion_heat / (dt * ratefact)) + eagle_cooling_rate(
          log(u_lower), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
          He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);
    }
  } else {
    u_lower /= bracket_factor;
    u_upper *= bracket_factor;
    
    LambdaNet = (He_reion_heat / (dt * ratefact)) + eagle_cooling_rate(
        log(u_upper), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
        He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);
    while(u_upper - u_ini - LambdaNet*ratefact*dt < 0) {
      u_lower *= bracket_factor;
      u_upper *= bracket_factor;
      LambdaNet = (He_reion_heat / (dt * ratefact)) + eagle_cooling_rate(
          log(u_upper), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
          He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);
    }
  }

  // bisection iterations
  i = 0;
  do {
    u_next = 0.5 * (u_lower + u_upper);
    LambdaNet = (He_reion_heat / (dt * ratefact)) + eagle_cooling_rate(
        log(u_next), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
        He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);
    if (u_next - u_ini - LambdaNet*ratefact*dt > 0.0) {
      u_upper = u_next;
    } else {
      u_lower = u_next;
    }
    i++;
  } while (fabs(u_upper - u_lower)/u_next > bisection_tolerance && i < 50*eagle_max_iterations);
  
  if (i >= 50*eagle_max_iterations) error("Particle id %llu failed to converge", p->id);

  return log(u_upper);
}


/*
 * @brief Newton Raphson integration scheme to calculate particle cooling over
 * timestep
 *
 * @param logu_init Initial guess for log(internal energy)
 * @param u_ini Internal energy at beginning of hydro step
 * @param z_index Redshift index of particle
 * @param dz Redshift offset of particle
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param He_reion_heat Heating due to helium reionization
 * (only depends on redshift, so passed as parameter)
 * @param p Particle structure
 * @param cosmo Cosmology structure
 * @param cooling Cooling data structure
 * @param phys_const Physical constants structure
 * @param abundance_ratio Array of ratios of metal abundance to solar
 * @param dt timestep
 * @param bisection_flag Flag to identify if scheme failed to converge
 */
__attribute__((always_inline)) INLINE float newton_iter(
    float logu_init, double u_ini, int z_index,
    float dz, int n_h_i, float d_n_h, int He_i, float d_He,
    float He_reion_heat, struct part *restrict p, 
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const, 
    float *abundance_ratio, float dt, int *bisection_flag) {
  /* this routine does the iteration scheme, call one and it iterates to
   * convergence */

  double logu, logu_old;
  double dLambdaNet_du, LambdaNet;
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction in Hydrogen number density */
  double inn_h =
      chemistry_get_number_density(p, cosmo, chemistry_element_H, phys_const) *
      cooling->number_density_scale;

  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  double ratefact = inn_h * (XH / eagle_proton_mass_cgs);

  logu_old = logu_init;
  logu = logu_old;
  int i = 0;
  float log_table_bound_high = (cooling->Therm[cooling->N_Temp-1]+1)/eagle_log_10_e; 
  float log_table_bound_low = (cooling->Therm[0]+1)/eagle_log_10_e; 

  do /* iterate to convergence */
  {
    logu_old = logu;
    LambdaNet = (He_reion_heat / (dt * ratefact)) +
        eagle_cooling_rate(
        logu_old, &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
        He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);

    // Newton iteration.
    logu = logu_old - (1.0 - u_ini * exp(-logu_old) -
	         LambdaNet * ratefact * dt * exp(-logu_old)) /
	    (1.0 - dLambdaNet_du * ratefact * dt);

    // check whether iterations go out of bounds of table,
    // if out of bounds, try again, guess average between old 
    // and table bound
    if (logu > log_table_bound_high) { 
      logu = (log_table_bound_high + logu_old)/2.0;
    } else if (logu < log_table_bound_low) {
      logu = (log_table_bound_low + logu_old)/2.0;
    }

    i++;
  } while (fabs(logu - logu_old) > newton_tolerance && i < eagle_max_iterations);
  if (i >= eagle_max_iterations) {
    // flag to trigger bisection scheme
    *bisection_flag = 1;
  }

  return logu;
}

/*
 * @brief Calculate ratio of particle element abundances
 * to solar abundance
 *
 * @param p Pointer to particle data
 * @param cooling Pointer to cooling data
 * @param ratio_solar Pointer to array or ratios
 */
__attribute__((always_inline)) INLINE void abundance_ratio_to_solar(
		const struct part *restrict p, 
		const struct cooling_function_data *restrict cooling,
		float *ratio_solar) {

  // compute ratios for all elements	
  for(enum chemistry_element elem = chemistry_element_H; elem < chemistry_element_count; elem++){
    if (elem == chemistry_element_Fe) {
      // NOTE: solar abundances have iron last with calcium and sulphur directly before, hence +2
      ratio_solar[elem] = p->chemistry_data.metal_mass_fraction[elem]/
		       cooling->SolarAbundances[elem+2]; 		
    } else {
      ratio_solar[elem] = p->chemistry_data.metal_mass_fraction[elem]/
		       cooling->SolarAbundances[elem];
    }
  }

  // assign ratios for Ca and S, note positions of these elements occur before Fe  
  ratio_solar[chemistry_element_count] = p->chemistry_data.metal_mass_fraction[chemistry_element_Si]*
               cooling->sulphur_over_silicon_ratio/
    	       cooling->SolarAbundances[chemistry_element_count-1];	
  ratio_solar[chemistry_element_count+1] = p->chemistry_data.metal_mass_fraction[chemistry_element_Si]*
               cooling->calcium_over_silicon_ratio/
    	       cooling->SolarAbundances[chemistry_element_count];	
}


/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The time-step of this particle.
 */
void cooling_cool_part(
		const struct phys_const *restrict phys_const,
		const struct unit_system *restrict us,
		const struct cosmology *restrict cosmo,
		const struct cooling_function_data *restrict cooling,
		struct part *restrict p, struct xpart *restrict xp, float dt) {

  float dz;
  int z_index;
  float XH, HeFrac;
  double inn_h;
  double ratefact, u, LambdaNet, LambdaTune = 0, dLambdaNet_du;
  int He_i, n_h_i;
  float d_He, d_n_h;
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);
  
  double u_old = (hydro_get_physical_internal_energy(p, cosmo) 
  		+ hydro_du_dt * dt) *
  	cooling->internal_energy_scale;
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  
  dt *= units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  
  u = u_old;
  double u_ini = u_old;
  
  float abundance_ratio[chemistry_element_count + 2];
  abundance_ratio_to_solar(p, cooling, abundance_ratio);
  
  XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
  	(XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  
  /* convert Hydrogen mass fraction in Hydrogen number density */
  inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H, phys_const) *
  	cooling->number_density_scale;
  
  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);

  /* set helium and hydrogen reheating term */
  if (cooling->he_reion == 1) LambdaTune = eagle_helium_reionization_extraheat(z_index, dz, cooling);
  
  // compute hydrogen number density and helium fraction table indices
  // and offsets (once per particle to save computation)
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  
  // Start cooling, compute cooling rate
  LambdaNet = LambdaTune/(dt*ratefact) + eagle_cooling_rate(
    log(u_ini), &dLambdaNet_du, z_index, dz, n_h_i, d_n_h,
    He_i, d_He, p, cooling, cosmo, phys_const, abundance_ratio);

  if (fabs(ratefact * LambdaNet * dt) < explicit_tolerance * u_old) {
    // if cooling rate is small, take the explicit solution
    u = u_old + ratefact * LambdaNet * dt;
  } else {
    float logu;
  
    double u_eq = u_ini;
  
    if (dt > 0) {
      // newton method
      double log_u_temp = log(u_eq);
      int bisection_flag = 0;
      logu = newton_iter(log_u_temp, u_ini, z_index, dz,
                n_h_i, d_n_h, He_i, d_He, LambdaTune, p, cosmo, cooling, phys_const, 
                abundance_ratio, dt, &bisection_flag);
      if (bisection_flag == 1) {
        // newton method didn't work, so revert to bisection
        logu = bisection_iter(log_u_temp,u_ini,z_index,dz,n_h_i,d_n_h,He_i,d_He,LambdaTune,p,cosmo,cooling,phys_const,abundance_ratio,dt);
      }
      u = exp(logu);
    }
  }

  // calculate du/dt
  float cooling_du_dt = 0.0;
  if (dt > 0) {
    cooling_du_dt = (u - u_ini) / dt / cooling->power_scale / cosmo->a2_inv;
  }

  /* Update the internal energy time derivative */

  hydro_set_physical_internal_energy_dt(p, cosmo, hydro_du_dt + cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy +=
      -hydro_get_mass(p) * cooling_du_dt *
      (dt / units_cgs_conversion_factor(us, UNIT_CONV_TIME));
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE void cooling_write_flavour(
    hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Cooling Model", "EAGLE");
}

/**
 * @brief Computes the cooling time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE float cooling_timestep(
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us, const struct part *restrict p, const struct xpart *restrict xp) {

  /* Remember to update when using an implicit integrator */
  // const float cooling_rate = cooling->cooling_rate;
  // const float internal_energy = hydro_get_comoving_internal_energy(p);
  // return cooling->cooling_tstep_mult * internal_energy / fabsf(cooling_rate);

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE void cooling_first_init_part(
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
__attribute__((always_inline)) INLINE float cooling_get_radiated_energy(
    const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * To do: incorporate into code so that not all tables are read at the beginning
 * @brief Checks the tables that are currently loaded in memory and read
 * new ones if necessary.
 *
 * @param cooling The #cooling_function_data we play with.
 * @param index_z The index along the redshift axis of the tables of the current
 * z.
 */
inline void eagle_check_cooling_tables(struct cooling_function_data* cooling,
                                int index_z) {

  /* Do we already have the right table in memory? */
  if (cooling->low_z_index == index_z) return;

  if (index_z >= cooling->N_Redshifts) {

    error("Missing implementation");
    // Add reading of high-z tables

    /* Record the table indices */
    cooling->low_z_index = index_z;
    cooling->high_z_index = index_z;
  } else {

    /* Record the table indices */
    cooling->low_z_index = index_z;
    cooling->high_z_index = index_z + 1;

    /* Load the damn thing */
    eagle_read_two_tables(cooling);
  }
}

/**
 * To do: incorporate into code so that not all tables are read at the beginning
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * Here we load the cooling tables corresponding to our redshift.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 */
INLINE void cooling_update(const struct phys_const* phys_const,
                                  const struct unit_system* us,
                                  const struct cosmology* cosmo,
                                  struct cooling_function_data* cooling) {
  /* Current redshift */
  const float redshift = cosmo->z;

  /* Get index along the redshift index of the tables */
  int index_z;
  float delta_z_table;
  get_redshift_index(redshift, &index_z, &delta_z_table, cooling);
  cooling->index_z = index_z;
  cooling->delta_z_table = delta_z_table;

  /* Load the current table (if different from what we have) */
  eagle_check_cooling_tables(cooling, index_z);
}


/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const Physical constants data structure 
 * @param cooling Cooling data structure containing properties to initialize
 */
INLINE void cooling_init_backend(
    struct swift_params *parameter_file, const struct unit_system *us,
    const struct phys_const *phys_const,
    struct cooling_function_data *cooling) {

  char fname[200];

  // read some parameters
  parser_get_param_string(parameter_file, "EagleCooling:filename",
                          cooling->cooling_table_path);
  cooling->reionisation_redshift = parser_get_param_float(
      parameter_file, "EagleCooling:reionisation_redshift");
  cooling->calcium_over_silicon_ratio = parser_get_param_float(
      parameter_file, "EAGLEChemistry:CalciumOverSilicon");
  cooling->sulphur_over_silicon_ratio = parser_get_param_float(
      parameter_file, "EAGLEChemistry:SulphurOverSilicon");
  cooling->he_reion = parser_get_param_float(
      parameter_file, "EagleCooling:he_reion");
  cooling->he_reion_z_center = parser_get_param_float(
      parameter_file, "EagleCooling:he_reion_z_center");
  cooling->he_reion_z_sigma = parser_get_param_float(
      parameter_file, "EagleCooling:he_reion_z_sigma");
  cooling->he_reion_ev_pH = parser_get_param_float(
      parameter_file, "EagleCooling:he_reion_ev_pH");

  // read in cooling tables
  GetCoolingRedshifts(cooling);
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  ReadCoolingHeader(fname, cooling);
  cooling->table = eagle_readtable(cooling->cooling_table_path, cooling);

  // allocate array for solar abundances
  cooling->solar_abundances = malloc(cooling->N_Elements * sizeof(double));

  // compute conversion factors
  cooling->internal_energy_scale =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) /
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->number_density_scale =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) /
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->power_scale = units_cgs_conversion_factor(us, UNIT_CONV_POWER) /
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->temperature_scale =
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
INLINE void cooling_print_backend(
    const struct cooling_function_data *cooling) {

  message("Cooling function is 'EAGLE'.");
}

