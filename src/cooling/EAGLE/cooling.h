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
#ifndef SWIFT_COOLING_EAGLE_H
#define SWIFT_COOLING_EAGLE_H

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <hdf5.h>

/* Local includes. */
// #include "cooling_struct.h"
#include "error.h"
#include "hydro.h"
#include "chemistry.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
// #include "physical_constants_cgs.h"
#include "units.h"
#include "eagle_cool_tables.h"

extern const int element_name_length;
extern const float cmb_temperature;
extern const double compton_rate;
extern const bool metal_cooling_on;
extern const int max_iterations;
extern const float proton_mass_cgs;

enum hdf5_allowed_types {
  hdf5_short,
  hdf5_int,
  hdf5_long,
  hdf5_float,
  hdf5_double,
  hdf5_char
};

__attribute__((always_inline)) INLINE int row_major_index_2d(int i, int j,
                                                               int nx, int ny) {
  return ((i % nx) * ny + (j % ny));
}

__attribute__((always_inline)) INLINE int row_major_index_3d(int i, int j,
                                                              int k, int nx, 
							      int ny, int nz) {
  return ((i % nx) * ny * nz + (j % ny) * nz + (k % nz));
}

__attribute__((always_inline)) INLINE int row_major_index_4d(int i, int j,
                                                              int k, int l, 
							      int nx, int ny, 
							      int nz, int nw) {
  return ((i % nx) * ny * nz * nw + (j % ny) * nz * nw + (k % nz) * nw + (l % nw));
}


/*
 * ----------------------------------------------------------------------
 * This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation.  The table is assumed to
 * be evenly spaced.
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE void get_index_1d(float *table, int ntable, double x, int *i, float *dx) {
  float dxm1;
  const float EPS = 1.e-4;

  dxm1 = (float)(ntable - 1) / (table[ntable - 1] - table[0]);

  if ((float)x <= table[0] + EPS) {
    *i = 0;
    *dx = 0;
  } else if ((float)x >= table[ntable - 1] - EPS) {
    *i = ntable - 2;
    *dx = 1;
  } else {
    *i = (int)floor(((float)x - table[0]) * dxm1);
    *dx = ((float)x - table[*i]) * dxm1;
  }
}

__attribute__((always_inline)) INLINE void get_index_1d_therm(float *table, int ntable, double x, int *i, float *dx) {
  float dxm1;
  const float EPS = 1.e-4;

  dxm1 = (float)(ntable - 1) / (table[ntable - 1] - table[0]);

  if ((float)x <= table[0] + EPS) {
    *i = 0;
    *dx = 0;
  } else if ((float)x >= table[ntable - 1] - EPS) {
    *i = ntable - 2;
    *dx = 1;
  } else {
    *i = (int)floor(((float)x - table[0]) * dxm1);
    //printf("i, x, dxm1: %d, %f, %f, %f\n", *i, (float) x, table[0], dxm1);
    *dx = ((float)x - table[*i]) * dxm1;
  }
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_1d(float *table, int i, float dx) {
  float result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE double interpol_1d_dbl(double *table, int i, float dx) {
  double result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_2d(float *table, int i, int j, float dx, float dy, int nx, int ny) {
  float result;
  int index[4];

  index[0] = row_major_index_2d(i,j,nx,ny);
  index[1] = row_major_index_2d(i,j+1,nx,ny);
  index[2] = row_major_index_2d(i+1,j,nx,ny);
  index[3] = row_major_index_2d(i+1,j+1,nx,ny);
  result = (1 - dx) * (1 - dy) * table[index[0]] + (1 - dx) * dy * table[index[1]] +
           dx * (1 - dy) * table[index[2]] + dx * dy * table[index[3]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE double interpol_2d_dbl(double *table, int i, int j, double dx, double dy, int nx, int ny) {
  double result;
  int index[4];

  index[0] = row_major_index_2d(i,j,nx,ny);
  index[1] = row_major_index_2d(i,j+1,nx,ny);
  index[2] = row_major_index_2d(i+1,j,nx,ny);
  index[3] = row_major_index_2d(i+1,j+1,nx,ny);
  result = (1 - dx) * (1 - dy) * table[index[0]] + (1 - dx) * dy * table[index[1]] +
           dx * (1 - dy) * table[index[2]] + dx * dy * table[index[3]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a tri-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_3d(float *table, int i, int j, int k, float dx, float dy,
                  float dz, int nx, int ny, int nz) {
  float result;
  int index[8];

  index[0] = row_major_index_3d(i,j,k,nx,ny,nz);
  index[1] = row_major_index_3d(i,j,k+1,nx,ny,nz);
  index[2] = row_major_index_3d(i,j+1,k,nx,ny,nz);
  index[3] = row_major_index_3d(i,j+1,k+1,nx,ny,nz);
  index[4] = row_major_index_3d(i+1,j,k,nx,ny,nz);
  index[5] = row_major_index_3d(i+1,j,k+1,nx,ny,nz);
  index[6] = row_major_index_3d(i+1,j+1,k,nx,ny,nz);
  index[7] = row_major_index_3d(i+1,j+1,k+1,nx,ny,nz);
  result = (1 - dx) * (1 - dy) * (1 - dz) * table[index[0]] +
           (1 - dx) * (1 - dy) * dz * table[index[1]] +
           (1 - dx) * dy * (1 - dz) * table[index[2]] +
           (1 - dx) * dy * dz * table[index[3]] +
           dx * (1 - dy) * (1 - dz) * table[index[4]] +
           dx * (1 - dy) * dz * table[index[5]] +
           dx * dy * (1 - dz) * table[index[6]] +
           dx * dy * dz * table[index[7]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_4d(float *table, int i, int j, int k, int l, float dx,
                        float dy, float dz, float dw, int nx, int ny, int nz, int nw) {
  float result;
  int index[16];

  index[0]  = row_major_index_4d(i,j,k,l,nx,ny,nz,nw);
  index[1]  = row_major_index_4d(i,j,k,l+1,nx,ny,nz,nw);
  index[2]  = row_major_index_4d(i,j,k+1,l,nx,ny,nz,nw);
  index[3]  = row_major_index_4d(i,j,k+1,l+1,nx,ny,nz,nw);
  index[4]  = row_major_index_4d(i,j+1,k,l,nx,ny,nz,nw);
  index[5]  = row_major_index_4d(i,j+1,k,l+1,nx,ny,nz,nw);
  index[6]  = row_major_index_4d(i,j+1,k+1,l,nx,ny,nz,nw);
  index[7]  = row_major_index_4d(i,j+1,k+1,l+1,nx,ny,nz,nw);
  index[8]  = row_major_index_4d(i+1,j,k,l,nx,ny,nz,nw);
  index[9]  = row_major_index_4d(i+1,j,k,l+1,nx,ny,nz,nw);
  index[10] = row_major_index_4d(i+1,j,k+1,l,nx,ny,nz,nw);
  index[11] = row_major_index_4d(i+1,j,k+1,l+1,nx,ny,nz,nw);
  index[12] = row_major_index_4d(i+1,j+1,k,l,nx,ny,nz,nw);
  index[13] = row_major_index_4d(i+1,j+1,k,l+1,nx,ny,nz,nw);
  index[14] = row_major_index_4d(i+1,j+1,k+1,l,nx,ny,nz,nw);
  index[15] = row_major_index_4d(i+1,j+1,k+1,l+1,nx,ny,nz,nw);

  result = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table[index[0]] +
           (1 - dx) * (1 - dy) * (1 - dz) * dw * table[index[1]] +
           (1 - dx) * (1 - dy) * dz * (1 - dw) * table[index[2]] +
           (1 - dx) * (1 - dy) * dz * dw * table[index[3]] +
           (1 - dx) * dy * (1 - dz) * (1 - dw) * table[index[4]] +
           (1 - dx) * dy * (1 - dz) * dw * table[index[5]] +
           (1 - dx) * dy * dz * (1 - dw) * table[index[6]] +
           (1 - dx) * dy * dz * dw * table[index[7]] +
           dx * (1 - dy) * (1 - dz) * (1 - dw) * table[index[8]] +
           dx * (1 - dy) * (1 - dz) * dw * table[index[9]] +
           dx * (1 - dy) * dz * (1 - dw) * table[index[10]] +
           dx * (1 - dy) * dz * dw * table[index[11]] +
           dx * dy * (1 - dz) * (1 - dw) * table[index[12]] +
           dx * dy * (1 - dz) * dw * table[index[13]] +
           dx * dy * dz * (1 - dw) * table[index[14]] +
           dx * dy * dz * dw * table[index[15]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * Get cooling table redshift index
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE void get_redshift_index(float z, int *z_index, float *dz, const struct cooling_function_data* restrict cooling) {
  int i, iz;

  if (get_redshift_index_first_call == 0) {
    get_redshift_index_first_call = 1;
    get_redshift_index_previous = cooling->N_Redshifts - 2;

    /* this routine assumes cooling_redshifts table is in increasing order. Test
     * this. */
    for (i = 0; i < cooling->N_Redshifts - 2; i++)
      if (cooling->Redshifts[i + 1] < cooling->Redshifts[i]) {
        error("[get_redshift_index]: table should be in increasing order\n");
      }
  }

  /* before the earliest redshift or before hydrogen reionization, flag for
   * collisional cooling */
  //printf("Eagle cooling.h z, n_redshifts-1, max redshift, reionisation_redshift: %.5e, %d, %.5e, %.5e\n",z,cooling->N_Redshifts-1,cooling->Redshifts[cooling->N_Redshifts - 1],cooling->reionisation_redshift);
  if (z > cooling->reionisation_redshift) {
    *z_index = cooling->N_Redshifts;
    *dz = 0.0;
  }
  /* from reionization use the cooling tables */
  else if (z > cooling->Redshifts[cooling->N_Redshifts - 1] &&
           z <= cooling->reionisation_redshift) {
    *z_index = cooling->N_Redshifts + 1;
    *dz = 0.0;
  }
  /* at the end, just use the last value */
  else if (z <= cooling->Redshifts[0]) {
    *z_index = 0;
    *dz = 0.0;
  } else {
    /* start at the previous index and search */
    for (iz = get_redshift_index_previous; iz >= 0; iz--) {
      if (z > cooling->Redshifts[iz]) {
        *dz = (z - cooling->Redshifts[iz]) /
              (cooling->Redshifts[iz + 1] - cooling->Redshifts[iz]);

        get_redshift_index_previous = *z_index = iz;

        break;
      }
    }
  }
}


/*
 * @brief interpolates temperature from internal energy based on table
 *
 */
__attribute__((always_inline)) INLINE static double eagle_convert_u_to_temp(const struct part* restrict p, const struct cooling_function_data* restrict cooling, const struct cosmology* restrict cosmo) {
  float d_z;
  float inn_h = p->chemistry_data.metal_mass_fraction[chemistry_element_H]*cooling->number_density_scale;
  float inHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He]*cooling->number_density_scale; // chemistry data
  float u = hydro_get_physical_internal_energy(p,cosmo)*cooling->internal_energy_scale;
  
  int u_i, He_i, n_h_i;

  float d_n_h, d_He, d_u, logT, T;
  float z = cosmo->z,dz;
  int z_index;

  get_redshift_index(z,&z_index,&dz,cooling);	// ADD OFFSET (d_z) CALCULATION INTO THIS FUNCTION
  get_index_1d(cooling->Redshifts, cooling->N_Redshifts, z, &z_index, &d_z);
  get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  get_index_1d_therm(cooling->Therm, cooling->N_Temp, log10(u), &u_i, &d_u);

  if (z_index == cooling->N_Redshifts+1){
  logT = interpol_4d(cooling->table.photoionisation_cooling.temperature, 0, He_i, n_h_i,
                     u_i, d_z, d_He, d_n_h, d_u, cooling->N_Redshifts,cooling->N_He,cooling->N_nH,cooling->N_Temp);
  } else if (z_index > cooling->N_Redshifts+1){
  logT = interpol_4d(cooling->table.collisional_cooling.temperature, 0, He_i, n_h_i,
                     u_i, d_z, d_He, d_n_h, d_u, cooling->N_Redshifts,cooling->N_He,cooling->N_nH,cooling->N_Temp);
  } else {
  logT = interpol_4d(cooling->table.element_cooling.temperature, 0, He_i, n_h_i,
                     u_i, d_z, d_He, d_n_h, d_u, cooling->N_Redshifts,cooling->N_He,cooling->N_nH,cooling->N_Temp);
  }

  T = pow(10.0, logT);

  if (u_i == 0 && d_u == 0) T *= u / pow(10.0, cooling->Therm[0]);

  return T;
}

/*
 * @brief calculates cooling rate
 *
 */
__attribute__((always_inline)) INLINE static double eagle_cooling_rate(const struct part* restrict p, const struct cooling_function_data* restrict cooling, const struct cosmology* restrict cosmo, const struct phys_const *internal_const) {
  double T_gam, LambdaNet = 0.0, solar_electron_abundance;
  float n_h = chemistry_get_number_density(p,chemistry_element_H,internal_const)*cooling->number_density_scale; // chemistry_data
  double z = cosmo->z;
  float d_z,dz; // ARE THESE ACTUALLY MEANT TO BE THE SAME VALUE???
  int z_index;
  float electron_abundance;
  float *cooling_solar_abundances = malloc(cooling->N_Elements*sizeof(float));

  int i;
  double temp;
  int n_h_i, He_i, temp_i;
  float d_n_h, d_He, d_temp;
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  
  get_redshift_index(z,&z_index,&dz,cooling);	// ADD OFFSET (d_z) CALCULATION INTO THIS FUNCTION
  get_index_1d(cooling->Redshifts, cooling->N_Redshifts, z, &z_index, &d_z);
  
  temp = eagle_convert_u_to_temp(p,cooling,cosmo);

  get_index_1d(cooling->Temp, cooling->N_Temp, log10(temp), &temp_i, &d_temp);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(n_h), &n_h_i, &d_n_h);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  if (z_index == cooling->N_Redshifts+1){
    LambdaNet =
        interpol_4d(cooling->table.photoionisation_cooling.H_plus_He_heating, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
    electron_abundance =
        interpol_4d(cooling->table.photoionisation_cooling.H_plus_He_electron_abundance, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
  } else if (z_index > cooling->N_Redshifts+1){
    LambdaNet =
        interpol_4d(cooling->table.collisional_cooling.H_plus_He_heating, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
    electron_abundance =
        interpol_4d(cooling->table.collisional_cooling.H_plus_He_electron_abundance, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
  } 

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */
  // Is this not handled in the above if statement?

  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {
    /* inverse Compton cooling is not in collisional table
       before reionisation so add now */

    T_gam = cmb_temperature * (1 + z);

    LambdaNet -= compton_rate * (temp - cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                 electron_abundance / n_h;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  if (metal_cooling_on) {
    /* for each element, find the abundance and multiply it
       by the interpolated heating-cooling */

    set_cooling_SolarAbundances(p->chemistry_data.metal_mass_fraction, cooling_solar_abundances, cooling);

    solar_electron_abundance =
        interpol_3d(cooling->SolarElectronAbundance, 0, n_h_i, temp_i, d_z,
                    d_n_h, d_temp,cooling->N_nH,cooling->N_Temp,cooling->N_Redshifts); /* ne/n_h */

    for (i = 0; i < cooling->N_Elements; i++)
      if (cooling_solar_abundances[i] > 0)
        LambdaNet +=
            interpol_4d(cooling->table.element_cooling.metal_heating, 0, i, n_h_i,
                        temp_i, d_z, 0.0, d_n_h, d_temp,cooling->N_Redshifts,cooling->N_Elements,cooling->N_nH,cooling->N_Temp) *
            (electron_abundance / solar_electron_abundance) *
            cooling_solar_abundances[i];
  }

  return LambdaNet;
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
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {
  
  double u_old = hydro_get_comoving_internal_energy(p)*cooling->internal_energy_scale;
  //double rho = hydro_get_comoving_density(p);
  float dz; 
  int z_index;
  get_redshift_index(cosmo->z,&z_index,&dz,cooling);

  float inn_h, XH, HeFrac;

  double du, ratefact, u, u_upper, u_lower, LambdaNet, LambdaTune = 0;
  int i;

  static double zold = 100, LambdaCumul = 0;

  u = u_old;
  u_lower = u;
  u_upper = u;
  //printf("eagle cooling.h internal energy: %.5e\n", u);

  XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] / (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  /* convert Hydrogen mass fraction in Hydrogen number density */
  //inn_h = rho * XH / PROTONMASS;
  inn_h = chemistry_get_number_density(p,chemistry_element_H,phys_const)*cooling->number_density_scale;
  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  ratefact = inn_h * (XH / proton_mass_cgs);
  /* set helium and hydrogen reheating term */
  //LambdaTune = eagle_helium_reionization_extraheat(z_index, dz); // INCLUDE HELIUM REIONIZATION????
  if (zold > z_index) {
    LambdaCumul += LambdaTune;
    printf(" EagleDoCooling %d %g %g %g\n", z_index, dz, LambdaTune, LambdaCumul);
    zold = z_index;
  }

  /* iterative, implicit cooling */
  if (dt > 0)
    {
      LambdaNet = (LambdaTune / (dt * ratefact)) + eagle_cooling_rate(p, cooling, cosmo, phys_const);
    }                                                                                                                  
  else
    {
      LambdaNet = eagle_cooling_rate(p, cooling, cosmo, phys_const);
    }

  if (fabs(ratefact * LambdaNet * dt) < 0.05 * u_old) {
    /* cooling rate is small, take the explicit solution */
    u = u_old + ratefact * LambdaNet * dt;
  }
  else
  {
    i = 0;

    /* bracketing  */
    if (u - u_old - ratefact * LambdaNet * dt < 0) /* heating  */
      {
        u_upper *= sqrt(1.1);
        u_lower /= sqrt(1.1);

        while (u_upper - u_old - LambdaTune - ratefact * eagle_cooling_rate(p, cooling, cosmo, phys_const) * dt < 0
               && i < max_iterations) {
          u_upper *= 1.1;
          u_lower *= 1.1;
          i++;
        }
      }

    if (i == max_iterations) printf("Problem with cooling finding upper bound\n");

    if (u - u_old - ratefact * LambdaNet * dt > 0) /* cooling */
      {
        u_lower /= sqrt(1.1);
        u_upper *= sqrt(1.1);

        i = 0;

        while(u_lower - u_old - LambdaTune - ratefact * eagle_cooling_rate(p, cooling, cosmo, phys_const) * dt > 0
              && i < max_iterations)
          {
            u_upper /= 1.1;
            u_lower /= 1.1;
            i++;
          }
      }

    if (i == max_iterations) printf("Problem with cooling finding lower bound\n");

    i = 0;

    do /* iterate to convergence */
      {
        u = 0.5 * (u_lower + u_upper);

        LambdaNet = (LambdaTune / (dt * ratefact)) + eagle_cooling_rate(p, cooling, cosmo, phys_const);

        if (u - u_old - ratefact * LambdaNet * dt > 0)
          u_upper = u;
        else
          u_lower = u;

        du = u_upper - u_lower;

        i++;

        if (i >= (max_iterations - 10)) printf("u = %g\n", u);
      } while (fabs(du / u) > 1.0e-6 && i < max_iterations);

    if (i >= max_iterations) printf("failed to converge in EagleDoCooling()\n");
  }

  float cooling_du_dt = 0.0;
  if (dt > 0) cooling_du_dt = (u - u_old)/dt/cooling->power_scale;
  const float hydro_du_dt = hydro_get_internal_energy_dt(p);

  /* Integrate cooling equation to enforce energy floor */
  if (u_old + hydro_du_dt * dt < cooling->min_energy) {
    cooling_du_dt = 0.f;
  } else if (u_old + (hydro_du_dt + cooling_du_dt) * dt < cooling->min_energy) {
    cooling_du_dt = (u_old + dt * hydro_du_dt - cooling->min_energy) / dt;
  }
  //printf("Eagle cooling.h new internal energy, old internal energy, dt, cooling_du_dt, cooling_rate: %.5e, %.5e, %.5e, %.5e, %.5e\n", u, u_old, dt, cooling_du_dt, eagle_cooling_rate(p, cooling, cosmo, phys_const));

  /* Update the internal energy time derivative */
  hydro_set_internal_energy_dt(p, hydro_du_dt + cooling_du_dt);
  //hydro_set_internal_energy_dt(p, hydro_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy += -hydro_get_mass(p) * cooling_du_dt * dt;
  //xp->cooling_data.radiated_energy += 0.0;

}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
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
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us, const struct part* restrict p) {

  /* Remember to update when using an implicit integrator!!!*/
  //const float cooling_rate = cooling->cooling_rate;
  //const float internal_energy = hydro_get_comoving_internal_energy(p);
  //return cooling->cooling_tstep_mult * internal_energy / fabsf(cooling_rate);

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct cooling_function_data* cooling) {

  xp->cooling_data.radiated_energy = 0.f; // Why is this zero???
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
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {
  
  char fname[200];

  parser_get_param_string(parameter_file, "EagleCooling:filename",cooling->cooling_table_path);
  GetCoolingRedshifts(cooling);
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  ReadCoolingHeader(fname,cooling);
  cooling->table = eagle_readtable(cooling->cooling_table_path,cooling);
  printf("Eagle cooling.h read table \n");

  cooling->internal_energy_scale = units_cgs_conversion_factor(us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
  cooling->number_density_scale = units_cgs_conversion_factor(us,UNIT_CONV_DENSITY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS)/(2.0*phys_const->const_proton_mass);
  cooling->power_scale = units_cgs_conversion_factor(us,UNIT_CONV_POWER)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
  cooling->temperature_scale = units_cgs_conversion_factor(us,UNIT_CONV_TEMPERATURE);

}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'EAGLE'.");
}

#endif /* SWIFT_COOLING_EAGLE_H */
