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

//static const int eagle_element_name_length = 20;
//static const float eagle_cmb_temperature = 2.728;
//static const double eagle_compton_rate = 1.0178e-37*2.728*2.728*2.728*2.728;
//static const bool eagle_metal_cooling_on = 1;
//static const int eagle_max_iterations = 150;
//static const float eagle_proton_mass_cgs = 1.6726e-24;

static int get_redshift_index_first_call = 0;
static int get_redshift_index_previous = -1;

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
  int index = (i % nx) * ny + (j % ny);
  if (index >= nx*ny) fprintf(stderr, "row_major_index_2d out of bounds, i, j, nx, ny, index, nx*ny %d, %d, %d, %d, %d, %d \n",i,j,nx,ny,index,nx*ny); 
  return index;
}

__attribute__((always_inline)) INLINE int row_major_index_3d(int i, int j,
                                                              int k, int nx, 
							      int ny, int nz) {
  int index = (i % nx) * ny * nz + (j % ny) * nz + (k % nz);
  if (index >= nx*ny*nz) fprintf(stderr, "row_major_index_3d out of bounds, i, j, k, nx, ny, nz, index, nx*ny*nz %d, %d, %d, %d, %d, %d, %d, %d \n",i,j,k,nx,ny,nz,index,nx*ny*nz); 
  return index;
}

__attribute__((always_inline)) INLINE int row_major_index_4d(int i, int j,
                                                              int k, int l, 
							      int nx, int ny, 
							      int nz, int nw) {
  int index = (i % nx) * ny * nz * nw + (j % ny) * nz * nw + (k % nz) * nw + (l % nw);
  if (index >= nx*ny*nz*nw) fprintf(stderr, "row_major_index_4d out of bounds, i, j, k, l, nx, ny, nz, nw, index, nx*ny*nz*nw %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n",i,j,k,l,nx,ny,nz,nw,index,nx*ny*nz*nw); 
  return index;
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
  if(index[0] >= nx*ny || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j %d, %d, table size %d\n", index[0],i,j,nx*ny);
  if(index[1] >= nx*ny || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j %d, %d, table size %d\n", index[1],i,j+1,nx*ny);
  if(index[2] >= nx*ny || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j %d, %d, table size %d\n", index[2],i+1,j,nx*ny);
  if(index[3] >= nx*ny || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j %d, %d, table size %d\n", index[3],i+1,j+1,nx*ny);

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
  if(index[0] >= nx*ny || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j %d, %d, table size %d\n", index[0],i,j,nx*ny);
  if(index[1] >= nx*ny || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j %d, %d, table size %d\n", index[1],i,j+1,nx*ny);
  if(index[2] >= nx*ny || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j %d, %d, table size %d\n", index[2],i+1,j,nx*ny);
  if(index[3] >= nx*ny || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j %d, %d, table size %d\n", index[3],i+1,j+1,nx*ny);

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
  if(index[0] >= nx*ny*nz || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[0],i,j,k,nx*ny*nz);
  if(index[1] >= nx*ny*nz || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[1],i,j,k+1,nx*ny*nz);
  if(index[2] >= nx*ny*nz || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[2],i,j+1,k,nx*ny*nz);
  if(index[3] >= nx*ny*nz || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[3],i,j+1,k+1,nx*ny*nz);
  if(index[4] >= nx*ny*nz || index[4] < 0) fprintf(stderr,"index 4 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[4],i+1,j,k,nx*ny*nz);
  if(index[5] >= nx*ny*nz || index[5] < 0) fprintf(stderr,"index 5 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[5],i+1,j,k+1,nx*ny*nz);
  if(index[6] >= nx*ny*nz || index[6] < 0) fprintf(stderr,"index 6 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[6],i+1,j+1,k,nx*ny*nz);
  if(index[7] >= nx*ny*nz || index[7] < 0) fprintf(stderr,"index 7 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[7],i+1,j+1,k+1,nx*ny*nz);

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
  if(index[0] >= nx*ny*nz*nw || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[0],i,j,k,l,nx*ny*nz*nw);
  if(index[1] >= nx*ny*nz*nw || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[1],i,j,k,l+1,nx*ny*nz*nw);
  if(index[2] >= nx*ny*nz*nw || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[2],i,j,k+1,l,nx*ny*nz*nw);
  if(index[3] >= nx*ny*nz*nw || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[3],i,j,k+1,l+1,nx*ny*nz*nw);
  if(index[4] >= nx*ny*nz*nw || index[4] < 0) fprintf(stderr,"index 4 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[4],i,j+1,k,l,nx*ny*nz*nw);
  if(index[5] >= nx*ny*nz*nw || index[5] < 0) fprintf(stderr,"index 5 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[5],i,j+1,k,l+1,nx*ny*nz*nw);
  if(index[6] >= nx*ny*nz*nw || index[6] < 0) fprintf(stderr,"index 6 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[6],i,j+1,k+1,l,nx*ny*nz*nw);
  if(index[7] >= nx*ny*nz*nw || index[7] < 0) fprintf(stderr,"index 7 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[7],i,j+1,k+1,l+1,nx*ny*nz*nw);
  if(index[8] >= nx*ny*nz*nw || index[8] < 0) fprintf(stderr,"index 8 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[8],i+1,j,k,l,nx*ny*nz*nw);
  if(index[9] >= nx*ny*nz*nw || index[9] < 0) fprintf(stderr,"index 9 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[9],i+1,j,k,l+1,nx*ny*nz*nw);
  if(index[10] >= nx*ny*nz*nw || index[10] < 0) fprintf(stderr,"index 10 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[10],i+1,j,k+1,l,nx*ny*nz*nw);
  if(index[11] >= nx*ny*nz*nw || index[11] < 0) fprintf(stderr,"index 11 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[11],i+1,j,k+1,l+1,nx*ny*nz*nw);
  if(index[12] >= nx*ny*nz*nw || index[12] < 0) fprintf(stderr,"index 12 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[12],i+1,j+1,k,l,nx*ny*nz*nw);
  if(index[13] >= nx*ny*nz*nw || index[13] < 0) fprintf(stderr,"index 13 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[13],i+1,j+1,k,l+1,nx*ny*nz*nw);
  if(index[14] >= nx*ny*nz*nw || index[14] < 0) fprintf(stderr,"index 14 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[14],i+1,j+1,k+1,l,nx*ny*nz*nw);
  if(index[15] >= nx*ny*nz*nw || index[15] < 0) fprintf(stderr,"index 15 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[15],i+1,j+1,k+1,l+1,nx*ny*nz*nw);

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
__attribute__((always_inline)) INLINE static double eagle_convert_u_to_temp(const struct part* restrict p, const struct cooling_function_data* restrict cooling, const struct cosmology* restrict cosmo, const struct phys_const *internal_const) {
  float d_z;
  double inn_h = chemistry_get_number_density(p,chemistry_element_H,internal_const)*cooling->number_density_scale;
  double inHe = chemistry_get_number_density(p,chemistry_element_He,internal_const)*cooling->number_density_scale; // chemistry data
  double u = hydro_get_physical_internal_energy(p,cosmo)*cooling->internal_energy_scale;
  
  int u_i, He_i, n_h_i;

  float d_n_h, d_He, d_u, logT, T;
  float z = cosmo->z,dz;
  int z_index;

  printf("Eagle cooling.h number density non-dim, dim, scale factor %.5e, %.5e, %.5e\n", chemistry_get_number_density(p,chemistry_element_H,internal_const), inn_h, cooling->number_density_scale);

  get_redshift_index(z,&z_index,&dz,cooling);	// ADD OFFSET (d_z) CALCULATION INTO THIS FUNCTION
  get_index_1d(cooling->Redshifts, cooling->N_Redshifts, z, &z_index, &d_z);
  get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  get_index_1d_therm(cooling->Therm, cooling->N_Temp, log10(u), &u_i, &d_u);

  //printf("Eagle cooling.h actual z, HeFrac, nH, internal energy %.5e, %.5e, %.5e, %.5e \n", z, inHe, log10(inn_h), log10(u));
  //printf("Eagle cooling.h interp z, HeFrac, nH, internal energy %.5e, %.5e, %.5e, %.5e \n", cooling->Redshifts[z_index], cooling->HeFrac[He_i], cooling->nH[n_h_i], cooling->Therm[u_i]);

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
  //T = 1.0e6;

  if (u_i == 0 && d_u == 0) T *= u / pow(10.0, cooling->Therm[0]);

  return T;
}

/*
 * @brief calculates cooling rate
 *
 */
__attribute__((always_inline)) INLINE static double eagle_cooling_rate(const struct part* restrict p, const struct cooling_function_data* restrict cooling, const struct cosmology* restrict cosmo, const struct phys_const *internal_const) {
  double T_gam, LambdaNet = 0.0, solar_electron_abundance;
  double n_h = chemistry_get_number_density(p,chemistry_element_H,internal_const)*cooling->number_density_scale; // chemistry_data
  double z = cosmo->z;
  float d_z,dz; // ARE THESE ACTUALLY MEANT TO BE THE SAME VALUE???
  int z_index;
  float h_plus_he_electron_abundance;
  float *cooling_solar_abundances = malloc(cooling->N_Elements*sizeof(float));

  int i;
  double temp;
  int n_h_i, He_i, temp_i;
  float d_n_h, d_He, d_temp;
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  //float HeFrac = 0.248;
  //n_h = pow(10,-4.4);
  //z = 0.0;

  
  get_redshift_index(z,&z_index,&dz,cooling);	// ADD OFFSET (d_z) CALCULATION INTO THIS FUNCTION
  get_index_1d(cooling->Redshifts, cooling->N_Redshifts, z, &z_index, &d_z);
  
  temp = eagle_convert_u_to_temp(p,cooling,cosmo,internal_const);

  get_index_1d(cooling->Temp, cooling->N_Temp, log10(temp), &temp_i, &d_temp);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(n_h), &n_h_i, &d_n_h);
  printf("Eagle cooling.h actual z, HeFrac, log nH, log temperature %.5e, %.5e, %.5e, %.5e \n", z, HeFrac, log10(n_h), log10(temp));
  printf("Eagle cooling.h interp z, HeFrac, log nH, log temperature %.5e, %.5e, %.5e, %.5e \n", cooling->Redshifts[z_index], cooling->HeFrac[He_i], cooling->nH[n_h_i], cooling->Temp[temp_i]);
    

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  if (z_index == cooling->N_Redshifts+1){
    LambdaNet =
        interpol_4d(cooling->table.photoionisation_cooling.H_plus_He_heating, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
    printf("Eagle cooling.h 1 LambdaNet = %.5e\n", LambdaNet);
    h_plus_he_electron_abundance =
        interpol_4d(cooling->table.photoionisation_cooling.H_plus_He_electron_abundance, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
  } else if (z_index > cooling->N_Redshifts+1){
    LambdaNet =
        interpol_4d(cooling->table.collisional_cooling.H_plus_He_heating, 0, He_i, n_h_i,
                     temp_i, 0, d_He, d_n_h, d_temp,1,cooling->N_He,cooling->N_nH,cooling->N_Temp);
    printf("Eagle cooling.h 2 LambdaNet = %.5e\n", LambdaNet);
    h_plus_he_electron_abundance =
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

    T_gam = eagle_cmb_temperature * (1 + z);

    LambdaNet -= eagle_compton_rate * (temp - eagle_cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                 h_plus_he_electron_abundance / n_h; // WATCH OUT WHERE h_plus_he_electron_abundance GETS DEFINED!!!
    printf("Eagle cooling.h 3 LambdaNet = %.5e, %.5e, %.5e, %.5e, %.5e, %.5e, %.5e\n", LambdaNet, eagle_compton_rate, temp, eagle_cmb_temperature*(1+z), pow((1+z),4), h_plus_he_electron_abundance, n_h);
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  if (eagle_metal_cooling_on) {
    /* for each element, find the abundance and multiply it
       by the interpolated heating-cooling */

    set_cooling_SolarAbundances(p->chemistry_data.metal_mass_fraction, cooling_solar_abundances, cooling);

    solar_electron_abundance =
        interpol_3d(cooling->table.element_cooling.electron_abundance, 0, n_h_i, temp_i, d_z,
                    d_n_h, d_temp,cooling->N_Redshifts,cooling->N_nH,cooling->N_Temp); /* ne/n_h */
    //printf("Eagle cooling.h solar electron abundance %.5e, %.5e\n",solar_electron_abundance, cooling->table.element_cooling.electron_abundance[n_h_i*cooling->N_Temp + temp_i]);
    //for(int ii = 0; ii < cooling->N_Redshifts; ii++){
    //  for(int ie = 0; ie < cooling->N_Elements; ie++){
    //    for(int ij = 0; ij < cooling->N_nH; ij++){
    //      for(int ik = 0; ik < cooling->N_Temp; ik++){
    //        double table_metal_heating = cooling->table.element_cooling.metal_heating[ii*cooling->N_Elements*cooling->N_nH*cooling->N_Temp + ie*cooling->N_nH*cooling->N_Temp + ij*cooling->N_Temp + ik];
    //        if (table_metal_heating != 0) printf("Eagle cooling.h cooling->table.element_cooling.metal_heating(z,elem,nh,temp),iz,elem,i_h,i_temp, %.5e, %.5e, %d, %.5e, %.5e\n",table_metal_heating,cooling->Redshifts[ii],ie,cooling->nH[ij],cooling->Temp[ik]);
    //        //printf("Eagle cooling.h cooling->table.element_cooling.metal_heating(iz,i_h,i_temp),%.5e\n",table_metal_heating);
    //      }
    //    }
    //  }
    //}

    for (i = 0; i < cooling->N_Elements; i++){
      if (cooling_solar_abundances[i] > 0 && h_plus_he_electron_abundance != 0 && solar_electron_abundance != 0){
        LambdaNet +=
            interpol_4d(cooling->table.element_cooling.metal_heating, 0, i, n_h_i,
                        temp_i, d_z, 0.0, d_n_h, d_temp,cooling->N_Redshifts,cooling->N_Elements,cooling->N_nH,cooling->N_Temp) *
            (h_plus_he_electron_abundance / solar_electron_abundance) *
            cooling_solar_abundances[i];
        printf("Eagle cooling.h 4 LambdaNet = %.5e\n", LambdaNet);
      } else if (cooling_solar_abundances[i] > 0 && (h_plus_he_electron_abundance == 0 || solar_electron_abundance == 0)){
        LambdaNet +=
            interpol_4d(cooling->table.element_cooling.metal_heating, 0, i, n_h_i,
                        temp_i, d_z, 0.0, d_n_h, d_temp,cooling->N_Redshifts,cooling->N_Elements,cooling->N_nH,cooling->N_Temp) *
            cooling_solar_abundances[i];
        printf("Eagle cooling.h 5 LambdaNet = %.5e\n", LambdaNet);
      }
    }
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

  float XH, HeFrac;
  double inn_h;

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
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);
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
               && i < eagle_max_iterations) {
          u_upper *= 1.1;
          u_lower *= 1.1;
          i++;
        }
      }

    if (i == eagle_max_iterations) printf("Problem with cooling finding upper bound\n");

    if (u - u_old - ratefact * LambdaNet * dt > 0) /* cooling */
      {
        u_lower /= sqrt(1.1);
        u_upper *= sqrt(1.1);

        i = 0;

        while(u_lower - u_old - LambdaTune - ratefact * eagle_cooling_rate(p, cooling, cosmo, phys_const) * dt > 0
              && i < eagle_max_iterations)
          {
            u_upper /= 1.1;
            u_lower /= 1.1;
            i++;
          }
      }

    if (i == eagle_max_iterations) printf("Problem with cooling finding lower bound\n");

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

        if (i >= (eagle_max_iterations - 10)) printf("u = %g\n", u);
      } while (fabs(du / u) > 1.0e-6 && i < eagle_max_iterations);

    if (i >= eagle_max_iterations) printf("failed to converge in EagleDoCooling()\n");
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
  cooling->number_density_scale = units_cgs_conversion_factor(us,UNIT_CONV_DENSITY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);// /(2.0*phys_const->const_proton_mass); // CHECK OUT NUMBER DENSITY SCALING!!!
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
