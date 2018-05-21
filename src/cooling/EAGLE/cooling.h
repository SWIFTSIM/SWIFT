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
#include <time.h>

/* Local includes. */
#include "cooling_struct.h"
#include "error.h"
#include "hydro.h"
#include "chemistry.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"
#include "eagle_cool_tables.h"

/* number of calls to eagle cooling rate */
extern int n_eagle_cooling_rate_calls_1;
extern int n_eagle_cooling_rate_calls_2;
extern int n_eagle_cooling_rate_calls_3;
extern int n_eagle_cooling_rate_calls_4;

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
  int index = i * ny + j;
#ifdef SWIFT_DEBUG_CHECKS
  assert(i < nx);
  assert(j < ny);
#endif
  return index;
}

__attribute__((always_inline)) INLINE int row_major_index_3d(int i, int j,
                                                              int k, int nx, 
							      int ny, int nz) {
  int index = i * ny * nz + j * nz + k;
#ifdef SWIFT_DEBUG_CHECKS
  assert(i < nx);
  assert(j < ny);
  assert(k < nz);
#endif
  return index;
}

__attribute__((always_inline)) INLINE int row_major_index_4d(int i, int j,
                                                              int k, int l, 
							      int nx, int ny, 
							      int nz, int nw) {
  int index = i * ny * nz * nw + j * nz * nw + k * nw + l;
#ifdef SWIFT_DEBUG_CHECKS
  //printf("Eagle cooling.h j, ny %d %d\n",j,ny);
  assert(i < nx);
  assert(j < ny);
  assert(k < nz);
  assert(l < nw);
#endif
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
    if (*i >= ntable || *i < 0){
      printf("Eagle cooling.h i, ntable, x, table[0], dxm1 %d %d %.5e %.5e %.5e \n", *i, ntable, x, table[0], dxm1);
      fflush(stdout);
    }
    *dx = ((float)x - table[*i]) * dxm1;
  }
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
    //printf("Eagle cooling.h z, z_index, z_index_previous, redshifts grid elem, dz %.5e %d %d %.5e %.5e %.5e\n", z,iz,get_redshift_index_previous,*dz, cooling->Redshifts[iz], cooling->Redshifts[iz+1]);
    }
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
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j %d, %d, table size %d\n", index[0],i,j,nx*ny);
  if(index[1] >= nx*ny || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j %d, %d, table size %d\n", index[1],i,j+1,nx*ny);
  if(index[2] >= nx*ny || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j %d, %d, table size %d\n", index[2],i+1,j,nx*ny);
  if(index[3] >= nx*ny || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j %d, %d, table size %d\n", index[3],i+1,j+1,nx*ny);
#endif

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
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j %d, %d, table size %d\n", index[0],i,j,nx*ny);
  if(index[1] >= nx*ny || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j %d, %d, table size %d\n", index[1],i,j+1,nx*ny);
  if(index[2] >= nx*ny || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j %d, %d, table size %d\n", index[2],i+1,j,nx*ny);
  if(index[3] >= nx*ny || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j %d, %d, table size %d\n", index[3],i+1,j+1,nx*ny);
#endif

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
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny*nz || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[0],i,j,k,nx*ny*nz);
  if(index[1] >= nx*ny*nz || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[1],i,j,k+1,nx*ny*nz);
  if(index[2] >= nx*ny*nz || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[2],i,j+1,k,nx*ny*nz);
  if(index[3] >= nx*ny*nz || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[3],i,j+1,k+1,nx*ny*nz);
  if(index[4] >= nx*ny*nz || index[4] < 0) fprintf(stderr,"index 4 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[4],i+1,j,k,nx*ny*nz);
  if(index[5] >= nx*ny*nz || index[5] < 0) fprintf(stderr,"index 5 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[5],i+1,j,k+1,nx*ny*nz);
  if(index[6] >= nx*ny*nz || index[6] < 0) fprintf(stderr,"index 6 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[6],i+1,j+1,k,nx*ny*nz);
  if(index[7] >= nx*ny*nz || index[7] < 0) fprintf(stderr,"index 7 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[7],i+1,j+1,k+1,nx*ny*nz);
#endif

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
#ifdef SWIFT_DEBUG_CHECKS
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
#endif

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

__attribute__((always_inline)) INLINE static void construct_1d_table_from_3d(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,float *result_table){
  int index[4];
  int He_i, n_h_i, z_i;
  float d_He, d_n_h, d_z;
  float inHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;
  
  get_redshift_index(cosmo->z,&z_i,&d_z,cooling);	
  get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  for(int i = 0; i < cooling->N_Temp; i++){
    index[0] = row_major_index_3d(z_i,   n_h_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    index[1] = row_major_index_3d(z_i,   n_h_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    index[2] = row_major_index_3d(z_i+1, n_h_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    index[3] = row_major_index_3d(z_i+1, n_h_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    
    result_table[i] = (1 - d_z) * (1 - d_n_h) * table[index[0]] +
                      (1 - d_z) * d_n_h       * table[index[1]] +
                      d_z       * (1 - d_n_h) * table[index[2]] +
                      d_z       * d_n_h       * table[index[3]];
  }
}

__attribute__((always_inline)) INLINE static void construct_1d_table_from_4d(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,float *result_table){
  int index[8];
  int He_i, n_h_i, z_i;
  float d_He, d_n_h, d_z;
  float inHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;
  
  get_redshift_index(cosmo->z,&z_i,&d_z,cooling);	
  get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  for(int i = 0; i < cooling->N_Temp; i++){
    index[0] = row_major_index_4d(z_i,   n_h_i,   He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[1] = row_major_index_4d(z_i,   n_h_i,   He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[2] = row_major_index_4d(z_i,   n_h_i+1, He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[3] = row_major_index_4d(z_i,   n_h_i+1, He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[4] = row_major_index_4d(z_i+1, n_h_i,   He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[5] = row_major_index_4d(z_i+1, n_h_i,   He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[6] = row_major_index_4d(z_i+1, n_h_i+1, He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[7] = row_major_index_4d(z_i+1, n_h_i+1, He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    
    result_table[i] = (1 - d_z) * (1 - d_n_h) * (1 - d_He) * table[index[0]] +
                      (1 - d_z) * (1 - d_n_h) * d_He       * table[index[1]] +
                      (1 - d_z) * d_n_h       * (1 - d_He) * table[index[2]] +
                      (1 - d_z) * d_n_h       * d_He       * table[index[3]] +
                      d_z       * (1 - d_n_h) * (1 - d_He) * table[index[4]] +
                      d_z       * (1 - d_n_h) * d_He       * table[index[5]] +
                      d_z       * d_n_h       * (1 - d_He) * table[index[6]] +
                      d_z       * d_n_h       * d_He       * table[index[7]];
  }
}

__attribute__((always_inline)) INLINE static void construct_1d_table_from_4d_elements(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,float *result_table){
  int index[4];
  int n_h_i, z_i;
  float d_n_h, d_z;
  float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;
  
  get_redshift_index(cosmo->z,&z_i,&d_z,cooling);	
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  for(int j = 0; j < cooling->N_Elements; j++){
    for(int i = 0; i < cooling->N_Temp; i++){
      index[0] = row_major_index_4d(z_i,   j, n_h_i,   i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      index[1] = row_major_index_4d(z_i,   j, n_h_i+1, i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      index[2] = row_major_index_4d(z_i+1, j, n_h_i,   i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      index[3] = row_major_index_4d(z_i+1, j, n_h_i+1, i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      
      result_table[i+j*cooling->N_Temp] = (1 - d_z) * (1 - d_n_h) * table[index[0]] +
                                          (1 - d_z) * d_n_h       * table[index[1]] +
                                          d_z       * (1 - d_n_h) * table[index[2]] +
                                          d_z       * d_n_h       * table[index[3]];
    }
  }
}

inline int set_cooling_SolarAbundances(const float *element_abundance,
                                double *cooling_element_abundance,
                                const struct cooling_function_data* restrict cooling,
				const struct part* restrict p) {
  int i, index;
  int Silicon_SPH_Index = -1;
  int Calcium_SPH_Index = -1;
  int Sulphur_SPH_Index = -1;
  
  int Silicon_CoolHeat_Index = -1;
  int Calcium_CoolHeat_Index = -1;
  int Sulphur_CoolHeat_Index = -1;

  //float *cooling->ElementAbundance_SOLARM1 = malloc(cooling->N_SolarAbundances*sizeof(float));

    /* determine (inverse of) solar abundance of these elements */
    for (i = 0; i < cooling->N_Elements; i++) {
      index =
          get_element_index(cooling->SolarAbundanceNames,
                            cooling->N_SolarAbundances, cooling->ElementNames[i]);

      if (index < 0) error("Eagle cooling.h index out of bounds");

      index = cooling->SolarAbundanceNamePointers[i];

      if(cooling->SolarAbundances[index] != 0) cooling->ElementAbundance_SOLARM1[i] = 1. / cooling->SolarAbundances[index];
      else cooling->ElementAbundance_SOLARM1[i] = 0.0;

    }

    /* Sulphur tracks Silicon: may choose not to follow Sulphur as SPH element
     */
    /* Same is true for Calcium */
    /* We will assume the code tracks Silicon, and may need to scale Calcium and
     * Sulphur accordingly */

    Silicon_SPH_Index = element_index("Silicon",cooling);
    Calcium_SPH_Index = element_index("Calcium",cooling);
    Sulphur_SPH_Index = element_index("Sulphur",cooling);

    Silicon_CoolHeat_Index =
        get_element_index(cooling->ElementNames, cooling->N_Elements, "Silicon");
    Calcium_CoolHeat_Index =
        get_element_index(cooling->ElementNames, cooling->N_Elements, "Calcium");
    Sulphur_CoolHeat_Index =
        get_element_index(cooling->ElementNames, cooling->N_Elements, "Sulphur");

    if (Silicon_CoolHeat_Index == -1 || Calcium_CoolHeat_Index == -1 ||
        Sulphur_CoolHeat_Index == -1) {
        error("Si, Ca, or S index out of bound\n");
    }

  int sili_index;
  for (i = 0; i < cooling->N_Elements; i++) {
    if (strcmp(chemistry_get_element_name((enum chemistry_element) i), "Silicon") == 0) sili_index = i;
  }

  // Eagle way of identifying and assigning element abundance with strange workaround for calcium and sulphur
  //for (i = 0; i < cooling->N_Elements; i++) {
  //  if (i == Calcium_CoolHeat_Index && Calcium_SPH_Index == -1)
  //    /* SPH does not track Calcium: use Si abundance */
  //    if (Silicon_SPH_Index == -1)
  //      cooling_element_abundance[i] = 0.0;
  //    else{
  //      cooling_element_abundance[i] =
  //          element_abundance[Silicon_SPH_Index] *
  //          cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
  //    }
  //  else if (i == Sulphur_CoolHeat_Index && Sulphur_SPH_Index == -1)
  //    /* SPH does not track Sulphur: use Si abundance */
  //    if (Silicon_SPH_Index == -1)
  //      cooling_element_abundance[i] = 0.0;
  //    else{
  //      cooling_element_abundance[i] =
  //          element_abundance[Silicon_SPH_Index] *
  //          cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
  //    }
  //  else{
  //    cooling_element_abundance[i] = element_abundance[cooling->ElementNamePointers[i]] *
  //                                   cooling_ElementAbundance_SOLARM1[i];
  //    //printf ("Eagle cooling.h element, name, abundance, solar abundance, solarm1, cooling abundance %d %s %.5e %.5e %.5e %.5e\n",cooling->ElementNamePointers[i],cooling->ElementNames[i], element_abundance[cooling->ElementNamePointers[i]],cooling->SolarAbundances[i], cooling_ElementAbundance_SOLARM1[i], cooling_element_abundance[i]);
  //  }
  //}
  
  for (i = 0; i < cooling->N_Elements; i++) {
    if (i == Calcium_CoolHeat_Index && Calcium_SPH_Index != -1)
      if (Silicon_SPH_Index == -1)
        cooling_element_abundance[i] = 0.0;
      else{
        cooling_element_abundance[i] =
            element_abundance[sili_index] * cooling->calcium_over_silicon_ratio *
            cooling->ElementAbundance_SOLARM1[Calcium_CoolHeat_Index];
      }
    else if (i == Sulphur_CoolHeat_Index && Sulphur_SPH_Index != -1)
      /* SPH does not track Sulphur: use Si abundance */
      if (Silicon_SPH_Index == -1)
        cooling_element_abundance[i] = 0.0;
      else{
        cooling_element_abundance[i] =
            element_abundance[sili_index] * cooling->sulphur_over_silicon_ratio *
            cooling->ElementAbundance_SOLARM1[Sulphur_CoolHeat_Index];
      }
    else{
      cooling_element_abundance[i] = element_abundance[cooling->ElementNamePointers[i]] *
                                     cooling->ElementAbundance_SOLARM1[i];
    }
    //printf ("Eagle cooling.h i, solar abundance name pointer, element, name, abundance, solarm1, cooling abundance %d %d %d %s %.5e %.5e %.5e\n", i, cooling->SolarAbundanceNamePointers[i],cooling->ElementNamePointers[i],cooling->ElementNames[i], element_abundance[cooling->ElementNamePointers[i]], cooling->ElementAbundance_SOLARM1[i], cooling_element_abundance[i]);
  }

  return 0;
}

__attribute__((always_inline)) INLINE float eagle_solar_abundance_factor(enum chemistry_element elem,
									 const struct part* restrict p,
									 const struct cooling_function_data* restrict cooling){
  float element_mass_fraction, solar_abundance;

  //if (elem == chemistry_element_S){
  //  element_mass_fraction = p->chemistry_data.metal_mass_fraction[chemistry_element_Si]*cooling->sulphur_over_silicon_ratio;
  //} else if (elem == chemistry_element_Ca){
  //  element_mass_fraction = p->chemistry_data.metal_mass_fraction[chemistry_element_Si]*cooling->calcium_over_silicon_ratio;
  //} else {
    element_mass_fraction = p->chemistry_data.metal_mass_fraction[elem];
    solar_abundance = cooling->SolarAbundances[elem];
  //}
  
  float element_abundance_factor = element_mass_fraction/solar_abundance;

  return element_abundance_factor;
}

/*
 * @brief interpolates internal energy from temperature based on table
 *
 */
__attribute__((always_inline)) INLINE static double eagle_convert_temp_to_u_1d_table(double temp,
										float *temperature_table,
										const struct part* restrict p,
										const struct cooling_function_data* restrict cooling,
										const struct cosmology* restrict cosmo,
										const struct phys_const *internal_const) {
  
  int temp_i;
  float d_temp, u;

  get_index_1d(temperature_table, cooling->N_Temp, log10(temp), &temp_i, &d_temp);

  u = pow(10.0,interpol_1d(cooling->Therm,temp_i,d_temp));

  return u;
}

/*
 * @brief interpolates temperature from internal energy based on table
 *
 */
__attribute__((always_inline)) INLINE static double eagle_convert_u_to_temp_1d_table(double u,
										float *delta_u,
										float *temperature_table,
										const struct part* restrict p,
										const struct cooling_function_data* restrict cooling,
										const struct cosmology* restrict cosmo,
										const struct phys_const *internal_const) {
  
  int u_i;
  float d_u, logT, T;

#ifdef SWIFT_DEBUG_CHECKS
      if (isnan(u)) printf("Eagle cooling.h convert u to temp u is nan");
      if (log10(u) <= cooling->Therm[0]) printf("Eagle cooling.h convert u to temp particle id, u, u_min %llu %.5e %.5e\n", p->id, u, cooling->Therm[0]);
      if (log10(u) >= cooling->Therm[cooling->N_Temp - 1]) printf("Eagle cooling.h convert u to temp particle id, u, u_max %llu %.5e %.5e\n", p->id, u, cooling->Therm[cooling->N_Temp - 1]);
#endif
  get_index_1d(cooling->Therm, cooling->N_Temp, log10(u), &u_i, &d_u);

  logT = interpol_1d(temperature_table,u_i,d_u);
  T = pow(10.0, logT);

  if (u_i == 0 && d_u == 0) T *= u / pow(10.0, cooling->Therm[0]);

  *delta_u = pow(10.0,cooling->Therm[u_i + 1]) - pow(10.0,cooling->Therm[u_i]);

  return T;
}


/*
 * @brief calculates cooling rate
 *
 */
__attribute__((always_inline)) INLINE static double eagle_metal_cooling_rate_1d_table(double u,
										      double *dlambda_du,
									              float *H_plus_He_heat_table,
										      float *H_plus_He_electron_abundance_table,
										      float *element_cooling_table,
										      float *element_electron_abundance_table,
										      float *temp_table,
										      const struct part* restrict p, 
										      const struct cooling_function_data* restrict cooling, 
										      const struct cosmology* restrict cosmo, 
										      const struct phys_const *internal_const, 
										      double* element_lambda) {
  double T_gam, solar_electron_abundance;
  double n_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale; // chemistry_data
  double z = cosmo->z;
  double cooling_rate = 0.0, temp_lambda;
  float dz; 
  float du;
  int z_index;
  float h_plus_he_electron_abundance;

  int i;
  double temp;
  int n_h_i, He_i, temp_i;
  float d_n_h, d_He, d_temp;
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  *dlambda_du = 0.0;
  
  get_redshift_index(z,&z_index,&dz,cooling);	
  
  temp = eagle_convert_u_to_temp_1d_table(u,&du,temp_table,p,cooling,cosmo,internal_const);

  get_index_1d(cooling->Temp, cooling->N_Temp, log10(temp), &temp_i, &d_temp);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(n_h), &n_h_i, &d_n_h);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

    /* Collisional cooling */
    //element_lambda[0] =
    //    interpol_2d(cooling->table.collisional_cooling.H_plus_He_heating, He_i,
    //                 temp_i, d_He, d_temp,cooling->N_He,cooling->N_Temp);
    //h_plus_he_electron_abundance =
    //    interpol_2d(cooling->table.collisional_cooling.H_plus_He_electron_abundance, He_i,
    //                 temp_i, d_He, d_temp,cooling->N_He,cooling->N_Temp);
    
    /* Photodissociation */
    //element_lambda[0] =
    //    interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_heating, He_i,
    //                 temp_i, n_h_i, d_He, d_temp, d_n_h,cooling->N_He,cooling->N_Temp,cooling->N_nH);
    //h_plus_he_electron_abundance =
    //    interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_electron_abundance, He_i,
    //                 temp_i, n_h_i, d_He, d_temp, d_n_h,cooling->N_He,cooling->N_Temp,cooling->N_nH);

    /* redshift tables */
    temp_lambda = interpol_1d(H_plus_He_heat_table, temp_i, d_temp);
    h_plus_he_electron_abundance = interpol_1d(H_plus_He_electron_abundance_table, temp_i, d_temp);
    cooling_rate += temp_lambda;
    *dlambda_du += (H_plus_He_heat_table[temp_i+1] - H_plus_He_heat_table[temp_i])/du;
    if (element_lambda != NULL) element_lambda[0] = temp_lambda;

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {
    /* inverse Compton cooling is not in collisional table
       before reionisation so add now */

    T_gam = eagle_cmb_temperature * (1 + z);

    temp_lambda = -eagle_compton_rate * (temp - eagle_cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                 h_plus_he_electron_abundance / n_h;
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

    /* for each element, find the abundance and multiply it
       by the interpolated heating-cooling */

    //set_cooling_SolarAbundances(p->chemistry_data.metal_mass_fraction, cooling->solar_abundances, cooling, p);

    /* Collisional cooling */
    //solar_electron_abundance =
    //    interpol_1d(cooling->table.collisional_cooling.electron_abundance, temp_i, d_temp); /* ne/n_h */

    //for (i = 0; i < cooling->N_Elements; i++){
    //    element_lambda[i+2] = interpol_2d(cooling->table.collisional_cooling.metal_heating, i,
    //                    temp_i, 0.0, d_temp,cooling->N_Elements,cooling->N_Temp) *
    //        (h_plus_he_electron_abundance / solar_electron_abundance) *
    //        cooling->solar_abundances[i];
    //}
    
    /* Photodissociation */
    //solar_electron_abundance =
    //    interpol_2d(cooling->table.photodissociation_cooling.electron_abundance, temp_i, n_h_i, d_temp, d_n_h, cooling->N_Temp, cooling->N_nH); /* ne/n_h */
      
    //for (i = 0; i < cooling->N_Elements; i++){
    //    element_lambda[i+2] = interpol_3d(cooling->table.photodissociation_cooling.metal_heating, i,
    //                    temp_i, n_h_i, 0.0, d_temp, d_n_h,cooling->N_Elements,cooling->N_Temp,cooling->N_nH) *
    //        (h_plus_he_electron_abundance / solar_electron_abundance) *
    //        cooling->solar_abundances[i];
    //}
    
    /* redshift tables */
    solar_electron_abundance = interpol_1d(element_electron_abundance_table, temp_i, d_temp);
    
    for (i = 0; i < cooling->N_Elements; i++){
	//temp_lambda = interpol_1d(element_cooling_table + i*cooling->N_Temp, temp_i, d_temp) *
        //    (h_plus_he_electron_abundance / solar_electron_abundance) *
        //    cooling->solar_abundances[i];
	temp_lambda = interpol_1d(element_cooling_table + i*cooling->N_Temp, temp_i, d_temp) *
            (h_plus_he_electron_abundance / solar_electron_abundance);
        cooling_rate += temp_lambda;
	*dlambda_du += (element_cooling_table[temp_i+1]*H_plus_He_electron_abundance_table[temp_i+1]/element_electron_abundance_table[temp_i+1] - 
			element_cooling_table[temp_i]*H_plus_He_electron_abundance_table[temp_i]/element_electron_abundance_table[temp_i])/du;
        if (element_lambda != NULL) element_lambda[i+2] = temp_lambda;
    }

    //for(enum chemistry_element k = chemistry_element_C; k < chemistry_element_count; k++){
    //    printf("Eagle cooling.h solar abundance factor %.5e %.5e\n", cooling->solar_abundances[k-2],eagle_solar_abundance_factor(k, p, cooling));
    //}

    return cooling_rate;
}

__attribute__((always_inline)) INLINE static double eagle_cooling_rate_1d_table(double u,
										double *dLambdaNet_du,
									        float *H_plus_He_heat_table,
										float *H_plus_He_electron_abundance_table,
										float *element_cooling_table,
										float *element_electron_abundance_table,
										float *temp_table,
										const struct part* restrict p, 
										const struct cooling_function_data* restrict cooling, 
										const struct cosmology* restrict cosmo, 
										const struct phys_const *internal_const) {
  double *element_lambda = NULL, lambda_net1 = 0.0, lambda_net2 = 0.0, delta;//, dLambdaNet_du_calc;
  float d_u;
  int u_i;
  get_index_1d(cooling->Therm, cooling->N_Temp, log10(u), &u_i, &d_u);
  if (u_i >= cooling->N_Temp-2){
    delta = pow(10.0,cooling->Therm[cooling->N_Temp-2]) - u;
  } else {
    delta = pow(10.0,cooling->Therm[u_i+1]) - pow(10.0,cooling->Therm[u_i]);
  }
  delta *= 0.1;

  lambda_net1 = eagle_metal_cooling_rate_1d_table(u, dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);
  lambda_net2 = eagle_metal_cooling_rate_1d_table(u + delta, dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);
  
  *dLambdaNet_du = (lambda_net2 - lambda_net1)/delta;
  
  //lambda_net1 = eagle_metal_cooling_rate_1d_table(u, dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);


  return lambda_net1;
}

__attribute__((always_inline)) INLINE static double eagle_print_metal_cooling_rate_1d_table(float *H_plus_He_heat_table,
											    float *H_plus_He_electron_abundance_table,
											    float *element_cooling_table,
											    float *element_electron_abundance_table,
											    float *temp_table,
											    const struct part* restrict p, 
											    const struct cooling_function_data* restrict cooling, 
											    const struct cosmology* restrict cosmo, 
											    const struct phys_const *internal_const) {
  double *element_lambda, lambda_net = 0.0;
  element_lambda = malloc((cooling->N_Elements+2)*sizeof(double));
  double u = hydro_get_physical_internal_energy(p,cosmo)*cooling->internal_energy_scale;
  double dLambdaNet_du;
  
  char output_filename[21];
  FILE** output_file = malloc((cooling->N_Elements+2)*sizeof(FILE*));
  for (int element = 0; element < cooling->N_Elements+2; element++){
    sprintf(output_filename, "%s%d%s", "cooling_output_",element,".dat");
    output_file[element] = fopen(output_filename, "a");
    if (output_file == NULL)
    {   
        printf("Error opening file!\n");
        exit(1);
    }
  }

  for (int j = 0; j < cooling->N_Elements+2; j++) element_lambda[j] = 0.0;
  lambda_net = eagle_metal_cooling_rate_1d_table(u,&dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);
  for (int j = 0; j < cooling->N_Elements+2; j++) {
    fprintf(output_file[j],"%.5e\n",element_lambda[j]);
  }
  
  for (int i = 0; i < cooling->N_Elements+2; i++) fclose(output_file[i]);

  return lambda_net;
}


__attribute__((always_inline)) INLINE static double newton_guess_iter(double u,
								      float *H_plus_He_heat_table,
                                                                      float *H_plus_He_electron_abundance_table,
                                                                      float *element_cooling_table,
                                                                      float *element_electron_abundance_table,
                                                                      float *temp_table,
								      const struct part* restrict p,
								      const struct cooling_function_data* restrict cooling,
								      const struct cosmology* restrict cosmo,
    							              const struct phys_const* restrict phys_const,
								      float dt){
  
  double LambdaNet, ratefact, inn_h, u_upper, u_lower, u_old, dLambdaNet_du;
  float XH, HeFrac;
  int i = 0;
  
  XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] / (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,phys_const)*cooling->number_density_scale;
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  //LambdaTune = eagle_helium_reionization_extraheat(z_index, dz); // INCLUDE HELIUM REIONIZATION????
  //if (zold > z_index) {
  //  LambdaCumul += LambdaTune;
  //  zold = z_index;
  //}

  u_upper = u;
  u_lower = u;
  u_old = u;

  LambdaNet = eagle_cooling_rate_1d_table(u, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const);

    /* bracketing  */
    if (LambdaNet < 0) /* heating  */
      {
        u_upper *= sqrt(1.1);
        u_lower /= sqrt(1.1);

        while (u_upper - u_old - ratefact * eagle_cooling_rate_1d_table(u_upper, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const) * dt < 0 && i < eagle_max_iterations) {
          u_upper *= 1.1;
          u_lower *= 1.1;
          i++;
          //n_eagle_cooling_rate_calls_2++;
        }
      }

    //if (i == eagle_max_iterations) printf("Problem with cooling finding upper bound\n");
#ifdef SWIFT_DEBUG_CHECKS
    //if (i < eagle_max_iterations) printf("Eagle cooling.h heating bound converged number of iterations, u_upper, u_lower, u, u_old, cooling %d %.8e %.8e %.8e %.8e %.8e\n", i, u_upper, u_lower, u, u_old,ratefact*eagle_cooling_rate(u_upper, p,cooling,cosmo,phys_const)*dt);
    if (i == eagle_max_iterations){
      printf("Problem with cooling finding upper bound, u_upper, u_lower, u, u_old, cooling %.5e %.5e %.5e %.5e %.5e\n", u_upper, u_lower, u, u_old,ratefact*eagle_cooling_rate(u_upper, p,cooling,cosmo,phys_const)*dt);
    }
#endif

    if (LambdaNet > 0) /* cooling */
      {
        u_lower /= sqrt(1.1);
        u_upper *= sqrt(1.1);

        i = 0;

        while(u_lower - u_old - ratefact * eagle_cooling_rate_1d_table(u_lower, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const) * dt > 0 && i < eagle_max_iterations)
          {
            u_upper /= 1.1;
            u_lower /= 1.1;
            i++;
            //n_eagle_cooling_rate_calls_3++;
          }
      }

    return log(0.5*(u_upper + u_lower));
}


__attribute__((always_inline)) INLINE static float newton_iter(float x_init,
							       double u_ini,
							       float *H_plus_He_heat_table,
							       float *H_plus_He_electron_abundance_table,
							       float *element_cooling_table,
							       float *element_electron_abundance_table,
							       float *temp_table,
    							       struct part* restrict p,
    							       const struct cosmology* restrict cosmo,
    							       const struct cooling_function_data* restrict cooling,
    							       const struct phys_const* restrict phys_const,
							       float dt){
  double x, x_old;
  double dLambdaNet_du, LambdaNet;
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  //float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] / (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  /* convert Hydrogen mass fraction in Hydrogen number density */
  double inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,phys_const)*cooling->number_density_scale;
  //inn_h = hydro_get_physical_density(p,cosmo)*units_cgs_conversion_factor(us,UNIT_CONV_DENSITY) * XH / eagle_proton_mass_cgs;
  
  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  double ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  
  /* set helium and hydrogen reheating term */
  //LambdaTune = eagle_helium_reionization_extraheat(z_index, dz); // INCLUDE HELIUM REIONIZATION????
  //if (zold > z_index) {
  //  LambdaCumul += LambdaTune;
  //  printf(" EagleDoCooling %d %g %g %g\n", z_index, dz, LambdaTune, LambdaCumul);
  //  zold = z_index;
  //}
  
  x_old = x_init ;
  x = x_old;
  //float du;
  int i = 0;
  
  do /* iterate to convergence */
    {
      x_old = x;

      //LambdaNet = (LambdaTune / (dt * ratefact)) + eagle_cooling_rate_1d_table(exp(x_old), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const);
      LambdaNet = eagle_cooling_rate_1d_table(exp(x_old), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const);
      n_eagle_cooling_rate_calls_1++;
      
      x = x_old - (1.0 - u_ini*exp(-x_old) - LambdaNet*ratefact*dt*exp(-x_old))/(1.0 - dLambdaNet_du*ratefact*dt);
      //if (dt > 0 && i >=10) printf("Eagle cooling.h particle id, log u, temperature, lambda_net, terms, error %llu %.5e %.5e %.5e %.5e %.5e %.5e %.5e %d\n", p->id, x, eagle_convert_u_to_temp_1d_table(exp(x),&du,temp_table,p,cooling,cosmo,phys_const), LambdaNet, (1.0 - u_ini*exp(-x_old) - LambdaNet*ratefact*dt*exp(-x_old)),(1.0 - dLambdaNet_du*ratefact*dt),(1.0 - u_ini*exp(-x_old) - LambdaNet*ratefact*dt*exp(-x_old))/(1.0 - dLambdaNet_du*ratefact*dt),fabs((exp(x) - exp(x_old)) / exp(x)), i);
      
      if (i == 10){
        x = newton_guess_iter(u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,p,cooling,cosmo,phys_const,dt);
        n_eagle_cooling_rate_calls_4++;
      }
      if (x > cooling->Therm[cooling->N_Temp - 1]*log(10) + 1) x = newton_guess_iter(u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,p,cooling,cosmo,phys_const,dt); 
      if (x < cooling->Therm[0]*log(10) - 1) x = newton_guess_iter(u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,p,cooling,cosmo,phys_const,dt); 

      if(dt > 0) printf("Eagle cooling.h u, u_old, error, step %.5e %.5e %.5e %d\n", exp(x), exp(x_old), fabs((x - x_old) / x), i);
      //if(dt > 0) printf("Eagle cooling.h %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %d\n", 1.0,  u_ini, exp(-x_old),  LambdaNet, ratefact, dt, exp(-x_old),1.0, dLambdaNet_du, ratefact, dt, i);

      i++;
    } while (fabs((x - x_old) / x) > 1.0e-5 && i < eagle_max_iterations);

  return x;

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
  
  double u_old = hydro_get_physical_internal_energy(p,cosmo)*cooling->internal_energy_scale;
  float dz; 
  int z_index;
  get_redshift_index(cosmo->z,&z_index,&dz,cooling);

  float XH, HeFrac;
  double inn_h;

  double ratefact, u, LambdaNet, LambdaTune = 0, dLambdaNet_du, LambdaNext;

  static double zold = 100, LambdaCumul = 0;
  dt *= units_cgs_conversion_factor(us,UNIT_CONV_TIME);

  u = u_old;
  double u_ini = u_old, u_temp;

  XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] / (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  //printf("Eagle cooling.h density %.5e\n", hydro_get_physical_density(p,cosmo)*units_cgs_conversion_factor(us,UNIT_CONV_DENSITY));

  /* convert Hydrogen mass fraction in Hydrogen number density */
  inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,phys_const)*cooling->number_density_scale;
  //inn_h = hydro_get_physical_density(p,cosmo)*units_cgs_conversion_factor(us,UNIT_CONV_DENSITY) * XH / eagle_proton_mass_cgs;
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

  // construct 1d table of cooling rates wrt temperature
  float H_plus_He_heat_table[176]; 			// WARNING sort out how it is declared/allocated
  float H_plus_He_electron_abundance_table[176]; 	// WARNING sort out how it is declared/allocated
  float temp_table[176]; 				// WARNING sort out how it is declared/allocated
  float element_cooling_table[9*176]; 			// WARNING sort out how it is declared/allocated
  float element_electron_abundance_table[176]; 		// WARNING sort out how it is declared/allocated
  construct_1d_table_from_4d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.H_plus_He_heating,H_plus_He_heat_table);
  construct_1d_table_from_4d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.H_plus_He_electron_abundance,H_plus_He_electron_abundance_table);
  construct_1d_table_from_4d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.temperature,temp_table);
  construct_1d_table_from_4d_elements(p,cooling,cosmo,phys_const,cooling->table.element_cooling.metal_heating,element_cooling_table);
  construct_1d_table_from_3d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.electron_abundance,element_electron_abundance_table);

  n_eagle_cooling_rate_calls_2++;

  LambdaNet = eagle_cooling_rate_1d_table(u_ini, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const);
  //if (fabs(ratefact * LambdaNet * dt) < 0.05 * u_old) {
  //  /* cooling rate is small, take the explicit solution */
  //  u = u_old + ratefact * LambdaNet * dt;
  //}

  u_temp = u_ini + LambdaNet*ratefact*dt;
  if (u_temp > 0) LambdaNext = eagle_cooling_rate_1d_table(u_temp, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const);
  if (fabs(LambdaNet - LambdaNext)/LambdaNet < 0.5) {
    u_temp = u_ini;
  } else {
    u_temp = eagle_convert_temp_to_u_1d_table(1.0e4,temp_table,p,cooling,cosmo,phys_const);
  }

  float x = newton_iter(log(u_temp),u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,p,cosmo,cooling,phys_const,dt);
  u = exp(x);
  //if (i >= eagle_max_iterations) n_eagle_cooling_rate_calls_3++;

  float cooling_du_dt = 0.0;
  if (dt > 0){ 
    cooling_du_dt = (u - u_ini)/dt/cooling->power_scale;
  }

  const float hydro_du_dt = hydro_get_internal_energy_dt(p);

  /* Update the internal energy time derivative */
  hydro_set_internal_energy_dt(p, hydro_du_dt + cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy += -hydro_get_mass(p) * cooling_du_dt * (dt/units_cgs_conversion_factor(us,UNIT_CONV_TIME));

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
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {
  
  char fname[200];

  parser_get_param_string(parameter_file, "EagleCooling:filename",cooling->cooling_table_path);
  cooling->reionisation_redshift = parser_get_param_float(parameter_file, "EagleCooling:reionisation_redshift");
  cooling->calcium_over_silicon_ratio = parser_get_param_float(parameter_file, "EAGLEChemistry:CalciumOverSilicon");
  cooling->sulphur_over_silicon_ratio = parser_get_param_float(parameter_file, "EAGLEChemistry:SulphurOverSilicon");
  GetCoolingRedshifts(cooling);
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  ReadCoolingHeader(fname,cooling);
  MakeNamePointers(cooling);
  cooling->table = eagle_readtable(cooling->cooling_table_path,cooling);
  printf("Eagle cooling.h read table \n");

  cooling->ElementAbundance_SOLARM1 = malloc(cooling->N_SolarAbundances*sizeof(float));
  cooling->solar_abundances = malloc(cooling->N_Elements*sizeof(double));

  cooling->delta_u = cooling->Therm[1] - cooling->Therm[0];

  cooling->internal_energy_scale = units_cgs_conversion_factor(us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
  cooling->number_density_scale = units_cgs_conversion_factor(us,UNIT_CONV_DENSITY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
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
