/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (jame.s.willis@durham.ac.uk)
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
#ifndef SWIFT_CACHE_H
#define SWIFT_CACHE_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "cell.h"
#include "error.h"
#include "part.h"
#include "vector.h"
#include "sort.h"

#define MAX_NO_OF_PARTS 1000
#define NUM_VEC_PROC 2
#define CACHE_ALIGN sizeof(float) * VEC_SIZE
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)
#define C2_CACHE_ALIGN sizeof(float) * VEC_SIZE

/* Cache struct to hold a local copy of a cells' particle
 * properties required for density/force calculations.*/
struct cache {

#ifdef DOPAIR1_AUTO_VEC
  float x[MAX_NO_OF_PARTS]      __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* x position*/
  float y[MAX_NO_OF_PARTS]      __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* y position*/
  float z[MAX_NO_OF_PARTS]      __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* z position*/
  float m[MAX_NO_OF_PARTS]      __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Mass */
  float vx[MAX_NO_OF_PARTS]    __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* x velocity */
  float vy[MAX_NO_OF_PARTS]    __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* y velocity */
  float vz[MAX_NO_OF_PARTS]    __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* z velocity */ 
  float h[MAX_NO_OF_PARTS]      __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Smoothing length */

  /*Cached arrays to hold particle updates*/
  float rho[MAX_NO_OF_PARTS]        __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Density */
  float rho_dh[MAX_NO_OF_PARTS]     __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Density gradient */
  float wcount[MAX_NO_OF_PARTS]     __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* No. of contributions*/
  float wcount_dh[MAX_NO_OF_PARTS]  __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Mass */
  float div_v[MAX_NO_OF_PARTS]      __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Velocity divergence */
  float curl_vx[MAX_NO_OF_PARTS]    __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Velocity curl x */
  float curl_vy[MAX_NO_OF_PARTS]    __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Velocity curl y */ 
  float curl_vz[MAX_NO_OF_PARTS]    __attribute__((aligned(sizeof(float) * VEC_SIZE))); /* Velocity curl z */

  int count;
#else
  /* Particle x position. */
  float *restrict x __attribute__((aligned(CACHE_ALIGN)));

  /* Particle y position. */
  float *restrict y __attribute__((aligned(CACHE_ALIGN)));

  /* Particle z position. */
  float *restrict z __attribute__((aligned(CACHE_ALIGN)));

  /* Particle smoothing length. */
  float *restrict h __attribute__((aligned(CACHE_ALIGN)));

  /* Particle mass. */
  float *restrict m __attribute__((aligned(CACHE_ALIGN)));

  /* Particle x velocity. */
  float *restrict vx __attribute__((aligned(CACHE_ALIGN)));

  /* Particle y velocity. */
  float *restrict vy __attribute__((aligned(CACHE_ALIGN)));

  /* Particle z velocity. */
  float *restrict vz __attribute__((aligned(CACHE_ALIGN)));

  /* Cache size. */
  int count;

#endif
  /* Particle x position. */
  //float *restrict rho __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle y position. */
  //float *restrict rho_dh __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle z position. */
  //float *restrict wcount __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle smoothing length. */
  //float *restrict wcount_dh __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle mass. */
  //float *restrict div_v __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle x velocity. */
  //float *restrict curl_vx __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle y velocity. */
  //float *restrict curl_vy __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  ///* Particle z velocity. */
  //float *restrict curl_vz __attribute__((aligned(sizeof(float) * VEC_SIZE)));

};

#ifdef DOPAIR1_AUTO_VEC
struct cache ci_cache, cj_cache;
#else
struct cache cj_cache;
#endif

/* Secondary cache struct to hold a list of interactions between two
 * particles.*/
struct c2_cache {

  /* Separation between two particles squared. */
  float r2q[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* x separation between two particles. */
  float dxq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* y separation between two particles. */
  float dyq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* z separation between two particles. */
  float dzq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* Mass of particle pj. */
  float mq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* x velocity of particle pj. */
  float vxq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* y velocity of particle pj. */
  float vyq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));

  /* z velocity of particle pj. */
  float vzq[C2_CACHE_SIZE] __attribute__((aligned(C2_CACHE_ALIGN)));
};

/**
 * @brief Allocate memory and initialise cache.
 *
 * @param c The cache.
 * @param count Number of particles to allocate space for.
 */
__attribute__((always_inline)) INLINE void cache_init(struct cache *c,
                                                      size_t count) {

  /* Align cache on correct byte boundary and pad cache size to include 2 vector
   * lengths for remainder operations. */
  unsigned long alignment = sizeof(float) * VEC_SIZE;
  unsigned int sizeBytes = (count + (2 * VEC_SIZE)) * sizeof(float);
  int error = 0;

  /* Free memory if cache has already been allocated. */
  if (c->count > 0) {
    free(c->x);
    free(c->y);
    free(c->z);
    free(c->m);
    free(c->vx);
    free(c->vy);
    free(c->vz);
    free(c->h);
    //free(c->rho);
    //free(c->rho_dh);
    //free(c->wcount);
    //free(c->wcount_dh);
    //free(c->div_v);
    //free(c->curl_vx);
    //free(c->curl_vy);
    //free(c->curl_vz);
  }

  error += posix_memalign((void **)&c->x, alignment, sizeBytes);
  error += posix_memalign((void **)&c->y, alignment, sizeBytes);
  error += posix_memalign((void **)&c->z, alignment, sizeBytes);
  error += posix_memalign((void **)&c->m, alignment, sizeBytes);
  error += posix_memalign((void **)&c->vx, alignment, sizeBytes);
  error += posix_memalign((void **)&c->vy, alignment, sizeBytes);
  error += posix_memalign((void **)&c->vz, alignment, sizeBytes);
  error += posix_memalign((void **)&c->h, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->rho, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->rho_dh, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->wcount, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->wcount_dh, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->div_v, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->curl_vx, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->curl_vy, alignment, sizeBytes);
  //error += posix_memalign((void **)&c->curl_vz, alignment, sizeBytes);

  if (error != 0)
    error("Couldn't allocate cache, no. of particles: %d", (int)count);
  c->count = count;
}

/**
 * @brief Populate cache by reading in the particles in unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 */
__attribute__((always_inline)) INLINE void cache_read_particles(
    const struct cell *const ci, struct cache *const ci_cache) {

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < ci->count; i++) {
    ci_cache->x[i] = ci->parts[i].x[0] - ci->loc[0];
    ci_cache->y[i] = ci->parts[i].x[1] - ci->loc[1];
    ci_cache->z[i] = ci->parts[i].x[2] - ci->loc[2];
    ci_cache->h[i] = ci->parts[i].h;

    ci_cache->m[i] = ci->parts[i].mass;
    ci_cache->vx[i] = ci->parts[i].v[0];
    ci_cache->vy[i] = ci->parts[i].v[1];
    ci_cache->vz[i] = ci->parts[i].v[2];
  }
}

/**
 * @brief Populate cache by reading in the particles from two cells in unsorted order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The cache for cell ci.
 * @param cj_cache The cache for cell cj.
 * @param shift The amount to shift the particle positions to account for BCs
 */
__attribute__((always_inline)) INLINE void cache_read_two_cells(
    const struct cell *const ci, const struct cell *const cj, struct cache *const ci_cache, struct cache *const cj_cache, const double *const shift) {

  /* Shift the particles positions to a local frame (ci frame) so single precision can be
   * used instead of double precision. Also shift the cell ci, particles positions due to BCs but leave cell cj. */
  for (int i = 0; i < ci->count; i++) {
    ci_cache->x[i] = ci->parts[i].x[0] - ci->loc[0] - shift[0];
    ci_cache->y[i] = ci->parts[i].x[1] - ci->loc[1] - shift[1];
    ci_cache->z[i] = ci->parts[i].x[2] - ci->loc[2] - shift[2];
    ci_cache->h[i] = ci->parts[i].h;

    ci_cache->m[i] = ci->parts[i].mass;
    ci_cache->vx[i] = ci->parts[i].v[0];
    ci_cache->vy[i] = ci->parts[i].v[1];
    ci_cache->vz[i] = ci->parts[i].v[2];
  }
  
  for (int i = 0; i < cj->count; i++) {
    cj_cache->x[i] = cj->parts[i].x[0] - ci->loc[0];
    cj_cache->y[i] = cj->parts[i].x[1] - ci->loc[1];
    cj_cache->z[i] = cj->parts[i].x[2] - ci->loc[2];
    cj_cache->h[i] = cj->parts[i].h;

    cj_cache->m[i] = cj->parts[i].mass;
    cj_cache->vx[i] = cj->parts[i].v[0];
    cj_cache->vy[i] = cj->parts[i].v[1];
    cj_cache->vz[i] = cj->parts[i].v[2];
  }
}

/**
 * @brief Populate cache by reading in the particles from two cells in sorted order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The cache for cell ci.
 * @param cj_cache The cache for cell cj.
 * @param shift The amount to shift the particle positions to account for BCs
 */
__attribute__((always_inline)) INLINE void cache_read_two_cells_sorted(
    const struct cell *const ci, const struct cell *const cj, struct cache *const ci_cache, struct cache *const cj_cache, const struct entry *restrict sort_i, const struct entry *restrict sort_j, const double *const shift) {

  int idx;
  /* Shift the particles positions to a local frame (ci frame) so single precision can be
   * used instead of double precision. Also shift the cell ci, particles positions due to BCs but leave cell cj. */
#ifdef WITH_VECTORIZATION
#pragma simd
#endif
  for (int i = 0; i < ci->count; i++) {
    idx = sort_i[i].i;
    ci_cache->x[i] = ci->parts[idx].x[0] - ci->loc[0] - shift[0];
    ci_cache->y[i] = ci->parts[idx].x[1] - ci->loc[1] - shift[1];
    ci_cache->z[i] = ci->parts[idx].x[2] - ci->loc[2] - shift[2];
    ci_cache->h[i] = ci->parts[idx].h;

    ci_cache->m[i] = ci->parts[idx].mass;
    ci_cache->vx[i] = ci->parts[idx].v[0];
    ci_cache->vy[i] = ci->parts[idx].v[1];
    ci_cache->vz[i] = ci->parts[idx].v[2];

#ifdef DOPAIR1_AUTO_VEC
    ci_cache->rho[i]         = 0.0f; 
    ci_cache->rho_dh[i]      = 0.0f; 
    ci_cache->wcount[i]      = 0.0f; 
    ci_cache->wcount_dh[i]   = 0.0f; 
    ci_cache->div_v[i]       = 0.0f; 
    ci_cache->curl_vx[i]     = 0.0f; 
    ci_cache->curl_vy[i]     = 0.0f; 
    ci_cache->curl_vz[i]     = 0.0f; 
#endif
  }
 
#ifdef WITH_VECTORIZATION
#pragma simd
#endif
  for (int i = 0; i < cj->count; i++) {
    idx = sort_j[i].i;
    cj_cache->x[i] = cj->parts[idx].x[0] - ci->loc[0];
    cj_cache->y[i] = cj->parts[idx].x[1] - ci->loc[1];
    cj_cache->z[i] = cj->parts[idx].x[2] - ci->loc[2];
    cj_cache->h[i] = cj->parts[idx].h;

    cj_cache->m[i] = cj->parts[idx].mass;
    cj_cache->vx[i] = cj->parts[idx].v[0];
    cj_cache->vy[i] = cj->parts[idx].v[1];
    cj_cache->vz[i] = cj->parts[idx].v[2];
#ifdef DOPAIR1_AUTO_VEC
    cj_cache->rho[i]         = 0.0f; 
    cj_cache->rho_dh[i]      = 0.0f; 
    cj_cache->wcount[i]      = 0.0f; 
    cj_cache->wcount_dh[i]   = 0.0f; 
    cj_cache->div_v[i]       = 0.0f; 
    cj_cache->curl_vx[i]     = 0.0f; 
    cj_cache->curl_vy[i]     = 0.0f; 
    cj_cache->curl_vz[i]     = 0.0f; 
#endif
  }
}

__attribute__((always_inline)) INLINE static void cache_write_sorted_particles(const struct cache *const ci_cache, const struct cache *const cj_cache, const struct cell *const ci, const struct cell *const cj, const struct entry *restrict sort_i, const struct entry *restrict sort_j) {

#ifdef DOPAIR1_AUTO_VEC
    struct part *restrict pi, *restrict pj;

    int idx = 0;
    for (int i=0; i<ci->count; i++) {
      idx = sort_i[i].i;  
      pi = &ci->parts[idx];
      pi->rho                 +=  ci_cache->rho[i];
      pi->density.rho_dh      +=  ci_cache->rho_dh[i];
      pi->density.wcount      +=  ci_cache->wcount[i];
      pi->density.wcount_dh   +=  ci_cache->wcount_dh[i];
      pi->density.div_v       +=  ci_cache->div_v[i];
      pi->density.rot_v[0]    +=  ci_cache->curl_vx[i];
      pi->density.rot_v[1]    +=  ci_cache->curl_vy[i];
      pi->density.rot_v[2]    +=  ci_cache->curl_vz[i];
    }

    for (int i=0; i<cj->count; i++) {
      idx = sort_j[i].i;  
      pj = &cj->parts[idx];
      pj->rho                 += cj_cache->rho[i];                     
      pj->density.rho_dh      += cj_cache->rho_dh[i];
      pj->density.wcount      += cj_cache->wcount[i];
      pj->density.wcount_dh   += cj_cache->wcount_dh[i];
      pj->density.div_v       += cj_cache->div_v[i];
      pj->density.rot_v[0]    += cj_cache->curl_vx[i];               
      pj->density.rot_v[1]    += cj_cache->curl_vy[i];               
      pj->density.rot_v[2]    += cj_cache->curl_vz[i];               
    }
#endif
}
#endif /* SWIFT_CACHE_H */
