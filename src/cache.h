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
#include "sort_part.h"
#include "vector.h"

#define NUM_VEC_PROC 2
#define CACHE_ALIGN 64
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)
#define C2_CACHE_ALIGN sizeof(float) * VEC_SIZE

#ifdef WITH_VECTORIZATION
/* Cache struct to hold a local copy of a cells' particle
 * properties required for density/force calculations.*/
struct cache {

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

  /* Maximum distance of particles into neighbouring cell. */
  float *restrict max_d __attribute__((aligned(CACHE_ALIGN)));

  /* Cache size. */
  int count;
};

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

  /* Align cache on correct byte boundary and pad cache size to be a multiple of
   * the vector size
   * and include 2 vector lengths for remainder operations. */
  unsigned int pad = 2 * VEC_SIZE, rem = count % VEC_SIZE;
  if (rem > 0) pad += VEC_SIZE - rem;
  unsigned int sizeBytes = (count + pad) * sizeof(float);
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
    free(c->max_d);
  }

  error += posix_memalign((void **)&c->x, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->y, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->z, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->m, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->vx, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->vy, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->vz, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->h, CACHE_ALIGN, sizeBytes);
  error += posix_memalign((void **)&c->max_d, CACHE_ALIGN, sizeBytes);

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

#if defined(GADGET2_SPH)

/* Shift the particles positions to a local frame so single precision can be
 * used instead of double precision. */
#if defined(WITH_VECTORIZATION) && defined(__ICC)
#pragma vector aligned
#endif
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

#endif
}

/**
 * @brief Populate cache by reading in the particles from two cells in unsorted
 * order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The cache for cell ci.
 * @param cj_cache The cache for cell cj.
 * @param shift The amount to shift the particle positions to account for BCs
 */
__attribute__((always_inline)) INLINE void cache_read_two_cells(
    const struct cell *const ci, const struct cell *const cj,
    struct cache *const ci_cache, struct cache *const cj_cache,
    const double *const shift) {

  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be
   * used instead of double precision. Also shift the cell ci, particles
   * positions due to BCs but leave cell cj. */
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

__attribute__((always_inline)) INLINE void cache_read_cell_sorted(
    const struct cell *const ci, struct cache *const ci_cache,
    const struct entry *restrict sort_i, double *const loc,
    double *const shift) {

  int idx;
/* Shift the particles positions to a local frame (ci frame) so single precision
 * can be
 * used instead of double precision. Also shift the cell ci, particles positions
 * due to BCs but leave cell cj. */
#if defined(WITH_VECTORIZATION) && defined(__ICC)
#pragma simd
#endif
  for (int i = 0; i < ci->count; i++) {
    idx = sort_i[i].i;

    ci_cache->x[i] = ci->parts[idx].x[0] - loc[0] - shift[0];
    ci_cache->y[i] = ci->parts[idx].x[1] - loc[1] - shift[1];
    ci_cache->z[i] = ci->parts[idx].x[2] - loc[2] - shift[2];
    ci_cache->h[i] = ci->parts[idx].h;

    ci_cache->m[i] = ci->parts[idx].mass;
    ci_cache->vx[i] = ci->parts[idx].v[0];
    ci_cache->vy[i] = ci->parts[idx].v[1];
    ci_cache->vz[i] = ci->parts[idx].v[2];
  }
}

/**
 * @brief Populate cache by reading in the particles from two cells in sorted
 * order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The #cache for cell ci.
 * @param cj_cache The #cache for cell cj.
 * @param sort_i The array of sorted particle indices for cell ci.
 * @param sort_j The array of sorted particle indices for cell ci.
 * @param shift The amount to shift the particle positions to account for BCs
 */
__attribute__((always_inline)) INLINE void cache_read_two_cells_sorted(
    const struct cell *const ci, const struct cell *const cj,
    struct cache *const ci_cache, struct cache *const cj_cache,
    const struct entry *restrict sort_i, const struct entry *restrict sort_j,
    const double *const shift) {

  int idx;
/* Shift the particles positions to a local frame (ci frame) so single precision
 * can be
 * used instead of double precision. Also shift the cell ci, particles positions
 * due to BCs but leave cell cj. */
#if defined(WITH_VECTORIZATION) && defined(__ICC)
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
  }

#if defined(WITH_VECTORIZATION) && defined(__ICC)
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
  }
}

/**
 * @brief Populate caches by only reading particles that are within range of
 * each other within the adjoining cell.Also read the particles into the cache
 * in sorted order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The #cache for cell ci.
 * @param cj_cache The #cache for cell cj.
 * @param sort_i The array of sorted particle indices for cell ci.
 * @param sort_j The array of sorted particle indices for cell ci.
 * @param shift The amount to shift the particle positions to account for BCs
 * @param first_pi The first particle in cell ci that is in range.
 * @param last_pj The last particle in cell cj that is in range.
 * @param num_vec_proc Number of vectors that will be used to process the
 * interaction.
 */
__attribute__((always_inline)) INLINE void cache_read_two_partial_cells_sorted(
    const struct cell *const ci, const struct cell *const cj,
    struct cache *const ci_cache, struct cache *const cj_cache,
    const struct entry *restrict sort_i, const struct entry *restrict sort_j,
    const double *const shift, int *first_pi, int *last_pj,
    const int num_vec_proc) {

  int idx, ci_cache_idx;
  /* Pad number of particles read to the vector size. */
  int rem = (ci->count - *first_pi) % (num_vec_proc * VEC_SIZE);
  if (rem != 0) {
    int pad = (num_vec_proc * VEC_SIZE) - rem;

    if (*first_pi - pad >= 0) *first_pi -= pad;
  }

  rem = *last_pj % (num_vec_proc * VEC_SIZE);
  if (rem != 0) {
    int pad = (num_vec_proc * VEC_SIZE) - rem;

    if (*last_pj + pad < cj->count) *last_pj += pad;
  }

  int first_pi_align = *first_pi;
  int last_pj_align = *last_pj;

/* Shift the particles positions to a local frame (ci frame) so single precision
 * can be
 * used instead of double precision. Also shift the cell ci, particles positions
 * due to BCs but leave cell cj. */
#if defined(WITH_VECTORIZATION) && defined(__ICC)
#pragma vector aligned
#endif
  for (int i = first_pi_align; i < ci->count; i++) {
    /* Make sure ci_cache is filled from the first element. */
    ci_cache_idx = i - first_pi_align;
    idx = sort_i[i].i;
    ci_cache->x[ci_cache_idx] = ci->parts[idx].x[0] - ci->loc[0] - shift[0];
    ci_cache->y[ci_cache_idx] = ci->parts[idx].x[1] - ci->loc[1] - shift[1];
    ci_cache->z[ci_cache_idx] = ci->parts[idx].x[2] - ci->loc[2] - shift[2];
    ci_cache->h[ci_cache_idx] = ci->parts[idx].h;

    ci_cache->m[ci_cache_idx] = ci->parts[idx].mass;
    ci_cache->vx[ci_cache_idx] = ci->parts[idx].v[0];
    ci_cache->vy[ci_cache_idx] = ci->parts[idx].v[1];
    ci_cache->vz[ci_cache_idx] = ci->parts[idx].v[2];
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact.*/
  float fake_pix = 2.0f * ci->parts[sort_i[ci->count - 1].i].x[0];
  for (int i = ci->count - first_pi_align;
       i < ci->count - first_pi_align + VEC_SIZE; i++) {
    ci_cache->x[i] = fake_pix;
    ci_cache->y[i] = 1.f;
    ci_cache->z[i] = 1.f;
    ci_cache->h[i] = 1.f;

    ci_cache->m[i] = 1.f;
    ci_cache->vx[i] = 1.f;
    ci_cache->vy[i] = 1.f;
    ci_cache->vz[i] = 1.f;
  }

#if defined(WITH_VECTORIZATION) && defined(__ICC)
#pragma vector aligned
#endif
  for (int i = 0; i <= last_pj_align; i++) {
    idx = sort_j[i].i;
    cj_cache->x[i] = cj->parts[idx].x[0] - ci->loc[0];
    cj_cache->y[i] = cj->parts[idx].x[1] - ci->loc[1];
    cj_cache->z[i] = cj->parts[idx].x[2] - ci->loc[2];
    cj_cache->h[i] = cj->parts[idx].h;

    cj_cache->m[i] = cj->parts[idx].mass;
    cj_cache->vx[i] = cj->parts[idx].v[0];
    cj_cache->vy[i] = cj->parts[idx].v[1];
    cj_cache->vz[i] = cj->parts[idx].v[2];
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact.*/
  float fake_pjx = 2.0f * cj->parts[sort_j[cj->count - 1].i].x[0];
  for (int i = last_pj_align + 1; i < last_pj_align + 1 + VEC_SIZE; i++) {
    cj_cache->x[i] = fake_pjx;
    cj_cache->y[i] = 1.f;
    cj_cache->z[i] = 1.f;
    cj_cache->h[i] = 1.f;

    cj_cache->m[i] = 1.f;
    cj_cache->vx[i] = 1.f;
    cj_cache->vy[i] = 1.f;
    cj_cache->vz[i] = 1.f;
  }
}

/* @brief Clean the memory allocated by a #cache object.
 *
 * @param c The #cache to clean.
 */
static INLINE void cache_clean(struct cache *c) {
  if (c->count > 0) {
    free(c->x);
    free(c->y);
    free(c->z);
    free(c->m);
    free(c->vx);
    free(c->vy);
    free(c->vz);
    free(c->h);
    free(c->max_d);
  }
}

#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_CACHE_H */
