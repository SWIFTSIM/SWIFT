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
#include "align.h"
#include "cell.h"
#include "error.h"
#include "part.h"
#include "sort_part.h"
#include "vector.h"

#define NUM_VEC_PROC 2
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)

#ifdef WITH_VECTORIZATION
/* Cache struct to hold a local copy of a cells' particle
 * properties required for density/force calculations.*/
struct cache {

  /* Particle x position. */
  float *restrict x SWIFT_CACHE_ALIGN;

  /* Particle y position. */
  float *restrict y SWIFT_CACHE_ALIGN;

  /* Particle z position. */
  float *restrict z SWIFT_CACHE_ALIGN;

  /* Particle smoothing length. */
  float *restrict h SWIFT_CACHE_ALIGN;

  /* Particle mass. */
  float *restrict m SWIFT_CACHE_ALIGN;

  /* Particle x velocity. */
  float *restrict vx SWIFT_CACHE_ALIGN;

  /* Particle y velocity. */
  float *restrict vy SWIFT_CACHE_ALIGN;

  /* Particle z velocity. */
  float *restrict vz SWIFT_CACHE_ALIGN;

  /* Maximum index into neighbouring cell for particles that are in range. */
  int *restrict max_index SWIFT_CACHE_ALIGN;

  /* Cache size. */
  int count;
};

/* Secondary cache struct to hold a list of interactions between two
 * particles.*/
struct c2_cache {

  /* Separation between two particles squared. */
  float r2q[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* x separation between two particles. */
  float dxq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* y separation between two particles. */
  float dyq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* z separation between two particles. */
  float dzq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* Mass of particle pj. */
  float mq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* x velocity of particle pj. */
  float vxq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* y velocity of particle pj. */
  float vyq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* z velocity of particle pj. */
  float vzq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;
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
  size_t pad = 2 * VEC_SIZE, rem = count % VEC_SIZE;
  if (rem > 0) pad += VEC_SIZE - rem;
  size_t sizeBytes = (count + pad) * sizeof(float);
  size_t sizeIntBytes = (count + pad) * sizeof(int);
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
    free(c->max_index);
  }

  error += posix_memalign((void **)&c->x, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->y, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->z, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->m, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->vx, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->vy, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->vz, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->h, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->max_index, SWIFT_CACHE_ALIGNMENT,
                          sizeIntBytes);

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
    const struct cell *restrict const ci,
    struct cache *restrict const ci_cache) {

#if defined(GADGET2_SPH)

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vx, ci_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vy, ci_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vz, ci_cache->vz, SWIFT_CACHE_ALIGNMENT);

  const struct part *restrict parts = ci->parts;
  double loc[3];
  loc[0] = ci->loc[0];
  loc[1] = ci->loc[1];
  loc[2] = ci->loc[2];

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < ci->count; i++) {
    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;

    m[i] = parts[i].mass;
    vx[i] = parts[i].v[0];
    vy[i] = parts[i].v[1];
    vz[i] = parts[i].v[2];
  }

#endif
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
    const struct cell *restrict const ci, const struct cell *restrict const cj,
    struct cache *restrict const ci_cache,
    struct cache *restrict const cj_cache, const struct entry *restrict sort_i,
    const struct entry *restrict sort_j, const double *restrict const shift,
    int *first_pi, int *last_pj, const int num_vec_proc) {

  int idx;
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
  const struct part *restrict parts_i = ci->parts;
  const struct part *restrict parts_j = cj->parts;
  double loc[3];
  loc[0] = ci->loc[0];
  loc[1] = ci->loc[1];
  loc[2] = ci->loc[2];

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vx, ci_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vy, ci_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vz, ci_cache->vz, SWIFT_CACHE_ALIGNMENT);

  int ci_cache_count = ci->count - first_pi_align;
  /* Shift the particles positions to a local frame (ci frame) so single
   * precision
   * can be
   * used instead of double precision. Also shift the cell ci, particles
   * positions
   * due to BCs but leave cell cj. */
  for (int i = 0; i < ci_cache_count; i++) {
    idx = sort_i[i + first_pi_align].i;
    x[i] = (float)(parts_i[idx].x[0] - loc[0] - shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - loc[1] - shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - loc[2] - shift[2]);
    h[i] = parts_i[idx].h;

    m[i] = parts_i[idx].mass;
    vx[i] = parts_i[idx].v[0];
    vy[i] = parts_i[idx].v[1];
    vz[i] = parts_i[idx].v[2];
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact.*/
  float fake_pix = 2.0f * parts_i[sort_i[ci->count - 1].i].x[0];
  for (int i = ci->count - first_pi_align;
       i < ci->count - first_pi_align + VEC_SIZE; i++) {
    x[i] = fake_pix;
    y[i] = 1.f;
    z[i] = 1.f;
    h[i] = 1.f;

    m[i] = 1.f;
    vx[i] = 1.f;
    vy[i] = 1.f;
    vz[i] = 1.f;
  }

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, xj, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, yj, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, zj, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, hj, cj_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, mj, cj_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vxj, cj_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vyj, cj_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vzj, cj_cache->vz, SWIFT_CACHE_ALIGNMENT);

  for (int i = 0; i <= last_pj_align; i++) {
    idx = sort_j[i].i;
    xj[i] = (float)(parts_j[idx].x[0] - loc[0]);
    yj[i] = (float)(parts_j[idx].x[1] - loc[1]);
    zj[i] = (float)(parts_j[idx].x[2] - loc[2]);
    hj[i] = parts_j[idx].h;

    mj[i] = parts_j[idx].mass;
    vxj[i] = parts_j[idx].v[0];
    vyj[i] = parts_j[idx].v[1];
    vzj[i] = parts_j[idx].v[2];
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact.*/
  float fake_pjx = 2.0f * cj->parts[sort_j[cj->count - 1].i].x[0];
  for (int i = last_pj_align + 1; i < last_pj_align + 1 + VEC_SIZE; i++) {
    xj[i] = fake_pjx;
    yj[i] = 1.f;
    zj[i] = 1.f;
    hj[i] = 1.f;

    mj[i] = 1.f;
    vxj[i] = 1.f;
    vyj[i] = 1.f;
    vzj[i] = 1.f;
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
    free(c->max_index);
  }
}

#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_CACHE_H */
