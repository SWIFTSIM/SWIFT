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

  /* Particle density. */
  float *restrict rho SWIFT_CACHE_ALIGN;

  /* Particle smoothing length gradient. */
  float *restrict grad_h SWIFT_CACHE_ALIGN;

  /* Pressure over density squared. */
  float *restrict pOrho2 SWIFT_CACHE_ALIGN;

  /* Balsara switch. */
  float *restrict balsara SWIFT_CACHE_ALIGN;

  /* Particle sound speed. */
  float *restrict soundspeed SWIFT_CACHE_ALIGN;

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
   * the vector size and include 2 vector lengths for remainder operations. */
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
    free(c->rho);
    free(c->grad_h);
    free(c->pOrho2);
    free(c->balsara);
    free(c->soundspeed);
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
  error += posix_memalign((void **)&c->rho, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->grad_h, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->pOrho2, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->balsara, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->soundspeed, SWIFT_CACHE_ALIGNMENT, sizeBytes);

  if (error != 0)
    error("Couldn't allocate cache, no. of particles: %d", (int)count);
  c->count = count;
}

/**
 * @brief Populate cache by reading in the particles in unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 * @return uninhibited_count The no. of uninhibited particles.
 */
__attribute__((always_inline)) INLINE int cache_read_particles(
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

  const int count = ci->hydro.count;
  const struct part *restrict parts = ci->hydro.parts;
  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const double max_dx = ci->hydro.dx_max_part;
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};
  const float h_padded = ci->hydro.h_max / 4.;

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < count; i++) {

    /* Pad inhibited particles. */
    if (parts[i].time_bin >= time_bin_inhibited) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
      h[i] = h_padded;

      continue;
    }

    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;
    m[i] = parts[i].mass;
    vx[i] = parts[i].v[0];
    vy[i] = parts[i].v[1];
    vz[i] = parts[i].v[2];
  }

  /* Pad cache if the no. of particles is not a multiple of double the vector
   * length. */
  int count_align = count;
  const int rem = count % (NUM_VEC_PROC * VEC_SIZE);
  if (rem != 0) {
    count_align += (NUM_VEC_PROC * VEC_SIZE) - rem;

    /* Set positions to something outside of the range of any particle */
    for (int i = count; i < count_align; i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
    }
  }

  return count_align;

#else
  error("Can't call the cache reading function with this flavour of SPH!");
  return 0;
#endif
}

/**
 * @brief Populate cache by reading in the particles in unsorted order for
 * doself_subset.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 * @return uninhibited_count The no. of uninhibited particles.
 */
__attribute__((always_inline)) INLINE int cache_read_particles_subset_self(
    const struct cell *restrict const ci,
    struct cache *restrict const ci_cache) {

#if defined(GADGET2_SPH)

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vx, ci_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vy, ci_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vz, ci_cache->vz, SWIFT_CACHE_ALIGNMENT);

  const int count = ci->hydro.count;
  const struct part *restrict parts = ci->hydro.parts;
  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const double max_dx = ci->hydro.dx_max_part;
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < count; i++) {

    /* Pad inhibited particles. */
    if (parts[i].time_bin >= time_bin_inhibited) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];

      continue;
    }

    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    m[i] = parts[i].mass;
    vx[i] = parts[i].v[0];
    vy[i] = parts[i].v[1];
    vz[i] = parts[i].v[2];
  }

  /* Pad cache if the no. of particles is not a multiple of double the vector
   * length. */
  int count_align = count;
  const int rem = count % (NUM_VEC_PROC * VEC_SIZE);
  if (rem != 0) {
    count_align += (NUM_VEC_PROC * VEC_SIZE) - rem;

    /* Set positions to something outside of the range of any particle */
    for (int i = count; i < count_align; i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
    }
  }

  return count_align;

#else
  error("Can't call the cache reading function with this flavour of SPH!");
  return 0;
#endif
}

/**
 * @brief Populate cache by only reading particles that are within range of
 * each other within the adjoining cell. Also read the particles into the cache
 * in sorted order.
 *
 * @param ci The i #cell.
 * @param ci_cache The #cache for cell ci.
 * @param sort_i The array of sorted particle indices for cell ci.
 * @param first_pi The first particle in cell ci that is in range.
 * @param last_pi The last particle in cell ci that is in range.
 * @param last_pi The last particle in cell ci that is in range.
 * @param loc The cell location to remove from the particle positions.
 * @param flipped Flag to check whether the cells have been flipped or not.
 */
__attribute__((always_inline)) INLINE void cache_read_particles_subset_pair(
    const struct cell *restrict const ci, struct cache *restrict const ci_cache,
    const struct entry *restrict sort_i, int *first_pi, int *last_pi,
    const double *loc, const int flipped) {

#if defined(GADGET2_SPH)

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vx, ci_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vy, ci_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vz, ci_cache->vz, SWIFT_CACHE_ALIGNMENT);

  const struct part *restrict parts = ci->hydro.parts;

  /* The cell is on the right so read the particles
   * into the cache from the start of the cell. */
  if (!flipped) {
    const int rem = (*last_pi + 1) % VEC_SIZE;
    if (rem != 0) {
      const int pad = VEC_SIZE - rem;

      /* Increase last_pi if there are particles in the cell left to read. */
      if (*last_pi + pad < ci->hydro.count) *last_pi += pad;
    }

    const double max_dx = ci->hydro.dx_max_part;
    const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                                 -(2. * ci->width[1] + max_dx),
                                 -(2. * ci->width[2] + max_dx)};

    /* Shift the particles positions to a local frame so single precision can be
     * used instead of double precision. */
    for (int i = 0; i < *last_pi; i++) {
      const int idx = sort_i[i].i;

      /* Put inhibited particles out of range. */
      if (parts[idx].time_bin >= time_bin_inhibited) {
        x[i] = pos_padded[0];
        y[i] = pos_padded[1];
        z[i] = pos_padded[2];
        m[i] = 1.f;
        vx[i] = 1.f;
        vy[i] = 1.f;
        vz[i] = 1.f;

        continue;
      }

      x[i] = (float)(parts[idx].x[0] - loc[0]);
      y[i] = (float)(parts[idx].x[1] - loc[1]);
      z[i] = (float)(parts[idx].x[2] - loc[2]);
      m[i] = parts[idx].mass;
      vx[i] = parts[idx].v[0];
      vy[i] = parts[idx].v[1];
      vz[i] = parts[idx].v[2];
    }

    /* Pad cache with fake particles that exist outside the cell so will not
     * interact. We use values of the same magnitude (but negative!) as the real
     * particles to avoid overflow problems. */
    for (int i = *last_pi; i < *last_pi + VEC_SIZE; i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];

      m[i] = 1.f;
      vx[i] = 1.f;
      vy[i] = 1.f;
      vz[i] = 1.f;
    }
  }
  /* The cell is on the left so read the particles
   * into the cache from the end of the cell. */
  else {
    const int rem = (ci->hydro.count - *first_pi) % VEC_SIZE;
    if (rem != 0) {
      const int pad = VEC_SIZE - rem;

      /* Decrease first_pi if there are particles in the cell left to read. */
      if (*first_pi - pad >= 0) *first_pi -= pad;
    }

    const int ci_cache_count = ci->hydro.count - *first_pi;
    const double max_dx = ci->hydro.dx_max_part;
    const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                                 -(2. * ci->width[1] + max_dx),
                                 -(2. * ci->width[2] + max_dx)};

    /* Shift the particles positions to a local frame so single precision can be
     * used instead of double precision. */
    for (int i = 0; i < ci_cache_count; i++) {
      const int idx = sort_i[i + *first_pi].i;

      /* Put inhibited particles out of range. */
      if (parts[idx].time_bin >= time_bin_inhibited) {
        x[i] = pos_padded[0];
        y[i] = pos_padded[1];
        z[i] = pos_padded[2];

        m[i] = 1.f;
        vx[i] = 1.f;
        vy[i] = 1.f;
        vz[i] = 1.f;

        continue;
      }

      x[i] = (float)(parts[idx].x[0] - loc[0]);
      y[i] = (float)(parts[idx].x[1] - loc[1]);
      z[i] = (float)(parts[idx].x[2] - loc[2]);
      m[i] = parts[idx].mass;
      vx[i] = parts[idx].v[0];
      vy[i] = parts[idx].v[1];
      vz[i] = parts[idx].v[2];
    }

    /* Pad cache with fake particles that exist outside the cell so will not
     * interact. We use values of the same magnitude (but negative!) as the real
     * particles to avoid overflow problems. */
    for (int i = ci->hydro.count - *first_pi;
         i < ci->hydro.count - *first_pi + VEC_SIZE; i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];

      m[i] = 1.f;
      vx[i] = 1.f;
      vy[i] = 1.f;
      vz[i] = 1.f;
    }
  }

#endif
}

/**
 * @brief Populate cache for force interactions by reading in the particles in
 * unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 * @return uninhibited_count The no. of uninhibited particles.
 */
__attribute__((always_inline)) INLINE int cache_read_force_particles(
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
  swift_declare_aligned_ptr(float, rho, ci_cache->rho, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, grad_h, ci_cache->grad_h,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pOrho2, ci_cache->pOrho2,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, balsara, ci_cache->balsara,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, soundspeed, ci_cache->soundspeed,
                            SWIFT_CACHE_ALIGNMENT);

  const int count = ci->hydro.count;
  const struct part *restrict parts = ci->hydro.parts;
  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const double max_dx = ci->hydro.dx_max_part;
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};
  const float h_padded = ci->hydro.h_max / 4.;

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < count; i++) {

    /* Skip inhibited particles. */
    if (parts[i].time_bin >= time_bin_inhibited) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
      h[i] = h_padded;
      rho[i] = 1.f;
      grad_h[i] = 1.f;
      pOrho2[i] = 1.f;
      balsara[i] = 1.f;
      soundspeed[i] = 1.f;

      continue;
    }

    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;
    m[i] = parts[i].mass;
    vx[i] = parts[i].v[0];
    vy[i] = parts[i].v[1];
    vz[i] = parts[i].v[2];
    rho[i] = parts[i].rho;
    grad_h[i] = parts[i].force.f;
    pOrho2[i] = parts[i].force.P_over_rho2;
    balsara[i] = parts[i].force.balsara;
    soundspeed[i] = parts[i].force.soundspeed;
  }

  /* Pad cache if there is a serial remainder. */
  int count_align = count;
  const int rem = count % VEC_SIZE;
  if (rem != 0) {
    count_align += VEC_SIZE - rem;

    /* Set positions to the same as particle pi so when the r2 > 0 mask is
     * applied these extra contributions are masked out.*/
    for (int i = count; i < count_align; i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
      h[i] = h_padded;
      rho[i] = 1.f;
      grad_h[i] = 1.f;
      pOrho2[i] = 1.f;
      balsara[i] = 1.f;
      soundspeed[i] = 1.f;
    }
  }

  return count_align;

#else
  error("Can't call the cache reading function with this flavour of SPH!");
  return 0;
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
 */
__attribute__((always_inline)) INLINE void cache_read_two_partial_cells_sorted(
    const struct cell *restrict const ci, const struct cell *restrict const cj,
    struct cache *restrict const ci_cache,
    struct cache *restrict const cj_cache, const struct entry *restrict sort_i,
    const struct entry *restrict sort_j, const double *restrict const shift,
    int *first_pi, int *last_pj) {

  /* Make the number of particles to be read a multiple of the vector size.
   * This eliminates serial remainder loops where possible when populating the
   * cache. */

  /* Is the number of particles to read a multiple of the vector size? */
  int rem = (ci->hydro.count - *first_pi) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Decrease first_pi if there are particles in the cell left to read. */
    if (*first_pi - pad >= 0) *first_pi -= pad;
  }

  rem = (*last_pj + 1) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Increase last_pj if there are particles in the cell left to read. */
    if (*last_pj + pad < cj->hydro.count) *last_pj += pad;
  }

  /* Get some local pointers */
  const int first_pi_align = *first_pi;
  const int last_pj_align = *last_pj;
  const struct part *restrict parts_i = ci->hydro.parts;
  const struct part *restrict parts_j = cj->hydro.parts;

  /* Shift particles to the local frame and account for boundary conditions.*/
  const double total_ci_shift[3] = {
      cj->loc[0] + shift[0], cj->loc[1] + shift[1], cj->loc[2] + shift[2]};
  const double total_cj_shift[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

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

  int ci_cache_count = ci->hydro.count - first_pi_align;
  const double max_dx = max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float pos_padded_i[3] = {-(2. * ci->width[0] + max_dx),
                                 -(2. * ci->width[1] + max_dx),
                                 -(2. * ci->width[2] + max_dx)};
  const float h_padded_i = ci->hydro.h_max / 4.;

  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be used instead of double precision.  */
  for (int i = 0; i < ci_cache_count; i++) {
    const int idx = sort_i[i + first_pi_align].i;

    /* Put inhibited particles out of range. */
    if (parts_i[idx].time_bin >= time_bin_inhibited) {
      x[i] = pos_padded_i[0];
      y[i] = pos_padded_i[1];
      z[i] = pos_padded_i[2];
      h[i] = h_padded_i;

      m[i] = 1.f;
      vx[i] = 1.f;
      vy[i] = 1.f;
      vz[i] = 1.f;

      continue;
    }

    x[i] = (float)(parts_i[idx].x[0] - total_ci_shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - total_ci_shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - total_ci_shift[2]);
    h[i] = parts_i[idx].h;
    vx[i] = parts_i[idx].v[0];
    vy[i] = parts_i[idx].v[1];
    vz[i] = parts_i[idx].v[2];
#ifdef GADGET2_SPH
    m[i] = parts_i[idx].mass;
#endif
  }

#ifdef SWIFT_DEBUG_CHECKS
  const float shift_threshold_x =
      2. * ci->width[0] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_y =
      2. * ci->width[1] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_z =
      2. * ci->width[2] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);

  /* Make sure that particle positions have been shifted correctly. */
  for (int i = 0; i < ci_cache_count; i++) {
    if (x[i] > shift_threshold_x || x[i] < -shift_threshold_x)
      error(
          "Error: ci->loc[%lf,%lf,%lf],cj->loc[%lf,%lf,%lf] Particle %d x pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. x=%f, ci->width[0]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, x[i], ci->width[0]);
    if (y[i] > shift_threshold_y || y[i] < -shift_threshold_y)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d y pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. y=%f, ci->width[1]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, y[i], ci->width[1]);
    if (z[i] > shift_threshold_z || z[i] < -shift_threshold_z)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d z pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. z=%f, ci->width[2]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, z[i], ci->width[2]);
  }
#endif

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  for (int i = ci->hydro.count - first_pi_align;
       i < ci->hydro.count - first_pi_align + VEC_SIZE; i++) {
    x[i] = pos_padded_i[0];
    y[i] = pos_padded_i[1];
    z[i] = pos_padded_i[2];
    h[i] = h_padded_i;

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

  const float pos_padded_j[3] = {-(2. * cj->width[0] + max_dx),
                                 -(2. * cj->width[1] + max_dx),
                                 -(2. * cj->width[2] + max_dx)};
  const float h_padded_j = cj->hydro.h_max / 4.;

  for (int i = 0; i <= last_pj_align; i++) {
    const int idx = sort_j[i].i;

    /* Put inhibited particles out of range. */
    if (parts_j[idx].time_bin >= time_bin_inhibited) {
      xj[i] = pos_padded_j[0];
      yj[i] = pos_padded_j[1];
      zj[i] = pos_padded_j[2];
      hj[i] = h_padded_j;

      mj[i] = 1.f;
      vxj[i] = 1.f;
      vyj[i] = 1.f;
      vzj[i] = 1.f;

      continue;
    }

    xj[i] = (float)(parts_j[idx].x[0] - total_cj_shift[0]);
    yj[i] = (float)(parts_j[idx].x[1] - total_cj_shift[1]);
    zj[i] = (float)(parts_j[idx].x[2] - total_cj_shift[2]);
    hj[i] = parts_j[idx].h;
    vxj[i] = parts_j[idx].v[0];
    vyj[i] = parts_j[idx].v[1];
    vzj[i] = parts_j[idx].v[2];
#ifdef GADGET2_SPH
    mj[i] = parts_j[idx].mass;
#endif
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure that particle positions have been shifted correctly. */
  for (int i = 0; i <= last_pj_align; i++) {
    if (xj[i] > shift_threshold_x || xj[i] < -shift_threshold_x)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d xj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. xj=%f, ci->width[0]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, xj[i], ci->width[0]);
    if (yj[i] > shift_threshold_y || yj[i] < -shift_threshold_y)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d yj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. yj=%f, ci->width[1]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, yj[i], ci->width[1]);
    if (zj[i] > shift_threshold_z || zj[i] < -shift_threshold_z)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d zj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. zj=%f, ci->width[2]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, zj[i], ci->width[2]);
  }
#endif

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  for (int i = last_pj_align + 1; i < last_pj_align + 1 + VEC_SIZE; i++) {
    xj[i] = pos_padded_j[0];
    yj[i] = pos_padded_j[1];
    zj[i] = pos_padded_j[2];
    hj[i] = h_padded_j;
    mj[i] = 1.f;
    vxj[i] = 1.f;
    vyj[i] = 1.f;
    vzj[i] = 1.f;
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
 */
__attribute__((always_inline)) INLINE void
cache_read_two_partial_cells_sorted_force(
    const struct cell *const ci, const struct cell *const cj,
    struct cache *const ci_cache, struct cache *const cj_cache,
    const struct entry *restrict sort_i, const struct entry *restrict sort_j,
    const double *const shift, int *first_pi, int *last_pj) {

  /* Make the number of particles to be read a multiple of the vector size.
   * This eliminates serial remainder loops where possible when populating the
   * cache. */

  /* Is the number of particles to read a multiple of the vector size? */
  int rem = (ci->hydro.count - *first_pi) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Decrease first_pi if there are particles in the cell left to read. */
    if (*first_pi - pad >= 0) *first_pi -= pad;
  }

  rem = (*last_pj + 1) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Increase last_pj if there are particles in the cell left to read. */
    if (*last_pj + pad < cj->hydro.count) *last_pj += pad;
  }

  /* Get some local pointers */
  const int first_pi_align = *first_pi;
  const int last_pj_align = *last_pj;
  const struct part *restrict parts_i = ci->hydro.parts;
  const struct part *restrict parts_j = cj->hydro.parts;

  /* Shift particles to the local frame and account for boundary conditions.*/
  const double total_ci_shift[3] = {
      cj->loc[0] + shift[0], cj->loc[1] + shift[1], cj->loc[2] + shift[2]};
  const double total_cj_shift[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

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
  swift_declare_aligned_ptr(float, rho, ci_cache->rho, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, grad_h, ci_cache->grad_h,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pOrho2, ci_cache->pOrho2,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, balsara, ci_cache->balsara,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, soundspeed, ci_cache->soundspeed,
                            SWIFT_CACHE_ALIGNMENT);

  int ci_cache_count = ci->hydro.count - first_pi_align;
  const double max_dx = max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float pos_padded_i[3] = {-(2. * ci->width[0] + max_dx),
                                 -(2. * ci->width[1] + max_dx),
                                 -(2. * ci->width[2] + max_dx)};
  const float h_padded_i = ci->hydro.h_max / 4.;

  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be  used instead of double precision.  */
  for (int i = 0; i < ci_cache_count; i++) {

    const int idx = sort_i[i + first_pi_align].i;

    /* Put inhibited particles out of range. */
    if (parts_i[idx].time_bin >= time_bin_inhibited) {
      x[i] = pos_padded_i[0];
      y[i] = pos_padded_i[1];
      z[i] = pos_padded_i[2];
      h[i] = h_padded_i;
      m[i] = 1.f;
      vx[i] = 1.f;
      vy[i] = 1.f;
      vz[i] = 1.f;
      rho[i] = 1.f;
      grad_h[i] = 1.f;
      pOrho2[i] = 1.f;
      balsara[i] = 1.f;
      soundspeed[i] = 1.f;

      continue;
    }

    x[i] = (float)(parts_i[idx].x[0] - total_ci_shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - total_ci_shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - total_ci_shift[2]);
    h[i] = parts_i[idx].h;
    vx[i] = parts_i[idx].v[0];
    vy[i] = parts_i[idx].v[1];
    vz[i] = parts_i[idx].v[2];
#ifdef GADGET2_SPH
    m[i] = parts_i[idx].mass;
    rho[i] = parts_i[idx].rho;
    grad_h[i] = parts_i[idx].force.f;
    pOrho2[i] = parts_i[idx].force.P_over_rho2;
    balsara[i] = parts_i[idx].force.balsara;
    soundspeed[i] = parts_i[idx].force.soundspeed;
#endif
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  for (int i = ci->hydro.count - first_pi_align;
       i < ci->hydro.count - first_pi_align + VEC_SIZE; i++) {
    x[i] = pos_padded_i[0];
    y[i] = pos_padded_i[1];
    z[i] = pos_padded_i[2];
    h[i] = h_padded_i;
    m[i] = 1.f;
    vx[i] = 1.f;
    vy[i] = 1.f;
    vz[i] = 1.f;
    rho[i] = 1.f;
    grad_h[i] = 1.f;
    pOrho2[i] = 1.f;
    balsara[i] = 1.f;
    soundspeed[i] = 1.f;
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
  swift_declare_aligned_ptr(float, rhoj, cj_cache->rho, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, grad_hj, cj_cache->grad_h,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pOrho2j, cj_cache->pOrho2,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, balsaraj, cj_cache->balsara,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, soundspeedj, cj_cache->soundspeed,
                            SWIFT_CACHE_ALIGNMENT);

  const float pos_padded_j[3] = {-(2. * cj->width[0] + max_dx),
                                 -(2. * cj->width[1] + max_dx),
                                 -(2. * cj->width[2] + max_dx)};
  const float h_padded_j = cj->hydro.h_max / 4.;

  for (int i = 0; i <= last_pj_align; i++) {
    const int idx = sort_j[i].i;

    /* Put inhibited particles out of range. */
    if (parts_j[idx].time_bin == time_bin_inhibited) {
      xj[i] = pos_padded_j[0];
      yj[i] = pos_padded_j[1];
      zj[i] = pos_padded_j[2];
      hj[i] = h_padded_j;
      mj[i] = 1.f;
      vxj[i] = 1.f;
      vyj[i] = 1.f;
      vzj[i] = 1.f;
      rhoj[i] = 1.f;
      grad_hj[i] = 1.f;
      pOrho2j[i] = 1.f;
      balsaraj[i] = 1.f;
      soundspeedj[i] = 1.f;

      continue;
    }

    xj[i] = (float)(parts_j[idx].x[0] - total_cj_shift[0]);
    yj[i] = (float)(parts_j[idx].x[1] - total_cj_shift[1]);
    zj[i] = (float)(parts_j[idx].x[2] - total_cj_shift[2]);
    hj[i] = parts_j[idx].h;
    vxj[i] = parts_j[idx].v[0];
    vyj[i] = parts_j[idx].v[1];
    vzj[i] = parts_j[idx].v[2];
#ifdef GADGET2_SPH
    mj[i] = parts_j[idx].mass;
    rhoj[i] = parts_j[idx].rho;
    grad_hj[i] = parts_j[idx].force.f;
    pOrho2j[i] = parts_j[idx].force.P_over_rho2;
    balsaraj[i] = parts_j[idx].force.balsara;
    soundspeedj[i] = parts_j[idx].force.soundspeed;
#endif
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  for (int i = last_pj_align + 1; i < last_pj_align + 1 + VEC_SIZE; i++) {
    xj[i] = pos_padded_j[0];
    yj[i] = pos_padded_j[1];
    zj[i] = pos_padded_j[2];
    hj[i] = h_padded_j;
    mj[i] = 1.f;
    vxj[i] = 1.f;
    vyj[i] = 1.f;
    vzj[i] = 1.f;
    rhoj[i] = 1.f;
    grad_hj[i] = 1.f;
    pOrho2j[i] = 1.f;
    balsaraj[i] = 1.f;
    soundspeedj[i] = 1.f;
  }
}

/**
 * @brief Clean the memory allocated by a #cache object.
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
    free(c->rho);
    free(c->grad_h);
    free(c->pOrho2);
    free(c->balsara);
    free(c->soundspeed);
  }
  c->count = 0;
}

#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_CACHE_H */
