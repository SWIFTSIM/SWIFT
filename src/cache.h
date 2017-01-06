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

#define NUM_VEC_PROC 2
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)
#define C2_CACHE_ALIGN sizeof(float) * VEC_SIZE

/* Cache struct to hold a local copy of a cells' particle
 * properties required for density/force calculations.*/
struct cache {

  /* Particle x position. */
  float *restrict x __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle y position. */
  float *restrict y __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle z position. */
  float *restrict z __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle smoothing length. */
  float *restrict h __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle mass. */
  float *restrict m __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle x velocity. */
  float *restrict vx __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle y velocity. */
  float *restrict vy __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  /* Particle z velocity. */
  float *restrict vz __attribute__((aligned(sizeof(float) * VEC_SIZE)));

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
  }

  error += posix_memalign((void **)&c->x, alignment, sizeBytes);
  error += posix_memalign((void **)&c->y, alignment, sizeBytes);
  error += posix_memalign((void **)&c->z, alignment, sizeBytes);
  error += posix_memalign((void **)&c->m, alignment, sizeBytes);
  error += posix_memalign((void **)&c->vx, alignment, sizeBytes);
  error += posix_memalign((void **)&c->vy, alignment, sizeBytes);
  error += posix_memalign((void **)&c->vz, alignment, sizeBytes);
  error += posix_memalign((void **)&c->h, alignment, sizeBytes);

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

#endif /* SWIFT_CACHE_H */
