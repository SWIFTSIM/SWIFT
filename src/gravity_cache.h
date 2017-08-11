/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_GRAVITY_CACHE_H
#define SWIFT_GRAVITY_CACHE_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "align.h"
#include "error.h"
#include "gravity.h"
#include "vector.h"

/**
 * @brief A SoA object for the #gpart of a cell.
 *
 * This is used to help vectorize the leaf-leaf gravity interactions.
 */
struct gravity_cache {

  /*! #gpart x position. */
  float *restrict x SWIFT_CACHE_ALIGN;

  /*! #gpart y position. */
  float *restrict y SWIFT_CACHE_ALIGN;

  /*! #gpart z position. */
  float *restrict z SWIFT_CACHE_ALIGN;

  /*! #gpart softening length. */
  float *restrict epsilon SWIFT_CACHE_ALIGN;

  /*! #gpart mass. */
  float *restrict m SWIFT_CACHE_ALIGN;

  /*! #gpart x acceleration. */
  float *restrict a_x SWIFT_CACHE_ALIGN;

  /*! #gpart y acceleration. */
  float *restrict a_y SWIFT_CACHE_ALIGN;

  /*! #gpart z acceleration. */
  float *restrict a_z SWIFT_CACHE_ALIGN;

  /*! Cache size */
  int count;
};

/**
 * @brief Frees the memory allocated in a #gravity_cache
 *
 * @param c The #gravity_cache to free.
 */
static INLINE void gravity_cache_clean(struct gravity_cache *c) {

  if (c->count > 0) {
    free(c->x);
    free(c->y);
    free(c->z);
    free(c->epsilon);
    free(c->m);
    free(c->a_x);
    free(c->a_y);
    free(c->a_z);
  }
  c->count = 0;
}

/**
 * @brief Allocates memory for the #gpart caches used in the leaf-leaf
 * interactions.
 *
 * The cache is padded for the vector size and aligned properly
 *
 * @param c The #gravity_cache to allocate.
 * @param count The number of #gpart to allocated for (space_splitsize is a good
 * choice).
 */
static INLINE void gravity_cache_init(struct gravity_cache *c, int count) {

  /* Size of the gravity cache */
  const int padded_count = count - (count % VEC_SIZE) + VEC_SIZE;
  const size_t sizeBytes = padded_count * sizeof(float);

  /* Delete old stuff if any */
  gravity_cache_clean(c);

  int error = 0;
  error += posix_memalign((void **)&c->x, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->y, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->z, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->epsilon, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->m, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->a_x, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->a_y, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->a_z, SWIFT_CACHE_ALIGNMENT, sizeBytes);

  if (error != 0)
    error("Couldn't allocate gravity cache, size: %d", padded_count);

  c->count = padded_count;
}

/**
 * @brief Fills a #gravity_cache structure with some #gpart and shift them.
 *
 * @param c The #gravity_cache to fill.
 * @param gparts The #gpart array to read from.
 * @param gcount The number of particles to read.
 * @param gcount_padded The number of particle to read padded to the next
 * multiple of the vector length.
 * @param shift A shift to apply to all the particles.
 */
__attribute__((always_inline)) INLINE void gravity_cache_populate(
    struct gravity_cache *c, const struct gpart *restrict gparts, int gcount,
    int gcount_padded, const double shift[3]) {

  /* Make the compiler understand we are in happy vectorization land */
  float *restrict x = c->x;
  float *restrict y = c->y;
  float *restrict z = c->z;
  float *restrict m = c->m;
  float *restrict epsilon = c->epsilon;
  swift_align_information(x, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(y, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(z, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(epsilon, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(m, SWIFT_CACHE_ALIGNMENT);
  swift_assume_size(gcount_padded, VEC_SIZE);

  /* Fill the input caches */
  for (int i = 0; i < gcount; ++i) {
    x[i] = (float)(gparts[i].x[0] - shift[0]);
    y[i] = (float)(gparts[i].x[1] - shift[1]);
    z[i] = (float)(gparts[i].x[2] - shift[2]);
    epsilon[i] = gparts[i].epsilon;
    m[i] = gparts[i].mass;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = 0.f;
    y[i] = 0.f;
    z[i] = 0.f;
    epsilon[i] = 0.f;
    m[i] = 0.f;
  }
}

/**
 * @brief Fills a #gravity_cache structure with some #gpart.
 *
 * @param c The #gravity_cache to fill.
 * @param gparts The #gpart array to read from.
 * @param gcount The number of particles to read.
 * @param gcount_padded The number of particle to read padded to the next
 * multiple of the vector length.
 */
__attribute__((always_inline)) INLINE void gravity_cache_populate_no_shift(
    struct gravity_cache *c, const struct gpart *restrict gparts, int gcount,
    int gcount_padded) {

  /* Make the compiler understand we are in happy vectorization land */
  float *restrict x = c->x;
  float *restrict y = c->y;
  float *restrict z = c->z;
  float *restrict m = c->m;
  float *restrict epsilon = c->epsilon;
  swift_align_information(x, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(y, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(z, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(epsilon, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(m, SWIFT_CACHE_ALIGNMENT);
  swift_assume_size(gcount_padded, VEC_SIZE);

  /* Fill the input caches */
  for (int i = 0; i < gcount; ++i) {
    x[i] = (float)(gparts[i].x[0]);
    y[i] = (float)(gparts[i].x[1]);
    z[i] = (float)(gparts[i].x[2]);
    epsilon[i] = gparts[i].epsilon;
    m[i] = gparts[i].mass;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = 0.f;
    y[i] = 0.f;
    z[i] = 0.f;
    epsilon[i] = 0.f;
    m[i] = 0.f;
  }
}

/**
 * @brief Write the output cache values back to the #gpart.
 *
 * @param c The #gravity_cache to read from.
 * @param gparts The #gpart array to write to.
 * @param gcount The number of particles to write.
 */
__attribute__((always_inline)) INLINE void gravity_cache_write_back(
    const struct gravity_cache *c, struct gpart *restrict gparts, int gcount) {

  /* Make the compiler understand we are in happy vectorization land */
  float *restrict a_x = c->a_x;
  float *restrict a_y = c->a_y;
  float *restrict a_z = c->a_z;
  swift_align_information(a_x, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(a_y, SWIFT_CACHE_ALIGNMENT);
  swift_align_information(a_z, SWIFT_CACHE_ALIGNMENT);

  /* Write stuff back to the particles */
  for (int i = 0; i < gcount; ++i) {
    gparts[i].a_grav[0] += a_x[i];
    gparts[i].a_grav[1] += a_y[i];
    gparts[i].a_grav[2] += a_z[i];
  }
}

#endif /* SWIFT_GRAVITY_CACHE_H */
