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

  /*! Is this #gpart active ? */
  int *restrict active SWIFT_CACHE_ALIGN;

  /*! Can this #gpart use a M2P interaction ? */
  int *restrict use_mpole SWIFT_CACHE_ALIGN;

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
    free(c->active);
    free(c->use_mpole);
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
  const size_t sizeBytesF = padded_count * sizeof(float);
  const size_t sizeBytesI = padded_count * sizeof(int);

  /* Delete old stuff if any */
  gravity_cache_clean(c);

  int e = 0;
  e += posix_memalign((void **)&c->x, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->y, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->z, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->epsilon, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->m, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->a_x, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->a_y, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->a_z, SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += posix_memalign((void **)&c->active, SWIFT_CACHE_ALIGNMENT, sizeBytesI);
  e +=
      posix_memalign((void **)&c->use_mpole, SWIFT_CACHE_ALIGNMENT, sizeBytesI);

  if (e != 0) error("Couldn't allocate gravity cache, size: %d", padded_count);

  c->count = padded_count;
}

/**
 * @brief Fills a #gravity_cache structure with some #gpart and shift them.
 *
 * Also checks whether the #gpart can use a M2P interaction instead of the
 * more expensive P2P.
 *
 * @param max_active_bin The largest active bin in the current time-step.
 * @param c The #gravity_cache to fill.
 * @param gparts The #gpart array to read from.
 * @param gcount The number of particles to read.
 * @param gcount_padded The number of particle to read padded to the next
 * multiple of the vector length.
 * @param shift A shift to apply to all the particles.
 * @param CoM The position of the multipole.
 * @param r_max2 The square of the multipole radius.
 * @param theta_crit2 The square of the opening angle.
 * @param cell The cell we play with (to get reasonable padding positions).
 */
__attribute__((always_inline)) INLINE static void gravity_cache_populate(
    timebin_t max_active_bin, struct gravity_cache *c,
    const struct gpart *restrict gparts, int gcount, int gcount_padded,
    const double shift[3], const float CoM[3], float r_max2, float theta_crit2,
    const struct cell *cell) {

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, x, c->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, c->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, c->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, epsilon, c->epsilon, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, c->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, active, c->active, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, use_mpole, c->use_mpole,
                            SWIFT_CACHE_ALIGNMENT);
  swift_assume_size(gcount_padded, VEC_SIZE);

  /* Fill the input caches */
  for (int i = 0; i < gcount; ++i) {
    x[i] = (float)(gparts[i].x[0] - shift[0]);
    y[i] = (float)(gparts[i].x[1] - shift[1]);
    z[i] = (float)(gparts[i].x[2] - shift[2]);
    epsilon[i] = gparts[i].epsilon;
    m[i] = gparts[i].mass;
    active[i] = (int)(gparts[i].time_bin <= max_active_bin);

    /* Check whether we can use the multipole instead of P-P */
    const float dx = x[i] - CoM[0];
    const float dy = y[i] - CoM[1];
    const float dz = z[i] - CoM[2];
    const float r2 = dx * dx + dy * dy + dz * dz;
    use_mpole[i] = gravity_M2P_accept(r_max2, theta_crit2, r2);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Particles used for padding should get impossible positions
   * that have a reasonable magnitude. We use the cell width for this */
  const float pos_padded[3] = {-2. * cell->width[0], -2. * cell->width[1],
                               -2. * cell->width[2]};

  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    epsilon[i] = 0.f;
    m[i] = 0.f;
    active[i] = 0;
    use_mpole[i] = 0;
  }
}

/**
 * @brief Fills a #gravity_cache structure with some #gpart and shift them.
 *
 * @param max_active_bin The largest active bin in the current time-step.
 * @param c The #gravity_cache to fill.
 * @param gparts The #gpart array to read from.
 * @param gcount The number of particles to read.
 * @param gcount_padded The number of particle to read padded to the next
 * multiple of the vector length.
 * @param shift A shift to apply to all the particles.
 * @param cell The cell we play with (to get reasonable padding positions).
 */
__attribute__((always_inline)) INLINE static void
gravity_cache_populate_no_mpole(timebin_t max_active_bin,
                                struct gravity_cache *c,
                                const struct gpart *restrict gparts, int gcount,
                                int gcount_padded, const double shift[3],
                                const struct cell *cell) {

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, x, c->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, c->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, c->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, epsilon, c->epsilon, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, c->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, active, c->active, SWIFT_CACHE_ALIGNMENT);
  swift_assume_size(gcount_padded, VEC_SIZE);

  /* Fill the input caches */
  for (int i = 0; i < gcount; ++i) {
    x[i] = (float)(gparts[i].x[0] - shift[0]);
    y[i] = (float)(gparts[i].x[1] - shift[1]);
    z[i] = (float)(gparts[i].x[2] - shift[2]);
    epsilon[i] = gparts[i].epsilon;
    m[i] = gparts[i].mass;
    active[i] = (int)(gparts[i].time_bin <= max_active_bin);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Particles used for padding should get impossible positions
   * that have a reasonable magnitude. We use the cell width for this */
  const float pos_padded[3] = {-2. * cell->width[0], -2. * cell->width[1],
                               -2. * cell->width[2]};
  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    epsilon[i] = 0.f;
    m[i] = 0.f;
    active[i] = 0;
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
  swift_declare_aligned_ptr(float, a_x, c->a_x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_y, c->a_y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_z, c->a_z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, active, c->active, SWIFT_CACHE_ALIGNMENT);

  /* Write stuff back to the particles */
  for (int i = 0; i < gcount; ++i) {
    if (active[i]) {
      gparts[i].a_grav[0] += a_x[i];
      gparts[i].a_grav[1] += a_y[i];
      gparts[i].a_grav[2] += a_z[i];
    }
  }
}

#endif /* SWIFT_GRAVITY_CACHE_H */
