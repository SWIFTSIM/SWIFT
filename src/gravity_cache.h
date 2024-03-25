/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* Local headers */
#include "accumulate.h"
#include "align.h"
#include "error.h"
#include "gravity.h"
#include "multipole_accept.h"
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

  /*! #gpart potential. */
  float *restrict pot SWIFT_CACHE_ALIGN;

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
    swift_free("gravity_cache", c->x);
    swift_free("gravity_cache", c->y);
    swift_free("gravity_cache", c->z);
    swift_free("gravity_cache", c->epsilon);
    swift_free("gravity_cache", c->m);
    swift_free("gravity_cache", c->a_x);
    swift_free("gravity_cache", c->a_y);
    swift_free("gravity_cache", c->a_z);
    swift_free("gravity_cache", c->pot);
    swift_free("gravity_cache", c->active);
    swift_free("gravity_cache", c->use_mpole);
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
static INLINE void gravity_cache_init(struct gravity_cache *c,
                                      const int count) {

  /* Size of the gravity cache */
  const int padded_count = count - (count % VEC_SIZE) + VEC_SIZE;
  const size_t sizeBytesF = padded_count * sizeof(float);
  const size_t sizeBytesI = padded_count * sizeof(int);

  /* Delete old stuff if any */
  gravity_cache_clean(c);

  int e = 0;
  e += swift_memalign("gravity_cache", (void **)&c->x, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->y, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->z, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->epsilon,
                      SWIFT_CACHE_ALIGNMENT, sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->m, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->a_x, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->a_y, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->a_z, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->pot, SWIFT_CACHE_ALIGNMENT,
                      sizeBytesF);
  e += swift_memalign("gravity_cache", (void **)&c->active,
                      SWIFT_CACHE_ALIGNMENT, sizeBytesI);
  e += swift_memalign("gravity_cache", (void **)&c->use_mpole,
                      SWIFT_CACHE_ALIGNMENT, sizeBytesI);

  if (e != 0) error("Couldn't allocate gravity cache, size: %d", padded_count);

  c->count = padded_count;
}

/**
 * @brief Zero all the output fields (acceleration and potential) of a
 * #gravity_cache.
 *
 * @param c The #gravity_cache to zero.
 * @param gcount_padded The padded size of the cache arrays.
 */
INLINE static void gravity_cache_zero_output(struct gravity_cache *c,
                                             const int gcount_padded) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded % VEC_SIZE != 0)
    error("Padded gcount size not a multiple of the vector length");
#endif

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, a_x, c->a_x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_y, c->a_y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_z, c->a_z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pot, c->pot, SWIFT_CACHE_ALIGNMENT);
  swift_assume_size(gcount_padded, VEC_SIZE);

  /* Zero everything */
  bzero(a_x, gcount_padded * sizeof(float));
  bzero(a_y, gcount_padded * sizeof(float));
  bzero(a_z, gcount_padded * sizeof(float));
  bzero(pot, gcount_padded * sizeof(float));
}

/**
 * @brief Fills a #gravity_cache structure with some #gpart and shift them.
 *
 * Also checks whether the #gpart can use a M2P interaction instead of the
 * more expensive P2P.
 *
 * @param max_active_bin The largest active bin in the current time-step.
 * @param allow_mpole Are we allowing the use of multipoles?
 * @param periodic Are we using periodic BCs ?
 * @param dim The size of the simulation volume along each dimension.
 * @param c The #gravity_cache to fill.
 * @param gparts The #gpart array to read from.
 * @param gcount The number of particles to read.
 * @param gcount_padded The number of particle to read padded to the next
 * multiple of the vector length.
 * @param shift A shift to apply to all the particles.
 * @param CoM The position of the multipole.
 * @param multipole The mulipole to check for.
 * @param cell The cell we play with (to get reasonable padding positions).
 * @param grav_props The global gravity properties.
 */
INLINE static void gravity_cache_populate(
    const timebin_t max_active_bin, const int allow_mpole, const int periodic,
    const float dim[3], struct gravity_cache *c,
    const struct gpart *restrict gparts, const int gcount,
    const int gcount_padded, const double shift[3], const float CoM[3],
    const struct gravity_tensors *multipole, const struct cell *cell,
    const struct gravity_props *grav_props) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Invalid padded cache size. Too small.");
  if (gcount_padded % VEC_SIZE != 0)
    error("Padded gravity cache size invalid. Not a multiple of SIMD length.");
#endif

  /* Do we need to grow the cache? */
  if (c->count < gcount_padded) gravity_cache_init(c, gcount_padded + VEC_SIZE);

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
#if !defined(SWIFT_DEBUG_CHECKS) && _OPENMP >= 201307
#pragma omp simd
#endif
  for (int i = 0; i < gcount; ++i) {

    x[i] = (float)(gparts[i].x[0] - shift[0]);
    y[i] = (float)(gparts[i].x[1] - shift[1]);
    z[i] = (float)(gparts[i].x[2] - shift[2]);
    epsilon[i] = gravity_get_softening(&gparts[i], grav_props);

#ifdef SWIFT_DEBUG_CHECKS
    if (gparts[i].time_bin == time_bin_not_created) {
      error("Found an extra gpart in the gravity cache");
    }
#endif

    /* Make a dummy particle out of the inhibted ones */
    if (gparts[i].time_bin == time_bin_inhibited) {
      m[i] = 0.f;
      active[i] = 0;
    } else {
      m[i] = gparts[i].mass;
      active[i] = (int)(gparts[i].time_bin <= max_active_bin);
    }

    /* Distance to the CoM of the other cell. */
    float dx = x[i] - CoM[0];
    float dy = y[i] - CoM[1];
    float dz = z[i] - CoM[2];

    /* Apply periodic BC */
    if (periodic) {
      dx = nearestf(dx, dim[0]);
      dy = nearestf(dy, dim[1]);
      dz = nearestf(dz, dim[2]);
    }
    const float r2 = dx * dx + dy * dy + dz * dz;

    /* Check whether we can use the multipole instead of P-P */
    use_mpole[i] = allow_mpole && gravity_M2P_accept(grav_props, &gparts[i],
                                                     multipole, r2, periodic);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Particles used for padding should get impossible positions
   * that have a reasonable magnitude. We use the cell width for this */
  const float pos_padded[3] = {-2.f * (float)cell->width[0],
                               -2.f * (float)cell->width[1],
                               -2.f * (float)cell->width[2]};
  const float eps_padded = epsilon[0];

  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    epsilon[i] = eps_padded;
    m[i] = 0.f;
    active[i] = 0;
    use_mpole[i] = 0;
  }

  /* Zero the output as well */
  gravity_cache_zero_output(c, gcount_padded);
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
 * @param grav_props The global gravity properties.
 */
INLINE static void gravity_cache_populate_no_mpole(
    const timebin_t max_active_bin, struct gravity_cache *c,
    const struct gpart *restrict gparts, const int gcount,
    const int gcount_padded, const double shift[3], const struct cell *cell,
    const struct gravity_props *grav_props) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Invalid padded cache size. Too small.");
  if (gcount_padded % VEC_SIZE != 0)
    error("Padded gravity cache size invalid. Not a multiple of SIMD length.");
#endif

  /* Do we need to grow the cache? */
  if (c->count < gcount_padded) gravity_cache_init(c, gcount_padded + VEC_SIZE);

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
    epsilon[i] = gravity_get_softening(&gparts[i], grav_props);

#ifdef SWIFT_DEBUG_CHECKS
    if (gparts[i].time_bin == time_bin_not_created) {
      error("Found an extra gpart in the gravity cache");
    }
#endif

    /* Make a dummy particle out of the inhibted ones */
    if (gparts[i].time_bin == time_bin_inhibited) {
      m[i] = 0.f;
      active[i] = 0;
    } else {
      m[i] = gparts[i].mass;
      active[i] = (int)(gparts[i].time_bin <= max_active_bin);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Particles used for padding should get impossible positions
   * that have a reasonable magnitude. We use the cell width for this */
  const float pos_padded[3] = {-2.f * (float)cell->width[0],
                               -2.f * (float)cell->width[1],
                               -2.f * (float)cell->width[2]};
  const float eps_padded = epsilon[0];

  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    epsilon[i] = eps_padded;
    m[i] = 0.f;
    active[i] = 0;
  }

  /* Zero the output as well */
  gravity_cache_zero_output(c, gcount_padded);
}

/**
 * @brief Fills a #gravity_cache structure with some #gpart and make them use
 * the multi-pole.
 *
 * @param max_active_bin The largest active bin in the current time-step.
 * @param periodic Are we using periodic BCs ?
 * @param dim The size of the simulation volume along each dimension.
 * @param c The #gravity_cache to fill.
 * @param gparts The #gpart array to read from.
 * @param gcount The number of particles to read.
 * @param gcount_padded The number of particle to read padded to the next
 * multiple of the vector length.
 * @param cell The cell we play with (to get reasonable padding positions).
 * @param CoM The position of the multipole.
 * @param multipole The mulipole to check for.
 * @param grav_props The global gravity properties.
 */
INLINE static void gravity_cache_populate_all_mpole(
    const timebin_t max_active_bin, const int periodic, const float dim[3],
    struct gravity_cache *c, const struct gpart *restrict gparts,
    const int gcount, const int gcount_padded, const struct cell *cell,
    const float CoM[3], const struct gravity_tensors *multipole,
    const struct gravity_props *grav_props) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Invalid padded cache size. Too small.");
  if (gcount_padded % VEC_SIZE != 0)
    error("Padded gravity cache size invalid. Not a multiple of SIMD length.");
#endif

  /* Do we need to grow the cache? */
  if (c->count < gcount_padded) gravity_cache_init(c, gcount_padded + VEC_SIZE);

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
    x[i] = (float)(gparts[i].x[0]);
    y[i] = (float)(gparts[i].x[1]);
    z[i] = (float)(gparts[i].x[2]);
    epsilon[i] = gravity_get_softening(&gparts[i], grav_props);
    m[i] = gparts[i].mass;
    active[i] = (int)(gparts[i].time_bin <= max_active_bin);
    use_mpole[i] = 1;

#ifdef SWIFT_DEBUG_CHECKS
    /* Distance to the CoM of the other cell. */
    float dx = x[i] - CoM[0];
    float dy = y[i] - CoM[1];
    float dz = z[i] - CoM[2];

    /* Apply periodic BC */
    if (periodic) {
      dx = nearestf(dx, dim[0]);
      dy = nearestf(dy, dim[1]);
      dz = nearestf(dz, dim[2]);
    }
    const float r2 = dx * dx + dy * dy + dz * dz;

    if (!gravity_M2P_accept(grav_props, &gparts[i], multipole, r2, periodic))
      error("Using m-pole where the test fails");
#endif
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (gcount_padded < gcount) error("Padded counter smaller than counter");
#endif

  /* Particles used for padding should get impossible positions
   * that have a reasonable magnitude. We use the cell width for this */
  const float pos_padded[3] = {-2.f * (float)cell->width[0],
                               -2.f * (float)cell->width[1],
                               -2.f * (float)cell->width[2]};
  const float eps_padded = epsilon[0];

  /* Pad the caches */
  for (int i = gcount; i < gcount_padded; ++i) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    epsilon[i] = eps_padded;
    m[i] = 0.f;
    active[i] = 0;
    use_mpole[i] = 0;
  }

  /* Zero the output as well */
  gravity_cache_zero_output(c, gcount_padded);
}

/**
 * @brief Write the output cache values back to the active #gpart.
 *
 * This function obviously omits the padded values in the cache.
 *
 * @param c The #gravity_cache to read from.
 * @param gparts The #gpart array to write to.
 * @param gcount The number of particles to write.
 */
INLINE static void gravity_cache_write_back(const struct gravity_cache *c,
                                            struct gpart *restrict gparts,
                                            const int gcount) {

  /* Make the compiler understand we are in happy vectorization land */
  swift_declare_aligned_ptr(float, a_x, c->a_x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_y, c->a_y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, a_z, c->a_z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pot, c->pot, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(int, active, c->active, SWIFT_CACHE_ALIGNMENT);

  /* Write stuff back to the particles */
#if !defined(SWIFT_DEBUG_CHECKS) && _OPENMP >= 201307
#pragma omp simd
#endif
  for (int i = 0; i < gcount; ++i) {
    if (active[i]) {
      gparts[i].a_grav[0] += a_x[i];
      gparts[i].a_grav[1] += a_y[i];
      gparts[i].a_grav[2] += a_z[i];
      gravity_add_comoving_potential(&gparts[i], pot[i]);
    }
  }
}

#endif /* SWIFT_GRAVITY_CACHE_H */
