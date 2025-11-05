/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_GHOST_STATS_H
#define SWIFT_GHOST_STATS_H

/* Config parameters. */
#include "minmax.h"
#include "part.h"

#include <config.h>
#include <float.h>
#include <stdio.h>

#ifdef SWIFT_GHOST_STATS

/*****************************
 * ghost_stats functionality *
 *****************************/

/**
 * @brief Entry for a single iteration in the ghost statistics table.
 */
struct ghost_stats_entry {
  /*! Number of particles processed during the iteration. */
  int count;
  /*! Number of particles in the iteration that had no neighbours. */
  int count_no_ngb;
  /*! Minimum initial smoothing length of the particles at the start of the
   *  iteration. */
  float hmin;
  /*! Maximum initial smoothing length of the particles at the start of the
   *  iteration. */
  float hmax;
  /*! Sum of the initial smoothing lengths, useful to compute averages in
   *  post-processing. */
  double hsum;
  /*! Sum of the initial smoothing lengths squared, useful to compute variances
   *  and standard deviations in post-processing. */
  double hsum2;
};

/**
 * @brief Ghost statistics stored in a cell.
 */
struct ghost_stats {
  /* Hydro ghost statistics. */
  struct ghost_stats_entry hydro[SWIFT_GHOST_STATS + 1];
  /* Stars ghost statistics. */
  struct ghost_stats_entry stars[SWIFT_GHOST_STATS + 1];
  /* Black holes ghost statistics. */
  struct ghost_stats_entry black_holes[SWIFT_GHOST_STATS + 1];
};

/* ghost_stats_entry struct functions */

/**
 * @brief Reset the given ghost stats bin.
 *
 * @param bin Ghost stats bin to reset.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_reset_entry(
    struct ghost_stats_entry *restrict bin) {

  bin->count = 0;
  bin->count_no_ngb = 0;
  bin->hmin = FLT_MAX;
  bin->hmax = 0.0f;
  bin->hsum = 0.;
  bin->hsum2 = 0.;
}

/**
 * @brief Write the given ghost stats bin to the given file.
 *
 * @param f File to write to.
 * @param bin Ghost stats bin to write.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_print_entry(
    FILE *f, const struct ghost_stats_entry *restrict bin) {

  if (bin->count > 0) {
    fprintf(f, "\t%i\t%i\t%g\t%g\t%g\t%g", bin->count, bin->count_no_ngb,
            bin->hmin, bin->hmax, bin->hsum, bin->hsum2);
  } else {
    fprintf(f, "\t0\t0\t0\t0\t0\t0");
  }
}

/* ghost_stats struct functions */

/**
 * @brief Reset all the entries in the ghost_stats struct.
 *
 * @param gstats Ghost stats struct.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_reset_entries(
    struct ghost_stats *restrict gstats) {

  for (int b = 0; b < SWIFT_GHOST_STATS + 1; ++b) {
    ghost_stats_reset_entry(&gstats->hydro[b]);
  }
  for (int b = 0; b < SWIFT_GHOST_STATS + 1; ++b) {
    ghost_stats_reset_entry(&gstats->stars[b]);
  }
  for (int b = 0; b < SWIFT_GHOST_STATS + 1; ++b) {
    ghost_stats_reset_entry(&gstats->black_holes[b]);
  }
}

/**
 * @brief Account for the star particles that are still under consideration at
 * the start of a ghost decision loop (so after the neighbour loop but before
 * the smoothing length is updated).
 *
 * @param gstats Ghost stats struct to update.
 * @param iteration_number Iteration number in the high-level ghost iteration
 * scheme.
 * @param scount Number of star particles still under consideration (smoothing
 * length has not converged yet or will do so during this iteration).
 * @param sparts Star particle array.
 * @param sid Indices of active star particles in the array.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_account_for_stars(
    struct ghost_stats *restrict gstats, int iteration_number, int scount,
    struct spart *restrict sparts, int *sid) {

  /* accumulate in highest bin, so we can spot out-of-bounds values */
  int binidx = min(iteration_number, SWIFT_GHOST_STATS - 1);
  struct ghost_stats_entry *restrict sbin = &gstats->stars[binidx];
  sbin->count += scount;
  for (int i = 0; i < scount; i++) {
    const float hi = sparts[sid[i]].h;
    sbin->hmin = min(sbin->hmin, hi);
    sbin->hmax = max(sbin->hmax, hi);
    sbin->hsum += hi;
    sbin->hsum2 += hi * hi;
  }
}

/**
 * @brief Account for the properties of the a converged star particle.
 *
 * @param gstats Ghost stats struct to update.
 * @param sp Star particle that has a converged smoothing length.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_converged_star(
    struct ghost_stats *restrict gstats, struct spart *restrict sp) {

  struct ghost_stats_entry *restrict sbin = &gstats->stars[SWIFT_GHOST_STATS];
  ++sbin->count;
  const float hi = sp->h;
  sbin->hmin = min(sbin->hmin, hi);
  sbin->hmax = max(sbin->hmax, hi);
  sbin->hsum += hi;
  sbin->hsum2 += hi * hi;
}

/**
 * @brief Register the occurrence of a "no neighbour" event during the current
 * star iteration step.
 *
 * @param gstats Ghost stats struct to update.
 * @param iteration_number Number of the current iteration.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_star_iteration(struct ghost_stats *restrict gstats,
                                  int iteration_number) {

  int binidx = min(iteration_number, SWIFT_GHOST_STATS - 1);
  struct ghost_stats_entry *restrict sbin = &gstats->stars[binidx];
  ++sbin->count_no_ngb;
}

/**
 * @brief Register the occurrence of a converged star particle without
 * neighbours.
 *
 * @param gstats Ghost stats struct to update.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_star_converged(struct ghost_stats *restrict gstats) {

  struct ghost_stats_entry *restrict sbin = &gstats->stars[SWIFT_GHOST_STATS];
  ++sbin->count_no_ngb;
}

/**
 * @brief Account for the black hole particles that are still under
 * consideration at the start of a ghost decision loop (so after the neighbour
 * loop but before the smoothing length is updated).
 *
 * @param gstats Ghost stats struct to update.
 * @param iteration_number Iteration number in the high-level ghost iteration
 * scheme.
 * @param bcount Number of black hole particles still under consideration
 * (smoothing length has not converged yet or will do so during this iteration).
 * @param bparts Black hole particle array.
 * @param sid Indices of active black hole particles in the array.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_account_for_black_holes(struct ghost_stats *restrict gstats,
                                    int iteration_number, int bcount,
                                    struct bpart *restrict bparts, int *sid) {

  /* accumulate in highest bin, so we can spot out-of-bounds values */
  int binidx = min(iteration_number, SWIFT_GHOST_STATS - 1);
  struct ghost_stats_entry *restrict bbin = &gstats->black_holes[binidx];
  bbin->count += bcount;
  for (int i = 0; i < bcount; i++) {
    const float hi = bparts[sid[i]].h;
    bbin->hmin = min(bbin->hmin, hi);
    bbin->hmax = max(bbin->hmax, hi);
    bbin->hsum += hi;
    bbin->hsum2 += hi * hi;
  }
}

/**
 * @brief Account for the properties of the a converged black hole particle.
 *
 * @param gstats Ghost stats struct to update.
 * @param bp Black hole particle that has a converged smoothing length.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_converged_black_hole(struct ghost_stats *restrict gstats,
                                 struct bpart *restrict bp) {

  struct ghost_stats_entry *restrict bbin =
      &gstats->black_holes[SWIFT_GHOST_STATS];
  ++bbin->count;
  const float hi = bp->h;
  bbin->hmin = min(bbin->hmin, hi);
  bbin->hmax = max(bbin->hmax, hi);
  bbin->hsum += hi;
  bbin->hsum2 += hi * hi;
}

/**
 * @brief Register the occurrence of a "no neighbour" event during the current
 * black hole iteration step.
 *
 * @param gstats Ghost stats struct to update.
 * @param iteration_number Number of the current iteration.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_black_hole_iteration(struct ghost_stats *restrict gstats,
                                        int iteration_number) {

  int binidx = min(iteration_number, SWIFT_GHOST_STATS - 1);
  struct ghost_stats_entry *restrict bbin = &gstats->black_holes[binidx];
  ++bbin->count_no_ngb;
}

/**
 * @brief Register the occurrence of a converged black hole particle without
 * neighbours.
 *
 * @param gstats Ghost stats struct to update.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_black_hole_converged(struct ghost_stats *restrict gstats) {

  struct ghost_stats_entry *restrict bbin =
      &gstats->black_holes[SWIFT_GHOST_STATS];
  ++bbin->count_no_ngb;
}

/**
 * @brief Account for the gas particles that are still under consideration at
 * the start of a ghost decision loop (so after the neighbour loop but before
 * the smoothing length is updated).
 *
 * @param gstats Ghost stats struct to update.
 * @param iteration_number Iteration number in the high-level ghost iteration
 * scheme.
 * @param count Number of gas particles still under consideration (smoothing
 * length has not converged yet or will do so during this iteration).
 * @param parts Gas particle array.
 * @param sid Indices of active gas particles in the array.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_account_for_hydro(
    struct ghost_stats *restrict gstats, int iteration_number, int count,
    struct part *restrict parts, int *pid) {

  /* accumulate in highest bin, so we can spot out-of-bounds values */
  int binidx = min(iteration_number, SWIFT_GHOST_STATS - 1);
  struct ghost_stats_entry *restrict hbin = &gstats->hydro[binidx];
  hbin->count += count;
  for (int i = 0; i < count; i++) {
    const float hi = parts[pid[i]].h;
    hbin->hmin = min(hbin->hmin, hi);
    hbin->hmax = max(hbin->hmax, hi);
    hbin->hsum += hi;
    hbin->hsum2 += hi * hi;
  }
}

/**
 * @brief Account for the properties of the a converged gas particle.
 *
 * @param gstats Ghost stats struct to update.
 * @param p Gas particle that has a converged smoothing length.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_converged_hydro(
    struct ghost_stats *restrict gstats, struct part *restrict p) {

  struct ghost_stats_entry *restrict hbin = &gstats->hydro[SWIFT_GHOST_STATS];
  ++hbin->count;
  const float hi = p->h;
  hbin->hmin = min(hbin->hmin, hi);
  hbin->hmax = max(hbin->hmax, hi);
  hbin->hsum += hi;
  hbin->hsum2 += hi * hi;
}

/**
 * @brief Register the occurrence of a "no neighbour" event during the current
 * hydro iteration step.
 *
 * @param gstats Ghost stats struct to update.
 * @param iteration_number Number of the current iteration.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_hydro_iteration(struct ghost_stats *restrict gstats,
                                   int iteration_number) {

  int binidx = min(iteration_number, SWIFT_GHOST_STATS - 1);
  struct ghost_stats_entry *restrict hbin = &gstats->hydro[binidx];
  ++hbin->count_no_ngb;
}

/**
 * @brief Register the occurrence of a converged gas particle without
 * neighbours.
 *
 * @param gstats Ghost stats struct to update.
 */
__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_hydro_converged(struct ghost_stats *restrict gstats) {

  struct ghost_stats_entry *restrict hbin = &gstats->hydro[SWIFT_GHOST_STATS];
  ++hbin->count_no_ngb;
}

/**
 * @brief Write the header of a ghost statistics file.
 *
 * @param f File to write to.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_write_header(
    FILE *f) {

  fprintf(f, "# Ghost statistics\n");
  fprintf(f, "# Values listed in blocks per particle type\n");
  fprintf(f, "# Order of types: hydro, stars, black holes\n");
  fprintf(f, "# Number of blocks per type: %i\n", SWIFT_GHOST_STATS + 1);
  fprintf(f, "# Number of values per block: 6\n");
  fprintf(f, "# Last block contains converged values\n");
  fprintf(f, "# Fields per block:\n");
  fprintf(f, "#  - count: i4\n");
  fprintf(f, "#  - no neighbour count: i4\n");
  fprintf(f, "#  - min h: f4\n");
  fprintf(f, "#  - max h: f4\n");
  fprintf(f, "#  - sum h: f8\n");
  fprintf(f, "#  - sum h^2: f8\n");
  fprintf(f, "# First column is cellID\n");
  fprintf(f, "# Cells with no values are omitted\n");
}

/**
 * @brief Write the ghost statistics for the given cell to the given file.
 *
 * @param f File to write to.
 * @param gstats Ghost statistics to write.
 * @param cellID Cell ID to use to identify the cell.
 */
__attribute__((always_inline)) INLINE static void ghost_stats_write_cell_stats(
    FILE *f, const struct ghost_stats *restrict gstats,
    const long long cellID) {

  if (gstats->hydro[0].count + gstats->stars[0].count > 0) {
    fprintf(f, "%lld", cellID);
    for (int b = 0; b < SWIFT_GHOST_STATS + 1; ++b) {
      ghost_stats_print_entry(f, &gstats->hydro[b]);
    }
    for (int b = 0; b < SWIFT_GHOST_STATS + 1; ++b) {
      ghost_stats_print_entry(f, &gstats->stars[b]);
    }
    for (int b = 0; b < SWIFT_GHOST_STATS + 1; ++b) {
      ghost_stats_print_entry(f, &gstats->black_holes[b]);
    }
    fprintf(f, "\n");
  }
}

/******************
 * cell interface *
 ******************/

struct cell;

void cell_reset_ghost_histograms(struct cell *c);
void cell_write_ghost_stats(FILE *f, const struct cell *c,
                            const long long cellID);
/*******************
 * space interface *
 *******************/

struct space;

void space_reset_ghost_histograms(struct space *s);
void space_write_ghost_stats(const struct space *s, int j);

#else

struct ghost_stats {};

/* stars */
__attribute__((always_inline)) INLINE static void ghost_stats_account_for_stars(
    struct ghost_stats *restrict gstats, int iteration_number, int scount,
    struct spart *restrict sparts, int *sid) {}

__attribute__((always_inline)) INLINE static void ghost_stats_converged_star(
    struct ghost_stats *restrict gstats, struct spart *restrict sp) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_star_iteration(struct ghost_stats *restrict gstats,
                                  int iteration_number) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_star_converged(struct ghost_stats *restrict gstats) {}

/* black holes */
__attribute__((always_inline)) INLINE static void
ghost_stats_account_for_black_holes(struct ghost_stats *restrict gstats,
                                    int iteration_number, int bcount,
                                    struct bpart *restrict bparts, int *sid) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_converged_black_hole(struct ghost_stats *restrict gstats,
                                 struct bpart *restrict bp) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_black_hole_iteration(struct ghost_stats *restrict gstats,
                                        int iteration_number) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_black_hole_converged(struct ghost_stats *restrict gstats) {}

/* hydro */
__attribute__((always_inline)) INLINE static void ghost_stats_account_for_hydro(
    struct ghost_stats *restrict gstats, int iteration_number, int count,
    struct part *restrict parts, int *pid) {}

__attribute__((always_inline)) INLINE static void ghost_stats_converged_hydro(
    struct ghost_stats *restrict gstats, struct part *restrict p) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_hydro_iteration(struct ghost_stats *restrict gstats,
                                   int iteration_number) {}

__attribute__((always_inline)) INLINE static void
ghost_stats_no_ngb_hydro_converged(struct ghost_stats *restrict gstats) {}

/// cell interface

struct cell;

__attribute__((always_inline)) INLINE static void cell_reset_ghost_histograms(
    struct cell *c) {}

/// space interface

struct space;

__attribute__((always_inline)) INLINE static void space_reset_ghost_histograms(
    struct space *s) {}

__attribute__((always_inline)) INLINE static void space_write_ghost_stats(
    const struct space *s, int j) {}
#endif

#endif /* SWIFT_GHOST_STATS */
