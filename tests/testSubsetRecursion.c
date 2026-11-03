/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/*
 * Checks the density loops that the ghosts rerun on a subset of particles
 * while their smoothing length converges. The cells keep the h_max_old of the
 * previous step, while the particle being redone already has a larger h. The
 * density sums of the subset loops must equal a brute-force sum over all the
 * gas of the two cells.
 *
 * Two adjacent top-level cells of width 1 are split once into octants of
 * width 0.5. The gas sits on a lattice. The particle being redone sits just
 * inside the face shared by the two cells.
 */

void runner_dosub_pair_subset_stars_density(struct runner *r, struct cell *ci,
                                            struct spart *sparts,
                                            const int *ind, const int scount,
                                            struct cell *cj,
                                            const int gettimer);
void runner_dosub_self_subset_stars_density(struct runner *r, struct cell *ci,
                                            struct spart *sparts,
                                            const int *ind, const int scount,
                                            const int gettimer);
void runner_dosub_pair_subset_density(struct runner *r, struct cell *ci,
                                      struct part *parts, const int *ind,
                                      const int count, struct cell *cj,
                                      const int gettimer);
void runner_dosub_self_subset_density(struct runner *r, struct cell *ci,
                                      struct part *parts, const int *ind,
                                      const int count, const int gettimer);
void runner_dosub_subset_bh_density(struct runner *r, struct cell *ci,
                                    struct bpart *bparts, int *ind,
                                    const int bcount, struct cell *cj,
                                    int gettimer);
void runner_dosub_subset_sinks_density(struct runner *r, struct cell *ci,
                                       struct sink *sinks, int *ind,
                                       const int scount, struct cell *cj,
                                       int gettimer);

#if !defined(BLACK_HOLES_NONE)
#define WITH_BH_LOOP 1
#endif
#if defined(SINK_GEAR)
#define WITH_SINK_LOOP 1
#endif

#define TI_CURRENT 8

/* Gas lattice points per axis in one octant, coarse to keep the sums small. */
#define N_SIDE 2
#define N_OCTANT (N_SIDE * N_SIDE * N_SIDE)

/* Smoothing length of the lattice gas, small enough to allow recursion. */
#define H_LATTICE 0.05f

/* Position of the particle being redone, just inside the shared face. */
static const double x_redo[3] = {0.999, 0.25, 0.25};

/* Gas planted just beyond the octant next to the shared face. */
#define N_PLANT 4
static const double x_plant[N_PLANT][3] = {{1.505, 0.25, 0.25},
                                           {1.505, 0.30, 0.25},
                                           {1.505, 0.25, 0.30},
                                           {1.510, 0.20, 0.20}};

/**
 * @brief Build a top-level cell of width 1 split into eight leaf octants.
 *
 * Octant 4 of the cell at the origin also holds the particle being redone,
 * as its first gas particle, its only star, black hole and sink. Octant 4 of
 * the other cell also holds the planted gas.
 *
 * @param c The #cell to fill.
 * @param progeny The array of 8 #cell used as octants.
 * @param loc_x The x position of the cell's lower corner.
 * @param with_redo Whether this cell holds the particle being redone.
 */
static void make_split_cell(struct cell *c, struct cell *progeny,
                            const double loc_x, const int with_redo) {

  const int count = 8 * N_OCTANT + (with_redo ? 1 : N_PLANT);
  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0)
    error("Couldn't allocate the gas");
  bzero(c->hydro.parts, count * sizeof(struct part));
  if (posix_memalign((void **)&c->hydro.xparts, xpart_align,
                     count * sizeof(struct xpart)) != 0)
    error("Couldn't allocate the gas extras");
  bzero(c->hydro.xparts, count * sizeof(struct xpart));
  if (posix_memalign((void **)&c->stars.parts, spart_align,
                     sizeof(struct spart)) != 0)
    error("Couldn't allocate the stars");
  bzero(c->stars.parts, sizeof(struct spart));
  /* Without the module these cell fields share a union: leave them unset. */
#ifdef WITH_BH_LOOP
  if (posix_memalign((void **)&c->black_holes.parts, bpart_align,
                     sizeof(struct bpart)) != 0)
    error("Couldn't allocate the black holes");
  bzero(c->black_holes.parts, sizeof(struct bpart));
#endif
#ifdef WITH_SINK_LOOP
  if (posix_memalign((void **)&c->sinks.parts, sink_align,
                     sizeof(struct sink)) != 0)
    error("Couldn't allocate the sinks");
  bzero(c->sinks.parts, sizeof(struct sink));
#endif

  int offset = 0;
  for (int k = 0; k < 8; k++) {
    struct cell *cp = &progeny[k];
    bzero(cp, sizeof(struct cell));
    cp->loc[0] = loc_x + 0.5 * ((k >> 2) & 1);
    cp->loc[1] = 0.5 * ((k >> 1) & 1);
    cp->loc[2] = 0.5 * (k & 1);
    cp->hydro.parts = &c->hydro.parts[offset];
    cp->hydro.xparts = &c->hydro.xparts[offset];
    cp->stars.parts = c->stars.parts;
#ifdef WITH_BH_LOOP
    cp->black_holes.parts = c->black_holes.parts;
#endif
#ifdef WITH_SINK_LOOP
    cp->sinks.parts = c->sinks.parts;
#endif

    int n = 0;
    if (with_redo && k == 4) {
      struct part *p = &c->hydro.parts[offset];
      for (int d = 0; d < 3; d++) p->x[d] = x_redo[d];
      p->id = 1;
      n = 1;
      cp->stars.count = 1;
      struct spart *sp = c->stars.parts;
      for (int d = 0; d < 3; d++) sp->x[d] = x_redo[d];
      sp->id = 1;
      sp->time_bin = 1;
#if defined(FEEDBACK_GEAR) || defined(FEEDBACK_AGORA)
      sp->feedback_data.will_do_feedback = 1;
#endif
#ifdef WITH_BH_LOOP
      cp->black_holes.count = 1;
      struct bpart *bp = c->black_holes.parts;
      for (int d = 0; d < 3; d++) bp->x[d] = x_redo[d];
      bp->id = 1;
      bp->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
      bp->ti_drift = TI_CURRENT;
      bp->ti_kick = TI_CURRENT;
#endif
#endif
#ifdef WITH_SINK_LOOP
      cp->sinks.count = 1;
      struct sink *si = c->sinks.parts;
      for (int d = 0; d < 3; d++) si->x[d] = x_redo[d];
      si->id = 1;
      si->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
      si->ti_drift = TI_CURRENT;
      si->ti_kick = TI_CURRENT;
#endif
#endif
#ifdef SWIFT_DEBUG_CHECKS
      sp->ti_drift = TI_CURRENT;
      sp->ti_kick = TI_CURRENT;
#endif
    } else if (k == 4) {
      for (n = 0; n < N_PLANT; n++) {
        struct part *p = &c->hydro.parts[offset + n];
        for (int d = 0; d < 3; d++) p->x[d] = x_plant[n][d];
        p->id = 200000 + n;
      }
    }
    for (int ix = 0; ix < N_SIDE; ix++) {
      for (int iy = 0; iy < N_SIDE; iy++) {
        for (int iz = 0; iz < N_SIDE; iz++) {
          struct part *p = &c->hydro.parts[offset + n];
          p->x[0] = cp->loc[0] + (ix + 0.5) * 0.5 / N_SIDE;
          p->x[1] = cp->loc[1] + (iy + 0.5) * 0.5 / N_SIDE;
          p->x[2] = cp->loc[2] + (iz + 0.5) * 0.5 / N_SIDE;
          p->id = 2 + offset + n + (with_redo ? 0 : 100000);
          n++;
        }
      }
    }
    for (int i = 0; i < n; i++) {
      struct part *p = &c->hydro.parts[offset + i];
      p->h = H_LATTICE;
      p->mass = 1.f;
      p->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
      p->ti_drift = TI_CURRENT;
      p->ti_kick = TI_CURRENT;
#endif
    }
    cp->hydro.count = n;
    offset += n;
  }

  bzero(c, sizeof(struct cell));
  c->hydro.parts = progeny[0].hydro.parts;
  c->hydro.xparts = progeny[0].hydro.xparts;
  c->stars.parts = progeny[0].stars.parts;
#ifdef WITH_BH_LOOP
  c->black_holes.parts = progeny[0].black_holes.parts;
  c->black_holes.count = with_redo ? 1 : 0;
#endif
#ifdef WITH_SINK_LOOP
  c->sinks.parts = progeny[0].sinks.parts;
  c->sinks.count = with_redo ? 1 : 0;
#endif
  c->hydro.count = count;
  c->stars.count = with_redo ? 1 : 0;
  c->loc[0] = loc_x;
  c->split = 1;

  for (int k = 0; k < 9; k++) {
    struct cell *cc = (k < 8) ? &progeny[k] : c;
    const double width = (k < 8) ? 0.5 : 1.;
    for (int d = 0; d < 3; d++) cc->width[d] = width;
    cc->dmin = width;
    cc->hydro.h_max_old = H_LATTICE;
    cc->hydro.ti_end_min = TI_CURRENT;
    cc->hydro.ti_old_part = TI_CURRENT;
    cc->stars.ti_end_min = TI_CURRENT;
    cc->stars.ti_old_part = TI_CURRENT;
#ifdef WITH_BH_LOOP
    cc->black_holes.ti_end_min = TI_CURRENT;
    cc->black_holes.ti_old_part = TI_CURRENT;
#endif
#ifdef WITH_SINK_LOOP
    cc->sinks.ti_end_min = TI_CURRENT;
    cc->sinks.ti_old_part = TI_CURRENT;
#endif
    cc->nodeID = 0;
    if (k < 8) {
      cc->parent = c;
      cc->top = c;
      cc->depth = 1;
      c->progeny[k] = cc;
    } else {
      cc->top = c;
    }
  }
}

/**
 * @brief Brute-force density sums of a particle at x_redo over the gas of
 * both cells.
 *
 * @param cells The two top-level #cell.
 * @param h The smoothing length of the particle.
 * @param skip_id Gas id to leave out (the particle itself), or -1.
 * @param far_only Whether to count only the gas in the far half of the
 * second cell.
 * @param wcount (return) The kernel-weighted neighbour count.
 * @return The number of neighbours.
 */
static int brute_force(struct cell *cells, const float h,
                       const long long skip_id, const int far_only,
                       double *wcount) {
  const float h_inv = 1.f / h;
  int n_ngb = 0;
  *wcount = 0.;
  for (int c = 0; c < 2; c++) {
    if (far_only && c == 0) continue;
    for (int i = 0; i < cells[c].hydro.count; i++) {
      const struct part *p = &cells[c].hydro.parts[i];
      if (p->id == skip_id) continue;
      if (far_only && p->x[0] < 1.5) continue;
      const float dx[3] = {(float)(x_redo[0] - p->x[0]),
                           (float)(x_redo[1] - p->x[1]),
                           (float)(x_redo[2] - p->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (r2 < h * h * kernel_gamma2) {
        float w, w_dx;
        kernel_deval(sqrtf(r2) * h_inv, &w, &w_dx);
        *wcount += w;
        n_ngb++;
      }
    }
  }
  return n_ngb;
}

int main(int argc, char *argv[]) {

  clocks_set_cpufreq(0);

  struct space space;
  bzero(&space, sizeof(struct space));
  space.periodic = 0;
  space.dim[0] = space.dim[1] = space.dim[2] = 10.;

  static struct hydro_props hp;
  static struct stars_props stars_p;
  static struct feedback_props feedback_properties;
  static struct phys_const physical_constants;
  static struct unit_system internal_units;
  static struct cooling_function_data cooling_function;
  static struct sink_props sink_properties;
  static struct chemistry_global_data chemistry_data;
  static struct gravity_props gravity_properties;
  static struct black_holes_props black_holes_properties;
  static struct entropy_floor_properties entropy_floor;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.feedback_props = &feedback_properties;
  engine.physical_constants = &physical_constants;
  engine.internal_units = &internal_units;
  engine.cooling_func = &cooling_function;
  engine.sink_properties = &sink_properties;
  engine.chemistry = &chemistry_data;
  engine.gravity_properties = &gravity_properties;
  engine.black_holes_properties = &black_holes_properties;
  engine.entropy_floor = &entropy_floor;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = TI_CURRENT;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;
  engine.stars_properties = &stars_p;
  engine.nodeID = 0;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  struct cell cells[2];
  struct cell progeny[2][8];
  make_split_cell(&cells[0], progeny[0], 0., /*with_redo=*/1);
  make_split_cell(&cells[1], progeny[1], 1., /*with_redo=*/0);
  space.cells_top = cells;

  struct part *p_redo = progeny[0][4].hydro.parts;
  struct spart *sp_redo = cells[0].stars.parts;
  const int ind_gas = (int)(p_redo - cells[0].hydro.parts);
  const int ind_star = 0;
  int ind_one = 0;
#ifdef WITH_BH_LOOP
  struct bpart *bp_redo = cells[0].black_holes.parts;
#endif
#ifdef WITH_SINK_LOOP
  struct sink *si_redo = cells[0].sinks.parts;
#endif
  (void)ind_one;
  static const char *loop_name[4] = {"stars", "hydro", "black-hole", "sink"};

  /* kernel_gamma * h in units of the octant width, before and after the
   * update. Only values below 1 allow recursion. The last case shrinks h. */
  const float r_old = 0.9f;
  const float r_new[3] = {1.1f, 1.8f, 0.8f};

  int n_fail = 0;
  for (int k = 0; k < 3; k++) {
    const float h_old = r_old * 0.5f / kernel_gamma;
    const float h_new = r_new[k] * 0.5f / kernel_gamma;

    for (int loop = 0; loop < 4; loop++) {
      const int is_stars = (loop == 0);
#ifndef WITH_BH_LOOP
      if (loop == 2) continue;
#endif
#ifndef WITH_SINK_LOOP
      if (loop == 3) continue;
#endif

      for (int c = 0; c < 2; c++) {
        for (int j = 0; j < 9; j++) {
          struct cell *cc = (j < 8) ? &progeny[c][j] : &cells[c];
          cc->stars.h_max_old = (c == 0) ? h_old : 0.f;
#ifdef WITH_BH_LOOP
          cc->black_holes.h_max_old = (c == 0) ? h_old : 0.f;
#endif
#ifdef WITH_SINK_LOOP
          cc->sinks.h_max_old = (c == 0) ? h_old : 0.f;
#endif
          cc->hydro.h_max_old = (c == 0) ? max(h_old, H_LATTICE) : H_LATTICE;
        }
      }

      float wcount_code;
      int n_code = -1;
      if (is_stars) {
        bzero(&sp_redo->density, sizeof(sp_redo->density));
#if defined(FEEDBACK_GEAR)
        bzero(&sp_redo->feedback_data, sizeof(sp_redo->feedback_data));
        sp_redo->feedback_data.will_do_feedback = 1;
#endif
        sp_redo->h = h_new;
        runner_dosub_self_subset_stars_density(&runner, &cells[0],
                                               cells[0].stars.parts, &ind_star,
                                               1, /*gettimer=*/0);
        runner_dosub_pair_subset_stars_density(
            &runner, &cells[0], cells[0].stars.parts, &ind_star, 1, &cells[1],
            /*gettimer=*/0);
        wcount_code = sp_redo->density.wcount;
#if defined(FEEDBACK_GEAR)
        n_code = sp_redo->feedback_data.num_ngbs;
#endif
      } else if (loop == 2) {
#ifdef WITH_BH_LOOP
        bzero(&bp_redo->density, sizeof(bp_redo->density));
        bp_redo->h = h_new;
#if defined(BLACK_HOLES_GEAR)
        bp_redo->num_ngbs = 0;
#endif
        runner_dosub_subset_bh_density(&runner, &cells[0], bp_redo, &ind_one, 1,
                                       NULL, /*gettimer=*/0);
        runner_dosub_subset_bh_density(&runner, &cells[0], bp_redo, &ind_one, 1,
                                       &cells[1], /*gettimer=*/0);
        wcount_code = bp_redo->density.wcount;
#if defined(BLACK_HOLES_GEAR)
        n_code = bp_redo->num_ngbs;
#endif
#else
        wcount_code = 0.f;
#endif
      } else if (loop == 3) {
#ifdef WITH_SINK_LOOP
        bzero(&si_redo->density, sizeof(si_redo->density));
        si_redo->h = h_new;
        si_redo->num_ngbs = 0;
        runner_dosub_subset_sinks_density(&runner, &cells[0], si_redo, &ind_one,
                                          1, NULL, /*gettimer=*/0);
        runner_dosub_subset_sinks_density(&runner, &cells[0], si_redo, &ind_one,
                                          1, &cells[1], /*gettimer=*/0);
        wcount_code = si_redo->density.wcount;
        n_code = si_redo->num_ngbs;
#else
        wcount_code = 0.f;
#endif
      } else {
        bzero(&p_redo->density, sizeof(p_redo->density));
        p_redo->h = h_new;
        runner_dosub_self_subset_density(&runner, &cells[0],
                                         cells[0].hydro.parts, &ind_gas, 1,
                                         /*gettimer=*/0);
        runner_dosub_pair_subset_density(&runner, &cells[0],
                                         cells[0].hydro.parts, &ind_gas, 1,
                                         &cells[1], /*gettimer=*/0);
        wcount_code = p_redo->density.wcount;
        p_redo->h = H_LATTICE;
      }

      const long long skip_id = (loop == 1) ? p_redo->id : -1;
      double wcount_ref, wcount_far;
      const int n_ref = brute_force(cells, h_new, skip_id, 0, &wcount_ref);
      const int n_far = brute_force(cells, h_new, skip_id, 1, &wcount_far);
      const double rel = fabs(wcount_code - wcount_ref) / wcount_ref;
      /* Rounding bound of a float sum of n_ref positive terms. */
      const double rel_tol = n_ref * FLT_EPSILON;
      const int ok = (rel < rel_tol) && (n_code < 0 || n_code == n_ref);
      if (!ok) n_fail++;

      message(
          "%s loop, kernel_gamma*h/octant: old %.2f new %.2f: neighbours %d "
          "(code %d, far half %d), wcount code %.6e brute force %.6e (far half "
          "%.6e), relative difference %.3e (bar %.3e): %s",
          loop_name[loop], r_old, r_new[k], n_ref, n_code, n_far, wcount_code,
          wcount_ref, wcount_far, rel, rel_tol, ok ? "agree" : "MISMATCH");
    }
  }

  for (int c = 0; c < 2; c++) {
    free(cells[c].hydro.parts);
    free(cells[c].hydro.xparts);
    free(cells[c].stars.parts);
#ifdef WITH_BH_LOOP
    free(cells[c].black_holes.parts);
#endif
#ifdef WITH_SINK_LOOP
    free(cells[c].sinks.parts);
#endif
  }

  if (n_fail) error("%d subset density sums differ from brute force", n_fail);
  return 0;
}
