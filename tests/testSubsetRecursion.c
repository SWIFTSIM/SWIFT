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
 * while their smoothing length converges. The density sums of the subset
 * loops must equal a brute-force sum over all the gas of the two cells.
 *
 * Two adjacent top-level cells of width 1 are split once into octants of
 * width 0.5. The gas sits on a lattice. The particle being redone sits just
 * inside the face shared by the two cells. Extra gas is planted in the octant
 * of the second cell that does not touch the shared face.
 *
 * Case 1: the cells keep the h of the last step while the particle being
 * redone already has a larger h.
 * Case 2: the planted gas has drifted out of its octant, inside the kernel.
 * Only the black hole and sink loops are checked: the stars and hydro loops
 * place each particle at a depth set from its h alone, so their reruns keep
 * the drift of the searching cell only. Pass --with-depth-drift to run these
 * loops too.
 * Case 3: the activation of a black hole density task must drift every gas
 * sub-cell that a rerun with a grown h can read.
 */

void runner_dosub_pair_subset_stars_density(struct runner *r, struct cell *ci,
                                            struct spart *sparts,
                                            const int *ind, const int scount,
                                            const float h_max_redo,
                                            struct cell *cj,
                                            const int gettimer);
void runner_dosub_self_subset_stars_density(struct runner *r, struct cell *ci,
                                            struct spart *sparts,
                                            const int *ind, const int scount,
                                            const float h_max_redo,
                                            const int gettimer);
void runner_dosub_pair_subset_density(struct runner *r, struct cell *ci,
                                      struct part *parts, const int *ind,
                                      const int count, const float h_max_redo,
                                      struct cell *cj, const int gettimer);
void runner_dosub_self_subset_density(struct runner *r, struct cell *ci,
                                      struct part *parts, const int *ind,
                                      const int count, const float h_max_redo,
                                      const int gettimer);
void runner_dosub_subset_bh_density(struct runner *r, struct cell *ci,
                                    struct bpart *bparts, int *ind,
                                    const int bcount, const float h_max_redo,
                                    struct cell *cj, int gettimer);
void runner_dosub_subset_sinks_density(struct runner *r, struct cell *ci,
                                       struct sink *sinks, int *ind,
                                       const int scount, const float h_max_redo,
                                       struct cell *cj, int gettimer);

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

/* Ids of the planted gas start here. */
#define ID_PLANT 200000

/* Position of the particle being redone, just inside the shared face. */
static const double x_redo[3] = {0.999, 0.25, 0.25};

/* Planted gas, in octant 4 of the second cell (x from 1.5 to 2). */
#define N_PLANT 4
static const double x_plant_inside[N_PLANT][3] = {{1.505, 0.25, 0.25},
                                                  {1.505, 0.30, 0.25},
                                                  {1.505, 0.25, 0.30},
                                                  {1.510, 0.20, 0.20}};

/* The same gas drifted by up to 0.04 out of its octant. */
static const double x_plant_drifted[N_PLANT][3] = {{1.46, 0.25, 0.25},
                                                   {1.46, 0.30, 0.25},
                                                   {1.46, 0.25, 0.30},
                                                   {1.47, 0.20, 0.20}};

/* Drift of the gas of the second cell in case 2. */
#define DX_DRIFT 0.05f

static const char *loop_name[4] = {"stars", "hydro", "black-hole", "sink"};

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
 * @param x_plant The positions of the planted gas.
 */
static void make_split_cell(struct cell *c, struct cell *progeny,
                            const double loc_x, const int with_redo,
                            const double x_plant[N_PLANT][3]) {

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
        p->id = ID_PLANT + n;
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
    cc->hydro.super = c;
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
    cc->top = c;
    if (k < 8) {
      cc->parent = c;
      cc->depth = 1;
      c->progeny[k] = cc;
    }
  }
}

/**
 * @brief Free the particle arrays of the two top-level cells.
 *
 * @param cells The two top-level #cell.
 */
static void free_cells(struct cell *cells) {
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
}

/**
 * @brief Brute-force density sums of a particle at x_redo over the gas of
 * both cells.
 *
 * @param cells The two top-level #cell.
 * @param h The smoothing length of the particle.
 * @param skip_id Gas id to leave out (the particle itself), or -1.
 * @param far_only Whether to count only the gas of the octants of the second
 * cell that do not touch the shared face.
 * @param wcount (return) The kernel-weighted neighbour count.
 * @return The number of neighbours.
 */
static int brute_force(struct cell *cells, const float h,
                       const long long skip_id, const int far_only,
                       double *wcount) {
  const float h_inv = 1.f / h;
  const struct part *far_start = cells[1].progeny[4]->hydro.parts;
  int n_ngb = 0;
  *wcount = 0.;
  for (int c = 0; c < 2; c++) {
    if (far_only && c == 0) continue;
    for (int i = 0; i < cells[c].hydro.count; i++) {
      const struct part *p = &cells[c].hydro.parts[i];
      if (p->id == skip_id) continue;
      if (far_only && p < far_start) continue;
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

/**
 * @brief Set the values of the last step on every cell, as the activation
 * of the tasks does.
 *
 * @param cells The two top-level #cell.
 * @param h_old The h of the last step of the particle being redone.
 * @param dx_j The drift of the gas of the second cell.
 */
static void set_old_values(struct cell *cells, const float h_old,
                           const float dx_j) {
  for (int c = 0; c < 2; c++) {
    for (int j = 0; j < 9; j++) {
      struct cell *cc = (j < 8) ? cells[c].progeny[j] : &cells[c];
      cc->stars.h_max_old = (c == 0) ? h_old : 0.f;
#ifdef WITH_BH_LOOP
      cc->black_holes.h_max_old = (c == 0) ? h_old : 0.f;
#endif
#ifdef WITH_SINK_LOOP
      cc->sinks.h_max_old = (c == 0) ? h_old : 0.f;
#endif
      cc->hydro.h_max_old = (c == 0) ? max(h_old, H_LATTICE) : H_LATTICE;
      cc->hydro.dx_max_part_old = (c == 1) ? dx_j : 0.f;
      cc->hydro.dx_max_part = cc->hydro.dx_max_part_old;
      cc->hydro.dx_outside_max = cc->hydro.dx_max_part_old;
    }
  }
}

/**
 * @brief Run the subset density loops of one particle type and compare the
 * result with the brute-force sum.
 *
 * @param r The #runner.
 * @param cells The two top-level #cell.
 * @param loop The loop: 0 stars, 1 hydro, 2 black holes, 3 sinks.
 * @param h_new The h of the particle being redone.
 * @param label The name of the case for the output.
 * @return 1 if the sums agree, 0 otherwise.
 */
static int run_subset_loop(struct runner *r, struct cell *cells, const int loop,
                           const float h_new, const char *label) {

  struct part *p_redo = cells[0].progeny[4]->hydro.parts;
  const int ind_gas = (int)(p_redo - cells[0].hydro.parts);
  const int ind_star = 0;
  int ind_one = 0;
  (void)ind_one;

  float wcount_code = 0.f;
  int n_code = -1;
  if (loop == 0) {
    struct spart *sp = cells[0].stars.parts;
    bzero(&sp->density, sizeof(sp->density));
#if defined(FEEDBACK_GEAR)
    bzero(&sp->feedback_data, sizeof(sp->feedback_data));
    sp->feedback_data.will_do_feedback = 1;
#endif
    sp->h = h_new;
    runner_dosub_self_subset_stars_density(r, &cells[0], cells[0].stars.parts,
                                           &ind_star, 1, /*h_max_redo=*/h_new,
                                           /*gettimer=*/0);
    runner_dosub_pair_subset_stars_density(r, &cells[0], cells[0].stars.parts,
                                           &ind_star, 1, /*h_max_redo=*/h_new,
                                           &cells[1], /*gettimer=*/0);
    wcount_code = sp->density.wcount;
  } else if (loop == 1) {
    bzero(&p_redo->density, sizeof(p_redo->density));
    p_redo->h = h_new;
    runner_dosub_self_subset_density(r, &cells[0], cells[0].hydro.parts,
                                     &ind_gas, 1, /*h_max_redo=*/h_new,
                                     /*gettimer=*/0);
    runner_dosub_pair_subset_density(r, &cells[0], cells[0].hydro.parts,
                                     &ind_gas, 1, /*h_max_redo=*/h_new,
                                     &cells[1], /*gettimer=*/0);
    wcount_code = p_redo->density.wcount;
    p_redo->h = H_LATTICE;
  } else if (loop == 2) {
#ifdef WITH_BH_LOOP
    struct bpart *bp = cells[0].black_holes.parts;
    bzero(&bp->density, sizeof(bp->density));
    bp->h = h_new;
#if defined(BLACK_HOLES_GEAR)
    bp->num_ngbs = 0;
#endif
    runner_dosub_subset_bh_density(r, &cells[0], bp, &ind_one, 1,
                                   /*h_max_redo=*/h_new, NULL,
                                   /*gettimer=*/0);
    runner_dosub_subset_bh_density(r, &cells[0], bp, &ind_one, 1,
                                   /*h_max_redo=*/h_new, &cells[1],
                                   /*gettimer=*/0);
    wcount_code = bp->density.wcount;
#if defined(BLACK_HOLES_GEAR)
    n_code = bp->num_ngbs;
#endif
#endif
  } else {
#ifdef WITH_SINK_LOOP
    struct sink *si = cells[0].sinks.parts;
    bzero(&si->density, sizeof(si->density));
    si->h = h_new;
    si->num_ngbs = 0;
    runner_dosub_subset_sinks_density(r, &cells[0], si, &ind_one, 1,
                                      /*h_max_redo=*/h_new, NULL,
                                      /*gettimer=*/0);
    runner_dosub_subset_sinks_density(r, &cells[0], si, &ind_one, 1,
                                      /*h_max_redo=*/h_new, &cells[1],
                                      /*gettimer=*/0);
    wcount_code = si->density.wcount;
    n_code = si->num_ngbs;
#endif
  }

  const long long skip_id = (loop == 1) ? p_redo->id : -1;
  double wcount_ref, wcount_far;
  const int n_ref = brute_force(cells, h_new, skip_id, 0, &wcount_ref);
  const int n_far = brute_force(cells, h_new, skip_id, 1, &wcount_far);
  const double rel = fabs(wcount_code - wcount_ref) / wcount_ref;
  /* Rounding bound of a float sum of n_ref positive terms. */
  const double rel_tol = n_ref * FLT_EPSILON;
  const int ok = (rel < rel_tol) && (n_code < 0 || n_code == n_ref);

  message(
      "%s, %s loop, kernel_gamma*h/octant %.2f: neighbours %d (code %d, far "
      "octants %d), wcount code %.6e brute force %.6e (far octants %.6e), "
      "relative difference %.3e (bar %.3e): %s",
      label, loop_name[loop], h_new * kernel_gamma / 0.5f, n_ref, n_code, n_far,
      wcount_code, wcount_ref, wcount_far, rel, rel_tol,
      ok ? "agree" : "MISMATCH");
  return ok;
}

#ifdef WITH_BH_LOOP

/**
 * @brief Is the gas of a leaf cell drifted this step, given the drift flags?
 *
 * A cell flagged for drift drifts its whole sub-tree (cell_drift_part()).
 *
 * @param c The leaf #cell.
 */
static int leaf_is_drifted(const struct cell *c) {
  for (const struct cell *f = c; f != NULL; f = f->parent)
    if (cell_get_flag(f, cell_flag_do_hydro_drift)) return 1;
  return 0;
}

/**
 * @brief Activate a black hole density pair task as the unskip does, then
 * count the gas inside the grown kernel that will not be drifted.
 *
 * @param cells The two top-level #cell.
 * @param s The #scheduler.
 * @param h_old The h of the last step of the black hole.
 * @param h_new The h of the black hole in the rerun.
 * @return 1 if every gas particle inside the kernel is drifted.
 */
static int check_bh_drift_activation(struct cell *cells, struct scheduler *s,
                                     const float h_old, const float h_new) {

  static struct task tasks[16];
  static int tid_active[16];
  bzero(tasks, sizeof(tasks));
  for (int t = 0; t < 16; t++) tasks[t].skip = 1;
  s->tasks = tasks;
  s->tid_active = tid_active;
  s->active_count = 0;

  struct task *density = &tasks[0];
  density->type = task_type_pair;
  density->subtype = task_subtype_bh_density;
  density->ci = &cells[0];
  density->cj = &cells[1];
  struct link link = {density, NULL};

  for (int c = 0; c < 2; c++) {
    cells[c].hydro.drift = &tasks[1 + 4 * c];
    cells[c].black_holes.drift = &tasks[2 + 4 * c];
    cells[c].black_holes.swallow_ghost_1 = &tasks[3 + 4 * c];
    cells[c].black_holes.black_holes_in = &tasks[4 + 4 * c];
    for (int j = 0; j < 9; j++) {
      struct cell *cc = (j < 8) ? cells[c].progeny[j] : &cells[c];
      cc->flags = 0;
      cc->hydro.h_max = H_LATTICE;
      cc->hydro.dx_max_part = 0.f;
      cc->black_holes.h_max = (c == 0) ? h_old : 0.f;
      cc->black_holes.dx_max_part = 0.f;
      /* Only the cells holding the black hole are active. */
      const int has_bh = (c == 0) && (j == 4 || j == 8);
      cc->black_holes.ti_end_min = has_bh ? TI_CURRENT : 2 * TI_CURRENT;
    }
  }
  cells[0].black_holes.density = &link;

  cell_unskip_black_holes_tasks(&cells[0], s);
  cells[0].black_holes.density = NULL;

  /* Gas of the second cell inside the grown kernel, by octant. */
  int n_missed = 0, n_in = 0;
  for (int k = 0; k < 8; k++) {
    const struct cell *cp = cells[1].progeny[k];
    int n_kernel = 0;
    for (int i = 0; i < cp->hydro.count; i++) {
      const struct part *p = &cp->hydro.parts[i];
      double r2 = 0.;
      for (int d = 0; d < 3; d++)
        r2 += (x_redo[d] - p->x[d]) * (x_redo[d] - p->x[d]);
      if (r2 < h_new * h_new * kernel_gamma2) n_kernel++;
    }
    n_in += n_kernel;
    if (!leaf_is_drifted(cp)) n_missed += n_kernel;
  }

  const int ok = (n_missed == 0);
  message(
      "black-hole drift activation, kernel_gamma*h/octant old %.2f new %.2f: "
      "%d gas of the second cell in the kernel, %d of them not drifted: %s",
      h_old * kernel_gamma / 0.5f, h_new * kernel_gamma / 0.5f, n_in, n_missed,
      ok ? "agree" : "MISMATCH");
  return ok;
}

#endif /* WITH_BH_LOOP */

int main(int argc, char *argv[]) {

  clocks_set_cpufreq(0);
  const int with_depth_drift =
      (argc > 1 && strcmp(argv[1], "--with-depth-drift") == 0);

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
  space.e = &engine;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  struct cell cells[2];
  struct cell progeny[2][8];

  /* kernel_gamma * h in units of the octant width. Only values below 1 allow
   * recursion. */
  const float r_old = 0.9f;
  const float h_old = r_old * 0.5f / kernel_gamma;

  int n_fail = 0;

  /* Case 1: h grows past the octant, or shrinks. */
  make_split_cell(&cells[0], progeny[0], 0., /*with_redo=*/1, x_plant_inside);
  make_split_cell(&cells[1], progeny[1], 1., /*with_redo=*/0, x_plant_inside);
  space.cells_top = cells;
  const float r_grow[3] = {1.1f, 1.8f, 0.8f};
  for (int k = 0; k < 3; k++) {
    const float h_new = r_grow[k] * 0.5f / kernel_gamma;
    for (int loop = 0; loop < 4; loop++) {
#ifndef WITH_BH_LOOP
      if (loop == 2) continue;
#endif
#ifndef WITH_SINK_LOOP
      if (loop == 3) continue;
#endif
      set_old_values(cells, h_old, 0.f);
      if (!run_subset_loop(&runner, cells, loop, h_new, "h growth")) n_fail++;
    }
  }
  free_cells(cells);

  /* Case 2: gas drifted out of the far octant into the kernel. */
  make_split_cell(&cells[0], progeny[0], 0., /*with_redo=*/1, x_plant_drifted);
  make_split_cell(&cells[1], progeny[1], 1., /*with_redo=*/0, x_plant_drifted);
  space.cells_top = cells;
  {
    const float h_new = 0.95f * 0.5f / kernel_gamma;
    for (int loop = 0; loop < 4; loop++) {
      if (loop < 2 && !with_depth_drift) continue;
#ifndef WITH_BH_LOOP
      if (loop == 2) continue;
#endif
#ifndef WITH_SINK_LOOP
      if (loop == 3) continue;
#endif
      set_old_values(cells, h_new, DX_DRIFT);
      if (!run_subset_loop(&runner, cells, loop, h_new, "neighbour drift"))
        n_fail++;
    }
  }
  free_cells(cells);

#ifdef WITH_BH_LOOP
  /* Case 3: the black hole task activation drifts the gas a rerun reads. */
  make_split_cell(&cells[0], progeny[0], 0., /*with_redo=*/1, x_plant_inside);
  make_split_cell(&cells[1], progeny[1], 1., /*with_redo=*/0, x_plant_inside);
  space.cells_top = cells;
  {
    struct scheduler sched;
    bzero(&sched, sizeof(sched));
    sched.space = &space;
    const float h_new = 1.1f * 0.5f / kernel_gamma;
    if (!check_bh_drift_activation(cells, &sched, h_old, h_new)) n_fail++;
  }
  free_cells(cells);
#endif

  if (n_fail) error("%d subset density checks failed", n_fail);
  return 0;
}
