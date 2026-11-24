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
#include "space_getsid.h"
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
 * Case 4: the stars rerun goes deeper than the task activation. The cells
 * are split twice, the stars h_max stops the activation at the top level and
 * the gas was sorted on an earlier step, then drifted. Two controls corrupt
 * the drift bounds of the cells and must give a wrong sum: one shows that the
 * rerun reads the sorted gas two levels down, the other that the drift margin
 * stops the descent.
 * Case 5: the stars density loop of a self task split in two, which
 * interacts the pairs of progeny directly, sorted and unsorted. A control
 * corrupts the sort drift of a progeny and must give a wrong sum.
 * Case 6: the gas density pair task and subset rerun with unsorted, sorted,
 * beyond-limit and stale sorts planted. With debugging checks, the sort state
 * counters must count the planted states, here and in cases 4 and 5.
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
void runner_dosub_self_stars_density(struct runner *r, struct cell *c,
                                     const int offset, const int ntasks,
                                     const float h_hi, const int gettimer);

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
      hydro_set_mass(p, 1.f);
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
#if defined(FEEDBACK_GEAR)
    n_code = sp->feedback_data.num_ngbs;
#endif
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

/**
 * @brief Compare the sort state counters of one particle type at one site
 * with the expected ones (debugging checks only).
 *
 * @param label The name of the case for the output.
 * @param type The particle type of the sort: gas (0) or stars (1).
 * @param site The #runner_sort_site.
 * @param want Whether the unsorted, stale and beyond counts must be positive.
 * @param sums_ok Whether the sums of the loop agree.
 * @return 1 if the counts and the sums are as expected, 0 otherwise.
 */
static int check_sort_counts(const char *label, const int type,
                             const enum runner_sort_site site,
                             const int want[3], const int sums_ok) {
#ifdef SWIFT_DEBUG_CHECKS
  struct runner_sort_counters counters;
  runner_sort_counters_get(&counters);
  const int *n = counters.count[type][site];
  int ok = sums_ok;
  for (int k = 0; k < 3; k++)
    if ((n[k] > 0) != want[k]) ok = 0;
  message("%s, %s sorts: unsorted %d stale %d beyond %d, sums %s: %s", label,
          type == 0 ? "gas" : "stars", n[0], n[1], n[2],
          sums_ok ? "agree" : "MISMATCH", ok ? "agree" : "MISMATCH");
  return ok;
#else
  return sums_ok;
#endif
}

/* Cells of case 4: index 0 is the top cell, 1 + k its octant k and
 * 9 + 8 k + l the octant l of octant k. */
#define N_DEEP 73

/* Positions of case 4 for a gas cell on the right of the star cell. */
static const double x_star_deep = 0.999;
static const double yz_deep = 0.125;

/* The first gas particle drifts towards the star after its sort. The second
 * drifts out of a leaf that does not touch the leaf of the star. */
#define N_EXTRA 2
static const double x_extra_old[N_EXTRA] = {1.17, 1.26};
static const double x_extra_new[N_EXTRA] = {1.12, 1.20};

/* kernel_gamma * h of the stars h_max, larger than half the top cell. */
#define H_STARS_MAX_DEEP 0.6f

/**
 * @brief Build a top-level cell of width 1 split twice, with a gas lattice in
 * each of its 64 leaves.
 *
 * @param cells The #N_DEEP cells of the hierarchy.
 * @param loc_x The x position of the cell's lower corner.
 * @param x_star The position of the star, or NULL for a cell without stars.
 * @param x_extra The positions of the extra gas, or NULL for none.
 * @param ind_extra (return) The indices of the extra gas in the cell.
 */
static void make_deep_cell(struct cell *cells, const double loc_x,
                           const double *x_star,
                           const double x_extra[N_EXTRA][3], int *ind_extra) {

  const int n_max = 64 * N_OCTANT + N_EXTRA;
  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  struct spart *sparts = NULL;
  if (posix_memalign((void **)&parts, part_align,
                     n_max * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&xparts, xpart_align,
                     n_max * sizeof(struct xpart)) != 0 ||
      posix_memalign((void **)&sparts, spart_align, sizeof(struct spart)) != 0)
    error("Couldn't allocate the particles");
  bzero(parts, n_max * sizeof(struct part));
  bzero(xparts, n_max * sizeof(struct xpart));
  bzero(sparts, sizeof(struct spart));
  bzero(cells, N_DEEP * sizeof(struct cell));

  if (x_star != NULL) {
    for (int d = 0; d < 3; d++) sparts->x[d] = x_star[d];
    sparts->id = 1;
    sparts->time_bin = 1;
#if defined(FEEDBACK_GEAR) || defined(FEEDBACK_AGORA)
    sparts->feedback_data.will_do_feedback = 1;
#endif
#ifdef SWIFT_DEBUG_CHECKS
    sparts->ti_drift = TI_CURRENT;
    sparts->ti_kick = TI_CURRENT;
#endif
  }

  for (int i = 0; i < N_DEEP; i++) {
    struct cell *c = &cells[i];
    const int depth = (i == 0) ? 0 : (i < 9) ? 1 : 2;
    const int k = (depth == 1) ? i - 1 : (depth == 2) ? (i - 9) / 8 : 0;
    const int l = (depth == 2) ? (i - 9) % 8 : 0;
    const double width = 1. / (1 << depth);
    c->loc[0] = loc_x;
    for (int d = 0; d < 3; d++) {
      if (depth >= 1) c->loc[d] += 0.5 * ((k >> (2 - d)) & 1);
      if (depth == 2) c->loc[d] += 0.25 * ((l >> (2 - d)) & 1);
      c->width[d] = width;
    }
    c->dmin = width;
    c->depth = depth;
    c->split = (depth < 2);
    c->top = &cells[0];
    c->hydro.super = &cells[0];
    c->nodeID = 0;
    c->hydro.h_max = H_LATTICE;
    c->hydro.h_max_old = H_LATTICE;
    c->hydro.ti_end_min = TI_CURRENT;
    c->hydro.ti_old_part = TI_CURRENT;
    c->stars.ti_end_min = (x_star != NULL) ? TI_CURRENT : 2 * TI_CURRENT;
    c->stars.ti_old_part = TI_CURRENT;
    c->stars.parts = sparts;
    if (x_star != NULL && cell_distance_outside_box(c, x_star) == 0.f) {
      c->stars.count = 1;
      c->stars.h_max = H_STARS_MAX_DEEP / kernel_gamma;
      c->stars.h_max_old = c->stars.h_max;
    }
    if (depth == 1) {
      c->parent = &cells[0];
      cells[0].progeny[k] = c;
    } else if (depth == 2) {
      c->parent = &cells[1 + k];
      cells[1 + k].progeny[l] = c;
    }
  }

  /* Fill the leaves in the order of the particle arrays. */
  int offset = 0;
  for (int k = 0; k < 8; k++) {
    struct cell *c1 = &cells[1 + k];
    c1->hydro.parts = &parts[offset];
    c1->hydro.xparts = &xparts[offset];
    for (int l = 0; l < 8; l++) {
      struct cell *leaf = &cells[9 + 8 * k + l];
      leaf->hydro.parts = &parts[offset];
      leaf->hydro.xparts = &xparts[offset];
      int n = 0;
      for (int e = 0; x_extra != NULL && e < N_EXTRA; e++) {
        if (cell_distance_outside_box(leaf, x_extra[e]) > 0.f) continue;
        for (int d = 0; d < 3; d++) parts[offset + n].x[d] = x_extra[e][d];
        parts[offset + n].id = ID_PLANT + e;
        ind_extra[e] = offset + n;
        n++;
      }
      for (int ix = 0; ix < N_SIDE; ix++) {
        for (int iy = 0; iy < N_SIDE; iy++) {
          for (int iz = 0; iz < N_SIDE; iz++) {
            struct part *p = &parts[offset + n];
            p->x[0] = leaf->loc[0] + (ix + 0.5) * 0.25 / N_SIDE;
            p->x[1] = leaf->loc[1] + (iy + 0.5) * 0.25 / N_SIDE;
            p->x[2] = leaf->loc[2] + (iz + 0.5) * 0.25 / N_SIDE;
            p->id = 2 + offset + n;
            n++;
          }
        }
      }
      for (int i = 0; i < n; i++) {
        struct part *p = &parts[offset + i];
        p->h = H_LATTICE;
        hydro_set_mass(p, 1.f);
        p->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
        p->ti_drift = TI_CURRENT;
        p->ti_kick = TI_CURRENT;
#endif
      }
      leaf->hydro.count = n;
      c1->hydro.count += n;
      offset += n;
    }
    cells[0].hydro.count += c1->hydro.count;
  }
  cells[0].hydro.parts = parts;
  cells[0].hydro.xparts = xparts;
}

/**
 * @brief Free the particle and sort arrays of a case 4 cell.
 *
 * @param cells The #N_DEEP cells of the hierarchy.
 */
static void free_deep_cell(struct cell *cells) {
  for (int i = 0; i < N_DEEP; i++) cell_free_hydro_sorts(&cells[i]);
  free(cells[0].hydro.parts);
  free(cells[0].hydro.xparts);
  free(cells[0].stars.parts);
}

/**
 * @brief Drift one gas particle to a new position and raise the drift bounds
 * of its leaf and of the cells above it, as cell_drift_part() does.
 *
 * @param cells The #N_DEEP cells of the hierarchy.
 * @param ind The index of the particle in the cell.
 * @param x_new The new position.
 */
static void drift_deep_part(struct cell *cells, const int ind,
                            const double x_new[3]) {
  struct part *p = &cells[0].hydro.parts[ind];
  struct xpart *xp = &cells[0].hydro.xparts[ind];
  float dx2 = 0.f;
  for (int d = 0; d < 3; d++) {
    xp->x_diff[d] -= (float)(x_new[d] - p->x[d]);
    xp->x_diff_sort[d] -= (float)(x_new[d] - p->x[d]);
    dx2 += xp->x_diff_sort[d] * xp->x_diff_sort[d];
    p->x[d] = x_new[d];
  }
  const float dx = sqrtf(dx2);

  struct cell *leaf = NULL;
  for (int i = 9; i < N_DEEP; i++)
    if (p >= cells[i].hydro.parts &&
        p < cells[i].hydro.parts + cells[i].hydro.count)
      leaf = &cells[i];
  const float dx_out = cell_distance_outside_box(leaf, x_new);
  for (struct cell *c = leaf; c != NULL; c = c->parent) {
    c->hydro.dx_max_part = max(c->hydro.dx_max_part, dx);
    c->hydro.dx_max_sort = max(c->hydro.dx_max_sort, dx);
    c->hydro.dx_outside_max = max(c->hydro.dx_outside_max, dx_out);
  }
}

/**
 * @brief Run the stars density loop of the star of case 4 or 5 and compare
 * it with a brute-force sum over the gas of the cells.
 *
 * @param r The #runner.
 * @param star The star #cell hierarchy.
 * @param gas The gas #cell hierarchy, or NULL for a self task only.
 * @param kernel_radius kernel_gamma * h of the star.
 * @param ntasks 0 for the subset rerun of the ghost, else the number of
 * split self tasks of the density loop.
 * @param label The name of the case for the output.
 * @return 1 if the sums agree, 0 otherwise.
 */
static int run_deep_stars(struct runner *r, struct cell *star, struct cell *gas,
                          const float kernel_radius, const int ntasks,
                          const char *label) {

  struct spart *sp = star[0].stars.parts;
  const float h = kernel_radius / kernel_gamma;
  const int ind_star = 0;
  bzero(&sp->density, sizeof(sp->density));
#if defined(FEEDBACK_GEAR)
  bzero(&sp->feedback_data, sizeof(sp->feedback_data));
  sp->feedback_data.will_do_feedback = 1;
#endif
  sp->h = h;
  if (ntasks == 0) {
    runner_dosub_self_subset_stars_density(r, &star[0], sp, &ind_star, 1,
                                           /*h_max_redo=*/h, /*gettimer=*/0);
    runner_dosub_pair_subset_stars_density(r, &star[0], sp, &ind_star, 1,
                                           /*h_max_redo=*/h, &gas[0],
                                           /*gettimer=*/0);
  } else {
    for (int i = 0; i < N_DEEP; i++)
      if (star[i].stars.count > 0) star[i].stars.h_max_active = h;
    for (int offset = 0; offset < ntasks; offset++)
      runner_dosub_self_stars_density(r, &star[0], offset, ntasks,
                                      /*h_hi=*/FLT_MAX, /*gettimer=*/0);
  }

  const float h_inv = 1.f / h;
  double wcount_ref = 0.;
  int n_ref = 0;
  for (int c = 0; c < 2; c++) {
    if (c == 1 && gas == NULL) continue;
    const struct cell *top = (c == 0) ? &star[0] : &gas[0];
    for (int i = 0; i < top->hydro.count; i++) {
      const struct part *p = &top->hydro.parts[i];
      const float dx[3] = {(float)(sp->x[0] - p->x[0]),
                           (float)(sp->x[1] - p->x[1]),
                           (float)(sp->x[2] - p->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (fabsf(sqrtf(r2) / kernel_radius - 1.f) < 1e-3f)
        error("Gas %lld is too close to the kernel edge", p->id);
      if (r2 < h * h * kernel_gamma2) {
        float w, w_dx;
        kernel_deval(sqrtf(r2) * h_inv, &w, &w_dx);
        wcount_ref += w;
        n_ref++;
      }
    }
  }
  const double wcount_code = sp->density.wcount;
  const double rel = fabs(wcount_code - wcount_ref) / wcount_ref;
  /* Rounding bound of a float sum of n_ref positive terms. */
  const double rel_tol = n_ref * FLT_EPSILON;
  const int ok = (rel < rel_tol);

  message(
      "%s, kernel_gamma*h/leaf %.2f: neighbours %d, wcount code %.6e brute "
      "force %.6e, relative difference %.3e (bar %.3e): %s",
      label, kernel_radius / 0.25f, n_ref, wcount_code, wcount_ref, rel,
      rel_tol, ok ? "agree" : "MISMATCH");
  return ok;
}

/**
 * @brief Case 4: the stars rerun descends below the level where the task
 * activation stopped.
 *
 * @param r The #runner.
 * @param s The #space.
 * @param mirror Whether the gas cell is on the left of the star cell.
 * @return The number of failed checks.
 */
static int check_deep_stars(struct runner *r, struct space *s,
                            const int mirror) {

  static struct cell star[N_DEEP], gas[N_DEEP];
  const double x_star[3] = {mirror ? 2. - x_star_deep : x_star_deep, yz_deep,
                            yz_deep};
  double x_old[N_EXTRA][3], x_new[N_EXTRA][3];
  for (int e = 0; e < N_EXTRA; e++) {
    x_old[e][0] = mirror ? 2. - x_extra_old[e] : x_extra_old[e];
    x_new[e][0] = mirror ? 2. - x_extra_new[e] : x_extra_new[e];
    x_old[e][1] = x_old[e][2] = x_new[e][1] = x_new[e][2] = yz_deep;
  }
  int ind_extra[N_EXTRA];
  make_deep_cell(star, mirror ? 1. : 0., x_star, NULL, NULL);
  make_deep_cell(gas, mirror ? 0. : 1., NULL, x_old, ind_extra);

  /* The gas was sorted on an earlier step along the axis of the pair only. */
  struct cell *ci = &star[0], *cj = &gas[0];
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);
  runner_do_hydro_sort(r, &gas[0], 1 << sid, /*cleanup=*/0, /*lock=*/0,
                       /*rt_requests_sort=*/0, /*clock=*/0);

  /* Then drifted. */
  for (int e = 0; e < N_EXTRA; e++)
    drift_deep_part(gas, ind_extra[e], x_new[e]);

  /* Activate the pair task as the unskip does. */
  static struct task tasks[8];
  static int tid_active[8];
  bzero(tasks, sizeof(tasks));
  for (int t = 0; t < 8; t++) tasks[t].skip = 1;
  struct scheduler sched;
  bzero(&sched, sizeof(sched));
  sched.space = s;
  sched.tasks = tasks;
  sched.tid_active = tid_active;
  for (int c = 0; c < 2; c++) {
    struct cell *top = (c == 0) ? &star[0] : &gas[0];
    top->hydro.drift = &tasks[4 * c];
    top->stars.drift = &tasks[4 * c + 1];
    top->hydro.sorts = &tasks[4 * c + 2];
    top->stars.sorts = &tasks[4 * c + 3];
  }
  cell_activate_subcell_stars_tasks(&star[0], &gas[0], &sched,
                                    /*with_star_formation=*/0,
                                    /*with_star_formation_sink=*/0,
                                    /*with_timestep_sync=*/0);

  /* The activation stops at the top level, drifts all the gas and asks no
   * new sort. */
  int n_deep_flags = 0;
  for (int i = 1; i < N_DEEP; i++)
    if (cell_get_flag(&gas[i], cell_flag_do_hydro_drift) ||
        gas[i].hydro.requires_sorts)
      n_deep_flags++;
  const int activation_ok = cell_get_flag(&gas[0], cell_flag_do_hydro_drift) &&
                            gas[0].hydro.requires_sorts == (1 << sid) &&
                            gas[0].hydro.do_sort == 0 && n_deep_flags == 0;
  message(
      "deep stars rerun%s, activation: drift of the top gas cell %d, sorts "
      "required %d asked %d, sub-cells flagged %d: %s",
      mirror ? " (mirrored)" : "",
      cell_get_flag(&gas[0], cell_flag_do_hydro_drift) != 0,
      gas[0].hydro.requires_sorts, gas[0].hydro.do_sort, n_deep_flags,
      activation_ok ? "agree" : "MISMATCH");
  int n_fail = activation_ok ? 0 : 1;

  /* The rerun reaches the leaves and reads the sort of the gas there. */
  const float r_deep = 0.15f, r_margin = 0.24f;
  runner_sort_counters_reset();
  int sums_ok = run_deep_stars(r, star, gas, r_deep, /*ntasks=*/0,
                               "deep stars rerun, leaves");
  if (!sums_ok) n_fail++;

  /* The gas of the star cell is not sorted, and the sort of the leaf of the
   * drifted gas is beyond its limit. */
  const int want[3] = {1, 0, 1};
  if (!check_sort_counts("deep stars rerun, leaves", /*type=*/0,
                         runner_sort_site_fallback, want, sums_ok))
    n_fail++;

  /* The gas planted on the sort axis drifted towards the star: a sorted
   * sweep without the drift of its leaf misses it. */
  /* A naive build never reads the sort: the control cannot mismatch there. */
#ifndef SWIFT_USE_NAIVE_INTERACTIONS
  struct cell *leaf_sorted = NULL;
  for (int i = 9; i < N_DEEP; i++)
    if (cell_distance_outside_box(&gas[i], x_new[0]) == 0.f)
      leaf_sorted = &gas[i];
  const float dx_sort_save = leaf_sorted->hydro.dx_max_sort;
  leaf_sorted->hydro.dx_max_sort = 0.f;
  if (run_deep_stars(r, star, gas, r_deep, /*ntasks=*/0,
                     "control without the sort drift of the leaf (must "
                     "mismatch)")) {
    message("The control did not mismatch: the leaf sort was not read");
    n_fail++;
  }
  leaf_sorted->hydro.dx_max_sort = dx_sort_save;
#endif

  /* The drift margin keeps the rerun of a larger star above the leaves. */
  if (!run_deep_stars(r, star, gas, r_margin, /*ntasks=*/0,
                      "deep stars rerun, margin"))
    n_fail++;

  float dx_out_save[N_DEEP];
  for (int i = 0; i < N_DEEP; i++) {
    dx_out_save[i] = gas[i].hydro.dx_outside_max;
    gas[i].hydro.dx_outside_max = 0.f;
  }
  if (run_deep_stars(r, star, gas, r_margin, /*ntasks=*/0,
                     "control without the drift margin (must mismatch)")) {
    message("The control did not mismatch: the margin is not tested");
    n_fail++;
  }
  for (int i = 0; i < N_DEEP; i++) gas[i].hydro.dx_outside_max = dx_out_save[i];

  free_deep_cell(star);
  free_deep_cell(gas);
  return n_fail;
}

/**
 * @brief Case 5: the stars density loop of a self task split in two.
 *
 * The star sits next to the face between two octants. Its kernel reaches
 * gas of the other octant that drifted towards it after the sort.
 *
 * @param r The #runner.
 * @param sorted Whether the gas and the stars are sorted.
 * @return The number of failed checks.
 */
static int check_split_self_stars(struct runner *r, const int sorted) {

  static struct cell star[N_DEEP];
  const double x_star[3] = {0.49, 0.26, 0.26};
  const double x_old[N_EXTRA][3] = {{0.72, 0.26, 0.26}, {0.1, 0.9, 0.9}};
  const double x_new[3] = {0.66, 0.26, 0.26};
  int ind_extra[N_EXTRA];
  make_deep_cell(star, 0., x_star, x_old, ind_extra);
  if (sorted) {
    runner_do_hydro_sort(r, &star[0], 0x1FFF, /*cleanup=*/0, /*lock=*/0,
                         /*rt_requests_sort=*/0, /*clock=*/0);
    runner_do_stars_sort(r, &star[0], 0x1FFF, /*cleanup=*/0, /*clock=*/0);
  }
  drift_deep_part(star, ind_extra[0], x_new);

  const float kernel_radius = 0.2f;
  int n_fail = 0;
  const char *label =
      sorted ? "split stars self, sorted" : "split stars self, unsorted";
  runner_sort_counters_reset();
  const int sums_ok =
      run_deep_stars(r, star, NULL, kernel_radius, /*ntasks=*/2, label);
  if (!sums_ok) n_fail++;

  /* Unsorted, the pairs fall back on naive interactions. Sorted, the octant
   * of the drifted gas is beyond its limit. */
  const int want_unsorted[3] = {1, 0, 0}, want_none[3] = {0, 0, 0};
  const int want_beyond[3] = {0, 0, 1};
  if (!check_sort_counts(label, /*type=*/0, runner_sort_site_fallback,
                         sorted ? want_beyond : want_unsorted, sums_ok) ||
      !check_sort_counts(label, /*type=*/1, runner_sort_site_fallback,
                         sorted ? want_none : want_unsorted, sums_ok))
    n_fail++;

  /* The octant of the drifted gas is a progeny of the top cell. */
  /* A naive stars build never reads the sort: the control cannot mismatch. */
#ifndef SWIFT_USE_NAIVE_INTERACTIONS_STARS
  if (sorted) {
    struct cell *octant = &star[1 + 4];
    const float dx_sort_save = octant->hydro.dx_max_sort;
    octant->hydro.dx_max_sort = 0.f;
    if (run_deep_stars(r, star, NULL, kernel_radius, /*ntasks=*/2,
                       "control without the sort drift of the octant (must "
                       "mismatch)")) {
      message("The control did not mismatch: the sorted pair was not run");
      n_fail++;
    }
    octant->hydro.dx_max_sort = dx_sort_save;
  }
#endif

  for (int i = 0; i < N_DEEP; i++) cell_free_stars_sorts(&star[i]);
  free_deep_cell(star);
  return n_fail;
}

#ifdef SWIFT_DEBUG_CHECKS

void runner_dosub_pair1_density(struct runner *r, struct cell *ci,
                                struct cell *cj, const float h_hi,
                                const int gettimer);

/* kernel_gamma * h of the gas of case 6, below the leaf width. */
#define R_GAS_SORT 0.2f

/**
 * @brief Run the gas density loop of case 6 and compare every density sum with
 * a brute-force sum over the gas of the other cell.
 *
 * @param r The #runner.
 * @param a The first gas #cell hierarchy.
 * @param b The second gas #cell hierarchy.
 * @param ind_subset The index in a of the particle of a subset rerun, or -1
 * for the pair task.
 * @return 1 if the sums agree, 0 otherwise.
 */
static int run_gas_pair(struct runner *r, struct cell *a, struct cell *b,
                        const int ind_subset) {

  const float h = R_GAS_SORT / kernel_gamma;
  for (int c = 0; c < 2; c++) {
    struct cell *top = (c == 0) ? a : b;
    for (int i = 0; i < top->hydro.count; i++)
      bzero(&top->hydro.parts[i].density, sizeof(top->hydro.parts[i].density));
  }

  runner_sort_counters_reset();
  if (ind_subset < 0)
    runner_dosub_pair1_density(r, &a[0], &b[0], /*h_hi=*/FLT_MAX,
                               /*gettimer=*/0);
  else
    runner_dosub_pair_subset_density(r, &a[0], a[0].hydro.parts, &ind_subset, 1,
                                     /*h_max_redo=*/h, &b[0],
                                     /*gettimer=*/0);

  const float h_inv = 1.f / h;
  int ok = 1;
  for (int c = 0; c < 2; c++) {
    const struct cell *top = (c == 0) ? a : b;
    const struct cell *other = (c == 0) ? b : a;
    for (int i = 0; i < top->hydro.count; i++) {
      if (ind_subset >= 0 && (c == 1 || i != ind_subset)) continue;
      const struct part *pi = &top->hydro.parts[i];
      double wcount_ref = 0.;
      int n_ref = 0;
      for (int j = 0; j < other->hydro.count; j++) {
        const struct part *pj = &other->hydro.parts[j];
        const float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                             (float)(pi->x[1] - pj->x[1]),
                             (float)(pi->x[2] - pj->x[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        if (r2 < h * h * kernel_gamma2) {
          float w, w_dx;
          kernel_deval(sqrtf(r2) * h_inv, &w, &w_dx);
          wcount_ref += w;
          n_ref++;
        }
      }
      const double diff = fabs(pi->density.wcount - wcount_ref);
      if (diff > (n_ref + 1) * FLT_EPSILON * max(wcount_ref, 1.)) ok = 0;
    }
  }
  return ok;
}

/**
 * @brief Case 6: the sort state counters of the gas density loops.
 *
 * @param r The #runner.
 * @param s The #space.
 * @return The number of failed checks.
 */
static int check_sort_counters(struct runner *r, struct space *s) {

  static struct cell a[N_DEEP], b[N_DEEP];
  make_deep_cell(a, 0., NULL, NULL, NULL);
  make_deep_cell(b, 1., NULL, NULL, NULL);
  const float h = R_GAS_SORT / kernel_gamma;
  for (int c = 0; c < 2; c++) {
    struct cell *top = (c == 0) ? a : b;
    for (int i = 0; i < top[0].hydro.count; i++) top[0].hydro.parts[i].h = h;
    for (int i = 0; i < N_DEEP; i++)
      top[i].hydro.h_max = top[i].hydro.h_max_old = top[i].hydro.h_max_active =
          h;
  }
  struct cell *ci = &a[0], *cj = &b[0];
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  const int recurse_size_save = space_recurse_size_pair_hydro;
  space_recurse_size_pair_hydro = 1;
  int n_fail = 0;

  /* The leaves are not sorted: the pair leaves sort them on demand. */
  const int want_unsorted[3] = {1, 0, 0}, want_none[3] = {0, 0, 0};
  const int want_stale[3] = {0, 1, 0}, want_beyond[3] = {0, 0, 1};
  int sums_ok = run_gas_pair(r, a, b, /*ind_subset=*/-1);
  if (!check_sort_counts("sort counters, unsorted pair leaves", /*type=*/0,
                         runner_sort_site_pair_leaf, want_unsorted, sums_ok))
    n_fail++;

  /* Sorted and valid at every level. */
  runner_do_hydro_sort(r, &a[0], 1 << sid, /*cleanup=*/0, /*lock=*/0,
                       /*rt_requests_sort=*/0, /*clock=*/0);
  runner_do_hydro_sort(r, &b[0], 1 << sid, /*cleanup=*/0, /*lock=*/0,
                       /*rt_requests_sort=*/0, /*clock=*/0);
  sums_ok = run_gas_pair(r, a, b, /*ind_subset=*/-1);
  if (!check_sort_counts("sort counters, sorted pair leaves", /*type=*/0,
                         runner_sort_site_pair_leaf, want_none, sums_ok))
    n_fail++;

  /* A leaf of b on the shared face moved beyond its sort limit after its
   * sort, but within the limit of its parents. */
  struct cell *leaf = NULL;
  for (int i = 9; i < N_DEEP && leaf == NULL; i++)
    if (b[i].loc[0] == 1.) leaf = &b[i];
  const float dx_beyond = 1.5f * space_maxreldx * leaf->dmin;
  for (struct cell *c = leaf; c != NULL; c = c->parent)
    c->hydro.dx_max_sort = max(c->hydro.dx_max_sort, dx_beyond);
  sums_ok = run_gas_pair(r, a, b, /*ind_subset=*/-1);
  if (!check_sort_counts("sort counters, pair leaf beyond its limit",
                         /*type=*/0, runner_sort_site_pair_leaf, want_beyond,
                         sums_ok))
    n_fail++;
  for (int i = 0; i < N_DEEP; i++) b[i].hydro.dx_max_sort = 0.f;

  /* The subset rerun of a particle of a on the shared face. */
  int ind = -1;
  for (int i = 0; i < a[0].hydro.count && ind < 0; i++)
    if (a[0].hydro.parts[i].x[0] > 0.9) ind = i;
  sums_ok = run_gas_pair(r, a, b, ind);
  if (!check_sort_counts("sort counters, sorted subset", /*type=*/0,
                         runner_sort_site_fallback, want_none, sums_ok))
    n_fail++;

  /* The sorts of b are stale: the subset falls back on naive interactions. */
  for (int i = 0; i < N_DEEP; i++)
    b[i].hydro.dx_max_sort_old = 2.f * space_maxreldx * b[i].dmin;
  sums_ok = run_gas_pair(r, a, b, ind);
  if (!check_sort_counts("sort counters, stale subset", /*type=*/0,
                         runner_sort_site_fallback, want_stale, sums_ok))
    n_fail++;
  for (int i = 0; i < N_DEEP; i++) b[i].hydro.dx_max_sort_old = 0.f;

  space_recurse_size_pair_hydro = recurse_size_save;
  free_deep_cell(a);
  free_deep_cell(b);
  return n_fail;
}
#endif /* SWIFT_DEBUG_CHECKS */

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

  /* Case 4: the stars rerun goes deeper than the task activation. */
  n_fail += check_deep_stars(&runner, &space, /*mirror=*/0);
  n_fail += check_deep_stars(&runner, &space, /*mirror=*/1);

  /* Case 5: a split stars self task, with a single star to recurse on. */
  space_recurse_size_self_stars = 1;
  n_fail += check_split_self_stars(&runner, /*sorted=*/0);
  n_fail += check_split_self_stars(&runner, /*sorted=*/1);

  /* Case 6: the sort state counters of the gas pair and subset loops. */
#ifdef SWIFT_DEBUG_CHECKS
  n_fail += check_sort_counters(&runner, &space);
#else
  message("Case 6 skipped: the sort state counters need debugging checks");
#endif

  if (n_fail) error("%d subset density checks failed", n_fail);
  return 0;
}
