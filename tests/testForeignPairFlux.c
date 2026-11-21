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
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "space_getsid.h"
#include "swift.h"

/*
 * Checks the force pair loop of a scheme that updates the inactive
 * neighbours of a particle, when one cell of the pair is foreign.
 *
 * Two adjacent top-level cells of width 1 are split once into octants. Some
 * gas has an h that puts its pairs at the top level, the rest at the level of
 * the octants. Half of the gas is active. The pair is run three times: with
 * both cells local, with the second cell foreign and with the first cell
 * foreign. The flux that each local particle receives must equal the flux it
 * receives when both cells are local, and a brute-force sum over the pairs of
 * the two cells. Each case runs with the cells split and unsplit.
 *
 * The bar is relative to the largest sum of the absolute pair terms of a
 * particle of the cell. A different order of the same float sums, or of the
 * two particles of a pair, gives about 1e-6 of it. A missing or a doubled
 * pair term gives the size of the term, so a term near the kernel edge can
 * stay below the bar.
 */

#if hydro_force_updates_inactive_neighbours

void runner_dosub_pair2_force(struct runner *r, struct cell *ci,
                              struct cell *cj, const float h_hi,
                              const int gettimer);

#define TI_CURRENT 8
#define ACTIVE_BIN 1
#define INACTIVE_BIN 4
#define MAX_ACTIVE_BIN 3

/* Gas particles per octant. */
#define N_OCTANT 27
#define N_CELL (8 * N_OCTANT)

/* Number of flux components compared: momentum and energy. */
#define N_FLUX 4

#define REL_TOL 1e-4

/* A flux of one particle and the sum of its absolute pair terms. */
struct flux {
  double f[N_FLUX];
  double abs_sum[N_FLUX];
};

/* Exponent test, which -ffast-math cannot remove. */
static int is_finite(const double x) {
  uint64_t bits;
  memcpy(&bits, &x, sizeof(bits));
  return ((bits >> 52) & 0x7ff) != 0x7ff;
}

static void get_flux(const struct part *p, double f[N_FLUX]) {
  for (int k = 0; k < 3; k++) f[k] = p->flux.momentum[k];
  f[3] = p->flux.energy;
}

/**
 * @brief Set up a Gizmo particle ready for the force loop.
 */
static void init_particle(struct part *p, struct xpart *xp, const double m,
                          const int active, const struct engine *e) {
  p->conserved.mass = m;
  p->conserved.energy = 1.f;
  hydro_first_init_part(p, xp);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) p->geometry.matrix_E[i][j] = (i == j);
  p->geometry.volume = 1.f;
  p->geometry.wcorr = 1.f;
  p->rho = random_uniform(0.5, 1.5);
  p->P = random_uniform(0.5, 1.5);
  hydro_prepare_force(p, xp, e->cosmology, e->hydro_properties,
                      e->pressure_floor_props, 0., 0.);
  hydro_reset_acceleration(p);
  for (int k = 0; k < 3; k++) {
    p->v[k] = random_uniform(-1., 1.);
    p->flux.momentum[k] = 0.f;
  }
  p->flux.energy = 0.f;
  p->flux.dt = active ? 1e-3f : -1.f;
  p->time_bin = active ? ACTIVE_BIN : INACTIVE_BIN;
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_drift = TI_CURRENT;
  p->ti_kick = TI_CURRENT;
#endif
}

static void set_cell_properties(struct cell *c, const struct engine *e) {
  c->hydro.h_max = 0.f;
  c->hydro.h_max_active = 0.f;
  int any_active = 0;
  for (int k = 0; k < c->hydro.count; k++) {
    const struct part *p = &c->hydro.parts[k];
    c->hydro.h_max = max(c->hydro.h_max, p->h);
    if (part_is_active(p, e)) {
      any_active = 1;
      c->hydro.h_max_active = max(c->hydro.h_max_active, p->h);
    }
  }
  c->hydro.ti_end_min = any_active ? TI_CURRENT : 2 * TI_CURRENT;
  c->hydro.ti_old_part = TI_CURRENT;
  c->dmin = c->width[0];
  c->h_min_allowed = c->dmin * 0.5 * (1. / kernel_gamma);
  c->h_max_allowed = c->dmin * (1. / kernel_gamma);
}

/**
 * @brief Build a top-level cell of width 1 split into eight leaf octants.
 *
 * @param loc The corner of the cell.
 * @param p_top The probability that a particle has its pairs at the top level.
 * @param all_active Are all the particles active?
 */
static struct cell *make_cell(const double loc[3], const double p_top,
                              const int all_active, const struct engine *e) {
  const double width = 1.;
  const double h_top_min = width * 0.5 / kernel_gamma;
  const double h_top_max = width / kernel_gamma;
  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, 9 * sizeof(struct cell)) != 0)
    error("Failed to allocate the cells.");
  bzero(c, 9 * sizeof(struct cell));
  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     N_CELL * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&c->hydro.xparts, xpart_align,
                     N_CELL * sizeof(struct xpart)) != 0)
    error("Failed to allocate the particles.");
  bzero(c->hydro.parts, N_CELL * sizeof(struct part));
  bzero(c->hydro.xparts, N_CELL * sizeof(struct xpart));
  for (int i = 0; i < 3; i++) {
    c->loc[i] = loc[i];
    c->width[i] = width;
  }
  c->hydro.count = N_CELL;
  c->top = c;
  c->hydro.super = c;
  c->split = 1;
  for (int k = 0; k < 8; k++) {
    const double sub_loc[3] = {loc[0] + 0.5 * width * ((k >> 2) & 1),
                               loc[1] + 0.5 * width * ((k >> 1) & 1),
                               loc[2] + 0.5 * width * (k & 1)};
    for (int n = 0; n < N_OCTANT; n++) {
      const int ind = k * N_OCTANT + n;
      struct part *p = &c->hydro.parts[ind];
      for (int i = 0; i < 3; i++)
        p->x[i] = sub_loc[i] + 0.5 * width * random_uniform(0.02, 0.98);
      if (random_uniform(0., 1.) < p_top)
        p->h = random_uniform(1.02 * h_top_min, 0.98 * h_top_max);
      else
        p->h = random_uniform(0.52 * h_top_min, 0.98 * h_top_min);
      p->id = ind;
      const int active = all_active || random_uniform(0., 1.) < 0.5;
      init_particle(p, &c->hydro.xparts[ind], 1.0 / N_CELL, active, e);
    }
    struct cell *cp = &c[1 + k];
    for (int i = 0; i < 3; i++) {
      cp->loc[i] = sub_loc[i];
      cp->width[i] = 0.5 * width;
    }
    cp->hydro.parts = c->hydro.parts + k * N_OCTANT;
    cp->hydro.xparts = c->hydro.xparts + k * N_OCTANT;
    cp->hydro.count = N_OCTANT;
    cp->depth = 1;
    cp->parent = c;
    cp->top = c;
    cp->hydro.super = c;
    c->progeny[k] = cp;
  }
  set_cell_properties(c, e);
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    set_cell_properties(cp, e);
    for (int n = 0; n < cp->hydro.count; n++)
      cell_set_part_h_depth(&cp->hydro.parts[n], cp);
  }
  return c;
}

static void free_cell(struct cell *c) {
  for (int k = 0; k < 8; k++) cell_free_hydro_sorts(c->progeny[k]);
  cell_free_hydro_sorts(c);
  free(c->hydro.parts);
  free(c->hydro.xparts);
  free(c);
}

static void set_node(struct cell *c, const int node_id) {
  c->nodeID = node_id;
  for (int k = 0; k < 8; k++) c->progeny[k]->nodeID = node_id;
}

/**
 * @brief Add the flux of each pair term that a particle of a receives from
 * the particles of b, by brute force over all the pairs.
 */
static void brute_force(const struct cell *a, const struct cell *b,
                        struct flux *out, const struct engine *e) {
  const float a_fac = e->cosmology->a;
  const float H = e->cosmology->H;
  struct part pk, pq;
  for (int i = 0; i < a->hydro.count; i++) {
    const struct part *k = &a->hydro.parts[i];
    for (int j = 0; j < b->hydro.count; j++) {
      const struct part *q = &b->hydro.parts[j];
      const int k_active = part_is_active(k, e);
      const int q_active = part_is_active(q, e);
      if (!k_active && !q_active) continue;
      float dx[3];
      for (int d = 0; d < 3; d++) dx[d] = k->x[d] - q->x[d];
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      const float h = max(k->h, q->h);
      if (r2 >= h * h * kernel_gamma2) continue;

      memcpy(&pk, k, sizeof(struct part));
      memcpy(&pq, q, sizeof(struct part));
      for (int d = 0; d < 3; d++) pk.flux.momentum[d] = 0.f;
      pk.flux.energy = 0.f;
      if (k_active && q_active) {
        runner_iact_force(r2, dx, k->h, q->h, &pk, &pq, a_fac, H);
      } else if (k_active) {
        runner_iact_nonsym_force(r2, dx, k->h, q->h, &pk, &pq, a_fac, H);
      } else {
        const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
        runner_iact_nonsym_force(r2, mdx, q->h, k->h, &pq, &pk, a_fac, H);
      }
      double f[N_FLUX];
      get_flux(&pk, f);
      for (int d = 0; d < N_FLUX; d++) {
        out[i].f[d] += f[d];
        out[i].abs_sum[d] += fabs(f[d]);
      }
    }
  }
}

/**
 * @brief Run the force pair task from fresh copies of the particles and
 * store the fluxes of the particles of both cells.
 */
static void run_pair(struct runner *r, struct cell *ci, struct cell *cj,
                     const struct part *si, const struct part *sj,
                     const int split, double (*fi)[N_FLUX],
                     double (*fj)[N_FLUX]) {
  memcpy(ci->hydro.parts, si, N_CELL * sizeof(struct part));
  memcpy(cj->hydro.parts, sj, N_CELL * sizeof(struct part));
  ci->split = cj->split = split;
  runner_dosub_pair2_force(r, ci, cj, /*h_hi=*/FLT_MAX, /*gettimer=*/0);
  for (int k = 0; k < N_CELL; k++) {
    get_flux(&ci->hydro.parts[k], fi[k]);
    get_flux(&cj->hydro.parts[k], fj[k]);
  }
}

/**
 * @brief Count the particles whose flux differs from the reference by more
 * than the bar, and return the largest difference in units of the bar scale.
 */
static int compare(const double (*f)[N_FLUX], const double (*f_ref)[N_FLUX],
                   const struct flux *bf, double *max_rel) {
  double scale[N_FLUX] = {0.};
  for (int k = 0; k < N_CELL; k++)
    for (int d = 0; d < N_FLUX; d++) scale[d] = max(scale[d], bf[k].abs_sum[d]);

  int n_bad = 0;
  for (int k = 0; k < N_CELL; k++) {
    int bad = 0;
    for (int d = 0; d < N_FLUX; d++) {
      if (!is_finite(f[k][d])) {
        bad = 1;
        continue;
      }
      const double diff = fabs(f[k][d] - f_ref[k][d]);
      if (scale[d] == 0.) {
        if (diff != 0.) bad = 1;
        continue;
      }
      const double rel = diff / scale[d];
      *max_rel = max(*max_rel, rel);
      if (rel > REL_TOL) bad = 1;
    }
    n_bad += bad;
  }
  return n_bad;
}

/**
 * @brief Count the inactive particles of a that have an active neighbour
 * in b, which only a foreign-pair path of the loop can reach.
 */
static int count_reached(const struct cell *a, const struct cell *b,
                         const struct engine *e) {
  int n = 0;
  for (int i = 0; i < a->hydro.count; i++) {
    const struct part *k = &a->hydro.parts[i];
    if (part_is_active(k, e)) continue;
    for (int j = 0; j < b->hydro.count; j++) {
      const struct part *q = &b->hydro.parts[j];
      if (!part_is_active(q, e)) continue;
      double r2 = 0.;
      for (int d = 0; d < 3; d++)
        r2 += (k->x[d] - q->x[d]) * (k->x[d] - q->x[d]);
      const double h = max(k->h, q->h);
      if (r2 < h * h * kernel_gamma2) {
        n++;
        break;
      }
    }
  }
  return n;
}

/**
 * @brief Run one configuration of the pair and check it.
 *
 * @return The number of failed checks.
 */
static int do_case(const char *label, struct runner *r, struct engine *e,
                   const double offset[3], const double p_top,
                   const int all_active) {

  const double origin[3] = {0., 0., 0.};
  struct cell *ci = make_cell(origin, p_top, all_active, e);
  struct cell *cj = make_cell(offset, p_top, all_active, e);
  runner_do_hydro_sort(r, ci, 0x1FFF, 0, 0, 0, 0);
  runner_do_hydro_sort(r, cj, 0x1FFF, 0, 0, 0, 0);

  struct part *si = malloc(N_CELL * sizeof(struct part));
  struct part *sj = malloc(N_CELL * sizeof(struct part));
  memcpy(si, ci->hydro.parts, N_CELL * sizeof(struct part));
  memcpy(sj, cj->hydro.parts, N_CELL * sizeof(struct part));

  struct flux *bf_i = calloc(N_CELL, sizeof(struct flux));
  struct flux *bf_j = calloc(N_CELL, sizeof(struct flux));
  brute_force(ci, cj, bf_i, e);
  brute_force(cj, ci, bf_j, e);
  double (*bfi)[N_FLUX] = malloc(N_CELL * sizeof(*bfi));
  double (*bfj)[N_FLUX] = malloc(N_CELL * sizeof(*bfj));
  for (int k = 0; k < N_CELL; k++)
    for (int d = 0; d < N_FLUX; d++) {
      bfi[k][d] = bf_i[k].f[d];
      bfj[k][d] = bf_j[k].f[d];
    }

  struct cell *first = ci, *second = cj;
  double shift[3];
  space_getsid_and_swap_cells(e->s, &first, &second, shift);

  const int reach_i = count_reached(ci, cj, e);
  const int reach_j = count_reached(cj, ci, e);
  message("%s: first cell %s, inactive with an active neighbour: ci %d cj %d",
          label, first == ci ? "ci" : "cj", reach_i, reach_j);

  double (*ref_i)[N_FLUX] = malloc(N_CELL * sizeof(*ref_i));
  double (*ref_j)[N_FLUX] = malloc(N_CELL * sizeof(*ref_j));
  double (*fi)[N_FLUX] = malloc(N_CELL * sizeof(*fi));
  double (*fj)[N_FLUX] = malloc(N_CELL * sizeof(*fj));

  int n_fail = 0;
  for (int split = 0; split < 2; split++) {

    /* Both cells local. */
    set_node(ci, 0);
    set_node(cj, 0);
    run_pair(r, ci, cj, si, sj, split, ref_i, ref_j);
    double rel_local = 0.;
    const int bad_local = compare(ref_i, bfi, bf_i, &rel_local) +
                          compare(ref_j, bfj, bf_j, &rel_local);

    /* cj foreign: the fluxes of ci are this rank's. */
    double rel_i = 0., rel_i_bf = 0.;
    set_node(cj, 1);
    run_pair(r, ci, cj, si, sj, split, fi, fj);
    const int bad_i = compare(fi, ref_i, bf_i, &rel_i);
    const int bad_i_bf = compare(fi, bfi, bf_i, &rel_i_bf);
    double sum[N_FLUX] = {0.}, sum_abs = 0.;
    for (int k = 0; k < N_CELL; k++)
      for (int d = 0; d < N_FLUX; d++) {
        sum[d] += fi[k][d];
        sum_abs += bf_i[k].abs_sum[d];
      }

    /* ci foreign: the fluxes of cj are this rank's. */
    double rel_j = 0., rel_j_bf = 0.;
    set_node(cj, 0);
    set_node(ci, 1);
    run_pair(r, ci, cj, si, sj, split, fi, fj);
    const int bad_j = compare(fj, ref_j, bf_j, &rel_j);
    const int bad_j_bf = compare(fj, bfj, bf_j, &rel_j_bf);
    for (int k = 0; k < N_CELL; k++)
      for (int d = 0; d < N_FLUX; d++) {
        sum[d] += fj[k][d];
        sum_abs += bf_j[k].abs_sum[d];
      }
    set_node(ci, 0);

    /* The momentum of the pair is conserved when the two ranks add up. */
    const double sum_mom =
        sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
    const int bad_sum = !is_finite(sum_mom) || sum_mom > REL_TOL * sum_abs;

    message(
        "  %-7s all local vs brute force: %d bad (max %.1e) | ci local vs all "
        "local: %d bad (max %.1e), vs brute force %d | cj local vs all local: "
        "%d bad (max %.1e), vs brute force %d | |sum momentum| / sum |term| "
        "%.1e",
        split ? "split" : "unsplit", bad_local, rel_local, bad_i, rel_i,
        bad_i_bf, bad_j, rel_j, bad_j_bf, sum_mom / sum_abs);

    n_fail += (bad_local > 0) + (bad_i > 0) + (bad_i_bf > 0) + (bad_j > 0) +
              (bad_j_bf > 0) + bad_sum;
  }

  free(si);
  free(sj);
  free(bf_i);
  free(bf_j);
  free(bfi);
  free(bfj);
  free(ref_i);
  free(ref_j);
  free(fi);
  free(fj);
  free_cell(ci);
  free_cell(cj);
  return n_fail;
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  const unsigned int seed = (argc > 1) ? atoi(argv[1]) : 1234;
  srand(seed);
  message("Seed: %u", seed);

  struct space space;
  bzero(&space, sizeof(struct space));
  space.dim[0] = space.dim[1] = space.dim[2] = 4.;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct hydro_props hydro_props;
  hydro_props_init_no_hydro(&hydro_props);
  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(pressure_floor));
  struct phys_const phys_const;
  bzero(&phys_const, sizeof(phys_const));
  phys_const.const_vacuum_permeability = 1.;

  struct engine engine;
  bzero(&engine, sizeof(engine));
  engine.s = &space;
  engine.nodeID = 0;
  engine.ti_current = TI_CURRENT;
  engine.max_active_bin = MAX_ACTIVE_BIN;
  engine.time_base = 1e-3;
  engine.cosmology = &cosmo;
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;
  engine.physical_constants = &phys_const;

  struct runner runner;
  bzero(&runner, sizeof(runner));
  runner.e = &engine;

  const double face[3] = {1., 0., 0.};
  const double edge[3] = {1., 1., 0.};
  const double corner[3] = {1., 1., 1.};
  const double face_below[3] = {-1., 0., 0.};

  int n_fail = 0;
  n_fail += do_case("all active, face", &runner, &engine, face, 0.3, 1);
  n_fail += do_case("half active, face", &runner, &engine, face, 0.3, 0);
  n_fail += do_case("half active, face, swapped", &runner, &engine, face_below,
                    0.3, 0);
  n_fail +=
      do_case("half active, face, all top", &runner, &engine, face, 1.0, 0);
  n_fail += do_case("half active, edge", &runner, &engine, edge, 0.3, 0);
  n_fail += do_case("half active, corner", &runner, &engine, corner, 0.3, 0);

  if (n_fail) {
    message("%d checks failed.", n_fail);
    return 1;
  }
  message("All checks passed.");
  return 0;
}

#else

int main(int argc, char *argv[]) {
  /* The hydro scheme does not update inactive neighbours: skip. */
  return 77;
}

#endif
