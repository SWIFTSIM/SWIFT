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
#include "periodic.h"
#include "space_getsid.h"
#include "swift.h"

/*
 * Checks the pair loops that update the inactive neighbours of a particle,
 * when one cell of the pair is foreign: the Gizmo force loop and the GEAR RT
 * transport loop.
 *
 * Two adjacent top-level cells of width 1 are split into a tree of depth 1
 * or 2. Each particle has an h that puts its pairs at one level of the tree.
 * The pair is run three times: with both cells local, with the second cell
 * foreign and with the first cell foreign. The flux that each local particle
 * receives must equal the flux it receives when both cells are local, and a
 * brute-force sum over the pairs of the two cells. Each case runs with the
 * top-level cells split and unsplit.
 *
 * The bar is relative to the largest sum of the absolute pair terms of a
 * particle of the cell. A different order of the same float sums, or of the
 * two particles of a pair, gives about 1e-6 of it. A missing or a doubled
 * pair term gives the size of the term, so a term near the kernel edge can
 * stay below the bar.
 *
 * The RT loops also run the GEAR RT gradient loop, which updates the active
 * particles only. All the gas is hydro-inactive there, as on an RT sub-cycle,
 * so the search of the RT loops cannot rely on the hydro h_max_active. The
 * gradient loop also runs as a self task on one split cell.
 *
 * A bare configure skips this test (exit 77), and so does a hydro scheme that
 * does not update inactive neighbours. Run it by hand in a debug build, with
 * the seed as the first argument (./testForeignPairFlux 1234), after one of:
 *   ./configure --with-hydro=gizmo-mfm --with-riemann-solver=hllc
 *     --enable-debug --enable-debugging-checks --disable-optimization
 *     --disable-hand-vec
 *   ./configure --with-hydro=gizmo-mfv --with-riemann-solver=hllc
 *     --with-rt=GEAR_1 --with-rt-riemann-solver=GLF --with-stars=GEAR
 *     --with-feedback=none --with-grackle=<path> --enable-debug
 *     --enable-debugging-checks --disable-optimization --disable-hand-vec
 * Only the second line checks the RT transport and gradient loops.
 */

#if hydro_force_updates_inactive_neighbours

void runner_dosub_pair2_force(struct runner *r, struct cell *ci,
                              struct cell *cj, const float h_hi,
                              const int gettimer);
void runner_dosub_pair2_rt_transport(struct runner *r, struct cell *ci,
                                     struct cell *cj, const float h_hi,
                                     const int gettimer);
void runner_dosub_pair1_rt_gradient(struct runner *r, struct cell *ci,
                                    struct cell *cj, const float h_hi,
                                    const int gettimer);
void runner_dosub_self1_rt_gradient(struct runner *r, struct cell *c,
                                    const float h_hi, const int gettimer);

#define TI_CURRENT 8
#define ACTIVE_BIN 1
#define INACTIVE_BIN 4
#define MAX_ACTIVE_BIN 2

/* Gas particles per leaf cell. */
#define N_LEAF 27
#define MAX_DEPTH 2

/* Number of flux components compared: mass, momentum and energy, the
 * radiation energy and flux of the first photon group, or the gradient of its
 * energy density and of the x component of its flux. */
#define N_FLUX 5

#define REL_TOL 1e-4

enum loop_type { loop_force, loop_rt_transport, loop_rt_gradient };

static const char *loop_name[3] = {"force", "rt_transport", "rt_gradient"};

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

static void *checked_malloc(const size_t size) {
  void *ptr = malloc(size);
  if (ptr == NULL) error("Failed to allocate %zu bytes.", size);
  return ptr;
}

static void *checked_calloc(const size_t count, const size_t size) {
  void *ptr = calloc(count, size);
  if (ptr == NULL) error("Failed to allocate %zu bytes.", count * size);
  return ptr;
}

static int is_active(const struct part *p, const struct engine *e,
                     const enum loop_type loop) {
  return loop == loop_force ? part_is_active(p, e) : part_is_rt_active(p, e);
}

static void get_flux(const struct part *p, const enum loop_type loop,
                     double f[N_FLUX]) {
  for (int d = 0; d < N_FLUX; d++) f[d] = 0.;
  if (loop == loop_force) {
#ifdef GIZMO_MFV_SPH
    f[0] = p->flux.mass;
#endif
    for (int k = 0; k < 3; k++) f[1 + k] = p->flux.momentum[k];
    f[4] = p->flux.energy;
  } else if (loop == loop_rt_transport) {
#if rt_transport_updates_inactive_neighbours
    f[0] = p->rt_data.flux[0].energy;
    for (int k = 0; k < 3; k++) f[1 + k] = p->rt_data.flux[0].flux[k];
#endif
  } else {
#if rt_transport_updates_inactive_neighbours
    for (int k = 0; k < 3; k++) f[k] = p->rt_data.gradient[0].energy_density[k];
    f[3] = p->rt_data.gradient[0].flux[0][0];
    f[4] = p->rt_data.gradient[0].flux[0][1];
#endif
  }
}

static void reset_flux(struct part *p) {
#ifdef GIZMO_MFV_SPH
  p->flux.mass = 0.f;
#endif
  for (int k = 0; k < 3; k++) p->flux.momentum[k] = 0.f;
  p->flux.energy = 0.f;
#if rt_transport_updates_inactive_neighbours
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.flux[g].energy = 0.f;
    for (int k = 0; k < 3; k++) p->rt_data.flux[g].flux[k] = 0.f;
  }
  rt_gradients_init(p);
#endif
}

static void interact(const enum loop_type loop, const int symmetric,
                     const float r2, const float dx[3], const float hi,
                     const float hj, struct part *pi, struct part *pj,
                     const float a, const float H) {
  if (loop == loop_force) {
    if (symmetric)
      runner_iact_force(r2, dx, hi, hj, pi, pj, a, H);
    else
      runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj, a, H);
  } else if (loop == loop_rt_gradient) {
#if rt_transport_updates_inactive_neighbours
    runner_iact_nonsym_rt_gradient(r2, dx, hi, hj, pi, pj, a, H);
#endif
  } else {
#if rt_transport_updates_inactive_neighbours
    if (symmetric)
      runner_iact_rt_transport(r2, dx, hi, hj, pi, pj, a, H);
    else
      runner_iact_nonsym_rt_transport(r2, dx, hi, hj, pi, pj, a, H);
#endif
  }
}

/**
 * @brief Set up a Gizmo particle ready for the force and RT transport loops.
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
  for (int k = 0; k < 3; k++) p->v[k] = random_uniform(-1., 1.);
  p->flux.dt = active ? 1e-3f : -1.f;
  p->time_bin = active ? ACTIVE_BIN : INACTIVE_BIN;
#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  /* FVPM chemistry: the fluxes were last applied at the start of the step */
  p->chemistry_data.flux.ti_flux_last =
      get_integer_time_begin(TI_CURRENT, active ? ACTIVE_BIN : INACTIVE_BIN);
#endif
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_drift = TI_CURRENT;
  p->ti_kick = TI_CURRENT;
#endif

  p->rt_time_data.time_bin = active ? ACTIVE_BIN : INACTIVE_BIN;

#if rt_transport_updates_inactive_neighbours
  p->rt_data.flux_dt = active ? 1e-3f : -1.f;
  for (int g = 0; g < RT_NGROUPS; g++) {
    const float E = random_uniform(0.5, 1.5);
    p->rt_data.radiation[g].energy_density = E;
    for (int k = 0; k < 3; k++)
      p->rt_data.radiation[g].flux[k] = E * random_uniform(-0.3, 0.3);
  }
#ifdef SWIFT_RT_DEBUG_CHECKS
  p->rt_data.debug_kicked = 1;
  p->rt_data.debug_injection_done = 1;
  p->rt_data.debug_gradients_done = 1;
#endif
#endif
  reset_flux(p);
}

static void set_cell_properties(struct cell *c, const struct engine *e) {
  c->hydro.h_max = 0.f;
  c->hydro.h_max_active = 0.f;
  int any_active = 0, any_rt_active = 0;
  for (int k = 0; k < c->hydro.count; k++) {
    const struct part *p = &c->hydro.parts[k];
    c->hydro.h_max = max(c->hydro.h_max, p->h);
    if (part_is_active(p, e)) {
      any_active = 1;
      c->hydro.h_max_active = max(c->hydro.h_max_active, p->h);
    }
    if (part_is_rt_active(p, e)) any_rt_active = 1;
  }
  c->hydro.ti_end_min = any_active ? TI_CURRENT : 2 * TI_CURRENT;
  c->rt.ti_rt_end_min = any_rt_active ? TI_CURRENT : 2 * TI_CURRENT;
  c->hydro.ti_old_part = TI_CURRENT;
}

/* A tree of cells and its particles. */
struct tree {
  struct cell *cells;
  int n_cells;
  int depth;
};

static int pow8(const int n) {
  int m = 1;
  for (int l = 0; l < n; l++) m *= 8;
  return m;
}

static int tree_size(const int depth) {
  int n = 0;
  for (int l = 0; l <= depth; l++) n += pow8(l);
  return n;
}

/**
 * @brief Fill the cell c of the tree and its progeny with the particles
 * starting at parts.
 *
 * @param p_level The probability that a particle has its pairs at each level.
 * @param active_fraction The probability that a particle is active.
 */
static void build_cell(struct tree *t, struct cell *c, const double loc[3],
                       const double width, const int depth, struct part *parts,
                       struct xpart *xparts, const double *p_level,
                       const double active_fraction, const int n_total,
                       const struct engine *e) {
  for (int i = 0; i < 3; i++) {
    c->loc[i] = loc[i];
    c->width[i] = width;
  }
  c->depth = depth;
  c->dmin = width;
  c->h_min_allowed = width * 0.5 * (1. / kernel_gamma);
  c->h_max_allowed = width * (1. / kernel_gamma);
  c->top = t->cells;
  c->hydro.super = t->cells;
  c->hydro.parts = parts;
  c->hydro.xparts = xparts;

  if (depth == t->depth) {
    c->hydro.count = N_LEAF;
    for (int n = 0; n < N_LEAF; n++) {
      struct part *p = &parts[n];
      for (int i = 0; i < 3; i++)
        p->x[i] = loc[i] + width * random_uniform(0.02, 0.98);
      const double u = random_uniform(0., 1.);
      int level = 0;
      double p_cumul = p_level[0];
      while (level < t->depth && u >= p_cumul) p_cumul += p_level[++level];
      const double h_band = pow(0.5, level + 1) / kernel_gamma;
      p->h = random_uniform(1.04 * h_band, 1.96 * h_band);
      p->id = (parts - t->cells->hydro.parts) + n;
      const int active = random_uniform(0., 1.) < active_fraction;
      init_particle(p, &xparts[n], 1.0 / n_total, active, e);
    }
    set_cell_properties(c, e);
    for (int n = 0; n < N_LEAF; n++) cell_set_part_h_depth(&parts[n], c);
    return;
  }

  c->split = 1;
  const int count_progeny = n_total / pow8(depth + 1);
  for (int k = 0; k < 8; k++) {
    struct cell *cp = &t->cells[t->n_cells];
    t->n_cells += 1;
    const double sub_loc[3] = {loc[0] + 0.5 * width * ((k >> 2) & 1),
                               loc[1] + 0.5 * width * ((k >> 1) & 1),
                               loc[2] + 0.5 * width * (k & 1)};
    cp->parent = c;
    c->progeny[k] = cp;
    build_cell(t, cp, sub_loc, 0.5 * width, depth + 1,
               parts + k * count_progeny, xparts + k * count_progeny, p_level,
               active_fraction, n_total, e);
  }
  c->hydro.count = 8 * count_progeny;
  set_cell_properties(c, e);
}

/**
 * @brief Build a top-level cell of width 1 split into a tree of leaves.
 *
 * @param loc The corner of the cell.
 * @param depth The depth of the leaves.
 * @param p_level The probability that a particle has its pairs at each level.
 * @param active_fraction The probability that a particle is active.
 */
static struct cell *make_cell(const double loc[3], const int depth,
                              const double *p_level,
                              const double active_fraction,
                              const struct engine *e) {
  struct tree t;
  t.depth = depth;
  const int n_cells = tree_size(depth);
  const int n_total = N_LEAF * pow8(depth);
  if (posix_memalign((void **)&t.cells, cell_align,
                     n_cells * sizeof(struct cell)) != 0)
    error("Failed to allocate the cells.");
  bzero(t.cells, n_cells * sizeof(struct cell));
  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  if (posix_memalign((void **)&parts, part_align,
                     n_total * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&xparts, xpart_align,
                     n_total * sizeof(struct xpart)) != 0)
    error("Failed to allocate the particles.");
  bzero(parts, n_total * sizeof(struct part));
  bzero(xparts, n_total * sizeof(struct xpart));
  t.cells->hydro.parts = parts;
  t.n_cells = 1;
  build_cell(&t, t.cells, loc, 1., 0, parts, xparts, p_level, active_fraction,
             n_total, e);
  if (t.n_cells != n_cells) error("Wrong number of cells in the tree.");
  return t.cells;
}

static void free_sorts(struct cell *c) {
  if (c->split)
    for (int k = 0; k < 8; k++) free_sorts(c->progeny[k]);
  cell_free_hydro_sorts(c);
}

static void free_cell(struct cell *c) {
  c->split = 1;
  free_sorts(c);
  free(c->hydro.parts);
  free(c->hydro.xparts);
  free(c);
}

static void set_node(struct cell *c, const int node_id) {
  c->nodeID = node_id;
  if (c->progeny[0] != NULL)
    for (int k = 0; k < 8; k++) set_node(c->progeny[k], node_id);
}

/**
 * @brief Make all the gas of a tree hydro-inactive, as on an RT sub-cycle.
 */
static void set_hydro_inactive(struct cell *c, const struct engine *e) {
  if (c->split) {
    for (int k = 0; k < 8; k++) set_hydro_inactive(c->progeny[k], e);
  } else {
    for (int k = 0; k < c->hydro.count; k++) {
      c->hydro.parts[k].time_bin = INACTIVE_BIN;
      c->hydro.parts[k].flux.dt = -1.f;
    }
  }
  set_cell_properties(c, e);
}

/**
 * @brief Add the flux of each pair term that a particle of a receives from
 * the particles of b, by brute force over all the pairs.
 */
static void brute_force(const struct cell *a, const struct cell *b,
                        struct flux *out, const enum loop_type loop,
                        const struct engine *e) {
  const struct space *s = e->s;
  const float a_fac = e->cosmology->a;
  const float H = e->cosmology->H;
  struct part pk, pq;
  for (int i = 0; i < a->hydro.count; i++) {
    const struct part *k = &a->hydro.parts[i];
    for (int j = 0; j < b->hydro.count; j++) {
      const struct part *q = &b->hydro.parts[j];
      const int k_active = is_active(k, e, loop);
      const int q_active = is_active(q, e, loop);
      if (k == q) continue;
      if (!k_active && !q_active) continue;
      if (loop == loop_rt_gradient && !k_active) continue;
      float dx[3];
      for (int d = 0; d < 3; d++) {
        const double dx_d = k->x[d] - q->x[d];
        dx[d] = s->periodic ? nearest(dx_d, s->dim[d]) : dx_d;
      }
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      const float h = loop == loop_rt_gradient ? k->h : max(k->h, q->h);
      if (r2 >= h * h * kernel_gamma2) continue;

      memcpy(&pk, k, sizeof(struct part));
      memcpy(&pq, q, sizeof(struct part));
      reset_flux(&pk);
      if (k_active && q_active) {
        interact(loop, 1, r2, dx, k->h, q->h, &pk, &pq, a_fac, H);
      } else if (k_active) {
        interact(loop, 0, r2, dx, k->h, q->h, &pk, &pq, a_fac, H);
      } else {
        const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
        interact(loop, 0, r2, mdx, q->h, k->h, &pq, &pk, a_fac, H);
      }
      double f[N_FLUX];
      get_flux(&pk, loop, f);
      for (int d = 0; d < N_FLUX; d++) {
        out[i].f[d] += f[d];
        out[i].abs_sum[d] += fabs(f[d]);
      }
    }
  }
}

/**
 * @brief Run the pair task from fresh copies of the particles and store the
 * fluxes of the particles of both cells.
 */
static void run_pair(struct runner *r, struct cell *ci, struct cell *cj,
                     const struct part *si, const struct part *sj,
                     const int split, const enum loop_type loop,
                     double (*fi)[N_FLUX], double (*fj)[N_FLUX]) {
  const int count_i = ci->hydro.count, count_j = cj->hydro.count;
  memcpy(ci->hydro.parts, si, count_i * sizeof(struct part));
  memcpy(cj->hydro.parts, sj, count_j * sizeof(struct part));
  ci->split = cj->split = split;
  if (loop == loop_force)
    runner_dosub_pair2_force(r, ci, cj, /*h_hi=*/FLT_MAX, /*gettimer=*/0);
  else if (loop == loop_rt_transport)
    runner_dosub_pair2_rt_transport(r, ci, cj, /*h_hi=*/FLT_MAX,
                                    /*gettimer=*/0);
  else
    runner_dosub_pair1_rt_gradient(r, ci, cj, /*h_hi=*/FLT_MAX,
                                   /*gettimer=*/0);
  for (int k = 0; k < count_i; k++) get_flux(&ci->hydro.parts[k], loop, fi[k]);
  for (int k = 0; k < count_j; k++) get_flux(&cj->hydro.parts[k], loop, fj[k]);
}

/**
 * @brief Count the particles whose flux differs from the reference by more
 * than the bar, and return the largest difference in units of the bar scale.
 */
static int compare(const int count, const double (*f)[N_FLUX],
                   const double (*f_ref)[N_FLUX], const struct flux *bf,
                   double *max_rel) {
  double scale[N_FLUX] = {0.};
  for (int k = 0; k < count; k++)
    for (int d = 0; d < N_FLUX; d++) scale[d] = max(scale[d], bf[k].abs_sum[d]);

  int n_bad = 0;
  for (int k = 0; k < count; k++) {
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
                         const enum loop_type loop, const struct engine *e) {
  const struct space *s = e->s;
  int n = 0;
  for (int i = 0; i < a->hydro.count; i++) {
    const struct part *k = &a->hydro.parts[i];
    if (is_active(k, e, loop)) continue;
    for (int j = 0; j < b->hydro.count; j++) {
      const struct part *q = &b->hydro.parts[j];
      if (!is_active(q, e, loop)) continue;
      double r2 = 0.;
      for (int d = 0; d < 3; d++) {
        const double dx_d = k->x[d] - q->x[d];
        const double dx = s->periodic ? nearest(dx_d, s->dim[d]) : dx_d;
        r2 += dx * dx;
      }
      const double h = loop == loop_rt_gradient ? q->h : max(k->h, q->h);
      if (r2 < h * h * kernel_gamma2) {
        n++;
        break;
      }
    }
  }
  return n;
}

/* One configuration of the pair. */
struct pair_case {
  const char *label;
  double offset[3];
  int periodic;
  int depth;
  double p_level[MAX_DEPTH + 1];
  double active_i, active_j;
};

/**
 * @brief Run one configuration of the pair and check it.
 *
 * @return The number of failed checks.
 */
static int do_case(const struct pair_case *pc, const enum loop_type loop,
                   struct runner *r, struct engine *e) {

  e->s->periodic = pc->periodic;
  const double origin[3] = {0., 0., 0.};
  struct cell *ci = make_cell(origin, pc->depth, pc->p_level, pc->active_i, e);
  struct cell *cj =
      make_cell(pc->offset, pc->depth, pc->p_level, pc->active_j, e);
  if (loop != loop_force) {
    set_hydro_inactive(ci, e);
    set_hydro_inactive(cj, e);
  }
  runner_do_hydro_sort(r, ci, 0x1FFF, 0, 0, 0, 0);
  runner_do_hydro_sort(r, cj, 0x1FFF, 0, 0, 0, 0);

  const int n = ci->hydro.count;
  struct part *si = checked_malloc(n * sizeof(struct part));
  struct part *sj = checked_malloc(n * sizeof(struct part));
  memcpy(si, ci->hydro.parts, n * sizeof(struct part));
  memcpy(sj, cj->hydro.parts, n * sizeof(struct part));

  struct flux *bf_i = checked_calloc(n, sizeof(struct flux));
  struct flux *bf_j = checked_calloc(n, sizeof(struct flux));
  brute_force(ci, cj, bf_i, loop, e);
  brute_force(cj, ci, bf_j, loop, e);
  double (*bfi)[N_FLUX] = checked_malloc(n * sizeof(*bfi));
  double (*bfj)[N_FLUX] = checked_malloc(n * sizeof(*bfj));
  for (int k = 0; k < n; k++)
    for (int d = 0; d < N_FLUX; d++) {
      bfi[k][d] = bf_i[k].f[d];
      bfj[k][d] = bf_j[k].f[d];
    }

  struct cell *first = ci, *second = cj;
  double shift[3];
  space_getsid_and_swap_cells(e->s, &first, &second, shift);

  const int reach_i = count_reached(ci, cj, loop, e);
  const int reach_j = count_reached(cj, ci, loop, e);
  int n_deepest = 0;
  for (int k = 0; k < n; k++)
    n_deepest += (si[k].h < pow(0.5, pc->depth) / kernel_gamma);
  message(
      "%s %s: first cell %s, inactive with an active neighbour: ci %d cj %d, "
      "pairs at depth %d: %d",
      loop_name[loop], pc->label, first == ci ? "ci" : "cj", reach_i, reach_j,
      pc->depth, n_deepest);

  double (*ref_i)[N_FLUX] = checked_malloc(n * sizeof(*ref_i));
  double (*ref_j)[N_FLUX] = checked_malloc(n * sizeof(*ref_j));
  double (*fi)[N_FLUX] = checked_malloc(n * sizeof(*fi));
  double (*fj)[N_FLUX] = checked_malloc(n * sizeof(*fj));

  int n_fail = 0;
  for (int split = 0; split < 2; split++) {

    /* Both cells local. */
    set_node(ci, 0);
    set_node(cj, 0);
    run_pair(r, ci, cj, si, sj, split, loop, ref_i, ref_j);
    double rel_local = 0.;
    const int bad_local = compare(n, ref_i, bfi, bf_i, &rel_local) +
                          compare(n, ref_j, bfj, bf_j, &rel_local);

    /* cj foreign: the fluxes of ci are this rank's. */
    double rel_i = 0., rel_i_bf = 0.;
    set_node(cj, 1);
    run_pair(r, ci, cj, si, sj, split, loop, fi, fj);
    const int bad_i = compare(n, fi, ref_i, bf_i, &rel_i);
    const int bad_i_bf = compare(n, fi, bfi, bf_i, &rel_i_bf);
    double sum[N_FLUX] = {0.}, sum_abs[N_FLUX] = {0.};
    for (int k = 0; k < n; k++)
      for (int d = 0; d < N_FLUX; d++) {
        sum[d] += fi[k][d];
        sum_abs[d] += bf_i[k].abs_sum[d];
      }

    /* ci foreign: the fluxes of cj are this rank's. */
    double rel_j = 0., rel_j_bf = 0.;
    set_node(cj, 0);
    set_node(ci, 1);
    run_pair(r, ci, cj, si, sj, split, loop, fi, fj);
    const int bad_j = compare(n, fj, ref_j, bf_j, &rel_j);
    const int bad_j_bf = compare(n, fj, bfj, bf_j, &rel_j_bf);
    for (int k = 0; k < n; k++)
      for (int d = 0; d < N_FLUX; d++) {
        sum[d] += fj[k][d];
        sum_abs[d] += bf_j[k].abs_sum[d];
      }
    set_node(ci, 0);

    /* The first four components of a flux are conserved when the two ranks
     * add up (the last one is the internal energy flux of the force loop). A
     * gradient is not conserved. */
    double sum_rel = 0.;
    int bad_sum = 0;
    const int n_conserved = loop == loop_rt_gradient ? 0 : 4;
    for (int d = 0; d < n_conserved; d++) {
      if (!is_finite(sum[d])) bad_sum = 1;
      if (sum_abs[d] == 0.) continue;
      sum_rel = max(sum_rel, fabs(sum[d]) / sum_abs[d]);
    }
    bad_sum |= sum_rel > REL_TOL;

    message(
        "  %-7s all local vs brute force: %d bad (max %.1e) | ci local vs all "
        "local: %d bad (max %.1e), vs brute force %d | cj local vs all local: "
        "%d bad (max %.1e), vs brute force %d | |sum| / sum |term| %.1e",
        split ? "split" : "unsplit", bad_local, rel_local, bad_i, rel_i,
        bad_i_bf, bad_j, rel_j, bad_j_bf, sum_rel);

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

#if rt_transport_updates_inactive_neighbours
/**
 * @brief Run the RT gradient self task on one cell with all the gas
 * hydro-inactive and check the gradients against a brute-force sum.
 *
 * Half of the particles are in the top h band, which the self task handles at
 * the level of the cell, where the search of an RT-active particle must not
 * rely on h_max_active.
 *
 * @return The number of failed checks.
 */
static int do_self_case(const int depth, const double *p_level,
                        struct runner *r, struct engine *e) {

  e->s->periodic = 0;
  const double origin[3] = {0., 0., 0.};
  struct cell *c = make_cell(origin, depth, p_level, 0.5, e);
  set_hydro_inactive(c, e);
  set_node(c, 0);
  runner_do_hydro_sort(r, c, 0x1FFF, 0, 0, 0, 0);

  const int n = c->hydro.count;
  struct part *s0 = checked_malloc(n * sizeof(struct part));
  memcpy(s0, c->hydro.parts, n * sizeof(struct part));

  struct flux *bf = checked_calloc(n, sizeof(struct flux));
  brute_force(c, c, bf, loop_rt_gradient, e);
  double (*ref)[N_FLUX] = checked_malloc(n * sizeof(*ref));
  double (*f)[N_FLUX] = checked_malloc(n * sizeof(*f));
  for (int k = 0; k < n; k++)
    for (int d = 0; d < N_FLUX; d++) ref[k][d] = bf[k].f[d];

  int n_fail = 0;
  for (int split = 0; split < 2; split++) {
    memcpy(c->hydro.parts, s0, n * sizeof(struct part));
    c->split = split;
    runner_dosub_self1_rt_gradient(r, c, /*h_hi=*/FLT_MAX, /*gettimer=*/0);
    for (int k = 0; k < n; k++)
      get_flux(&c->hydro.parts[k], loop_rt_gradient, f[k]);
    double max_rel = 0.;
    const int n_bad = compare(n, f, ref, bf, &max_rel);
    message(
        "rt_gradient self, depth %d, %-7s vs brute force: %d bad (max %.1e)",
        depth, split ? "split" : "unsplit", n_bad, max_rel);
    n_fail += n_bad > 0;
  }

  free(s0);
  free(bf);
  free(ref);
  free(f);
  free_cell(c);
  return n_fail;
}
#endif

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

#if rt_transport_updates_inactive_neighbours
  rt_params.reduced_speed_of_light = 1.f;
  rt_params.reduced_speed_of_light_inverse = 1.f;
#endif

  struct engine engine;
  bzero(&engine, sizeof(engine));
  engine.s = &space;
  engine.nodeID = 0;
  engine.ti_current = TI_CURRENT;
  engine.ti_current_subcycle = TI_CURRENT;
  engine.max_active_bin = MAX_ACTIVE_BIN;
  engine.max_active_bin_subcycle = MAX_ACTIVE_BIN;
  engine.time_base = 1e-3;
  engine.cosmology = &cosmo;
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;
  engine.physical_constants = &phys_const;

  struct runner runner;
  bzero(&runner, sizeof(runner));
  runner.e = &engine;

  const struct pair_case cases[] = {
      {"all active, face", {1., 0., 0.}, 0, 1, {0.3, 0.7}, 1., 1.},
      {"half active, face", {1., 0., 0.}, 0, 1, {0.3, 0.7}, 0.5, 0.5},
      {"half active, face, swapped", {-1., 0., 0.}, 0, 1, {0.3, 0.7}, 0.5, 0.5},
      {"half active, face, all top", {1., 0., 0.}, 0, 1, {1.0, 0.}, 0.5, 0.5},
      {"half active, edge", {1., 1., 0.}, 0, 1, {0.3, 0.7}, 0.5, 0.5},
      {"half active, corner", {1., 1., 1.}, 0, 1, {0.3, 0.7}, 0.5, 0.5},
      {"ci all inactive, face", {1., 0., 0.}, 0, 1, {0.3, 0.7}, 0., 0.5},
      {"half active, periodic face", {3., 0., 0.}, 1, 1, {0.3, 0.7}, 0.5, 0.5},
      {"half active, periodic corner",
       {3., 3., 3.},
       1,
       1,
       {0.3, 0.7},
       0.5,
       0.5},
      {"half active, depth 2, face",
       {1., 0., 0.},
       0,
       2,
       {0.15, 0.3, 0.55},
       0.5,
       0.5},
  };
  const int n_cases = sizeof(cases) / sizeof(cases[0]);

  int n_fail = 0;
  for (int c = 0; c < n_cases; c++)
    n_fail += do_case(&cases[c], loop_force, &runner, &engine);
#if rt_transport_updates_inactive_neighbours
  for (int c = 0; c < n_cases; c++)
    n_fail += do_case(&cases[c], loop_rt_transport, &runner, &engine);
  for (int c = 0; c < n_cases; c++)
    n_fail += do_case(&cases[c], loop_rt_gradient, &runner, &engine);
  const double p_self[MAX_DEPTH + 1] = {0.5, 0.5};
  n_fail += do_self_case(1, p_self, &runner, &engine);
  const double p_self2[MAX_DEPTH + 1] = {0.3, 0.3, 0.4};
  n_fail += do_self_case(2, p_self2, &runner, &engine);
#endif

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
