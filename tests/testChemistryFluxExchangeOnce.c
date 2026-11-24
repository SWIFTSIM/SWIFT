/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#include <config.h>

/* Local headers. */
#include "space_getsid.h"
#include "swift.h"

/* System headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file testChemistryFluxExchangeOnce.c
 * @brief Check that the force and the chemistry prep loops apply the FVPM
 * chemistry flux of every interacting pair exactly once to each particle.
 *
 * The per-particle exchange counters (SWIFT_CHEMISTRY_DEBUG_CHECKS) are
 * compared to a brute-force count: one per pair with r < gamma max(hi, hj)
 * and at least one active particle. Cells are split up to twice, with h in
 * the depth band of each level, so that the recursion of the sub-cell tasks
 * and the active pairs in different depth bands are exercised. Further
 * cases pair a split with an unsplit cell, pair with a fully inactive cell,
 * wrap across the periodic box, put particles up to the cell faces, and
 * move the particles after the sort.
 *
 * A pair with a foreign cell misses exchanges in the pair enumeration itself
 * when the force loop does not update inactive neighbours (SPH). There, each
 * particle's count must equal the reference minus the exact number of misses
 * predicted for it (predict_foreign_gap()). The prep and apply counts of
 * every particle must be equal.
 */

#if (defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) ||             \
     defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)) && \
    defined(SWIFT_CHEMISTRY_DEBUG_CHECKS)

#define TI_CURRENT 8
#define ACTIVE_BIN 1
#define INACTIVE_BIN 4
#define MAX_ACTIVE_BIN 2
#define PARTS_PER_OCTANT 27
#define MAX_DEPTH 2
#define PARTS_PER_LEAF_DEPTH2 13
#define MAX_FRAMES (2 * (MAX_DEPTH + 1))

void runner_dosub_self2_force(struct runner *r, struct cell *c,
                              const float h_hi, const int gettimer);
void runner_dosub_pair2_force(struct runner *r, struct cell *ci,
                              struct cell *cj, const float h_hi,
                              const int gettimer);
void runner_dosub_self2_chemistry_flux_prep(struct runner *r, struct cell *c,
                                            const float h_hi,
                                            const int gettimer);
void runner_dosub_pair2_chemistry_flux_prep(struct runner *r, struct cell *ci,
                                            struct cell *cj, const float h_hi,
                                            const int gettimer);

/**
 * @brief Set the hydro and chemistry state of a particle for a force loop.
 *
 * @param p The #part.
 * @param xp Its #xpart.
 * @param m The mass of the particle.
 * @param active Is the particle active in this step?
 * @param e The #engine.
 */
static void init_particle(struct part *p, struct xpart *xp, const double m,
                          const int active, const struct engine *e) {

#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
  p->conserved.mass = m;
  p->conserved.energy = 1.f;
  hydro_first_init_part(p, xp);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) p->geometry.matrix_E[i][j] = (i == j);
  p->geometry.volume = 1.f;
  p->geometry.wcorr = 1.f;
  p->rho = 1.f;
  p->P = 1.f;
#else
  p->mass = m;
  p->u = 1.f;
  p->rho = 1.f;
#endif
  hydro_prepare_force(p, xp, e->cosmology, e->hydro_properties,
                      e->pressure_floor_props, 0., 0.);
  hydro_reset_acceleration(p);

  p->time_bin = active ? ACTIVE_BIN : INACTIVE_BIN;
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_drift = TI_CURRENT;
  p->ti_kick = TI_CURRENT;
#endif

  struct chemistry_part_data *ch = &p->chemistry_data;
  ch->kappa = 1.0;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  ch->tau = 1.0;
#endif
  ch->limiter.maxr = 1.f;
  for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
    ch->limiter.rhoZ[k][0] = -1e10;
    ch->limiter.rhoZ[k][1] = 1e10;
    ch->limiter.Z[k][0] = -1e10;
    ch->limiter.Z[k][1] = 1e10;
    ch->fct_theta[k] = 1.0;
    ch->metal_mass[k] = m * 1e-3 * random_uniform(0., 1.);
  }

  /* Each particle applied its fluxes last at the start of its step. */
  ch->flux.ti_flux_last =
      get_integer_time_begin(TI_CURRENT, active ? ACTIVE_BIN : INACTIVE_BIN);
}

/**
 * @brief Fill the cell properties derived from its particles.
 */
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
  c->hydro.dx_max_part = 0.f;
  c->hydro.dx_max_sort = 0.f;
  c->hydro.dx_max_sort_old = 0.f;
  c->hydro.sorted = 0;
  c->hydro.sort = NULL;
  c->dmin = c->width[0];
  c->h_min_allowed = c->dmin * 0.5 * (1. / kernel_gamma);
  c->h_max_allowed = c->dmin * (1. / kernel_gamma);
}

/**
 * @brief Set the geometry and the links of a progeny cell.
 *
 * @param cp The progeny cell.
 * @param parent Its parent.
 * @param k Its octant in the parent.
 * @param first The index of its first particle in the parent.
 * @param count Its number of particles.
 */
static void link_progeny(struct cell *cp, struct cell *parent, const int k,
                         const int first, const int count) {
  const double w = 0.5 * parent->width[0];
  cp->loc[0] = parent->loc[0] + w * ((k >> 2) & 1);
  cp->loc[1] = parent->loc[1] + w * ((k >> 1) & 1);
  cp->loc[2] = parent->loc[2] + w * (k & 1);
  for (int i = 0; i < 3; i++) cp->width[i] = w;
  cp->hydro.parts = parent->hydro.parts + first;
  cp->hydro.xparts = parent->hydro.xparts + first;
  cp->hydro.count = count;
  cp->depth = parent->depth + 1;
  cp->parent = parent;
  cp->top = parent->top;
  cp->hydro.super = parent->top;
  cp->nodeID = parent->nodeID;
  parent->progeny[k] = cp;
}

/**
 * @brief Fill the properties of a cell and of its progeny, and the depth_h
 * of the particles of its leaves.
 */
static void set_tree_properties(struct cell *c, const struct engine *e) {
  set_cell_properties(c, e);
  if (!c->split) {
    for (int n = 0; n < c->hydro.count; n++)
      cell_set_part_h_depth(&c->hydro.parts[n], c);
    return;
  }
  for (int k = 0; k < 8; k++) set_tree_properties(c->progeny[k], e);
}

/**
 * @brief Build a unit cell split @p depth times.
 *
 * A cell of depth 0 or 1 holds 8 * #PARTS_PER_OCTANT particles, a cell of
 * depth 2 holds 64 * #PARTS_PER_LEAF_DEPTH2. h is drawn from the band of one
 * of the levels of the cell, so the particles sit at several depths.
 *
 * @param loc The lower corner of the cell.
 * @param depth The number of splits, at most #MAX_DEPTH.
 * @param margin The smallest distance of a particle to the faces of its
 * octant (of its leaf at depth 2), in units of the octant width.
 * @param node_id The rank that owns the cell.
 * @param frac_active The probability that a particle is active.
 * @param top_band_active Are all particles in the band of the top cell
 * active?
 * @param e The #engine.
 * @return The cell; the cells below it follow it in the same allocation.
 */
static struct cell *make_cell(const double loc[3], const int depth,
                              const double margin, const int node_id,
                              const double frac_active,
                              const int top_band_active,
                              const struct engine *e) {

  const int n_units = (depth == 2) ? 64 : 8;
  const int parts_per_unit =
      (depth == 2) ? PARTS_PER_LEAF_DEPTH2 : PARTS_PER_OCTANT;
  const int count = n_units * parts_per_unit;
  const int n_cells = (depth == 2) ? 73 : 9;
  const double width = 1.;
  const double unit_width = (depth == 2) ? 0.25 * width : 0.5 * width;
  const double h_top_min = width * 0.5 / kernel_gamma;
  const double h_top_max = width / kernel_gamma;
  const double h_band[MAX_DEPTH + 1][2] = {
      {1.02 * h_top_min, 0.98 * h_top_max},
      {0.52 * h_top_min, 0.98 * h_top_min},
      {0.26 * h_top_min, 0.49 * h_top_min}};

  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, n_cells * sizeof(struct cell)) !=
      0)
    error("Could not allocate the cells.");
  bzero(c, n_cells * sizeof(struct cell));
  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&c->hydro.xparts, xpart_align,
                     count * sizeof(struct xpart)) != 0)
    error("Could not allocate the particles.");
  bzero(c->hydro.parts, count * sizeof(struct part));
  bzero(c->hydro.xparts, count * sizeof(struct xpart));

  for (int i = 0; i < 3; i++) {
    c->loc[i] = loc[i];
    c->width[i] = width;
  }
  c->hydro.count = count;
  c->nodeID = node_id;
  c->top = c;
  c->hydro.super = c;
  c->split = (depth > 0);

  for (int u = 0; u < n_units; u++) {
    const int k = (depth == 2) ? u / 8 : u;
    const int l = u % 8;
    double unit_loc[3] = {loc[0] + 0.5 * width * ((k >> 2) & 1),
                          loc[1] + 0.5 * width * ((k >> 1) & 1),
                          loc[2] + 0.5 * width * (k & 1)};
    if (depth == 2) {
      unit_loc[0] += unit_width * ((l >> 2) & 1);
      unit_loc[1] += unit_width * ((l >> 1) & 1);
      unit_loc[2] += unit_width * (l & 1);
    }
    for (int n = 0; n < parts_per_unit; n++) {
      const int ind = u * parts_per_unit + n;
      struct part *p = &c->hydro.parts[ind];
      for (int i = 0; i < 3; i++)
        p->x[i] =
            unit_loc[i] + unit_width * random_uniform(margin, 1. - margin);
      int level = 0;
      if (depth == 1 && random_uniform(0., 1.) >= 0.3) level = 1;
      if (depth == 2) {
        const double draw = random_uniform(0., 1.);
        if (draw >= 0.1) level = (draw >= 0.45) ? 2 : 1;
      }
      p->h = random_uniform(h_band[level][0], h_band[level][1]);
      p->id = ind + 1 + 10000 * node_id;
      const int active = (top_band_active && level == 0) ||
                         random_uniform(0., 1.) < frac_active;
      init_particle(p, &c->hydro.xparts[ind], 1.0 / count, active, e);
    }
  }

  if (depth > 0) {
    const int per_octant = count / 8;
    for (int k = 0; k < 8; k++) {
      struct cell *cp = &c[1 + k];
      link_progeny(cp, c, k, k * per_octant, per_octant);
      cp->split = (depth == 2);
      if (depth == 2)
        for (int l = 0; l < 8; l++)
          link_progeny(&c[9 + 8 * k + l], cp, l, l * PARTS_PER_LEAF_DEPTH2,
                       PARTS_PER_LEAF_DEPTH2);
    }
  }

  /* Foreign cells hold no xparts. */
  if (node_id != e->nodeID) {
    free(c->hydro.xparts);
    for (int k = 0; k < n_cells; k++) c[k].hydro.xparts = NULL;
  }

  set_tree_properties(c, e);
  return c;
}

/**
 * @brief Free the sorts of a cell and of its progeny.
 */
static void free_tree_sorts(struct cell *c) {
  if (c->split)
    for (int k = 0; k < 8; k++) free_tree_sorts(c->progeny[k]);
  cell_free_hydro_sorts(c);
}

static void free_cell(struct cell *c) {
  free_tree_sorts(c);
  free(c->hydro.parts);
  free(c->hydro.xparts);
  free(c);
}

static void reset_counts(struct cell *c) {
  for (int k = 0; k < c->hydro.count; k++) {
    c->hydro.parts[k].chemistry_data.flux_exchange_count[0] = 0;
    c->hydro.parts[k].chemistry_data.flux_exchange_count[1] = 0;
  }
}

/* Number of particle pairs whose in-kernel test depends on the frame or on
   the precision of the distance. */
static int n_ambiguous_pairs = 0;

/* The naive pair loop visits every pair from both sides, and a force loop
   that updates inactive neighbours lets a foreign active particle reach
   them: neither has a foreign gap. */
#if defined(SWIFT_USE_NAIVE_INTERACTIONS) || \
    hydro_force_updates_inactive_neighbours
static const int foreign_gap_exists = 0;
#else
static const int foreign_gap_exists = 1;
#endif

/**
 * @brief The progeny of a split cell that holds a particle.
 */
static const struct cell *progeny_of(const struct cell *c,
                                     const struct part *p) {
  for (int k = 0; k < 8; k++) {
    const struct cell *cp = c->progeny[k];
    if (p >= cp->hydro.parts && p < cp->hydro.parts + cp->hydro.count)
      return cp;
  }
  error("Particle not found in the progeny.");
  return NULL;
}

/* The frames in which the loops can measure the distance of a pair: the
   origins of the cells that hold either particle, at every level, and the
   periodic image of the second particle as seen from the first. */
struct pair_frames {
  const struct cell *a, *b;
  const double *loc[MAX_FRAMES];
  int of_second[MAX_FRAMES];
  int count;
  double image[3];
};

/**
 * @brief Prepare the frames of the pairs of a particle of a with a particle
 * of b. The origins are filled by within_kernel() when needed.
 */
static void pair_frames(const struct cell *a, const struct cell *b,
                        const struct space *s, struct pair_frames *f) {
  f->a = a;
  f->b = b;
  f->count = 0;
  for (int k = 0; k < 3; k++) {
    const double d = b->loc[k] - a->loc[k];
    f->image[k] = 0.;
    if (s->periodic && d > 0.5 * s->dim[k]) f->image[k] = -s->dim[k];
    if (s->periodic && d < -0.5 * s->dim[k]) f->image[k] = s->dim[k];
  }
}

/**
 * @brief Is pj within the kernel support of radius h around pi, tested in
 * float as the interaction loops do?
 *
 * The loops compute the distance from positions relative to the origin of
 * one of the two cells at the level of the visit, the other particle moved
 * by the periodic shift. The test is evaluated in all candidate frames and
 * in double; a disagreement is counted in #n_ambiguous_pairs. The verdict in
 * the first frame is returned. A pair beyond the support by 1 percent is out
 * in every frame.
 */
static int within_kernel(const struct part *pi, const struct part *pj,
                         const float h, struct pair_frames *f) {
  double r2d = 0.;
  for (int k = 0; k < 3; k++) {
    const double d = pi->x[k] - (pj->x[k] + f->image[k]);
    r2d += d * d;
  }
  if (r2d > 1.01 * h * h * kernel_gamma2) return 0;

  f->count = 0;
  for (const struct cell *c = f->b; c != NULL;
       c = c->split ? progeny_of(c, pj) : NULL) {
    f->loc[f->count] = c->loc;
    f->of_second[f->count++] = 1;
  }
  for (const struct cell *c = f->a; c != NULL;
       c = c->split ? progeny_of(c, pi) : NULL) {
    f->loc[f->count] = c->loc;
    f->of_second[f->count++] = 0;
  }

  const float hg2 = h * h * kernel_gamma2;
  int verdict = -1, ambiguous = 0;
  for (int n = 0; n < f->count; n++) {
    float dx[3];
    for (int k = 0; k < 3; k++) {
      const double shift_i = f->of_second[n] ? f->image[k] : 0.;
      const double shift_j = f->of_second[n] ? 0. : f->image[k];
      const float xi = pi->x[k] - (f->loc[n][k] + shift_i);
      const float xj = pj->x[k] - (f->loc[n][k] - shift_j);
      dx[k] = xj - xi;
    }
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    const int v = r2 < hg2;
    if (verdict < 0) verdict = v;
    if (v != verdict) ambiguous = 1;
  }
  if ((r2d < (double)h * h * kernel_gamma2) != verdict) ambiguous = 1;
  n_ambiguous_pairs += ambiguous;
  return verdict;
}

/**
 * @brief Add the brute-force number of interacting pairs of the particles of
 * a to the reference counts: one per pair with r < gamma max(hi, hj) and at
 * least one active particle.
 *
 * @return The number of pairs with both particles active in different depth
 * bands.
 */
static int brute_force_count(const struct cell *a, const struct cell *b,
                             int *ref_a, int *ref_b, const struct engine *e) {
  int mixed_band_pairs = 0;
  struct pair_frames frames;
  pair_frames(a, b, e->s, &frames);
  for (int i = 0; i < a->hydro.count; i++) {
    const struct part *pi = &a->hydro.parts[i];
    const int j_start = (a == b) ? i + 1 : 0;
    for (int j = j_start; j < b->hydro.count; j++) {
      const struct part *pj = &b->hydro.parts[j];
      const int active_i = part_is_active(pi, e);
      const int active_j = part_is_active(pj, e);
      if (within_kernel(pi, pj, max(pi->h, pj->h), &frames) &&
          (active_i || active_j)) {
        ref_a[i]++;
        ref_b[j]++;
        if (active_i && active_j && pi->depth_h != pj->depth_h)
          mixed_band_pairs++;
      }
    }
  }
  return mixed_band_pairs;
}

/**
 * @brief Is the cell of the local particle the first of the pair at the
 * level where the recursion of the pair task visits the pair?
 *
 * This follows runner_dosub_pair2: the pair is visited at the first level
 * where a cell is not split or too small to recurse, or where max(hi, hj)
 * reaches the split h of the cell pair.
 */
static int local_cell_is_first(const struct cell *local,
                               const struct cell *foreign,
                               const struct part *pi, const struct part *pj,
                               struct space *s) {
  const struct cell *a = local, *b = foreign;
  while (1) {
    struct cell *first = (struct cell *)a, *second = (struct cell *)b;
    double shift[3];
    space_getsid_and_swap_cells(s, &first, &second, shift);
    if (!a->split || !b->split ||
        a->hydro.count < space_recurse_size_pair_hydro ||
        b->hydro.count < space_recurse_size_pair_hydro ||
        max(pi->h, pj->h) >= cell_pair_hydro_h_split(first, second))
      return first == a;
    a = progeny_of(a, pi);
    b = progeny_of(b, pj);
  }
}

/**
 * @brief Predict the missed exchanges of the local cell of a pair with a
 * foreign cell.
 *
 * A local inactive particle is reached only from its own side, which loops
 * over its neighbours within its own kernel if its cell is the first of the
 * pair after the canonical swap, and over those beyond the kernel of the
 * neighbour otherwise. It therefore misses every foreign active neighbour
 * whose kernel reaches it, except those within its own kernel when its cell
 * is first.
 */
static void predict_foreign_gap(const struct cell *local,
                                const struct cell *foreign, int *miss,
                                const struct engine *e) {
  struct pair_frames frames;
  pair_frames(local, foreign, e->s, &frames);
  for (int i = 0; i < local->hydro.count; i++) {
    const struct part *pi = &local->hydro.parts[i];
    if (part_is_active(pi, e)) continue;
    for (int j = 0; j < foreign->hydro.count; j++) {
      const struct part *pj = &foreign->hydro.parts[j];
      if (!part_is_active(pj, e)) continue;
      if (!within_kernel(pi, pj, pj->h, &frames)) continue;
      if (local_cell_is_first(local, foreign, pi, pj, e->s) &&
          within_kernel(pi, pj, pi->h, &frames))
        continue;
      miss[i]++;
    }
  }
}

/**
 * @brief Move every particle of a sorted cell by less than @p d, and record
 * the move in the drift bounds of the cell and of its progeny.
 */
static void move_after_sort(struct cell *c, const double d) {
  if (c == c->top)
    for (int n = 0; n < c->hydro.count; n++)
      for (int k = 0; k < 3; k++)
        c->hydro.parts[n].x[k] += random_uniform(-d, d) / sqrt(3.);
  c->hydro.dx_max_part = d;
  c->hydro.dx_max_sort = d;
  c->hydro.dx_outside_max = d;
  if (c->split)
    for (int k = 0; k < 8; k++) move_after_sort(c->progeny[k], d);
}

/**
 * @brief Compare the exchange counters of a local cell with the reference
 * minus the predicted misses, in both passes.
 *
 * @return The number of (particle, pass) cases that violate it.
 */
static int count_mismatches(const struct cell *c, const int *ref,
                            const int *miss, int *doubles, int *drops) {
  int bad = 0;
  for (int k = 0; k < c->hydro.count; k++) {
    for (int pass = 0; pass < 2; pass++) {
      const int n = c->hydro.parts[k].chemistry_data.flux_exchange_count[pass];
      if (n > ref[k]) (*doubles)++;
      if (n < ref[k] - miss[k]) (*drops)++;
      if (n != ref[k] - miss[k]) bad++;
    }
  }
  return bad;
}

/**
 * @brief Count the particles whose prep and apply counters differ.
 */
static int count_unequal_passes(const struct cell *c) {
  int unequal = 0;
  for (int k = 0; k < c->hydro.count; k++) {
    const int *n = c->hydro.parts[k].chemistry_data.flux_exchange_count;
    if (n[0] != n[1]) unequal++;
  }
  return unequal;
}

/**
 * @brief Run the force and the prep loop on a cell pair (or a single cell
 * if cj is NULL) and compare the counters with the reference.
 *
 * With a foreign cj, the pair is run once with each cell local. Each run
 * checks the counts of its local cell, and the equality of the prep and
 * apply counts of every particle.
 *
 * @param expect_misses Must the prediction hold missed exchanges?
 * @param expect_mixed Must the pairs hold active pairs in different depth
 * bands?
 * @param d_move The largest move of the particles after the sort.
 * @return 1 on success.
 */
static int run_case(const char *label, struct runner *r, struct engine *e,
                    struct cell *ci, struct cell *cj, const int foreign,
                    const int expect_misses, const int expect_mixed,
                    const double d_move) {

  runner_do_hydro_sort(r, ci, 0x1FFF, 0, 0, 0, 0);
  if (cj != NULL) runner_do_hydro_sort(r, cj, 0x1FFF, 0, 0, 0, 0);
  if (d_move > 0.) {
    move_after_sort(ci, d_move);
    if (cj != NULL) move_after_sort(cj, d_move);
  }

  const int ambiguous_before = n_ambiguous_pairs;
  int *ref_i = calloc(ci->hydro.count, sizeof(int));
  int *ref_j = cj ? calloc(cj->hydro.count, sizeof(int)) : NULL;
  if (ref_i == NULL || (cj && ref_j == NULL))
    error("Could not allocate the reference counts.");
  const int mixed = cj ? brute_force_count(ci, cj, ref_i, ref_j, e)
                       : brute_force_count(ci, ci, ref_i, ref_i, e);

  int ref_total = 0;
  for (int k = 0; k < ci->hydro.count; k++) ref_total += ref_i[k];
  if (cj)
    for (int k = 0; k < cj->hydro.count; k++) ref_total += ref_j[k];

  /* Predicted misses, per rank for a foreign pair. */
  int *miss_i[2], *miss_j[2];
  for (int rank = 0; rank < 2; rank++) {
    miss_i[rank] = calloc(ci->hydro.count, sizeof(int));
    miss_j[rank] = cj ? calloc(cj->hydro.count, sizeof(int)) : NULL;
    if (miss_i[rank] == NULL || (cj && miss_j[rank] == NULL))
      error("Could not allocate the predicted misses.");
  }
  if (cj && foreign && foreign_gap_exists) {
    predict_foreign_gap(ci, cj, miss_i[0], e);
    predict_foreign_gap(cj, ci, miss_j[1], e);
  }

  int predicted = 0;
  for (int rank = 0; rank < 2; rank++) {
    for (int k = 0; k < ci->hydro.count; k++) predicted += miss_i[rank][k];
    if (cj)
      for (int k = 0; k < cj->hydro.count; k++) predicted += miss_j[rank][k];
  }

  int bad = 0, unequal = 0, doubles = 0, drops = 0;
  const int n_ranks = foreign ? 2 : 1;
  for (int rank = 0; rank < n_ranks; rank++) {
    e->nodeID = rank;
    reset_counts(ci);
    if (cj) reset_counts(cj);

    if (cj) {
      runner_dosub_pair2_chemistry_flux_prep(r, ci, cj, /*h_hi=*/FLT_MAX,
                                             /*gettimer=*/0);
      runner_dosub_pair2_force(r, ci, cj, /*h_hi=*/FLT_MAX, /*gettimer=*/0);
    } else {
      runner_dosub_self2_chemistry_flux_prep(r, ci, /*h_hi=*/FLT_MAX,
                                             /*gettimer=*/0);
      runner_dosub_self2_force(r, ci, /*h_hi=*/FLT_MAX, /*gettimer=*/0);
    }

    unequal += count_unequal_passes(ci);
    if (cj) unequal += count_unequal_passes(cj);
    if (ci->nodeID == rank)
      bad += count_mismatches(ci, ref_i, miss_i[rank], &doubles, &drops);
    if (cj && cj->nodeID == rank)
      bad += count_mismatches(cj, ref_j, miss_j[rank], &doubles, &drops);
  }
  e->nodeID = 0;

  const int ambiguous = n_ambiguous_pairs - ambiguous_before;
  const int reach_ok = (expect_misses ? (predicted > 0) : (predicted == 0)) &&
                       (!expect_mixed || mixed > 0);
  message(
      "%-42s reference exchanges=%6d mixed-band pairs=%5d predicted "
      "misses=%4d doubles=%d unpredicted drops=%d mismatches=%d "
      "unequal passes=%d ambiguous pairs=%d",
      label, ref_total, mixed, predicted, doubles, drops, bad, unequal,
      ambiguous);

  for (int rank = 0; rank < 2; rank++) {
    free(miss_i[rank]);
    free(miss_j[rank]);
  }
  free(ref_i);
  free(ref_j);
  return bad == 0 && unequal == 0 && ambiguous == 0 && reach_ok;
}

int main(int argc, char *argv[]) {

  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  const unsigned int seed = (argc > 1) ? atoi(argv[1]) : 1234;
  message("Seed used for RNG: %u", seed);
  srand(seed);

  struct space space;
  bzero(&space, sizeof(struct space));
  space.periodic = 0;
  space.dim[0] = 4.;
  space.dim[1] = 4.;
  space.dim[2] = 4.;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct hydro_props hydro_props;
  hydro_props_init_no_hydro(&hydro_props);
  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  struct phys_const phys_const;
  bzero(&phys_const, sizeof(struct phys_const));
  phys_const.const_vacuum_permeability = 1.;

  struct chemistry_global_data chem_data;
  bzero(&chem_data, sizeof(struct chemistry_global_data));
  chem_data.diffusion_coefficient = 1.f;
  chem_data.diffusion_mode = isotropic_constant;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  chem_data.tau = 1.;
  chem_data.relaxation_time_mode = constant_mode;
  chem_data.hyperbolic_limiter_scope = limiter_all_components;
#endif
  chem_data.hll_riemann_solver_psi = 0.1f;
  chem_data.C_CFL_chemistry = 0.5f;
  chem_data.flux_limiter_noise_gate = 1e-15;
  chem_data.flux_limiter_safety = 0.5;
  chem_data.flux_limiter_sink_stability = 0.25;
  chem_data.flux_limiter_startup = 0.1;
  chem_data.use_positivity_limiter = 1;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.ti_current = TI_CURRENT;
  engine.max_active_bin = MAX_ACTIVE_BIN;
  engine.time_base = 1e-3;
  engine.nodeID = 0;
  engine.cosmology = &cosmo;
  engine.hydro_properties = &hydro_props;
  engine.pressure_floor_props = &pressure_floor;
  engine.physical_constants = &phys_const;
  engine.chemistry = &chem_data;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  const double loc_i[3] = {0., 0., 0.};
  const double offsets[3][3] = {{1., 0., 0.}, {1., 1., 0.}, {1., 1., 1.}};
  const char *pair_names[3] = {"face", "edge", "corner"};

  int pass = 1;
  char label[64];

  for (int split = 0; split < 2; split++) {
    struct cell *c = make_cell(loc_i, split, 0.02, 0, 0.5, 0, &engine);
    sprintf(label, "self, %s", split ? "split" : "unsplit");
    pass &= run_case(label, &runner, &engine, c, NULL, 0, 0, split, 0.);
    free_cell(c);
  }

  for (int o = 0; o < 3; o++) {
    for (int variant = 0; variant < 4; variant++) {
      const int split = (variant >= 2);
      const int foreign = (variant == 1);
      const int top_band_active = (variant == 2);
      const int expect_misses = foreign && foreign_gap_exists;

      struct cell *ci =
          make_cell(loc_i, split, 0.02, 0, 0.5, top_band_active, &engine);
      struct cell *cj = make_cell(offsets[o], split, 0.02, foreign, 0.5,
                                  top_band_active, &engine);
      sprintf(label, "pair %s, %s%s%s", pair_names[o],
              split ? "split" : "unsplit", foreign ? ", foreign" : "",
              top_band_active ? ", top band active" : "");
      pass &= run_case(label, &runner, &engine, ci, cj, foreign, expect_misses,
                       split, 0.);
      free_cell(ci);
      free_cell(cj);
    }
  }

  /* The cases below have no margins, in a periodic space. */
  space.periodic = 1;
  if (8 * PARTS_PER_LEAF_DEPTH2 < space_recurse_size_pair_hydro ||
      8 * PARTS_PER_LEAF_DEPTH2 < space_recurse_size_self_hydro)
    error("The tasks of a cell of depth 2 do not recurse to its leaves.");

  struct cell *c = make_cell(loc_i, 2, 0., 0, 0.5, 0, &engine);
  pass &= run_case("self, depth 2", &runner, &engine, c, NULL, 0, 0, 1, 0.);
  free_cell(c);

  c = make_cell(loc_i, 2, 0., 0, 0.5, 0, &engine);
  pass &= run_case("self, depth 2, moved after sort", &runner, &engine, c, NULL,
                   0, 0, 1, 0.02);
  free_cell(c);

  /* Pair cases: depth of ci and cj, rank of cj, fraction of active particles
     in cj, position of cj and move after the sort. */
  const struct {
    const char *label;
    int depth_i, depth_j, foreign;
    double frac_active_j, loc_j[3], d_move;
  } cases[] = {
      {"pair face, split and unsplit", 1, 0, 0, 0.5, {1., 0., 0.}, 0.},
      {"pair face, split, inactive cj", 1, 1, 0, 0., {1., 0., 0.}, 0.},
      {"pair face, split, foreign", 1, 1, 1, 0.5, {1., 0., 0.}, 0.},
      {"pair face, depth 2", 2, 2, 0, 0.5, {1., 0., 0.}, 0.},
      {"pair face, depth 2, foreign", 2, 2, 1, 0.5, {1., 0., 0.}, 0.},
      {"pair face, split, periodic", 1, 1, 0, 0.5, {3., 0., 0.}, 0.},
      {"pair face, split, moved after sort", 1, 1, 0, 0.5, {1., 0., 0.}, 0.02},
      {"pair corner, depth 2, periodic, moved, foreign",
       2,
       2,
       1,
       0.5,
       {3., 3., 3.},
       0.02}};
  const int n_cases = sizeof(cases) / sizeof(cases[0]);

  for (int n = 0; n < n_cases; n++) {
    struct cell *ci =
        make_cell(loc_i, cases[n].depth_i, 0., 0, 0.5, 0, &engine);
    struct cell *cj =
        make_cell(cases[n].loc_j, cases[n].depth_j, 0., cases[n].foreign,
                  cases[n].frac_active_j, 0, &engine);
    const int expect_mixed = cases[n].frac_active_j > 0.;
    pass &= run_case(cases[n].label, &runner, &engine, ci, cj, cases[n].foreign,
                     cases[n].foreign && foreign_gap_exists, expect_mixed,
                     cases[n].d_move);
    free_cell(ci);
    free_cell(cj);
  }

  if (!pass) error("testChemistryFluxExchangeOnce: FAILED");
  message("testChemistryFluxExchangeOnce: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
