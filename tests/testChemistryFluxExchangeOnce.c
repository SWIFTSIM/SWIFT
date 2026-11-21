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
 * and at least one active particle. Cells are optionally split once, with h
 * in two depth bands, so that the recursion of the sub-cell tasks and the
 * ownership rule of mixed-band pairs are exercised.
 *
 * Two configurations miss exchanges in the pair enumeration itself. There,
 * each particle's count must equal the reference minus the exact number of
 * misses predicted for it (predict_split_pair_gap(), predict_foreign_gap()).
 * The prep and apply counts of every particle must be equal.
 */

#if (defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) ||             \
     defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)) && \
    defined(SWIFT_CHEMISTRY_DEBUG_CHECKS)

#define TI_CURRENT 8
#define ACTIVE_BIN 1
#define INACTIVE_BIN 4
#define MAX_ACTIVE_BIN 3
#define PARTS_PER_OCTANT 27

void runner_dosub_self2_force(struct runner *r, struct cell *c,
                              int recurse_below_h_max, const int gettimer);
void runner_dosub_pair2_force(struct runner *r, struct cell *ci,
                              struct cell *cj, int recurse_below_h_max,
                              const int gettimer);
void runner_dosub_self2_chemistry_flux_prep(struct runner *r, struct cell *c,
                                            int recurse_below_h_max,
                                            const int gettimer);
void runner_dosub_pair2_chemistry_flux_prep(struct runner *r, struct cell *ci,
                                            struct cell *cj,
                                            int recurse_below_h_max,
                                            const int gettimer);

/**
 * @brief Set the hydro and chemistry state of a particle for a force loop.
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

  /* Inactive particles carry a negative flux time-step. */
  ch->flux.dt = active ? 1e-3f : -1.f;
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
 * @brief Build a unit cell of 8 * #PARTS_PER_OCTANT particles, split once if
 * requested.
 *
 * With split, h is drawn from the band of the top cell or from the band of
 * its progeny, so the particles sit at two depths. Without split, h is in
 * the band of the top cell.
 */
static struct cell *make_cell(const double loc[3], const int split,
                              const int node_id, const double frac_active,
                              const int top_band_active,
                              const struct engine *e) {

  const int count = 8 * PARTS_PER_OCTANT;
  const double width = 1.;
  const double h_top_min = width * 0.5 / kernel_gamma;
  const double h_top_max = width / kernel_gamma;

  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, 9 * sizeof(struct cell)) != 0)
    error("Could not allocate the cells.");
  bzero(c, 9 * sizeof(struct cell));
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
  c->split = split;

  for (int k = 0; k < 8; k++) {
    const double sub_loc[3] = {loc[0] + 0.5 * width * ((k >> 2) & 1),
                               loc[1] + 0.5 * width * ((k >> 1) & 1),
                               loc[2] + 0.5 * width * (k & 1)};
    for (int n = 0; n < PARTS_PER_OCTANT; n++) {
      const int ind = k * PARTS_PER_OCTANT + n;
      struct part *p = &c->hydro.parts[ind];
      for (int i = 0; i < 3; i++)
        p->x[i] = sub_loc[i] + 0.5 * width * random_uniform(0.02, 0.98);
      const int in_top_band = !split || random_uniform(0., 1.) < 0.3;
      if (in_top_band)
        p->h = random_uniform(1.02 * h_top_min, 0.98 * h_top_max);
      else
        p->h = random_uniform(0.52 * h_top_min, 0.98 * h_top_min);
      p->id = ind + 1 + 1000 * node_id;
      const int active = (top_band_active && in_top_band) ||
                         random_uniform(0., 1.) < frac_active;
      init_particle(p, &c->hydro.xparts[ind], 1.0 / count, active, e);
    }

    if (split) {
      struct cell *cp = &c[1 + k];
      for (int i = 0; i < 3; i++) {
        cp->loc[i] = sub_loc[i];
        cp->width[i] = 0.5 * width;
      }
      cp->hydro.parts = c->hydro.parts + k * PARTS_PER_OCTANT;
      cp->hydro.xparts = c->hydro.xparts + k * PARTS_PER_OCTANT;
      cp->hydro.count = PARTS_PER_OCTANT;
      cp->depth = 1;
      cp->parent = c;
      cp->top = c;
      cp->hydro.super = c;
      cp->nodeID = node_id;
      c->progeny[k] = cp;
    }
  }

  /* Foreign cells hold no xparts. */
  if (node_id != e->nodeID) {
    free(c->hydro.xparts);
    c->hydro.xparts = NULL;
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) c->progeny[k]->hydro.xparts = NULL;
  }

  set_cell_properties(c, e);
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp == NULL) continue;
    set_cell_properties(cp, e);
    for (int n = 0; n < cp->hydro.count; n++)
      cell_set_part_h_depth(&cp->hydro.parts[n], cp);
  }
  if (!split)
    for (int n = 0; n < count; n++)
      cell_set_part_h_depth(&c->hydro.parts[n], c);

  return c;
}

static void free_cell(struct cell *c) {
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) cell_free_hydro_sorts(c->progeny[k]);
  cell_free_hydro_sorts(c);
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

/* The naive pair loop visits every pair from both sides, so it has no
   foreign gap. */
#ifdef SWIFT_USE_NAIVE_INTERACTIONS
static const int foreign_gap_exists = 0;
#else
static const int foreign_gap_exists = 1;
#endif

/**
 * @brief Origin of the progeny cell holding the particle of index k, or of
 * the cell itself if it is not split.
 */
static const double *octant_loc(const struct cell *c, const int k) {
  return c->split ? c->progeny[k / PARTS_PER_OCTANT]->loc : c->loc;
}

/**
 * @brief Is pj within the kernel support of radius h around pi, tested in
 * float as the interaction loops do?
 *
 * The loops compute the distance from positions relative to the origin of
 * one of the two cells at the level of the visit. The test is evaluated in
 * all candidate frames and in double; a disagreement is counted in
 * #n_ambiguous_pairs. The verdict in the first frame is returned.
 */
static int within_kernel(const struct part *pi, const struct part *pj,
                         const float h, const double *frames[4]) {
  const float hg2 = h * h * kernel_gamma2;
  int verdict = -1, ambiguous = 0;
  for (int f = 0; f < 4; f++) {
    float dx[3];
    for (int k = 0; k < 3; k++) {
      const float xi = pi->x[k] - frames[f][k];
      const float xj = pj->x[k] - frames[f][k];
      dx[k] = xj - xi;
    }
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    const int v = r2 < hg2;
    if (verdict < 0) verdict = v;
    if (v != verdict) ambiguous = 1;
  }
  double r2d = 0.;
  for (int k = 0; k < 3; k++)
    r2d += (pi->x[k] - pj->x[k]) * (pi->x[k] - pj->x[k]);
  if ((r2d < (double)h * h * kernel_gamma2) != verdict) ambiguous = 1;
  n_ambiguous_pairs += ambiguous;
  return verdict;
}

/**
 * @brief Fill the candidate frames of the pair (a[i], b[j]).
 */
static void pair_frames(const struct cell *a, const int i, const struct cell *b,
                        const int j, const double *frames[4]) {
  frames[0] = b->loc;
  frames[1] = a->loc;
  frames[2] = octant_loc(b, j);
  frames[3] = octant_loc(a, i);
}

/**
 * @brief Add the brute-force number of interacting pairs of the particles of
 * a to the reference counts: one per pair with r < gamma max(hi, hj) and at
 * least one active particle.
 *
 * @return The number of pairs with both particles active in different depth
 * bands, which only the ownership rule settles.
 */
static int brute_force_count(const struct cell *a, const struct cell *b,
                             int *ref_a, int *ref_b, const struct engine *e) {
  int mixed_band_pairs = 0;
  for (int i = 0; i < a->hydro.count; i++) {
    const struct part *pi = &a->hydro.parts[i];
    const int j_start = (a == b) ? i + 1 : 0;
    for (int j = j_start; j < b->hydro.count; j++) {
      const struct part *pj = &b->hydro.parts[j];
      const double *frames[4];
      pair_frames(a, i, b, j, frames);
      const int active_i = part_is_active(pi, e);
      const int active_j = part_is_active(pj, e);
      if (within_kernel(pi, pj, max(pi->h, pj->h), frames) &&
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
 * @brief Do the progeny holding a[i] and b[j] share a face, an edge or a
 * corner?
 */
static int octants_touch(const struct cell *a, const int i,
                         const struct cell *b, const int j) {
  const struct cell *pa = a->progeny[i / PARTS_PER_OCTANT];
  const struct cell *pb = b->progeny[j / PARTS_PER_OCTANT];
  for (int k = 0; k < 3; k++) {
    if (pa->loc[k] > pb->loc[k] + pb->width[k] + 1e-9) return 0;
    if (pb->loc[k] > pa->loc[k] + pa->width[k] + 1e-9) return 0;
  }
  return 1;
}

/**
 * @brief Predict the missed exchanges of a pair of split local cells.
 *
 * The sub-cell recursion only visits pairs of progeny that touch. A pair of
 * particles in progeny that do not touch is visited only at the top level,
 * where only the particles in the top depth band are taken. An inactive
 * top-band particle and an active progeny-band partner within its kernel
 * therefore never exchange: both miss one count.
 */
static void predict_split_pair_gap(const struct cell *ci, const struct cell *cj,
                                   int *miss_i, int *miss_j,
                                   const struct engine *e) {
  for (int i = 0; i < ci->hydro.count; i++) {
    const struct part *pi = &ci->hydro.parts[i];
    for (int j = 0; j < cj->hydro.count; j++) {
      const struct part *pj = &cj->hydro.parts[j];
      if (octants_touch(ci, i, cj, j)) continue;
      const int active_i = part_is_active(pi, e);
      const int active_j = part_is_active(pj, e);
      const int top_i = (pi->depth_h == 0) && (pj->depth_h == 1);
      const int top_j = (pj->depth_h == 0) && (pi->depth_h == 1);
      if (!(top_i && !active_i && active_j) &&
          !(top_j && !active_j && active_i))
        continue;
      const double *frames[4];
      pair_frames(ci, i, cj, j, frames);
      if (within_kernel(pi, pj, max(pi->h, pj->h), frames)) {
        miss_i[i]++;
        miss_j[j]++;
      }
    }
  }
}

/**
 * @brief Predict the missed exchanges of the local cell of a pair with a
 * foreign cell.
 *
 * A local inactive particle is reached only from its own side, which loops
 * over its neighbours within its own kernel, and only if the cell is the
 * first of the pair after the canonical swap. It therefore misses every
 * foreign active neighbour whose kernel reaches it, except those within its
 * own kernel when its cell is first.
 */
static void predict_foreign_gap(const struct cell *local,
                                const struct cell *foreign,
                                const int local_first, int *miss,
                                const struct engine *e) {
  for (int i = 0; i < local->hydro.count; i++) {
    const struct part *pi = &local->hydro.parts[i];
    if (part_is_active(pi, e)) continue;
    for (int j = 0; j < foreign->hydro.count; j++) {
      const struct part *pj = &foreign->hydro.parts[j];
      if (!part_is_active(pj, e)) continue;
      const double *frames[4];
      pair_frames(local, i, foreign, j, frames);
      if (!within_kernel(pi, pj, pj->h, frames)) continue;
      if (local_first && within_kernel(pi, pj, pi->h, frames)) continue;
      miss[i]++;
    }
  }
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
 * @return 1 on success.
 */
static int run_case(const char *label, struct runner *r, struct engine *e,
                    struct cell *ci, struct cell *cj, const int foreign,
                    const int expect_misses) {

  runner_do_hydro_sort(r, ci, 0x1FFF, 0, 0, 0, 0);
  if (cj != NULL) runner_do_hydro_sort(r, cj, 0x1FFF, 0, 0, 0, 0);

  const int ambiguous_before = n_ambiguous_pairs;
  int *ref_i = calloc(ci->hydro.count, sizeof(int));
  int *ref_j = cj ? calloc(cj->hydro.count, sizeof(int)) : NULL;
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
  }
  if (cj && foreign && foreign_gap_exists) {
    struct cell *first = ci, *second = cj;
    double shift[3];
    space_getsid_and_swap_cells(e->s, &first, &second, shift);
    predict_foreign_gap(ci, cj, first == ci, miss_i[0], e);
    predict_foreign_gap(cj, ci, first == cj, miss_j[1], e);
  } else if (cj && ci->split && cj->split) {
    predict_split_pair_gap(ci, cj, miss_i[0], miss_j[0], e);
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
      runner_dosub_pair2_chemistry_flux_prep(r, ci, cj, 0, 0);
      runner_dosub_pair2_force(r, ci, cj, 0, 0);
    } else {
      runner_dosub_self2_chemistry_flux_prep(r, ci, 0, 0);
      runner_dosub_self2_force(r, ci, 0, 0);
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
  const int reach_ok = expect_misses ? (predicted > 0) : (predicted == 0);
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
    struct cell *c = make_cell(loc_i, split, 0, 0.5, 0, &engine);
    sprintf(label, "self, %s", split ? "split" : "unsplit");
    pass &= run_case(label, &runner, &engine, c, NULL, 0, 0);
    free_cell(c);
  }

  for (int o = 0; o < 3; o++) {
    for (int variant = 0; variant < 4; variant++) {
      const int split = (variant >= 2);
      const int foreign = (variant == 1);
      const int top_band_active = (variant == 2);
      const int expect_misses =
          (foreign && foreign_gap_exists) || (variant == 3);

      struct cell *ci =
          make_cell(loc_i, split, 0, 0.5, top_band_active, &engine);
      struct cell *cj =
          make_cell(offsets[o], split, foreign, 0.5, top_band_active, &engine);
      sprintf(label, "pair %s, %s%s%s", pair_names[o],
              split ? "split" : "unsplit", foreign ? ", foreign" : "",
              top_band_active ? ", top band active" : "");
      pass &= run_case(label, &runner, &engine, ci, cj, foreign, expect_misses);
      free_cell(ci);
      free_cell(cj);
    }
  }

  /* Split foreign pairs are not run: see the TODO(MPI) in
     engine_maketasks.c. */

  if (!pass) error("testChemistryFluxExchangeOnce: FAILED");
  message("testChemistryFluxExchangeOnce: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
