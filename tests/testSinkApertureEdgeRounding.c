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

/* The fixed-aperture gas-gas loop and the brute-force count check of sink
 * formation must classify a pair at the edge of the aperture identically.
 * The check calls sink_formation_gas_loop_r2() and
 * sink_formation_gas_loop_in_aperture(); this test places isolated pairs
 * within a few hundred float ULP of r_cut^2 and requires every loop variant
 * (sorted and naive pair, optimised and naive self, with and without a
 * periodic wrap, all activity combinations) to count exactly the pairs those
 * functions accept. */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

#define NODE_ID 0

/* Top-level cell width of the sink box (32 x 0.027248 / 3) and its aperture. */
#define CELL_WIDTH 0.29065
#define R_CUT 5e-3f

/* Half-width of the band of r^2 / r_cut^2 - 1 in which pairs are placed. */
#define EDGE_BAND 2.5e-5

/* Active particles are in bin 1, passive ones in bin 4, at ti_current = 8. */
#define TI_CURRENT 8
#define BIN_ACTIVE 1
#define BIN_PASSIVE 4

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_hydro_aperture_edge_rounding(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H, const int with_self_gravity, const struct cosmology *cosmo,
    const struct sink_props *sink_props) {
  pi->density.wcount += 1.f;
}

__attribute__((always_inline)) INLINE static void
runner_iact_hydro_aperture_edge_rounding(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const int with_self_gravity, const struct cosmology *cosmo,
    const struct sink_props *sink_props) {
  pi->density.wcount += 1.f;
  pj->density.wcount += 1.f;
}

#include "space_getsid.h"
#define FUNCTION edge_rounding
#define FUNCTION_TASK_LOOP TASK_LOOP_PREP_SINK_FORMATION
#include "runner_doiact_functions_hydro_aperture.h"
#include "runner_doiact_undef.h"

static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;

/** @brief Uniform deviate in [0, 1) from a fixed-seed xorshift generator. */
static double rng_uniform(void) {
  rng_state ^= rng_state << 13;
  rng_state ^= rng_state >> 7;
  rng_state ^= rng_state << 17;
  return (double)(rng_state >> 11) * (1.0 / 9007199254740992.0);
}

/**
 * @brief Allocate an empty cell of width CELL_WIDTH at @p loc with room for
 * @p capacity particles.
 */
static struct cell *make_empty_cell(const double loc[3], const int capacity) {
  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, sizeof(struct cell)) != 0)
    error("Couldn't allocate cell.");
  bzero(c, sizeof(struct cell));
  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     capacity * sizeof(struct part)) != 0)
    error("Couldn't allocate particles.");
  bzero(c->hydro.parts, capacity * sizeof(struct part));
  for (int k = 0; k < 3; k++) {
    c->loc[k] = loc[k];
    c->width[k] = CELL_WIDTH;
  }
  c->dmin = CELL_WIDTH;
  c->hydro.super = c;
  c->hydro.ti_old_part = TI_CURRENT;
  c->hydro.ti_end_min = TI_CURRENT;
  c->grav.ti_old_part = TI_CURRENT;
  c->grav.ti_end_min = TI_CURRENT;
  c->nodeID = NODE_ID;
  return c;
}

/** @brief Append a particle at @p x to cell @p c. */
static struct part *add_part(struct cell *c, const double x[3],
                             const long long id, const int active) {
  struct part *p = &c->hydro.parts[c->hydro.count++];
  p->x[0] = x[0];
  p->x[1] = x[1];
  p->x[2] = x[2];
  p->h = R_CUT;
  p->id = id;
  p->time_bin = active ? BIN_ACTIVE : BIN_PASSIVE;
  hydro_set_mass(p, 1.f);
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_drift = TI_CURRENT;
  p->ti_kick = TI_CURRENT;
#endif
  if (p->h > c->hydro.h_max) c->hydro.h_max = p->h;
  c->hydro.h_max_active = c->hydro.h_max;
  return p;
}

static void free_cell(struct cell *c) {
  free(c->hydro.parts);
  free(c->hydro.sort);
  free(c);
}

/** @brief Wrap a coordinate into [0, dim). */
static double wrap(double x, const double dim) {
  if (x >= dim) x -= dim;
  if (x < 0.) x += dim;
  return x;
}

/**
 * @brief Separation vector of length r_cut * sqrt(1 + delta), with delta
 * uniform in [-EDGE_BAND, EDGE_BAND] and a component of at least 0.2 r_cut
 * along +axis (zero axis: any direction).
 */
static void edge_vector(const int axis, double d[3]) {
  double u[3], norm;
  do {
    for (int k = 0; k < 3; k++) u[k] = 2. * rng_uniform() - 1.;
    norm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  } while (norm < 0.1 || norm > 1. || (axis >= 0 && u[axis] / norm < 0.2));
  const double delta = EDGE_BAND * (2. * rng_uniform() - 1.);
  const double r = (double)R_CUT * sqrt(1. + delta);
  for (int k = 0; k < 3; k++) d[k] = r * u[k] / norm;
}

struct tally {
  int counted;
  int disagree;
  int inside;
  int outside;
  int one_ulp;
};

/**
 * @brief Compare every active particle's loop count with the shared
 * classification of its single candidate neighbour.
 */
static void compare(struct part *const *first, struct part *const *second,
                    const int n_pairs, const double dim[3], const int periodic,
                    const char *label, struct tally *t) {

  const double r_cut2 = (double)R_CUT * (double)R_CUT;
  const double ulp_r_cut2 = nextafterf(R_CUT * R_CUT, 1.f) - R_CUT * R_CUT;

  for (int i = 0; i < n_pairs; i++) {
    const struct part *a = first[i];
    const struct part *b = second[i];
    float dx_ab[3], dx_ba[3];
    const double r2_ab =
        sink_formation_gas_loop_r2(a->x, b->x, dim, periodic, dx_ab);
    const double r2_ba =
        sink_formation_gas_loop_r2(b->x, a->x, dim, periodic, dx_ba);
    if (r2_ab != r2_ba)
      error("%s: r2 depends on the particle order (%.17g, %.17g)", label, r2_ab,
            r2_ba);
    const int expect = sink_formation_gas_loop_in_aperture(r2_ab, R_CUT);
    if (expect)
      t->inside++;
    else
      t->outside++;
    if (fabs(r2_ab - r_cut2) <= ulp_r_cut2) t->one_ulp++;

    const struct part *pair[2] = {a, b};
    for (int k = 0; k < 2; k++) {
      const struct part *p = pair[k];
      if (p->time_bin != BIN_ACTIVE) continue;
      t->counted++;
      if ((int)p->density.wcount != expect) {
        t->disagree++;
        if (t->disagree <= 5)
          message(
              "%s: particle %lld counted %d, shared test says %d "
              "(r2 / r_cut2 - 1 = %.3e)",
              label, p->id, (int)p->density.wcount, expect,
              r2_ab / r_cut2 - 1.);
      }
    }
  }
}

static void reset_counts(struct cell *c) {
  for (int i = 0; i < c->hydro.count; i++)
    c->hydro.parts[i].density.wcount = 0.f;
}

/**
 * @brief Pairs straddling the face of @p ci and @p cj normal to @p axis.
 *
 * @param face Coordinate of the face along @p axis in the frame of ci.
 */
static void run_pair_case(struct runner *r, const double loc_i[3],
                          const double loc_j[3], const int axis,
                          const double face, const char *label,
                          long long *next_id, struct tally *t_sorted,
                          struct tally *t_naive) {

  const struct space *s = r->e->s;
  const double spacing = 3. * R_CUT;
  const int n_side = (int)((CELL_WIDTH - 4. * R_CUT) / spacing);
  const int capacity = n_side * n_side;

  struct cell *ci = make_empty_cell(loc_i, capacity);
  struct cell *cj = make_empty_cell(loc_j, capacity);
  struct part **first = malloc(capacity * sizeof(struct part *));
  struct part **second = malloc(capacity * sizeof(struct part *));
  const int u_axis = (axis + 1) % 3, v_axis = (axis + 2) % 3;

  int n_pairs = 0;
  for (int iu = 0; iu < n_side; iu++) {
    for (int iv = 0; iv < n_side; iv++) {
      double d[3];
      edge_vector(axis, d);
      const double f = 0.1 + 0.8 * rng_uniform();
      double xa[3], xb[3];
      xa[axis] = face - f * d[axis];
      xa[u_axis] = loc_i[u_axis] + 2. * R_CUT + iu * spacing;
      xa[v_axis] = loc_i[v_axis] + 2. * R_CUT + iv * spacing;
      for (int k = 0; k < 3; k++) xb[k] = wrap(xa[k] + d[k], s->dim[k]);
      const int mode = (int)(3. * rng_uniform());
      first[n_pairs] = add_part(ci, xa, ++(*next_id), mode != 2);
      second[n_pairs] = add_part(cj, xb, ++(*next_id), mode != 1);
      n_pairs++;
    }
  }

  runner_do_drift_part(r, ci, 0);
  runner_do_drift_part(r, cj, 0);
  runner_do_hydro_sort(r, ci, 0x1FFF, 0, 0, 0, 0);
  runner_do_hydro_sort(r, cj, 0x1FFF, 0, 0, 0, 0);

  reset_counts(ci);
  reset_counts(cj);
  runner_dopair1_branch_hydro_aperture_edge_rounding(r, ci, cj, R_CUT);
  compare(first, second, n_pairs, s->dim, s->periodic, label, t_sorted);

  reset_counts(ci);
  reset_counts(cj);
  runner_dopair1_naive_hydro_aperture_edge_rounding(r, ci, cj, R_CUT);
  compare(first, second, n_pairs, s->dim, s->periodic, label, t_naive);

  free(first);
  free(second);
  free_cell(ci);
  free_cell(cj);
}

/** @brief Pairs inside one cell, in random array order. */
static void run_self_case(struct runner *r, long long *next_id,
                          struct tally *t_opt, struct tally *t_naive) {

  const struct space *s = r->e->s;
  const double spacing = 3. * R_CUT;
  const int n_side = (int)((CELL_WIDTH - 4. * R_CUT) / spacing);
  const int n_max = n_side * n_side * n_side;
  const double loc[3] = {CELL_WIDTH, CELL_WIDTH, CELL_WIDTH};

  struct cell *c = make_empty_cell(loc, 2 * n_max);
  struct part **first = malloc(n_max * sizeof(struct part *));
  struct part **second = malloc(n_max * sizeof(struct part *));

  int n_pairs = 0;
  for (int ix = 0; ix < n_side; ix++) {
    for (int iy = 0; iy < n_side; iy++) {
      for (int iz = 0; iz < n_side; iz++) {
        double d[3];
        edge_vector(-1, d);
        double xa[3], xb[3];
        const int idx[3] = {ix, iy, iz};
        for (int k = 0; k < 3; k++) {
          xa[k] = loc[k] + 2. * R_CUT + idx[k] * spacing;
          xb[k] = xa[k] + d[k];
        }
        const int mode = (int)(3. * rng_uniform());
        const int a_first = rng_uniform() < 0.5;
        struct part *pa, *pb;
        if (a_first) {
          pa = add_part(c, xa, ++(*next_id), mode != 2);
          pb = add_part(c, xb, ++(*next_id), mode != 1);
        } else {
          pb = add_part(c, xb, ++(*next_id), mode != 1);
          pa = add_part(c, xa, ++(*next_id), mode != 2);
        }
        first[n_pairs] = pa;
        second[n_pairs] = pb;
        n_pairs++;
      }
    }
  }

  runner_do_drift_part(r, c, 0);

  reset_counts(c);
  runner_doself1_branch_hydro_aperture_edge_rounding(r, c, R_CUT);
  compare(first, second, n_pairs, s->dim, s->periodic, "self", t_opt);

  reset_counts(c);
  runner_doself1_naive_hydro_aperture_edge_rounding(r, c, R_CUT);
  compare(first, second, n_pairs, s->dim, s->periodic, "self naive", t_naive);

  free(first);
  free(second);
  free_cell(c);
}

static void report(const char *label, const struct tally *t) {
  message(
      "%-24s active particles %6d, disagreements %4d, pairs inside %5d, "
      "outside %5d, within 1 float ULP of r_cut^2 %4d",
      label, t->counted, t->disagree, t->inside, t->outside, t->one_ulp);
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  struct space space;
  bzero(&space, sizeof(struct space));
  space.periodic = 1;
  space.dim[0] = 3. * CELL_WIDTH;
  space.dim[1] = 3. * CELL_WIDTH;
  space.dim[2] = 3. * CELL_WIDTH;

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  hp.h_max = FLT_MAX;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = TI_CURRENT;
  engine.max_active_bin = 3;
  engine.hydro_properties = &hp;
  engine.nodeID = NODE_ID;

  struct phys_const prog_const;
  bzero(&prog_const, sizeof(struct phys_const));
  prog_const.const_vacuum_permeability = 1.0;
  engine.physical_constants = &prog_const;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct sink_props sink_props;
  bzero(&sink_props, sizeof(struct sink_props));
  engine.sink_properties = &sink_props;

  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  engine.pressure_floor_props = &pressure_floor;

  struct lightcone_array_props lightcone_array_properties;
  lightcone_array_properties.nr_lightcones = 0;
  engine.lightcone_array_properties = &lightcone_array_properties;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  long long next_id = 0;
  struct tally sorted = {0}, naive = {0}, self = {0}, self_naive = {0};
  const double W = CELL_WIDTH;

  /* Faces inside the box, then the face across the periodic boundary. */
  for (int axis = 0; axis < 3; axis++) {
    const double loc_i[3] = {0., 0., 0.};
    double loc_j[3] = {0., 0., 0.};
    loc_j[axis] = W;
    run_pair_case(&runner, loc_i, loc_j, axis, W, "pair", &next_id, &sorted,
                  &naive);
  }
  for (int axis = 0; axis < 3; axis++) {
    double loc_i[3] = {0., 0., 0.};
    const double loc_j[3] = {0., 0., 0.};
    loc_i[axis] = 2. * W;
    run_pair_case(&runner, loc_i, loc_j, axis, 3. * W, "pair periodic",
                  &next_id, &sorted, &naive);
  }
  run_self_case(&runner, &next_id, &self, &self_naive);

  report("pair, sorted", &sorted);
  report("pair, naive", &naive);
  report("self, optimised", &self);
  report("self, naive", &self_naive);

  const struct tally *all[4] = {&sorted, &naive, &self, &self_naive};
  int failed = 0;
  for (int k = 0; k < 4; k++) {
    if (all[k]->disagree > 0) failed = 1;
    if (all[k]->inside == 0 || all[k]->outside == 0 || all[k]->one_ulp == 0) {
      message("No pair on one side of the edge, or none within 1 ULP of it.");
      failed = 1;
    }
  }

  if (failed) {
    message("FAIL: a loop and the shared aperture test disagree.");
    return 1;
  }
  message("PASS: every loop counts exactly the pairs the shared test accepts.");
  return 0;
}
