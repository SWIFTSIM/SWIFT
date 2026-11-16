/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

/**
 * @file testSortShiftBound.c
 *
 * @brief Checks that c->hydro.dx_max_sort (and c->stars.dx_max_sort) remains
 * an upper bound on the particle motion since the sort of *every* sorted
 * direction, including directions added to a cell after its last reset.
 *
 * The scenario, applied with the real drift and sort routines:
 *  1. Sort a leaf cell in direction a (this resets the x_diff_sort).
 *  2. Drift all the particles by +D.
 *  3. Sort the cell in a new direction b without resetting.
 *  4. Drift all the particles by -D.
 * The particles are back where they were at the reset, so x_diff_sort = 0,
 * but they moved by D since the sort of direction b. D is kept below
 * space_maxreldx * dmin so that the cell is never considered stale.
 *
 * Test A checks the invariant |x . s_j - d_j| <= dx_max_sort for every sorted
 * direction j and every particle (hydro and stars).
 *
 * Test B shows the consequence in the density loop: a gas particle in ci
 * misses a neighbour in cj whose sort was added late, compared to the brute
 * force answer.
 *
 * Each test is also run as a control where all the directions are sorted at
 * the reset (no late addition), which must always pass.
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/* Just a forward declaration... */
void runner_dopair1_branch_density(struct runner *r, struct cell *ci,
                                   struct cell *cj, int limit_h_min,
                                   int limit_h_max);

/* Late-added direction: the x face, i.e. the direction of the pair in test B */
#define SID_LATE 4
/* Direction sorted at the reset */
#define SID_EARLY 12

/* Displacement applied (in units of dmin = 1) */
#define DISP 0.08

/* Number of particles per dimension in the cells of test A */
#define N_SIDE 6

/* Tolerance on the invariant (float sort entries) */
#define TOL 1e-5

/**
 * @brief Allocates an empty leaf cell of size 1 at a given location.
 */
static struct cell *make_empty_cell(const double loc[3]) {

  struct cell *c = NULL;
  if (posix_memalign((void **)&c, cell_align, sizeof(struct cell)) != 0)
    error("Couldn't allocate the cell");
  bzero(c, sizeof(struct cell));

  for (int k = 0; k < 3; k++) {
    c->loc[k] = loc[k];
    c->width[k] = 1.;
  }
  c->dmin = 1.;
  c->split = 0;
  c->depth = 0;
  c->nodeID = 0;
  c->h_min_allowed = c->dmin * 0.5 * (1. / kernel_gamma);
  c->h_max_allowed = c->dmin * (1. / kernel_gamma);
  c->hydro.super = c;
  c->hydro.ti_old_part = 8;
  c->hydro.ti_end_min = 8;
  c->stars.ti_old_part = 8;
  c->stars.ti_end_min = 8;
  return c;
}

/**
 * @brief Adds gas particles to an empty cell.
 */
static void add_parts(struct cell *c, const int count, const double (*x)[3],
                      const float *h, long long *id) {

  if (posix_memalign((void **)&c->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0)
    error("Couldn't allocate the parts");
  if (posix_memalign((void **)&c->hydro.xparts, xpart_align,
                     count * sizeof(struct xpart)) != 0)
    error("Couldn't allocate the xparts");
  bzero(c->hydro.parts, count * sizeof(struct part));
  bzero(c->hydro.xparts, count * sizeof(struct xpart));

  float h_max = 0.f;
  for (int i = 0; i < count; i++) {
    struct part *p = &c->hydro.parts[i];
    for (int k = 0; k < 3; k++) p->x[k] = x[i][k];
    p->h = h[i];
    hydro_set_mass(p, 1.f);

    /* Give the particles a valid thermodynamic state: the drift predicts
     * the pressure and the sound speed, which would otherwise be 0/0. The
     * density is set last as hydro_first_init_part() zeroes it. */
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
    p->conserved.energy = 1.f;
#else
    p->rho = 1.f;
    hydro_set_init_internal_energy(p, 1.f);
#endif
#if defined(PLANETARY_SPH) || defined(REMIX_SPH)
    set_idg_def(&eos.idg_def, 0);
    p->mat_id = 0;
#endif
    hydro_first_init_part(p, &c->hydro.xparts[i]);
#if !defined(GIZMO_MFV_SPH) && !defined(GIZMO_MFM_SPH)
    p->rho = 1.f;
#endif

    p->id = ++(*id);
    p->time_bin = 1;
    p->depth_h = 0;
#ifdef SWIFT_DEBUG_CHECKS
    p->ti_drift = 8;
    p->ti_kick = 8;
#endif
    h_max = max(h_max, p->h);
  }
  c->hydro.count = count;
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max;
}

/**
 * @brief Adds star particles to an empty cell.
 */
static void add_sparts(struct cell *c, const int count, const double (*x)[3],
                       long long *id) {

  if (posix_memalign((void **)&c->stars.parts, spart_align,
                     count * sizeof(struct spart)) != 0)
    error("Couldn't allocate the sparts");
  bzero(c->stars.parts, count * sizeof(struct spart));

  for (int i = 0; i < count; i++) {
    struct spart *sp = &c->stars.parts[i];
    for (int k = 0; k < 3; k++) sp->x[k] = x[i][k];
    sp->h = 0.01f;
    sp->id = ++(*id);
    sp->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
    sp->ti_drift = 8;
    sp->ti_kick = 8;
#endif
  }
  c->stars.count = count;
  c->stars.h_max = 0.01f;
  c->stars.h_max_active = 0.01f;
}

static void free_cell(struct cell *c) {
  free(c->hydro.parts);
  free(c->hydro.xparts);
  free(c->hydro.sort);
  free(c->stars.parts);
  free(c->stars.sort);
  free(c);
}

/**
 * @brief Sets the drift velocity of all the particles in a cell.
 */
static void set_velocity(struct cell *c, const double v[3]) {
  for (int i = 0; i < c->hydro.count; i++)
    for (int k = 0; k < 3; k++) {
      c->hydro.xparts[i].v_full[k] = v[k];
      c->hydro.parts[i].v[k] = v[k];
    }
  for (int i = 0; i < c->stars.count; i++)
    for (int k = 0; k < 3; k++) c->stars.parts[i].v[k] = v[k];
}

/**
 * @brief Moves the engine forward by one unit of time and drifts the cells.
 */
static void step(struct engine *e, struct cell **cells, const int n) {
  e->ti_current += 8;
  e->time += 1.;
  for (int i = 0; i < n; i++) {
    if (cells[i]->hydro.count)
      cell_drift_part(cells[i], e, /*force=*/1, /*init=*/0, NULL);
    if (cells[i]->stars.count)
      cell_drift_spart(cells[i], e, /*force=*/1, /*init=*/0, NULL);
  }
}

/**
 * @brief Largest violation of the sort invariant over all sorted directions.
 *
 * @return max_{j, k} (|x_k . s_j - d_jk| - dx_max_sort), positive if violated.
 */
static double hydro_violation(const struct cell *c) {
  double worst = -FLT_MAX;
  for (int j = 0; j < 13; j++) {
    if (!(c->hydro.sorted & (1 << j))) continue;
    const struct sort_entry *sort = cell_get_hydro_sorts(c, j);
    for (int k = 0; k < c->hydro.count; k++) {
      const struct part *p = &c->hydro.parts[sort[k].i];
      const double d = p->x[0] * runner_shift[j][0] +
                       p->x[1] * runner_shift[j][1] +
                       p->x[2] * runner_shift[j][2];
      worst = max(worst, fabs(d - sort[k].d) - c->hydro.dx_max_sort);
    }
  }
  return worst;
}

static double stars_violation(const struct cell *c) {
  double worst = -FLT_MAX;
  for (int j = 0; j < 13; j++) {
    if (!(c->stars.sorted & (1 << j))) continue;
    const struct sort_entry *sort = cell_get_stars_sorts(c, j);
    for (int k = 0; k < c->stars.count; k++) {
      const struct spart *sp = &c->stars.parts[sort[k].i];
      const double d = sp->x[0] * runner_shift[j][0] +
                       sp->x[1] * runner_shift[j][1] +
                       sp->x[2] * runner_shift[j][2];
      worst = max(worst, fabs(d - sort[k].d) - c->stars.dx_max_sort);
    }
  }
  return worst;
}

/**
 * @brief Test A: the invariant on a single cell of gas and star particles.
 *
 * @return 1 if the invariant holds, 0 otherwise.
 */
static int test_invariant(struct runner *r, struct engine *e, const int late) {

  static long long id = 0;
  const double loc[3] = {0., 0., 0.};
  struct cell *c = make_empty_cell(loc);

  /* A perturbed lattice of particles */
  const int count = N_SIDE * N_SIDE * N_SIDE;
  double (*x)[3] = malloc(count * sizeof(double[3]));
  float *h = malloc(count * sizeof(float));
  for (int i = 0; i < count; i++) {
    const int ix = i / (N_SIDE * N_SIDE);
    const int iy = (i / N_SIDE) % N_SIDE;
    const int iz = i % N_SIDE;
    x[i][0] = (ix + 0.5 + random_uniform(-0.3, 0.3)) / N_SIDE;
    x[i][1] = (iy + 0.5 + random_uniform(-0.3, 0.3)) / N_SIDE;
    x[i][2] = (iz + 0.5 + random_uniform(-0.3, 0.3)) / N_SIDE;
    h[i] = 0.01f;
  }
  add_parts(c, count, (const double (*)[3])x, h, &id);
  add_sparts(c, count, (const double (*)[3])x, &id);
  free(x);
  free(h);
  e->ti_current = 8;

  /* 1. Sort at the reset */
  const int flags_reset =
      late ? (1 << SID_EARLY) : ((1 << SID_EARLY) | (1 << SID_LATE));
  runner_do_hydro_sort(r, c, flags_reset, /*cleanup=*/1, /*lock=*/0,
                       /*rt_requests_sort=*/0, /*clock=*/0);
  runner_do_stars_sort(r, c, flags_reset, /*cleanup=*/1, /*clock=*/0);

  /* 2. Move by +D along the late direction */
  const double v_plus[3] = {DISP * runner_shift[SID_LATE][0],
                            DISP * runner_shift[SID_LATE][1],
                            DISP * runner_shift[SID_LATE][2]};
  set_velocity(c, v_plus);
  step(e, &c, 1);

  /* 3. Add the late direction without resetting */
  if (late) {
    runner_do_hydro_sort(r, c, 1 << SID_LATE, /*cleanup=*/0, /*lock=*/0,
                         /*rt_requests_sort=*/0, /*clock=*/0);
    runner_do_stars_sort(r, c, 1 << SID_LATE, /*cleanup=*/0, /*clock=*/0);
  }

  /* 4. Move back by -D */
  const double v_minus[3] = {-v_plus[0], -v_plus[1], -v_plus[2]};
  set_velocity(c, v_minus);
  step(e, &c, 1);

  /* The cell must not be stale, otherwise the normal machinery re-sorts it */
  if (c->hydro.dx_max_sort > space_maxreldx * c->dmin ||
      c->stars.dx_max_sort > space_maxreldx * c->dmin)
    error("Test set-up error: the cell went stale.");

  const double viol_h = hydro_violation(c);
  const double viol_s = stars_violation(c);

  message("[%s] hydro: dx_max_sort=%.4e worst |d - d_sort| - dx_max_sort=%+.4e",
          late ? "late " : "ctrl ", c->hydro.dx_max_sort, viol_h);
  message("[%s] stars: dx_max_sort=%.4e worst |d - d_sort| - dx_max_sort=%+.4e",
          late ? "late " : "ctrl ", c->stars.dx_max_sort, viol_s);

  free_cell(c);
  return (viol_h <= TOL) && (viol_s <= TOL);
}

/**
 * @brief Test B: a density pair where cj had its pair direction added late.
 *
 * ci holds one gas particle at x = 0.9 with h * kernel_gamma = 0.2. cj holds
 * a particle at x = 1.05 (a true neighbour at distance 0.15) and some
 * particles far from the interface. cj is moved away from ci by D, sorted
 * along x, and moved back. Its sort array then places the neighbour at
 * x = 1.13, beyond the search radius unless dx_max_sort accounts for D.
 *
 * @return 1 if the pair finds the same number of neighbours as brute force.
 */
static int test_pair(struct runner *r, struct engine *e, const int late) {

  static long long id = 1000000;
  const double loc_i[3] = {0., 0., 0.};
  const double loc_j[3] = {1., 0., 0.};
  struct cell *ci = make_empty_cell(loc_i);
  struct cell *cj = make_empty_cell(loc_j);
  struct cell *cells[2] = {ci, cj};

  const double xi[1][3] = {{0.9, 0.5, 0.5}};
  const float hi[1] = {0.2f / kernel_gamma};
  add_parts(ci, 1, xi, hi, &id);

  const double xj[4][3] = {
      {1.05, 0.5, 0.5}, {1.6, 0.2, 0.3}, {1.7, 0.8, 0.6}, {1.8, 0.4, 0.9}};
  const float hj[4] = {0.01f, 0.01f, 0.01f, 0.01f};
  add_parts(cj, 4, xj, hj, &id);
  e->ti_current = 8;

  /* 1. Sort at the reset */
  runner_do_hydro_sort(r, ci, 0x1FFF, /*cleanup=*/1, 0, 0, 0);
  const int flags_reset =
      late ? (1 << SID_EARLY) : ((1 << SID_EARLY) | (1 << SID_LATE));
  runner_do_hydro_sort(r, cj, flags_reset, /*cleanup=*/1, 0, 0, 0);

  /* 2. Move cj away from ci by D */
  const double v_zero[3] = {0., 0., 0.};
  const double v_plus[3] = {DISP, 0., 0.};
  set_velocity(ci, v_zero);
  set_velocity(cj, v_plus);
  step(e, cells, 2);

  /* 3. Add the pair direction to cj without resetting */
  if (late) runner_do_hydro_sort(r, cj, 1 << SID_LATE, /*cleanup=*/0, 0, 0, 0);

  /* 4. Move cj back */
  const double v_minus[3] = {-DISP, 0., 0.};
  set_velocity(cj, v_minus);
  step(e, cells, 2);

  /* Prepare the cells as the unskip would */
  for (int i = 0; i < 2; i++) {
    cells[i]->hydro.dx_max_sort_old = cells[i]->hydro.dx_max_sort;
    cells[i]->hydro.ti_end_min = e->ti_current;
    for (int k = 0; k < cells[i]->hydro.count; k++)
      hydro_init_part(&cells[i]->hydro.parts[k], NULL);
  }

  /* Brute-force neighbour count of the particle in ci */
  const struct part *pi = &ci->hydro.parts[0];
  int n_brute = 0;
  for (int k = 0; k < cj->hydro.count; k++) {
    const struct part *pj = &cj->hydro.parts[k];
    double r2 = 0.;
    for (int d = 0; d < 3; d++)
      r2 += (pi->x[d] - pj->x[d]) * (pi->x[d] - pj->x[d]);
    if (r2 < pi->h * pi->h * kernel_gamma2) n_brute++;
  }

  /* The sorted pair interaction */
  runner_dopair1_branch_density(r, ci, cj, /*limit_min_h=*/0,
                                /*limit_max_h=*/0);

  /* Neighbour count of the SPH loop, from the un-normalised weight */
  const int n_pair = pi->density.wcount > 0.f ? 1 : 0;

  message(
      "[%s] pair: dx_max_sort(cj)=%.4e brute-force neighbours=%d "
      "found by the pair=%d (wcount=%e)",
      late ? "late " : "ctrl ", cj->hydro.dx_max_sort, n_brute, n_pair,
      pi->density.wcount);

  free_cell(ci);
  free_cell(cj);
  return n_pair == n_brute;
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  srand(1234);

  /* Build the infrastructure */
  struct space space;
  bzero(&space, sizeof(struct space));
  space.periodic = 1;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  hp.h_max = FLT_MAX;
  hp.h_min = 0.f;

  struct phys_const prog_const;
  bzero(&prog_const, sizeof(struct phys_const));
  prog_const.const_vacuum_permeability = 1.0;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  struct entropy_floor_properties entropy_floor;
  bzero(&entropy_floor, sizeof(struct entropy_floor_properties));

  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));

  struct lightcone_array_props lightcone_array_properties;
  bzero(&lightcone_array_properties, sizeof(struct lightcone_array_props));

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.time = 0.;
  engine.time_base = 1. / 8.;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;
  engine.physical_constants = &prog_const;
  engine.cosmology = &cosmo;
  engine.entropy_floor = &entropy_floor;
  engine.pressure_floor_props = &pressure_floor;
  engine.lightcone_array_properties = &lightcone_array_properties;
  engine.nodeID = 0;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  message("Displacement D = %.3f, staleness threshold = %.3f", DISP,
          space_maxreldx);

  const int ctrl_inv = test_invariant(&runner, &engine, /*late=*/0);
  const int late_inv = test_invariant(&runner, &engine, /*late=*/1);
  const int ctrl_pair = test_pair(&runner, &engine, /*late=*/0);
  const int late_pair = test_pair(&runner, &engine, /*late=*/1);

  message("Test A (invariant)  control: %s  late sort: %s",
          ctrl_inv ? "PASS" : "FAIL", late_inv ? "PASS" : "FAIL");
  message("Test B (pair)       control: %s  late sort: %s",
          ctrl_pair ? "PASS" : "FAIL", late_pair ? "PASS" : "FAIL");

  return (ctrl_inv && late_inv && ctrl_pair && late_pair) ? 0 : 1;
}
