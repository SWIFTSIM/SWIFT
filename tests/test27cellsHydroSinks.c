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

/* Test for the gas-gas neighbour loop used in sink particle formation.
 *
 * Generates 27 cells arranged in a 3x3x3 grid, fills them with gas particles,
 * and verifies that the optimised pair and self interaction functions
 *   runner_dopair1_branch_hydro_aperture_test_formation()
 *   runner_doself1_branch_hydro_aperture_test_formation()
 * produce exactly the same neighbour count (density.wcount) as the reference
 * brute-force implementation that loops over all pairs within the fixed
 * aperture radius r_cut.
 *
 * The loop template is instantiated locally in this translation unit with a
 * model-agnostic test iact that counts each pair found within r_cut by
 * incrementing density.wcount.  This makes the test independent of any
 * particular sink model and sensitive to any error in the geometric
 * neighbour-finding logic (wrong pairs included or excluded).
 *
 * Gas particles may have any smoothing length — the loop's sole cutoff
 * criterion is the fixed geometric aperture r_cut, so no h-based constraint
 * is needed.
 *
 * Usage: test27cellsHydroSinks -n N -r RUNS [-s SIZE] [-d PERT] [-f SUFFIX]
 *   -n N      Cube root of particles per cell (e.g. 3 → 27 particles/cell)
 *   -r RUNS   Number of optimised-version runs (brute force runs once)
 *   -s SIZE   Physical size of each cell (default 1.0)
 *   -d PERT   Position perturbation in units of inter-particle spacing [0,1)
 *   -f SUFFIX Output file name suffix
 */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

/* The fixed aperture radius as a fraction of the cell size. */
#define R_CUT_FRACTION 0.25f

/* Smoothing length in units of inter-particle spacing.  The loop uses a fixed
 * geometric r_cut, so h may be anything — use a value similar to standard
 * hydro tests. */
#define H_FRAC 1.2348f

#define NODE_ID 0

/* ============================================================
 * Test-local iact: count neighbours within r_cut.
 *
 * These functions are model-agnostic — they accumulate density.wcount so
 * that both the optimised loop and the brute-force reference produce an
 * integer neighbour count per particle, making disagreements detectable.
 * ============================================================ */

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_hydro_aperture_test_formation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H, const int with_self_gravity, const struct cosmology *cosmo,
    const struct sink_props *sink_props) {
  pi->density.wcount += 1.0f;
}

__attribute__((always_inline)) INLINE static void
runner_iact_hydro_aperture_test_formation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const int with_self_gravity, const struct cosmology *cosmo,
    const struct sink_props *sink_props) {
  pi->density.wcount += 1.0f;
  pj->density.wcount += 1.0f;
}

/* ============================================================
 * Instantiate the loop template with the test iact.
 *
 * IACT_NONSYM_HYDRO_APERTURE / IACT_HYDRO_APERTURE are defined in
 * runner_doiact_hydro_aperture.h (included by
 * runner_doiact_functions_hydro_aperture.h) and expand to
 * runner_iact_nonsym_hydro_aperture_FUNCTION / runner_iact_hydro_aperture_FUNCTION.
 * With FUNCTION=test_formation they resolve to the test functions above, so the
 * generated loops count neighbours.
 *
 * space_getsid_and_swap_cells is used internally by the pair functions; pull
 * in its declaration here just as runner_doiact_hydro_aperture.c does.
 * ============================================================ */
#include "space_getsid.h"
#define FUNCTION test_formation
#define FUNCTION_TASK_LOOP TASK_LOOP_PREP_SINK_FORMATION
#include "runner_doiact_functions_hydro_aperture.h"
#include "runner_doiact_undef.h"

/* The locally instantiated function names (used directly below). */
#define DOSELF1_NAME "runner_doself1_branch_hydro_aperture_test_formation"
#define DOPAIR1_NAME "runner_dopair1_branch_hydro_aperture_test_formation"

/* ============================================================
 * Cell construction helper.
 * ============================================================ */

/**
 * @brief Construct a cell and all of its particles in a valid state.
 *
 * Particles are placed on a regular cubic grid with optional position
 * perturbations.
 *
 * @param n Cube root of the number of particles.
 * @param offset Position of the cell corner relative to (0,0,0).
 * @param size Physical side length of the cell.
 * @param h_frac Smoothing length in units of the inter-particle spacing.
 * @param partId Running counter for unique particle IDs.
 * @param pert Position perturbation amplitude in units of spacing.
 */
struct cell *make_cell(size_t n, double *offset, double size, double h_frac,
                       long long *partId, double pert) {

  const size_t count = n * n * n;
  float h_max = 0.f;

  struct cell *cell = NULL;
  if (posix_memalign((void **)&cell, cell_align, sizeof(struct cell)) != 0)
    error("Couldn't allocate cell.");
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0)
    error("Couldn't allocate particles (%d).", (int)count);
  bzero(cell->hydro.parts, count * sizeof(struct part));

  /* Build the particle grid. */
  struct part *part = cell->hydro.parts;
  for (size_t x = 0; x < n; ++x) {
    for (size_t y = 0; y < n; ++y) {
      for (size_t z = 0; z < n; ++z) {
        part->x[0] =
            offset[0] +
            size * (x + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;
        part->x[1] =
            offset[1] +
            size * (y + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;
        part->x[2] =
            offset[2] +
            size * (z + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;

        part->v[0] = 0.f;
        part->v[1] = 0.f;
        part->v[2] = 0.f;

        /* Smoothing length: h_frac * spacing, where spacing = size/n. */
        part->h = (float)(size * h_frac / (float)n);
        h_max = fmaxf(h_max, part->h);

        part->id = ++(*partId);
        part->time_bin = 1;

        /* Set mass so that the density loop has a non-trivial input. */
        hydro_set_mass(part, 1.0f);

#ifdef SWIFT_DEBUG_CHECKS
        part->ti_drift = 8;
        part->ti_kick = 8;
#endif
        ++part;
      }
    }
  }

  /* Cell metadata. */
  cell->split = 0;
  cell->hydro.h_max = h_max;
  cell->hydro.h_max_active = h_max;
  cell->hydro.count = count;
  cell->hydro.dx_max_part = 0.f;
  cell->hydro.dx_max_sort = 0.f;
  cell->hydro.dx_max_sort_old = 0.f;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->dmin = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];
  cell->h_min_allowed = size * 0.5f * (1.f / kernel_gamma);
  cell->h_max_allowed = size * (1.f / kernel_gamma);
  cell->hydro.super = cell;
  cell->hydro.ti_old_part = 8;
  cell->hydro.ti_end_min = 8;
  cell->grav.ti_old_part = 8;
  cell->grav.ti_end_min = 8;
  cell->nodeID = NODE_ID;

  shuffle_particles(cell->hydro.parts, cell->hydro.count);

  cell->hydro.sorted = 0;
  cell->hydro.sort = NULL;

  return cell;
}

/**
 * @brief Free a cell and its particle arrays.
 */
void clean_up(struct cell *ci) {
  free(ci->hydro.parts);
  free(ci->hydro.sort);
  free(ci);
}

/* ============================================================
 * Field reset and dump helpers.
 * ============================================================ */

/**
 * @brief Reset density.wcount to zero for all particles in a cell.
 */
void zero_particle_fields(struct cell *c) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    c->hydro.parts[pid].density.wcount = 0.f;
  }
}

/**
 * @brief Dump density.wcount and particle data for the main cell and all 26
 *        neighbouring cells to @p fileName.
 *
 * The wcount column holds the neighbour count accumulated by the test iact.
 * The output is compared between the optimised and brute-force runs to verify
 * that both find exactly the same set of pairs within r_cut.
 */
void dump_particle_fields(const char *fileName, struct cell *main_cell,
                          struct cell **cells) {
  FILE *file = fopen(fileName, "w");
  if (!file) error("Could not open output file '%s'.", fileName);

  /* Header. */
  fprintf(file, "# %6s %10s %10s %10s %10s\n", "ID", "pos_x", "pos_y", "pos_z",
          "wcount");

  fprintf(file,
          "# Main cell ------------------------------------------------\n");

  /* Main cell particles. */
  for (int pid = 0; pid < main_cell->hydro.count; pid++) {
    const struct part *p = &main_cell->hydro.parts[pid];
    fprintf(file, "%8llu %10.6f %10.6f %10.6f %10.6f\n", p->id, p->x[0],
            p->x[1], p->x[2], p->density.wcount);
  }

  /* Neighbouring cells. */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        struct cell *cj = cells[i * 9 + j * 3 + k];
        if (cj == main_cell) continue;
        fprintf(
            file,
            "# Offset: [%2d %2d %2d] ------------------------------------\n",
            i - 1, j - 1, k - 1);
        for (int pjd = 0; pjd < cj->hydro.count; pjd++) {
          const struct part *p = &cj->hydro.parts[pjd];
          fprintf(file, "%8llu %10.6f %10.6f %10.6f %10.6f\n", p->id, p->x[0],
                  p->x[1], p->x[2], p->density.wcount);
        }
      }
    }
  }
  fclose(file);
}

/* ============================================================
 * Brute-force reference implementations.
 *
 * These O(N²) functions are guaranteed correct by construction and serve as
 * the ground truth against which the optimised loops are verified.  They call
 * the same test iact as the loop template, so the comparison is exact.
 * ============================================================ */

/**
 * @brief Brute-force pair interaction between all gas particles in @p ci and
 *        @p cj within the fixed aperture radius @p r_cut.
 *
 * Applies the non-symmetric interaction to (active pi, any pj) pairs in the
 * ci→cj direction and to (active pj, any pi) pairs in the cj→ci direction.
 * Periodic shift between cell centres is computed internally.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius.
 */
void pairs_all_hydro_sinks(struct runner *r, struct cell *ci, struct cell *cj,
                           const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;
  const float r_cut2 = r_cut * r_cut;

  /* Periodic shift between cell centres. */
  double shift[3] = {0., 0., 0.};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Forward direction: active pi from ci interacts with all pj in cj. */
  for (int i = 0; i < ci->hydro.count; i++) {
    struct part *pi = &ci->hydro.parts[i];
    if (!part_is_active(pi, e)) continue;
    if (part_is_inhibited(pi, e)) continue;

    for (int j = 0; j < cj->hydro.count; j++) {
      struct part *pj = &cj->hydro.parts[j];
      if (part_is_inhibited(pj, e)) continue;

      const float dx[3] = {(float)(pi->x[0] - pj->x[0] - shift[0]),
                           (float)(pi->x[1] - pj->x[1] - shift[1]),
                           (float)(pi->x[2] - pj->x[2] - shift[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (r2 < r_cut2) {
        runner_iact_nonsym_hydro_aperture_test_formation(r2, dx, pi->h, pj->h, pi,
                                                      pj, a, H,
                                                      /*with_self_gravity=*/0,
                                                      /*cosmo=*/NULL,
                                                      /*sink_props=*/NULL);
      }
    }
  }

  /* Reverse direction: active pj from cj interacts with all pi in ci. */
  for (int j = 0; j < cj->hydro.count; j++) {
    struct part *pj = &cj->hydro.parts[j];
    if (!part_is_active(pj, e)) continue;
    if (part_is_inhibited(pj, e)) continue;

    for (int i = 0; i < ci->hydro.count; i++) {
      struct part *pi = &ci->hydro.parts[i];
      if (part_is_inhibited(pi, e)) continue;

      /* dx points from pi to pj (pj is the receiver). */
      const float mdx[3] = {(float)(pj->x[0] - pi->x[0] + shift[0]),
                            (float)(pj->x[1] - pi->x[1] + shift[1]),
                            (float)(pj->x[2] - pi->x[2] + shift[2])};
      const float r2 = mdx[0] * mdx[0] + mdx[1] * mdx[1] + mdx[2] * mdx[2];

      if (r2 < r_cut2) {
        runner_iact_nonsym_hydro_aperture_test_formation(r2, mdx, pj->h, pi->h, pj,
                                                      pi, a, H,
                                                      /*with_self_gravity=*/0,
                                                      /*cosmo=*/NULL,
                                                      /*sink_props=*/NULL);
      }
    }
  }
}

/**
 * @brief Brute-force self interaction for all gas particles in @p c within the
 *        fixed aperture radius @p r_cut.
 *
 * Iterates over all unordered pairs (i < j) and applies the non-symmetric
 * interaction to each active particle independently.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius.
 */
void self_all_hydro_sinks(struct runner *r, struct cell *c, const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;
  const float r_cut2 = r_cut * r_cut;

  for (int i = 0; i < c->hydro.count; i++) {
    struct part *pi = &c->hydro.parts[i];
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_active(pi, e);

    /* Only visit j > i to cover each pair once. */
    for (int j = i + 1; j < c->hydro.count; j++) {
      struct part *pj = &c->hydro.parts[j];
      if (part_is_inhibited(pj, e)) continue;

      const int pj_active = part_is_active(pj, e);

      const float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                           (float)(pi->x[1] - pj->x[1]),
                           (float)(pi->x[2] - pj->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (r2 < r_cut2) {
        /* Update active pi with contribution from pj. */
        if (pi_active) {
          runner_iact_nonsym_hydro_aperture_test_formation(
              r2, dx, pi->h, pj->h, pi, pj, a, H, /*with_self_gravity=*/0,
              /*cosmo=*/NULL, /*sink_props=*/NULL);
        }
        /* Update active pj with contribution from pi. */
        if (pj_active) {
          const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
          runner_iact_nonsym_hydro_aperture_test_formation(
              r2, mdx, pj->h, pi->h, pj, pi, a, H, /*with_self_gravity=*/0,
              /*cosmo=*/NULL, /*sink_props=*/NULL);
        }
      }
    }
  }
}

/* ============================================================
 * Main test driver.
 * ============================================================ */

int main(int argc, char *argv[]) {

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  size_t runs = 0, particles = 0;
  double size = 1.;
  double perturbation = 0.;
  char outputFileNameExtension[100] = "";
  char outputFileName[200] = "";

  /* Initialise CPU frequency counter. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Enable FP exception trapping. */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Fixed random seed for reproducibility. */
  srand(42);

  int c;
  while ((c = getopt(argc, argv, "s:n:r:d:f:")) != -1) {
    switch (c) {
      case 's':
        sscanf(optarg, "%lf", &size);
        break;
      case 'n':
        sscanf(optarg, "%zu", &particles);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case 'd':
        sscanf(optarg, "%lf", &perturbation);
        break;
      case 'f':
        strcpy(outputFileNameExtension, optarg);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (particles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -n PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 27 cells of gas particles and tests the hydro-sinks "
        "neighbour loop\n(%s and\n%s).\n"
        "\nOptions:\n"
        "  -n N     Cube root of particles per cell\n"
        "  -r RUNS  Number of runs for the optimised version\n"
        "  -s SIZE  Cell size (default 1.0)\n"
        "  -d PERT  Position perturbation in units of spacing [0,1)\n"
        "  -f SUFF  Output file name suffix\n",
        argv[0], DOPAIR1_NAME, DOSELF1_NAME);
    exit(1);
  }

  /* The fixed aperture radius and the smoothing length. */
  const float r_cut = (float)(R_CUT_FRACTION * size);
  const double h_frac = H_FRAC;

  message("DOSELF1 function: %s", DOSELF1_NAME);
  message("DOPAIR1 function: %s", DOPAIR1_NAME);
  message("Cell size:        %.3f", size);
  message("Aperture r_cut:   %.3f  (%.3f * cell_size)", r_cut, R_CUT_FRACTION);
  message("Smoothing length: %.3f  (h_frac = %.4f)", h_frac * size, h_frac);
  message("Kernel:           %s", kernel_name);
  printf("\n");

  /* ---- Infrastructure setup ---- */

  struct space space;
  space.periodic = 1;
  space.dim[0] = 3. * size;
  space.dim[1] = 3. * size;
  space.dim[2] = 3. * size;

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  hp.eta_neighbours = h_frac;
  hp.h_tolerance = 1e0;
  hp.h_max = FLT_MAX;
  hp.max_smoothing_iterations = 1;
  hp.CFL_condition = 0.1f;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
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
  sink_props.cut_off_radius = r_cut;
  engine.sink_properties = &sink_props;

  struct pressure_floor_props pressure_floor;
  bzero(&pressure_floor, sizeof(struct pressure_floor_props));
  engine.pressure_floor_props = &pressure_floor;

  struct lightcone_array_props lightcone_array_properties;
  lightcone_array_properties.nr_lightcones = 0;
  engine.lightcone_array_properties = &lightcone_array_properties;

  struct runner runner;
  runner.e = &engine;

  /* ---- Construct the 3x3x3 cell grid ---- */

  struct cell *cells[27];
  struct cell *main_cell;
  static long long partId = 0;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        double offset[3] = {i * size, j * size, k * size};
        cells[i * 9 + j * 3 + k] =
            make_cell(particles, offset, size, h_frac, &partId, perturbation);

        runner_do_drift_part(&runner, cells[i * 9 + j * 3 + k], 0);
        runner_do_hydro_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0, 0, 0,
                             0);
      }
    }
  }

  /* Centre cell is our main cell. */
  main_cell = cells[13];

  /* ---- Run the optimised interaction ---- */

  ticks timings[27];
  for (int i = 0; i < 27; i++) timings[i] = 0;
  ticks time = 0;

  for (size_t i = 0; i < runs; ++i) {

    /* Reset wcount = 0 for all particles. */
    for (int j = 0; j < 27; ++j) zero_particle_fields(cells[j]);

    const ticks tic = getticks();

    /* Pair interactions: main_cell against each of its 26 neighbours. */
    for (int j = 0; j < 27; ++j) {
      if (cells[j] != main_cell) {
        const ticks sub_tic = getticks();
        runner_dopair1_branch_hydro_aperture_test_formation(&runner, main_cell,
                                                         cells[j], r_cut);
        timings[j] += getticks() - sub_tic;
      }
    }

    /* Self interaction. */
    const ticks self_tic = getticks();
    runner_doself1_branch_hydro_aperture_test_formation(&runner, main_cell, r_cut);
    timings[13] += getticks() - self_tic;

    time += getticks() - tic;

    /* Dump the last run. */
    if (i == runs - 1) {
      sprintf(outputFileName, "swift_hydro_sinks_dopair_27_%.150s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell, cells);
    }
  }

  /* Timing report. */
  ticks corner_time = timings[0] + timings[2] + timings[6] + timings[8] +
                      timings[18] + timings[20] + timings[24] + timings[26];
  ticks edge_time = timings[1] + timings[3] + timings[5] + timings[7] +
                    timings[9] + timings[11] + timings[15] + timings[17] +
                    timings[19] + timings[21] + timings[23] + timings[25];
  ticks face_time = timings[4] + timings[10] + timings[12] + timings[14] +
                    timings[16] + timings[22];
  ticks self_time = timings[13];

  message("Corner calculations took:  %.3f %s.",
          clocks_from_ticks(corner_time / runs), clocks_getunit());
  message("Edge calculations took:    %.3f %s.",
          clocks_from_ticks(edge_time / runs), clocks_getunit());
  message("Face calculations took:    %.3f %s.",
          clocks_from_ticks(face_time / runs), clocks_getunit());
  message("Self calculations took:    %.3f %s.",
          clocks_from_ticks(self_time / runs), clocks_getunit());
  message("SWIFT total took:          %.3f %s.", clocks_from_ticks(time / runs),
          clocks_getunit());

  /* ---- Brute-force reference run ---- */

  /* Reset fields. */
  for (int i = 0; i < 27; ++i) zero_particle_fields(cells[i]);

  const ticks tic = getticks();

  /* Brute-force pair interactions. */
  for (int j = 0; j < 27; ++j)
    if (cells[j] != main_cell)
      pairs_all_hydro_sinks(&runner, main_cell, cells[j], r_cut);

  /* Brute-force self interaction. */
  self_all_hydro_sinks(&runner, main_cell, r_cut);

  const ticks toc = getticks();

  /* Dump brute-force results. */
  sprintf(outputFileName, "brute_force_hydro_sinks_27_%.150s.dat",
          outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell, cells);

  message("Brute-force calculation took: %.3f %s.",
          clocks_from_ticks(toc - tic), clocks_getunit());

  /* ---- Clean up ---- */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return 0;
}
