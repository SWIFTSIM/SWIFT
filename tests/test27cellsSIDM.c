/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2025 Katy Proctor (katy.proctor@fysik.su.se).
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
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

#if defined(WITH_VECTORIZATION)
#define DOSELF1 runner_doself1_branch_sidm_density
#define DOSELF1_SUBSET runner_doself_subset_branch_sidm_density
#define DOPAIR1_SUBSET runner_dopair_subset_branch_sidm_density
#define DOPAIR1 runner_dopair1_branch_sidm_density
#ifdef TEST_DOSELF_SUBSET
#define DOSELF1_NAME "runner_doself_subset_branch_sidm_density"
#else
#define DOSELF1_NAME "runner_doself1_branch_sidm_density"
#endif
#ifdef TEST_DOPAIR_SUBSET
#define DOPAIR1_NAME "runner_dopair_subset_branch_sidm_density"
#else
#define DOPAIR1_NAME "runner_dopair1_branch_sidm_density"
#endif
#endif

#ifndef DOSELF1
#define DOSELF1 runner_doself1_branch_sidm_density
#define DOSELF1_SUBSET runner_doself_subset_branch_sidm_density
#ifdef TEST_DOSELF_SUBSET
#define DOSELF1_NAME "runner_doself_subset_branch_sidm_density"
#else
#define DOSELF1_NAME "runner_doself1_branch_sidm_density"
#endif
#endif

#ifndef DOPAIR1
#define DOPAIR1 runner_dopair1_branch_sidm_density
#define DOPAIR1_SUBSET runner_dopair_subset_branch_sidm_density
#ifdef TEST_DOPAIR_SUBSET
#define DOPAIR1_NAME "runner_dopair1_subset_branch_sidm_density"
#else
#define DOPAIR1_NAME "runner_dopair1_branch_sidm_density"
#endif
#endif

#define NODE_ID 0

/**
 * @brief Constructs a cell and all of its particle in a valid state prior to
 * a DOPAIR or DOSELF calcuation.
 *
 * @param n_sidm The cube root of the number of siparticles.
 * @param offset The position of the cell offset from (0,0,0).
 * @param size The cell size.
 * @param h The smoothing length of the particles in units of the inter-particle
 * separation.
 * @param density The density of the fluid.
 * @param sipartId The running counter of IDs.
 * @param pert The perturbation to apply to the particles in the cell in units
 * of the inter-particle separation.
 * @param h_pert The perturbation to apply to the smoothing length.
 */
struct cell *make_cell(size_t n_sidm, double *offset, double size, double h,
                       double density, long long *sipartId, double pert,
                       double h_pert) {
  const size_t sicount = n_sidm * n_sidm * n_sidm;
  const double volume = size * size * size;
  float h_max = 0.f;
  struct cell *cell = NULL;
  if (posix_memalign((void **)&cell, cell_align, sizeof(struct cell)) != 0) {
    error("Couldn't allocate the cell");
  }
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->sidm.parts, sipart_align,
                     sicount * sizeof(struct sipart)) != 0) {
    error("couldn't allocate siparticles, no. of siparticles: %d",
          (int)sicount);
  }
  bzero(cell->sidm.parts, sicount * sizeof(struct sipart));

  /* Construct the parts */
  struct sipart *sipart = cell->sidm.parts;
  for (size_t x = 0; x < n_sidm; ++x) {
    for (size_t y = 0; y < n_sidm; ++y) {
      for (size_t z = 0; z < n_sidm; ++z) {
        sipart->x[0] =
            offset[0] +
            size * (x + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n_sidm;
        sipart->x[1] =
            offset[1] +
            size * (y + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n_sidm;
        sipart->x[2] =
            offset[2] +
            size * (z + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n_sidm;

        sipart->v[0] = 0.f;
        sipart->v[1] = 0.f;
        sipart->v[2] = 0.f;

        if (h_pert)
          sipart->h = size * h * random_uniform(1.f, h_pert) / (float)n_sidm;
        else
          sipart->h = size * h / (float)n_sidm;
        h_max = fmaxf(h_max, sipart->h);
        sipart->id = ++(*sipartId);
        // sipart->depth_h = 0;
        sipart->mass = density * volume / sicount;
        sipart->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        sipart->ti_drift = 8;
        sipart->ti_kick = 8;
#endif

        ++sipart;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  // cell->depth = 0;
  cell->sidm.h_max = h_max;
  cell->sidm.h_max_active = h_max;
  cell->sidm.count = sicount;
  cell->sidm.dx_max_part = 0.;
  cell->sidm.dx_max_sort = 0.;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->dmin = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];
  cell->h_min_allowed = cell->dmin * 0.5 * (1. / kernel_gamma);
  cell->h_max_allowed = cell->dmin * (1. / kernel_gamma);

  cell->sidm.super = cell;
  cell->sidm.ti_old_part = 8;
  cell->sidm.ti_end_min = 8;
  cell->nodeID = NODE_ID;

  shuffle_siparticles(cell->sidm.parts, cell->sidm.count);

  cell->sidm.sorted = 0;
  cell->sidm.sort = NULL;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->sidm.parts);
  free(ci->sidm.sort);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields(struct cell *c) {

  for (int siid = 0; siid < c->sidm.count; siid++) {
    sidm_init_sipart(&c->sidm.parts[siid]);
  }
}

/**
 * @brief Ends the loop by adding the appropriate coefficients
 */
void end_calculation(struct cell *c) {

  for (int siid = 0; siid < c->sidm.count; siid++) {
    sidm_end_density(&c->sidm.parts[siid]);

    /* Recover the common "Neighbour number" definition */
    c->sidm.parts[siid].density.wcount *= pow_dimension(c->sidm.parts[siid].h);
    c->sidm.parts[siid].density.wcount *= kernel_norm;
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *main_cell,
                          struct cell **cells) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file, "# %4s %10s %10s %10s %13s %13s %13s %13s\n", "ID", "pos_x",
          "pos_y", "pos_z", "rho", "rho_dh", "wcount", "wcount_dh");

  fprintf(file, "# Main cell --------------------------------------------\n");

  /* Write main cell */
  for (int siid = 0; siid < main_cell->sidm.count; siid++) {
    fprintf(file, "%6llu %10f %10f %10f %13e %13e %13e %13e\n",
            main_cell->sidm.parts[siid].id, main_cell->sidm.parts[siid].x[0],
            main_cell->sidm.parts[siid].x[1], main_cell->sidm.parts[siid].x[2],
            sidm_get_comoving_density(&main_cell->sidm.parts[siid]),
            main_cell->sidm.parts[siid].density.rho_dh,
            main_cell->sidm.parts[siid].density.wcount,
            main_cell->sidm.parts[siid].density.wcount_dh);
  }

  /* Write all other cells */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        struct cell *cj = cells[i * 9 + j * 3 + k];
        if (cj == main_cell) continue;

        fprintf(file,
                "# Offset: [%2d %2d %2d] -----------------------------------\n",
                i - 1, j - 1, k - 1);

        for (int sipjd = 0; sipjd < cj->sidm.count; sipjd++) {
          fprintf(file, "%6llu %10f %10f %10f %13e %13e %13e %13e\n",
                  cj->sidm.parts[sipjd].id, cj->sidm.parts[sipjd].x[0],
                  cj->sidm.parts[sipjd].x[1], cj->sidm.parts[sipjd].x[2],
                  sidm_get_comoving_density(&main_cell->sidm.parts[sipjd]),
                  main_cell->sidm.parts[sipjd].density.rho_dh,
                  cj->sidm.parts[sipjd].density.wcount,
                  cj->sidm.parts[sipjd].density.wcount_dh);
        }
      }
    }
  }
  fclose(file);
}

/* Just a forward declaration... */
void runner_dopair1_branch_sidm_density(struct runner *r, struct cell *ci,
                                        struct cell *cj, int limit_h_min,
                                        int limit_h_max);
void runner_doself1_branch_sidm_density(struct runner *r, struct cell *c,
                                        int limit_h_min, int limit_h_max);
void runner_dopair_subset_branch_sidm_density(struct runner *r,
                                              struct cell *restrict ci,
                                              struct sipart *restrict siparts_i,
                                              int *restrict ind, int sicount,
                                              struct cell *restrict cj);
void runner_doself_subset_branch_sidm_density(struct runner *r,
                                              struct cell *restrict ci,
                                              struct sipart *restrict siparts,
                                              int *restrict ind, int sicount);

/* And go... */
int main(int argc, char *argv[]) {

  /* Do not run test if code is not compiled in SIDM mode  */
#ifdef SIDM_NONE
  FILE *f = fopen("SIDM_NONE.dat", "w");
  fclose(f);
  return 0;
#endif

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  size_t runs = 0, siparticles = 0;
  double h = 1.23485, size = 1., rho = 1.;
  double perturbation = 0., h_pert = 0.;
  char outputFileNameExtension[100] = "";
  char outputFileName[200] = "";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  srand(0);

  int c;
  while ((c = getopt(argc, argv, "m:s:h:p:n:r:t:d:f:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 'p':
        sscanf(optarg, "%lf", &h_pert);
        break;
      case 's':
        sscanf(optarg, "%lf", &size);
        break;
      case 'n':
        sscanf(optarg, "%zu", &siparticles);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case 'd':
        sscanf(optarg, "%lf", &perturbation);
        break;
      case 'm':
        sscanf(optarg, "%lf", &rho);
        break;
      case 'f':
        strcpy(outputFileNameExtension, optarg);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || siparticles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -n PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 27 cells, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair_sidm_density() and "
        "runner_doself_sidm_density()."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-p                 - Random fractional change in h, h=h*random(1,p)"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("DOSELF1 function called: %s", DOSELF1_NAME);
  message("DOPAIR1 function called: %s", DOPAIR1_NAME);
  message("Vector size: %d", VEC_SIZE);
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f", pow_dimension(h) * kernel_norm);
  message("Density target: rho = %f", rho);

  printf("\n");

  /* Build the infrastructure */
  struct space space;
  space.periodic = 1;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  struct sidm_props sidm_p;
  sidm_p.eta_neighbours = h;
  sidm_p.h_tolerance = 1e0;
  sidm_p.h_max = FLT_MAX;
  sidm_p.max_smoothing_iterations = 1;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.sidm_properties = &sidm_p;
  engine.nodeID = NODE_ID;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct runner runner;
  runner.e = &engine;

  /* Construct some cells */
  struct cell *cells[27];
  struct cell *main_cell;
  static long long sipartId = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        double offset[3] = {i * size, j * size, k * size};
        cells[i * 9 + j * 3 + k] = make_cell(siparticles, offset, size, h, rho,
                                             &sipartId, perturbation, h_pert);

        runner_do_drift_sipart(&runner, cells[i * 9 + j * 3 + k], 0);

        // runner_do_hydro_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0, 0,
        // 0,
        //                      0); // TODO:SIDM sorting
      }
    }
  }

  /* Store the main cell for future use */
  main_cell = cells[13];

  ticks timings[27];
  for (int i = 0; i < 27; i++) timings[i] = 0;

  ticks time = 0;
  for (size_t i = 0; i < runs; ++i) {
    /* Zero the fields */
    for (int j = 0; j < 27; ++j) zero_particle_fields(cells[j]);

    const ticks tic = getticks();

#if defined(TEST_DOSELF_SUBSET) || defined(TEST_DOPAIR_SUBSET)
    int *siid = NULL;
    int sicount = 0;
    if ((siid = (int *)malloc(sizeof(int) * main_cell->sidm.count)) == NULL)
      error("Can't allocate memory for siid.");
    for (int k = 0; k < main_cell->sidm.count; k++) {
      siid[sicount] = k;
      ++sicount;
    }
#endif

    /* Run all the pairs */
    for (int j = 0; j < 27; ++j) {
      if (cells[j] != main_cell) {
        const ticks sub_tic = getticks();

#ifdef TEST_DOPAIR_SUBSET
        DOPAIR1_SUBSET(&runner, main_cell, main_cell->sidm.parts, siid, sicount,
                       cells[j]);
#else
        DOPAIR1(&runner, main_cell, cells[j], /*limit_h_min=*/0,
                /*limit_h_max=*/0);
#endif

        timings[j] += getticks() - sub_tic;
      }
    }

    /* And now the self-interaction */
    const ticks self_tic = getticks();

#ifdef TEST_DOSELF_SUBSET
    DOSELF1_SUBSET(&runner, main_cell, main_cell->sidm.parts, siid, sicount);
#else
    DOSELF1(&runner, main_cell, /*limit_h_min=*/0, /*limit_h_max=*/0);
#endif

    timings[13] += getticks() - self_tic;

    const ticks toc = getticks();
    time += toc - tic;

    /* Let's get physical ! */
    end_calculation(main_cell);

    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_sidm_dopair_27_%.150s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell, cells);
    }
  }

  /* Output timing */
  ticks corner_time = timings[0] + timings[2] + timings[6] + timings[8] +
                      timings[18] + timings[20] + timings[24] + timings[26];

  ticks edge_time = timings[1] + timings[3] + timings[5] + timings[7] +
                    timings[9] + timings[11] + timings[15] + timings[17] +
                    timings[19] + timings[21] + timings[23] + timings[25];

  ticks face_time = timings[4] + timings[10] + timings[12] + timings[14] +
                    timings[16] + timings[22];

  ticks self_time = timings[13];

  message("Corner calculations took:     %.3f %s.",
          clocks_from_ticks(corner_time / runs), clocks_getunit());
  message("Edge calculations took:       %.3f %s.",
          clocks_from_ticks(edge_time / runs), clocks_getunit());
  message("Face calculations took:       %.3f %s.",
          clocks_from_ticks(face_time / runs), clocks_getunit());
  message("Self calculations took:       %.3f %s.",
          clocks_from_ticks(self_time / runs), clocks_getunit());
  message("SWIFT calculation took:       %.3f %s.",
          clocks_from_ticks(time / runs), clocks_getunit());

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  for (int i = 0; i < 27; ++i) zero_particle_fields(cells[i]);

  const ticks tic = getticks();

  /* Run all the brute-force pairs */
  for (int j = 0; j < 27; ++j)
    if (cells[j] != main_cell)
      pairs_all_sidm_density(&runner, main_cell, cells[j]);

  /* And now the self-interaction */
  self_all_sidm_density(&runner, main_cell);

  const ticks toc = getticks();

  /* Let's get physical ! */
  end_calculation(main_cell);

  /* Dump */
  sprintf(outputFileName, "sidm_brute_force_27_%.150s.dat",
          outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell, cells);

  /* Output timing */
  message("Brute force calculation took : %.3f %s.",
          clocks_from_ticks(toc - tic), clocks_getunit());

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return 0;
}
