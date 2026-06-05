/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Katy Proctor (katy.proctor@fysik.su.se).
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

#define NODE_ID 0

enum velocity_field {
  velocity_zero,
  velocity_const,
  velocity_divergent,
  velocity_rotating
};

void set_velocity(struct sipart *sipart, enum velocity_field vel, float size) {
  switch (vel) {
    case velocity_zero:
      sipart->v[0] = 0.f;
      sipart->v[1] = 0.f;
      sipart->v[2] = 0.f;
      break;
    case velocity_const:
      sipart->v[0] = 1.f;
      sipart->v[1] = 0.f;
      sipart->v[2] = 0.f;
      break;
    case velocity_divergent:
      sipart->v[0] = sipart->x[0] - 2.5 * size;
      sipart->v[1] = sipart->x[1] - 2.5 * size;
      sipart->v[2] = sipart->x[2] - 2.5 * size;
      break;
    case velocity_rotating:
      sipart->v[0] = sipart->x[1];
      sipart->v[1] = -sipart->x[0];
      sipart->v[2] = 0.f;
      break;
  }
}

void reset_siparticles(struct cell *c, enum velocity_field vel, float size) {
  for (int i = 0; i < c->sidm.count; ++i) {
    set_velocity(&c->sidm.parts[i], vel, size);
    sidm_init_sipart(&c->sidm.parts[i]);
  }
}

struct cell *make_cell(size_t n, const double offset[3], double size, double h,
                       double density, long long *sipartId, double pert,
                       double h_pert, struct sidm_props *sidm_properties,
                       enum velocity_field vel) {

  const size_t sicount = n * n * n;
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
  for (size_t x = 0; x < n; ++x) {
    for (size_t y = 0; y < n; ++y) {
      for (size_t z = 0; z < n; ++z) {
        sipart->x[0] =
            offset[0] +
            size * (x + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;
        sipart->x[1] =
            offset[1] +
            size * (y + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;
        sipart->x[2] =
            offset[2] +
            size * (z + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;

        set_velocity(sipart, vel, size);

        if (h_pert)
          sipart->h = size * h * random_uniform(1.f, h_pert) / (float)n;
        else
          sipart->h = size * h / (float)n;
        h_max = fmaxf(h_max, sipart->h);

        sidm_first_init_sipart(sipart, sidm_properties);

        sipart->id = ++(*sipartId);
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
  cell->sidm.h_max = h_max;
  cell->sidm.h_max_active = h_max;
  cell->sidm.count = sicount;
  // cell->grav.count = 0;
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
 * @brief Initializes all particle fields to be ready for a density calculation
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

void dump_particle_fields(char *fileName, struct cell *main_cell) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file, "# %4s %10s %10s %10s %13s %13s %13s %13s %13s\n", "ID",
          "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho", "SIDM_rate");

  fprintf(file, "# Main cell --------------------------------------------\n");

  /* Write main cell */
  for (int siid = 0; siid < main_cell->sidm.count; siid++) {
    fprintf(file, "%6llu %10f %10f %10f %13e %13e %13e %13e %13e\n",
            main_cell->sidm.parts[siid].id, main_cell->sidm.parts[siid].x[0],
            main_cell->sidm.parts[siid].x[1], main_cell->sidm.parts[siid].x[2],
            main_cell->sidm.parts[siid].v[0], main_cell->sidm.parts[siid].v[1],
            main_cell->sidm.parts[siid].v[2],
            sidm_get_comoving_density(&main_cell->sidm.parts[siid]),
            main_cell->sidm.parts[siid].SIDM_rate);
  }

  fclose(file);
}

/* Just forward declarations... */
void runner_dopair1_branch_sidm_density(struct runner *r, struct cell *ci,
                                        struct cell *cj, int limit_h_min,
                                        int limit_h_max);
void runner_doself1_branch_sidm_density(struct runner *r, struct cell *c,
                                        int limit_h_min, int limit_h_max);

void runner_dopair1_branch_sidm_force(struct runner *r, struct cell *ci,
                                      struct cell *cj, int limit_h_min,
                                      int limit_h_max);
void runner_doself1_branch_sidm_force(struct runner *r, struct cell *ci,
                                      int limit_h_min, int limit_h_max);
void runner_do_sidm_density_ghost(struct runner *r, struct cell *c, int timer);

/* And go... */
int main(int argc, char *argv[]) {

  /* Do not run test if code is not compiled in SIDM mode */
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
  enum velocity_field vel = velocity_zero;

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
  while ((c = getopt(argc, argv, "m:s:h:p:n:r:d:f:v:")) != -1) {
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
      case 'v':
        sscanf(optarg, "%d", (int *)&vel);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (optind >= argc) {
    error("Missing parameter file.");
  }

  if (h < 0 || siparticles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -n PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 125 cells, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_branch_sidm_density()"
        " and runner_doself1_branch_sidm_density() followed by "
        "runner_dopair1_sidm_force() and "
        "runner_doself1_sidm_force()"
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-p                 - Random fractional change in h, h=h*random(1,p)"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-v type (0,1,2,3)  - Velocity field: (zero, constant, divergent, "
        "rotating)"
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f", pow_dimension(h) * kernel_norm);
  message("Density target: rho = %f", rho);
  if (vel == velocity_zero)
    message("Velocity field: zero");
  else if (vel == velocity_const)
    message("Velocity field: constant");
  else if (vel == velocity_divergent)
    message("Velocity field: divergent");
  else if (vel == velocity_rotating)
    message("Velocity field: rotating");

  printf("\n");

  /* parse parameters */
  message("Reading parameters.");
  struct swift_params param_file;
  parser_read_file(argv[optind], &param_file);

  /* Default unit system */
  message("Initialization of the unit system.");
  struct unit_system us;
  units_init_cgs(&us);

  /* Default physical constants */
  message("Initialization of the physical constants.");
  struct phys_const prog_const;
  phys_const_init(&us, &param_file, &prog_const);

  struct output_options *output_options =
      (struct output_options *)malloc(sizeof(struct output_options));
  output_options_init(&param_file, 0, output_options);

  struct hydro_props hydro_properties;
  hydro_props_init(&hydro_properties, &prog_const, &us, &param_file);

  /* Build the infrastructure */
  struct space space;
  space.periodic = 1;
  space.dim[0] = 5.;
  space.dim[1] = 5.;
  space.dim[2] = 5.;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.nodeID = NODE_ID;

  struct runner runner;
  runner.e = &engine;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct sidm_props sidm_p;
  sidm_props_init(&sidm_p, &prog_const, &us, &param_file, &hydro_properties,
                  &cosmo);
  sidm_p.eta_neighbours = h;
  sidm_p.h_tolerance = 1e0;
  sidm_p.h_max = FLT_MAX;
  sidm_p.h_min = 0.f;
  sidm_p.h_min_ratio = 0.f;
  sidm_p.max_smoothing_iterations = 10;

  engine.sidm_properties = &sidm_p;

  /* Construct some cells */
  struct cell *cells[125];
  struct cell *inner_cells[27];
  struct cell *main_cell;
  int count = 0;
  static long long sipartId = 0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 5; ++k) {

        /* Position of the cell */
        const double offset[3] = {i * size, j * size, k * size};

        /* Construct it */
        cells[i * 25 + j * 5 + k] =
            make_cell(siparticles, offset, size, h, rho, &sipartId,
                      perturbation, h_pert, &sidm_p, vel);
        /* Store the inner cells */
        if (i > 0 && i < 4 && j > 0 && j < 4 && k > 0 && k < 4) {
          inner_cells[count] = cells[i * 25 + j * 5 + k];
          count++;
        }
      }
    }
  }

  /* Store the main cell for future use */
  main_cell = cells[62];

  ticks timings[27];
  for (int i = 0; i < 27; i++) timings[i] = 0;

  /* Start the test */
  ticks time = 0;
  for (size_t n = 0; n < runs; ++n) {
    const ticks tic = getticks();

    /* Initialise the particles */
    for (int j = 0; j < 125; ++j) runner_do_drift_sipart(&runner, cells[j], 0);

    /* Zero the density and SIDM_rate fields */
    for (int j = 0; j < 125; ++j) zero_particle_fields(cells[j]);

    /* Do the density calculation */

    /* Run all the pairs (only once !) */
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        for (int k = 0; k < 5; k++) {

          struct cell *ci = cells[i * 25 + j * 5 + k];

          for (int ii = -1; ii < 2; ii++) {
            int iii = i + ii;
            if (iii < 0 || iii >= 5) continue;
            iii = (iii + 5) % 5;
            for (int jj = -1; jj < 2; jj++) {
              int jjj = j + jj;
              if (jjj < 0 || jjj >= 5) continue;
              jjj = (jjj + 5) % 5;
              for (int kk = -1; kk < 2; kk++) {
                int kkk = k + kk;
                if (kkk < 0 || kkk >= 5) continue;
                kkk = (kkk + 5) % 5;

                struct cell *cj = cells[iii * 25 + jjj * 5 + kkk];

                if (cj > ci)
                  runner_dopair1_branch_sidm_density(
                      &runner, ci, cj, /*limit_h_min=*/0, /*limit_h_max=*/0);
              }
            }
          }
        }
      }
    }

    /* And now the self-interaction for the central cells */
    for (int j = 0; j < 27; ++j)
      runner_doself1_branch_sidm_density(&runner, inner_cells[j],
                                         /*limit_h_min=*/0, /*limit_h_max=*/0);

    /* Ghost to finish everything on the central cells */
    for (int j = 0; j < 27; ++j)
      runner_do_sidm_density_ghost(&runner, inner_cells[j], /*timer=*/0);

    /* Do the force calculation */

    int ctr = 0;
    /* Do the pairs (for the central 27 cells) */
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        for (int k = 1; k < 4; k++) {

          struct cell *cj = cells[i * 25 + j * 5 + k];

          if (main_cell != cj) {
            const ticks sub_tic = getticks();

            runner_dopair1_branch_sidm_force(&runner, main_cell, cj,
                                             /*limit_h_min=*/0,
                                             /*limit_h_max=*/0);

            timings[ctr++] += getticks() - sub_tic;
          }
        }
      }
    }

    const ticks self_tic = getticks();

    /* And now the self-interaction for the main cell */
    runner_doself1_branch_sidm_force(&runner, main_cell, /*limit_h_min=*/0,
                                     /*limit_h_max=*/0);

    timings[26] += getticks() - self_tic;

    /* Finally, end the force loop */
    runner_do_end_sidm_force(&runner, main_cell, 0);
    const ticks toc = getticks();
    time += toc - tic;

    /* Dump if necessary */
    if (n == 0) {
      sprintf(outputFileName, "swift_sidm_dopair_125_%.150s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell);
    }

    for (int i = 0; i < 125; ++i) {
      for (int pid = 0; pid < cells[i]->sidm.count; ++pid) {
        sidm_init_sipart(&cells[i]->sidm.parts[pid]);
      }
    }
  }

  /* Output timing */
  ticks corner_time = timings[0] + timings[2] + timings[6] + timings[8] +
                      timings[17] + timings[19] + timings[23] + timings[25];

  ticks edge_time = timings[1] + timings[3] + timings[5] + timings[7] +
                    timings[9] + timings[11] + timings[14] + timings[16] +
                    timings[18] + timings[20] + timings[22] + timings[24];

  ticks face_time = timings[4] + timings[10] + timings[12] + timings[13] +
                    timings[15] + timings[21];

  ticks self_time = timings[26];

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

  /* NOW BRUTE-FORCE CALCULATION */

  /* Reset particles to the same velocity field as the optimised run */
  for (int j = 0; j < 125; ++j) reset_siparticles(cells[j], vel, size);

  const ticks tic = getticks();

  /* Do the density calculation */

  /* Run all the pairs (only once !) */
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      for (int k = 0; k < 5; k++) {

        struct cell *ci = cells[i * 25 + j * 5 + k];

        for (int ii = -1; ii < 2; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= 5) continue;
          iii = (iii + 5) % 5;
          for (int jj = -1; jj < 2; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= 5) continue;
            jjj = (jjj + 5) % 5;
            for (int kk = -1; kk < 2; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= 5) continue;
              kkk = (kkk + 5) % 5;

              struct cell *cj = cells[iii * 25 + jjj * 5 + kkk];

              if (cj > ci) pairs_all_sidm_density(&runner, ci, cj);
            }
          }
        }
      }
    }
  }

  /* And now the self-interaction for the central cells */
  for (int j = 0; j < 27; ++j) self_all_sidm_density(&runner, inner_cells[j]);

  /* Ghost to finish everything on the central cells */
  for (int j = 0; j < 27; ++j)
    runner_do_sidm_density_ghost(&runner, inner_cells[j], /*timer=*/0);

  /* Do the force calculation */

  /* Do the pairs (for the central 27 cells) */
  for (int i = 1; i < 4; i++) {
    for (int j = 1; j < 4; j++) {
      for (int k = 1; k < 4; k++) {

        struct cell *cj = cells[i * 25 + j * 5 + k];

        if (main_cell != cj) pairs_all_sidm_force(&runner, main_cell, cj);
      }
    }
  }

  /* And now the self-interaction for the main cell */
  self_all_sidm_force(&runner, main_cell);

  /* Finally, end the force loop */
  runner_do_end_sidm_force(&runner, main_cell, 0);

  const ticks toc = getticks();

  /* Output timing */
  message("Brute force calculation took : %.3f %s.",
          clocks_from_ticks(toc - tic), clocks_getunit());

  sprintf(outputFileName, "sidm_brute_force_125_%.150s.dat",
          outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 125; ++i) clean_up(cells[i]);

  return 0;
}