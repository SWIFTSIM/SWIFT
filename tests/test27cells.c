/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

#define ACC_THRESHOLD 1e-5

#if defined(WITH_VECTORIZATION)
#define DOSELF1 runner_doself1_density_vec
#define DOPAIR1 runner_dopair1_density_vec
#define DOSELF1_NAME "runner_doself1_density_vec"
#define DOPAIR1_NAME "runner_dopair1_density_vec"
#endif

#ifndef DOSELF1
#define DOSELF1 runner_doself1_density
#define DOSELF1_NAME "runner_doself1_density"
#endif

#ifndef DOPAIR1
#define DOPAIR1 runner_dopair1_density
#define DOPAIR1_NAME "runner_dopair1_density"
#endif

enum velocity_types {
  velocity_zero,
  velocity_random,
  velocity_divergent,
  velocity_rotating
};

/**
 * @brief Constructs a cell and all of its particle in a valid state prior to
 * a DOPAIR or DOSELF calcuation.
 *
 * @param n The cube root of the number of particles.
 * @param offset The position of the cell offset from (0,0,0).
 * @param size The cell size.
 * @param h The smoothing length of the particles in units of the inter-particle
 *separation.
 * @param density The density of the fluid.
 * @param partId The running counter of IDs.
 * @param pert The perturbation to apply to the particles in the cell in units
 *of the inter-particle separation.
 * @param vel The type of velocity field (0, random, divergent, rotating)
 */
struct cell *make_cell(size_t n, double *offset, double size, double h,
                       double density, long long *partId, double pert,
                       enum velocity_types vel) {
  const size_t count = n * n * n;
  const double volume = size * size * size;
  struct cell *cell = malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->parts, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->parts, count * sizeof(struct part));

  /* Construct the parts */
  struct part *part = cell->parts;
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
        switch (vel) {
          case velocity_zero:
            part->v[0] = 0.f;
            part->v[1] = 0.f;
            part->v[2] = 0.f;
            break;
          case velocity_random:
            part->v[0] = random_uniform(-0.05, 0.05);
            part->v[1] = random_uniform(-0.05, 0.05);
            part->v[2] = random_uniform(-0.05, 0.05);
            break;
          case velocity_divergent:
            part->v[0] = part->x[0] - 1.5 * size;
            part->v[1] = part->x[1] - 1.5 * size;
            part->v[2] = part->x[2] - 1.5 * size;
            break;
          case velocity_rotating:
            part->v[0] = part->x[1];
            part->v[1] = -part->x[0];
            part->v[2] = 0.f;
            break;
        }
        part->h = size * h / (float)n;
        part->id = ++(*partId);

#if defined(GIZMO_SPH) || defined(SHADOWFAX_SPH)
        part->conserved.mass = density * volume / count;

#ifdef SHADOWFAX_SPH
        double anchor[3] = {0., 0., 0.};
        double side[3] = {1., 1., 1.};
        voronoi_cell_init(&part->cell, part->x, anchor, side);
#endif

#else
        part->mass = density * volume / count;
#endif

#if defined(HOPKINS_PE_SPH)
        part->entropy = 1.f;
        part->entropy_one_over_gamma = 1.f;
#endif

        part->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        part->ti_drift = 8;
        part->ti_kick = 8;
#endif

        ++part;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->dx_max = 0.;
  cell->dx_max_sort = 0.;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->ti_old = 8;
  cell->ti_end_min = 8;
  cell->ti_end_max = 8;
  cell->ti_sort = 8;

  shuffle_particles(cell->parts, cell->count);

  cell->sorted = 0;
  cell->sort = NULL;
  cell->sortsize = 0;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->parts);
  free(ci->sort);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields(struct cell *c) {
  for (int pid = 0; pid < c->count; pid++) {
    hydro_init_part(&c->parts[pid], NULL);
  }
}

/**
 * @brief Ends the loop by adding the appropriate coefficients
 */
void end_calculation(struct cell *c) {
  for (int pid = 0; pid < c->count; pid++) {
    hydro_end_density(&c->parts[pid]);
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *main_cell,
                          struct cell **cells) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s %10s %13s %13s %13s %13s %13s "
          "%13s %13s %13s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho", "rho_dh",
          "wcount", "wcount_dh", "div_v", "curl_vx", "curl_vy", "curl_vz");

  fprintf(file, "# Main cell --------------------------------------------\n");

  /* Write main cell */
  for (int pid = 0; pid < main_cell->count; pid++) {
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
            "%13e %13e %13e\n",
            main_cell->parts[pid].id, main_cell->parts[pid].x[0],
            main_cell->parts[pid].x[1], main_cell->parts[pid].x[2],
            main_cell->parts[pid].v[0], main_cell->parts[pid].v[1],
            main_cell->parts[pid].v[2],
            hydro_get_density(&main_cell->parts[pid]),
#if defined(GIZMO_SPH) || defined(SHADOWFAX_SPH)
            0.f,
#else
            main_cell->parts[pid].density.rho_dh,
#endif
            main_cell->parts[pid].density.wcount,
            main_cell->parts[pid].density.wcount_dh,
#if defined(GADGET2_SPH) || defined(DEFAULT_SPH) || defined(HOPKINS_PE_SPH)
            main_cell->parts[pid].density.div_v,
            main_cell->parts[pid].density.rot_v[0],
            main_cell->parts[pid].density.rot_v[1],
            main_cell->parts[pid].density.rot_v[2]
#else
            0., 0., 0., 0.
#endif
            );
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

        for (int pjd = 0; pjd < cj->count; pjd++) {
          fprintf(
              file,
              "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
              "%13e %13e %13e\n",
              cj->parts[pjd].id, cj->parts[pjd].x[0], cj->parts[pjd].x[1],
              cj->parts[pjd].x[2], cj->parts[pjd].v[0], cj->parts[pjd].v[1],
              cj->parts[pjd].v[2], hydro_get_density(&cj->parts[pjd]),
#if defined(GIZMO_SPH) || defined(SHADOWFAX_SPH)
              0.f,
#else
              main_cell->parts[pjd].density.rho_dh,
#endif
              cj->parts[pjd].density.wcount, cj->parts[pjd].density.wcount_dh,
#if defined(GADGET2_SPH) || defined(DEFAULT_SPH) || defined(HOPKINS_PE_SPH)
              cj->parts[pjd].density.div_v, cj->parts[pjd].density.rot_v[0],
              cj->parts[pjd].density.rot_v[1], cj->parts[pjd].density.rot_v[2]
#else
              0., 0., 0., 0.
#endif
              );
        }
      }
    }
  }
  fclose(file);
}

/**
 * @brief Compares the vectorised result against
 * the serial result of the interaction.
 *
 * @param serial_parts Particle array that has been interacted serially
 * @param vec_parts Particle array to be interacted using vectors
 * @param count No. of particles that have been interacted
 * @param threshold Level of accuracy needed
 *
 * @return Non-zero value if difference found, 0 otherwise
 */
int check_results(struct part *serial_parts, struct part *vec_parts, int count,
                  double threshold) {
  int result = 0;

  for (int i = 0; i < count; i++)
    result += compare_particles(serial_parts[i], vec_parts[i], threshold);

  return result;
}

/* Just a forward declaration... */
void runner_doself1_density(struct runner *r, struct cell *ci);
void runner_doself1_density_vec(struct runner *r, struct cell *ci);
void runner_dopair1_density(struct runner *r, struct cell *ci, struct cell *cj);
void runner_dopair1_density_vec(struct runner *r, struct cell *ci,
                                struct cell *cj);

/* And go... */
int main(int argc, char *argv[]) {

  engine_pin();
  size_t runs = 0, particles = 0;
  double h = 1.23485, size = 1., rho = 1.;
  double perturbation = 0.;
  double threshold = ACC_THRESHOLD;
  char outputFileNameExtension[200] = "";
  char outputFileName[200] = "";
  enum velocity_types vel = velocity_zero;

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FP-exceptions */
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "m:s:h:n:r:t:d:f:v:a:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
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
      case 'm':
        sscanf(optarg, "%lf", &rho);
        break;
      case 'f':
        strcpy(outputFileNameExtension, optarg);
        break;
      case 'v':
        sscanf(optarg, "%d", (int *)&vel);
        break;
      case 'a':
        sscanf(optarg, "%lf", &threshold);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || particles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -n PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 27 cells, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density() and "
        "runner_doself1_density()."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-v type (0,1,2,3)  - Velocity field: (zero, random, divergent, "
        "rotating)"
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("DOSELF1 function called: %s", DOSELF1_NAME);
  message("DOPAIR1 function called: %s", DOPAIR1_NAME);
  message("Vector size: %d", VEC_SIZE);
  message("Adiabatic index: ga = %f", hydro_gamma);
  message("Hydro implementation: %s", SPH_IMPLEMENTATION);
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f", pow_dimension(h) * kernel_norm);
  message("Density target: rho = %f", rho);
  message("div_v target:   div = %f", vel == 2 ? 3.f : 0.f);
  message("curl_v target: curl = [0., 0., %f]", vel == 3 ? -2.f : 0.f);

  printf("\n");

  /* Build the infrastructure */
  struct space space;
  space.periodic = 1;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  struct hydro_props hp;
  hp.h_max = FLT_MAX;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;

  struct runner runner;
  runner.e = &engine;

  /* Construct some cells */
  struct cell *cells[27];
  struct cell *main_cell;
  static long long partId = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        double offset[3] = {i * size, j * size, k * size};
        cells[i * 9 + j * 3 + k] = make_cell(particles, offset, size, h, rho,
                                             &partId, perturbation, vel);

        runner_do_drift_particles(&runner, cells[i * 9 + j * 3 + k], 0);

        runner_do_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0);
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

#if !(defined(MINIMAL_SPH) && defined(WITH_VECTORIZATION))

#ifdef WITH_VECTORIZATION
    runner.ci_cache.count = 0;
    cache_init(&runner.ci_cache, 512);
    runner.cj_cache.count = 0;
    cache_init(&runner.cj_cache, 512);
#endif

    /* Run all the pairs */
    for (int j = 0; j < 27; ++j) {
      if (cells[j] != main_cell) {
        const ticks sub_tic = getticks();

        DOPAIR1(&runner, main_cell, cells[j]);

        const ticks sub_toc = getticks();
        timings[j] += sub_toc - sub_tic;
      }
    }

    /* And now the self-interaction */
    const ticks self_tic = getticks();

    DOSELF1(&runner, main_cell);

    const ticks self_toc = getticks();

    timings[13] += self_toc - self_tic;

#endif

    const ticks toc = getticks();
    time += toc - tic;

    /* Let's get physical ! */
    end_calculation(main_cell);

    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_dopair_27_%s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell, cells);
    }
  }

  /* Store the vectorised particle results. */
  struct part vec_parts[main_cell->count];
  for (int i = 0; i < main_cell->count; i++) vec_parts[i] = main_cell->parts[i];

  /* Output timing */
  ticks corner_time = timings[0] + timings[2] + timings[6] + timings[8] +
                      timings[18] + timings[20] + timings[24] + timings[26];

  ticks edge_time = timings[1] + timings[3] + timings[5] + timings[7] +
                    timings[9] + timings[11] + timings[15] + timings[17] +
                    timings[19] + timings[21] + timings[23] + timings[25];

  ticks face_time = timings[4] + timings[10] + timings[12] + timings[14] +
                    timings[16] + timings[22];

  message("Corner calculations took       : %15lli ticks.", corner_time / runs);
  message("Edge calculations took         : %15lli ticks.", edge_time / runs);
  message("Face calculations took         : %15lli ticks.", face_time / runs);
  message("Self calculations took         : %15lli ticks.", timings[13] / runs);
  message("SWIFT calculation took         : %15lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  for (int i = 0; i < 27; ++i) zero_particle_fields(cells[i]);

  const ticks tic = getticks();

#if !(defined(MINIMAL_SPH) && defined(WITH_VECTORIZATION))

  /* Run all the brute-force pairs */
  for (int j = 0; j < 27; ++j)
    if (cells[j] != main_cell) pairs_all_density(&runner, main_cell, cells[j]);

  /* And now the self-interaction */
  self_all_density(&runner, main_cell);

#endif

  const ticks toc = getticks();

  /* Let's get physical ! */
  end_calculation(main_cell);

  /* Dump */
  sprintf(outputFileName, "brute_force_27_%s.dat", outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell, cells);

  /* Check serial results against the vectorised results. */
  if (check_results(main_cell->parts, vec_parts, main_cell->count, threshold))
    message("Differences found...");

  /* Output timing */
  message("Brute force calculation took : %15lli ticks.", toc - tic);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return 0;
}
