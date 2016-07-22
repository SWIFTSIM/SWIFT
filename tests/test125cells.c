/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/**
 * @brief Constructs a cell and all of its particle in a valid state prior to
 * a SPH time-step.
 *
 * @param n The cube root of the number of particles.
 * @param offset The position of the cell offset from (0,0,0).
 * @param size The cell size.
 * @param h The smoothing length of the particles in units of the inter-particle
 * separation.
 * @param density The density of the fluid.
 * @param partId The running counter of IDs.
 * @param pert The perturbation to apply to the particles in the cell in units
 * of the inter-particle separation.
 */
struct cell *make_cell(size_t n, const double offset[3], double size, double h,
                       double density, long long *partId, double pert) {

  const size_t count = n * n * n;
  const double volume = size * size * size;
  struct cell *cell = malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->parts, part_align,
                     count * sizeof(struct part)) != 0)
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  if (posix_memalign((void **)&cell->xparts, xpart_align,
                     count * sizeof(struct xpart)) != 0)
    error("couldn't allocate particles, no. of x-particles: %d", (int)count);
  bzero(cell->parts, count * sizeof(struct part));
  bzero(cell->xparts, count * sizeof(struct xpart));

  /* Construct the parts */
  struct part *part = cell->parts;
  struct xpart *xpart = cell->xparts;
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
        part->h = size * h / (float)n;
        part->entropy = 1.f;
        part->id = ++(*partId);
        part->mass = density * volume / count;
        part->ti_begin = 0;
        part->ti_end = 1;
        xpart->v_full[0] = part->v[0];
        xpart->v_full[1] = part->v[1];
        xpart->v_full[2] = part->v[2];
        ++part;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->gcount = 0;
  cell->dx_max = 0.;
  cell->h[0] = size;
  cell->h[1] = size;
  cell->h[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->ti_end_min = 1;
  cell->ti_end_max = 1;

  // shuffle_particles(cell->parts, cell->count);

  cell->sorted = 0;
  cell->sort = NULL;
  cell->sortsize = 0;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->parts);
  free(ci->xparts);
  free(ci->sort);
  free(ci);
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
  for (size_t pid = 0; pid < main_cell->count; pid++) {
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
            "%13e %13e %13e\n",
            main_cell->parts[pid].id, main_cell->parts[pid].x[0],
            main_cell->parts[pid].x[1], main_cell->parts[pid].x[2],
            main_cell->parts[pid].v[0], main_cell->parts[pid].v[1],
            main_cell->parts[pid].v[2], main_cell->parts[pid].rho,
            main_cell->parts[pid].rho_dh, main_cell->parts[pid].density.wcount,
            main_cell->parts[pid].density.wcount_dh,
#if defined(GADGET2_SPH)
            main_cell->parts[pid].div_v, main_cell->parts[pid].density.rot_v[0],
            main_cell->parts[pid].density.rot_v[1],
            main_cell->parts[pid].density.rot_v[2]
#elif defined(DEFAULT_SPH)
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
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 5; ++k) {
        struct cell *cj = cells[i * 25 + j * 5 + k];
        if (cj == main_cell) continue;

        fprintf(file,
                "# Offset: [%2d %2d %2d] -----------------------------------\n",
                i - 2, j - 2, k - 2);

        for (size_t pjd = 0; pjd < cj->count; pjd++) {
          fprintf(
              file,
              "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
              "%13e %13e %13e\n",
              cj->parts[pjd].id, cj->parts[pjd].x[0], cj->parts[pjd].x[1],
              cj->parts[pjd].x[2], cj->parts[pjd].v[0], cj->parts[pjd].v[1],
              cj->parts[pjd].v[2], cj->parts[pjd].rho, cj->parts[pjd].rho_dh,
              cj->parts[pjd].density.wcount, cj->parts[pjd].density.wcount_dh,
#if defined(GADGET2_SPH)
              cj->parts[pjd].div_v, cj->parts[pjd].density.rot_v[0],
              cj->parts[pjd].density.rot_v[1], cj->parts[pjd].density.rot_v[2]
#elif defined(DEFAULT_SPH)
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

/* Just a forward declaration... */
void runner_dopair1_density(struct runner *r, struct cell *ci, struct cell *cj);
void runner_doself1_density(struct runner *r, struct cell *ci);
void runner_dopair2_force(struct runner *r, struct cell *ci, struct cell *cj);
void runner_doself2_force(struct runner *r, struct cell *ci);

/* And go... */
int main(int argc, char *argv[]) {

  size_t runs = 0, particles = 0;
  double h = 1.2348, size = 1., rho = 1.;
  double perturbation = 0.;
  char outputFileNameExtension[200] = "";
  char outputFileName[200] = "";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FP-exceptions */
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "m:s:h:p:r:t:d:f:v:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 's':
        sscanf(optarg, "%lf", &size);
        break;
      case 'p':
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
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || particles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -p PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates a cell pair, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "rotating)"
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f",
          h * h * h * 4.0 * M_PI * kernel_gamma3 / 3.0);
  message("Density target: rho = %f", rho);
  // message("div_v target:   div = %f", vel == 2 ? 3.f : 0.f);
  // message("curl_v target: curl = [0., 0., %f]", vel == 3 ? -2.f : 0.f);
  printf("\n");

  /* Build the infrastructure */
  struct space space;
  space.periodic = 0;
  space.h_max = h;

  struct hydro_props hp;
  hp.target_neighbours = h * h * h * 4.0 * M_PI * kernel_gamma3 / 3.0;
  hp.delta_neighbours = 1;
  hp.max_smoothing_iterations =
      1; /* We construct correct h values, 1 is enough */

  struct engine engine;
  engine.hydro_properties = &hp;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 1;

  struct runner runner;
  runner.e = &engine;

  /* Construct some cells */
  struct cell *cells[125];
  struct cell *inner_cells[27];
  struct cell *main_cell;
  int count = 0;
  static long long partId = 0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 5; ++k) {

        /* Position of the cell */
        const double offset[3] = {i * size, j * size, k * size};

        /* Construct it */
        cells[i * 25 + j * 5 + k] =
            make_cell(particles, offset, size, h, rho, &partId, perturbation);

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

  ticks time = 0;
  for (size_t i = 0; i < runs; ++i) {

    const ticks tic = getticks();

    /* First, sort stuff */
    for (int j = 0; j < 125; ++j) runner_do_sort(&runner, cells[j], 0x1FFF, 0);

    /* Initialise the particles */
    for (int j = 0; j < 125; ++j) runner_do_init(&runner, cells[j], 0);

/* Do the density calculation */
#if defined(DEFAULT_SPH) || !defined(WITH_VECTORIZATION)

    /* Run all the pairs (only once !)*/
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

                if (cj > ci) runner_dopair1_density(&runner, ci, cj);
              }
            }
          }
        }
      }
    }

    /* And now the self-interaction for the central cells*/
    for (int j = 0; j < 27; ++j)
      runner_doself1_density(&runner, inner_cells[j]);

#endif

    /* Ghost to finish everything on the central cells */
    for (int j = 0; j < 27; ++j) runner_do_ghost(&runner, inner_cells[j]);

    message("N_ngb = %f", main_cell->parts[0].density.wcount);

/* Do the force calculation */
#if defined(DEFAULT_SPH) || !defined(WITH_VECTORIZATION)

    /* Do the pairs (for the central 27 cells) */
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        for (int k = 1; k < 4; k++) {

          struct cell *cj = cells[i * 25 + j * 5 + k];

          if (main_cell != cj) runner_dopair2_force(&runner, main_cell, cj);
        }
      }
    }

    /* And now the self-interaction for the main cell */
    runner_doself2_force(&runner, main_cell);
#endif

    /* Finally, give a gentle kick */
    runner_do_kick(&runner, main_cell, 0);

    const ticks toc = getticks();
    time += toc - tic;

    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_dopair_125_%s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell, cells);
    }
  }

  /* Output timing */
  message("SWIFT calculation took       : %15lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /*   /\* Zero the fields *\/ */
  /*   for (int i = 0; i < 125; ++i) zero_particle_fields(cells[i]); */

  /*   const ticks tic = getticks(); */

  /* #if defined(DEFAULT_SPH) || !defined(WITH_VECTORIZATION) */

  /*   /\* Run all the brute-force pairs *\/ */
  /*   for (int j = 0; j < 125; ++j) */
  /*     if (cells[j] != main_cell) pairs_all_density(&runner, main_cell,
   * cells[j]); */

  /*   /\* And now the self-interaction *\/ */
  /*   self_all_density(&runner, main_cell); */

  /* #endif */

  /*   const ticks toc = getticks(); */

  /*   /\* Let's get physical ! *\/ */
  /*   end_calculation(main_cell); */

  /*   /\* Dump *\/ */
  /*   sprintf(outputFileName, "brute_force_125_%s.dat",
   * outputFileNameExtension); */
  /*   dump_particle_fields(outputFileName, main_cell, cells); */

  /*   /\* Output timing *\/ */
  /*   message("Brute force calculation took : %15lli ticks.", toc - tic); */

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 125; ++i) clean_up(cells[i]);

  return 0;
}
