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

#include <fenv.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "swift.h"

/**
 * Returns a random number (uniformly distributed) in [a,b[
 */
double random_uniform(double a, double b) {
  return (rand() / (double)RAND_MAX) * (b - a) + a;
}

/* n is both particles per axis and box size:
 * particles are generated on a mesh with unit spacing
 */
struct cell *make_cell(size_t n, double *offset, double size, double h,
                       double density, long long *partId, double pert) {
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
        // part->v[0] = part->x[0] - 1.5;
        // part->v[1] = part->x[1] - 1.5;
        // part->v[2] = part->x[2] - 1.5;
        part->v[0] = random_uniform(-0.05, 0.05);
        part->v[1] = random_uniform(-0.05, 0.05);
        part->v[2] = random_uniform(-0.05, 0.05);
        part->h = size * h / (float)n;
        part->id = ++(*partId);
        part->mass = density * volume / count;
        part->ti_begin = 0;
        part->ti_end = 1;
        ++part;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->dx_max = 0.;
  cell->h[0] = size;
  cell->h[1] = size;
  cell->h[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->ti_end_min = 1;
  cell->ti_end_max = 1;

  cell->sorted = 0;
  cell->sort = NULL;
  cell->sortsize = 0;
  runner_dosort(NULL, cell, 0x1FFF, 0);

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

  for (size_t pid = 0; pid < c->count; pid++) {
    c->parts[pid].rho = 0.f;
    c->parts[pid].rho_dh = 0.f;
    hydro_init_part(&c->parts[pid]);
  }
}

/**
 * @brief Ends the loop by adding the appropriate coefficients
 */
void end_calculation(struct cell *c) {

  for (size_t pid = 0; pid < c->count; pid++) {
    hydro_end_density(&c->parts[pid], 1);
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
#ifdef GADGET2_SPH
            main_cell->parts[pid].div_v, main_cell->parts[pid].density.rot_v[0],
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

        for (size_t pjd = 0; pjd < cj->count; pjd++) {
          fprintf(
              file,
              "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
              "%13e %13e %13e\n",
              cj->parts[pjd].id, cj->parts[pjd].x[0], cj->parts[pjd].x[1],
              cj->parts[pjd].x[2], cj->parts[pjd].v[0], cj->parts[pjd].v[1],
              cj->parts[pjd].v[2], cj->parts[pjd].rho, cj->parts[pjd].rho_dh,
              cj->parts[pjd].density.wcount, cj->parts[pjd].density.wcount_dh,
#ifdef GADGET2_SPH
              cj->parts[pjd].div_v, cj->parts[pjd].density.rot_v[0],
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

/* And go... */
int main(int argc, char *argv[]) {

  size_t runs = 0, particles = 0;
  double h = 1.12575, size = 1., rho = 1.;
  double perturbation = 0.;
  char outputFileNameExtension[200] = "";
  char outputFileName[200] = "";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "m:s:h:p:r:t:d:f:")) != -1) {
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
        "\n-h DISTANCE=1.1255 - Smoothing length"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("Smoothing length: h = %f", h);
  message("Neighbour target: N = %f", kernel_nwneigh);

  /* Build the infrastructure */
  struct space space;
  space.periodic = 0;
  space.h_max = h;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 1;

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
        cells[i * 9 + j * 3 + k] =
            make_cell(particles, offset, size, h, rho, &partId, perturbation);
      }
    }
  }

  /* Store the main cell for future use */
  main_cell = cells[13];

  ticks time = 0;
  for (size_t i = 0; i < runs; ++i) {

    /* Zero the fields */
    for (int j = 0; j < 27; ++j) zero_particle_fields(cells[j]);

    const ticks tic = getticks();

    /* Run all the pairs */
    for (int j = 0; j < 27; ++j)
      if (cells[j] != main_cell)
        runner_dopair1_density(&runner, main_cell, cells[j]);

    /* And now the self-interaction */
    runner_doself1_density(&runner, main_cell);

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

  /* Output timing */
  message("SWIFT calculation took       : %15lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  for (int i = 0; i < 27; ++i) zero_particle_fields(cells[i]);

  const ticks tic = getticks();

  /* Run all the brute-force pairs */
  for (int j = 0; j < 27; ++j)
    if (cells[j] != main_cell) pairs_all_density(&runner, main_cell, cells[j]);

  /* And now the self-interaction */
  self_all_density(&runner, main_cell);

  const ticks toc = getticks();

  /* Let's get physical ! */
  end_calculation(main_cell);

  /* Dump */
  sprintf(outputFileName, "brute_force_27_%s.dat", outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell, cells);

  /* Output timing */
  message("Brute force calculation took : %15lli ticks.", toc - tic);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return 0;
}
