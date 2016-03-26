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
                       double density, unsigned long long *partId,
                       double pert) {
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
  cell->h[0] = n;
  cell->h[1] = n;
  cell->h[2] = n;
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
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *ci, struct cell *cj) {

  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s %10s %13s %13s %13s %13s %13s "
          "%13s %13s %13s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho", "rho_dh",
          "wcount", "wcount_dh", "div_v", "curl_vx", "curl_vy", "curl_vz");

  fprintf(file, "# ci --------------------------------------------\n");

  for (size_t pid = 0; pid < ci->count; pid++) {
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
            "%13e %13e %13e\n",
            ci->parts[pid].id, ci->parts[pid].x[0], ci->parts[pid].x[1],
            ci->parts[pid].x[2], ci->parts[pid].v[0], ci->parts[pid].v[1],
            ci->parts[pid].v[2], ci->parts[pid].rho, ci->parts[pid].rho_dh,
            ci->parts[pid].density.wcount, ci->parts[pid].density.wcount_dh,
#ifdef GADGET2_SPH
            ci->parts[pid].div_v, ci->parts[pid].density.rot_v[0],
            ci->parts[pid].density.rot_v[1], ci->parts[pid].density.rot_v[2]
#else
            0., 0., 0., 0.
#endif
            );
  }

  fprintf(file, "# cj --------------------------------------------\n");

  for (size_t pjd = 0; pjd < cj->count; pjd++) {
    fprintf(file,
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

  fclose(file);
}

/* Just a forward declaration... */
void runner_dopair1_density(struct runner *r, struct cell *ci, struct cell *cj);

int main(int argc, char *argv[]) {
  size_t particles = 0, runs = 0, volume, type = 0;
  double offset[3] = {0, 0, 0}, h = 1.1255, size = 1., rho = 1.;
  double perturbation = 0.1;
  struct cell *ci, *cj;
  struct space space;
  struct engine engine;
  struct runner runner;
  char c;
  static unsigned long long partId = 0;
  char outputFileNameExtension[200] = "";
  char outputFileName[200] = "";
  ticks tic, toc, time;

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  srand(0);

  while ((c = getopt(argc, argv, "h:p:r:t:d:f:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 'p':
        sscanf(optarg, "%zu", &particles);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case 't':
        sscanf(optarg, "%zu", &type);
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

  if (h < 0 || particles == 0 || runs == 0 || type > 2) {
    printf(
        "\nUsage: %s -p PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates a cell pair, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density."
        "\n\nOptions:"
        "\n-t TYPE=0          - cells share face (0), edge (1) or corner (2)"
        "\n-h DISTANCE=1.1255 - smoothing length"
        "\n-d pert            - perturbation to apply to the particles [0,1["
        "\n-f fileName        - part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  space.periodic = 0;
  space.h_max = h;

  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 1;
  runner.e = &engine;

  volume = particles * particles * particles;
  message("particles: %zu B\npositions: 0 B", 2 * volume * sizeof(struct part));

  ci = make_cell(particles, offset, size, h, rho, &partId, perturbation);
  for (size_t i = 0; i < type + 1; ++i) offset[i] = 1.;
  cj = make_cell(particles, offset, size, h, rho, &partId, perturbation);

  time = 0;
  for (size_t i = 0; i < runs; ++i) {

    /* Zero the fields */
    zero_particle_fields(ci);
    zero_particle_fields(cj);

    tic = getticks();

    /* Run the test */
    runner_dopair1_density(&runner, ci, cj);

    toc = getticks();
    time += toc - tic;

    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_dopair_%s.dat", outputFileNameExtension);
      dump_particle_fields(outputFileName, ci, cj);
    }
  }

  /* Output timing */
  message("SWIFT calculation took       %lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  zero_particle_fields(ci);
  zero_particle_fields(cj);

  tic = getticks();

  /* Run the brute-force test */
  pairs_all_density(&runner, ci, cj);

  toc = getticks();

  /* Dump */
  sprintf(outputFileName, "brute_force_%s.dat", outputFileNameExtension);
  dump_particle_fields(outputFileName, ci, cj);

  /* Output timing */
  message("Brute force calculation took %lli ticks.", toc - tic);

  /* Clean things to make the sanitizer happy ... */
  clean_up(ci);
  clean_up(cj);

  return 0;
}
