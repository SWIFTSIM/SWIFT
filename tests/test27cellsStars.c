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

#define DOSELF1 runner_doself_branch_stars_density
#define DOSELF1_SUBSET runner_doself_subset_branch_stars_density
#ifdef TEST_DOSELF_SUBSET
#define DOSELF1_NAME "runner_doself_subset_branch_stars_density"
#else
#define DOSELF1_NAME "runner_doself1_branch_stars_density"
#endif

#define DOPAIR1 runner_dopair_branch_stars_density
#define DOPAIR1_SUBSET runner_dopair_subset_branch_stars_density
#ifdef TEST_DOPAIR_SUBSET
#define DOPAIR1_NAME "runner_dopair_subset_branch_stars_density"
#else
#define DOPAIR1_NAME "runner_dopair_branch_stars_density"
#endif

#define NODE_ID 0

/**
 * @brief Constructs a cell and all of its particle in a valid state prior to
 * a DOPAIR or DOSELF calcuation.
 *
 * @param n The cube root of the number of particles.
 * @param n_stars The cube root of the number of star particles.
 * @param offset The position of the cell offset from (0,0,0).
 * @param size The cell size.
 * @param h The smoothing length of the particles in units of the inter-particle
 * separation.
 * @param partId The running counter of IDs for gas.
 * @param spartId The running counter of IDs for stars.
 * @param pert The perturbation to apply to the particles in the cell in units
 * of the inter-particle separation.
 * @param h_pert The perturbation to apply to the smoothing length.
 */
struct cell *make_cell(size_t n, size_t n_stars, double *offset, double size,
                       double h, long long *partId, long long *spartId,
                       double pert, double h_pert) {
  const size_t count = n * n * n;
  const size_t scount = n_stars * n_stars * n_stars;
  float h_max = 0.f;
  struct cell *cell = (struct cell *)malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->hydro.parts, count * sizeof(struct part));

  /* Construct the parts */
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

        part->v[0] = 0;
        part->v[1] = 0;
        part->v[2] = 0;
        if (h_pert)
          part->h = size * h * random_uniform(1.f, h_pert) / (float)n;
        else
          part->h = size * h / (float)n;
        h_max = fmaxf(h_max, part->h);
        part->id = ++(*partId);

        part->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        part->ti_drift = 8;
        part->ti_kick = 8;
#endif
        ++part;
      }
    }
  }

  /* Construct the sparts */
  if (posix_memalign((void **)&cell->stars.parts, spart_align,
                     scount * sizeof(struct spart)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)scount);
  }
  bzero(cell->stars.parts, scount * sizeof(struct spart));

  struct spart *spart = cell->stars.parts;
  for (size_t x = 0; x < n_stars; ++x) {
    for (size_t y = 0; y < n_stars; ++y) {
      for (size_t z = 0; z < n_stars; ++z) {
        spart->x[0] =
            offset[0] + size * (x + 0.5 + random_uniform(-0.5, 0.5) * pert) /
                            (float)n_stars;
        spart->x[1] =
            offset[1] + size * (y + 0.5 + random_uniform(-0.5, 0.5) * pert) /
                            (float)n_stars;
        spart->x[2] =
            offset[2] + size * (z + 0.5 + random_uniform(-0.5, 0.5) * pert) /
                            (float)n_stars;

        spart->v[0] = 0;
        spart->v[1] = 0;
        spart->v[2] = 0;
        if (h_pert)
          spart->h = size * h * random_uniform(1.f, h_pert) / (float)n_stars;
        else
          spart->h = size * h / (float)n_stars;
        h_max = fmaxf(h_max, spart->h);
        spart->id = ++(*spartId);

        spart->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        spart->ti_drift = 8;
        spart->ti_kick = 8;
#endif
        ++spart;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->hydro.h_max = h_max;
  cell->hydro.count = count;
  cell->stars.count = scount;
  cell->hydro.dx_max_part = 0.;
  cell->hydro.dx_max_sort = 0.;
  cell->stars.dx_max_sort = 0.;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->hydro.ti_old_part = 8;
  cell->hydro.ti_end_min = 8;
  cell->hydro.ti_end_max = 8;
  cell->grav.ti_old_part = 8;
  cell->grav.ti_end_min = 8;
  cell->grav.ti_end_max = 8;
  cell->stars.ti_end_min = 8;
  cell->nodeID = NODE_ID;

  shuffle_particles(cell->hydro.parts, cell->hydro.count);
  shuffle_sparticles(cell->stars.parts, cell->stars.count);

  cell->hydro.sorted = 0;
  for (int k = 0; k < 13; k++) cell->hydro.sort[k] = NULL;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->hydro.parts);
  free(ci->stars.parts);
  for (int k = 0; k < 13; k++)
    if (ci->hydro.sort[k] != NULL) free(ci->hydro.sort[k]);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields(struct cell *c) {
  for (int pid = 0; pid < c->stars.count; pid++) {
    stars_init_spart(&c->stars.parts[pid]);
  }
}

/**
 * @brief Ends the loop by adding the appropriate coefficients
 */
void end_calculation(struct cell *c, const struct cosmology *cosmo) {
  for (int pid = 0; pid < c->stars.count; pid++) {
    stars_end_density(&c->stars.parts[pid], cosmo);

    /* Recover the common "Neighbour number" definition */
    c->stars.parts[pid].density.wcount *= pow_dimension(c->stars.parts[pid].h);
    c->stars.parts[pid].density.wcount *= kernel_norm;
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *main_cell,
                          struct cell **cells) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file, "# %4s %10s %10s %10s %13s %13s\n", "ID", "pos_x", "pos_y",
          "pos_z", "wcount", "wcount_dh");

  fprintf(file, "# Main cell --------------------------------------------\n");

  /* Write main cell */
  for (int pid = 0; pid < main_cell->stars.count; pid++) {
    fprintf(file, "%6llu %10f %10f %10f %13e %13e\n",
            main_cell->stars.parts[pid].id, main_cell->stars.parts[pid].x[0],
            main_cell->stars.parts[pid].x[1], main_cell->stars.parts[pid].x[2],
            main_cell->stars.parts[pid].density.wcount,
            main_cell->stars.parts[pid].density.wcount_dh);
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

        for (int pjd = 0; pjd < cj->stars.count; pjd++) {
          fprintf(file, "%6llu %10f %10f %10f %13e %13e\n",
                  cj->stars.parts[pjd].id, cj->stars.parts[pjd].x[0],
                  cj->stars.parts[pjd].x[1], cj->stars.parts[pjd].x[2],
                  cj->stars.parts[pjd].density.wcount,
                  cj->stars.parts[pjd].density.wcount_dh);
        }
      }
    }
  }
  fclose(file);
}

/* Just a forward declaration... */
void runner_dopair_branch_stars_density(struct runner *r, struct cell *ci,
                                        struct cell *cj);
void runner_doself_branch_stars_density(struct runner *r, struct cell *c);
void runner_dopair_subset_branch_stars_density(struct runner *r,
                                               struct cell *restrict ci,
                                               struct spart *restrict sparts_i,
                                               int *restrict ind, int scount,
                                               struct cell *restrict cj);
void runner_doself_subset_branch_stars_density(struct runner *r,
                                               struct cell *restrict ci,
                                               struct spart *restrict sparts,
                                               int *restrict ind, int scount);

/* And go... */
int main(int argc, char *argv[]) {

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  size_t runs = 0, particles = 0;
  size_t sparticles = 0;
  double h = 1.23485, size = 1.;
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

  char c;
  while ((c = getopt(argc, argv, "s:h:p:n:N:r:t:d:f:")) != -1) {
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
        sscanf(optarg, "%zu", &particles);
        break;
      case 'N':
        sscanf(optarg, "%zu", &sparticles);
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

  if (h < 0 || particles == 0 || runs == 0 || sparticles == 0) {
    printf(
        "\nUsage: %s -n PARTICLES_PER_AXIS -N SPARTICLES_PER_AXIS -r "
        "NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 27 cells, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair_stars_density() and "
        "runner_doself_stars_density()."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-p                 - Random fractional change in h, h=h*random(1,p)"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("DOSELF1 function called: %s", DOSELF1_NAME);
  message("DOPAIR1 function called: %s", DOPAIR1_NAME);
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f", pow_dimension(h) * kernel_norm);

  printf("\n");

  /* Build the infrastructure */
  struct space space;
  space.periodic = 1;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  struct hydro_props hp;
  hp.eta_neighbours = h;
  hp.h_tolerance = 1e0;
  hp.h_max = FLT_MAX;
  hp.max_smoothing_iterations = 1;
  hp.CFL_condition = 0.1;

  struct stars_props stars_p;
  stars_p.eta_neighbours = h;
  stars_p.h_tolerance = 1e0;
  stars_p.h_max = FLT_MAX;
  stars_p.max_smoothing_iterations = 1;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;
  engine.stars_properties = &stars_p;
  engine.nodeID = NODE_ID;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct runner runner;
  runner.e = &engine;

  /* Construct some cells */
  struct cell *cells[27];
  struct cell *main_cell;
  static long long partId = 0;
  long long spartId = particles * particles * particles * 27;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        double offset[3] = {i * size, j * size, k * size};
        cells[i * 9 + j * 3 + k] =
            make_cell(particles, sparticles, offset, size, h, &partId, &spartId,
                      perturbation, h_pert);

        runner_do_drift_part(&runner, cells[i * 9 + j * 3 + k], 0);

        runner_do_hydro_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0, 0);
        runner_do_stars_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0, 0);
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
    int *pid = NULL;
    int scount = 0;
    if ((pid = (int *)malloc(sizeof(int) * main_cell->stars.count)) == NULL)
      error("Can't allocate memory for pid.");
    for (int k = 0; k < main_cell->stars.count; k++)
      if (spart_is_active(&main_cell->stars.parts[k], &engine)) {
        pid[scount] = k;
        ++scount;
      }
#endif

    /* Run all the pairs */
    for (int j = 0; j < 27; ++j) {
      if (cells[j] != main_cell) {
        const ticks sub_tic = getticks();

#ifdef TEST_DOPAIR_SUBSET
        DOPAIR1_SUBSET(&runner, main_cell, main_cell->stars.parts, pid, scount,
                       cells[j]);
#else
        DOPAIR1(&runner, main_cell, cells[j]);
#endif

        timings[j] += getticks() - sub_tic;
      }
    }

    /* And now the self-interaction */
    const ticks self_tic = getticks();

#ifdef TEST_DOSELF_SUBSET
    DOSELF1_SUBSET(&runner, main_cell, main_cell->stars.parts, pid, scount);
#else
    DOSELF1(&runner, main_cell);
#endif

    timings[13] += getticks() - self_tic;

    const ticks toc = getticks();
    time += toc - tic;

    /* Let's get physical ! */
    end_calculation(main_cell, &cosmo);

    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_star_dopair_27_%.150s.dat",
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

  message("Corner calculations took       : %15lli ticks.", corner_time / runs);
  message("Edge calculations took         : %15lli ticks.", edge_time / runs);
  message("Face calculations took         : %15lli ticks.", face_time / runs);
  message("Self calculations took         : %15lli ticks.", timings[13] / runs);
  message("SWIFT calculation took         : %15lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  for (int i = 0; i < 27; ++i) zero_particle_fields(cells[i]);

  const ticks tic = getticks();

  /* Run all the brute-force pairs */
  for (int j = 0; j < 27; ++j)
    if (cells[j] != main_cell)
      pairs_all_stars_density(&runner, main_cell, cells[j]);

  /* And now the self-interaction */
  self_all_stars_density(&runner, main_cell);

  const ticks toc = getticks();

  /* Let's get physical ! */
  end_calculation(main_cell, &cosmo);

  /* Dump */
  sprintf(outputFileName, "star_brute_force_27_%.150s.dat",
          outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell, cells);

  /* Output timing */
  message("Brute force calculation took : %15lli ticks.", toc - tic);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return 0;
}
