/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Rob McGibbon (mcgibbon@strw.leidenuniv.nl)
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

/* This test requires a black hole model with star-density support (e.g.
 * EAGLE). With the 'none' model the BH-related cell fields are compacted
 * into a union and cannot be used. */
#ifdef BLACK_HOLES_HAVE_STAR_DENSITY

#include "black_holes_iact.h"
#include "runner_doiact_bh_stars.h"

#define DOSELF1 runner_doself_branch_bh_stars_density
#define DOSELF1_SUBSET runner_doself_subset_branch_bh_stars_density
#ifdef TEST_DOSELF_SUBSET
#define DOSELF1_NAME "runner_doself_subset_branch_bh_stars_density"
#else
#define DOSELF1_NAME "runner_doself_branch_bh_stars_density"
#endif

#define DOPAIR1 runner_dopair_branch_bh_stars_density
#define DOPAIR1_SUBSET runner_dopair_subset_branch_bh_stars_density
#ifdef TEST_DOPAIR_SUBSET
#define DOPAIR1_NAME "runner_dopair_subset_branch_bh_stars_density"
#else
#define DOPAIR1_NAME "runner_dopair_branch_bh_stars_density"
#endif

#define NODE_ID 0

/**
 * @brief Constructs a cell with star particles and black holes in a valid
 * state prior to a DOPAIR or DOSELF calculation.
 *
 * @param n_stars The cube root of the number of star particles.
 * @param n_bh The cube root of the number of black holes.
 * @param offset The position of the cell offset from (0,0,0).
 * @param size The cell size.
 * @param h The search radius of the black holes in units of the star
 * inter-particle separation.
 * @param spartId The running counter of IDs for stars.
 * @param bpartId The running counter of IDs for black holes.
 * @param pert The perturbation to apply to the particles in the cell in units
 * of the inter-particle separation.
 * @param h_pert The perturbation to apply to the search radius.
 */
struct cell *make_cell(size_t n_stars, size_t n_bh, double *offset, double size,
                       double h, long long *spartId, long long *bpartId,
                       double pert, double h_pert) {
  const size_t scount = n_stars * n_stars * n_stars;
  const size_t bcount = n_bh * n_bh * n_bh;
  float stars_h_max = 0.f;
  float h_star_max = 0.f;
  struct cell *cell = NULL;
  if (posix_memalign((void **)&cell, cell_align, sizeof(struct cell)) != 0) {
    error("Couldn't allocate the cell");
  }
  bzero(cell, sizeof(struct cell));

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
        spart->h = size * h / (float)n_stars;
        spart->mass = 1.f;
        stars_h_max = fmaxf(stars_h_max, spart->h);
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

  /* Construct the bparts */
  if (posix_memalign((void **)&cell->black_holes.parts, bpart_align,
                     bcount * sizeof(struct bpart)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)bcount);
  }
  bzero(cell->black_holes.parts, bcount * sizeof(struct bpart));

  struct bpart *bpart = cell->black_holes.parts;
  for (size_t x = 0; x < n_bh; ++x) {
    for (size_t y = 0; y < n_bh; ++y) {
      for (size_t z = 0; z < n_bh; ++z) {
        bpart->x[0] = offset[0] +
                      size * (x + 0.5 + random_uniform(-0.5, 0.5) * pert) /
                          (float)n_bh;
        bpart->x[1] = offset[1] +
                      size * (y + 0.5 + random_uniform(-0.5, 0.5) * pert) /
                          (float)n_bh;
        bpart->x[2] = offset[2] +
                      size * (z + 0.5 + random_uniform(-0.5, 0.5) * pert) /
                          (float)n_bh;

        bpart->v[0] = 0;
        bpart->v[1] = 0;
        bpart->v[2] = 0;

        /* Search radius set relative to the *star* spacing */
        if (h_pert)
          bpart->h_star = size * h * random_uniform(1.f, h_pert) /
                          (float)n_stars;
        else
          bpart->h_star = size * h / (float)n_stars;
        bpart->h = bpart->h_star;
        h_star_max = fmaxf(h_star_max, bpart->h_star);
        bpart->id = ++(*bpartId);

        bpart->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        bpart->ti_drift = 8;
        bpart->ti_kick = 8;
#endif
        ++bpart;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->stars.h_max = stars_h_max;
  cell->stars.h_max_active = stars_h_max;
  cell->stars.h_max_old = stars_h_max;
  cell->stars.count = scount;
  cell->stars.dx_max_part = 0.;
  cell->stars.dx_max_sort = 0.;
  cell->black_holes.count = bcount;
  cell->black_holes.h_max = h_star_max;
  cell->black_holes.h_max_active = h_star_max;
  cell->black_holes.h_max_old = h_star_max;
  cell->black_holes.h_star_max = h_star_max;
  cell->black_holes.h_star_max_active = h_star_max;
  cell->black_holes.h_star_max_old = h_star_max;
  cell->black_holes.dx_max_part = 0.;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->dmin = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->hydro.super = cell;
  cell->stars.ti_old_part = 8;
  cell->stars.ti_end_min = 8;
  cell->black_holes.ti_old_part = 8;
  cell->black_holes.ti_end_min = 8;
  cell->hydro.ti_old_part = 8;
  cell->hydro.ti_end_min = 8;
  cell->nodeID = NODE_ID;

  shuffle_sparticles(cell->stars.parts, cell->stars.count);

  cell->stars.sorted = 0;
  cell->stars.sort = NULL;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->stars.parts);
  free(ci->black_holes.parts);
  free(ci->stars.sort);
  free(ci);
}

/**
 * @brief Initializes all black holes to be ready for a star-density
 * calculation
 */
void zero_particle_fields(struct cell *c) {
  for (int bid = 0; bid < c->black_holes.count; bid++) {
    black_holes_init_stars_density(&c->black_holes.parts[bid]);
  }
}

/**
 * @brief Ends the loop by adding the appropriate coefficients
 */
void end_calculation(struct cell *c, const struct cosmology *cosmo) {
  for (int bid = 0; bid < c->black_holes.count; bid++) {
    struct bpart *bp = &c->black_holes.parts[bid];
    black_holes_stars_end_density(bp, cosmo);

    /* Recover the common "Neighbour number" definition */
    bp->stars_density.wcount *= pow_dimension(bp->h_star);
    bp->stars_density.wcount *= kernel_norm;
  }
}

/**
 * @brief Dump the black holes of the main cell to a file
 */
void dump_particle_fields(char *fileName, struct cell *main_cell) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file, "# %4s %10s %10s %10s %13s %13s %13s\n", "ID", "pos_x", "pos_y",
          "pos_z", "wcount", "wcount_dh", "rho");

  fprintf(file, "# Main cell --------------------------------------------\n");

  /* Write main cell */
  for (int bid = 0; bid < main_cell->black_holes.count; bid++) {
    const struct bpart *bp = &main_cell->black_holes.parts[bid];
    fprintf(file, "%6llu %10f %10f %10f %13e %13e %13e\n", bp->id, bp->x[0],
            bp->x[1], bp->x[2], bp->stars_density.wcount,
            bp->stars_density.wcount_dh, bp->stars_density.rho);
  }
  fclose(file);
}

/**
 * @brief Brute-force version of the self interaction.
 */
void self_all_bh_stars_density(struct runner *r, struct cell *ci) {

  for (int bid = 0; bid < ci->black_holes.count; bid++) {
    struct bpart *bi = &ci->black_holes.parts[bid];
    const float hi = bi->h_star;
    const float hig2 = hi * hi * kernel_gamma2;

    for (int sjd = 0; sjd < ci->stars.count; sjd++) {
      const struct spart *sj = &ci->stars.parts[sjd];

      const float dx[3] = {(float)(bi->x[0] - sj->x[0]),
                           (float)(bi->x[1] - sj->x[1]),
                           (float)(bi->x[2] - sj->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (r2 < hig2) {
        runner_iact_nonsym_bh_stars_density(r2, dx, hi, bi, sj);
      }
    }
  }
}

/**
 * @brief Brute-force version of the pair interaction (periodic wrapping).
 */
void pairs_all_bh_stars_density(struct runner *r, struct cell *ci,
                                struct cell *cj) {

  const double *dim = r->e->s->dim;

  for (int bid = 0; bid < ci->black_holes.count; bid++) {
    struct bpart *bi = &ci->black_holes.parts[bid];
    const float hi = bi->h_star;
    const float hig2 = hi * hi * kernel_gamma2;

    for (int sjd = 0; sjd < cj->stars.count; sjd++) {
      const struct spart *sj = &cj->stars.parts[sjd];

      double dxd[3] = {bi->x[0] - sj->x[0], bi->x[1] - sj->x[1],
                       bi->x[2] - sj->x[2]};
      for (int k = 0; k < 3; k++) {
        if (dxd[k] > dim[k] / 2)
          dxd[k] -= dim[k];
        else if (dxd[k] < -dim[k] / 2)
          dxd[k] += dim[k];
      }
      const float dx[3] = {(float)dxd[0], (float)dxd[1], (float)dxd[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (r2 < hig2) {
        runner_iact_nonsym_bh_stars_density(r2, dx, hi, bi, sj);
      }
    }
  }
}

/* And go... */
int main(int argc, char *argv[]) {

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  size_t runs = 0;
  size_t sparticles = 0, bparticles = 0;
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

  int c;
  while ((c = getopt(argc, argv, "s:h:p:N:B:r:d:f:")) != -1) {
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
      case 'N':
        sscanf(optarg, "%zu", &sparticles);
        break;
      case 'B':
        sscanf(optarg, "%zu", &bparticles);
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

  if (h < 0 || sparticles == 0 || bparticles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -N SPARTICLES_PER_AXIS -B BPARTICLES_PER_AXIS -r "
        "NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 27 cells, filled with stars and black holes on a "
        "Cartesian grid."
        "\nThese are then interacted using "
        "runner_dopair_branch_bh_stars_density() and "
        "runner_doself_branch_bh_stars_density()."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Search radius in units of <star spacing>"
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
  message("Search radius: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f", pow_dimension(h) * kernel_norm);

  printf("\n");

  /* Build the infrastructure */
  struct space space;
  bzero(&space, sizeof(struct space));
  space.periodic = 1;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  struct hydro_props hp;
  bzero(&hp, sizeof(struct hydro_props));
  hp.eta_neighbours = h;
  hp.h_tolerance = 1e0;
  hp.h_max = FLT_MAX;
  hp.max_smoothing_iterations = 1;
  hp.CFL_condition = 0.1;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;
  engine.nodeID = NODE_ID;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &engine;

  struct lightcone_array_props lightcone_array_properties;
  lightcone_array_properties.nr_lightcones = 0;
  engine.lightcone_array_properties = &lightcone_array_properties;

  /* Construct some cells */
  struct cell *cells[27];
  struct cell *main_cell;
  long long spartId = 0;
  long long bpartId = sparticles * sparticles * sparticles * 27;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        double offset[3] = {i * size, j * size, k * size};
        cells[i * 9 + j * 3 + k] =
            make_cell(sparticles, bparticles, offset, size, h, &spartId,
                      &bpartId, perturbation, h_pert);

        runner_do_drift_spart(&runner, cells[i * 9 + j * 3 + k], 0);

        runner_do_stars_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0, 0);
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

#if defined(TEST_DOSELF_SUBSET) || defined(TEST_DOPAIR_SUBSET)
    int *pid = NULL;
    int bcount = 0;
    if ((pid = (int *)malloc(sizeof(int) * main_cell->black_holes.count)) ==
        NULL)
      error("Can't allocate memory for pid.");
    for (int k = 0; k < main_cell->black_holes.count; k++)
      if (bpart_is_active(&main_cell->black_holes.parts[k], &engine)) {
        pid[bcount] = k;
        ++bcount;
      }
#endif

    /* Run all the pairs */
    for (int j = 0; j < 27; ++j) {
      if (cells[j] != main_cell) {

#ifdef TEST_DOPAIR_SUBSET
        DOPAIR1_SUBSET(&runner, main_cell, main_cell->black_holes.parts, pid,
                       bcount, cells[j]);
#else
        DOPAIR1(&runner, main_cell, cells[j]);
#endif
      }
    }

    /* And now the self-interaction */
#ifdef TEST_DOSELF_SUBSET
    DOSELF1_SUBSET(&runner, main_cell, main_cell->black_holes.parts, pid,
                   bcount);
#else
    DOSELF1(&runner, main_cell);
#endif

#if defined(TEST_DOSELF_SUBSET) || defined(TEST_DOPAIR_SUBSET)
    free(pid);
#endif

    const ticks toc = getticks();
    time += toc - tic;

    /* Let's get physical ! */
    end_calculation(main_cell, &cosmo);

    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_bh_stars_dopair_27_%.150s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell);
    }
  }

  /* Output timing */
  message("SWIFT calculation took:       %.3f %s.",
          clocks_from_ticks(time / runs), clocks_getunit());

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  for (int i = 0; i < 27; ++i) zero_particle_fields(cells[i]);

  const ticks tic = getticks();

  /* Run all the brute-force pairs */
  for (int j = 0; j < 27; ++j)
    if (cells[j] != main_cell)
      pairs_all_bh_stars_density(&runner, main_cell, cells[j]);

  /* And now the self-interaction */
  self_all_bh_stars_density(&runner, main_cell);

  const ticks toc = getticks();

  /* Let's get physical ! */
  end_calculation(main_cell, &cosmo);

  /* Dump */
  sprintf(outputFileName, "bh_stars_brute_force_27_%.150s.dat",
          outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell);

  /* Output timing */
  message("Brute force calculation took : %.3f %s.",
          clocks_from_ticks(toc - tic), clocks_getunit());

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return 0;
}

#else /* BLACK_HOLES_HAVE_STAR_DENSITY */

int main(int argc, char *argv[]) {
  printf(
      "Test not run: the selected black hole model does not support the "
      "star-density loops.\n");
  return 0;
}

#endif /* BLACK_HOLES_HAVE_STAR_DENSITY */
