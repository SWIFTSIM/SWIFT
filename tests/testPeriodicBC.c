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
#define DOSELF1 runner_doself1_branch_density
#define DOPAIR1 runner_dopair1_branch_density
#define DOSELF1_NAME "runner_doself1_density_vec"
#define DOPAIR1_NAME "runner_dopair1_density_vec"
#endif

#ifndef DOSELF1
#define DOSELF1 runner_doself1_branch_density
#define DOSELF1_NAME "runner_doself1_density"
#endif

#ifndef DOPAIR1
#define DOPAIR1 runner_dopair1_branch_density
#define DOPAIR1_NAME "runner_dopair1_density"
#endif

#define NODE_ID 0

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
  struct cell *cell = (struct cell *)malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->hydro.parts, count * sizeof(struct part));

  float h_max = 0.f;

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
        h_max = fmax(h_max, part->h);
        part->id = ++(*partId);

#if defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
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
  cell->hydro.h_max = h_max;
  cell->hydro.count = count;
  cell->hydro.dx_max_part = 0.;
  cell->hydro.dx_max_sort = 0.;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->hydro.ti_old_part = 8;
  cell->hydro.ti_end_min = 8;
  cell->hydro.ti_end_max = 8;
  cell->nodeID = NODE_ID;

  shuffle_particles(cell->hydro.parts, cell->hydro.count);

  cell->hydro.sorted = 0;
  for (int k = 0; k < 13; k++) cell->hydro.sort[k] = NULL;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->hydro.parts);
  for (int k = 0; k < 13; k++)
    if (ci->hydro.sort[k] != NULL) free(ci->hydro.sort[k]);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields(struct cell *c) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    hydro_init_part(&c->hydro.parts[pid], NULL);
  }
}

/**
 * @brief Ends the loop by adding the appropriate coefficients
 */
void end_calculation(struct cell *c, const struct cosmology *cosmo) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    hydro_end_density(&c->hydro.parts[pid], cosmo);
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *main_cell, int i, int j,
                          int k) {
  FILE *file = fopen(fileName, "a");

  /* Write header */
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s %10s %13s %13s %13s %13s %13s "
          "%13s %13s %13s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho", "rho_dh",
          "wcount", "wcount_dh", "div_v", "curl_vx", "curl_vy", "curl_vz");

  fprintf(file, "# Centre cell at (i,j,k)=(%d, %d, %d) ---------------------\n",
          i, j, k);

  /* Write main cell */
  for (int pid = 0; pid < main_cell->hydro.count; pid++) {
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
            "%13e %13e %13e\n",
            main_cell->hydro.parts[pid].id, main_cell->hydro.parts[pid].x[0],
            main_cell->hydro.parts[pid].x[1], main_cell->hydro.parts[pid].x[2],
            main_cell->hydro.parts[pid].v[0], main_cell->hydro.parts[pid].v[1],
            main_cell->hydro.parts[pid].v[2],
            hydro_get_comoving_density(&main_cell->hydro.parts[pid]),
#if defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
            0.f,
#else
            main_cell->hydro.parts[pid].density.rho_dh,
#endif
            main_cell->hydro.parts[pid].density.wcount,
            main_cell->hydro.parts[pid].density.wcount_dh,
#if defined(GADGET2_SPH) || defined(DEFAULT_SPH) || defined(HOPKINS_PE_SPH)
            main_cell->hydro.parts[pid].density.div_v,
            main_cell->hydro.parts[pid].density.rot_v[0],
            main_cell->hydro.parts[pid].density.rot_v[1],
            main_cell->hydro.parts[pid].density.rot_v[2]
#else
            0., 0., 0., 0.
#endif
    );
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
    result += compare_particles(&serial_parts[i], &vec_parts[i], threshold);

  return result;
}

/* Just a forward declaration... */
void runner_doself1_density(struct runner *r, struct cell *ci);
void runner_doself1_density_vec(struct runner *r, struct cell *ci);
void runner_dopair1_branch_density(struct runner *r, struct cell *ci,
                                   struct cell *cj);
void runner_doself1_branch_density(struct runner *r, struct cell *c);

void test_boundary_conditions(struct cell **cells, struct runner runner,
                              const int loc_i, const int loc_j, const int loc_k,
                              const int dim, char *swiftOutputFileName,
                              char *bruteForceOutputFileName) {

  /* Store the main cell for future use */
  struct cell *main_cell = cells[loc_i * (dim * dim) + loc_j * dim + loc_k];

  /* Zero the fields */
  for (int j = 0; j < dim * dim * dim; ++j) zero_particle_fields(cells[j]);

/* Run all the pairs */
#ifdef WITH_VECTORIZATION
  runner.ci_cache.count = 0;
  cache_init(&runner.ci_cache, 512);
  runner.cj_cache.count = 0;
  cache_init(&runner.cj_cache, 512);
#endif

  /* Now loop over all the neighbours of this cell
   * and perform the pair interactions. */
  for (int ii = -1; ii < 2; ii++) {
    int iii = loc_i + ii;
    iii = (iii + dim) % dim;
    for (int jj = -1; jj < 2; jj++) {
      int jjj = loc_j + jj;
      jjj = (jjj + dim) % dim;
      for (int kk = -1; kk < 2; kk++) {
        int kkk = loc_k + kk;
        kkk = (kkk + dim) % dim;

        /* Get the neighbouring cell */
        struct cell *cj = cells[iii * (dim * dim) + jjj * dim + kkk];

        if (cj != main_cell) DOPAIR1(&runner, main_cell, cj);
      }
    }
  }

  /* And now the self-interaction */

  DOSELF1(&runner, main_cell);

  /* Let's get physical ! */
  end_calculation(main_cell, runner.e->cosmology);

  /* Dump particles from the main cell. */
  dump_particle_fields(swiftOutputFileName, main_cell, loc_i, loc_j, loc_k);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  for (int i = 0; i < dim * dim * dim; ++i) zero_particle_fields(cells[i]);

  /* Now loop over all the neighbours of this cell
   * and perform the pair interactions. */
  for (int ii = -1; ii < 2; ii++) {
    int iii = loc_i + ii;
    iii = (iii + dim) % dim;
    for (int jj = -1; jj < 2; jj++) {
      int jjj = loc_j + jj;
      jjj = (jjj + dim) % dim;
      for (int kk = -1; kk < 2; kk++) {
        int kkk = loc_k + kk;
        kkk = (kkk + dim) % dim;

        /* Get the neighbouring cell */
        struct cell *cj = cells[iii * (dim * dim) + jjj * dim + kkk];

        if (cj != main_cell) pairs_all_density(&runner, main_cell, cj);
      }
    }
  }

  /* And now the self-interaction */
  self_all_density(&runner, main_cell);

  /* Let's get physical ! */
  end_calculation(main_cell, runner.e->cosmology);

  /* Dump */
  dump_particle_fields(bruteForceOutputFileName, main_cell, loc_i, loc_j,
                       loc_k);
}

/* And go... */
int main(int argc, char *argv[]) {

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  size_t runs = 0, particles = 0;
  double h = 1.23485, size = 1., rho = 1.;
  double perturbation = 0.;
  double threshold = ACC_THRESHOLD;
  char outputFileNameExtension[100] = "";
  char swiftOutputFileName[200] = "";
  char bruteForceOutputFileName[200] = "";
  enum velocity_types vel = velocity_zero;

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
  const int dim = 8;
  struct space space;
  space.periodic = 1;
  space.dim[0] = dim;
  space.dim[1] = dim;
  space.dim[2] = dim;

  struct hydro_props hp;
  hp.h_max = FLT_MAX;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;
  engine.nodeID = NODE_ID;

  struct runner runner;
  runner.e = &engine;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  /* Construct some cells */
  struct cell *cells[dim * dim * dim];
  static long long partId = 0;
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      for (int k = 0; k < dim; ++k) {
        double offset[3] = {i * size, j * size, k * size};
        cells[i * (dim * dim) + j * dim + k] = make_cell(
            particles, offset, size, h, rho, &partId, perturbation, vel);

        runner_do_drift_part(&runner, cells[i * (dim * dim) + j * dim + k], 0);

        runner_do_hydro_sort(&runner, cells[i * (dim * dim) + j * dim + k],
                             0x1FFF, 0, 0);
      }
    }
  }

  /* Create output file names. */
  sprintf(swiftOutputFileName, "swift_periodic_BC_%.150s.dat",
          outputFileNameExtension);
  sprintf(bruteForceOutputFileName, "brute_force_periodic_BC_%.150s.dat",
          outputFileNameExtension);

  /* Delete files if they already exist. */
  remove(swiftOutputFileName);
  remove(bruteForceOutputFileName);

  const int half_dim = (dim - 1) / 2;

  /* Test the periodic boundary conditions for each of the 8 corners. Interact
   * each corner with all of its 26 neighbours.*/
  test_boundary_conditions(cells, runner, 0, 0, 0, dim, swiftOutputFileName,
                           bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, 0, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, 0, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, 0, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, dim - 1, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, dim - 1, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, dim - 1, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, dim - 1, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);

  /* Test the boundary conditions for cells at the centre of each face of the
   * box. */
  test_boundary_conditions(cells, runner, half_dim, half_dim, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, half_dim, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, half_dim, half_dim, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, half_dim, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, half_dim, 0, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, half_dim, dim - 1, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);

  /* Test the boundary conditions for cells at the centre of each edge of the
   * box. */
  test_boundary_conditions(cells, runner, half_dim, dim - 1, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, dim - 1, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, half_dim, dim - 1, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, dim - 1, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);

  test_boundary_conditions(cells, runner, 0, half_dim, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, half_dim, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, half_dim, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, half_dim, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);

  test_boundary_conditions(cells, runner, half_dim, 0, 0, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, dim - 1, 0, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, half_dim, 0, dim - 1, dim,
                           swiftOutputFileName, bruteForceOutputFileName);
  test_boundary_conditions(cells, runner, 0, 0, half_dim, dim,
                           swiftOutputFileName, bruteForceOutputFileName);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < dim * dim * dim; ++i) clean_up(cells[i]);

  return 0;
}
