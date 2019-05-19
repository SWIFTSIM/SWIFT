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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

#define NODE_ID 1

/* Typdef function pointer for interaction function. */
typedef void (*interaction_func)(struct runner *, struct cell *, struct cell *);
typedef void (*init_func)(struct cell *, const struct cosmology *,
                          const struct hydro_props *);
typedef void (*finalise_func)(struct cell *, const struct cosmology *);

/**
 * @brief Constructs a cell and all of its particle in a valid state prior to
 * a DOPAIR or DOSELF calcuation.
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
 * @param h_pert The perturbation to apply to the smoothing length.
 * @param fraction_active The fraction of particles that should be active in the
 * cell.
 */
struct cell *make_cell(size_t n, double *offset, double size, double h,
                       double density, long long *partId, double pert,
                       double h_pert, double fraction_active) {
  const size_t count = n * n * n;
  const double volume = size * size * size;
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
        part->v[0] = random_uniform(-0.05, 0.05);
        part->v[1] = random_uniform(-0.05, 0.05);
        part->v[2] = random_uniform(-0.05, 0.05);

        if (h_pert)
          part->h = size * h * random_uniform(1.f, h_pert) / (float)n;
        else
          part->h = size * h / (float)n;
        h_max = fmaxf(h_max, part->h);
        part->id = ++(*partId);

/* Set the mass */
#if defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
        part->conserved.mass = density * volume / count;

#ifdef SHADOWFAX_SPH
        double anchor[3] = {0., 0., 0.};
        double side[3] = {1., 1., 1.};
        voronoi_cell_init(&part->cell, part->x, anchor, side);
#endif /* SHADOWFAX_SPH */

#else
        part->mass = density * volume / count;
#endif /* GIZMO_MFV_SPH */

/* Set the thermodynamic variable */
#if defined(GADGET2_SPH)
        part->entropy = 1.f;
#elif defined(MINIMAL_SPH) || defined(HOPKINS_PU_SPH) ||           \
    defined(HOPKINS_PU_SPH_MONAGHAN) || defined(ANARCHY_PU_SPH) || \
    defined(ANARCHY_DU_SPH)
        part->u = 1.f;
#elif defined(HOPKINS_PE_SPH)
        part->entropy = 1.f;
        part->entropy_one_over_gamma = 1.f;
#endif

        /* Set the time-bin */
        if (random_uniform(0, 1.f) < fraction_active)
          part->time_bin = 1;
        else
          part->time_bin = num_time_bins + 1;

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
  cell->hydro.ti_end_max = 10;
  cell->nodeID = NODE_ID;

  shuffle_particles(cell->hydro.parts, cell->hydro.count);

  cell->hydro.sorted = 0;
  for (int k = 0; k < 13; k++) cell->hydro.sort[k] = NULL;

  return cell;
}

void clean_up(struct cell *ci) {
  cell_free_hydro_sorts(ci);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields_density(struct cell *c, const struct cosmology *cosmo,
                                  const struct hydro_props *hydro_props) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    hydro_init_part(&c->hydro.parts[pid], NULL);
  }
}

/**
 * @brief Initializes all particles field to be ready for a force calculation
 */
void zero_particle_fields_force(struct cell *c, const struct cosmology *cosmo,
                                const struct hydro_props *hydro_props) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    struct part *p = &c->hydro.parts[pid];
    struct xpart *xp = &c->hydro.xparts[pid];

/* Mimic the result of a density calculation */
#ifdef GADGET2_SPH
    p->rho = 1.f;
    p->density.rho_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
    p->density.rot_v[0] = 0.f;
    p->density.rot_v[1] = 0.f;
    p->density.rot_v[2] = 0.f;
    p->density.div_v = 0.f;
#endif /* GADGET-2 */
#if defined(MINIMAL_SPH) || defined(ANARCHY_DU_SPH)
    p->rho = 1.f;
    p->density.rho_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
#endif /* MINIMAL */
#ifdef HOPKINS_PE_SPH
    p->rho = 1.f;
    p->rho_bar = 1.f;
    p->density.rho_dh = 0.f;
    p->density.pressure_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
#endif /* PRESSURE-ENTROPY */
#if defined(HOPKINS_PU_SPH) || defined(HOPKINS_PU_SPH_MONAGHAN) || \
    defined(ANARCHY_PU_SPH)
    p->rho = 1.f;
    p->pressure_bar = 0.6666666;
    p->density.rho_dh = 0.f;
    p->density.pressure_bar_dh = 0.f;
    p->density.wcount = 48.f / (kernel_norm * pow_dimension(p->h));
    p->density.wcount_dh = 0.f;
#endif /* PRESSURE-ENERGY */
#if defined(ANARCHY_PU_SPH) || defined(ANARCHY_DU_SPH)
    /* Initialise viscosity variables */
    p->viscosity.alpha = 0.8;
    p->viscosity.div_v = 0.f;
    p->viscosity.div_v_previous_step = 0.f;
    p->viscosity.v_sig = hydro_get_comoving_soundspeed(p);
#endif /* ANARCHY_PU_SPH viscosity variables */

    /* And prepare for a round of force tasks. */
    hydro_prepare_force(p, xp, cosmo, hydro_props, 0.);
    hydro_reset_acceleration(p);
  }
}

/**
 * @brief Ends the density loop by adding the appropriate coefficients
 */
void end_calculation_density(struct cell *c, const struct cosmology *cosmo) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    hydro_end_density(&c->hydro.parts[pid], cosmo);

    /* Recover the common "Neighbour number" definition */
    c->hydro.parts[pid].density.wcount *= pow_dimension(c->hydro.parts[pid].h);
    c->hydro.parts[pid].density.wcount *= kernel_norm;
  }
}

/**
 * @brief Ends the force loop by adding the appropriate coefficients
 */
void end_calculation_force(struct cell *c, const struct cosmology *cosmo) {
  for (int pid = 0; pid < c->hydro.count; pid++) {
    hydro_end_force(&c->hydro.parts[pid], cosmo);
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *ci, struct cell *cj) {
  FILE *file = fopen(fileName, "a");

  /* Write header */
  fprintf(file, "# %4s %13s %13s\n", "ID", "wcount", "h_dt");

  fprintf(file, "# ci --------------------------------------------\n");

  for (int pid = 0; pid < ci->hydro.count; pid++) {
    fprintf(file, "%6llu %13e %13e\n", ci->hydro.parts[pid].id,
            ci->hydro.parts[pid].density.wcount,
            ci->hydro.parts[pid].force.h_dt);
  }

  fprintf(file, "# cj --------------------------------------------\n");

  for (int pjd = 0; pjd < cj->hydro.count; pjd++) {
    fprintf(file, "%6llu %13e %13e\n", cj->hydro.parts[pjd].id,
            cj->hydro.parts[pjd].density.wcount,
            cj->hydro.parts[pjd].force.h_dt);
  }

  fclose(file);
}

/* Just a forward declaration... */
void runner_dopair1_density(struct runner *r, struct cell *ci, struct cell *cj);
void runner_dopair2_force_vec(struct runner *r, struct cell *ci,
                              struct cell *cj);
void runner_doself1_density_vec(struct runner *r, struct cell *ci);
void runner_dopair1_branch_density(struct runner *r, struct cell *ci,
                                   struct cell *cj);
void runner_dopair2_branch_force(struct runner *r, struct cell *ci,
                                 struct cell *cj);

/**
 * @brief Computes the pair interactions of two cells using SWIFT and a brute
 * force implementation.
 */
void test_pair_interactions(struct runner *runner, struct cell **ci,
                            struct cell **cj, char *swiftOutputFileName,
                            char *bruteForceOutputFileName,
                            interaction_func serial_interaction,
                            interaction_func vec_interaction, init_func init,
                            finalise_func finalise) {

  runner_do_hydro_sort(runner, *ci, 0x1FFF, 0, 0);
  runner_do_hydro_sort(runner, *cj, 0x1FFF, 0, 0);

  /* Zero the fields */
  init(*ci, runner->e->cosmology, runner->e->hydro_properties);
  init(*cj, runner->e->cosmology, runner->e->hydro_properties);

  /* Run the test */
  vec_interaction(runner, *ci, *cj);

  /* Let's get physical ! */
  finalise(*ci, runner->e->cosmology);
  finalise(*cj, runner->e->cosmology);

  /* Dump if necessary */
  dump_particle_fields(swiftOutputFileName, *ci, *cj);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  init(*ci, runner->e->cosmology, runner->e->hydro_properties);
  init(*cj, runner->e->cosmology, runner->e->hydro_properties);

  /* Run the brute-force test */
  serial_interaction(runner, *ci, *cj);

  /* Let's get physical ! */
  finalise(*ci, runner->e->cosmology);
  finalise(*cj, runner->e->cosmology);

  dump_particle_fields(bruteForceOutputFileName, *ci, *cj);
}

/**
 * @brief Computes the pair interactions of two cells in various configurations.
 */
void test_all_pair_interactions(
    struct runner *runner, double *offset2, size_t particles, double size,
    double h, double rho, long long *partId, double perturbation, double h_pert,
    char *swiftOutputFileName, char *bruteForceOutputFileName,
    interaction_func serial_interaction, interaction_func vec_interaction,
    init_func init, finalise_func finalise) {

  double offset1[3] = {0, 0, 0};
  struct cell *ci, *cj;

  /* Only one particle in each cell. */
  ci = make_cell(1, offset1, size, h, rho, partId, perturbation, h_pert, 1.);
  cj = make_cell(1, offset2, size, h, rho, partId, perturbation, h_pert, 1.);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* All active particles. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 1.);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 1.);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* Half particles are active. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 0.5);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 0.5);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* All particles inactive. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 0.);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 0.);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* 10% of particles active. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 0.1);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 0.1);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* One active cell one inactive cell. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 1.0);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 0.);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* One active cell one inactive cell. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 0.);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 1.0);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* Smaller cells, all active. */
  ci = make_cell(2, offset1, size, h, rho, partId, perturbation, h_pert, 1.0);
  cj = make_cell(2, offset2, size, h, rho, partId, perturbation, h_pert, 1.0);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* Different numbers of particles in each cell. */
  ci = make_cell(10, offset1, size, h, rho, partId, perturbation, h_pert, 0.5);
  cj = make_cell(3, offset2, size, h, rho, partId, perturbation, h_pert, 0.75);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* One cell inactive and the other only half active. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 0.5);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 0.);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  clean_up(ci);
  clean_up(cj);

  /* One cell inactive and the other only half active. */
  ci = make_cell(particles, offset1, size, h, rho, partId, perturbation, h_pert,
                 0.);
  cj = make_cell(particles, offset2, size, h, rho, partId, perturbation, h_pert,
                 0.5);

  test_pair_interactions(runner, &ci, &cj, swiftOutputFileName,
                         bruteForceOutputFileName, serial_interaction,
                         vec_interaction, init, finalise);

  /* Clean things to make the sanitizer happy ... */
  clean_up(ci);
  clean_up(cj);
}

int main(int argc, char *argv[]) {
  size_t particles = 0, runs = 0, type = 0;
  double h = 1.23485, size = 1., rho = 1.;
  double perturbation = 0.1, h_pert = 1.1;
  struct space space;
  struct engine engine;
  struct cosmology cosmo;
  struct hydro_props hydro_props;
  struct runner *runner;
  char c;
  static long long partId = 0;
  char outputFileNameExtension[100] = "";
  char swiftOutputFileName[200] = "";
  char bruteForceOutputFileName[200] = "";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Generate a RNG seed from time. */
  unsigned int seed = time(NULL);

  while ((c = getopt(argc, argv, "h:p:n:r:t:d:s:f:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 'p':
        sscanf(optarg, "%lf", &h_pert);
        break;
      case 'n':
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
      case 's':
        sscanf(optarg, "%u", &seed);
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
        "\nUsage: %s -n PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates a cell pair, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density."
        "\n\nOptions:"
        "\n-t TYPE=0          - cells share face (0), edge (1) or corner (2)"
        "\n-h DISTANCE=1.2348 - smoothing length"
        "\n-p                 - Random fractional change in h, h=h*random(1,p)"
        "\n-d pert            - perturbation to apply to the particles [0,1["
        "\n-s seed            - seed for RNG"
        "\n-f fileName        - part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Seed RNG. */
  message("Seed used for RNG: %d", seed);
  srand(seed);

  space.periodic = 0;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.nodeID = NODE_ID;

  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;
  hydro_props_init_no_hydro(&hydro_props);
  engine.hydro_properties = &hydro_props;

  if (posix_memalign((void **)&runner, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct runner)) != 0) {
    error("couldn't allocate runner");
  }

  runner->e = &engine;

  /* Create output file names. */
  sprintf(swiftOutputFileName, "swift_dopair_%.150s.dat",
          outputFileNameExtension);
  sprintf(bruteForceOutputFileName, "brute_force_pair_%.150s.dat",
          outputFileNameExtension);

  /* Delete files if they already exist. */
  remove(swiftOutputFileName);
  remove(bruteForceOutputFileName);

#ifdef WITH_VECTORIZATION
  runner->ci_cache.count = 0;
  cache_init(&runner->ci_cache, 512);
  runner->cj_cache.count = 0;
  cache_init(&runner->cj_cache, 512);
#endif

  double offset[3] = {1., 0., 0.};

  /* Define which interactions to call */
  interaction_func serial_inter_func = &pairs_all_density;
  interaction_func vec_inter_func = &runner_dopair1_branch_density;
  init_func init = &zero_particle_fields_density;
  finalise_func finalise = &end_calculation_density;

  /* Test a pair of cells face-on. */
  test_all_pair_interactions(runner, offset, particles, size, h, rho, &partId,
                             perturbation, h_pert, swiftOutputFileName,
                             bruteForceOutputFileName, serial_inter_func,
                             vec_inter_func, init, finalise);

  /* Test a pair of cells edge-on. */
  offset[0] = 1.;
  offset[1] = 1.;
  offset[2] = 0.;
  test_all_pair_interactions(runner, offset, particles, size, h, rho, &partId,
                             perturbation, h_pert, swiftOutputFileName,
                             bruteForceOutputFileName, serial_inter_func,
                             vec_inter_func, init, finalise);

  /* Test a pair of cells corner-on. */
  offset[0] = 1.;
  offset[1] = 1.;
  offset[2] = 1.;
  test_all_pair_interactions(runner, offset, particles, size, h, rho, &partId,
                             perturbation, h_pert, swiftOutputFileName,
                             bruteForceOutputFileName, serial_inter_func,
                             vec_inter_func, init, finalise);

  /* Re-assign function pointers. */
  serial_inter_func = &pairs_all_force;
  vec_inter_func = &runner_dopair2_branch_force;
  init = &zero_particle_fields_force;
  finalise = &end_calculation_force;

  /* Create new output file names. */
  sprintf(swiftOutputFileName, "swift_dopair2_force_%.150s.dat",
          outputFileNameExtension);
  sprintf(bruteForceOutputFileName, "brute_force_dopair2_%.150s.dat",
          outputFileNameExtension);

  /* Delete files if they already exist. */
  remove(swiftOutputFileName);
  remove(bruteForceOutputFileName);

  /* Test a pair of cells face-on. */
  offset[0] = 1.;
  offset[1] = 0.;
  offset[2] = 0.;
  test_all_pair_interactions(runner, offset, particles, size, h, rho, &partId,
                             perturbation, h_pert, swiftOutputFileName,
                             bruteForceOutputFileName, serial_inter_func,
                             vec_inter_func, init, finalise);

  /* Test a pair of cells edge-on. */
  offset[0] = 1.;
  offset[1] = 1.;
  offset[2] = 0.;
  test_all_pair_interactions(runner, offset, particles, size, h, rho, &partId,
                             perturbation, h_pert, swiftOutputFileName,
                             bruteForceOutputFileName, serial_inter_func,
                             vec_inter_func, init, finalise);

  /* Test a pair of cells corner-on. */
  offset[0] = 1.;
  offset[1] = 1.;
  offset[2] = 1.;
  test_all_pair_interactions(runner, offset, particles, size, h, rho, &partId,
                             perturbation, h_pert, swiftOutputFileName,
                             bruteForceOutputFileName, serial_inter_func,
                             vec_inter_func, init, finalise);
  return 0;
}
