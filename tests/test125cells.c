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

#if defined(WITH_VECTORIZATION)
#define DOSELF2_NAME "runner_doself2_force_vec"
#define DOPAIR2_NAME "runner_dopair2_force_vec"
#endif

#ifndef DOSELF2_NAME
#define DOSELF2_NAME "runner_doself2_density"
#define DOPAIR2_NAME "runner_dopair2_force"
#endif

#define NODE_ID 0

enum velocity_field {
  velocity_zero,
  velocity_const,
  velocity_divergent,
  velocity_rotating
};

enum pressure_field { pressure_const, pressure_gradient, pressure_divergent };

void set_velocity(struct part *part, enum velocity_field vel, float size) {

  switch (vel) {
    case velocity_zero:
      part->v[0] = 0.f;
      part->v[1] = 0.f;
      part->v[2] = 0.f;
      break;
    case velocity_const:
      part->v[0] = 1.f;
      part->v[1] = 0.f;
      part->v[2] = 0.f;
      break;
    case velocity_divergent:
      part->v[0] = part->x[0] - 2.5 * size;
      part->v[1] = part->x[1] - 2.5 * size;
      part->v[2] = part->x[2] - 2.5 * size;
      break;
    case velocity_rotating:
      part->v[0] = part->x[1];
      part->v[1] = -part->x[0];
      part->v[2] = 0.f;
      break;
  }
}

float get_pressure(double x[3], enum pressure_field press, float size) {

  float r2 = 0.;
  float dx[3] = {0.f};

  switch (press) {
    case pressure_const:
      return 1.5f;
      break;
    case pressure_gradient:
      return 1.5f * x[0]; /* gradient along x */
      break;
    case pressure_divergent:
      dx[0] = x[0] - 2.5 * size;
      dx[1] = x[1] - 2.5 * size;
      dx[2] = x[2] - 2.5 * size;
      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      return sqrt(r2) + 1.5f;
      break;
  }
  return 0.f;
}

void set_energy_state(struct part *part, enum pressure_field press, float size,
                      float density) {

  const float pressure = get_pressure(part->x, press, size);

#if defined(GADGET2_SPH)
  part->entropy = pressure / pow_gamma(density);
#elif defined(HOPKINS_PE_SPH)
  part->entropy = pressure / pow_gamma(density);
#elif defined(DEFAULT_SPH)
  part->u = pressure / (hydro_gamma_minus_one * density);
#elif defined(MINIMAL_SPH) || defined(HOPKINS_PU_SPH) || \
    defined(HOPKINS_PU_SPH_MONAGHAN) || defined(ANARCHY_PU_SPH)
  part->u = pressure / (hydro_gamma_minus_one * density);
#elif defined(PLANETARY_SPH)
  part->u = pressure / (hydro_gamma_minus_one * density);
#elif defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
  part->primitives.P = pressure;
#else
  error("Need to define pressure here !");
#endif
}

struct solution_part {

  long long id;
  double x[3];
  float v[3];
  float a_hydro[3];
  float h;
  float rho;
  float div_v;
  float S;
  float u;
  float P;
  float c;
  float h_dt;
  float v_sig;
  float S_dt;
  float u_dt;
};

void get_solution(const struct cell *main_cell, struct solution_part *solution,
                  float density, enum velocity_field vel,
                  enum pressure_field press, float size) {

  for (int i = 0; i < main_cell->hydro.count; ++i) {

    solution[i].id = main_cell->hydro.parts[i].id;

    solution[i].x[0] = main_cell->hydro.parts[i].x[0];
    solution[i].x[1] = main_cell->hydro.parts[i].x[1];
    solution[i].x[2] = main_cell->hydro.parts[i].x[2];

    solution[i].v[0] = main_cell->hydro.parts[i].v[0];
    solution[i].v[1] = main_cell->hydro.parts[i].v[1];
    solution[i].v[2] = main_cell->hydro.parts[i].v[2];

    solution[i].h = main_cell->hydro.parts[i].h;

    solution[i].rho = density;

    solution[i].P = get_pressure(solution[i].x, press, size);
    solution[i].u = solution[i].P / (solution[i].rho * hydro_gamma_minus_one);
    solution[i].S = solution[i].P / pow_gamma(solution[i].rho);
    solution[i].c = sqrt(hydro_gamma * solution[i].P / solution[i].rho);

    if (vel == velocity_divergent)
      solution[i].div_v = 3.f;
    else
      solution[i].div_v = 0.f;

    solution[i].h_dt = solution[i].h * solution[i].div_v / 3.;

    float gradP[3] = {0.f};
    if (press == pressure_gradient) {
      gradP[0] = 1.5f;
      gradP[1] = 0.f;
      gradP[2] = 0.f;
    } else if (press == pressure_divergent) {
      float dx[3];
      dx[0] = solution[i].x[0] - 2.5 * size;
      dx[1] = solution[i].x[1] - 2.5 * size;
      dx[2] = solution[i].x[2] - 2.5 * size;
      float r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
      if (r > 0.) {
        gradP[0] = dx[0] / r;
        gradP[1] = dx[1] / r;
        gradP[2] = dx[2] / r;
      }
    }

    solution[i].a_hydro[0] = -gradP[0] / solution[i].rho;
    solution[i].a_hydro[1] = -gradP[1] / solution[i].rho;
    solution[i].a_hydro[2] = -gradP[2] / solution[i].rho;

    solution[i].v_sig = 2.f * solution[i].c;

    solution[i].S_dt = 0.f;
    solution[i].u_dt = -(solution[i].P / solution[i].rho) * solution[i].div_v;
  }
}

void reset_particles(struct cell *c, struct hydro_space *hs,
                     enum velocity_field vel, enum pressure_field press,
                     float size, float density) {

  for (int i = 0; i < c->hydro.count; ++i) {

    struct part *p = &c->hydro.parts[i];

    set_velocity(p, vel, size);
    set_energy_state(p, press, size, density);

    hydro_init_part(p, hs);

#if defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
    float volume = p->conserved.mass / density;
#if defined(GIZMO_MFV_SPH)
    p->geometry.volume = volume;
#else
    p->cell.volume = volume;
#endif
    p->primitives.rho = density;
    p->primitives.v[0] = p->v[0];
    p->primitives.v[1] = p->v[1];
    p->primitives.v[2] = p->v[2];
    p->conserved.momentum[0] = p->conserved.mass * p->v[0];
    p->conserved.momentum[1] = p->conserved.mass * p->v[1];
    p->conserved.momentum[2] = p->conserved.mass * p->v[2];
    p->conserved.energy =
        p->primitives.P / hydro_gamma_minus_one * volume +
        0.5f *
            (p->conserved.momentum[0] * p->conserved.momentum[0] +
             p->conserved.momentum[1] * p->conserved.momentum[1] +
             p->conserved.momentum[2] * p->conserved.momentum[2]) /
            p->conserved.mass;
#endif
  }
}

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
 *of the inter-particle separation.
 * @param vel The type of velocity field.
 * @param press The type of pressure field.
 */
struct cell *make_cell(size_t n, const double offset[3], double size, double h,
                       double density, long long *partId, double pert,
                       enum velocity_field vel, enum pressure_field press) {

  const size_t count = n * n * n;
  const double volume = size * size * size;
  struct cell *cell = (struct cell *)malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0)
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  if (posix_memalign((void **)&cell->hydro.xparts, xpart_align,
                     count * sizeof(struct xpart)) != 0)
    error("couldn't allocate particles, no. of x-particles: %d", (int)count);
  bzero(cell->hydro.parts, count * sizeof(struct part));
  bzero(cell->hydro.xparts, count * sizeof(struct xpart));

  float h_max = 0.f;

  /* Construct the parts */
  struct part *part = cell->hydro.parts;
  struct xpart *xpart = cell->hydro.xparts;
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
        part->h = size * h / (float)n;
        h_max = fmax(h_max, part->h);

#if defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
        part->conserved.mass = density * volume / count;
#else
        part->mass = density * volume / count;
#endif

        set_velocity(part, vel, size);
        set_energy_state(part, press, size, density);

        hydro_first_init_part(part, xpart);

        part->id = ++(*partId);
        part->time_bin = 1;

#if defined(GIZMO_MFV_SPH)
        part->geometry.volume = part->conserved.mass / density;
        part->primitives.rho = density;
        part->primitives.v[0] = part->v[0];
        part->primitives.v[1] = part->v[1];
        part->primitives.v[2] = part->v[2];
        part->conserved.momentum[0] = part->conserved.mass * part->v[0];
        part->conserved.momentum[1] = part->conserved.mass * part->v[1];
        part->conserved.momentum[2] = part->conserved.mass * part->v[2];
        part->conserved.energy =
            part->primitives.P / hydro_gamma_minus_one * volume +
            0.5f *
                (part->conserved.momentum[0] * part->conserved.momentum[0] +
                 part->conserved.momentum[1] * part->conserved.momentum[1] +
                 part->conserved.momentum[2] * part->conserved.momentum[2]) /
                part->conserved.mass;
#endif

#ifdef SWIFT_DEBUG_CHECKS
        part->ti_drift = 8;
        part->ti_kick = 8;
#endif

        ++part;
        ++xpart;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->hydro.h_max = h_max;
  cell->hydro.count = count;
  cell->grav.count = 0;
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

  // shuffle_particles(cell->hydro.parts, cell->hydro.count);

  cell->hydro.sorted = 0;
  for (int k = 0; k < 13; k++) cell->hydro.sort[k] = NULL;

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->hydro.parts);
  free(ci->hydro.xparts);
  for (int k = 0; k < 13; k++)
    if (ci->hydro.sort[k] != NULL) free(ci->hydro.sort[k]);
  free(ci);
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *main_cell,
                          struct solution_part *solution, int with_solution) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
  fprintf(file,
          "# %4s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %13s %13s "
          "%13s %13s %13s %8s %8s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "h", "rho",
          "div_v", "S", "u", "P", "c", "a_x", "a_y", "a_z", "h_dt", "v_sig",
          "dS/dt", "du/dt");

  fprintf(file, "# Main cell --------------------------------------------\n");

  /* Write main cell */
  for (int pid = 0; pid < main_cell->hydro.count; pid++) {
    fprintf(file,
            "%6llu %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f "
            "%8.5f "
            "%8.5f %8.5f %13e %13e %13e %13e %13e %8.5f %8.5f\n",
            main_cell->hydro.parts[pid].id, main_cell->hydro.parts[pid].x[0],
            main_cell->hydro.parts[pid].x[1], main_cell->hydro.parts[pid].x[2],
            main_cell->hydro.parts[pid].v[0], main_cell->hydro.parts[pid].v[1],
            main_cell->hydro.parts[pid].v[2], main_cell->hydro.parts[pid].h,
            hydro_get_comoving_density(&main_cell->hydro.parts[pid]),
#if defined(MINIMAL_SPH) || defined(PLANETARY_SPH) ||              \
    defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH) ||            \
    defined(HOPKINS_PU_SPH) || defined(HOPKINS_PU_SPH_MONAGHAN) || \
    defined(ANARCHY_PU_SPH)
            0.f,
#elif defined(ANARCHY_PU_SPH)
            main_cell->hydro.parts[pid].viscosity.div_v,
#else
            main_cell->hydro.parts[pid].density.div_v,
#endif
            hydro_get_drifted_comoving_entropy(&main_cell->hydro.parts[pid]),
            hydro_get_drifted_comoving_internal_energy(
                &main_cell->hydro.parts[pid]),
            hydro_get_comoving_pressure(&main_cell->hydro.parts[pid]),
            hydro_get_comoving_soundspeed(&main_cell->hydro.parts[pid]),
            main_cell->hydro.parts[pid].a_hydro[0],
            main_cell->hydro.parts[pid].a_hydro[1],
            main_cell->hydro.parts[pid].a_hydro[2],
            main_cell->hydro.parts[pid].force.h_dt,
#if defined(GADGET2_SPH)
            main_cell->hydro.parts[pid].force.v_sig,
            main_cell->hydro.parts[pid].entropy_dt, 0.f
#elif defined(DEFAULT_SPH)
            main_cell->hydro.parts[pid].force.v_sig, 0.f,
            main_cell->hydro.parts[pid].force.u_dt
#elif defined(MINIMAL_SPH) || defined(HOPKINS_PU_SPH) || \
    defined(HOPKINS_PU_SPH_MONAGHAN)
            main_cell->hydro.parts[pid].force.v_sig, 0.f,
            main_cell->hydro.parts[pid].u_dt
#elif defined(ANARCHY_PU_SPH)
            main_cell->hydro.parts[pid].viscosity.v_sig, 0.f,
            main_cell->hydro.parts[pid].u_dt
#else
            0.f, 0.f, 0.f
#endif
    );
  }

  if (with_solution) {

    fprintf(file, "# Solution ---------------------------------------------\n");

    for (int pid = 0; pid < main_cell->hydro.count; pid++) {
      fprintf(file,
              "%6llu %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f "
              "%8.5f %8.5f "
              "%8.5f %8.5f %13f %13f %13f %13f %13f %8.5f %8.5f\n",
              solution[pid].id, solution[pid].x[0], solution[pid].x[1],
              solution[pid].x[2], solution[pid].v[0], solution[pid].v[1],
              solution[pid].v[2], solution[pid].h, solution[pid].rho,
              solution[pid].div_v, solution[pid].S, solution[pid].u,
              solution[pid].P, solution[pid].c, solution[pid].a_hydro[0],
              solution[pid].a_hydro[1], solution[pid].a_hydro[2],
              solution[pid].h_dt, solution[pid].v_sig, solution[pid].S_dt,
              solution[pid].u_dt);
    }
  }

  fclose(file);
}

/* Just a forward declaration... */
void runner_dopair1_branch_density(struct runner *r, struct cell *ci,
                                   struct cell *cj);
void runner_doself1_branch_density(struct runner *r, struct cell *ci);
#ifdef EXTRA_HYDRO_LOOP
void runner_dopair1_branch_gradient(struct runner *r, struct cell *ci,
                                    struct cell *cj);
void runner_doself1_branch_gradient(struct runner *r, struct cell *ci);
#endif /* EXTRA_HYDRO LOOP */
void runner_dopair2_branch_force(struct runner *r, struct cell *ci,
                                 struct cell *cj);
void runner_doself2_branch_force(struct runner *r, struct cell *ci);
void runner_doself2_force(struct runner *r, struct cell *ci);
void runner_doself2_force_vec(struct runner *r, struct cell *ci);

/* And go... */
int main(int argc, char *argv[]) {

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  size_t runs = 0, particles = 0;
  double h = 1.23485, size = 1., rho = 2.5;
  double perturbation = 0.;
  char outputFileNameExtension[100] = "";
  char outputFileName[200] = "";
  enum velocity_field vel = velocity_zero;
  enum pressure_field press = pressure_const;

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
  while ((c = getopt(argc, argv, "m:s:h:n:r:t:d:f:v:p:")) != -1) {
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
      case 'p':
        sscanf(optarg, "%d", (int *)&press);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || particles == 0 || runs == 0) {
    printf(
        "\nUsage: %s -n PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates 125 cells, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density() and "
        "runner_doself1_density() followed by runner_dopair2_force() and "
        "runner_doself2_force()"
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-m rho             - Physical density in the cell"
        "\n-s size            - Physical size of the cell"
        "\n-d pert            - Perturbation to apply to the particles [0,1["
        "\n-v type (0,1,2,3)  - Velocity field: (zero, constant, divergent, "
        "rotating)"
        "\n-p type (0,1,2)    - Pressure field: (constant, gradient divergent)"
        "\n-f fileName        - Part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  /* Help users... */
  message("DOSELF2 function called: %s", DOSELF2_NAME);
  message("DOPAIR2 function called: %s", DOPAIR2_NAME);
  message("Adiabatic index: ga = %f", hydro_gamma);
  message("Hydro implementation: %s", SPH_IMPLEMENTATION);
  message("Smoothing length: h = %f", h * size);
  message("Kernel:               %s", kernel_name);
  message("Neighbour target: N = %f", pow_dimension(h) * kernel_norm);
  message("Density target: rho = %f", rho);
  message("div_v target:   div = %f", vel == 2 ? 3.f : 0.f);
  message("curl_v target: curl = [0., 0., %f]", vel == 3 ? -2.f : 0.f);
  if (press == pressure_const)
    message("P field constant");
  else if (press == pressure_gradient)
    message("P field gradient");
  else
    message("P field divergent");

  printf("\n");

#if !defined(HYDRO_DIMENSION_3D)
  message("test125cells only useful in 3D. Change parameters in const.h !");
  return 1;
#endif

  /* Build the infrastructure */
  struct space space;
  space.periodic = 1;
  space.dim[0] = 5.;
  space.dim[1] = 5.;
  space.dim[2] = 5.;
  hydro_space_init(&space.hs, &space);

  struct phys_const prog_const;
  prog_const.const_newton_G = 1.f;

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  hp.eta_neighbours = h;
  hp.h_tolerance = 1e0;
  hp.h_max = FLT_MAX;
  hp.h_min = 0.f;
  hp.h_min_ratio = 0.f;
  hp.max_smoothing_iterations = 10;
  hp.CFL_condition = 0.1;

  struct engine engine;
  bzero(&engine, sizeof(struct engine));
  engine.hydro_properties = &hp;
  engine.physical_constants = &prog_const;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.nodeID = NODE_ID;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

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
        cells[i * 25 + j * 5 + k] = make_cell(
            particles, offset, size, h, rho, &partId, perturbation, vel, press);

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

  /* Construct the real solution */
  struct solution_part *solution = (struct solution_part *)malloc(
      main_cell->hydro.count * sizeof(struct solution_part));
  get_solution(main_cell, solution, rho, vel, press, size);

  ticks timings[27];
  for (int i = 0; i < 27; i++) timings[i] = 0;

  /* Start the test */
  ticks time = 0;
  for (size_t n = 0; n < runs; ++n) {

    const ticks tic = getticks();

    /* Initialise the particles */
    for (int j = 0; j < 125; ++j) runner_do_drift_part(&runner, cells[j], 0);

    /* Reset particles. */
    for (int i = 0; i < 125; ++i) {
      for (int pid = 0; pid < cells[i]->hydro.count; ++pid)
        hydro_init_part(&cells[i]->hydro.parts[pid], &space.hs);
    }

    /* First, sort stuff */
    for (int j = 0; j < 125; ++j)
      runner_do_hydro_sort(&runner, cells[j], 0x1FFF, 0, 0);

      /* Do the density calculation */

/* Initialise the particle cache. */
#ifdef WITH_VECTORIZATION
    runner.ci_cache.count = 0;
    runner.cj_cache.count = 0;
    cache_init(&runner.ci_cache, 512);
    cache_init(&runner.cj_cache, 512);
#endif

    /* Run all the  (only once !)*/
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

                if (cj > ci) runner_dopair1_branch_density(&runner, ci, cj);
              }
            }
          }
        }
      }
    }

    /* And now the self-interaction for the central cells*/
    for (int j = 0; j < 27; ++j)
      runner_doself1_branch_density(&runner, inner_cells[j]);

    /* Ghost to finish everything on the central cells */
    for (int j = 0; j < 27; ++j) runner_do_ghost(&runner, inner_cells[j], 0);

#ifdef EXTRA_HYDRO_LOOP
    /* We need to do the gradient loop and the extra ghost! */
    message(
        "Extra hydro loop detected, running gradient loop in test125cells.");

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

                if (cj > ci) runner_dopair1_branch_gradient(&runner, ci, cj);
              }
            }
          }
        }
      }
    }

    /* And now the self-interaction for the central cells */
    for (int j = 0; j < 27; ++j)
      runner_doself1_branch_gradient(&runner, inner_cells[j]);

    /* Extra ghost to finish everything on the central cells */
    for (int j = 0; j < 27; ++j)
      runner_do_extra_ghost(&runner, inner_cells[j], 0);

#endif /* EXTRA_HYDRO_LOOP */

      /* Do the force calculation */

#ifdef WITH_VECTORIZATION
    /* Initialise the cache. */
    cache_clean(&runner.ci_cache);
    cache_clean(&runner.cj_cache);
    cache_init(&runner.ci_cache, 512);
    cache_init(&runner.cj_cache, 512);
#endif

    int ctr = 0;
    /* Do the pairs (for the central 27 cells) */
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        for (int k = 1; k < 4; k++) {

          struct cell *cj = cells[i * 25 + j * 5 + k];

          if (main_cell != cj) {

            const ticks sub_tic = getticks();

            runner_dopair2_branch_force(&runner, main_cell, cj);

            timings[ctr++] += getticks() - sub_tic;
          }
        }
      }
    }

    ticks self_tic = getticks();

    /* And now the self-interaction for the main cell */
    runner_doself2_branch_force(&runner, main_cell);

    timings[26] += getticks() - self_tic;

    /* Finally, give a gentle kick */
    runner_do_end_force(&runner, main_cell, 0);
    const ticks toc = getticks();
    time += toc - tic;

    /* Dump if necessary */
    if (n == 0) {
      sprintf(outputFileName, "swift_dopair_125_%.150s.dat",
              outputFileNameExtension);
      dump_particle_fields(outputFileName, main_cell, solution, 0);
    }

    for (int i = 0; i < 125; ++i) {
      for (int pid = 0; pid < cells[i]->hydro.count; ++pid)
        hydro_init_part(&cells[i]->hydro.parts[pid], &space.hs);
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

  message("Corner calculations took:     %15lli ticks.", corner_time / runs);
  message("Edge calculations took:       %15lli ticks.", edge_time / runs);
  message("Face calculations took:       %15lli ticks.", face_time / runs);
  message("Self calculations took:       %15lli ticks.", self_time / runs);
  message("SWIFT calculation took:       %15lli ticks.", time / runs);

  for (int j = 0; j < 125; ++j)
    reset_particles(cells[j], &space.hs, vel, press, size, rho);

  /* NOW BRUTE-FORCE CALCULATION */

  const ticks tic = getticks();

  /* Kick the central cell */
  // runner_do_kick1(&runner, main_cell, 0);

  /* And drift it */
  // runner_do_drift_particles(&runner, main_cell, 0);

  /* Initialise the particles */
  // for (int j = 0; j < 125; ++j) runner_do_drift_particles(&runner, cells[j],
  // 0);

  /* Do the density calculation */

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

              if (cj > ci) pairs_all_density(&runner, ci, cj);
            }
          }
        }
      }
    }
  }

  /* And now the self-interaction for the central cells*/
  for (int j = 0; j < 27; ++j) self_all_density(&runner, inner_cells[j]);

  /* Ghost to finish everything on the central cells */
  for (int j = 0; j < 27; ++j) runner_do_ghost(&runner, inner_cells[j], 0);

#ifdef EXTRA_HYDRO_LOOP
  /* We need to do the gradient loop and the extra ghost! */

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

              if (cj > ci) pairs_all_gradient(&runner, ci, cj);
            }
          }
        }
      }
    }
  }

  /* And now the self-interaction for the central cells */
  for (int j = 0; j < 27; ++j) self_all_gradient(&runner, inner_cells[j]);

  /* Extra ghost to finish everything on the central cells */
  for (int j = 0; j < 27; ++j)
    runner_do_extra_ghost(&runner, inner_cells[j], 0);

#endif /* EXTRA_HYDRO_LOOP */

  /* Do the force calculation */

  /* Do the pairs (for the central 27 cells) */
  for (int i = 1; i < 4; i++) {
    for (int j = 1; j < 4; j++) {
      for (int k = 1; k < 4; k++) {

        struct cell *cj = cells[i * 25 + j * 5 + k];

        if (main_cell != cj) pairs_all_force(&runner, main_cell, cj);
      }
    }
  }

  /* And now the self-interaction for the main cell */
  self_all_force(&runner, main_cell);

  /* Finally, give a gentle kick */
  runner_do_end_force(&runner, main_cell, 0);
  // runner_do_kick2(&runner, main_cell, 0);

  const ticks toc = getticks();

  /* Output timing */
  message("Brute force calculation took: %15lli ticks.", toc - tic);

  sprintf(outputFileName, "brute_force_125_%.150s.dat",
          outputFileNameExtension);
  dump_particle_fields(outputFileName, main_cell, solution, 0);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 125; ++i) clean_up(cells[i]);
  free(solution);

#ifdef WITH_VECTORIZATION
  cache_clean(&runner.ci_cache);
  cache_clean(&runner.cj_cache);
#endif

  return 0;
}
