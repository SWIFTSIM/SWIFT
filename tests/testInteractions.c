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

#ifndef WITH_VECTORIZATION
int main() { return 0; }
#else

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "swift.h"

#define array_align sizeof(float) * VEC_SIZE
#define ACC_THRESHOLD 1e-5

/* Typdef function pointers for serial and vectorised versions of the
 * interaction functions. */
typedef void (*serial_interaction)(float, float *, float, float, struct part *,
                                   struct part *);
typedef void (*vec_interaction)(float *, float *, float *, float *,
                                struct part **, struct part **);

/**
 * @brief Constructs an array of particles in a valid state prior to
 * a IACT_NONSYM and IACT_NONSYM_VEC call.
 *
 * @param count No. of particles to create
 * @param offset The position of the particle offset from (0,0,0).
 * @param spacing Particle spacing.
 * @param h The smoothing length of the particles in units of the inter-particle
 *separation.
 * @param partId The running counter of IDs.
 */
struct part *make_particles(size_t count, double *offset, double spacing,
                            double h, long long *partId) {

  struct part *particles;
  if (posix_memalign((void **)&particles, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(particles, count * sizeof(struct part));

  /* Construct the particles */
  struct part *p;

  /* Set test particle at centre of unit sphere. */
  p = &particles[0];

  /* Place the test particle at the centre of a unit sphere. */
  p->x[0] = 0.0f;
  p->x[1] = 0.0f;
  p->x[2] = 0.0f;

  p->h = h;
  p->id = ++(*partId);
  p->mass = 1.0f;

  /* Place rest of particles around the test particle
   * with random position within a unit sphere. */
  for (size_t i = 1; i < count; ++i) {
    p = &particles[i];

    /* Randomise positions within a unit sphere. */
    p->x[0] = random_uniform(-1.0, 1.0);
    p->x[1] = random_uniform(-1.0, 1.0);
    p->x[2] = random_uniform(-1.0, 1.0);

    /* Randomise velocities. */
    p->v[0] = random_uniform(-0.05, 0.05);
    p->v[1] = random_uniform(-0.05, 0.05);
    p->v[2] = random_uniform(-0.05, 0.05);

    p->h = h;
    p->id = ++(*partId);
    p->mass = 1.0f;
  }
  return particles;
}

/**
 * @brief Populates particle properties needed for the force calculation.
 */
void prepare_force(struct part *parts, size_t count) {

  struct part *p;
  for (size_t i = 0; i < count; ++i) {
    p = &parts[i];
    p->rho = i + 1;
    p->force.balsara = random_uniform(0.0, 1.0);
    p->force.P_over_rho2 = i + 1;
    p->force.soundspeed = random_uniform(2.0, 3.0);
    p->force.v_sig = 0.0f;
    p->force.h_dt = 0.0f;
  }
}

/**
 * @brief Dumps all particle information to a file
 */
void dump_indv_particle_fields(char *fileName, struct part *p) {

  FILE *file = fopen(fileName, "a");

  fprintf(file,
          "%6llu %10f %10f %10f %10f %10f %10f %10e %10e %10e %13e %13e %13e "
          "%13e %13e %13e %13e "
          "%13e %13e %13e %10f\n",
          p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
          p->a_hydro[0], p->a_hydro[1], p->a_hydro[2], p->rho,
          p->density.rho_dh, p->density.wcount, p->density.wcount_dh,
          p->force.h_dt, p->force.v_sig,
#if defined(GADGET2_SPH)
          p->density.div_v, p->density.rot_v[0], p->density.rot_v[1],
          p->density.rot_v[2], p->entropy_dt
#elif defined(DEFAULT_SPH)
          p->density.div_v, p->density.rot_v[0], p->density.rot_v[1],
          p->density.rot_v[2], 0.
#else
          p->density.div_v, p->density.rot_v[0], p->density.rot_v[1],
          p->density.rot_v[2]
#endif
          );
  fclose(file);
}

/**
 * @brief Creates a header for the output file
 */
void write_header(char *fileName) {

  FILE *file = fopen(fileName, "w");
  /* Write header */
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s %10s %10s %10s %10s %13s %13s %13s "
          "%13s %13s %13s %13s"
          "%13s %13s %13s %13s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "a_x", "a_y",
          "a_z", "rho", "rho_dh", "wcount", "wcount_dh", "dh/dt", "v_sig",
          "div_v", "curl_vx", "curl_vy", "curl_vz", "dS/dt");
  fprintf(file, "\n# PARTICLES BEFORE INTERACTION:\n");
  fclose(file);
}

/**
 * @brief Compares the vectorised result against
 * the serial result of the interaction.
 *
 * @param serial_test_part Particle that has been updated serially
 * @param serial_parts Particle array that has been interacted serially
 * @param vec_test_part Particle that has been updated using vectors
 * @param vec_parts Particle array to be interacted using vectors
 * @param count No. of particles that have been interacted
 *
 * @return Non-zero value if difference found, 0 otherwise
 */
int check_results(struct part serial_test_part, struct part *serial_parts,
                  struct part vec_test_part, struct part *vec_parts,
                  int count) {
  int result = 0;
  result += compare_particles(serial_test_part, vec_test_part, ACC_THRESHOLD);

  for (int i = 0; i < count; i++)
    result += compare_particles(serial_parts[i], vec_parts[i], ACC_THRESHOLD);

  return result;
}

/*
 * @brief Calls the serial and vectorised version of an interaction
 * function given by the function pointers.
 *
 * @param test_part Particle that will be updated
 * @param parts Particle array to be interacted
 * @param count No. of particles to be interacted
 * @param serial_inter_func Serial interaction function to be called
 * @param vec_inter_func Vectorised interaction function to be called
 * @param runs No. of times to call interactions
 *
 */
void test_interactions(struct part test_part, struct part *parts, size_t count,
                       serial_interaction serial_inter_func,
                       vec_interaction vec_inter_func, char *filePrefix,
                       size_t runs) {

  ticks serial_time = 0, vec_time = 0;

  FILE *file;
  char serial_filename[200] = "";
  char vec_filename[200] = "";

  strcpy(serial_filename, filePrefix);
  strcpy(vec_filename, filePrefix);
  sprintf(serial_filename + strlen(serial_filename), "_serial.dat");
  sprintf(vec_filename + strlen(vec_filename), "_vec.dat");

  write_header(serial_filename);
  write_header(vec_filename);

  /* Test particle at the center of a unit sphere. */
  struct part pi_serial, pi_vec;

  /* Remaining particles in the sphere that will interact with test particle. */
  struct part pj_serial[count], pj_vec[count];

  /* Stores the separation, smoothing length and pointers to particles
   * needed for the vectorised interaction. */
  float r2q[count] __attribute__((aligned(array_align)));
  float hiq[count] __attribute__((aligned(array_align)));
  float hjq[count] __attribute__((aligned(array_align)));
  float dxq[3 * count] __attribute__((aligned(array_align)));
  struct part *piq[count], *pjq[count];

  /* Call serial interaction a set number of times. */
  for (size_t k = 0; k < runs; k++) {
    /* Reset particle to initial setup */
    pi_serial = test_part;
    for (size_t i = 0; i < count; i++) pj_serial[i] = parts[i];

    /* Only dump data on first run. */
    if (k == 0) {
      /* Dump state of particles before serial interaction. */
      dump_indv_particle_fields(serial_filename, &pi_serial);
      for (size_t i = 0; i < count; i++)
        dump_indv_particle_fields(serial_filename, &pj_serial[i]);
    }

    /* Perform serial interaction */
    for (size_t i = 0; i < count; i++) {
      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (size_t k = 0; k < 3; k++) {
        dx[k] = pi_serial.x[k] - pj_serial[i].x[k];
        r2 += dx[k] * dx[k];
      }

      const ticks tic = getticks();

      serial_inter_func(r2, dx, pi_serial.h, pj_serial[i].h, &pi_serial,
                        &pj_serial[i]);

      serial_time += getticks() - tic;
    }
  }

  file = fopen(serial_filename, "a");
  fprintf(file, "\n# PARTICLES AFTER INTERACTION:\n");
  fclose(file);

  /* Dump result of serial interaction. */
  dump_indv_particle_fields(serial_filename, &pi_serial);
  for (size_t i = 0; i < count; i++)
    dump_indv_particle_fields(serial_filename, &pj_serial[i]);

  /* Call vector interaction a set number of times. */
  for (size_t k = 0; k < runs; k++) {
    /* Reset particle to initial setup */
    pi_vec = test_part;
    for (size_t i = 0; i < count; i++) pj_vec[i] = parts[i];

    /* Setup arrays for vector interaction. */
    for (size_t i = 0; i < count; i++) {
      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (size_t k = 0; k < 3; k++) {
        dx[k] = pi_vec.x[k] - pj_vec[i].x[k];
        r2 += dx[k] * dx[k];
      }

      r2q[i] = r2;
      dxq[3 * i + 0] = dx[0];
      dxq[3 * i + 1] = dx[1];
      dxq[3 * i + 2] = dx[2];
      hiq[i] = pi_vec.h;
      hjq[i] = pj_vec[i].h;
      piq[i] = &pi_vec;
      pjq[i] = &pj_vec[i];
    }

    /* Only dump data on first run. */
    if (k == 0) {
      /* Dump state of particles before vector interaction. */
      dump_indv_particle_fields(vec_filename, piq[0]);
      for (size_t i = 0; i < count; i++)
        dump_indv_particle_fields(vec_filename, pjq[i]);
    }

    const ticks vec_tic = getticks();

    /* Perform vector interaction. */
    for (size_t i = 0; i < count; i += VEC_SIZE) {
      vec_inter_func(&(r2q[i]), &(dxq[3 * i]), &(hiq[i]), &(hjq[i]), &(piq[i]),
                     &(pjq[i]));
    }

    vec_time += getticks() - vec_tic;
  }

  file = fopen(vec_filename, "a");
  fprintf(file, "\n# PARTICLES AFTER INTERACTION:\n");
  fclose(file);

  /* Dump result of vector interaction. */
  dump_indv_particle_fields(vec_filename, piq[0]);
  for (size_t i = 0; i < count; i++)
    dump_indv_particle_fields(vec_filename, pjq[i]);

  /* Check serial results against the vectorised results. */
  if (check_results(pi_serial, pj_serial, pi_vec, pj_vec, count))
    message("Differences found...");

  message("The serial interactions took     : %15lli ticks.",
          serial_time / runs);
  message("The vectorised interactions took : %15lli ticks.", vec_time / runs);
}

/* And go... */
int main(int argc, char *argv[]) {
  size_t runs = 10000;
  double h = 1.0, spacing = 0.5;
  double offset[3] = {0.0, 0.0, 0.0};
  size_t count = 256;

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "h:s:n:r:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 's':
        sscanf(optarg, "%lf", &spacing);
      case 'n':
        sscanf(optarg, "%zu", &count);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || spacing < 0) {
    printf(
        "\nUsage: %s [OPTIONS...]\n"
        "\nGenerates a particle array with equal particle separation."
        "\nThese are then interacted using runner_iact_density and "
        "runner_iact_vec_density."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-s SPACING=0.5     - Spacing between particles"
        "\n-n NUMBER=9        - No. of particles",
        argv[0]);
    exit(1);
  }

  /* Correct count so that VEC_SIZE of particles interact with the test
   * particle. */
  count = count - (count % VEC_SIZE) + 1;

  /* Build the infrastructure */
  static long long partId = 0;
  struct part density_test_particle, force_test_particle;
  struct part *density_particles =
      make_particles(count, offset, spacing, h, &partId);
  struct part *force_particles =
      make_particles(count, offset, spacing, h, &partId);
  prepare_force(force_particles, count);

  /* Define which interactions to call */
  serial_interaction serial_inter_func = &runner_iact_nonsym_density;
  vec_interaction vec_inter_func = &runner_iact_nonsym_vec_density;

  density_test_particle = density_particles[0];
  /* Call the non-sym density test. */
  message("Testing non-symmetrical density interaction...");
  test_interactions(density_test_particle, &density_particles[1], count - 1,
                    serial_inter_func, vec_inter_func, "test_nonsym_density",
                    runs);

  density_particles = make_particles(count, offset, spacing, h, &partId);

  /* Re-assign function pointers. */
  serial_inter_func = &runner_iact_density;
  vec_inter_func = &runner_iact_vec_density;

  density_test_particle = density_particles[0];
  /* Call the symmetrical density test. */
  message("Testing symmetrical density interaction...");
  test_interactions(density_test_particle, &density_particles[1], count - 1,
                    serial_inter_func, vec_inter_func, "test_sym_density",
                    runs);

  /* Re-assign function pointers. */
  serial_inter_func = &runner_iact_nonsym_force;
  vec_inter_func = &runner_iact_nonsym_vec_force;

  force_test_particle = force_particles[0];
  /* Call the test non-sym force test. */
  message("Testing non-symmetrical force interaction...");
  test_interactions(force_test_particle, &force_particles[1], count - 1,
                    serial_inter_func, vec_inter_func, "test_nonsym_force",
                    runs);

  force_particles = make_particles(count, offset, spacing, h, &partId);
  prepare_force(force_particles, count);

  /* Re-assign function pointers. */
  serial_inter_func = &runner_iact_force;
  vec_inter_func = &runner_iact_vec_force;

  force_test_particle = force_particles[0];
  /* Call the test symmetrical force test. */
  message("Testing symmetrical force interaction...");
  test_interactions(force_test_particle, &force_particles[1], count - 1,
                    serial_inter_func, vec_inter_func, "test_sym_force", runs);

  return 0;
}

#endif /* WITH_VECTORIZATION */
