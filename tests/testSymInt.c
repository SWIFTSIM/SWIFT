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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "swift.h"

char *serial_filename = "test_sym_serial.dat";
char *vec_filename = "test_sym_vec.dat";

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
 * @param vel The type of velocity field (0, random, divergent, rotating)
 */
struct part *make_particles(int count, double *offset, double spacing, double h, 
                            long long *partId) {
 
  struct part *particles;
  if (posix_memalign((void **)&particles, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
   
  /* Construct the particles */
  struct part *p;
  for (size_t i = 0; i < VEC_SIZE + 1; ++i) {
    p = &particles[i];
    p->x[0] = offset[0] + spacing * i;
    p->x[1] = offset[1] + spacing * i;
    p->x[2] = offset[2] + spacing * i;
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
 * @brief Dumps all particle information to a file
 */
void dump_indv_particle_fields(char *fileName, struct part *p) {

  FILE *file = fopen(fileName, "a");
  
  fprintf(
      file,
      "%6llu %10f %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
      "%13e %13e %13e\n",
      p->id, p->x[0], p->x[1],
      p->x[2], p->v[0], p->v[1],
      p->v[2], p->rho, p->rho_dh,
      p->density.wcount, p->density.wcount_dh,
#if defined(GADGET2_SPH)
      p->div_v, p->density.rot_v[0],
      p->density.rot_v[1], p->density.rot_v[2]
#elif defined(DEFAULT_SPH)
      p->density.div_v, p->density.curl_v[0],
      p->density.curl_v[1], p->density.curl_v[2]
#else
      0., 0., 0., 0.
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
            "# %4s %10s %10s %10s %10s %10s %10s %13s %13s %13s %13s %13s "
            "%13s %13s %13s\n",
            "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho", "rho_dh",
            "wcount", "wcount_dh", "div_v", "curl_vx", "curl_vy", "curl_vz");
    fprintf(file,"\nPARTICLES BEFORE INTERACTION:\n");
    fclose(file);
}

/*
 * @brief Calls the serial and vectorised version of the non-symmetrical density interaction.
 *
 * @param parts Particle array to be interacted
 * @param count No. of particles to be interacted
 * 
 */
void test_sym_density_interaction(struct part *parts, int count) {
  
  /* Use the first particle in the array as the one that gets updated. */
  struct part pi = parts[0];
  //const float hig2 = hi * hi * kernel_gamma2;

  FILE *file;
  write_header(serial_filename);
  write_header(vec_filename); 

  /* Dump state of particles before serial interaction. */
  dump_indv_particle_fields(serial_filename,&pi); 
  for(size_t i = 1; i < count; i++)
    dump_indv_particle_fields(serial_filename,&parts[i]); 

  /* Make copy of pi to be used in vectorised version. */
  struct part pi_vec = pi;
  struct part pj_vec[VEC_SIZE];
  for(int i=0; i<VEC_SIZE; i++)
    pj_vec[i] = parts[i + 1];
  
  float r2q[VEC_SIZE] __attribute__((aligned(sizeof(float) * VEC_SIZE)));
  float hiq[VEC_SIZE] __attribute__((aligned(sizeof(float) * VEC_SIZE)));
  float hjq[VEC_SIZE] __attribute__((aligned(sizeof(float) * VEC_SIZE)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(sizeof(float) * VEC_SIZE)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
  
  /* Perform serial interaction */
  for(int i=1; i<count; i++) {
    /* Compute the pairwise distance. */
    float r2 = 0.0f;
    float dx[3];
    for (int k = 0; k < 3; k++) {
      dx[k] = pi.x[k] - parts[i].x[k];
      r2 += dx[k] * dx[k];
    }

    runner_iact_density(r2, dx, pi.h, parts[i].h, &pi, &parts[i]);
  }    

  file = fopen(serial_filename, "a");
  fprintf(file,"\nPARTICLES AFTER INTERACTION:\n");
  fclose(file);

  /* Dump result of serial interaction. */
  dump_indv_particle_fields(serial_filename,&pi); 
  for(size_t i = 1; i < count; i++)
    dump_indv_particle_fields(serial_filename,&parts[i]); 

  /* Setup arrays for vector interaction. */
  for(int i=0; i<VEC_SIZE; i++) {
    /* Compute the pairwise distance. */
    float r2 = 0.0f;
    float dx[3];
    for (int k = 0; k < 3; k++) {
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

  /* Dump state of particles before vector interaction. */
  dump_indv_particle_fields(vec_filename,piq[0]); 
  for(size_t i = 0; i < VEC_SIZE; i++)
    dump_indv_particle_fields(vec_filename,pjq[i]);
  
  /* Perform vector interaction. */
  runner_iact_vec_density(r2q, dxq, hiq, hjq, piq, pjq);
  
  file = fopen(vec_filename, "a");
  fprintf(file,"\nPARTICLES AFTER INTERACTION:\n");
  fclose(file);
  
  /* Dump result of serial interaction. */
  dump_indv_particle_fields(vec_filename,piq[0]); 
  for(size_t i = 0; i < VEC_SIZE; i++)
    dump_indv_particle_fields(vec_filename,pjq[i]);
}

/* And go... */
int main(int argc, char *argv[]) {
  double h = 1.2348, spacing = 0.5;
  double offset[3] = {0.0,0.0,0.0};
  int count = VEC_SIZE + 1;

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "s:h:v:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 's':
        sscanf(optarg, "%lf", &spacing);
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
        "\nThese are then interacted using runner_iact_density and runner_iact_vec_density."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-s spacing         - Spacing between particles"
        "\n-v type (0,1,2,3)  - Velocity field: (zero, random, divergent, "
        "rotating)",
        argv[0]);
    exit(1);
  }
  
  /* Build the infrastructure */
  static long long partId = 0;
  struct part *particles = make_particles(count,offset,spacing,h,&partId);

  /* Call the test non-sym density test. */
  test_sym_density_interaction(particles,count);
  
  return 0;
}
