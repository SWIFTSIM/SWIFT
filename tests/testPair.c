#include <fenv.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "swift.h"

/* n is both particles per axis and box size:
 * particles are generated on a mesh with unit spacing
 */
struct cell *make_cell(size_t n, double *offset, double h,
                       unsigned long long *partId) {
  size_t count = n * n * n;
  struct cell *cell = malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));
  struct part *part;
  size_t x, y, z, size;

  size = count * sizeof(struct part);
  if (posix_memalign((void **)&cell->parts, part_align, size) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->parts, count * sizeof(struct part));

  part = cell->parts;
  for (x = 0; x < n; ++x) {
    for (y = 0; y < n; ++y) {
      for (z = 0; z < n; ++z) {
        // Add .5 for symmetry: 0.5, 1.5, 2.5 vs. 0, 1, 2
        part->x[0] = x + offset[0] + 0.5;
        part->x[1] = y + offset[1] + 0.5;
        part->x[2] = z + offset[2] + 0.5;
        part->v[0] = 1.0f;
        part->v[1] = 1.0f;
        part->v[2] = 1.0f;
        part->h = h;
        part->id = ++(*partId);
        part->mass = 1.0f;
        part->ti_begin = 0;
        part->ti_end = 1;
        ++part;
      }
    }
  }

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

  fprintf(file,
          "# ID  rho  rho_dh  wcount  wcount_dh  div_v  curl_v:[x y z]\n");

  for (size_t pid = 0; pid < ci->count; pid++) {
    fprintf(file, "%6llu %f %f %f %f %f %f %f %f\n", ci->parts[pid].id,
            ci->parts[pid].rho, ci->parts[pid].rho_dh,
            ci->parts[pid].density.wcount, ci->parts[pid].density.wcount_dh,
            ci->parts[pid].div_v, ci->parts[pid].density.rot_v[0],
            ci->parts[pid].density.rot_v[1], ci->parts[pid].density.rot_v[2]);
  }

  fprintf(file, "# -----------------------------------\n");

  for (size_t pjd = 0; pjd < cj->count; pjd++) {
    fprintf(file, "%6llu %f %f %f %f %f %f %f %f\n", cj->parts[pjd].id,
            cj->parts[pjd].rho, cj->parts[pjd].rho_dh,
            cj->parts[pjd].density.wcount, cj->parts[pjd].density.wcount_dh,
            cj->parts[pjd].div_v, cj->parts[pjd].density.rot_v[0],
            cj->parts[pjd].density.rot_v[1], cj->parts[pjd].density.rot_v[2]);
  }

  fclose(file);
}

/* Just a forward declaration... */
void runner_dopair1_density(struct runner *r, struct cell *ci, struct cell *cj);

int main(int argc, char *argv[]) {
  size_t particles = 0, runs = 0, volume, type = 0;
  double offset[3] = {0, 0, 0}, h = 1.1255;  // * DIM/PARTS_PER_AXIS == * 1
  struct cell *ci, *cj;
  struct space space;
  struct engine engine;
  struct runner runner;
  char c;
  static unsigned long long partId = 0;
  ticks tic, toc, time;

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  while ((c = getopt(argc, argv, "h:p:r:t:")) != -1) {
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
    }
  }

  if (h < 0 || particles == 0 || runs == 0 || type > 2) {
    printf(
        "\nUsage: %s -p PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates a cell pair, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density."
        "\n\nOptions:"
        "\n-t TYPE=0          - cells share face (0), edge (1) or corner (2)"
        "\n-h DISTANCE=1.1255 - smoothing length\n",
        argv[0]);
    exit(1);
  }

  space.periodic = 0;
  volume = particles * particles * particles;
  message("particles: %zu B\npositions: 0 B", 2 * volume * sizeof(struct part));

  ci = make_cell(particles, offset, h, &partId);
  for (size_t i = 0; i < type + 1; ++i) offset[i] = particles;
  cj = make_cell(particles, offset, h, &partId);

  for (int i = 0; i < 3; ++i) {
    space.h_max = h;
    space.dt_step = 0.1;
  }

  engine.s = &space;
  engine.time = 0.1f;
  runner.e = &engine;

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
    if (i % 50 == 0) dump_particle_fields("swift_dopair.dat", ci, cj);
  }

  /* Output timing */
  message("SWIFT calculation took       %lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  zero_particle_fields(ci);
  zero_particle_fields(cj);

  tic = getticks();

  /* Run the test */
  pairs_all_density(&runner, ci, cj);

  toc = getticks();

  /* Dump */
  dump_particle_fields("brute_force.dat", ci, cj);

  /* Output timing */
  message("Brute force calculation took %lli ticks.", toc - tic);

  /* Clean things to make the sanitizer happy ... */
  clean_up(ci);
  clean_up(cj);

  return 0;
}
