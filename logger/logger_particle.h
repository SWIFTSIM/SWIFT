#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "logger_header.h"
#include "logger_tools.h"
#include "logger_time.h"

#include <stdio.h>
#include <stdlib.h>

#define DIM 3

struct particle {
  /* position */
  double pos[DIM];

  /* velocity */
  float vel[DIM];

  /* acceleration */
  float acc[DIM];

  /* entropy */
  float entropy;

  /* cutoff radius */
  float h;

  /* density */
  float density;

  /* mass */
  float mass;

  /* id */
  size_t id;

  /* time */
  double time;
};

enum reader_type {
  reader_const,
  reader_lin,
};

void particle_print(const struct particle *p);

void particle_read(struct particle *part, const struct header *h, void *map,
                  size_t *offset, const double time, const int reader,
                  struct time_array *times);

void particle_init(struct particle *part);

void particle_read_field(struct particle *part, void *map, size_t *offset,
                        const char *field, const size_t size);

void particle_interpolate(struct particle *part_curr,
                         const struct particle *part_next, const double time);

#endif  //__PARTICLE_H__
