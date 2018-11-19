#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "header.h"
#include "timeline.h"
#include "tools.h"

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

/**
 * @brief print a particle
 *
 * @param p #particle particle to print
 */
void particle_print(const struct particle *p);

/**
 * @brief read a particle in the dump file
 *
 * @param part #particle particle to update
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset offset of the chunk to read (update it to the end of the chunk)
 * @param time time to interpolate (if #reader_type is an interpolating one)
 * @param reader #reader_type
 * @param times #time_array times in the dump
 *
 * @return error code
 */
int particle_read(struct particle *part, const struct header *h, void *map,
                  size_t *offset, const double time, const int reader,
                  struct time_array *times);

/**
 * @brief Check if dump data type are compatible with the particle type
 *
 * @param h #header structure of the file
 *
 * @return error code
 */
int particle_check_data_type(const struct header *h);

/**
 * @brief initialize a particle
 *
 * @param part #particle particle to initialize
 */
void particle_init(struct particle *part);

/**
 * @brief read a single field for a particle
 *
 * @param part @particle particle to update
 * @param map file mapping
 * @param offset In: read position, Out: input shifted by the required amount of
 * data
 * @param field field to read
 * @param size number of bits to read
 *
 * @return error code
 */
int particle_read_field(struct particle *part, void *map, size_t *offset,
                        const char *field, const size_t size);

/**
 * @brief interpolate two particles at a given time
 *
 * @param part_curr #particle In: current particle (before time), Out:
 * interpolated particle
 * @param part_next #particle next particle (after time)
 * @param time interpolation time
 *
 * @return error code
 */
int particle_interpolate(struct particle *part_curr,
                         const struct particle *part_next,
                         const double time);

#endif  //__PARTICLE_H__
