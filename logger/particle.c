#include "particle.h"
#include "header.h"
#include "io.h"
#include "timeline.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/**
 * @brief print a particle
 *
 * @param p #particle particle to print
 */
void particle_print(const struct particle *p) {
  printf("ID:            %lu\n", p->id);
  printf("Mass:          %g\n", p->mass);
  printf("Time:          %g\n", p->time);
  printf("Cutoff Radius: %g\n", p->h);
  printf("Position:      (%g, %g, %g)\n", p->pos[0], p->pos[1], p->pos[2]);
  printf("Velocity:      (%g, %g, %g)\n", p->vel[0], p->vel[1], p->vel[2]);
  printf("Acceleration:  (%g, %g, %g)\n", p->acc[0], p->acc[1], p->acc[2]);
  printf("Entropy:       %g\n", p->entropy);
  printf("Density:       %g\n", p->density);
}

/**
 * @brief Check if dump data type are compatible with the particle type
 *
 * @param h #header structure of the file
 *
 * @return error code
 */
int particle_check_data_type(__attribute__((unused)) const struct header *h) {
  printf("TODO check_data_type\n");
  return 1;
}

/**
 * @brief initialize a particle
 *
 * @param part #particle particle to initialize
 */
void particle_init(struct particle *part) {
  for (size_t k = 0; k < DIM; k++) {
    part->pos[k] = 0;
    part->vel[k] = 0;
    part->acc[k] = 0;
  }

  part->entropy = -1;
  part->density = -1;
  part->h = -1;
  part->mass = -1;
  part->id = SIZE_MAX;
}

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
                        const char *field, const size_t size) {
  void *p = NULL;

  if (strcmp("position", field) == 0) {
    p = &part->pos;
  } else if (strcmp("velocity", field) == 0) {
    p = &part->vel;
  } else if (strcmp("acceleration", field) == 0) {
    p = &part->acc;
  } else if (strcmp("entropy", field) == 0) {
    p = &part->entropy;
  } else if (strcmp("cutoff radius", field) == 0) {
    p = &part->h;
  } else if (strcmp("density", field) == 0) {
    p = &part->density;
  } else if (strcmp("consts", field) == 0) {
    p = malloc(size);
  } else {
    error(ENOTSUP, "Type %s not defined", field);
  }

  int error_code = io_read_data(map, size, p, offset);

  if (error_code != 0) return error_code;

  if (strcmp("consts", field) == 0) {
    part->mass = 0;
    part->id = 0;
    memcpy(&part->mass, p, sizeof(float));
    p += sizeof(float);
    memcpy(&part->id, p, sizeof(size_t));
    p -= sizeof(float);
    free(p);
  }
  return error_code;
}

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
                  struct time_array *times) {
  int error_code = 0;
  size_t mask = 0;
  size_t h_offset = 0;

  particle_init(part);

  error_code = io_read_mask(h, map, offset, &mask, &h_offset);
  if (error_code != 0) return error_code;

  if (mask != 127) error(EIO, "Unexpected mask: %lu", mask);

  for (size_t i = 0; i < h->nber_mask; i++) {
    if (mask & h->masks[i]) {
      error_code = particle_read_field(part, map, offset, h->masks_name[i],
                                       h->masks_size[i]);
      if (error_code != 0) return error_code;
    }
  }

  if (times) /* move offset by 1 in order to be in the required chunk */
    part->time = time_array_get_time(times, *offset - 1);
  else
    part->time = -1;

  /* end of const case */
  if (reader == reader_const) return 0;

  /* read next particle */
  struct particle part_next;

  if (!h->forward_offset) error(ENOSYS, "TODO");

  if (h_offset == 0) return 0;
  /* get absolute offset of next particle */
  h_offset += *offset - header_get_mask_size(h, mask) - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;

  part_next.time = time_array_get_time(times, h_offset);

  /* previous part exists */
  error_code = particle_read(&part_next, h, map, &h_offset, part_next.time,
                             reader_const, times);
  if (error_code != 0) return error_code;

  error_code = particle_interpolate(part, &part_next, time);
  if (error_code != 0) return error_code;

  return 0;
}

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
                         const double time) {

  if (!part_curr) error(EFAULT, "part_curr is NULL");
  if (!part_next) error(EFAULT, "part_next is NULL");

#ifdef SWIFT_DEBUG_CHECKS
  if (part_next->time <= part_curr->time)
    error(EIO, "Wrong particle order (next before current)");
  if ((time < part_curr->time) || (part_next->time < time))
    error(EIO,
          "Interpolating, not extrapolating (particle time: %f, "
          "interpolating time: %f, next particle time: %f)",
          part_curr->time, time, part_next->time);
#endif

  double scaling = part_next->time - part_curr->time;

  scaling = (time - part_curr->time) / scaling;

  double tmp;
  float ftmp;

  /* interpolate vectors */
  for (size_t i = 0; i < DIM; i++) {
    tmp = (part_next->pos[i] - part_curr->pos[i]);
    part_curr->pos[i] += tmp * scaling;

    ftmp = (part_next->vel[i] - part_curr->vel[i]);
    part_curr->vel[i] += ftmp * scaling;

    ftmp = (part_next->acc[i] - part_curr->acc[i]);
    part_curr->acc[i] += ftmp * scaling;
  }

  /* interpolate scalars */
  ftmp = (part_next->entropy - part_curr->entropy);
  part_curr->entropy += ftmp * scaling;

  /* set time */
  part_curr->time = time;

  return 0;
}
