/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *                    Florian Cabot (florian.cabot@epfl.ch)
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

/* Include configuration */
#include "config.h"

/* Standard include */
#include <stdlib.h>

/* Local include */
#include "generate_log.h"
#include "hydro.h"
#include "logger_reader.h"

#define number_steps 10.
#define number_parts 100

/**
 * This function test the logger in the VR mode.
 * The idea is to simply read a snapshot at a given time and
 * then simply advance in time the particles.
 */
int main(int argc, char *argv[]) {
  /* Create required structures. */
  struct swift_params params;
  char filename[200] = "testVR.yml";

  /* Read parameters. */
  parser_read_file(filename, &params);

  /* Initialize the particles. */
  struct part *parts;
  if ((parts = (struct part *)malloc(sizeof(struct part) * number_parts)) ==
      NULL)
    error("Failed to allocate particles array.");

  struct xpart *xparts;
  if ((xparts = (struct xpart *)malloc(sizeof(struct xpart) * number_parts)) ==
      NULL)
    error("Failed to allocate xparticles array.");

  /* Write a 'simulation' */
  generate_log(&params, parts, xparts, number_parts);

  /* Initialize the reader */
  struct logger_reader reader;
  char basename[200];
  parser_get_param_string(&params, "Logger:basename", basename);
  logger_reader_init(&reader, basename,
                     /* Verbose */ 0);

  /* Read the time limits */
  double begin = logger_reader_get_time_begin(&reader);
  double end = logger_reader_get_time_end(&reader);

  /* Set the time */
  message("Time begin: %f end: %f", begin, end);
  logger_reader_set_time(&reader, begin);

  /* Get the number of particles */
  int n_type = 0;
  uint64_t n_tot = 0;
  const uint64_t *n_parts =
      logger_reader_get_number_particles(&reader, &n_type);
  for (int i = 0; i < n_type; i++) {
    n_tot += n_parts[i];
  }

  /* Allocate the particles memory */
  struct logger_particle *particles =
      malloc(n_tot * sizeof(struct logger_particle));

  logger_reader_read_all_particles(&reader, begin, logger_reader_const,
                                   particles, n_tot);

  /* Loop over time for a single particle */
  size_t id = 0;
  struct logger_particle p = particles[id];
  for (double t = begin; t < end; t += (end - begin) / number_steps) {
    /* Get the offset of the given time */
    size_t o = logger_reader_get_next_offset_from_time(&reader, t);
    message("time: %f offset: %ld", t, o);

    /* Read the next particle */
    struct logger_particle n;
    logger_reader_get_next_particle(&reader, &p, &n, o);

    message("Particle %zi: %f %f %f %f", id, p.pos[0], p.pos[1], p.pos[2],
            p.time);

    /* Now you can interpolate */
    logger_particle_interpolate(&p, &n, t);
  }

  /* Cleanup the memory */
  free(particles);
  logger_reader_free(&reader);
  return 0;
}
