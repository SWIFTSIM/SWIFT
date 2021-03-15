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
#include "logger.h"
#include "logger_reader.h"

#define number_steps 10.
#define number_parts 100

/**
 * This function test the logger in the Virtual Reality mode.
 * The idea is to simply read a snapshot at a given time and
 * then simply advance in time the particles.
 */
int main(int argc, char *argv[]) {
  /* Create required structures. */
  struct swift_params params;
  char filename[200] = "testVirtualReality.yml";

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
  strcat(basename, "_0000");
  logger_reader_init(&reader, basename,
                     /* Verbose */ 0, /* number_threads */1);

  /* Read the time limits */
  double begin = logger_reader_get_time_begin(&reader);
  double end = logger_reader_get_time_end(&reader);

  /* Set the time */
  message("Time begin: %f end: %f", begin, end);
  logger_reader_set_time(&reader, begin);

  /* Create the variables for the number of particles */
  const int n_type = swift_type_count;
  uint64_t *n_parts = (uint64_t *)malloc(n_type * sizeof(uint64_t));
  int *read_types = (int *)malloc(n_type * sizeof(int));
  if (read_types == NULL || n_parts == NULL) {
    error("Failed to allocate arrays.");
  }

  /* Set the flags in order to read everything */
  for (int i = 0; i < n_type; i++) {
    read_types[i] = 1;
  }

  /* Get the number of particles */
  logger_reader_get_number_particles(&reader, n_parts, read_types);

  uint64_t n_tot = 0;
  for (int i = 0; i < n_type; i++) {
    n_tot += n_parts[i];
  }

  /* Allocate the particles memory */
  double *pos = (double *)malloc(n_tot * 3 * sizeof(double));
  long long *ids = (long long *)malloc(n_tot * sizeof(long long));

  /* Create the list of fields. */
  const int n_fields = 2;
  int *required_fields = (int *)malloc(n_fields * sizeof(int));
  const struct header *h = &reader.log.header;
  for (int i = 0; i < n_fields; i++) {
    required_fields[i] = -1;
  }
  for (int j = 0; j < h->masks_count; j++) {
    if (strcmp(h->masks[j].name, "Coordinates") == 0) {
      required_fields[0] = j;
    } else if (strcmp(h->masks[j].name, "ParticleIDs") == 0) {
      required_fields[1] = j;
    }
  }
  if (required_fields[0] == -1) {
    error("Coordinates not found");
  }
  if (required_fields[1] == -1) {
    error("ParticleIDs not found.");
  }

  /* Create the output */
  void **output = malloc(n_fields * sizeof(void *));
  output[0] = (void *)pos;
  output[1] = (void *)ids;

  /* Loop over time for a single particle */
  int part_ind = 0;
  for (double t = begin; t < end; t += (end - begin) / number_steps) {
    /* Set the time of the next reading */
    logger_reader_set_time(&reader, t);

    /* Read the next time */
    logger_reader_read_all_particles(&reader, t, logger_reader_lin,
                                     required_fields, n_fields, output,
                                     n_parts);

    message("Particle %lli: %f %f %f %f", ids[part_ind], pos[3 * part_ind + 0],
            pos[3 * part_ind + 1], pos[3 * part_ind + 2], t);
  }

  /* Cleanup the memory */
  free(required_fields);
  free(ids);
  free(pos);
  free(parts);
  free(xparts);
  logger_reader_free(&reader);
  free(output);
  return 0;
}
