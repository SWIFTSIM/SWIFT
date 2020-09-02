/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include config */
#include "../../config.h"

/* Local headers */
#include "engine.h"
#include "hydro.h"
#include "logger.h"
#include "logger_io.h"

/* Not all the fields are written at every step.
 * Here we define how often a few fields are written.
 */
#define period_rho 2
#define period_h 4
#define const_time_base 1e-4

/**
 * @brief Generate the data of a bunch of particles.
 *
 * @param parts The list of particles.
 * @param xparts The list of extra particles.
 * @param nparts The number of particles.
 */
void generate_particles(struct part *parts, struct xpart *xparts,
                        size_t nparts) {
  struct hydro_space hs;

  for (size_t i = 0; i < nparts; i++) {
    /* Set internal energy. */
    hydro_set_init_internal_energy(&parts[i], 100);

    /* Initialize particle. */
    hydro_first_init_part(&parts[i], &xparts[i]);
    hydro_init_part(&parts[i], &hs);
    logger_part_data_init(&xparts[i].logger_data);

    for (int j = 0; j < 3; j++) {
      parts[i].x[j] = 0;
      parts[i].v[j] = (j == 0) ? -1 : 0;
      parts[i].a_hydro[j] = (j == 1) ? 1e-2 : 0;
      xparts[i].a_grav[j] = 0;
    }
    parts[i].h = 15;
    parts[i].rho = 50;
    parts[i].id = i;
    hydro_set_mass(&parts[i], 1.5);

    /* Add time bin in order to skip particles. */
    parts[i].time_bin = (i % 10) + 1;
  }
}

/** Provides a integer time given the step number.*/
integertime_t get_integer_time(int step) { return step; }

/** Provides a double time given the step number. */
double get_double_time(int step) { return step * const_time_base; }

/**
 * @brief Write a few particles during multiple time steps.
 *
 * As only the logger is tested, there is no need to really
 * evolve the particles.
 *
 * @param log The #logger_writer.
 * @param e The #engine.
 */
void write_particles(struct logger_writer *log, struct engine *e) {

  size_t nparts = e->total_nr_parts;
  struct part *parts = e->s->parts;
  struct xpart *xparts = e->s->xparts;

  const int number_steps = 100;
  const int number_index = 5;

  /* Loop over all the steps. */
  for (int i = 0; i < number_steps; i++) {
    e->time = get_double_time(i);
    e->ti_current = get_integer_time(i);
    /* Dump an index file if required */
    if (i % (number_steps / number_index) == number_index - 1) {
      engine_dump_index(e);
    }
    integertime_t ti_int = get_integer_time(i);
    double ti_double = get_double_time(i);

    /* Mark the current time step in the particle logger file. */
    logger_log_timestamp(log, ti_int, ti_double, &log->timestamp_offset);
    /* Make sure that we have enough space in the particle logger file
     * to store the particles in current time step. */
    logger_ensure_size(log, nparts, /* number gpart */ 0, 0);

    /* Loop over all the particles. */
    for (size_t j = 0; j < nparts; j++) {

      /* Skip some particles. */
      if (i % parts[j].time_bin != 0) continue;

      /* Write a time information to check that the correct particle is read. */
      parts[j].x[0] = i;

      // TODO write only a few masks at the time

      logger_log_part(log, &parts[j], &xparts[j], e, /* log_all */ 0,
                      /* special flags */ 0);
    }
  }
}

void generate_log(struct swift_params *params, struct part *parts,
                  struct xpart *xparts, size_t nparts) {
  /* Initialize the particles */
  generate_particles(parts, xparts, nparts);

  /* initialize the engine */
  struct engine e;
  e.policy = engine_policy_hydro;
  e.total_nr_parts = nparts;
  e.total_nr_gparts = 0;
  e.total_nr_sparts = 0;
  e.total_nr_bparts = 0;
  e.verbose = 1;
  e.ti_current = 0;
  e.time = 0;
  e.time_base = const_time_base;
  e.time_begin = 0;
  threadpool_init(&e.threadpool, 1);
  struct space s;
  e.s = &s;
  s.xparts = xparts;
  s.parts = parts;
  s.gparts = NULL;
  s.nr_parts = nparts;
  s.nr_gparts = 0;
  s.nr_sparts = 0;
  s.nr_bparts = 0;
  s.nr_inhibited_parts = 0;
  s.nr_inhibited_gparts = 0;
  s.nr_inhibited_sparts = 0;
  s.nr_inhibited_bparts = 0;
  s.nr_extra_gparts = 0;
  s.nr_extra_parts = 0;
  s.nr_extra_sparts = 0;
  s.nr_extra_bparts = 0;
  struct logger_writer log;
  e.logger = &log;

  /* Initialize the writer */
  logger_init(&log, &e, params);

  /* Write file header */
  logger_write_file_header(&log);

  /* Mark the current time step in the particle logger file. */
  logger_log_timestamp(&log, e.ti_current, e.time, &log.timestamp_offset);
  /* Make sure that we have enough space in the particle logger file
   * to store the particles in current time step. */
  logger_ensure_size(&log, nparts, /* number gpart */ 0, 0);

  /* Log all the particles before starting */
  logger_log_all_particles(&log, &e);
  engine_dump_index(&e);

  /* Write particles */
  write_particles(&log, &e);

  /* Write a sentinel timestamp */
  logger_log_timestamp(e.logger, e.ti_current, e.time,
                       &e.logger->timestamp_offset);

  /* Write all the particles at the end */
  logger_log_all_particles(e.logger, &e);

  /* Write a sentinel timestamp */
  logger_log_timestamp(e.logger, e.ti_current, e.time,
                       &e.logger->timestamp_offset);

  /* Cleanup the memory */
  logger_free(&log);
}
