/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Loic Hausammann (loic.hausammann@epfl.ch).
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

#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_particle.h"
#include "logger_reader.h"
#include "swift.h"

#define number_parts 100
/* Not all the fields are written at every step.
 * Here we define how often a few fields are written.
 */
#define period_rho 2
#define period_h 4

/**
 * @brief Initialize the particles.
 *
 * @param p The array of #part.
 * @param xp The array of #xpart.
 */
void init_particles(struct part *p, struct xpart *xp) {
  struct hydro_space hs;

  for (int i = 0; i < number_parts; i++) {
    /* Set internal energy. */
    hydro_set_init_internal_energy(&p[i], 100);

    /* Initialize particle. */
    hydro_first_init_part(&p[i], &xp[i]);
    hydro_init_part(&p[i], &hs);

    for (int j = 0; j < 3; j++) {
      p[i].x[j] = i;
      p[i].v[j] = (j == 0) ? -1 : 0;
      p[i].a_hydro[j] = (j == 1) ? 1e-2 : 0;
    }
    p[i].h = 15;
    p[i].rho = 50;
    p[i].id = i;
    hydro_set_mass(&p[i], 1.5);
    xp[i].logger_data.last_offset = 0;

    /* Add time bin in order to skip particles. */
    p[i].time_bin = (i % 10) + 1;
  }
}

/** Provides a integer time given the step number.*/
integertime_t get_integer_time(int step) { return step; }

/** Provides a double time given the step number. */
double get_double_time(int step) {
  const double time_base = 1e-4;
  return step * time_base;
}

/**
 * @brief Write a few particles during multiple time steps.
 *
 * As only the logger is tested, there is no need to really
 * evolve the particles.
 */
void write_particles(struct logger_writer *log, struct part *parts,
                     struct xpart *xparts) {

  const int number_steps = 100;

  /* Loop over all the steps. */
  for (int i = 0; i < number_steps; i++) {
    integertime_t ti_int = get_integer_time(i);
    double ti_double = get_double_time(i);

    /* Mark the current time step in the particle logger file. */
    logger_log_timestamp(log, ti_int, ti_double, &log->timestamp_offset);
    /* Make sure that we have enough space in the particle logger file
     * to store the particles in current time step. */
    logger_ensure_size(log, number_parts, /* number gpart */ 0, 0);

    /* Loop over all the particles. */
    for (int j = 0; j < number_parts; j++) {

      /* Skip some particles. */
      if (i % parts[j].time_bin != 0) continue;

      /* Write a time information to check that the correct particle is read. */
      parts[j].x[0] = i;

      /* Write this particle. */
      unsigned int mask =
          logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
          logger_mask_data[logger_a].mask | logger_mask_data[logger_u].mask |
          logger_mask_data[logger_consts].mask;

      int number_particle_step = i / parts[j].time_bin;

      if (number_particle_step % period_h == 0)
        mask |= logger_mask_data[logger_h].mask;
      if (number_particle_step % period_rho == 0)
        mask |= logger_mask_data[logger_rho].mask;

      logger_log_part(log, &parts[j], mask, &xparts[j].logger_data.last_offset);
    }

    // TODO write index files.
  }

  /* Mark the current time step in the particle logger file. */
  integertime_t ti_int = get_integer_time(number_steps);
  double ti_double = get_double_time(number_steps);
  logger_log_timestamp(log, ti_int, ti_double, &log->timestamp_offset);
}

/** Count the number of active particles. */
int get_number_active_particles(int step, struct part *p) {
  int count = 0;
  for (int i = 0; i < number_parts; i++) {
    if (step % p[i].time_bin == 0) count += 1;
  }
  return count;
}
/**
 * @brief Check that the reader contains the correct data
 *
 * @param reader The #logger_reader.
 */
void check_data(struct logger_reader *reader, struct part *parts,
                struct xpart *xparts) {

  /* No need to check the header, this is already done in testHeader.c */

  /* Get required structures. */
  struct logger_logfile *logfile = &reader->log;

  struct logger_particle lp;
  logger_particle_init(&lp);

  /* Define a few variables */
  double time = get_double_time(0);
  int is_particle = 0;
  int step = -1;

  /* Number of particle found during this time step. */
  int count = 0;
  /* Set it to an impossible value in order to flag it. */
  const size_t id_flag = 5 * number_parts;
  size_t previous_id = id_flag;

  /* Loop over each record. */
  for (size_t offset = reader_read_record(reader, &lp, &time, &is_particle,
                                          logfile->header.offset_first_record);
       offset < logfile->log.file_size;
       offset = reader_read_record(reader, &lp, &time, &is_particle, offset)) {

    /* Do the particle case */
    if (is_particle) {
      count += 1;

      /*
        Check that we are really increasing the id in the logfile.
        See the writing part to see that we are always increasing the id.
      */
      if (previous_id != id_flag && previous_id >= lp.id) {
        error("Wrong particle found");
        previous_id = lp.id;
      }

      /* Get the corresponding particle */
      if (lp.id >= number_parts) error("Wrong id %zi", lp.id);

      struct part *p = &parts[lp.id];

      /* Check the record's data. */
      for (int i = 0; i < 3; i++) {
        /* in the first index, we are storing the step information. */
        if (i == 0)
          assert(step == lp.pos[i]);
        else
          assert(p->x[i] == lp.pos[i]);
        assert(p->v[i] == lp.vel[i]);
        assert(p->a_hydro[i] == lp.acc[i]);
      }

      assert(p->entropy == lp.entropy);
      assert(p->mass == lp.mass);

      /* Check optional fields. */
      int number_steps = step / p->time_bin;
      if (number_steps % period_h == 0) {
        assert(p->h == lp.h);
      } else {
        assert(-1 == lp.h);
      }
      if (number_steps % period_rho == 0) {
        assert(p->rho == lp.density);
      } else {
        assert(-1 == lp.density);
      }
    }
    /* Time stamp case. */
    else {

      /* Check if we have the current amount of particles in previous step. */
      if (step != -1 && count != get_number_active_particles(step, parts))
        error(
            "The reader did not find the correct number of particles during "
            "step %i",
            step);

      step += 1;

      /* Reset some variables. */
      previous_id = id_flag;
      count = 0;

      /* Check the record's data. */
      assert(time == get_double_time(step));
    }
  }
}

int main(int argc, char *argv[]) {

  /*
    First generate the file.
  */

  message("Generating the dump.");

  /* Create required structures. */
  struct logger_writer log;
  struct swift_params params;
  char filename[200] = "testLogfileReader.yml";

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

  init_particles(parts, xparts);

  /* Initialize the logger. */
  logger_init(&log, &params);

  /* get dump filename. */
  char dump_filename[PARSER_MAX_LINE_SIZE];
  message("%s", log.base_name);
  strcpy(dump_filename, log.base_name);
  strcat(dump_filename, ".dump");

  /* Write file header. */
  logger_write_file_header(&log);

  /* Write particles. */
  write_particles(&log, parts, xparts);

  /* clean memory */
  logger_free(&log);
  /*
    Then read the file.
  */

  message("Reading the header.");

  /* Generate required structure for reading. */
  struct logger_reader reader;

  /* Set verbose level. */
  reader.verbose = 1;

  /* Read the header. */
  logger_reader_init(&reader, dump_filename, /* verbose */ 1);

  /*
    Finally check everything.
  */

  check_data(&reader, parts, xparts);

  /* Do some cleanup. */
  free(parts);
  free(xparts);

  return 0;
}
