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

/* Local header */
#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_particle.h"
#include "logger_reader.h"
#include "swift.h"

/* Tests header */
#include "generate_log.h"

#define number_parts 100
#define max_step 99

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
  int step = 0;
  int init_log_all_done = -1;

  /* Number of particle found during this time step. */
  int count = 0;
  /* Set it to an impossible value in order to flag it. */
  const size_t id_flag = 5 * number_parts;
  long long previous_id = id_flag;

  /* Loop over each record. */
  for (size_t offset = reader_read_record(reader, &lp, &time, &is_particle,
                                          logfile->header.offset_first_record);
       offset < logfile->log.mmap_size;
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
      if (lp.id >= number_parts) error("Wrong id %lli", lp.id);

      struct part *p = &parts[lp.id];

      /* Check the record's data. */
      for (int i = 0; i < 3; i++) {
        /* in the first index, we are storing the step information. */
        if (i == 0) {
          double tmp = step;
          /* At the end, we are not updating the particle */
          if (step >= max_step) {
            tmp = max_step - max_step % p->time_bin;
          }
          assert(tmp == lp.pos[i]);
        } else
          assert(p->x[i] == lp.pos[i]);
        assert(p->v[i] == lp.vel[i]);
        assert(p->a_hydro[i] == lp.acc[i]);
      }

      assert(p->entropy == lp.entropy);
      assert(p->mass == lp.mass);

      /* Check optional fields. */
      int number_steps = step / p->time_bin;
      if (number_steps % period_h == 0 || step > max_step) {
        assert(p->h == lp.h);
      } else {
        assert(-1 == lp.h);
      }
      if (number_steps % period_rho == 0 || step > max_step) {
        assert(p->rho == lp.density);
      } else {
        assert(-1 == lp.density);
      }
    }
    /* Time stamp case. */
    else {
      message("Step: %i", step);
      /* Check if we have the current amount of particles in previous step. */
      if (step != 0 && count != get_number_active_particles(step, parts))
        error(
            "The reader did not find the correct number of particles during "
            "step %i: %i != %i",
            step, count, get_number_active_particles(step, parts));

      /* Avoid the initial log */
      if (init_log_all_done > 0) {
        step += 1;
      }

      init_log_all_done += 1;

      /* Reset some variables. */
      previous_id = id_flag;
      count = 0;

      /* Check the record's data. */
      const int tmp_step = step >= max_step ? max_step : step;
      assert(time == get_double_time(tmp_step));
    }
  }
}

int main(int argc, char *argv[]) {

  /*
    First generate the file.
  */

  message("Generating the dump.");

  /* Create required structures. */
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

  /* Write a 'simulation' */
  generate_log(&params, parts, xparts, number_parts);

  /*
    Then read the file.
  */

  message("Reading the header.");

  /* Generate required structure for reading. */
  struct logger_reader reader;

  /* Set verbose level. */
  reader.verbose = 1;

  /* Read the header. */
  char basename[200];
  parser_get_param_string(&params, "Logger:basename", basename);
  logger_reader_init(&reader, basename, /* verbose */ 1);

  /*
    Finally check everything.
  */

  check_data(&reader, parts, xparts);

  /* Do some cleanup. */
  free(parts);
  free(xparts);

  return 0;
}
