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

#include "swift.h"
#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_particle.h"
#include "logger_reader.h"

#define number_parts 100

/**
 * @brief Initialize the particles.
 *
 * @param p The array of #part.
 * @param xp The array of #xpart.
 */
void init_particles(struct part *p, struct xpart *xp) {
  struct hydro_space hs;

  for(int i = 0; i < number_parts; i++) {
    /* Set internal energy. */
    hydro_set_init_internal_energy(&p[i], 100);

    /* Initialize particle. */
    hydro_first_init_part(&p[i], &xp[i]);
    hydro_init_part(&p[i], &hs);

    for(int j = 0; j < 3; j++) {
      p[i].x[j] = i;
      p[i].v[j] = j == 0 ? -1 : 0;
      p[i].a_hydro[j] = j == 1 ? 1e-2 : 0;
    }
    p[i].h = 15;
    p[i].rho = 50;
    p[i].id = i;
    hydro_set_mass(&p[i], 1.5);
    xp[i].logger_data.last_offset = 0;
  }
}

/**
 * @brief Provides a time given the step number.
 *
 * @param step The required step.
 */
integertime_t get_integer_time(int step) {
  return step;
}

/**
 * @brief Provides a time given the step number.
 *
 * @param step The required step.
 */
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
void write_particles(struct logger *log) {


  const int number_steps = 100;

  /* Create particles and initialize them. */
  struct part *parts;
  if ((parts = (struct part *)malloc(sizeof(struct part) * number_parts)) == NULL)
    error("Failed to allocate particles array.");

  struct xpart *xparts;
  if ((xparts = (struct xpart *)malloc(sizeof(struct xpart) * number_parts)) == NULL)
    error("Failed to allocate xparticles array.");

  init_particles(parts, xparts);

  /*
    Init logger
  */
    
  /* Do step */
  for(int i = 0; i < number_steps; i++) {
    integertime_t ti_int = get_integer_time(i);
    double ti_double = get_double_time(i);

    /* Mark the current time step in the particle logger file. */
    logger_log_timestamp(log, ti_int, ti_double,
			 &log->timestamp_offset);
    /* Make sure that we have enough space in the particle logger file
     * to store the particles in current time step. */
    logger_ensure_size(log, number_parts, /* number gpart */0, 0);

    /* Run step */
    for(int j = 0; j < number_parts; j++) {

      /* Write particle */
      /* TODO Currently writing everything, should adapt it through time */
      logger_log_part(log, &parts[j],
		      logger_mask_data[logger_x].mask |
		      logger_mask_data[logger_v].mask |
		      logger_mask_data[logger_a].mask |
		      logger_mask_data[logger_u].mask |
		      logger_mask_data[logger_h].mask |
		      logger_mask_data[logger_rho].mask |
		      logger_mask_data[logger_consts].mask,
		      &xparts[j].logger_data.last_offset);      
    }

    // TODO write index files
  }

  /* Mark the current time step in the particle logger file. */
  integertime_t ti_int = get_integer_time(number_steps);
  double ti_double = get_double_time(number_steps);
  logger_log_timestamp(log, ti_int, ti_double,
		       &log->timestamp_offset);


  /* Cleanup */
  free(parts);
  free(xparts);
}

/**
 * @brief Read a record (timestamp or particle)
 *
 * @param reader The #reader.
 * @param lp (out) The #logger_particle (if the record is a particle).
 * @param time (out) The time read (if the record is a timestamp).
 * @param is_particle Is the record a particle (or a timestamp)?
 * @param offset The offset in the file.
 *
 * @return The offset after this record.
 */
size_t read_record(struct logger_reader *reader, struct logger_particle *lp,
		   double *time, int *is_particle, size_t offset) {

  struct logger_logfile *log = &reader->log;

  /* Read mask to find out if timestamp or particle */
  size_t mask = 0;
  logger_loader_io_read_mask(&log->header, log->log.map + offset, &mask, NULL);

  /* Check if timestamp or not */
  int ind = header_get_field_index(&log->header, "timestamp");
  if (ind == -1) {
    error("File header does not contain a mask for time");
  }
  if (log->header.masks[ind].mask == mask) {
    *is_particle = 0;
    integertime_t int_time = 0;
    offset = time_read(&int_time, time, reader, offset);
  }
  else {
    *is_particle = 1;
    offset = logger_particle_read(lp, reader, offset, *time, logger_reader_const);
  }

  return offset;
}

/**
 * @brief Check that the reader contains the correct data
 *
 * @param reader The #logger_reader.
 */
void check_data(struct logger_reader *reader) {

  /* No need to check the header, this is already done in testHeader.c */

  /* Get required structures */
  struct logger_logfile *logfile = &reader->log;
  
  struct logger_particle lp;
  logger_particle_init(&lp);

  /* Generate the particles again */
  struct part *parts;
  if ((parts = (struct part *)malloc(sizeof(struct part) * number_parts)) == NULL)
    error("Failed to allocate particles array.");

  struct xpart *xparts;
  if ((xparts = (struct xpart *)malloc(sizeof(struct xpart) * number_parts)) == NULL)
    error("Failed to allocate xparticles array.");

  init_particles(parts, xparts);

  double time = get_double_time(0);
  int is_particle = 0;
  int step = 0;

  /* Loop over each record */
  for(size_t offset = read_record(reader, &lp, &time, &is_particle, logfile->header.offset_first_record);
      offset < logfile->log.file_size;
      offset = read_record(reader, &lp, &time, &is_particle, offset)) {

    if (is_particle) {
      if (lp.id >= number_parts)
	error("Wrong id %zi", lp.id);

      struct part *p = &parts[lp.id];

      for(int i = 0; i < 3; i++) {
	assert(p->x[i] == lp.pos[i]);
	assert(p->v[i] == lp.vel[i]);
	assert(p->a_hydro[i] == lp.acc[i]);
      }

      assert(p->entropy == lp.entropy);
      assert(p->h == lp.h);
      assert(p->rho == lp.density);
      assert(p->mass == lp.mass);
    }
    else {
      assert(time == get_double_time(step));

      step += 1;
    }
  }
  
}


int main(int argc, char *argv[]) {

  /*
    First generate the file.
  */

  message("Generating the dump.");

  /* Create required structures. */
  struct logger log;
  struct swift_params params;
  char filename[200] = "testReader.yml";

  /* Read parameters. */
  parser_read_file(filename, &params);

  /* Initialize the logger. */
  logger_init(&log, &params);

  /* get dump filename */
  char dump_filename[PARSER_MAX_LINE_SIZE];
  message("%s", log.base_name);
  strcpy(dump_filename, log.base_name);
  strcat(dump_filename, ".dump");

  /* Write file header. */
  logger_write_file_header(&log);

  /* Write particles. */
  write_particles(&log);

  /* clean memory */
  logger_clean(&log);
  /*
    Then read the file.
  */

  message("Reading the header.");

  /* Generate required structure for reading. */
  struct logger_reader reader;

  /* Set verbose level */
  reader.verbose = 1;
  
  /* Read the header */
  logger_reader_init(&reader, dump_filename, /* verbose */ 1);

  /*
    Finally check everything
  */

  check_data(&reader);

  return 0;
}
