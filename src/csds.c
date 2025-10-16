/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Config parameters. */
#include "part_type.h"
#include <config.h>

#ifdef HAVE_POSIX_FALLOCATE /* Are we on a sensible platform? */
#ifdef WITH_CSDS

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Define the particles first */
#include "part.h"

/* This object's header. */
#include "csds.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "chemistry_csds.h"
#include "engine.h"
#include "error.h"
#include "gravity_csds.h"
#include "hydro_csds.h"
#include "star_formation_csds.h"
#include "stars_csds.h"
#include "units.h"


/* Max number of entries that can be written for all given particle type */
static const int csds_max_size_output_list = 200;

/**
 * @brief log all particles in the engine.
 *
 * If this is the first log of all the particles,
 * we include a flag and write the type of particle.
 * This will be used by the reader to generate the index files.
 *
 * TODO use threadpool + csds function for multiple particles.
 * @param log The #csds_writer.
 * @param e The #engine.
 * @param flag The flag to use when writing the particles
 */
void csds_log_all_particles(struct csds_writer *log, const struct engine *e,
			    const enum csds_special_flags flag) {

  /* Ensure that enough space is available. */
  csds_ensure_size(log, e);

  /* csds_ensure_size is tracked separately. */
  ticks tic = getticks();

  /* some constants. */
  const struct space *s = e->s;

  /* log the parts. */
  for (size_t i = 0; i < s->nr_parts; i++) {
    struct part *p = &s->parts[i];
    struct xpart *xp = &s->xparts[i];
    if (!part_is_inhibited(p, e) && p->time_bin != time_bin_not_created) {
      csds_log_part(log, p, xp, e, /* log_all_fields */ 1, flag,
		    /* flag_data */ 0);
    }
  }

  /* log the gparts */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    struct gpart *gp = &s->gparts[i];
    if (!gpart_is_inhibited(gp, e) && gp->time_bin != time_bin_not_created &&
	(gp->type == swift_type_dark_matter ||
	 gp->type == swift_type_dark_matter_background)) {
      csds_log_gpart(log, gp, e, /* log_all_fields */ 1, flag,
		     /* flag_data */ 0);
    }
  }

  /* log the parts */
  for (size_t i = 0; i < s->nr_sparts; i++) {
    struct spart *sp = &s->sparts[i];
    if (!spart_is_inhibited(sp, e) && sp->time_bin != time_bin_not_created) {
      csds_log_spart(log, sp, e, /* log_all_fields */ 1, flag,
		     /* flag_data */ 0);
    }
  }

  if (s->nr_bparts > 0) error("Not implemented");
  if (s->nr_sinks > 0) error("Not implemented");

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
	    clocks_getunit());
}

/**
 * @brief Copy the particle fields into a given buffer.
 *
 * @param log The #csds_writer
 * @param p The #part to copy.
 * @param xp The #xpart to copy.
 * @param e The #engine.
 * @param mask The mask for the fields to write.
 * @param offset The offset to the previous log.
 * @param offset_new The offset of the current record.
 * @param buff The buffer to use when writing.
 * @param special_flags The data for the special flags.
 */
void csds_copy_part_fields(const struct csds_writer *log, const struct part *p,
			   const struct xpart *xp, const struct engine *e,
			   unsigned int mask, size_t *offset, size_t offset_new,
			   char *buff, const uint32_t special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mask == 0) {
    error("You should always log at least one field.");
  }
#endif

  /* Write the header. */
  buff = csds_write_record_header(buff, &mask, offset, offset_new);

  /* --- Fixed Fields --- */
  /* Special flags (CSDS_SPECIAL_FLAGS_INDEX = 0) */
  if (mask & log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask) {
    memcpy(buff, &special_flags, log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size);
    buff += log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask;
  }

  /* Positions (CSDS_POS_INDEX = 2) */
  if (mask & log->fixed_fields[CSDS_POS_INDEX].mask) {
    memcpy(buff, ((char *)p) + offsetof(struct part, x), log->fixed_fields[CSDS_POS_INDEX].size);
    buff += log->fixed_fields[CSDS_POS_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_POS_INDEX].mask;
  }

  /* Velocities (CSDS_VEL_INDEX = 3) */
  if (mask & log->fixed_fields[CSDS_VEL_INDEX].mask) {
    memcpy(buff, ((char *)p) + offsetof(struct part, v), log->fixed_fields[CSDS_VEL_INDEX].size);
    buff += log->fixed_fields[CSDS_VEL_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_VEL_INDEX].mask;
  }

  /* Accelerations (CSDS_ACCEL_INDEX = 4) */
  if (mask & log->fixed_fields[CSDS_ACCEL_INDEX].mask) {
    /* Since Acceleration is a calculated field for hydro, we use the
       conversion function. The `fixed_fields` only has the size/mask, not the
       conversion. We must ensure the `csds_hydro_define_fields` correctly
       registers the conversion. */
    /* **CRITICAL:** The macro approach used here is fragile for fixed fields
       that need conversion. For now, we manually call the function for the
       known 'Accelerations' field.  This assumes `csds_hydro_convert_acc` is
       the function to use. */
    csds_hydro_convert_acc(p, xp, e, buff);
    buff += log->fixed_fields[CSDS_ACCEL_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_ACCEL_INDEX].mask;
  }

  /* --- Particle-Type-Specific Fields (Hydro/Gas) --- */

  const unsigned int type_mask = log->part_type_masks[swift_type_gas];
  if (mask & type_mask) {
    for (int i = 0; i < log->number_fields[swift_type_gas]; i++) {
      struct csds_field *field = log->part_type_fields[swift_type_gas] + i;

      // The particle mask is set, so we write ALL fields in the list,
      // as they all share the same mask bit.

      /* Do we have a conversion function? */
      if (field->conversion_hydro) {
	char *tmp_buff = field->conversion_hydro(p, xp, e, buff);
	if ((tmp_buff - buff) != (long int)field->size) {
	  error("The field %s wrote an unexpected number of bits", field->name);
	}
      }
      /* Write it manually */
      else {
	if (field->use_xpart == 1)
	  memcpy(buff, ((char *)xp) + field->offset, field->size);
	else if (field->use_xpart == 0)
	  memcpy(buff, ((char *)p) + field->offset, field->size);
	else
	  error(
		"It seems that you are using the wrong CSDS function in the hydro."
		" You need to use csds_define_hydro_standard_field and not"
		" the general one.");
      }

      buff += field->size;
    }
    mask &= ~type_mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in parts. Remaining mask: %u", mask);
  }
#endif
}

/**
 * @brief Dump a #part to the log.
 *
 * @param log The #csds_writer
 * @param p The #part to dump.
 * @param xp The #xpart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void csds_log_part(struct csds_writer *log, const struct part *p,
		   struct xpart *xp, const struct engine *e,
		   const int log_all_fields, const enum csds_special_flags flag,
		   const int flag_data) {

  csds_log_parts(log, p, xp, /* count= */ 1, e, log_all_fields, flag,
		 flag_data);
}

/**
 * @brief Compute the total size and mask of a complete particle record.
 *
 * This function calculates the size and mask by combining:
 * 1. Fixed fields (Pos, Vel, Accel).
 * 2. Particle-type specific fields (the group fields).
 * 3. The Special Flag (if requested).
 * 4. The Record Header.
 *
 * @param log The #csds_writer.
 * @param type The swift particle type index (e.g., swift_type_gas).
 * @param flag The value of the special flags.
 * @param size (output) The total size of the record in bytes.
 * @param mask (output) The total mask to use in the header.
 */
void csds_compute_total_record_size_and_mask(const struct csds_writer *log,
					     const int type,
					     const enum csds_special_flags flag,
					     size_t *size,
					     unsigned int *mask) {
    *size = 0;
    *mask = 0;

    /* Add Fixed Fields (assuming we log them all for a particle record) */
    *mask |= log->fixed_fields[CSDS_POS_INDEX].mask;
    *size += log->fixed_fields[CSDS_POS_INDEX].size;
    *mask |= log->fixed_fields[CSDS_VEL_INDEX].mask;
    *size += log->fixed_fields[CSDS_VEL_INDEX].size;
    *mask |= log->fixed_fields[CSDS_ACCEL_INDEX].mask;
    *size += log->fixed_fields[CSDS_ACCEL_INDEX].size;

    /* Add Particle-Type-Specific Fields (using pre-calculated totals) */
    if (log->number_fields[type] > 0) {
	*mask |= log->part_type_masks[type];
	*size += log->part_type_total_size[type];
    }

    /* Add Special Flag */
    if (flag != csds_flag_none) {
	*size += log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size;
	*mask |= log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask;
    }

    /* Add Header Size */
    *size += CSDS_HEADER_SIZE;
}

/**
 * @brief Dump a group of #part to the log.
 *
 * @param log The #csds_writer.
 * @param p The #part to dump.
 * @param xp The #xpart to dump.
 * @param count The number of particle to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void csds_log_parts(struct csds_writer *log, const struct part *p,
		    struct xpart *xp, int count, const struct engine *e,
		    const int log_all_fields,
		    const enum csds_special_flags flag, const int flag_data) {

  /* Build the special flag */
  const uint32_t special_flags =
      csds_pack_flags_and_data(flag, flag_data, swift_type_gas);

  /* Compute the size of the buffer and the mask. */
  size_t size = 0;
  unsigned int mask = 0;
  csds_compute_total_record_size_and_mask(log, swift_type_gas, flag, &size, &mask);
  size_t size_total = count * size;

  /* Allocate a chunk of memory in the logfile of the right size. */
  size_t offset_new;
  char *buff =
      (char *)csds_logfile_writer_get(&log->logfile, size_total, &offset_new);

#ifdef SWIFT_DEBUG_CHECKS
  /* Save the buffer position in order to test if the requested buffer was
   * really used */
  const char *buff_before = buff;
#endif

  /* Write the particles */
  for (int i = 0; i < count; i++) {
    /* reset the offset of the previous log */
    if (flag == csds_flag_create || flag == csds_flag_mpi_enter) {
      xp[i].csds_data.last_offset = 0;
    }

    /* Copy everything into the buffer */
    csds_copy_part_fields(log, &p[i], &xp[i], e, mask,
			  &xp[i].csds_data.last_offset, offset_new, buff,
			  special_flags);

    /* Update the pointers */
    xp[i].csds_data.last_offset = offset_new;
    xp[i].csds_data.steps_since_last_output = 0;
    buff += size;
    offset_new += size;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure that the buffer was fully used */
  const int diff = buff - buff_before;
  if (diff != (int)size_total) {
    error("The requested buffer was not totally used: %i != %zi", diff,
	  size_total);
  }
#endif
}

/**
 * @brief Copy the particle fields into a given buffer.
 *
 * @param log The #csds_writer.
 * @param sp The #spart to copy.
 * @param e The #engine.
 * @param mask The mask for the fields to write.
 * @param offset The offset to the previous log.
 * @param offset_new The offset of the current record.
 * @param buff The buffer to use when writing.
 * @param special_flags The data for the special flags.
 */
void csds_copy_spart_fields(const struct csds_writer *log,
			    const struct spart *sp, const struct engine *e,
			    unsigned int mask, size_t *offset,
			    size_t offset_new, char *buff,
			    const uint32_t special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mask == 0) {
    error("You should always log at least one field.");
  }
#endif

  /* Write the header. */
  buff = csds_write_record_header(buff, &mask, offset, offset_new);

  /* --- Fixed Fields --- */
  /* Special flags (CSDS_SPECIAL_FLAGS_INDEX = 0) */
  if (mask & log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask) {
    memcpy(buff, &special_flags, log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size);
    buff += log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask;
  }

  /* Positions (CSDS_POS_INDEX = 2) */
  if (mask & log->fixed_fields[CSDS_POS_INDEX].mask) {
    memcpy(buff, ((char *)sp) + offsetof(struct spart, x), log->fixed_fields[CSDS_POS_INDEX].size);
    buff += log->fixed_fields[CSDS_POS_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_POS_INDEX].mask;
  }

  /* Velocities (CSDS_VEL_INDEX = 3) */
  if (mask & log->fixed_fields[CSDS_VEL_INDEX].mask) {
    memcpy(buff, ((char *)sp) + offsetof(struct spart, v), log->fixed_fields[CSDS_VEL_INDEX].size);
    buff += log->fixed_fields[CSDS_VEL_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_VEL_INDEX].mask;
  }

  /* Accelerations (CSDS_ACCEL_INDEX = 4) */
  if (mask & log->fixed_fields[CSDS_ACCEL_INDEX].mask) {
    // Accelerations requires a conversion function for stars
    csds_stars_convert_acc(sp, e, buff);
    buff += log->fixed_fields[CSDS_ACCEL_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_ACCEL_INDEX].mask;
  }

  /* --- Particle-Type-Specific Fields (Stars) --- */

  const unsigned int type_mask = log->part_type_masks[swift_type_stars];
  if (mask & type_mask) {
    for (int i = 0; i < log->number_fields[swift_type_stars]; i++) {
      struct csds_field *field = log->part_type_fields[swift_type_stars] + i;

      /* Do we have a conversion function? */
      if (field->conversion_stars) {
	char *tmp_buff = field->conversion_stars(sp, e, buff);
	if ((tmp_buff - buff) != (long int)field->size) {
	  error("The field %s wrote an unexpected number of bits", field->name);
	}
      }
      /* Write it manually */
      else {
	memcpy(buff, ((char *)sp) + field->offset, field->size);
      }

      buff += field->size;
    }
    mask &= ~type_mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in sparts. Remaining mask: %u", mask);
  }
#endif
}

/**
 * @brief Dump a #spart to the log.
 *
 * @param log The #csds_writer
 * @param sp The #spart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void csds_log_spart(struct csds_writer *log, struct spart *sp,
		    const struct engine *e, const int log_all_fields,
		    const enum csds_special_flags flag, const int flag_data) {

  csds_log_sparts(log, sp, /* count */ 1, e, log_all_fields, flag, flag_data);
}

/**
 * @brief Dump a group of #spart to the log.
 *
 * @param log The #csds_writer
 * @param sp The #spart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param count The number of particle to dump.
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void csds_log_sparts(struct csds_writer *log, struct spart *sp, int count,
		     const struct engine *e, const int log_all_fields,
		     const enum csds_special_flags flag, const int flag_data) {
  /* Build the special flag */
  const uint32_t special_flags =
      csds_pack_flags_and_data(flag, flag_data, swift_type_stars);

  /* Compute the size of the buffer and the mask. */
  unsigned int mask = 0;
  size_t size = 0;
  csds_compute_total_record_size_and_mask(log, swift_type_stars, flag, &size, &mask);
  size_t size_total = count * size;

  /* Allocate a chunk of memory in the logfile of the right size. */
  size_t offset_new;
  char *buff =
      (char *)csds_logfile_writer_get(&log->logfile, size_total, &offset_new);
#ifdef SWIFT_DEBUG_CHECKS
  /* Save the buffer position in order to test if the requested buffer was
   * really used */
  const char *buff_before = buff;
#endif

  for (int i = 0; i < count; i++) {
    /* reset the offset of the previous log */
    if (flag == csds_flag_create || flag == csds_flag_mpi_enter) {
      sp[i].csds_data.last_offset = 0;
    }

    /* Copy everything into the buffer */
    csds_copy_spart_fields(log, &sp[i], e, mask, &sp[i].csds_data.last_offset,
			   offset_new, buff, special_flags);

    /* Update the pointers */
    sp[i].csds_data.last_offset = offset_new;
    sp[i].csds_data.steps_since_last_output = 0;
    buff += size;
    offset_new += size;
  }
#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure that the buffer was fully used */
  const int diff = buff - buff_before;
  if (diff != (int)size_total) {
    error("It seems that the requested buffer was not totally used: %i != %zi",
	  diff, size_total);
  }
#endif
}

/**
 * @brief Copy the particle fields into a given buffer.
 *
 * @param log The #csds_writer.
 * @param gp The #gpart to copy.
 * @param e The #engine.
 * @param mask The mask for the fields to write.
 * @param offset The offset to the previous log.
 * @param offset_new The offset of the current record.
 * @param buff The buffer to use when writing.
 * @param special_flags The data of the special flag.
 */
void csds_copy_gpart_fields(const struct csds_writer *log,
			    const struct gpart *gp, const struct engine *e,
			    unsigned int mask, size_t *offset,
			    size_t offset_new, char *buff,
			    const uint32_t special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mask == 0) {
    error("You should always log at least one field.");
  }
#endif

  /* Write the header. */
  buff = csds_write_record_header(buff, &mask, offset, offset_new);

  /* --- Fixed Fields --- */
  /* Special flags (CSDS_SPECIAL_FLAGS_INDEX = 0) */
  if (mask & log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask) {
    memcpy(buff, &special_flags, log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size);
    buff += log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_SPECIAL_FLAGS_INDEX].mask;
  }

  /* Positions (CSDS_POS_INDEX = 2) */
  if (mask & log->fixed_fields[CSDS_POS_INDEX].mask) {
    memcpy(buff, ((char *)gp) + offsetof(struct gpart, x), log->fixed_fields[CSDS_POS_INDEX].size);
    buff += log->fixed_fields[CSDS_POS_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_POS_INDEX].mask;
  }

  /* Velocities (CSDS_VEL_INDEX = 3) */
  if (mask & log->fixed_fields[CSDS_VEL_INDEX].mask) {
    memcpy(buff, ((char *)gp) + offsetof(struct gpart, v_full), log->fixed_fields[CSDS_VEL_INDEX].size);
    buff += log->fixed_fields[CSDS_VEL_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_VEL_INDEX].mask;
  }

  /* Accelerations (CSDS_ACCEL_INDEX = 4) */
  if (mask & log->fixed_fields[CSDS_ACCEL_INDEX].mask) {
    /* Accelerations requires a conversion function for gravity */
    csds_gravity_convert_acc(gp, e, buff);
    buff += log->fixed_fields[CSDS_ACCEL_INDEX].size;
    mask &= ~log->fixed_fields[CSDS_ACCEL_INDEX].mask;
  }

  /* --- Particle-Type-Specific Fields (Gravity/DM) --- */

  const unsigned int type_mask = log->part_type_masks[gp->type];
  if (mask & type_mask) {
    for (int i = 0; i < log->number_fields[gp->type]; i++) {
      struct csds_field *field = log->part_type_fields[gp->type] + i;

      /* Do we have a conversion function? */
      if (field->conversion_grav) {
	char *tmp_buff = field->conversion_grav(gp, e, buff);
	if ((tmp_buff - buff) != (long int)field->size) {
	  error("The field %s wrote an unexpected number of bits", field->name);
	}
      }
      /* Write it manually */
      else {
	memcpy(buff, ((char *)gp) + field->offset, field->size);
      }

      buff += field->size;
    }
    mask &= ~type_mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in gparts. Remaining mask: %u", mask);
  }
#endif
}

/**
 * @brief Dump a #gpart to the log.
 *
 * @param log The #csds_writer
 * @param p The #gpart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void csds_log_gpart(struct csds_writer *log, struct gpart *p,
		    const struct engine *e, const int log_all_fields,
		    const enum csds_special_flags flag, const int flag_data) {
  csds_log_gparts(log, p, /* count */ 1, e, log_all_fields, flag, flag_data);
}

/**
 * @brief Dump a group of #gpart to the log.
 *
 * @param log The #csds_writer
 * @param p The #gpart to dump.
 * @param count The number of particle to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void csds_log_gparts(struct csds_writer *log, struct gpart *p, int count,
		     const struct engine *e, const int log_all_fields,
		     const enum csds_special_flags flag, const int flag_data) {
  /* Build the special flag */
  const uint32_t special_flags =
      csds_pack_flags_and_data(flag, flag_data, swift_type_dark_matter);

  /* Compute the size of the buffer. */
  /* As we might have some non DM particles, we cannot log_all_fields blindly */
  int count_dm = 0;
  for (int i = 0; i < count; i++) {
    /* Log only the dark matter */
    if (p[i].type != swift_type_dark_matter &&
	p[i].type != swift_type_dark_matter_background)
      continue;

    count_dm += 1;
  }
  unsigned int mask = 0;
  size_t size = 0;
  csds_compute_total_record_size_and_mask(log, swift_type_dark_matter, flag, &size, &mask);
  size_t size_total = size * count_dm;

  /* Allocate a chunk of memory in the logfile of the right size. */
  size_t offset_new;
  char *buff =
      (char *)csds_logfile_writer_get(&log->logfile, size_total, &offset_new);
#ifdef SWIFT_DEBUG_CHECKS
  /* Save the buffer position in order to test if the requested buffer was
   * really used */
  const char *buff_before = buff;
#endif

  for (int i = 0; i < count; i++) {
    /* Log only the dark matter */
    if (p[i].type != swift_type_dark_matter &&
	p[i].type != swift_type_dark_matter_background)
      continue;

    /* reset the offset of the previous log */
    if (flag == csds_flag_create || flag == csds_flag_mpi_enter) {
      p[i].csds_data.last_offset = 0;
    }

    /* Copy everything into the buffer */
    csds_copy_gpart_fields(log, &p[i], e, mask, &p[i].csds_data.last_offset,
			   offset_new, buff, special_flags);

    /* Update the pointers */
    p[i].csds_data.last_offset = offset_new;
    p[i].csds_data.steps_since_last_output = 0;
    buff += size;
    offset_new += size;
  }
#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure that the buffer was fully used */
  const int diff = buff - buff_before;
  if (diff != (int)size_total) {
    error("It seems that the requested buffer was not totally used: %i != %zi",
	  diff, size_total);
  }
#endif
}

/**
 * @brief write a timestamp
 *
 * @param log The #csds_writer
 * @param timestamp time to write
 * @param time time or scale factor
 * @param offset Pointer to the offset of the previous log of this particle;
 * (return) offset of this log.
 */
void csds_log_timestamp(struct csds_writer *log, integertime_t timestamp,
                        double time, size_t *offset) {
  struct csds_logfile_writer *logfile = &log->logfile;
  /* Start by computing the size of the message. */
  const int size =
      log->list_fields[CSDS_TIMESTAMP_INDEX].size + CSDS_HEADER_SIZE;

  /* Allocate a chunk of memory in the logfile of the right size. */
  size_t offset_new;
  char *buff = (char *)csds_logfile_writer_get(logfile, size, &offset_new);

  /* Write the header. */
  unsigned int mask = log->list_fields[CSDS_TIMESTAMP_INDEX].mask;
  buff = csds_write_record_header(buff, &mask, offset, offset_new);

  /* Store the timestamp. */
  memcpy(buff, &timestamp, sizeof(integertime_t));
  buff += sizeof(integertime_t);

  /* Store the time. */
  memcpy(buff, &time, sizeof(double));

  /* Update the log message offset. */
  *offset = offset_new;
}

/**
 * @brief Ensure that the buffer is large enough for a step.
 *
 * Check if csds parameters are large enough to write all particles
 * and ensure that enough space is available in the buffer.
 *
 * @param log The #csds_writer
 * @param e The #engine.
 */
void csds_ensure_size(struct csds_writer *log, const struct engine *e) {

  ticks tic = getticks();

  /* count part memory */
  size_t limit = 0;

  /* count part memory */
  limit += e->s->nr_parts;

  /* count gpart memory */
  limit += e->s->nr_gparts;

  /* count spart memory. */
  limit += e->s->nr_sparts;

  // TODO improve estimate with the size of each particle
  limit *= log->max_record_size;

  /* ensure enough space in logfile */
  csds_logfile_writer_ensure(&log->logfile, limit, log->buffer_scale * limit);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/** @brief Generate the name of the logfile files
 *
 * @param log The #csds_writer.
 * @param filename The filename of the logfile file.
 */
void csds_get_logfile_name(struct csds_writer *log, char *filename) {
  sprintf(filename, "%s_%04i.dump", log->base_name, engine_rank);
}

/**
 * @brief Initialize the variable list_fields.
 *
 * @param log The #csds_writer.
 * @param e The #engine.
 */
void csds_init_masks(struct csds_writer *log, const struct engine *e) {
  /* Set the pointers to 0 */
  for (int i = 0; i < swift_type_count; i++) {
    log->field_pointers[i] = NULL;
    log->number_fields[i] = 0;
  }

  struct csds_field list[100];
  int num_fields = 0;

  /* The next fields must be the two first ones. */
  /* Add the special flags (written manually => no need of offset) */
  if (CSDS_SPECIAL_FLAGS_INDEX != 0) {
    error("Expecting the special flags to be the first element.");
  }
  csds_define_common_field(list[CSDS_SPECIAL_FLAGS_INDEX], "SpecialFlags",
                           sizeof(uint32_t));
  num_fields += 1;

  /* Add the timestamp */
  if (CSDS_TIMESTAMP_INDEX != 1) {
    error("Expecting the timestamp to be the first element.");
  }
  csds_define_common_field(list[CSDS_TIMESTAMP_INDEX], "Timestamp",
                           sizeof(integertime_t) + sizeof(double));
  list[num_fields].type = mask_for_timestep;  // flag it as timestamp
  num_fields += 1;

  /* Initialize all the particles types */
  for (int i = 0; i < swift_type_count; i++) {
    int tmp_num_fields = 0;
    struct csds_field *current = &list[num_fields];
    enum mask_for_type mask_for_type;

    /* Set the pointer */
    log->field_pointers[i] = current;

    switch (i) {
      /* Hydro */
      case swift_type_gas:
        /* Set the mask type */
        mask_for_type = mask_for_gas;

        /* Set the masks */
        tmp_num_fields = csds_hydro_define_fields(current);
        tmp_num_fields +=
            csds_chemistry_define_fields_parts(current + tmp_num_fields);
        // TODO add cooling + SF
        break;

      /* Stars */
      case swift_type_stars:
        /* Set the mask type */
        mask_for_type = mask_for_stars;

        /* Set the masks */
        tmp_num_fields = csds_stars_define_fields(current);
        tmp_num_fields +=
            csds_chemistry_define_fields_sparts(current + tmp_num_fields);
        tmp_num_fields +=
            csds_star_formation_define_fields(current + tmp_num_fields);
        break;

      case swift_type_dark_matter:
        /* Set the mask type */
        mask_for_type = mask_for_dark_matter;

        /* Set the masks */
        tmp_num_fields = csds_gravity_define_fields(current);
        break;

      default:
        log->field_pointers[i] = NULL;
        break;
    }

    /* Set the particle type */
    for (int j = 0; j < tmp_num_fields; j++) {
      current[j].type = mask_for_type;
    }

    /* Update the number of fields */
    num_fields += tmp_num_fields;
    log->number_fields[i] = tmp_num_fields;
  }

  /* Set the counter */
  log->total_number_fields = num_fields;

  /* Set the masks and ensure to have only one for the common fields
     (e.g. Coordinates).
     Initially we have (Name, mask, part_type):
      - Coordinates, 0, 0
      - Velocity, 0, 0
      - Coordinates, 0, 1

      And get:
      - Coordinates, 1, 0
      - Velocity, 2, 0
      - Coordinates, 1, 1
  */
  int mask = 0;
  for (int i = 0; i < num_fields; i++) {
    /* Skip the elements already processed. */
    if (list[i].mask != 0) {
      continue;
    }
    const char *name = list[i].name;
    list[i].mask = 1 << mask;

    /* Check if the field exists in the other particle type. */
    for (int j = i + 1; j < num_fields; j++) {
      /* Check if the name is the same */
      if (strcmp(name, list[j].name) == 0) {
        /* Check if the data are the same */
        if (list[i].size != list[j].size) {
          error("Found two same fields but with different data size (%s).",
                name);
        }

        list[j].mask = 1 << mask;
      }
    }
    mask += 1;
  }

  /* Check that we have enough available flags. */
  if (mask >= 8 * CSDS_MASK_SIZE) {
    error(
        "Not enough available flags for all the fields. "
        "Please reduce the number of output fields.");
  }

  /* Save the data */
  size_t size_list = sizeof(struct csds_field) * num_fields;
  log->list_fields = (struct csds_field *)malloc(size_list);
  memcpy(log->list_fields, list, size_list);

  /* Update the pointers */
  for (int i = 0; i < swift_type_count; i++) {
    if (log->field_pointers[i] != NULL) {
      log->field_pointers[i] =
          log->list_fields + (log->field_pointers[i] - list);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) {
    message("The CSDS contains the following masks:");
    for (int i = 0; i < log->total_number_fields; i++) {
      message("%20s:\t mask=%03u\t size=%zi", log->list_fields[i].name,
              log->list_fields[i].mask, log->list_fields[i].size);
    }
  }
#endif
}

/**
 * @brief intialize the csds structure
 *
 * @param log The #csds_writer
 * @param e The #engine.
 * @param params The #swift_params
 */
void csds_init(struct csds_writer *log, const struct engine *e,
               struct swift_params *params) {

  ticks tic = getticks();

#ifdef WITH_MPI
  /* Should be safe, but better to check */
  if (e->nr_nodes >= 1 << 16)
    error(
        "The special flag does not contain enough bits"
        "to store the information about the ranks.");
#endif

  /* read parameters. */
  log->delta_step = parser_get_param_int(params, "CSDS:delta_step");
  size_t buffer_size =
      parser_get_opt_param_float(params, "CSDS:initial_buffer_size", 0.5) * 1e9;
  log->buffer_scale =
      parser_get_opt_param_float(params, "CSDS:buffer_scale", 10);
  parser_get_param_string(params, "CSDS:basename", log->base_name);

  /* Initialize the list_fields */
  csds_init_masks(log, e);

  /* set initial value of parameters. */
  log->timestamp_offset = 0;

  /* generate logfile filename. */
  char csds_name_file[PARSER_MAX_LINE_SIZE];
  csds_get_logfile_name(log, csds_name_file);

  /* Compute max size for a particle record. */
  int max_size = CSDS_OFFSET_SIZE + CSDS_MASK_SIZE;

  /* Loop over all fields except timestamp. */
  for (int i = 0; i < log->total_number_fields; i++) {
    /* Skip the timestamp */
    if (i == CSDS_TIMESTAMP_INDEX) continue;

    max_size += log->list_fields[i].size;
  }
  log->max_record_size = max_size;

  /* init logfile. */
  csds_logfile_writer_init(&log->logfile, csds_name_file, buffer_size);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Close logfile file and desallocate memory
 *
 * @param log The #csds_writer
 */
void csds_free(struct csds_writer *log) {
  csds_logfile_writer_close(&log->logfile);

  free(log->list_fields);
  log->list_fields = NULL;
  log->total_number_fields = 0;
}

/**
 * @brief Write a file header to a csds file
 *
 * @param log The #csds_writer
 *
 */
void csds_write_file_header(struct csds_writer *log) {

  /* get required variables. */
  struct csds_logfile_writer *logfile = &log->logfile;

  /* Write the beginning of the header */
  char *offset_first_record =
      csds_logfile_writer_write_begining_header(logfile);

  /* placeholder to write the number of unique masks. */
  size_t file_offset = 0;
  char *skip_unique_masks =
      csds_logfile_writer_get(logfile, sizeof(unsigned int), &file_offset);

  /* write masks. */
  // loop over all mask type.
  unsigned int unique_mask = 0;
  for (int i = 0; i < log->total_number_fields; i++) {
    /* Check if the mask was not already written */
    int is_written = 0;
    for (int j = 0; j < i; j++) {
      if (log->list_fields[i].mask == log->list_fields[j].mask) {
        is_written = 1;
        break;
      }
    }

    if (is_written) {
      continue;
    }

    unique_mask += 1;

    // mask name.
    csds_write_data(logfile, &file_offset, CSDS_STRING_SIZE,
                    (const char *)&log->list_fields[i].name);

    // mask size.
    csds_write_data(logfile, &file_offset, sizeof(unsigned int),
                    (const char *)&log->list_fields[i].size);
  }
  memcpy(skip_unique_masks, &unique_mask, sizeof(unsigned int));

  /* Write the number of fields per particle */
  csds_write_data(logfile, &file_offset, sizeof(log->number_fields),
                  (const char *)log->number_fields);

  /* Now write the order for each particle type */
  for (int i = 0; i < swift_type_count; i++) {
    int number_fields = log->number_fields[i];
    if (number_fields == 0) continue;

    struct csds_field *field = log->field_pointers[i];

    /* Loop over all the fields from this particle type */
    for (int k = 0; k < number_fields; k++) {
      unsigned int current_mask = field[k].mask;

      /* Find the index among all the fields. */
      for (int m = 0; m < log->total_number_fields; m++) {
        if (log->list_fields[m].mask == current_mask) {
          csds_write_data(logfile, &file_offset, sizeof(int), (const char *)&m);
          break;
        }
      }
    }
  }

  /* Write the end of the header. */
  csds_logfile_writer_write_end_header(logfile, offset_first_record);
}

/**
 * @brief read record header
 *
 * @param buff The reading buffer
 * @param mask The mask to read
 * @param offset (return) the offset pointed by this record (absolute)
 * @param cur_offset The current record offset
 *
 * @return Number of bytes read
 */
__attribute__((always_inline)) INLINE static int csds_read_record_header(
    const char *buff, unsigned int *mask, size_t *offset, size_t cur_offset) {
  memcpy(mask, buff, CSDS_MASK_SIZE);
  buff += CSDS_MASK_SIZE;

  *offset = 0;
  memcpy(offset, buff, CSDS_OFFSET_SIZE);
  *offset = cur_offset - *offset;

  return CSDS_MASK_SIZE + CSDS_OFFSET_SIZE;
}

/**
 * @brief Read a csds message and store the data in a #part.
 *
 * @param p The #part in which to store the values.
 * @param offset Pointer to the offset of the csds message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded csds message.
 *
 * @return The mask containing the values read.
 */
int csds_read_part(const struct csds_writer *log, struct part *p,
                   size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the csds mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += csds_read_record_header(buff, &mask, offset, cur_offset);

  for (int i = 0; i < log->total_number_fields; i++) {
    if ((mask & log->list_fields[i].mask) &&
        (log->list_fields[i].type == mask_for_gas)) {

      const char *name = log->list_fields[i].name;
      if (strcmp("Coordinates", name) == 0) {
        memcpy(p->x, buff, 3 * sizeof(double));
        buff += 3 * sizeof(double);
      } else if (strcmp("Velocities", name) == 0) {
        memcpy(p->v, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("Accelerations", name) == 0) {
        memcpy(p->a_hydro, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("SmoothingLengths", name) == 0) {
        memcpy(&p->h, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("InternalEnergies", name) == 0) {
        memcpy(&p->u, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("Masses", name) == 0) {
        memcpy(&p->mass, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("Densities", name) == 0) {
        memcpy(&p->rho, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("ParticleIDs", name) == 0) {
        memcpy(&p->id, buff, sizeof(long long));
        buff += sizeof(long long);
      } else if (strcmp("SPHENIXSecondaryFields", name) == 0) {
        // No need to read it for testing
        buff += 7 * sizeof(float);
      } else {
        error("Field '%s' not found", name);
      }
    }
  }

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Read a csds message and store the data in a #gpart.
 *
 * @param p The #gpart in which to store the values.
 * @param offset Pointer to the offset of the csds message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded csds message.
 *
 * @return The mask containing the values read.
 */
int csds_read_gpart(const struct csds_writer *log, struct gpart *p,
                    size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the csds mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += csds_read_record_header(buff, &mask, offset, cur_offset);

  for (int i = 0; i < log->total_number_fields; i++) {
    if ((mask & log->list_fields[i].mask) &&
        (log->list_fields[i].type == mask_for_dark_matter)) {

      const char *name = log->list_fields[i].name;
      if (strcmp("Coordinates", name) == 0) {
        memcpy(p->x, buff, 3 * sizeof(double));
        buff += 3 * sizeof(double);
      } else if (strcmp("Velocities", name) == 0) {
        memcpy(p->v_full, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("Accelerations", name) == 0) {
        memcpy(p->a_grav, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("ParticleIDs", name) == 0) {
        memcpy(&p->id_or_neg_offset, buff, sizeof(long long));
        buff += sizeof(long long);
      } else if (strcmp("Masses", name) == 0) {
        memcpy(&p->mass, buff, sizeof(float));
        buff += sizeof(float);
      } else {
        error("Field '%s' not found", name);
      }
    }
  }

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Read a csds message for a timestamp.
 *
 * @param log The #csds_writer.
 * @param t The timestamp in which to store the value.
 * @param time The time in which to store the value.
 * @param offset Pointer to the offset of the csds message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded csds message.
 *
 * @return The mask containing the values read.
 */
int csds_read_timestamp(const struct csds_writer *log, integertime_t *t,
                        double *time, size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the csds mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += csds_read_record_header(buff, &mask, offset, cur_offset);

  /* We are only interested in timestamps. */
  if (!(mask & log->list_fields[CSDS_TIMESTAMP_INDEX].mask))
    error("Trying to read timestamp from a particle.");

  /* Make sure we don't have extra fields. */
  if (mask != log->list_fields[CSDS_TIMESTAMP_INDEX].mask)
    error("Timestamp message contains extra fields.");

  /* Copy the timestamp value from the buffer. */
  memcpy(t, buff, sizeof(integertime_t));
  buff += sizeof(integertime_t);

  /* Copy the timestamp value from the buffer. */
  memcpy(time, buff, sizeof(double));

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Write a swift_params struct to the given FILE as a stream of bytes.
 *
 * @param log the struct
 * @param stream the file stream
 */
void csds_struct_dump(const struct csds_writer *log, FILE *stream) {
  restart_write_blocks((void *)log, sizeof(struct csds_writer), 1, stream,
                       "csds", "csds");

  /* Write the masks */
  restart_write_blocks((void *)log->list_fields, sizeof(struct csds_field),
                       log->total_number_fields, stream, "csds_masks",
                       "csds_masks");
}

/**
 * @brief Restore a csds struct from the given FILE as a stream of
 * bytes.
 *
 * @param csds the struct
 * @param stream the file stream
 */
void csds_struct_restore(struct csds_writer *log, FILE *stream) {
  /* Read the block */
  restart_read_blocks((void *)log, sizeof(struct csds_writer), 1, stream, NULL,
                      "csds");

  /* Read the masks */
  const struct csds_field *old_list_fields = log->list_fields;
  log->list_fields = (struct csds_field *)malloc(sizeof(struct csds_field) *
                                                 log->total_number_fields);

  restart_read_blocks((void *)log->list_fields, sizeof(struct csds_field),
                      log->total_number_fields, stream, NULL, "csds_masks");

  /* Restore the pointers */
  for (int i = 0; i < swift_type_count; i++) {
    if (log->field_pointers[i] == NULL) continue;

    log->field_pointers[i] =
        log->list_fields + (log->field_pointers[i] - old_list_fields);
  }

  /* Restart the logfile. */
  char csds_name_file[PARSER_MAX_LINE_SIZE];
  csds_get_logfile_name(log, csds_name_file);

  csds_logfile_writer_restart(&log->logfile, csds_name_file);
}

#endif /* WITH_CSDS */

#endif /* HAVE_POSIX_FALLOCATE */
