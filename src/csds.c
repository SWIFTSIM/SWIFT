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
#include "../config.h"

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
#include "dump.h"
#include "engine.h"
#include "error.h"
#include "gravity_csds.h"
#include "hydro_csds.h"
#include "star_formation_csds.h"
#include "stars_csds.h"
#include "units.h"

/*
 * Thoses are definitions from the format and therefore should not be changed!
 */
/* Number of bytes for a mask. */
// TODO change this to number of bits
#define csds_mask_size 2

/* Number of bits for record header. */
#define csds_header_bytes 8

/* Number bytes for an offset. */
#define csds_offset_size csds_header_bytes - csds_mask_size

/* Number of bytes for the file format information. */
#define csds_format_size 20

char csds_file_format[csds_format_size] = "SWIFT_CSDS";

/*
 * The two following defines need to correspond to the list's order
 * in csds_init_masks.
 */
/* Index of the special flags in the list of masks */
#define csds_index_special_flags 0
/* Index of the timestamp in the list of masks */
#define csds_index_timestamp 1

/**
 * @brief Print the current size used by the logger in GB (not the allocated
 * one).
 *
 * @param log The #csds_writer.
 * @param e The #engine.
 */
float csds_get_current_filesize_used_gb(const struct csds_writer *log,
                                        const struct engine *e) {
  return log->dump.count / (1024.f * 1024.f * 1024.f);
}

/**
 * @brief Write the header of a record (offset + mask).
 *
 * This is maybe broken for big(?) endian.
 *
 * @param buff The buffer where to write the mask and offset.
 * @param mask The mask to write inside the buffer.
 * @param offset The offset of the previous record.
 * @param offset_new The offset of the current record.
 *
 * @return updated buff
 */
char *csds_write_record_header(char *buff, const unsigned int *mask,
                               const size_t *offset, const size_t offset_new) {
  /* write mask. */
  memcpy(buff, mask, csds_mask_size);
  buff += csds_mask_size;

  /* write offset. */
  uint64_t diff_offset = offset_new - *offset;
  memcpy(buff, &diff_offset, csds_offset_size);
  buff += csds_offset_size;

  return buff;
}

/**
 * @brief Write to the dump.
 *
 * @param d #dump file
 * @param offset (return) offset of the data
 * @param size number of bytes to write
 * @param p pointer to the data
 */
void csds_write_data(struct dump *d, size_t *offset, size_t size,
                     const void *p) {
  /* get buffer. */
  char *buff = dump_get(d, size, offset);

  /* write data to the buffer. */
  memcpy(buff, p, size);

  /* Update offset to end of record. */
  *offset += size;
}

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
 * @param first_log Is it the first log of the particles?
 */
void csds_log_all_particles(struct csds_writer *log, const struct engine *e,
                            int first_log) {

  /* Ensure that enough space is available. */
  csds_ensure_size(log, e);

  /* csds_ensure_size is tracked separately. */
  ticks tic = getticks();

  /* some constants. */
  const struct space *s = e->s;

  /* Create the flags */
  const enum csds_special_flags flag =
      first_log ? csds_flag_create : csds_flag_none;

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

  /* Special flags */
  if (mask & log->list_fields[csds_index_special_flags].mask) {
    memcpy(buff, &special_flags,
           log->list_fields[csds_index_special_flags].size);
    buff += log->list_fields[csds_index_special_flags].size;
    mask &= ~log->list_fields[csds_index_special_flags].mask;
  }

  /* Write the hydro fields */
  for (int i = 0; i < log->number_fields[swift_type_gas]; i++) {
    struct csds_field *field = &log->field_pointers[swift_type_gas][i];

    /* Skip the fields that are not required. */
    if (!(mask & field->mask)) continue;

    /* Do we have a conversion function? */
    if (field->conversion_hydro) {
      char *tmp_buff = field->conversion_hydro(p, xp, e, buff);
      /* Check that the correct number of bits are written */
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

    /* Update the variables */
    buff += field->size;
    mask &= ~field->mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in parts. %u", mask);
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
 * @brief Compute the size and the mask of all the fields that will be written.
 *
 * @param fields The list of fields to write.
 * @param n_fields The number of fields to write.
 * @param size (output) The size of all the fields.
 * @param mask (output) The mask to use.
 */
void csds_compute_size_and_mask(struct csds_field *fields, int n_fields,
                                size_t *size, unsigned int *mask) {
  *size = 0;
  *mask = 0;
  // TODO: write only some fields
  for (int i = 0; i < n_fields; i++) {
    *size += fields[i].size;
    *mask |= fields[i].mask;
  }
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
  const int size_special_flag = log->list_fields[csds_index_special_flags].size;
  const uint32_t special_flags =
      csds_pack_flags_and_data(flag, flag_data, swift_type_gas);

  /* Compute the size of the buffer. */
  size_t size = 0;
  unsigned int mask = 0;
  // TODO: write only some fields
  csds_compute_size_and_mask(log->field_pointers[swift_type_gas],
                             log->number_fields[swift_type_gas], &size, &mask);
  /* Add the flag */
  if (flag != csds_flag_none) {
    size += size_special_flag;
    mask |= log->list_fields[csds_index_special_flags].mask;
  }
  size += csds_header_bytes;
  size_t size_total = count * size;

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size_total, &offset_new);

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

  /* Special flags */
  if (mask & log->list_fields[csds_index_special_flags].mask) {
    memcpy(buff, &special_flags,
           log->list_fields[csds_index_special_flags].size);
    buff += log->list_fields[csds_index_special_flags].size;
    mask &= ~log->list_fields[csds_index_special_flags].mask;
  }

  /* Write the stellar fields */
  for (int i = 0; i < log->number_fields[swift_type_stars]; i++) {
    struct csds_field *field = &log->field_pointers[swift_type_stars][i];

    /* Skip the fields that are not required. */
    if (!(mask & field->mask)) continue;

    /* Do we have a conversion function? */
    if (field->conversion_stars) {
      char *tmp_buff = field->conversion_stars(sp, e, buff);
      /* Check that the correct number of bits are written */
      if ((tmp_buff - buff) != (long int)field->size) {
        error("The field %s wrote an unexpected number of bits", field->name);
      }
    }
    /* Write it manually */
    else {
      memcpy(buff, ((char *)sp) + field->offset, field->size);
    }

    /* Update the variables */
    buff += field->size;
    mask &= ~field->mask;
  }
#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in sparts. %u", mask);
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
  const int size_special_flag = log->list_fields[csds_index_special_flags].size;
  const uint32_t special_flags =
      csds_pack_flags_and_data(flag, flag_data, swift_type_stars);

  /* Compute the size of the buffer. */
  // TODO: write only some fields
  unsigned int mask = 0;
  size_t size = 0;
  csds_compute_size_and_mask(log->field_pointers[swift_type_stars],
                             log->number_fields[swift_type_stars], &size,
                             &mask);

  /* Add the flag */
  if (flag != csds_flag_none) {
    mask |= log->list_fields[csds_index_special_flags].mask;
    size += size_special_flag;
  }
  size += csds_header_bytes;
  size_t size_total = count * size;

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size_total, &offset_new);
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

  /* Special flags */
  if (mask & log->list_fields[csds_index_special_flags].mask) {
    memcpy(buff, &special_flags,
           log->list_fields[csds_index_special_flags].size);
    buff += log->list_fields[csds_index_special_flags].size;
    mask &= ~log->list_fields[csds_index_special_flags].mask;
  }

  /* Write the gravity fields */
  for (int i = 0; i < log->number_fields[swift_type_dark_matter]; i++) {
    struct csds_field *field = &log->field_pointers[swift_type_dark_matter][i];

    /* Skip the fields that are not required. */
    if (!(mask & field->mask)) continue;

    /* Do we have a conversion function? */
    if (field->conversion_grav) {
      char *tmp_buff = field->conversion_grav(gp, e, buff);
      /* Check that the correct number of bits are written */
      if ((tmp_buff - buff) != (long int)field->size) {
        error("The field %s wrote an unexpected number of bits", field->name);
      }
    }
    /* Write it manually */
    else {
      memcpy(buff, ((char *)gp) + field->offset, field->size);
    }

    /* Update the variables */
    buff += field->size;
    mask &= ~field->mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in gparts. %u", mask);
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
  const int size_special_flag = log->list_fields[csds_index_special_flags].size;
  const uint32_t special_flags =
      csds_pack_flags_and_data(flag, flag_data, swift_type_dark_matter);

  /* Compute the size of the buffer. */
  /* As we might have some non DM particles, we cannot log_all_fields blindly */
  int count_dm = 0;
  // TODO: write only some fields
  for (int i = 0; i < count; i++) {
    /* Log only the dark matter */
    if (p[i].type != swift_type_dark_matter &&
        p[i].type != swift_type_dark_matter_background)
      continue;

    count_dm += 1;
  }
  unsigned int mask = 0;
  size_t size = 0;
  csds_compute_size_and_mask(log->field_pointers[swift_type_dark_matter],
                             log->number_fields[swift_type_dark_matter], &size,
                             &mask);

  /* Add the flag */
  if (flag != csds_flag_none) {
    mask |= log->list_fields[csds_index_special_flags].mask;
    size += size_special_flag;
  }
  size += csds_header_bytes;
  size_t size_total = size * count_dm;

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size_total, &offset_new);
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
  struct dump *dump = &log->dump;
  /* Start by computing the size of the message. */
  const int size =
      log->list_fields[csds_index_timestamp].size + csds_header_bytes;

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  unsigned int mask = log->list_fields[csds_index_timestamp].mask;
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

  /* ensure enough space in dump */
  dump_ensure(&log->dump, limit, log->buffer_scale * limit);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/** @brief Generate the name of the dump files
 *
 * @param log The #csds_writer.
 * @param filename The filename of the dump file.
 */
void csds_get_dump_name(struct csds_writer *log, char *filename) {
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
  if (csds_index_special_flags != 0) {
    error("Expecting the special flags to be the first element.");
  }
  csds_define_common_field(list[csds_index_special_flags], "SpecialFlags",
                           sizeof(uint32_t));
  num_fields += 1;

  /* Add the timestamp */
  if (csds_index_timestamp != 1) {
    error("Expecting the timestamp to be the first element.");
  }
  csds_define_common_field(list[csds_index_timestamp], "Timestamp",
                           sizeof(integertime_t) + sizeof(double));
  list[num_fields].type = mask_type_timestep;  // flag it as timestamp
  num_fields += 1;

  /* Initialize all the particles types */
  for (int i = 0; i < swift_type_count; i++) {
    int tmp_num_fields = 0;
    struct csds_field *current = &list[num_fields];
    enum mask_type mask_type;

    /* Set the pointer */
    log->field_pointers[i] = current;

    switch (i) {
      /* Hydro */
      case swift_type_gas:
        /* Set the mask type */
        mask_type = mask_type_gas;

        /* Set the masks */
        tmp_num_fields = csds_hydro_define_fields(current);
        tmp_num_fields +=
            csds_chemistry_define_fields_parts(current + tmp_num_fields);
        // TODO add cooling + SF
        break;

      /* Stars */
      case swift_type_stars:
        /* Set the mask type */
        mask_type = mask_type_stars;

        /* Set the masks */
        tmp_num_fields = csds_stars_define_fields(current);
        tmp_num_fields +=
            csds_chemistry_define_fields_sparts(current + tmp_num_fields);
        tmp_num_fields +=
            csds_star_formation_define_fields(current + tmp_num_fields);
        break;

      case swift_type_dark_matter:
        /* Set the mask type */
        mask_type = mask_type_dark_matter;

        /* Set the masks */
        tmp_num_fields = csds_gravity_define_fields(current);
        break;

      default:
        log->field_pointers[i] = NULL;
        break;
    }

    /* Set the particle type */
    for (int j = 0; j < tmp_num_fields; j++) {
      current[j].type = mask_type;
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
  if (mask >= 8 * csds_mask_size) {
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

  /* generate dump filename. */
  char csds_name_file[PARSER_MAX_LINE_SIZE];
  csds_get_dump_name(log, csds_name_file);

  /* Compute max size for a particle record. */
  int max_size = csds_offset_size + csds_mask_size;

  /* Loop over all fields except timestamp. */
  for (int i = 0; i < log->total_number_fields; i++) {
    /* Skip the timestamp */
    if (i == csds_index_timestamp) continue;

    max_size += log->list_fields[i].size;
  }
  log->max_record_size = max_size;

  /* init dump. */
  dump_init(&log->dump, csds_name_file, buffer_size);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Close dump file and desallocate memory
 *
 * @param log The #csds_writer
 */
void csds_free(struct csds_writer *log) {
  dump_close(&log->dump);

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
  struct dump *dump = &log->dump;

  uint64_t file_offset = dump->file_offset;

  if (file_offset != 0)
    error(
        "The CSDS is not empty."
        "This function should be called before writing anything in the CSDS");

  /* Write format information. */
  csds_write_data(dump, &file_offset, csds_format_size, &csds_file_format);

  /* Write the major version number. */
  int major = csds_major_version;
  csds_write_data(dump, &file_offset, sizeof(int), &major);

  /* Write the minor version number. */
  int minor = csds_minor_version;
  csds_write_data(dump, &file_offset, sizeof(int), &minor);

  /* write offset direction. */
  const int reversed = 0;
  csds_write_data(dump, &file_offset, sizeof(int), &reversed);

  /* placeholder to write the offset of the first log here. */
  char *skip_header = dump_get(dump, csds_offset_size, &file_offset);

  /* write number of bytes used for names. */
  const unsigned int label_size = csds_string_length;
  csds_write_data(dump, &file_offset, sizeof(unsigned int), &label_size);

  /* placeholder to write the number of unique masks. */
  char *skip_unique_masks = dump_get(dump, sizeof(unsigned int), &file_offset);

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
    csds_write_data(dump, &file_offset, csds_string_length,
                    &log->list_fields[i].name);

    // mask size.
    csds_write_data(dump, &file_offset, sizeof(unsigned int),
                    &log->list_fields[i].size);
  }
  memcpy(skip_unique_masks, &unique_mask, sizeof(unsigned int));

  /* Write the number of fields per particle */
  csds_write_data(dump, &file_offset, sizeof(log->number_fields),
                  log->number_fields);

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
          csds_write_data(dump, &file_offset, sizeof(int), &m);
          break;
        }
      }
    }
  }

  /* last step: write first offset. */
  memcpy(skip_header, &file_offset, csds_offset_size);
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
  memcpy(mask, buff, csds_mask_size);
  buff += csds_mask_size;

  *offset = 0;
  memcpy(offset, buff, csds_offset_size);
  *offset = cur_offset - *offset;

  return csds_mask_size + csds_offset_size;
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
        (log->list_fields[i].type == mask_type_gas)) {

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
        (log->list_fields[i].type == mask_type_dark_matter)) {

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
  if (!(mask & log->list_fields[csds_index_timestamp].mask))
    error("Trying to read timestamp from a particle.");

  /* Make sure we don't have extra fields. */
  if (mask != log->list_fields[csds_index_timestamp].mask)
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
    if (log->field_pointers == NULL) continue;

    log->field_pointers[i] =
        log->list_fields + (log->field_pointers[i] - old_list_fields);
  }

  /* Restart the dump file. */
  char csds_name_file[PARSER_MAX_LINE_SIZE];
  csds_get_dump_name(log, csds_name_file);

  dump_restart(&log->dump, csds_name_file);
}

#endif /* WITH_CSDS */

#endif /* HAVE_POSIX_FALLOCATE */
