/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

#if defined(HAVE_POSIX_FALLOCATE) && \
    defined(WITH_LOGGER) /* Are we on a sensible platform? */

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

void test_log_parts(struct logger_writer *log) {
  struct dump *d = &log->dump;

  /* Write several copies of a part to the dump. */
  struct part p;
  bzero(&p, sizeof(struct part));
  p.x[0] = 1.0;
  p.v[0] = 0.1;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Write the full part. */
  logger_log_part(log, &p,
                  logger_mask_data[logger_x].mask |
                              logger_mask_data[logger_v].mask |
                              logger_mask_data[logger_a].mask |
                              logger_mask_data[logger_u].mask |
                              logger_mask_data[logger_h].mask |
                              logger_mask_data[logger_rho].mask |
                              logger_mask_data[logger_consts].mask,
                  &offset);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Write only the position. */
  p.x[0] = 2.0;
  logger_log_part(log, &p, logger_mask_data[logger_x].mask, &offset);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Write the position and velocity. */
  p.x[0] = 3.0;
  p.v[0] = 0.3;
  logger_log_part(log, &p, logger_mask_data[logger_x].mask |
		  logger_mask_data[logger_v].mask
		  , &offset);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Recover the last part from the dump. */
  bzero(&p, sizeof(struct part));
  size_t offset_old = offset;
  int mask = logger_read_part(&p, &offset, (const char *)d->data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 3.0 || p.v[0] != 0.3f) {
    printf("FAIL: could not read position and velocity of stored particle.\n");
    abort();
  }

  /* Recover the second part from the dump (only position). */
  bzero(&p, sizeof(struct part));
  offset_old = offset;
  mask = logger_read_part(&p, &offset, (const char *)d->data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 2.0 || p.v[0] != 0.0) {
    printf("FAIL: could not read position and velocity of stored particle.\n");
    abort();
  }

  /* Recover the first part from the dump. */
  bzero(&p, sizeof(struct part));
  offset_old = offset;
  mask = logger_read_part(&p, &offset, (const char *)d->data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 1.0 || p.v[0] != 0.1f) {
    printf("FAIL: could not read position and velocity of stored particle.\n");
    abort();
  }
}

void test_log_gparts(struct logger_writer *log) {
  struct dump *d = &log->dump;

  /* Write several copies of a part to the dump. */
  struct gpart p;
  bzero(&p, sizeof(struct gpart));
  p.x[0] = 1.0;
  p.v_full[0] = 0.1;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Write the full part. */
  logger_log_gpart(log, &p,
                   logger_mask_data[logger_x].mask |
		   logger_mask_data[logger_v].mask |
		   logger_mask_data[logger_a].mask |
		   logger_mask_data[logger_h].mask |
		   logger_mask_data[logger_consts].mask,
                   &offset);
  printf("Wrote gpart at offset %#016zx.\n", offset);

  /* Write only the position. */
  p.x[0] = 2.0;
  logger_log_gpart(log, &p, logger_mask_data[logger_x].mask, &offset);
  printf("Wrote gpart at offset %#016zx.\n", offset);

  /* Write the position and velocity. */
  p.x[0] = 3.0;
  p.v_full[0] = 0.3;
  logger_log_gpart(log, &p, logger_mask_data[logger_x].mask |
		   logger_mask_data[logger_v].mask, &offset);
  printf("Wrote gpart at offset %#016zx.\n", offset);

  /* Recover the last part from the dump. */
  bzero(&p, sizeof(struct gpart));
  size_t offset_old = offset;
  int mask = logger_read_gpart(&p, &offset, (const char *)d->data);
  printf(
      "Recovered gpart at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v_full[0]);
  if (p.x[0] != 3.0 || p.v_full[0] != 0.3f) {
    printf("FAIL: could not read position and velocity of stored gpart.\n");
    abort();
  }

  /* Recover the second part from the dump. */
  bzero(&p, sizeof(struct gpart));
  offset_old = offset;
  mask = logger_read_gpart(&p, &offset, (const char *)d->data);
  printf(
      "Recovered gpart at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v_full[0]);
  if (p.x[0] != 2.0 || p.v_full[0] != 0.0) {
    printf("FAIL: could not read position and velocity of stored gpart.\n");
    abort();
  }

  /* Recover the first part from the dump. */
  bzero(&p, sizeof(struct gpart));
  offset_old = offset;
  mask = logger_read_gpart(&p, &offset, (const char *)d->data);
  printf(
      "Recovered gpart at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v_full[0]);
  if (p.x[0] != 1.0 || p.v_full[0] != 0.1f) {
    printf("FAIL: could not read position and velocity of stored gpart.\n");
    abort();
  }
}

void test_log_timestamps(struct logger_writer *log) {
  struct dump *d = &log->dump;

  /* The timestamp to log. */
  unsigned long long int t = 10;
  double time = 0.1;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Log three consecutive timestamps. */
  logger_log_timestamp(log, t, time, &offset);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);
  t += 10;
  time = 0.2;
  logger_log_timestamp(log, t, time, &offset);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);
  t += 10;
  time = 0.3;
  logger_log_timestamp(log, t, time, &offset);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);

  /* Recover the three timestamps. */
  size_t offset_old = offset;
  t = 0;
  time = 0;
  int mask = logger_read_timestamp(&t, &time, &offset, (const char *)d->data);
  printf("Recovered timestamp %020llu at offset %#016zx with mask %#04x.\n", t,
         offset_old, mask);
  if (t != 30) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }
  if (time != 0.3) {
    printf("FAIL: could not recover correct time %g.\n", time);
    abort();
  }

  offset_old = offset;
  t = 0;
  time = 0;
  mask = logger_read_timestamp(&t, &time, &offset, (const char *)d->data);
  printf("Recovered timestamp %020llu at offset %#016zx with mask %#04x.\n", t,
         offset_old, mask);
  if (t != 20) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }
  if (time != 0.2) {
    printf("FAIL: could not recover correct time.\n");
    abort();
  }

  offset_old = offset;
  t = 0;
  time = 0;
  mask = logger_read_timestamp(&t, &time, &offset, (const char *)d->data);
  printf("Recovered timestamp %020llu at offset %#016zx with mask %#04x.\n", t,
         offset_old, mask);
  if (t != 10) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }
  if (time != 0.1) {
    printf("FAIL: could not recover correct time.\n");
    abort();
  }
}

int main(int argc, char *argv[]) {

  /* Prepare a logger. */
  struct logger_writer log;
  struct swift_params params;
  parser_read_file("logger.yml", &params);
  logger_init(&log, &params);

  /* Test writing/reading parts. */
  test_log_parts(&log);

  /* Test writing/reading gparts. */
  test_log_gparts(&log);

  /* Test writing/reading timestamps. */
  test_log_timestamps(&log);

  /* Be clean */
  char filename[256];
  sprintf(filename, "%s.dump", log.base_name);
  remove(filename);

  /* Clean the logger. */
  logger_clean(&log);

  /* Return a happy number. */
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* HAVE_POSIX_FALLOCATE */
