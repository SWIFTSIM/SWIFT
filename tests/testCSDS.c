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
    defined(WITH_CSDS) /* Are we on a sensible platform? */

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

void test_log_parts(struct csds_writer *log) {
  struct dump *d = &log->dump;
  struct engine e;

  /* Write several copies of a part to the dump. */
  struct part p;
  struct xpart xp;
  bzero(&p, sizeof(struct part));
  bzero(&xp, sizeof(struct xpart));
  p.x[0] = 1.0;
  p.v[0] = 0.1;
  xp.csds_data.last_offset = 0;

  /* Write the full part. */
  csds_log_part(log, &p, &xp, &e, /* log_all */ 1, csds_flag_none,
                /* flag_data */ 0);
  printf("Wrote part at offset %#016zx.\n", xp.csds_data.last_offset);

  /* Write only the position. */
  p.x[0] = 2.0;
  p.v[0] = 0.;
  csds_log_part(log, &p, &xp, &e, /* log_all */ 0, csds_flag_none,
                /* flag_data */ 0);
  printf("Wrote part at offset %#016zx.\n", xp.csds_data.last_offset);

  /* Write the position and velocity. */
  p.x[0] = 3.0;
  p.v[0] = 0.3;
  csds_log_part(log, &p, &xp, &e, /* log_all */ 0, csds_flag_none,
                /* flag_data */ 0);
  printf("Wrote part at offset %#016zx.\n", xp.csds_data.last_offset);

  /* Recover the last part from the dump. */
  bzero(&p, sizeof(struct part));
  size_t offset = xp.csds_data.last_offset;
  size_t offset_old = offset;
  unsigned int mask = csds_read_part(log, &p, &offset, (const char *)d->data);
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
  mask = csds_read_part(log, &p, &offset, (const char *)d->data);
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
  mask = csds_read_part(log, &p, &offset, (const char *)d->data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 1.0 || p.v[0] != 0.1f) {
    printf("FAIL: could not read position and velocity of stored particle.\n");
    abort();
  }
}

void test_log_gparts(struct csds_writer *log) {
  struct dump *d = &log->dump;
  struct engine e;

  /* Write several copies of a part to the dump. */
  struct gpart p;
  bzero(&p, sizeof(struct gpart));
  p.x[0] = 1.0;
  p.v_full[0] = 0.1;
  p.type = swift_type_dark_matter;
  p.csds_data.last_offset = 0;

  /* Write the full part. */
  csds_log_gpart(log, &p, &e, /* log_all */ 1, csds_flag_none,
                 /* flag_data */ 0);
  printf("Wrote gpart at offset %#016zx.\n", p.csds_data.last_offset);

  /* Write only the position. */
  p.x[0] = 2.0;
  p.v_full[0] = 0.;
  csds_log_gpart(log, &p, &e, /* log_all */ 0, csds_flag_none,
                 /* flag_data */ 0);
  printf("Wrote gpart at offset %#016zx.\n", p.csds_data.last_offset);

  /* Write the position and velocity. */
  p.x[0] = 3.0;
  p.v_full[0] = 0.3;
  csds_log_gpart(log, &p, &e, /* log_all */ 0, csds_flag_none,
                 /* flag_data */ 0);
  printf("Wrote gpart at offset %#016zx.\n", p.csds_data.last_offset);

  /* Recover the last part from the dump. */
  size_t offset = p.csds_data.last_offset;
  bzero(&p, sizeof(struct gpart));
  size_t offset_old = offset;
  int mask = csds_read_gpart(log, &p, &offset, (const char *)d->data);
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
  mask = csds_read_gpart(log, &p, &offset, (const char *)d->data);
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
  mask = csds_read_gpart(log, &p, &offset, (const char *)d->data);
  printf(
      "Recovered gpart at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v_full[0]);
  if (p.x[0] != 1.0 || p.v_full[0] != 0.1f) {
    printf("FAIL: could not read position and velocity of stored gpart.\n");
    abort();
  }
}

void test_log_timestamps(struct csds_writer *log) {
  struct dump *d = &log->dump;

  /* The timestamp to log. */
  integertime_t t = 10;
  double time = 0.1;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Log three consecutive timestamps. */
  csds_log_timestamp(log, t, time, &offset);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);
  t += 10;
  time = 0.2;
  csds_log_timestamp(log, t, time, &offset);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);
  t += 10;
  time = 0.3;
  csds_log_timestamp(log, t, time, &offset);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);

  /* Recover the three timestamps. */
  size_t offset_old = offset;
  t = 0;
  time = 0;
  int mask =
      csds_read_timestamp(log, &t, &time, &offset, (const char *)d->data);
  printf(
      "Recovered timestamp %020llu with time %g at offset %#016zx with mask "
      "%#04x.\n",
      t, time, offset_old, mask);
  if (t != 30) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }
  if (time != 0.3) {
    printf("FAIL: could not recover correct time.\n");
    abort();
  }

  offset_old = offset;
  t = 0;
  time = 0;
  mask = csds_read_timestamp(log, &t, &time, &offset, (const char *)d->data);
  printf(
      "Recovered timestamp %020llu with time %g at offset %#016zx with mask "
      "%#04x.\n",
      t, time, offset_old, mask);
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
  mask = csds_read_timestamp(log, &t, &time, &offset, (const char *)d->data);
  printf(
      "Recovered timestamp %020llu with time %g at offset %#016zx with mask "
      "%#04x.\n",
      t, time, offset_old, mask);
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

  /* Prepare a csds. */
  struct csds_writer log;
  struct swift_params params;
  struct engine e;
  e.policy = engine_policy_hydro | engine_policy_self_gravity;

  parser_read_file("csds.yml", &params);
  csds_init(&log, &e, &params);

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

  /* Clean the csds. */
  csds_free(&log);

  /* Return a happy number. */
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* HAVE_POSIX_FALLOCATE */
