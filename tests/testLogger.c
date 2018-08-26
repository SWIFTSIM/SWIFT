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

#ifdef HAVE_POSIX_FALLOCATE /* Are we on a sensible platform? */

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

void test_log_parts(struct dump *d) {

  /* Write several copies of a part to the dump. */
  struct part p;
  bzero(&p, sizeof(struct part));
  p.x[0] = 1.0;
  p.v[0] = 0.1;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Write the full part. */
  logger_log_part(&p,
                  logger_mask_x | logger_mask_v | logger_mask_a |
                      logger_mask_u | logger_mask_h | logger_mask_rho |
                      logger_mask_consts,
                  &offset, d);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Write only the position. */
  p.x[0] = 2.0;
  logger_log_part(&p, logger_mask_x, &offset, d);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Write the position and velocity. */
  p.x[0] = 3.0;
  p.v[0] = 0.3;
  logger_log_part(&p, logger_mask_x | logger_mask_v, &offset, d);
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

void test_log_gparts(struct dump *d) {

  /* Write several copies of a part to the dump. */
  struct gpart p;
  bzero(&p, sizeof(struct gpart));
  p.x[0] = 1.0;
  p.v_full[0] = 0.1;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Write the full part. */
  logger_log_gpart(&p,
                   logger_mask_x | logger_mask_v | logger_mask_a |
                       logger_mask_h | logger_mask_consts,
                   &offset, d);
  printf("Wrote gpart at offset %#016zx.\n", offset);

  /* Write only the position. */
  p.x[0] = 2.0;
  logger_log_gpart(&p, logger_mask_x, &offset, d);
  printf("Wrote gpart at offset %#016zx.\n", offset);

  /* Write the position and velocity. */
  p.x[0] = 3.0;
  p.v_full[0] = 0.3;
  logger_log_gpart(&p, logger_mask_x | logger_mask_v, &offset, d);
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

void test_log_timestamps(struct dump *d) {

  /* The timestamp to log. */
  unsigned long long int t = 10;

  /* Start with an offset at the end of the dump. */
  size_t offset = d->count;

  /* Log three consecutive timestamps. */
  logger_log_timestamp(t, &offset, d);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);
  t += 10;
  logger_log_timestamp(t, &offset, d);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);
  t += 10;
  logger_log_timestamp(t, &offset, d);
  printf("Logged timestamp %020llu at offset %#016zx.\n", t, offset);

  /* Recover the three timestamps. */
  size_t offset_old = offset;
  t = 0;
  int mask = logger_read_timestamp(&t, &offset, (const char *)d->data);
  printf("Recovered timestamp %020llu at offset %#016zx with mask %#04x.\n", t,
         offset_old, mask);
  if (t != 30) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }

  offset_old = offset;
  t = 0;
  mask = logger_read_timestamp(&t, &offset, (const char *)d->data);
  printf("Recovered timestamp %020llu at offset %#016zx with mask %#04x.\n", t,
         offset_old, mask);
  if (t != 20) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }

  offset_old = offset;
  t = 0;
  mask = logger_read_timestamp(&t, &offset, (const char *)d->data);
  printf("Recovered timestamp %020llu at offset %#016zx with mask %#04x.\n", t,
         offset_old, mask);
  if (t != 10) {
    printf("FAIL: could not recover correct timestamp.\n");
    abort();
  }
}

int main(int argc, char *argv[]) {

  /* Some constants. */
  char filename[256];
  const int now = time(NULL);
  sprintf(filename, "/tmp/SWIFT_logger_test_%d.out", now);

  /* Prepare a dump. */
  struct dump d;
  dump_init(&d, filename, 1024 * 1024);

  /* Test writing/reading parts. */
  test_log_parts(&d);

  /* Test writing/reading gparts. */
  test_log_gparts(&d);

  /* Test writing/reading timestamps. */
  test_log_timestamps(&d);

  /* Finalize the dump. */
  dump_close(&d);

  /* Be clean */
  remove(filename);

  /* Return a happy number. */
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* HAVE_POSIX_FALLOCATE */
