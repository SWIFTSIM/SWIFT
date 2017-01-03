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

/* Some standard headers. */
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* This object's header. */
#include "../src/logger.h"

/* Local headers. */
#include "../src/dump.h"
#include "../src/part.h"

int main(int argc, char *argv[]) {

  /* Some constants. */
  const char *filename = "/tmp/dump_test.out";

  /* Prepare a dump. */
  struct dump d;
  dump_init(&d, filename, 1024 * 1024);

  /* Write several copies of a part to the dump. */
  struct part p;
  bzero(&p, sizeof(struct part));
  size_t offset = 0;

  /* Write the full part. */
  logger_log_part(&p, logger_mask_x | logger_mask_v | logger_mask_a |
                          logger_mask_u | logger_mask_h | logger_mask_rho |
                          logger_mask_consts,
                  &offset, &d);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Write only the position. */
  p.x[0] = 1.0;
  logger_log_part(&p, logger_mask_x, &offset, &d);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Write the position and velocity. */
  p.x[0] = 2.0;
  p.v[0] = 2.0;
  logger_log_part(&p, logger_mask_x | logger_mask_v, &offset, &d);
  printf("Wrote part at offset %#016zx.\n", offset);

  /* Recover the last part from the dump. */
  bzero(&p, sizeof(struct part));
  size_t offset_old = offset;
  int mask = logger_read_part(&p, &offset, d.data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 2.0 || p.v[0] != 2.0) {
    printf("FAIL: could not read position and velocity of stored particle.");
    abort();
  }

  /* Recover the second part from the dump. */
  bzero(&p, sizeof(struct part));
  offset_old = offset;
  mask = logger_read_part(&p, &offset, d.data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 1.0 || p.v[0] != 0.0) {
    printf("FAIL: could not read position and velocity of stored particle.");
    abort();
  }

  /* Recover the first part from the dump. */
  bzero(&p, sizeof(struct part));
  offset_old = offset;
  mask = logger_read_part(&p, &offset, d.data);
  printf(
      "Recovered part at offset %#016zx with mask %#04x: p.x[0]=%e, "
      "p.v[0]=%e.\n",
      offset_old, mask, p.x[0], p.v[0]);
  if (p.x[0] != 0.0 || p.v[0] != 0.0) {
    printf("FAIL: could not read position and velocity of stored particle.");
    abort();
  }

  /* Finalize the dump. */
  dump_close(&d);

  /* Return a happy number. */
  printf("PASS\n");
  return 0;
}
