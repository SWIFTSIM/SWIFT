/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
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
#include <config.h>

/* Some standard headers. */
#include <stdio.h>

/* Local headers. */
#include "runner_radiation_feedback.h"
#include "swift.h"

/**
 * @brief Compare the HII search shift of one cell pair with the expected one.
 *
 * @param label Name of the case, printed with the result.
 * @param periodic Is the box periodic?
 * @param loc_c Lower corner of the star's cell.
 * @param loc_cj Lower corner of the neighbour cell.
 * @param expected The expected shift.
 * @return 1 on a mismatch, 0 otherwise.
 */
static int check_shift(const char *label, const int periodic,
                       const double loc_c[3], const double loc_cj[3],
                       const double expected[3]) {

  /* A box of 3 top-level cells of unit width along each axis. */
  const double dim[3] = {3.0, 3.0, 3.0};
  double shift[3];
  runner_hii_get_pair_shift(dim, periodic, loc_c, loc_cj, shift);

  int fail = 0;
  for (int k = 0; k < 3; k++)
    if (shift[k] != expected[k]) fail = 1;

  printf("%s %s: shift = (%g, %g, %g), expected (%g, %g, %g)\n",
         fail ? "FAIL" : "ok  ", label, shift[0], shift[1], shift[2],
         expected[0], expected[1], expected[2]);
  return fail;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  int fails = 0;
  const double zero[3] = {0.0, 0.0, 0.0};

  /* Non-periodic box. A leaf in the upper half of the middle top-level cell
   * is 1.75 widths from the lower corner of its left neighbour, more than
   * half the box. The neighbour must be searched where it is. */
  {
    const double loc_c[3] = {1.75, 1.0, 1.0};
    const double loc_cj[3] = {0.0, 1.0, 1.0};
    fails +=
        check_shift("non-periodic, separation -1.75", 0, loc_c, loc_cj, zero);
  }

  /* Same in the positive direction along every axis at once. */
  {
    const double loc_c[3] = {0.0, 0.0, 0.0};
    const double loc_cj[3] = {2.0, 2.0, 2.0};
    fails += check_shift("non-periodic, separation +2", 0, loc_c, loc_cj, zero);
  }

  /* Periodic box. A leaf at the lower edge of top-level cell 0 meets top-level
   * cell 2 across the boundary: the neighbour moves down by one box. */
  {
    const double loc_c[3] = {0.0, 1.0, 1.0};
    const double loc_cj[3] = {2.0, 1.0, 1.0};
    const double expected[3] = {-3.0, 0.0, 0.0};
    fails += check_shift("periodic, wrap down", 1, loc_c, loc_cj, expected);
  }

  /* Periodic box, the mirror case: the neighbour moves up by one box. */
  {
    const double loc_c[3] = {2.75, 2.75, 1.0};
    const double loc_cj[3] = {0.0, 0.0, 1.0};
    const double expected[3] = {3.0, 3.0, 0.0};
    fails += check_shift("periodic, wrap up", 1, loc_c, loc_cj, expected);
  }

  /* Periodic box, adjacent cells inside the box: no shift. */
  {
    const double loc_c[3] = {1.0, 1.0, 1.0};
    const double loc_cj[3] = {0.0, 2.0, 1.0};
    fails += check_shift("periodic, direct neighbour", 1, loc_c, loc_cj, zero);
  }

  if (fails) {
    printf("%d case(s) failed.\n", fails);
    return 1;
  }
  return 0;
}
