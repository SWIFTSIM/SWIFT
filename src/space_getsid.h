/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_SPACE_GETSID_H
#define SWIFT_SPACE_GETSID_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stddef.h>

/* Includes. */
#include "cell.h"
#include "runner.h"
#include "space.h"

/**
 * @brief Get the shift-id of the given pair of cells, swapping them
 *      if need be.
 *
 * WARNING: This function may swap the cells ci and cj.
 *
 * @param s The space
 * @param ci Pointer to first #cell.
 * @param cj Pointer second #cell.
 * @param shift (return) Vector from ci to cj.
 *
 * @return The shift ID and set shift, may or may not swap ci and cj.
 */
__attribute__((always_inline)) INLINE static int space_getsid(
    const struct space *s, struct cell **ci, struct cell **cj,
    double shift[3]) {

  /* Get the relative distance between the pairs, wrapping. */
  const int periodic = s->periodic;
  double dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = (*cj)->loc[k] - (*ci)->loc[k];
    if (periodic && dx[k] < -s->dim[k] / 2)
      shift[k] = s->dim[k];
    else if (periodic && dx[k] > s->dim[k] / 2)
      shift[k] = -s->dim[k];
    else
      shift[k] = 0.0;
    dx[k] += shift[k];
  }

  /* Get the sorting index. */
  int sid = 0;
  for (int k = 0; k < 3; k++)
    sid = 3 * sid + ((dx[k] < 0.0) ? 0 : ((dx[k] > 0.0) ? 2 : 1));

  /* Switch the cells around? */
  if (runner_flip[sid]) {
    struct cell *temp = *ci;
    *ci = *cj;
    *cj = temp;
    for (int k = 0; k < 3; k++) shift[k] = -shift[k];
  }
  sid = sortlistID[sid];

  /* Return the sort ID. */
  return sid;
}

#endif /* SWIFT_SPACE_GETSID_H */
