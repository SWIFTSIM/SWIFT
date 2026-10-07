/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Pure-function unit test for cell_need_rebuild_for_hydro_aperture_pair()
 * (src/cell.h), the rebuild criterion for the fixed-aperture gas-gas sink
 * formation preparation loop.
 *
 * The predicate is: rebuild if r_cut + ci->hydro.dx_max_part +
 * cj->hydro.dx_max_part > cj->dmin. This test constructs two bare cells (no
 * particles, no space/engine needed) and checks the boolean strictly below,
 * exactly at, and strictly above that boundary, for both self-consistent
 * argument orders.
 */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <string.h>

/* Local headers. */
#include "cell.h"
#include "error.h"

static struct cell make_bare_cell(double dmin, float dx_max_part) {
  struct cell c;
  memset(&c, 0, sizeof(struct cell));
  c.dmin = dmin;
  c.hydro.dx_max_part = dx_max_part;
  return c;
}

int main(int argc, char *argv[]) {

  const double dmin = 1.0;

  /* Split the margin between the two cells asymmetrically to also exercise
   * the non-symmetric dx_max_part_i + dx_max_part_j sum, not just 2*dx. */
  const float dx_i = 0.06f;
  const float dx_j = 0.02f;
  const float dx_sum = dx_i + dx_j;

  /* --- Strictly below the boundary: r_cut + dx_sum < dmin --- */
  {
    const float r_cut = (float)dmin - dx_sum - 0.05f;
    struct cell ci = make_bare_cell(dmin, dx_i);
    struct cell cj = make_bare_cell(dmin, dx_j);

    if (cell_need_rebuild_for_hydro_aperture_pair(&ci, &cj, r_cut))
      error("Below the boundary: rebuild predicate incorrectly returned 1.");
    if (cell_need_rebuild_for_hydro_aperture_pair(&cj, &ci, r_cut))
      error(
          "Below the boundary (swapped args): rebuild predicate incorrectly "
          "returned 1.");
  }

  /* --- Exactly at the boundary: r_cut + dx_sum == dmin (not > dmin, so no
   * rebuild -- the predicate is a strict inequality). --- */
  {
    const float r_cut = (float)dmin - dx_sum;
    struct cell ci = make_bare_cell(dmin, dx_i);
    struct cell cj = make_bare_cell(dmin, dx_j);

    if (cell_need_rebuild_for_hydro_aperture_pair(&ci, &cj, r_cut))
      error("At the boundary: rebuild predicate incorrectly returned 1.");
    if (cell_need_rebuild_for_hydro_aperture_pair(&cj, &ci, r_cut))
      error(
          "At the boundary (swapped args): rebuild predicate incorrectly "
          "returned 1.");
  }

  /* --- Strictly above the boundary: r_cut + dx_sum > dmin --- */
  {
    const float r_cut = (float)dmin - dx_sum + 0.05f;
    struct cell ci = make_bare_cell(dmin, dx_i);
    struct cell cj = make_bare_cell(dmin, dx_j);

    if (!cell_need_rebuild_for_hydro_aperture_pair(&ci, &cj, r_cut))
      error("Above the boundary: rebuild predicate incorrectly returned 0.");
    if (!cell_need_rebuild_for_hydro_aperture_pair(&cj, &ci, r_cut))
      error(
          "Above the boundary (swapped args): rebuild predicate incorrectly "
          "returned 0.");
  }

  message("All cell_need_rebuild_for_hydro_aperture_pair() checks passed.");

  return 0;
}
