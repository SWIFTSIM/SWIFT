/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#include "hydro/Shadowswift/voronoi2d_algorithm.h"

VORONOI_DECLARE_GLOBAL_VARIABLES()

int main() {

  float anchor[3] = {-0.5f, -0.5f, -0.5f};
  float side[3] = {2.0f, 2.0f, 2.0f};
  voronoi_set_box(anchor, side);

  struct voronoi_cell cell;
  double x[3] = {0.5, 0.5, 0.5};

  voronoi_cell_init(&cell, x);

  float maxradius = voronoi_cell_finalize(&cell);

  assert(maxradius == 2.0f * sqrtf(2.0f));

  assert(cell.volume == 4.0f);

  assert(cell.centroid[0] == 0.5f);
  assert(cell.centroid[1] == 0.5f);

  /* reinitialize cell */
  voronoi_cell_init(&cell, x);

  float dx[2] = {0.25, 0.25};
  voronoi_cell_interact(&cell, dx, 0);

  voronoi_print_cell(&cell);

  return 0;
}
