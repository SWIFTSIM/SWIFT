/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#include "hydro/Shadowswift/voronoi1d_algorithm.h"

int main() {

  /* Create a Voronoi cell */
  double x[1] = {0.5f};
  struct voronoi_cell cell;
  voronoi_cell_init(&cell, x);

  /* Interact with a left and right neighbour */
  float xL[1] = {0.5f};
  float xR[1] = {-0.5f};
  voronoi_cell_interact(&cell, xL, 1);
  voronoi_cell_interact(&cell, xR, 2);

  /* Interact with some more neighbours to check if they are properly ignored */
  float x0[1] = {0.6f};
  float x1[1] = {-0.7f};
  voronoi_cell_interact(&cell, x0, 3);
  voronoi_cell_interact(&cell, x1, 4);

  /* Finalize cell and check results */
  voronoi_cell_finalize(&cell);

  if (cell.volume != 0.5f) {
    error("Wrong volume: %g!", cell.volume);
  }
  if (cell.centroid != 0.5f) {
    error("Wrong centroid: %g!", cell.centroid);
  }
  if (cell.idL != 1) {
    error("Wrong left neighbour: %llu!", cell.idL);
  }
  if (cell.idR != 2) {
    error("Wrong right neighbour: %llu!", cell.idR);
  }

  return 0;
}
