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

#include <stdlib.h>
#include "error.h"
#include "hydro/Shadowswift/voronoi3d_algorithm.h"

/**
 * @brief Check if voronoi_volume_tetrahedron() works
 */
void test_voronoi_volume_tetrahedron() {
  float v1[3] = {0., 0., 0.};
  float v2[3] = {0., 0., 1.};
  float v3[3] = {0., 1., 0.};
  float v4[3] = {1., 0., 0.};

  float V = voronoi_volume_tetrahedron(v1, v2, v3, v4);
  assert(V == 1.0f / 6.0f);
}

/**
 * @brief Check if voronoi_centroid_tetrahedron() works
 */
void test_voronoi_centroid_tetrahedron() {
  float v1[3] = {0., 0., 0.};
  float v2[3] = {0., 0., 1.};
  float v3[3] = {0., 1., 0.};
  float v4[3] = {1., 0., 0.};

  float centroid[3];
  voronoi_centroid_tetrahedron(centroid, v1, v2, v3, v4);
  assert(centroid[0] == 0.25f);
  assert(centroid[1] == 0.25f);
  assert(centroid[2] == 0.25f);
}

/**
 * @brief Check if voronoi_calculate_cell() works
 */
void test_calculate_cell() {
  struct voronoi_cell cell;

  cell.x[0] = 0.5f;
  cell.x[1] = 0.5f;
  cell.x[2] = 0.5f;

  /* Initialize the cell to a large cube. */
  voronoi_initialize(&cell);
  /* Calculate the volume and centroid of the large cube. */
  voronoi_calculate_cell(&cell);
  /* Calculate the faces. */
  voronoi_calculate_faces(&cell);

  /* Update these values if you ever change to another large cube! */
  assert(cell.volume == 27.0f);
  assert(cell.centroid[0] = 0.5f);
  assert(cell.centroid[1] = 0.5f);
  assert(cell.centroid[2] = 0.5f);

  /* Check cell neighbours. */
  assert(cell.nface == 6);
  assert(cell.ngbs[0] == VORONOI3D_BOX_FRONT);
  assert(cell.ngbs[1] == VORONOI3D_BOX_LEFT);
  assert(cell.ngbs[2] == VORONOI3D_BOX_BOTTOM);
  assert(cell.ngbs[3] == VORONOI3D_BOX_TOP);
  assert(cell.ngbs[4] == VORONOI3D_BOX_BACK);
  assert(cell.ngbs[5] == VORONOI3D_BOX_RIGHT);

  /* Check cell faces */
  assert(cell.face_areas[0] == 9.0f);
  assert(cell.face_midpoints[0][0] == 0.5f);
  assert(cell.face_midpoints[0][1] == -1.0f);
  assert(cell.face_midpoints[0][2] == 0.5f);

  assert(cell.face_areas[1] == 9.0f);
  assert(cell.face_midpoints[1][0] == -1.0f);
  assert(cell.face_midpoints[1][1] == 0.5f);
  assert(cell.face_midpoints[1][2] == 0.5f);

  assert(cell.face_areas[2] == 9.0f);
  assert(cell.face_midpoints[2][0] == 0.5f);
  assert(cell.face_midpoints[2][1] == 0.5f);
  assert(cell.face_midpoints[2][2] == -1.0f);

  assert(cell.face_areas[3] == 9.0f);
  assert(cell.face_midpoints[3][0] == 0.5f);
  assert(cell.face_midpoints[3][1] == 0.5f);
  assert(cell.face_midpoints[3][2] == 2.0f);

  assert(cell.face_areas[4] == 9.0f);
  assert(cell.face_midpoints[4][0] == 0.5f);
  assert(cell.face_midpoints[4][1] == 2.0f);
  assert(cell.face_midpoints[4][2] == 0.5f);

  assert(cell.face_areas[5] == 9.0f);
  assert(cell.face_midpoints[5][0] == 2.0f);
  assert(cell.face_midpoints[5][1] == 0.5f);
  assert(cell.face_midpoints[5][2] == 0.5f);
}

int main() {

  /* Check basic Voronoi cell functions */
  test_voronoi_volume_tetrahedron();
  test_voronoi_centroid_tetrahedron();
  test_calculate_cell();

  /* Create a Voronoi cell */
  double x[3] = {0.5f, 0.5f, 0.5f};
  struct voronoi_cell cell;
  voronoi_cell_init(&cell, x);

  /* Interact with neighbours */
  float x0[3] = {0.5f, 0.0f, 0.0f};
  float x1[3] = {-0.5f, 0.0f, 0.0f};
  float x2[3] = {0.0f, 0.5f, 0.0f};
  float x3[3] = {0.0f, -0.5f, 0.0f};
  float x4[3] = {0.0f, 0.0f, 0.5f};
  float x5[3] = {0.0f, 0.0f, -0.5f};
  voronoi_cell_interact(&cell, x0, 1);
  voronoi_cell_interact(&cell, x1, 2);
  voronoi_cell_interact(&cell, x2, 3);
  voronoi_cell_interact(&cell, x3, 4);
  voronoi_cell_interact(&cell, x4, 5);
  voronoi_cell_interact(&cell, x5, 6);

  /* Interact with some more neighbours to check if they are properly ignored */
  float xE0[3] = {0.6f, 0.0f, 0.1f};
  float xE1[3] = {-0.7f, 0.2f, 0.04f};
  voronoi_cell_interact(&cell, xE0, 7);
  voronoi_cell_interact(&cell, xE1, 8);

  /* Finalize cell and check results */
  voronoi_cell_finalize(&cell);

  if (cell.volume != 0.125f) {
    error("Wrong volume: %g!", cell.volume);
  }
  if (fabs(cell.centroid[0] - 0.5f) > 1.e-5f ||
      fabs(cell.centroid[1] - 0.5f) > 1.e-5f ||
      fabs(cell.centroid[2] - 0.5f) > 1.e-5f) {
    error("Wrong centroid: %g %g %g!", cell.centroid[0], cell.centroid[1],
          cell.centroid[2]);
  }

  /* Check faces. */
  float A, midpoint[3];
  A = voronoi_get_face(&cell, 1, midpoint);
  if (A) {
    if (A != 0.25f) {
      error("Wrong surface area: %g!", A);
    }
    if (fabs(midpoint[0] - 0.25f) > 1.e-5 || fabs(midpoint[1] - 0.5f) > 1.e-5 ||
        fabs(midpoint[2] - 0.5f) > 1.e-5) {
      error("Wrong face midpoint: %g %g %g!", midpoint[0], midpoint[1],
            midpoint[2]);
    }
  } else {
    error("Neighbour not found!");
  }

  return 0;
}
