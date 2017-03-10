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
#include "part.h"

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
  assert(cell.volume == voronoi_get_box_volume());
  float box_centroid[3];
  voronoi_get_box_centroid(box_centroid);
  assert(cell.centroid[0] = box_centroid[0]);
  assert(cell.centroid[1] = box_centroid[1]);
  assert(cell.centroid[2] = box_centroid[2]);

  /* Check cell neighbours. */
  assert(cell.nface == 6);
  assert(cell.ngbs[0] == VORONOI3D_BOX_FRONT);
  assert(cell.ngbs[1] == VORONOI3D_BOX_LEFT);
  assert(cell.ngbs[2] == VORONOI3D_BOX_BOTTOM);
  assert(cell.ngbs[3] == VORONOI3D_BOX_TOP);
  assert(cell.ngbs[4] == VORONOI3D_BOX_BACK);
  assert(cell.ngbs[5] == VORONOI3D_BOX_RIGHT);

  /* Check cell faces */
  float face_midpoint[3], face_area;
  face_area = voronoi_get_box_face(VORONOI3D_BOX_FRONT, face_midpoint);
  assert(cell.face_areas[0] == face_area);
  assert(cell.face_midpoints[0][0] == face_midpoint[0] - cell.x[0]);
  assert(cell.face_midpoints[0][1] == face_midpoint[1] - cell.x[1]);
  assert(cell.face_midpoints[0][2] == face_midpoint[2] - cell.x[2]);

  face_area = voronoi_get_box_face(VORONOI3D_BOX_LEFT, face_midpoint);
  assert(cell.face_areas[1] == face_area);
  assert(cell.face_midpoints[1][0] == face_midpoint[0] - cell.x[0]);
  assert(cell.face_midpoints[1][1] == face_midpoint[1] - cell.x[1]);
  assert(cell.face_midpoints[1][2] == face_midpoint[2] - cell.x[2]);

  face_area = voronoi_get_box_face(VORONOI3D_BOX_BOTTOM, face_midpoint);
  assert(cell.face_areas[2] == face_area);
  assert(cell.face_midpoints[2][0] == face_midpoint[0] - cell.x[0]);
  assert(cell.face_midpoints[2][1] == face_midpoint[1] - cell.x[1]);
  assert(cell.face_midpoints[2][2] == face_midpoint[2] - cell.x[2]);

  face_area = voronoi_get_box_face(VORONOI3D_BOX_TOP, face_midpoint);
  assert(cell.face_areas[3] == face_area);
  assert(cell.face_midpoints[3][0] == face_midpoint[0] - cell.x[0]);
  assert(cell.face_midpoints[3][1] == face_midpoint[1] - cell.x[1]);
  assert(cell.face_midpoints[3][2] == face_midpoint[2] - cell.x[2]);

  face_area = voronoi_get_box_face(VORONOI3D_BOX_BACK, face_midpoint);
  assert(cell.face_areas[4] == face_area);
  assert(cell.face_midpoints[4][0] == face_midpoint[0] - cell.x[0]);
  assert(cell.face_midpoints[4][1] == face_midpoint[1] - cell.x[1]);
  assert(cell.face_midpoints[4][2] == face_midpoint[2] - cell.x[2]);

  face_area = voronoi_get_box_face(VORONOI3D_BOX_RIGHT, face_midpoint);
  assert(cell.face_areas[5] == face_area);
  assert(cell.face_midpoints[5][0] == face_midpoint[0] - cell.x[0]);
  assert(cell.face_midpoints[5][1] == face_midpoint[1] - cell.x[1]);
  assert(cell.face_midpoints[5][2] == face_midpoint[2] - cell.x[2]);
}

void test_paths() {
  float u, l, q;
  int up, us, uw, lp, ls, lw, qp, qs, qw;
  float r2, dx[3];
  struct voronoi_cell cell;

  /* PATH 1.0 */
  // the first vertex is above the cutting plane and its first edge is below the
  // plane
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -1.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.nvert = 2;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.edges[0] = 1;
    cell.edgeindices[0] = 0;
    cell.edges[3] = 0;
    cell.edgeindices[3] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 0);
    assert(us == 0);
    assert(uw == 1);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(ls == 0);
    assert(lw == -1);
    assert(l == -0.75f);
  }

  /* PATH 1.1 */
  // the first vertex is above the cutting plane and its second edge is below
  // the plane
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 2.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.edges[0] = 1;
    cell.edges[1] = 2;
    cell.edges[6] = 0;
    cell.edgeindices[1] = 0;
    cell.edgeindices[6] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 0);
    assert(us == 1);
    assert(uw == 1);
    assert(u == 0.25f);
    assert(lp == 2);
    assert(ls == 0);
    assert(lw == -1);
    assert(l == -0.75f);
  }

  /* PATH 1.2 */
  // the first vertex is above the cutting plane and its second and last edge
  // is below the plane
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 2.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 2;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 2;
    cell.offsets[2] = 5;
    cell.edges[0] = 1;
    cell.edges[1] = 2;
    cell.edges[6] = 0;
    cell.edgeindices[1] = 0;
    cell.edgeindices[5] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 0);
    assert(us == 1);
    assert(uw == 1);
    assert(u == 0.25f);
    assert(lp == 2);
    assert(ls == 0);
    assert(lw == -1);
    assert(l == -0.75f);
  }

  /* PATH 1.3 */
  // the first vertex is above the cutting plane and is the closest vertex to
  // the plane. The code should crash.
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 2.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = 2.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.offsets[0] = 0;
    cell.edges[0] = 1;
    cell.edges[1] = 2;
    cell.edges[2] = 3;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == -1);
  }

  /* PATH 1.4.0 */
  // first vertex is above the plane, second vertex is closer and third vertex
  // lies below
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.edges[0] = 1;
    cell.edges[3] = 2;
    cell.edges[6] = 1;
    cell.edgeindices[0] = 2;
    cell.edgeindices[3] = 0;
    cell.edgeindices[6] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 0);
    assert(u == 0.125f);
    assert(lp == 2);
    assert(ls == 0);
    assert(l == -0.75f);
  }

  /* PATH 1.4.1 */
  // first vertex is above the plane, second vertex is closer and fourth vertex
  // is below
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = -1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.offsets[3] = 9;
    cell.edges[0] = 1;
    cell.edges[3] = 2;
    cell.edges[4] = 3;
    cell.edges[5] = 0;
    cell.edges[6] = 1;
    cell.edges[9] = 1;
    // this is the only difference between PATH 1.4.1 and PATH 1.4.2
    cell.edgeindices[0] = 3;
    cell.edgeindices[3] = 0;
    cell.edgeindices[4] = 0;
    cell.edgeindices[9] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 1);
    assert(u == 0.125f);
    assert(lp == 3);
    assert(ls == 0);
    assert(l == -0.75f);
  }

  /* PATH 1.4.2 */
  // first vertex is above the plane, second is closer, fourth is below
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = -1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.offsets[3] = 9;
    cell.edges[0] = 1;
    cell.edges[3] = 2;
    cell.edges[4] = 3;
    cell.edges[5] = 0;
    cell.edges[6] = 1;
    cell.edges[9] = 1;
    // this is the only difference between PATH 1.4.1 and PATH 1.4.2
    cell.edgeindices[0] = 2;
    cell.edgeindices[3] = 1;
    cell.edgeindices[4] = 0;
    cell.edgeindices[9] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 1);
    assert(u == 0.125f);
    assert(lp == 3);
    assert(ls == 0);
    assert(l == -0.75f);
  }

  /* PATH 1.4.3 */
  // first vertex is above the plane, second is closer, third is below
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.edges[0] = 1;
    // this is the difference between PATH 1.4.0 and this path
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edges[6] = 1;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[6] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 1);
    assert(u == 0.125f);
    assert(lp == 2);
    assert(ls == 0);
    assert(l == -0.75f);
  }

  /* PATH 1.4.4 */
  // first vertex is above the plane, second is closer, fourth is below
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = -1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 4;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 7;
    cell.offsets[3] = 10;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edges[5] = 3;
    cell.edges[7] = 1;
    cell.edges[10] = 1;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[5] = 0;
    cell.edgeindices[10] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 2);
    assert(u == 0.125f);
    assert(lp == 3);
    assert(ls == 0);
    assert(l == -0.75f);
  }

  /* PATH 1.4.5 */
  // same as 1.4.4, but with an order 3 second vertex
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = -1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.offsets[3] = 9;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edges[5] = 3;
    cell.edges[6] = 1;
    cell.edges[9] = 1;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[5] = 0;
    cell.edgeindices[9] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 2);
    assert(u == 0.125f);
    assert(lp == 3);
    assert(ls == 0);
    assert(l == -0.75f);
  }

  /* PATH 1.4.6 */
  // first vertex is above the plane, second is closer and is the closest
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.75f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 2;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 5;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == -1);
  }

  /* PATH 1.5 */
  // first vertex is above the plane, second vertex is too close to call
  {
    cell.vertices[0] = 1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.nvert = 2;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 2);
    assert(up == 1);
  }

  /* PATH 2.0 */
  // the first vertex is below the plane and its first edge is above the plane
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 1.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.nvert = 2;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.edges[0] = 1;
    cell.edgeindices[0] = 0;
    cell.edges[3] = 0;
    cell.edgeindices[3] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 1);
    assert(us == 0);
    assert(uw == -1);
    assert(u == 0.25f);
    assert(lp == 0);
    assert(ls == 0);
    assert(qw == 1);
    assert(l == -0.75f);
  }

  /* PATH 2.1 */
  // the first vertex is below the plane and its second edge is above the plane
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -2.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.edges[0] = 1;
    cell.edges[1] = 2;
    cell.edges[6] = 0;
    cell.edgeindices[1] = 0;
    cell.edgeindices[6] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 2);
    assert(us == 0);
    assert(uw == -1);
    assert(u == 0.25f);
    assert(lp == 0);
    assert(ls == 1);
    assert(qw == 1);
    assert(l == -0.75f);
  }

  /* PATH 2.2 */
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -2.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 2;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 2;
    cell.offsets[2] = 5;
    cell.edges[0] = 1;
    cell.edges[1] = 2;
    cell.edges[5] = 0;
    cell.edgeindices[1] = 0;
    cell.edgeindices[5] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 2);
    assert(us == 0);
    assert(uw == -1);
    assert(u == 0.25f);
    assert(lp == 0);
    assert(ls == 1);
    assert(qw == 1);
    assert(l == -0.75f);
  }

  /* PATH 2.3 */
  // the first vertex is below the plane and is the closest vertex to the plane
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -2.0f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = -2.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.offsets[0] = 0;
    cell.edges[0] = 1;
    cell.edges[1] = 2;
    cell.edges[2] = 3;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 0);
  }

  /* PATH 2.4.0 */
  // first vertex is below the plane, second is closer and third is above
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.edges[0] = 1;
    cell.edges[3] = 2;
    cell.edges[5] = 0;
    cell.edges[6] = 1;
    cell.edgeindices[0] = 2;
    cell.edgeindices[3] = 0;
    cell.edgeindices[5] = 0;
    cell.edgeindices[6] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 2);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(l == -0.5f);
  }

  /* PATH 2.4.1 */
  // first vertex is below, second is closer and fourth is above
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = 1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 4;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 7;
    cell.offsets[3] = 10;
    cell.edges[0] = 1;
    cell.edges[3] = 2;
    cell.edges[4] = 3;
    cell.edges[6] = 0;
    cell.edges[10] = 1;
    cell.edgeindices[0] = 3;
    cell.edgeindices[3] = 0;
    cell.edgeindices[4] = 0;
    cell.edgeindices[6] = 0;
    cell.edgeindices[10] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 3);
    assert(us == 0);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(ls == 1);
    assert(l == -0.5f);
  }

  /* PATH 2.4.2 */
  // first vertex is below, second is closer and fourth is above
  // same as 2.4.1, but with order 3 second vertex
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = 1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.offsets[3] = 9;
    cell.edges[0] = 1;
    cell.edges[3] = 2;
    cell.edges[4] = 3;
    cell.edges[5] = 0;
    cell.edges[9] = 1;
    cell.edgeindices[0] = 3;
    cell.edgeindices[3] = 0;
    cell.edgeindices[4] = 0;
    cell.edgeindices[5] = 0;
    cell.edgeindices[9] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 3);
    assert(us == 0);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(ls == 1);
    assert(l == -0.5f);
  }

  /* PATH 2.4.3 */
  // first vertex is below, second is closer, third is above
  // first vertex is first edge of second
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = 1.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edges[6] = 1;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[4] = 0;
    cell.edgeindices[6] = 1;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 2);
    assert(us == 0);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(ls == 1);
    assert(l == -0.5f);
  }

  /* PATH 2.4.4 */
  // first vertex is below, second is closer, fourth is above
  // first vertex is first edge of second
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = 1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 4;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 7;
    cell.offsets[3] = 10;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edges[5] = 3;
    cell.edges[10] = 1;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[5] = 0;
    cell.edgeindices[10] = 2;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 3);
    assert(us == 0);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(ls == 2);
    assert(l == -0.5f);
  }

  /* PATH 2.4.5 */
  // first vertex is below, second is closer, fourth is above
  // first vertex is first edge of second
  // second vertex is order 3 vertex (and not order 4 like 2.4.4)
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.vertices[9] = 1.0f;
    cell.vertices[10] = 0.0f;
    cell.vertices[11] = 0.0f;
    cell.nvert = 4;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.orders[2] = 3;
    cell.orders[3] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 6;
    cell.offsets[3] = 9;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edges[5] = 3;
    cell.edges[9] = 1;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[5] = 0;
    cell.edgeindices[9] = 2;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 1);
    assert(up == 3);
    assert(us == 0);
    assert(u == 0.25f);
    assert(lp == 1);
    assert(ls == 2);
    assert(l == -0.5f);
  }

  /* PATH 2.4.6 */
  // first vertex is below, second is closer and is closest
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = -0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.vertices[6] = -2.0f;
    cell.vertices[7] = 0.0f;
    cell.vertices[8] = 0.0f;
    cell.nvert = 3;
    cell.orders[0] = 3;
    cell.orders[1] = 2;
    cell.orders[2] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.offsets[2] = 5;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edges[4] = 2;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    cell.edgeindices[4] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 0);
  }

  /* PATH 2.5 */
  // first vertex is below, second is too close to call
  {
    cell.vertices[0] = -1.0f;
    cell.vertices[1] = 0.0f;
    cell.vertices[2] = 0.0f;
    cell.vertices[3] = 0.5f;
    cell.vertices[4] = 0.0f;
    cell.vertices[5] = 0.0f;
    cell.nvert = 2;
    cell.orders[0] = 3;
    cell.orders[1] = 3;
    cell.offsets[0] = 0;
    cell.offsets[1] = 3;
    cell.edges[0] = 1;
    cell.edges[3] = 0;
    cell.edgeindices[0] = 0;
    cell.edgeindices[3] = 0;
    dx[0] = 0.5f;
    dx[1] = 0.0f;
    dx[2] = 0.0f;
    r2 = 0.25f;
    int result = voronoi_intersect_find_closest_vertex(
        &cell, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
    assert(result == 2);
  }
}

#ifdef SHADOWFAX_SPH
void set_coordinates(struct part *p, double x, double y, double z,
                     unsigned int id) {
  p->x[0] = x;
  p->x[1] = y;
  p->x[2] = z;
  p->id = id;
  voronoi_cell_init(&p->cell, p->x);
}
#endif

void test_degeneracies() {
#ifdef SHADOWFAX_SPH
  int idx = 0;
  /* make a small cube */
  struct part particles[100];
  set_coordinates(&particles[idx], 0.1, 0.1, 0.1, idx);
  idx++;
  set_coordinates(&particles[idx], 0.2, 0.1, 0.1, idx);
  idx++;
  set_coordinates(&particles[idx], 0.1, 0.2, 0.1, idx);
  idx++;
  set_coordinates(&particles[idx], 0.1, 0.1, 0.2, idx);
  idx++;
  /* corner on cutting plane */
  set_coordinates(&particles[idx], 0.2, 0.2, 0.2, idx);
  idx++;
  /* edge on cutting plane */
  set_coordinates(&particles[idx], 0.2, 0.1, 0.2, idx);
  idx++;
  set_coordinates(&particles[idx], 0.2, 0.2, 0.1, idx);
  idx++;
  /* cutting plane is diagonal */
  set_coordinates(&particles[idx], 0.05, 0.1, 0.05, idx);
  idx++;
  /* order 4 vertex (found after an impressive display of analytical geometry
     of which I'm rather proud) */
  float t = 0.5 / 0.0475;
  set_coordinates(&particles[idx], 0.0075 * t + 0.1, 0.0075 * t + 0.1,
                  0.1 - 0.0025 * t, idx);
  idx++;
  /* order 4 vertex with float edge */
  t = 0.35 / 0.06125;
  set_coordinates(&particles[idx], 0.0075 * t + 0.1, 0.015 * t + 0.1,
                  0.1 - 0.005 * t, idx);
  idx++;
  /* plane that was already encountered */
  t = 0.5 / 0.0475;
  set_coordinates(&particles[idx], 0.0075 * t + 0.1, 0.0075 * t + 0.1,
                  0.1 - 0.0025 * t, idx);
  idx++;
  /* no intersection (just to cover all code) */
  set_coordinates(&particles[idx], 0.3, 0.3, 0.3, idx);
  idx++;
  set_coordinates(&particles[idx], 0.3, 0.1, 0.3, idx);
  idx++;
  /* order 5 vertex */
  t = 0.04 / 0.0175;
  set_coordinates(&particles[idx], 0.1 - 0.0075 * t, 0.1 + 0.00375 * t,
                  0.1 + 0.00625 * t, idx);
  idx++;
  /* plane with order 5 vertex */
  set_coordinates(&particles[idx], 0.1, 0.2, 0.1, idx);
  idx++;
  /* edge with order 5 vertex that looses an edge */
  t = -0.1 / 0.095;
  set_coordinates(&particles[idx], 0.1 - 0.015 * t, 0.1 + 0.015 * t,
                  0.1 - 0.005 * t, idx);
  idx++;
  for (int i = 1; i < idx; i++) {
    float dx[3];
    dx[0] = particles[0].x[0] - particles[i].x[0];
    dx[1] = particles[0].x[1] - particles[i].x[1];
    dx[2] = particles[0].x[2] - particles[i].x[2];
    voronoi_cell_interact(&particles[0].cell, dx, particles[i].id);
  }
#endif
}

VORONOI3D_DECLARE_GLOBAL_VARIABLES()

int main() {

  /* Set the all enclosing simulation box dimensions */
  float box_anchor[3] = {0.0f, 0.0f, 0.0f};
  float box_side[3] = {1.0f, 1.0f, 1.0f};
  voronoi_set_box(box_anchor, box_side);

  /* Check basic Voronoi cell functions */
  test_voronoi_volume_tetrahedron();
  test_voronoi_centroid_tetrahedron();
  test_calculate_cell();

  /* Test the different paths */
  test_paths();

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
  float expected_midpoints[6][3], expected_areas[6];
  expected_areas[0] = 0.25f;
  expected_midpoints[0][0] = 0.25f;
  expected_midpoints[0][1] = 0.5f;
  expected_midpoints[0][2] = 0.5f;
  expected_areas[1] = 0.25f;
  expected_midpoints[1][0] = 0.75f;
  expected_midpoints[1][1] = 0.5f;
  expected_midpoints[1][2] = 0.5f;
  expected_areas[2] = 0.25f;
  expected_midpoints[2][0] = 0.5f;
  expected_midpoints[2][1] = 0.25f;
  expected_midpoints[2][2] = 0.5f;
  expected_areas[3] = 0.25f;
  expected_midpoints[3][0] = 0.5f;
  expected_midpoints[3][1] = 0.75f;
  expected_midpoints[3][2] = 0.5f;
  expected_areas[4] = 0.25f;
  expected_midpoints[4][0] = 0.5f;
  expected_midpoints[4][1] = 0.5f;
  expected_midpoints[4][2] = 0.25f;
  expected_areas[5] = 0.25f;
  expected_midpoints[5][0] = 0.5f;
  expected_midpoints[5][1] = 0.5f;
  expected_midpoints[5][2] = 0.75f;

  /* Interact with some more neighbours to check if they are properly ignored */
  float xE0[3] = {0.6f, 0.0f, 0.1f};
  float xE1[3] = {-0.7f, 0.2f, 0.04f};
  voronoi_cell_interact(&cell, xE0, 7);
  voronoi_cell_interact(&cell, xE1, 8);

  /* Finalize cell and check results */
  voronoi_cell_finalize(&cell);

  if (fabs(cell.volume - 0.125f) > 1.e-5) {
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
  for (int i = 0; i < 6; ++i) {
    A = voronoi_get_face(&cell, i + 1, midpoint);
    if (A) {
      if (fabs(A - expected_areas[i]) > 1.e-5) {
        error("Wrong surface area: %g!", A);
      }
      if (fabs(midpoint[0] - expected_midpoints[i][0] + cell.x[0]) > 1.e-5 ||
          fabs(midpoint[1] - expected_midpoints[i][1] + cell.x[1]) > 1.e-5 ||
          fabs(midpoint[2] - expected_midpoints[i][2] + cell.x[2]) > 1.e-5) {
        error("Wrong face midpoint: %g %g %g (should be %g %g %g)!",
              midpoint[0], midpoint[1], midpoint[2], expected_midpoints[i][0],
              expected_midpoints[i][1], expected_midpoints[i][2]);
      }
    } else {
      error("Neighbour %i not found!", i);
    }
  }

  test_degeneracies();

  return 0;
}
