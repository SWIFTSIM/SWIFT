/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 SWIFT developers.
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
#include "cell_boundary.h"
#include "swift.h"

#include <assert.h>
#include <string.h>

int main(void) {
  struct space s;
  struct engine e;
  struct cell c;
  struct gpart gparts[3];
  bzero(&s, sizeof(s));
  bzero(&e, sizeof(e));
  bzero(&c, sizeof(c));
  bzero(gparts, sizeof(gparts));

  e.s = &s;
  e.nodeID = 0;
  s.dim[0] = s.dim[1] = s.dim[2] = 1.;
  const int lock_error = lock_init(&s.lock);
  assert(lock_error == 0);
  c.nodeID = 0;
  c.grav.count = 3;
  c.grav.parts = gparts;
  for (int i = 0; i < 3; i++) {
    gparts[i].type = swift_type_dark_matter;
    gparts[i].time_bin = 1;
  }

  gparts[0].x[0] = -0.1;
  gparts[1].x[0] = 1.; /* Exact upper face remains inside. */
  gparts[2].x[0] = 0.5;
  cell_boundary_gpart(&c, &e);
  assert(gparts[0].time_bin == time_bin_inhibited);
  assert(gparts[1].time_bin != time_bin_inhibited);
  assert(s.nr_inhibited_gparts == 1);

  gparts[1].x[0] = 1.1;
  cell_boundary_gpart(&c, &e);
  assert(gparts[1].time_bin == time_bin_inhibited);
  assert(s.nr_inhibited_gparts == 2); /* No double removal. */

  /* Reflective lower and upper x walls; all other faces stay open. */
  s.boundary_conditions.type = boundary_condition_reflective;
  s.boundary_conditions.walls = boundary_wall_x_min | boundary_wall_x_max;
  bzero(gparts, sizeof(gparts));
  for (int i = 0; i < 3; i++) {
    gparts[i].type = swift_type_dark_matter;
    gparts[i].time_bin = 1;
    gparts[i].x[1] = gparts[i].x[2] = 0.5;
  }

  /* Crossed the lower wall: mirrored with its normal velocity flipped. */
  gparts[0].x[0] = -0.125;
  gparts[0].v_full[0] = -1.f;
  gparts[0].v_full[1] = 2.f;
  /* Crossed the upper wall. */
  gparts[1].x[0] = 1.25;
  gparts[1].v_full[0] = 1.f;
  /* Crossed a wall and an open face: removed. */
  gparts[2].x[0] = -0.125;
  gparts[2].x[1] = -0.125;
  cell_boundary_gpart(&c, &e);
  assert(gparts[0].x[0] == 0.125 && gparts[0].v_full[0] == 1.f);
  assert(gparts[0].x[1] == 0.5 && gparts[0].v_full[1] == 2.f);
  assert(gparts[0].time_bin != time_bin_inhibited);
  assert(gparts[1].x[0] == 0.75 && gparts[1].v_full[0] == -1.f);
  assert(gparts[1].time_bin != time_bin_inhibited);
  assert(gparts[2].time_bin == time_bin_inhibited);
  assert(s.nr_inhibited_gparts == 3);

  /* A second pass leaves the reflected particles alone. */
  cell_boundary_gpart(&c, &e);
  assert(gparts[0].x[0] == 0.125 && gparts[0].v_full[0] == 1.f);
  assert(s.nr_inhibited_gparts == 3);
  return 0;
}
