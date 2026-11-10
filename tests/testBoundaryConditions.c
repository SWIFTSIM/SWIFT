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
#include "boundary_conditions.h"
#include "parser.h"

#include <assert.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[]) {
  struct swift_params params;
  struct boundary_conditions bc;

  assert(fabs(boundary_conditions_wrap_periodic(-0.1, 1.) - 0.9) < 1e-15);
  assert(boundary_conditions_wrap_periodic(1.25, 1.) == 0.25);
  assert(boundary_conditions_wrap_periodic(1., 1.) == 0.);

  bzero(&params, sizeof(params));
  parser_set_param(&params, "InitialConditions:periodic:1");
  boundary_conditions_init(&bc, &params);
  assert(bc.type == boundary_condition_periodic && bc.walls == 0);

  bzero(&params, sizeof(params));
  parser_set_param(&params, "InitialConditions:periodic:0");
  boundary_conditions_init(&bc, &params);
  assert(bc.type == boundary_condition_open && bc.walls == 0);

  bzero(&params, sizeof(params));
  parser_set_param(&params, "BoundaryConditions:type:periodic");
  boundary_conditions_init(&bc, &params);
  assert(bc.type == boundary_condition_periodic && bc.walls == 0);

  bzero(&params, sizeof(params));
  parser_set_param(&params, "BoundaryConditions:type:open");
  boundary_conditions_init(&bc, &params);
  assert(bc.type == boundary_condition_open && bc.walls == 0);

  bzero(&params, sizeof(params));
  parser_set_param(&params, "InitialConditions:periodic:0");
  parser_set_param(&params, "BoundaryConditions:type:reflective");
  parser_set_param(&params, "BoundaryConditions:walls:[x_min, x_max, z_min]");
  boundary_conditions_init(&bc, &params);
  assert(bc.type == boundary_condition_reflective);
  assert(bc.walls ==
         (boundary_wall_x_min | boundary_wall_x_max | boundary_wall_z_min));

  if (argc > 1) {
    bzero(&params, sizeof(params));
    for (int i = 1; i < argc; i++) parser_set_param(&params, argv[i]);
    boundary_conditions_init(&bc, &params);
  }

  return 0;
}
