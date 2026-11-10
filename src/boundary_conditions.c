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
#include <config.h>

#include "boundary_conditions.h"

#include <string.h>

#include "error.h"
#include "parser.h"

/**
 * @brief Read the box boundary conditions from the parameter file.
 *
 * @param bc Boundary conditions to initialise.
 * @param params Simulation parameters.
 */
void boundary_conditions_init(struct boundary_conditions *bc,
                              struct swift_params *params) {

  const int has_legacy =
      parser_does_param_exist(params, "InitialConditions:periodic");
  const int has_type = parser_does_param_exist(params, "BoundaryConditions:type");
  const int has_walls =
      parser_does_param_exist(params, "BoundaryConditions:walls");

  if (!has_legacy && !has_type)
    error("Specify BoundaryConditions:type or InitialConditions:periodic.");

  int legacy_periodic = 0;
  if (has_legacy) {
    legacy_periodic = parser_get_param_int(params, "InitialConditions:periodic");
    if (legacy_periodic != 0 && legacy_periodic != 1)
      error("InitialConditions:periodic must be 0 or 1.");
  }

  char type[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(params, "BoundaryConditions:type", type,
                              legacy_periodic ? "periodic" : "open");

  if (!strcmp(type, "periodic"))
    bc->type = boundary_condition_periodic;
  else if (!strcmp(type, "open"))
    bc->type = boundary_condition_open;
  else if (!strcmp(type, "reflective"))
    bc->type = boundary_condition_reflective;
  else
    error("Unknown BoundaryConditions:type '%s'.", type);

  if (has_legacy && legacy_periodic != (bc->type == boundary_condition_periodic))
    error("BoundaryConditions:type conflicts with InitialConditions:periodic.");

  bc->walls = 0;
  if (bc->type != boundary_condition_reflective) {
    if (has_walls)
      error("BoundaryConditions:walls requires type: reflective.");
    return;
  }

  if (!has_walls) error("Reflective boundaries require BoundaryConditions:walls.");

  static const char *const names[6] = {"x_min", "x_max", "y_min",
                                       "y_max", "z_min", "z_max"};
  char **walls;
  int nr_walls;
  parser_get_param_string_array(params, "BoundaryConditions:walls", &nr_walls,
                                &walls);
  if (nr_walls == 0) error("BoundaryConditions:walls cannot be empty.");

  for (int i = 0; i < nr_walls; i++) {
    int j = 0;
    for (; j < 6; j++)
      if (!strcmp(walls[i], names[j])) break;
    if (j == 6) error("Unknown boundary wall '%s'.", walls[i]);
    if (bc->walls & (1u << j)) error("Duplicate boundary wall '%s'.", walls[i]);
    bc->walls |= 1u << j;
  }

  parser_free_param_string_array(nr_walls, walls);
}
