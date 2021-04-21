/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#ifndef SWIFT_LIGHTCONE_MAP_TYPES_H
#define SWIFT_LIGHTCONE_MAP_TYPES_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "parser.h"
#include "part_type.h"

/* Avoid cyclic inclusions */
struct engine;
struct lightcone_map;
struct gpart;

/* Type to store pointer to function for updating a healpix map */
typedef void (*map_update_function_t)(struct lightcone_map *map, const struct engine *e,
                                      const struct gpart *gp, const double a_cross,
                                      const double x_cross[3]);

/**
 * @brief Struct to store information on one type of lightcone map
 */
struct lightcone_map_type {
  char name[PARSER_MAX_LINE_SIZE];
  map_update_function_t update_map;
};


void lightcone_map_total_mass(struct lightcone_map *map, const struct engine *e,
                              const struct gpart *gp, const double a_cross,
                              const double x_cross[3]);


/* This associates map names to the appropriate update function */
static const struct lightcone_map_type lightcone_map_types[] = {
  {"TotalMass", lightcone_map_total_mass},
  {"",          NULL}, /* NULL function indicates end of array */
};

#endif
