/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_PARAM_H
#define SWIFT_PARAM_H


/* Initial partition grid struct. Defines type of partitioning to use and any
 * related parameters. */
enum grid_types {
  GRID_GRID = 0,
  GRID_RANDOM,
  GRID_VECTORIZE,
  GRID_METIS_WEIGHT,
  GRID_METIS_NOWEIGHT
};

struct pgrid {
  enum grid_types type;
  int grid[3];
};


#endif /* SWIFT_PARAM_H */
