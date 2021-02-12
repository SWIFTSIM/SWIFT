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

#ifndef SWIFT_PERIODIC_REPLICATIONS_H
#define SWIFT_PERIODIC_REPLICATIONS_H

/* Config parameters. */
#include "../config.h"


/* Struct to store information about one periodic replication of the simulation box */
struct replication {

  /* Minimum distance from the observer to any point in the replication */
  double rmin;

  /* Maximum distance from the observer to any point in the replication */
  double rmax;

  /* Integer coordinates of the replication */
  int coord[3];
};


/* Struct to store an array of periodic replications  */
struct replication_list {

  /* Number of replications*/
  int nrep;

  /* Array of replications with nrep elements */
  struct replication *replication;
}

#endif
