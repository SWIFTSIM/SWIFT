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

#include <stdio.h>

#ifndef SWIFT_PERIODIC_REPLICATIONS_H
#define SWIFT_PERIODIC_REPLICATIONS_H

/* Forward declarations */
struct cell;

/* Struct to store information about one periodic replication of the simulation box */
struct replication {

  /* Minimum distance squared from the observer to any point in the replication */
  double rmin2;

  /* Maximum distance squared from the observer to any point in the replication */
  double rmax2;

  /* Coordinates of the replication */
  double coord[3];
};


/* Struct to store an array of periodic replications  */
struct replication_list {

  /* Number of replications*/
  int nrep;

  /* Distance limits used to make this replication list */
  double lightcone_rmin;
  double lightcone_rmax;

  /* Array of replications with nrep elements */
  struct replication *replication;
};

void replication_list_init(struct replication_list *replication_list,
                           double boxsize, double cell_width,
                           double observer_position[3],
                           double lightcone_rmin, double lightcone_rmax);

void replication_list_clean(struct replication_list *replication_list);

void replication_list_write(struct replication_list *replication_list, FILE *fd);

void replication_list_subset_for_cell(const struct replication_list *rep_in,
                                      const struct cell *cell,
                                      const double observer_position[3],
                                      struct replication_list *rep_out);

#endif
