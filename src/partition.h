/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_PARTITION_H
#define SWIFT_PARTITION_H

#include "space.h"
#include "task.h"

/* Initial partitioning types. */
enum initpart_type {
  INITPART_GRID = 0,
  INITPART_VECTORIZE,
  INITPART_METIS_WEIGHT,
  INITPART_METIS_NOWEIGHT
};

/* Simple descriptions of types for reports. */
extern const char *initpart_name[];

/* The selected initial partition type and any related metadata. */
struct initpart {
  enum initpart_type type;
  int grid[3];
};

/* Repartition type to use. */
enum repart_type {
  REPART_NONE = 0,
  REPART_METIS_BOTH,
  REPART_METIS_VERTEX,
  REPART_METIS_EDGE,
  REPART_METIS_VERTEX_EDGE
};

/* Simple descriptions of types for reports. */
extern const char *repart_name[];

void part_pick_metis(struct space *s, int nregions, int *vertexw, int *edgew,
                     int *celllist);

void part_repart(enum repart_type reparttype, int nodeID, int nr_nodes, 
                 struct space *s, struct task *tasks, int nr_tasks);
void part_part(struct initpart *ipart, int nodeID, int nr_nodes,
               struct space *s);
int part_check_complete(struct space *s, int verbose, int nregions);

#endif /* SWIFT_PARTITION_H */
