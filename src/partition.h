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

#include "parser.h"
#include "space.h"
#include "task.h"

/* Initial partitioning types. */
enum partition_type {
  INITPART_GRID = 0,
  INITPART_VECTORIZE,
  INITPART_METIS_WEIGHT,
  INITPART_METIS_NOWEIGHT
};

/* Simple descriptions of types for reports. */
extern const char *initial_partition_name[];

/* The selected initial partition type and any related metadata. */
struct partition {
  enum partition_type type;
  int grid[3];
  int usemetis;
};

/* Repartition type to use. */
enum repartition_type {
  REPART_NONE = 0,
  REPART_METIS_VERTEX_EDGE_COSTS,
  REPART_METIS_EDGE_COSTS,
  REPART_METIS_VERTEX_COUNTS,
  REPART_METIS_VERTEX_COSTS_TIMEBINS
};

/* Repartition preferences. */
struct repartition {
  enum repartition_type type;
  float trigger;
  float minfrac;
  float itr;
  int usemetis;
  int adaptive;

  int use_fixed_costs;
  int use_ticks;

  /* The partition as a cell-list. */
  int ncelllist;
  int *celllist;
};

/* Simple descriptions of types for reports. */
extern const char *repartition_name[];

void partition_repartition(struct repartition *reparttype, int nodeID,
                           int nr_nodes, struct space *s, struct task *tasks,
                           int nr_tasks);
void partition_initial_partition(struct partition *initial_partition,
                                 int nodeID, int nr_nodes, struct space *s);

int partition_space_to_space(double *oldh, double *oldcdim, int *oldnodeID,
                             struct space *s);
void partition_init(struct partition *partition,
                    struct repartition *repartition,
                    struct swift_params *params, int nr_nodes);

/* Dump/restore. */
void partition_store_celllist(struct space *s, struct repartition *reparttype);
void partition_restore_celllist(struct space *s,
                                struct repartition *reparttype);
void partition_struct_dump(struct repartition *reparttype, FILE *stream);
void partition_struct_restore(struct repartition *reparttype, FILE *stream);

#endif /* SWIFT_PARTITION_H */
