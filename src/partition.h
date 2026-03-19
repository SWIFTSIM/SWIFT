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

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <limits.h>

/* Local headers. */
#include "parser.h"
#include "space.h"
#include "task.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>

/* METIS/ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif

#endif

#ifndef IDX_MAX
/* When compiled without METIS/ParMETIS support, we need to define IDX_MAX
 * for accumulate_sizes which can be used by the wedge decomposition. */
#define IDX_MAX LLONG_MAX
#endif

/* Forward declarations. */
struct cell;

/* Initial partitioning types. */
enum partition_type {
  INITPART_GRID = 0,
  INITPART_VECTORIZE,
  INITPART_METIS_WEIGHT,
  INITPART_METIS_NOWEIGHT,
  INITPART_METIS_WEIGHT_EDGE,
  INITPART_WEDGE,
};

/* Simple descriptions of types for reports. */
extern const char *initial_partition_name[];

/* The selected initial partition type and any related metadata. */
struct partition {
  enum partition_type type;
  int grid[3];
  int usemetis;
  int nr_wedges;
  int nr_theta_slices;
  int nr_phi_slices;
  double theta_width;
  double phi_width;
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

void partition_clean(struct partition *partition,
                     struct repartition *repartition);

/* Partition type specific functions (common between zoom and uniform
 * partitions). */
void pick_vector(const int cdim[3], const int nregions, int *samplecells);
void split_vector(struct cell *cells_top, const int cdim[3], const int nregions,
                  int *samplecells);

/* Accumulate the total particle counts in each cell for weighting. */
void partition_accumulate_sizes(struct space *s, int verbose, double *counts);

/* Dump/restore. */
void partition_store_celllist(struct space *s, struct repartition *reparttype);
void partition_restore_celllist(struct space *s,
                                struct repartition *reparttype);
void partition_struct_dump(struct repartition *reparttype, FILE *stream);
void partition_struct_restore(struct repartition *reparttype, FILE *stream);

/* Metis/ParMetis support. */
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
struct weights_mapper_data {
  double *weights_e;
  double *weights_v;
  idx_t *inds;
  idx_t *xadj;
  int nr_edges;
  int eweights;
  int nodeID;
  int timebins;
  int vweights;
  int nr_cells;
  int use_ticks;
  struct cell *cells;
};

extern double repartition_costs[task_type_count][task_subtype_count];

void partition_sizes_to_edges(struct space *s, double *counts, double *edges,
                              const int *cell_edge_offsets);
int partition_count_edges(struct space *s, int periodic, int verbose,
                          int *cell_edge_offsets);
void partition_graph_init(struct space *s, int periodic, idx_t *adjncy,
                          int *nadjcny, idx_t *xadj, int *nxadj,
                          const int *cell_edge_offsets);
void partition_pick_metis(int nodeID, struct space *s, int nregions,
                          double *vertexw, double *edgew, int *celllist,
                          const int *cell_edge_offsets, int nedges);
void partition_pick_parmetis(int nodeID, struct space *s, int nregions,
                             double *vertexw, double *edgew, int refine,
                             int adaptive, float itr, int *celllist,
                             const int *cell_edge_offsets, int nedges);
void partition_split_metis(struct space *s, int nregions, int *celllist);
#endif

/* Debugging. */
#ifdef SWIFT_DEBUG_CHECKS
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
void partition_check_weights(struct task *tasks, int nr_tasks,
                             struct weights_mapper_data *mydata,
                             double *ref_weights_v, double *ref_weights_e);
#endif
#endif

#endif /* SWIFT_PARTITION_H */
