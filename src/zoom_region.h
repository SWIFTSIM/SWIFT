/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "cell.h"
#include "gravity_properties.h"
#include "partition.h"
#include "space.h"

#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif

#ifndef SWIFT_ZOOM_H
#define SWIFT_ZOOM_H

void zoom_region_init(struct swift_params *params, struct space *s,
                      struct gravity_props *gravity_properties,
                      int verbose);
int cell_getid_zoom(const struct space *s, const double x, const double y,
                    const double z);
void construct_zoom_region(struct space *s, int nr_nodes, int verbose);
void construct_tl_cells_with_zoom_region(
    struct space *s, const int *cdim, const float dmin,
    const integertime_t ti_current, struct gravity_props *gravity_properties,
    int verbose);
void find_void_cells(struct space *s, const int verbose);
void find_neighbouring_cells(struct space *s,
                             struct gravity_props *gravity_properties,
                             const int verbose);
void link_zoom_to_void(struct space *s, struct cell *c);
void find_vertex_edges(struct space *s, const int verbose);
double cell_min_dist2_diff_size(const struct cell *restrict ci,
                                const struct cell *restrict cj,
                                int periodic, const double dim[3]);
double cell_min_dist2(const struct cell *restrict ci,
                      const struct cell *restrict cj, int periodic,
                      const double dim[3]);
void engine_makeproxies_with_zoom_region(struct engine *e);
void engine_make_self_gravity_tasks_mapper_natural_cells(void *map_data,
                                                         int num_elements,
                                                         void *extra_data);
void engine_make_self_gravity_tasks_mapper_buffer_cells(void *map_data,
                                                        int num_elements,
                                                        void *extra_data);
void engine_make_self_gravity_tasks_mapper_zoom_cells(void *map_data,
                                                      int num_elements,
                                                      void *extra_data);
void engine_make_self_gravity_tasks_mapper_zoom_bkg(void *map_data,
                                                    int num_elements,
                                                    void *extra_data);
void engine_make_self_gravity_tasks_mapper_buffer_bkg(void *map_data,
                                                      int num_elements,
                                                      void *extra_data);
void engine_make_hydroloop_tasks_mapper_with_zoom(void *map_data,
                                                  int num_elements,
                                                  void *extra_data);
void engine_make_drift_tasks_for_wanderers(struct engine *e, struct cell *c);
void engine_make_drift_tasks_for_wanderers_mapper(void *map_data,
                                                  int num_elements,
                                                  void *extra_data);
void engine_make_fofloop_tasks_mapper_with_zoom(void *map_data,
                                                int num_elements,
                                                void *extra_data);
/* Parition prototypes */
void partition_initial_partition_zoom(struct partition *initial_partition,
                                      int nodeID, int nr_nodes,
                                      struct space *s);
void partition_repartition_zoom(struct repartition *reparttype, int nodeID,
                                int nr_nodes, struct space *s,
                                struct task *tasks,
                                int nr_tasks);
int partition_space_to_space_zoom(double *oldh, double *oldcdim,
                                  double *oldzoomh, double *oldzoomcdim,
                                  int *oldnodeIDs, struct space *s);

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
void edge_loop(const int *cdim, int offset, struct space *s,
               idx_t *adjncy, idx_t *xadj, double *counts, double *edges,
               int *iedge);
int get_wedge_index(struct space *s, struct cell *c);
void graph_init_zoom(struct space *s, int periodic, idx_t *weights_e,
                     idx_t *adjncy, int *nadjcny, idx_t *xadj,
                     int *nxadj, int nverts, int offset, int *cdim);
void sizes_to_edges_zoom(struct space *s, double *counts, double *edges,
                         int offset, int *cdim);
void split_metis_zoom(struct space *s, int nregions, int *celllist, int ncells,
                      int offset);
#endif
#ifdef WITH_MPI
void engine_addtasks_recv_void(struct engine *e);
void engine_addtasks_send_void(struct engine *e);
void engine_addtasks_recv_zoom_gravity(struct engine *e, struct cell *c,
                                       struct cell *zoom_c,
                                       struct task *t_grav,
                                       struct task *const tend);
void engine_addtasks_send_zoom_gravity(struct engine *e, struct cell *ci,
                                       struct cell *cj,
                                       struct task *t_grav,
                                       int tag);
#endif
void activate_void_tasks(struct engine *e);

void void_count_send_gparts(struct cell *c, struct engine *e, size_t *count,
                           int nodeID);
void void_count_recv_gparts(struct cell *c, struct engine *e, size_t *count,
                           int nodeID);
int void_attach_send_gparts(struct cell *c, struct engine *e, int count,
                            struct gpart *buff, int nodeID);
/* Regrid prototypes */

void space_regrid_zoom(struct space *s, struct gravity_props *p,
                       int nr_nodes, int verbose);

#endif /* SWIFT_ZOOM_H */
