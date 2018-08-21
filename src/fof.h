/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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

#ifndef SWIFT_FOF_H
#define SWIFT_FOF_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "cell.h"
#include "space.h"

/* MPI message required for FOF. */
struct fof_mpi {
  
  /* The local particle's root ID.*/
  int group_i;

  /* The local group's size.*/
  int group_i_size;

  /* The local group's mass.*/
  float group_i_mass;
  
  /* The foreign particle's root ID.*/
  int group_j;

} SWIFT_STRUCT_ALIGN; 

/* Function prototypes. */
void fof_init(struct space *s, long long Ngas, long long Ngparts);
void fof_search_serial(struct space *s);
void fof_search_cell(struct space *s, struct cell *c);
void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj);
void fof_search_pair_cells_foreign(struct space *s, struct cell *ci, struct cell *cj, size_t *link_count, struct fof_mpi *part_links);
void fof_search_tree_serial(struct space *s);
void fof_search_tree(struct space *s);
void fof_dump_group_data(char *out_file, const size_t nr_gparts, int *group_index, int *num_in_groups, long long *group_id, double *group_mass, const int min_group_size);

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype fof_mpi_type;
#endif

#endif /* SWIFT_FOF_H */
