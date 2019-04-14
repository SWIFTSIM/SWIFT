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
#include "align.h"
#include "parser.h"

/* Constants. */
#define UNION_BY_SIZE_OVER_MPI (1)
#define FOF_COMPRESS_PATHS_MIN_LENGTH (2)
#define FOF_NO_GAS (-1)
#define FOF_BLACK_HOLE (-2)
#define FOF_LOW_HALO_MASS (-3)

/* Avoid cyclic inclusions */
struct space;

/* MPI message required for FOF. */
struct fof_mpi {

  /* The local particle's root ID.*/
  size_t group_i;

  /* The local group's size.*/
  size_t group_i_size;

  /* The foreign particle's root ID.*/
  size_t group_j;

  /* The local group's size.*/
  size_t group_j_size;

} SWIFT_STRUCT_ALIGN;

struct fof {

  size_t *group_index;
  size_t *group_size;
  double *group_mass;
  long long *max_part_density_index;
  float *max_part_density;
  
  /*! The extra no. of black holes to seed locally. */
  int extra_bh_seed_count;
  
  /*! The FOF linking length squared. */
  double l_x2;
  
  /*! The minimum halo mass for black hole seeding. */
  double seed_halo_mass;
  
  /*! The no. of steps between each FOF search. */
  int run_freq;
  
  int num_groups;
  size_t min_group_size;
  size_t group_id_default;
  size_t group_id_offset;
  int group_links_size_default;

  int group_link_count;
  struct fof_mpi *group_links;
  int group_links_size;

  char base_name[PARSER_MAX_LINE_SIZE];

} SWIFT_STRUCT_ALIGN;

/* Store group size and offset into array. */
struct group_length {

  size_t index, size;

} SWIFT_STRUCT_ALIGN;

#ifdef WITH_MPI
/* Struct used to find final group ID when using MPI */
struct fof_final_index {
  size_t local_root;
  size_t global_root;
} SWIFT_STRUCT_ALIGN;

/* Struct used to find the total mass of a group when using MPI */
struct fof_final_mass {
  size_t global_root;
  double group_mass;
  long long max_part_density_index;
  float max_part_density;
} SWIFT_STRUCT_ALIGN;

/* Struct used to iterate over the hash table and unpack the mass fragments of a group when using MPI */
struct fof_mass_send_hashmap {
  struct fof_final_mass *mass_send;
  size_t nsend; 
} SWIFT_STRUCT_ALIGN;
#endif

/* Store local and foreign cell indices that touch. */
struct cell_pair_indices {

  struct cell *local, *foreign;

} SWIFT_STRUCT_ALIGN;

/* Function prototypes. */
void fof_init(struct space *s);
void fof_search_cell(struct space *s, struct cell *c);
void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj);
void fof_search_pair_cells_foreign(struct space *s, struct cell *ci,
                                   struct cell *cj, int *link_count,
                                   struct fof_mpi **group_links,
                                   int *group_links_size);
void fof_search_tree(struct space *s);
void fof_dump_group_data(char *out_file, struct space *s,
                         int num_groups, struct group_length *group_sizes);
void rec_fof_search_self(struct cell *c, struct space *s, const double dim[3],
                         const double search_r2);
void rec_fof_search_pair(struct cell *restrict ci, struct cell *restrict cj,
                         struct space *s, const double dim[3],
                         const double search_r2);

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype fof_mpi_type;
extern MPI_Datatype group_length_mpi_type;
#endif

#endif /* SWIFT_FOF_H */
