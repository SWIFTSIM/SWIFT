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

/* Avoid cyclic inclusions */
struct gpart;
struct space;
struct engine;
struct unit_system;
struct phys_const;
struct black_holes_props;
struct cosmology;

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

struct fof_props {

  /* ----------- Parameters of the FOF search ------- */

  /*! The linking length in units of the mean DM inter-particle separation. */
  double l_x_ratio;

  /*! The absolute linking length in internal units if the user wants to
   *   overwrite the one based on the mean inter-particle separation. */
  double l_x_absolute;

  /*! The square of the linking length. */
  double l_x2;

  /*! The minimum halo mass for black hole seeding. */
  double seed_halo_mass;

  /*! Minimal number of particles in a group */
  size_t min_group_size;

  /*! Default group ID to give to particles not in a group */
  size_t group_id_default;

  /*! ID of the first (largest) group. */
  size_t group_id_offset;

  /*! The base name of the output file */
  char base_name[PARSER_MAX_LINE_SIZE];

  /* ------------  Group properties ----------------- */

  /*! Number of groups */
  int num_groups;

  /*! Number of local black holes that belong to groups whose roots are on a
   * different node. */
  int extra_bh_seed_count;

  /*! Index of the root particle of the group a given gpart belongs to. */
  size_t *group_index;

  /*! Size of the group a given gpart belongs to. */
  size_t *group_size;

  /*! Mass of the group a given gpart belongs to. */
  double *group_mass;

  /*! Index of the part with the maximal density of each group. */
  long long *max_part_density_index;

  /*! Maximal density of all parts of each group. */
  float *max_part_density;

  /* ------------ MPI-related arrays --------------- */

  /*! The number of links between pairs of particles on this node and
   * a foreign node */
  int group_link_count;

  /*! The allocated size of the links array */
  int group_links_size;

  /*! The links between pairs of particles on this node and a foreign
   * node */
  struct fof_mpi *group_links;

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

/* Struct used to iterate over the hash table and unpack the mass fragments of a
 * group when using MPI */
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
void fof_init(struct fof_props *props, struct swift_params *params,
              const struct phys_const *phys_const, const struct unit_system *us,
              const int stand_alone_fof);
void fof_create_mpi_types(void);
void fof_allocate(const struct space *s, const long long total_nr_DM_particles,
                  struct fof_props *props);
void fof_search_tree(struct fof_props *props,
                     const struct black_holes_props *bh_props,
                     const struct phys_const *constants,
                     const struct cosmology *cosmo, struct space *s,
                     const int dump_results, const int seed_black_holes);
void rec_fof_search_self(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *c);
void rec_fof_search_pair(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *restrict ci, struct cell *restrict cj);
void fof_struct_dump(const struct fof_props *props, FILE *stream);
void fof_struct_restore(struct fof_props *props, FILE *stream);
#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype fof_mpi_type;
extern MPI_Datatype group_length_mpi_type;
#endif

#endif /* SWIFT_FOF_H */
