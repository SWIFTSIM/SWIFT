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
#include <config.h>

/* Local headers */
#include "align.h"
#include "parser.h"
#include "part_type.h"

/* Avoid cyclic inclusions */
struct cell;
struct gpart;
struct space;
struct engine;
struct unit_system;
struct phys_const;
struct black_holes_props;
struct cosmology;

/* Store group size and offset into array. */
struct group_length {

  size_t index, size;

} SWIFT_STRUCT_ALIGN;

struct fof_props {

  /*! Whether we're doing periodic FoF calls to seed black holes. */
  int seed_black_holes_enabled;

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

  /*! The types of particles to use for linking */
  int fof_linking_types[swift_type_count];

  /*! The types of particles to use for attaching */
  int fof_attach_types[swift_type_count];

  /* ------------  Group properties ----------------- */

  /*! Number of groups */
  long long num_groups;

  /*! Index of the root particle of the group a given gpart belongs to. */
  size_t *group_index;

  /*! Index of the root particle of the group a given gpart is attached to. */
  size_t *attach_index;

  /*! Has the particle found a linkable to attach to? */
  char *found_attachable_link;

  /*! For attachable particles: distance to the current nearest linkable part */
  float *distance_to_link;

  /*! Size of the group a given gpart belongs to. */
  size_t *group_size;

  /*! Size of the local groups a given gpart belongs to. */
  struct group_length *high_group_sizes;

  /*! Final size of the group a given gpart belongs to. */
  long long *final_group_size;

  /*! Final index of the group a given gpart belongs to. */
  long long *final_group_index;

  /*! Mass of the group a given gpart belongs to. */
  double *group_mass;

  /*! Gas mass of the group a given gpart belongs to. */
  float *group_gas_mass;

  /*! Stellar mass of the group a given gpart belongs to. */
  float *group_stellar_mass;

  /*! SFR of the group a given gpart belongs to. */
  float *group_star_formation_rate;

  /*! Does the group have a black hole? */
  char *has_black_hole;

  /*! Centre of mass of the group a given gpart belongs to. */
  double *group_centre_of_mass;

  /*! Distance to the particle furthest from the centre of mass in
   * the group a given gpart belongs to. */
  float *group_radii;

  /*! Maximal density of all parts of each group. */
  float *max_part_density;

  /* Number of groups on each node */
  size_t *num_on_node;

  /* First group on each node */
  size_t *first_on_node;

  /* Total number of groups on lower numbered MPI ranks */
  size_t num_groups_prev;

  /* ------------ MPI-related arrays --------------- */

  /*! The number of links between pairs of particles on this node and
   * a foreign node */
  int group_link_count;

  /*! The allocated size of the links array */
  int group_links_size;

  /*! The links between pairs of particles on this node and a foreign
   * node */
  struct fof_mpi *group_links;
};

#ifdef WITH_MPI

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
};

/* Struct used to find final group ID when using MPI */
struct fof_final_index {
  size_t local_root;
  size_t global_root;
};

/* Store local and foreign cell indices that touch. */
struct cell_pair_indices {
  struct cell *local, *foreign;
};
#endif

/* Function prototypes. */
void fof_init(struct fof_props *props, struct swift_params *params,
              const struct phys_const *phys_const, const struct unit_system *us,
              const int stand_alone_fof);
void fof_create_mpi_types(void);
void fof_allocate(const struct space *s, struct fof_props *props);
void fof_compute_local_sizes(struct fof_props *props, struct space *s);
void fof_search_foreign_cells(struct fof_props *props, const struct space *s);
void fof_link_attachable_particles(struct fof_props *props,
                                   const struct space *s);
void fof_finalise_attachables(struct fof_props *props, struct space *s);
void fof_link_foreign_fragments(struct fof_props *props, const struct space *s);
void fof_assign_group_ids(struct fof_props *props, struct space *s);
void fof_compute_group_props(struct fof_props *props,
                             const struct black_holes_props *bh_props,
                             const struct phys_const *constants,
                             const struct cosmology *cosmo, struct space *s,
                             const int dump_results,
                             const int dump_debug_results,
                             const int seed_black_holes);
void rec_fof_search_self(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *c);
void rec_fof_search_pair(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *restrict ci, struct cell *restrict cj);
void rec_fof_attach_self(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         const size_t nr_gparts, struct cell *c);
void rec_fof_attach_pair(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         const size_t nr_gparts, struct cell *restrict ci,
                         struct cell *restrict cj, const int ci_local,
                         const int cj_local);
void fof_free_arrays(struct fof_props *props);
void fof_struct_dump(const struct fof_props *props, FILE *stream);
void fof_struct_restore(struct fof_props *props, FILE *stream);

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype fof_mpi_type;
extern MPI_Datatype group_length_mpi_type;
#endif

#endif /* SWIFT_FOF_H */
