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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <errno.h>
#include <libgen.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "fof.h"

/* Local headers. */
#include "black_holes.h"
#include "common_io.h"
#include "engine.h"
#include "hashmap.h"
#include "memuse.h"
#include "proxy.h"
#include "threadpool.h"

#define fof_props_default_run_freq 2000
#define fof_props_default_group_id 2147483647
#define fof_props_default_group_id_offset 1
#define fof_props_default_group_link_size 20000

/* Constants. */
#define UNION_BY_SIZE_OVER_MPI (1)
#define FOF_COMPRESS_PATHS_MIN_LENGTH (2)

/**
 * @brief Properties of a group used for black hole seeding
 */
enum fof_halo_seeding_props {
  fof_halo_has_no_gas = -1LL,
  fof_halo_has_black_hole = -2LL,
  fof_halo_has_too_low_mass = -3LL
};

#ifdef WITH_MPI

/* MPI types used for communications */
MPI_Datatype fof_mpi_type;
MPI_Datatype group_length_mpi_type;
MPI_Datatype fof_final_index_type;
MPI_Datatype fof_final_mass_type;

/*! Offset between the first particle on this MPI rank and the first particle in
 * the global order */
size_t node_offset;
#endif

#ifdef SWIFT_DEBUG_CHECKS
static integertime_t ti_current;
#endif

/**
 * @brief Initialise the properties of the FOF code.
 *
 * @param props the #fof_props structure to fill.
 * @param params the parameter file parser.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 */
void fof_init(struct fof_props *props, struct swift_params *params,
              const struct phys_const *phys_const,
              const struct unit_system *us) {

  /* Base name for the FOF output file */
  parser_get_param_string(params, "FOF:basename", props->base_name);

  /* Check that we can write outputs by testing if the output
   * directory exists and is searchable and writable. */
  const char *dirp = dirname(props->base_name);
  if (access(dirp, W_OK | X_OK) != 0) {
    error("Cannot write FOF outputs in directory %s (%s)", dirp,
          strerror(errno));
  }

  /* Read the FOF search frequency. */
  props->run_freq = parser_get_opt_param_int(params, "FOF:run_freq",
                                             fof_props_default_run_freq);

  /* Read the minimum group size. */
  props->min_group_size = parser_get_param_int(params, "FOF:min_group_size");

  /* Read the default group ID of particles in groups below the minimum group
   * size. */
  props->group_id_default = parser_get_opt_param_int(
      params, "FOF:group_id_default", fof_props_default_group_id);

  /* Read the starting group ID. */
  props->group_id_offset = parser_get_opt_param_int(
      params, "FOF:group_id_offset", fof_props_default_group_id_offset);

  /* Read the linking length ratio to the mean inter-particle separation. */
  props->l_x_ratio =
      parser_get_param_double(params, "FOF:linking_length_ratio");

  if (props->l_x_ratio <= 0.)
    error("The FOF linking length can't be negative!");

  /* Read value of absolute linking length aksed by the user */
  props->l_x_absolute =
      parser_get_opt_param_double(params, "FOF:absolute_linking_length", -1.);

  if (props->l_x_ratio <= 0. && props->l_x_ratio != -1.)
    error("The FOF linking length can't be negative!");

  /* Read the minimal halo mass for black hole seeding */
  props->seed_halo_mass =
      parser_get_param_double(params, "FOF:black_hole_seed_halo_mass_Msun");

  /* Convert to internal units */
  props->seed_halo_mass *= phys_const->const_solar_mass;

#if defined(WITH_MPI) && defined(UNION_BY_SIZE_OVER_MPI)
  if (engine_rank == 0)
    message(
        "Performing FOF over MPI using union by size and union by rank "
        "locally.");
#else
  message("Performing FOF using union by rank.");
#endif
}

/**
 * @brief Registers MPI types used by FOF.
 */
void fof_create_mpi_types() {

#ifdef WITH_MPI
  if (MPI_Type_contiguous(sizeof(struct fof_mpi) / sizeof(unsigned char),
                          MPI_BYTE, &fof_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof.");
  }
  if (MPI_Type_contiguous(sizeof(struct group_length) / sizeof(unsigned char),
                          MPI_BYTE, &group_length_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&group_length_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for group_length.");
  }
  /* Define type for sending fof_final_index struct */
  if (MPI_Type_contiguous(sizeof(struct fof_final_index), MPI_BYTE,
                          &fof_final_index_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_final_index_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof_final_index.");
  }
  /* Define type for sending fof_final_mass struct */
  if (MPI_Type_contiguous(sizeof(struct fof_final_mass), MPI_BYTE,
                          &fof_final_mass_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_final_mass_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof_final_mass.");
  }
#else
  error("Calling an MPI function in non-MPI code.");
#endif
}

/**
 * @brief Allocate the memory and initialise the arrays for a FOF calculation.
 *
 * @param s The #space to act on.
 * @param total_nr_DM_particles The total number of DM particles in the
 * simulation.
 * @param props The properties of the FOF structure.
 */
void fof_allocate(const struct space *s, const long long total_nr_DM_particles,
                  struct fof_props *props) {

  /* Calculate the particle linking length based upon the mean inter-particle
   * spacing of the DM particles. */
  const double mean_inter_particle_sep =
      s->dim[0] / cbrt((double)total_nr_DM_particles);
  const double l_x = props->l_x_ratio * mean_inter_particle_sep;

  /* Are we using the aboslute value or the one derived from the mean
     inter-particle sepration? */
  if (props->l_x_absolute != -1.) {
    props->l_x2 = props->l_x_absolute * props->l_x_absolute;
  } else {
    props->l_x2 = l_x * l_x;
  }

#ifdef WITH_MPI
  /* Check size of linking length against the top-level cell dimensions. */
  if (props->l_x2 > s->width[0] * s->width[0])
    error(
        "Linking length greater than the width of a top-level cell. Need to "
        "check more than one layer of top-level cells for links.");
#endif

  const size_t nr_local_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;

  /* Allocate and initialise a group index array. */
  if (swift_memalign("fof_group_index", (void **)&props->group_index, 64,
                     nr_local_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle group indices for FOF search.");

  /* Allocate and initialise a group size array. */
  if (swift_memalign("fof_group_size", (void **)&props->group_size, 64,
                     nr_local_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of group size for FOF search.");

  /* Set initial group ID of the gparts */
  const size_t group_id_default = props->group_id_default;
  for (size_t i = 0; i < nr_local_gparts; i++) {
    gparts[i].group_id = group_id_default;
  }

  /* Set initial group index and group size */
  size_t *group_index = props->group_index;
  size_t *group_size = props->group_size;
  for (size_t i = 0; i < nr_local_gparts; i++) {
    group_index[i] = i;
    group_size[i] = 1;
  }

#ifdef SWIFT_DEBUG_CHECKS
  ti_current = s->e->ti_current;
#endif
}

/**
 * @brief Comparison function for qsort call comparing group sizes.
 *
 * @param a The first #group_length object.
 * @param b The second #group_length object.
 * @return 1 if the size of the group b is larger than the size of group a, -1
 * if a is the largest and 0 if they are equal.
 */
int cmp_func_group_size(const void *a, const void *b) {
  struct group_length *a_group_size = (struct group_length *)a;
  struct group_length *b_group_size = (struct group_length *)b;
  if (b_group_size->size > a_group_size->size)
    return 1;
  else if (b_group_size->size < a_group_size->size)
    return -1;
  else
    return 0;
}

#ifdef WITH_MPI

/**
 * @brief Comparison function for qsort call comparing group global roots.
 *
 * @param a The first #fof_final_index object.
 * @param b The second #fof_final_index object.
 * @return 1 if the global of the group b is *smaller* than the global group of
 * group a, -1 if a is the smaller one and 0 if they are equal.
 */
int compare_fof_final_index_global_root(const void *a, const void *b) {
  struct fof_final_index *fof_final_index_a = (struct fof_final_index *)a;
  struct fof_final_index *fof_final_index_b = (struct fof_final_index *)b;
  if (fof_final_index_b->global_root < fof_final_index_a->global_root)
    return 1;
  else if (fof_final_index_b->global_root > fof_final_index_a->global_root)
    return -1;
  else
    return 0;
}

/**
 * @brief Comparison function for qsort call comparing group global roots
 *
 * @param a The first #fof_final_mass object.
 * @param b The second #fof_final_mass object.
 * @return 1 if the global of the group b is *smaller* than the global group of
 * group a, -1 if a is the smaller one and 0 if they are equal.
 */
int compare_fof_final_mass_global_root(const void *a, const void *b) {
  struct fof_final_mass *fof_final_mass_a = (struct fof_final_mass *)a;
  struct fof_final_mass *fof_final_mass_b = (struct fof_final_mass *)b;
  if (fof_final_mass_b->global_root < fof_final_mass_a->global_root)
    return 1;
  else if (fof_final_mass_b->global_root > fof_final_mass_a->global_root)
    return -1;
  else
    return 0;
}
#endif

/**
 * @brief Check whether a given group ID is on the local node.
 *
 * This function only makes sense in MPI mode.
 *
 * @param group_id The ID to check.
 * @param nr_gparts The number of gparts on this node.
 */
__attribute__((always_inline)) INLINE static int is_local(
    const size_t group_id, const size_t nr_gparts) {
#ifdef WITH_MPI
  return (group_id >= node_offset && group_id < node_offset + nr_gparts);
#else
  error("Calling MPI function in non-MPI mode");
  return 1;
#endif
}

/**
 * @brief Find the global root ID of a given particle
 *
 * This function only makes sense in MPI mode.
 *
 * @param i Index of the particle.
 * @param group_index Array of group root indices.
 * @param nr_gparts The number of g-particles on this node.
 */
__attribute__((always_inline)) INLINE static size_t fof_find_global(
    const size_t i, const size_t *group_index, const size_t nr_gparts) {

#ifdef WITH_MPI
  size_t root = node_offset + i;
  if (!is_local(root, nr_gparts)) {

    /* Non local --> This is the root */
    return root;
  } else {

    /* Local --> Follow the links until we find the root */
    while (root != group_index[root - node_offset]) {
      root = group_index[root - node_offset];
      if (!is_local(root, nr_gparts)) break;
    }
  }

  /* Perform path compression. */
  // int index = i;
  // while(index != root) {
  //  int next = group_index[index];
  //  group_index[index] = root;
  //  index = next;
  //}

  return root;
#else
  error("Calling MPI function in non-MPI mode");
  return -1;
#endif
}

/**
 * @brief   Finds the local root ID of the group a particle exists in
 * when group_index contains globally unique identifiers -
 * i.e. we stop *before* we advance to a foreign root.
 *
 * Here we assume that the input i is a local index and we
 * return the local index of the root.
 *
 * @param i Index of the particle.
 * @param nr_gparts The number of g-particles on this node.
 * @param group_index Array of group root indices.
 */
__attribute__((always_inline)) INLINE static size_t fof_find_local(
    const size_t i, const size_t nr_gparts, const size_t *group_index) {
#ifdef WITH_MPI
  size_t root = node_offset + i;

  while ((group_index[root - node_offset] != root) &&
         (group_index[root - node_offset] >= node_offset) &&
         (group_index[root - node_offset] < node_offset + nr_gparts)) {
    root = group_index[root - node_offset];
  }

  return root - node_offset;
#else
  size_t root = i;

  while ((group_index[root] != root) && (group_index[root] < nr_gparts)) {
    root = group_index[root];
  }

  return root;
#endif
}

/**
 * @brief Finds the local root ID of the group a particle exists in.
 *
 * We follow the group_index array until reaching the root of the group.
 *
 * Also performs path compression if the path is long.
 *
 * @param i The index of the particle.
 * @param group_index Array of group root indices.
 */
__attribute__((always_inline)) INLINE static size_t fof_find(
    const size_t i, size_t *group_index) {

  size_t root = i;
  int tree_depth = 0;

  while (root != group_index[root]) {
#ifdef PATH_HALVING
    atomic_cas(&group_index[root], group_index[root],
               group_index[group_index[root]]);
#endif
    root = group_index[root];
    tree_depth++;
  }

  /* Only perform path compression on trees with a depth of
   * FOF_COMPRESS_PATHS_MIN_LENGTH or higher. */
  if (tree_depth >= FOF_COMPRESS_PATHS_MIN_LENGTH)
    atomic_cas(&group_index[i], group_index[i], root);

  return root;
}

/**
 * @brief Atomically update the root of a group
 *
 * @param address The address of the value to update.
 * @param y The new value to write.
 *
 * @return 1 If successful, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int atomic_update_root(
    volatile size_t *address, const size_t y) {

  size_t *size_t_ptr = (size_t *)address;

  size_t old_val = *address;
  size_t test_val = old_val;
  size_t new_val = y;

  /* atomic_cas returns old_val if *size_t_ptr has not changed since being
   * read.*/
  old_val = atomic_cas(size_t_ptr, test_val, new_val);

  if (test_val == old_val)
    return 1;
  else
    return 0;
}

/**
 * @brief Unifies two groups by setting them to the same root.
 *
 * @param root_i The root of the first group. Will be updated.
 * @param root_j The root of the second group.
 * @param group_index The list of group roots.
 */
__attribute__((always_inline)) INLINE static void fof_union(
    size_t *root_i, const size_t root_j, size_t *group_index) {

  int result = 0;

  /* Loop until the root can be set to a new value. */
  do {
    size_t root_i_new = fof_find(*root_i, group_index);
    const size_t root_j_new = fof_find(root_j, group_index);

    /* Skip particles in the same group. */
    if (root_i_new == root_j_new) return;

    /* If the root ID of pj is lower than pi's root ID set pi's root to point to
     * pj's. Otherwise set pj's root to point to pi's.*/
    if (root_j_new < root_i_new) {

      /* Updates the root and checks that its value has not been changed since
       * being read. */
      result = atomic_update_root(&group_index[root_i_new], root_j_new);

      /* Update root_i on the fly. */
      *root_i = root_j_new;
    } else {

      /* Updates the root and checks that its value has not been changed since
       * being read. */
      result = atomic_update_root(&group_index[root_j_new], root_i_new);

      /* Update root_i on the fly. */
      *root_i = root_i_new;
    }
  } while (result != 1);
}

/**
 * @brief Compute th minimal distance between any two points in two cells.
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param dim The size of the simulation domain.
 */
__attribute__((always_inline)) INLINE static double cell_min_dist(
    const struct cell *restrict ci, const struct cell *restrict cj,
    const double dim[3]) {

  /* Get cell locations. */
  const double cix_min = ci->loc[0];
  const double ciy_min = ci->loc[1];
  const double ciz_min = ci->loc[2];
  const double cjx_min = cj->loc[0];
  const double cjy_min = cj->loc[1];
  const double cjz_min = cj->loc[2];

  const double cix_max = ci->loc[0] + ci->width[0];
  const double ciy_max = ci->loc[1] + ci->width[1];
  const double ciz_max = ci->loc[2] + ci->width[2];
  const double cjx_max = cj->loc[0] + cj->width[0];
  const double cjy_max = cj->loc[1] + cj->width[1];
  const double cjz_max = cj->loc[2] + cj->width[2];

  double not_same_range[3];

  /* If two cells are in the same range of coordinates along
     any of the 3 axis, the distance along this axis is 0 */
  if (ci->width[0] > cj->width[0]) {
    if ((cix_min <= cjx_min) && (cjx_max <= cix_max))
      not_same_range[0] = 0.;
    else
      not_same_range[0] = 1.;
  } else {
    if ((cjx_min <= cix_min) && (cix_max <= cjx_max))
      not_same_range[0] = 0.;
    else
      not_same_range[0] = 1.;
  }
  if (ci->width[1] > cj->width[1]) {
    if ((ciy_min <= cjy_min) && (cjy_max <= ciy_max))
      not_same_range[1] = 0.;
    else
      not_same_range[1] = 1.;
  } else {
    if ((cjy_min <= ciy_min) && (ciy_max <= cjy_max))
      not_same_range[1] = 0.;
    else
      not_same_range[1] = 1.;
  }
  if (ci->width[2] > cj->width[2]) {
    if ((ciz_min <= cjz_min) && (cjz_max <= ciz_max))
      not_same_range[2] = 0.;
    else
      not_same_range[2] = 1.;
  } else {
    if ((cjz_min <= ciz_min) && (ciz_max <= cjz_max))
      not_same_range[2] = 0.;
    else
      not_same_range[2] = 1.;
  }

  /* Find the shortest distance between cells, remembering to account for
   * periodic boundary conditions. */
  double dx[3];
  dx[0] = min4(fabs(nearest(cix_min - cjx_min, dim[0])),
               fabs(nearest(cix_min - cjx_max, dim[0])),
               fabs(nearest(cix_max - cjx_min, dim[0])),
               fabs(nearest(cix_max - cjx_max, dim[0])));

  dx[1] = min4(fabs(nearest(ciy_min - cjy_min, dim[1])),
               fabs(nearest(ciy_min - cjy_max, dim[1])),
               fabs(nearest(ciy_max - cjy_min, dim[1])),
               fabs(nearest(ciy_max - cjy_max, dim[1])));

  dx[2] = min4(fabs(nearest(ciz_min - cjz_min, dim[2])),
               fabs(nearest(ciz_min - cjz_max, dim[2])),
               fabs(nearest(ciz_max - cjz_min, dim[2])),
               fabs(nearest(ciz_max - cjz_max, dim[2])));

  double r2 = 0.;
  for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k] * not_same_range[k];

  return r2;
}

#ifdef WITH_MPI
/* Add a group to the hash table. */
__attribute__((always_inline)) INLINE static void hashmap_add_group(
    const size_t group_id, const size_t group_offset, hashmap_t *map) {

  int created_new_element = 0;
  hashmap_value_t *offset =
      hashmap_get_new(map, group_id, &created_new_element);

  if (offset != NULL) {

    /* If the element is a new entry set its value. */
    if (created_new_element) {
      (*offset).value_st = group_offset;
    }
  } else
    error("Couldn't find key (%zu) or create new one.", group_id);
}

/* Find a group in the hash table. */
__attribute__((always_inline)) INLINE static size_t hashmap_find_group_offset(
    const size_t group_id, hashmap_t *map) {

  hashmap_value_t *group_offset = hashmap_get(map, group_id);

  if (group_offset == NULL)
    error("Couldn't find key (%zu) or create new one.", group_id);

  return (size_t)(*group_offset).value_st;
}

/* Compute send/recv offsets for MPI communication. */
__attribute__((always_inline)) INLINE static void fof_compute_send_recv_offsets(
    const int nr_nodes, int *sendcount, int **recvcount, int **sendoffset,
    int **recvoffset, size_t *nrecv) {

  /* Determine number of entries to receive */
  *recvcount = malloc(nr_nodes * sizeof(int));
  MPI_Alltoall(sendcount, 1, MPI_INT, *recvcount, 1, MPI_INT, MPI_COMM_WORLD);

  /* Compute send/recv offsets */
  *sendoffset = malloc(nr_nodes * sizeof(int));

  (*sendoffset)[0] = 0;
  for (int i = 1; i < nr_nodes; i++)
    (*sendoffset)[i] = (*sendoffset)[i - 1] + sendcount[i - 1];

  *recvoffset = malloc(nr_nodes * sizeof(int));

  (*recvoffset)[0] = 0;
  for (int i = 1; i < nr_nodes; i++)
    (*recvoffset)[i] = (*recvoffset)[i - 1] + (*recvcount)[i - 1];

  /* Allocate receive buffer */
  *nrecv = 0;
  for (int i = 0; i < nr_nodes; i++) (*nrecv) += (*recvcount)[i];
}
#endif

/**
 * @brief Perform a FOF search using union-find on a given leaf-cell
 *
 * @param props The properties fof the FOF scheme.
 * @param l_x2 The square of the FOF linking length.
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param c The #cell in which to perform FOF.
 */
void fof_search_self_cell(const struct fof_props *props, const double l_x2,
                          const struct gpart *const space_gparts,
                          struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->split) error("Performing the FOF search at a non-leaf level!");
#endif

  const size_t count = c->grav.count;
  struct gpart *gparts = c->grav.parts;

  /* Index of particles in the global group list */
  size_t *group_index = props->group_index;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset = group_index + (ptrdiff_t)(gparts - space_gparts);

  if (c->nodeID != engine_rank)
    error("Performing self FOF search on foreign cell.");

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    struct gpart *pi = &gparts[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Find the root of pi. */
    size_t root_i = fof_find(offset[i], group_index);

    for (size_t j = i + 1; j < count; j++) {

      struct gpart *pj = &gparts[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Compute the pairwise distance */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Merge the groups */
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

/**
 * @brief Perform a FOF search using union-find between two cells
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the simulation volume.
 * @param l_x2 The square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 */
void fof_search_pair_cells(const struct fof_props *props, const double dim[3],
                           const double l_x2, const int periodic,
                           const struct gpart *const space_gparts,
                           struct cell *restrict ci, struct cell *restrict cj) {

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart *gparts_i = ci->grav.parts;
  struct gpart *gparts_j = cj->grav.parts;

  /* Index of particles in the global group list */
  size_t *group_index = props->group_index;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset_i = group_index + (ptrdiff_t)(gparts_i - space_gparts);
  size_t *const offset_j = group_index + (ptrdiff_t)(gparts_j - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS
  if (offset_j > offset_i && (offset_j < offset_i + count_i))
    error("Overlapping cells");
  if (offset_i > offset_j && (offset_i < offset_j + count_j))
    error("Overlapping cells");
#endif

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  double diff[3];
  for (int k = 0; k < 3; k++) {
    diff[k] = cj->loc[k] - ci->loc[k];
    if (periodic && diff[k] < -dim[k] * 0.5)
      shift[k] = dim[k];
    else if (periodic && diff[k] > dim[k] * 0.5)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
    diff[k] += shift[k];
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    struct gpart *pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Find the root of pi. */
    size_t root_i = fof_find(offset_i[i], group_index);

    for (size_t j = 0; j < count_j; j++) {

      struct gpart *pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset_j[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance, remembering to account for boundary
       * conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Merge the groups */
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

/* Perform a FOF search between a local and foreign cell using the Union-Find
 * algorithm. Store any links found between particles.*/
void fof_search_pair_cells_foreign(
    const struct fof_props *props, const double dim[3], const double l_x2,
    const int periodic, const struct gpart *const space_gparts,
    const size_t nr_gparts, const struct cell *restrict ci,
    const struct cell *restrict cj, int *restrict link_count,
    struct fof_mpi **group_links, int *restrict group_links_size) {

#ifdef WITH_MPI
  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  const struct gpart *gparts_i = ci->grav.parts;
  const struct gpart *gparts_j = cj->grav.parts;

  /* Get local pointers */
  size_t *group_index = props->group_index;
  size_t *group_size = props->group_size;

  /* Values local to this function to avoid dereferencing */
  struct fof_mpi *local_group_links = *group_links;
  int local_link_count = *link_count;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset_i = group_index + (ptrdiff_t)(gparts_i - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS

  /* Check whether cells are local to the node. */
  const int ci_local = (ci->nodeID == engine_rank);
  const int cj_local = (cj->nodeID == engine_rank);

  if ((ci_local && cj_local) || (!ci_local && !cj_local))
    error(
        "FOF search of foreign cells called on two local cells or two foreign "
        "cells.");

  if (!ci_local) {
    error("Cell ci, is not local.");
  }
#endif

  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  for (int k = 0; k < 3; k++) {
    const double diff = cj->loc[k] - ci->loc[k];
    if (periodic && diff < -dim[k] / 2)
      shift[k] = dim[k];
    else if (periodic && diff > dim[k] / 2)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    const struct gpart *pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Find the root of pi. */
    const size_t root_i =
        fof_find_global(offset_i[i] - node_offset, group_index, nr_gparts);

    for (size_t j = 0; j < count_j; j++) {

      const struct gpart *pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance, remembering to account for boundary
       * conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        int found = 0;

        /* Check that the links have not already been added to the list. */
        for (int l = 0; l < local_link_count; l++) {
          if ((local_group_links)[l].group_i == root_i &&
              (local_group_links)[l].group_j == pj->group_id) {
            found = 1;
            break;
          }
        }

        if (!found) {

          /* If the group_links array is not big enough re-allocate it. */
          if (local_link_count + 1 > *group_links_size) {

            const int new_size = 2 * (*group_links_size);

            *group_links_size = new_size;

            (*group_links) = (struct fof_mpi *)realloc(
                *group_links, new_size * sizeof(struct fof_mpi));

            /* Reset the local pointer */
            local_group_links = *group_links;

            message("Re-allocating local group links from %d to %d elements.",
                    local_link_count, new_size);
          }

          /* Store the particle group properties for communication. */
          local_group_links[local_link_count].group_i = root_i;
          local_group_links[local_link_count].group_i_size =
              group_size[root_i - node_offset];

          local_group_links[local_link_count].group_j = pj->group_id;
          local_group_links[local_link_count].group_j_size = pj->group_size;

          local_link_count++;
        }
      }
    }
  }

  /* Update the returned values */
  *link_count = local_link_count;

#else
  error("Calling MPI function in non-MPI mode.");
#endif
}

/**
 * @brief Recursively perform a union-find FOF between two cells.
 *
 * If cells are more distant than the linking length, we abort early.
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the space.
 * @param search_r2 the square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 */
void rec_fof_search_pair(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *restrict ci, struct cell *restrict cj) {

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Pair FOF called on same cell!!!");
#endif

  /* Return if cells are out of range of each other. */
  if (r2 > search_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts,
                                ci->progeny[k], cj->progeny[l]);
      }
    }
  } else if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts,
                            ci->progeny[k], cj);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts, ci,
                            cj->progeny[k]);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells(props, dim, search_r2, periodic, space_gparts, ci,
                          cj);
  }
}
#ifdef WITH_MPI

/* Recurse on a pair of cells (one local, one foreign) and perform a FOF search
 * between cells that are within range. */
void rec_fof_search_pair_foreign(
    const struct fof_props *props, const double dim[3], const double search_r2,
    const int periodic, const struct gpart *const space_gparts,
    const size_t nr_gparts, const struct cell *ci, const struct cell *cj,
    int *restrict link_count, struct fof_mpi **group_links,
    int *restrict group_links_size) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Pair FOF called on same cell!!!");
#endif

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

  /* Return if cells are out of range of each other. */
  if (r2 > search_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                        space_gparts, nr_gparts, ci->progeny[k],
                                        cj->progeny[l], link_count, group_links,
                                        group_links_size);
      }
    }
  } else if (ci->split) {

    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                    space_gparts, nr_gparts, ci->progeny[k], cj,
                                    link_count, group_links, group_links_size);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                    space_gparts, nr_gparts, ci, cj->progeny[k],
                                    link_count, group_links, group_links_size);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells_foreign(props, dim, search_r2, periodic, space_gparts,
                                  nr_gparts, ci, cj, link_count, group_links,
                                  group_links_size);
  }
}

#endif

/**
 * @brief Recursively perform a union-find FOF on a cell.
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the space.
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param search_r2 the square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param c The #cell in which to perform FOF.
 */
void rec_fof_search_self(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *c) {

  /* Recurse? */
  if (c->split) {

    /* Loop over all progeny. Perform pair and self recursion on progenies.*/
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        rec_fof_search_self(props, dim, search_r2, periodic, space_gparts,
                            c->progeny[k]);

        for (int l = k + 1; l < 8; l++)
          if (c->progeny[l] != NULL)
            rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts,
                                c->progeny[k], c->progeny[l]);
      }
    }
  }
  /* Otherwise, compute self-interaction. */
  else
    fof_search_self_cell(props, search_r2, space_gparts, c);
}

/* Mapper function to atomically update the group size array. */
void fof_update_group_size_mapper(hashmap_key_t key, hashmap_value_t *value,
                                  void *data) {

  size_t *group_size = (size_t *)data;

  /* Use key to index into group size array. */
  atomic_add(&group_size[key], value->value_st);
}

/**
 * @brief Mapper function to calculate the group sizes.
 *
 * @param map_data An array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_calc_group_size_mapper(void *map_data, int num_elements,
                                void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct gpart *gparts = (struct gpart *)map_data;
  size_t *group_index = s->e->fof_properties->group_index;
  size_t *group_size = s->e->fof_properties->group_size;

  /* Offset into gparts array. */
  ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);
  size_t *const group_index_offset = group_index + gparts_offset;

  /* Create hash table. */
  hashmap_t map;
  hashmap_init(&map);

  /* Loop over particles and find which cells are in range of each other to
   * perform the FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    hashmap_key_t root =
        (hashmap_key_t)fof_find(group_index_offset[ind], group_index);
    const size_t gpart_index = gparts_offset + ind;

    /* Only add particles which aren't the root of a group. Stops groups of size
     * 1 being added to the hash table. */
    if (root != gpart_index) {
      hashmap_value_t *size = hashmap_get(&map, root);

      if (size != NULL)
        (*size).value_st++;
      else
        error("Couldn't find key (%zu) or create new one.", root);
    }
  }

  /* Update the group size array. */
  if (map.size > 0)
    hashmap_iterate(&map, fof_update_group_size_mapper, group_size);

  hashmap_free(&map);
}

/* Mapper function to atomically update the group mass array. */
static INLINE void fof_update_group_mass_mapper(hashmap_key_t key,
                                                hashmap_value_t *value,
                                                void *data) {

  double *group_mass = (double *)data;

  /* Use key to index into group mass array. */
  atomic_add_d(&group_mass[key], value->value_dbl);
}

/**
 * @brief Mapper function to calculate the group masses.
 *
 * @param map_data An array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_calc_group_mass_mapper(void *map_data, int num_elements,
                                void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct gpart *gparts = (struct gpart *)map_data;
  double *group_mass = s->e->fof_properties->group_mass;
  const size_t group_id_default = s->e->fof_properties->group_id_default;
  const size_t group_id_offset = s->e->fof_properties->group_id_offset;

  /* Create hash table. */
  hashmap_t map;
  hashmap_init(&map);

  /* Loop over particles and increment the group mass for groups above
   * min_group_size. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Only check groups above the minimum size. */
    if (gparts[ind].group_id != group_id_default) {

      hashmap_key_t index = gparts[ind].group_id - group_id_offset;
      hashmap_value_t *data = hashmap_get(&map, index);

      /* Update group mass */
      if (data != NULL)
        (*data).value_dbl += gparts[ind].mass;
      else
        error("Couldn't find key (%zu) or create new one.", index);
    }
  }

  /* Update the group mass array. */
  if (map.size > 0)
    hashmap_iterate(&map, fof_update_group_mass_mapper, group_mass);

  hashmap_free(&map);
}

#ifdef WITH_MPI
/* Mapper function to unpack hash table into array. */
void fof_unpack_group_mass_mapper(hashmap_key_t key, hashmap_value_t *value,
                                  void *data) {

  struct fof_mass_send_hashmap *fof_mass_send =
      (struct fof_mass_send_hashmap *)data;
  struct fof_final_mass *mass_send = fof_mass_send->mass_send;
  size_t *nsend = &fof_mass_send->nsend;

  /* Store elements from hash table in array. */
  mass_send[*nsend].global_root = key;
  mass_send[(*nsend)].group_mass = value->value_dbl;
  mass_send[(*nsend)].max_part_density_index = value->value_st;
  mass_send[(*nsend)++].max_part_density = value->value_flt;
}
#endif

/**
 * @brief Calculates the total mass of each group above min_group_size and finds
 * the densest particle for black hole seeding.
 */
void fof_calc_group_mass(struct fof_props *props, const struct space *s,
                         const size_t num_groups_local,
                         const size_t num_groups_prev,
                         size_t *restrict num_on_node,
                         size_t *restrict first_on_node, double *group_mass) {

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;
  const struct part *parts = s->parts;
  const size_t group_id_offset = props->group_id_offset;
  const size_t group_id_default = props->group_id_default;
  const double seed_halo_mass = props->seed_halo_mass;

#ifdef WITH_MPI
  size_t *group_index = props->group_index;
  const int nr_nodes = s->e->nr_nodes;

  /* Allocate and initialise the densest particle array. */
  if (swift_memalign("max_part_density_index",
                     (void **)&props->max_part_density_index, 32,
                     num_groups_local * sizeof(long long)) != 0)
    error(
        "Failed to allocate list of max group density indices for FOF search.");

  if (swift_memalign("max_part_density", (void **)&props->max_part_density, 32,
                     num_groups_local * sizeof(float)) != 0)
    error("Failed to allocate list of max group densities for FOF search.");

  /* Direct pointers to the arrays */
  long long *max_part_density_index = props->max_part_density_index;
  float *max_part_density = props->max_part_density;

  /* No densest particle found so far */
  bzero(max_part_density, num_groups_local * sizeof(float));

  /* Start by assuming that the haloes have no gas */
  for (size_t i = 0; i < num_groups_local; i++) {
    max_part_density_index[i] = fof_halo_has_no_gas;
  }

  /* Start the hash map */
  hashmap_t map;
  hashmap_init(&map);

  /* JSW TODO: Parallelise with threadpool*/
  for (size_t i = 0; i < nr_gparts; i++) {

    /* Check if the particle is in a group above the threshold. */
    if (gparts[i].group_id != group_id_default) {

      const size_t root = fof_find_global(i, group_index, nr_gparts);

      /* Increment the mass of groups that are local */
      if (is_local(root, nr_gparts)) {

        const size_t index =
            gparts[i].group_id - group_id_offset - num_groups_prev;

        /* Update group mass */
        group_mass[index] += gparts[i].mass;

      }
      /* Add mass fragments of groups that have a foreign root to a hash table.
       */
      else {

        hashmap_value_t *data = hashmap_get(&map, (hashmap_key_t)root);

        if (data != NULL) {
          data->value_dbl += gparts[i].mass;

          /* Find the densest gas particle.
           * Account for groups that already have a black hole and groups that
           * contain no gas. */
          if (gparts[i].type == swift_type_gas &&
              data->value_st != fof_halo_has_black_hole) {

            /* Update index if a denser gas particle is found. */
            if (parts[-gparts[i].id_or_neg_offset].rho > data->value_flt) {
              data->value_flt = parts[-gparts[i].id_or_neg_offset].rho;
              data->value_st = -gparts[i].id_or_neg_offset;
            }
          }
          /* If there is already a black hole in the group we don't need to
             create a new one. */
          else if (gparts[i].type == swift_type_black_hole) {
            data->value_st = fof_halo_has_black_hole;
            data->value_flt = 0.f;
          }
        } else
          error("Couldn't find key (%zu) or create new one.", root);
      }
    }
  }

  /* Loop over particles and find the densest particle in each group. */
  /* JSW TODO: Parallelise with threadpool*/
  for (size_t i = 0; i < nr_gparts; i++) {

    /* Only check groups above the minimum size and mass threshold. */
    if (gparts[i].group_id != group_id_default) {

      size_t root = fof_find_global(i, group_index, nr_gparts);

      /* Increment the mass of groups that are local */
      if (is_local(root, nr_gparts)) {

        size_t index = gparts[i].group_id - group_id_offset - num_groups_prev;

        /* Only seed groups above the mass threshold. */
        if (group_mass[index] > seed_halo_mass) {

          /* Find the densest gas particle.
           * Account for groups that already have a black hole and groups that
           * contain no gas. */
          if (gparts[i].type == swift_type_gas &&
              max_part_density_index[index] != fof_halo_has_black_hole) {

            /* Update index if a denser gas particle is found. */
            if (parts[-gparts[i].id_or_neg_offset].rho >
                max_part_density[index]) {
              max_part_density_index[index] = -gparts[i].id_or_neg_offset;
              max_part_density[index] = parts[-gparts[i].id_or_neg_offset].rho;
            }
          }
          /* If there is already a black hole in the group we don't need to
             create a new one. */
          else if (gparts[i].type == swift_type_black_hole) {
            max_part_density_index[index] = fof_halo_has_black_hole;
          }
        }
      }
    }
  }

  size_t nsend = map.size;
  struct fof_mass_send_hashmap hashmap_mass_send;

  /* Allocate and initialise a mass array. */
  if (posix_memalign((void **)&hashmap_mass_send.mass_send, 32,
                     nsend * sizeof(struct fof_mass_send_hashmap)) != 0)
    error("Failed to allocate list of group masses for FOF search.");

  hashmap_mass_send.nsend = 0;

  struct fof_final_mass *fof_mass_send = hashmap_mass_send.mass_send;

  /* Unpack mass fragments and roots from hash table. */
  if (map.size > 0)
    hashmap_iterate(&map, fof_unpack_group_mass_mapper, &hashmap_mass_send);

  nsend = hashmap_mass_send.nsend;

  if (nsend != map.size)
    error("No. of mass fragments to send != elements in hash table.");

  hashmap_free(&map);

  /* Sort by global root - this puts the groups in order of which node they're
   * stored on */
  qsort(fof_mass_send, nsend, sizeof(struct fof_final_mass),
        compare_fof_final_mass_global_root);

  /* Determine how many entries go to each node */
  int *sendcount = malloc(nr_nodes * sizeof(int));
  for (int i = 0; i < nr_nodes; i += 1) sendcount[i] = 0;
  int dest = 0;
  for (size_t i = 0; i < nsend; i += 1) {
    while ((fof_mass_send[i].global_root >=
            first_on_node[dest] + num_on_node[dest]) ||
           (num_on_node[dest] == 0))
      dest += 1;
    if (dest >= nr_nodes) error("Node index out of range!");
    sendcount[dest] += 1;
  }

  int *recvcount = NULL, *sendoffset = NULL, *recvoffset = NULL;
  size_t nrecv = 0;

  fof_compute_send_recv_offsets(nr_nodes, sendcount, &recvcount, &sendoffset,
                                &recvoffset, &nrecv);

  struct fof_final_mass *fof_mass_recv =
      malloc(nrecv * sizeof(struct fof_final_mass));

  /* Exchange group mass */
  MPI_Alltoallv(fof_mass_send, sendcount, sendoffset, fof_final_mass_type,
                fof_mass_recv, recvcount, recvoffset, fof_final_mass_type,
                MPI_COMM_WORLD);

  /* For each received global root, look up the group ID we assigned and
   * increment the group mass */
  for (size_t i = 0; i < nrecv; i += 1) {
    if ((fof_mass_recv[i].global_root < node_offset) ||
        (fof_mass_recv[i].global_root >= node_offset + nr_gparts)) {
      error("Received global root index out of range!");
    }
    group_mass[gparts[fof_mass_recv[i].global_root - node_offset].group_id -
               group_id_offset - num_groups_prev] +=
        fof_mass_recv[i].group_mass;
  }

  /* For each received global root, look up the group ID we assigned and find
   * the global maximum gas density */
  for (size_t i = 0; i < nrecv; i++) {

    const int offset =
        gparts[fof_mass_recv[i].global_root - node_offset].group_id -
        group_id_offset - num_groups_prev;

    /* Only seed groups above the mass threshold. */
    if (group_mass[offset] > seed_halo_mass) {

      /* Only check groups that don't already contain a black hole. */
      if (max_part_density_index[offset] != fof_halo_has_black_hole) {

        /* Find the densest particle in each group using the densest particle
         * from each group fragment. */
        if (fof_mass_recv[i].max_part_density > max_part_density[offset]) {
          max_part_density[offset] = fof_mass_recv[i].max_part_density;
          max_part_density_index[offset] =
              fof_mass_recv[i].max_part_density_index;
        }
      }
      /* If there is already a black hole in the group we don't need to create a
         new one. */
      else if (fof_mass_recv[i].max_part_density_index ==
               fof_halo_has_black_hole) {
        max_part_density_index[offset] = fof_halo_has_black_hole;
      }
    } else {
      max_part_density_index[offset] = fof_halo_has_no_gas;
    }
  }

  /* For each received global root, look up the group ID we assigned and send
   * the global maximum gas density index back */
  for (size_t i = 0; i < nrecv; i++) {
    if ((fof_mass_recv[i].global_root < node_offset) ||
        (fof_mass_recv[i].global_root >= node_offset + nr_gparts)) {
      error("Received global root index out of range!");
    }

    const int offset =
        gparts[fof_mass_recv[i].global_root - node_offset].group_id -
        group_id_offset - num_groups_prev;

    /* If the densest particle found locally is not the global max, make sure we
     * don't seed two black holes. */
    /* If the local index has been set to a foreign index then we don't need to
     * seed a black hole locally. */
    if (max_part_density_index[offset] ==
        fof_mass_recv[i].max_part_density_index) {
      max_part_density_index[offset] = fof_halo_has_black_hole;
    }
    /* The densest particle is on the same node as the global root so we don't
       need to seed a black hole on the other node. */
    else {
      fof_mass_recv[i].max_part_density_index = fof_halo_has_black_hole;
    }
  }

  /* Send the result back */
  MPI_Alltoallv(fof_mass_recv, recvcount, recvoffset, fof_final_mass_type,
                fof_mass_send, sendcount, sendoffset, fof_final_mass_type,
                MPI_COMM_WORLD);

  int extra_seed_count = 0;
  size_t density_index_size = num_groups_local;

  /* Add the index of the densest particle to the local list if the global root
   * is not on this node. */
  for (size_t i = 0; i < nsend; i++) {

    /* Only add the index if:
     * 1) there is not already a black hole in the group
     * AND
     * 2) there is gas in the group. */
    if (fof_mass_send[i].max_part_density_index >= 0) {

      /* Re-allocate the list if it's needed. */
      if (num_groups_local + extra_seed_count >= density_index_size) {
        const size_t new_size = 2 * density_index_size;

        max_part_density_index = (long long *)realloc(
            max_part_density_index, new_size * sizeof(long long));

        message(
            "Re-allocating max_part_density_index from %zu to %zu elements.",
            density_index_size, new_size);

        density_index_size = new_size;

        props->max_part_density_index = max_part_density_index;
      }

      /* Add particle index onto the end of the array. */
      max_part_density_index[num_groups_local + extra_seed_count] =
          fof_mass_send[i].max_part_density_index;
      extra_seed_count++;
    }
  }

  props->extra_bh_seed_count = extra_seed_count;

  free(sendcount);
  free(recvcount);
  free(sendoffset);
  free(recvoffset);
  free(fof_mass_send);
  free(fof_mass_recv);
#else

  /* Allocate and initialise the densest particle array. */
  if (swift_memalign("max_part_density_index",
                     (void **)&props->max_part_density_index, 32,
                     num_groups_local * sizeof(long long)) != 0)
    error(
        "Failed to allocate list of max group density indices for FOF search.");

  if (swift_memalign("max_part_density", (void **)&props->max_part_density, 32,
                     num_groups_local * sizeof(float)) != 0)
    error("Failed to allocate list of max group densities for FOF search.");

  /* Direct pointers to the arrays */
  long long *max_part_density_index = props->max_part_density_index;
  float *max_part_density = props->max_part_density;

  /* No densest particle found so far */
  bzero(max_part_density, num_groups_local * sizeof(float));

  /* Start by assuming that the haloes have no gas */
  for (size_t i = 0; i < num_groups_local; i++) {
    max_part_density_index[i] = fof_halo_has_no_gas;
  }

  /* Increment the group mass for groups above min_group_size. */
  threadpool_map(&s->e->threadpool, fof_calc_group_mass_mapper, gparts,
                 nr_gparts, sizeof(struct gpart), 0, (struct space *)s);

  /* Loop over particles and find the densest particle in each group. */
  /* JSW TODO: Parallelise with threadpool*/
  for (size_t i = 0; i < nr_gparts; i++) {

    const size_t index = gparts[i].group_id - group_id_offset;

    /* Only check groups above the minimum mass threshold. */
    if (gparts[i].group_id != group_id_default) {

      if (group_mass[index] > seed_halo_mass) {

        /* Find the densest gas particle.
         * Account for groups that already have a black hole and groups that
         * contain no gas. */
        if (gparts[i].type == swift_type_gas &&
            max_part_density_index[index] != fof_halo_has_black_hole) {

          const size_t gas_index = -gparts[i].id_or_neg_offset;
          const float rho_com = hydro_get_comoving_density(&parts[gas_index]);

          /* Update index if a denser gas particle is found. */
          if (rho_com > max_part_density[index]) {
            max_part_density_index[index] = gas_index;
            max_part_density[index] = rho_com;
          }
        }
        /* If there is already a black hole in the group we don't need to create
           a new one. */
        else if (gparts[i].type == swift_type_black_hole) {
          max_part_density_index[index] = fof_halo_has_black_hole;
        }

      } else {
        max_part_density_index[index] = fof_halo_has_too_low_mass;
      }
    }
  }

  props->extra_bh_seed_count = 0;
#endif
}

/**
 * @brief Mapper function to perform FOF search.
 *
 * @param map_data An array of #cell pair indices.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_find_foreign_links_mapper(void *map_data, int num_elements,
                                   void *extra_data) {

#ifdef WITH_MPI

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  const int periodic = s->periodic;
  const size_t nr_gparts = s->nr_gparts;
  const struct gpart *const gparts = s->gparts;
  const struct engine *e = s->e;
  struct fof_props *props = e->fof_properties;
  struct cell_pair_indices *cell_pairs = (struct cell_pair_indices *)map_data;

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = props->l_x2;

  /* Store links in an array local to this thread. */
  int local_link_count = 0;
  int local_group_links_size = props->group_links_size / e->nr_threads;

  /* Init the local group links buffer. */
  struct fof_mpi *local_group_links = (struct fof_mpi *)swift_calloc(
      "fof_local_group_links", sizeof(struct fof_mpi), local_group_links_size);
  if (local_group_links == NULL)
    error("Failed to allocate temporary group links buffer.");

  /* Loop over all pairs of local and foreign cells, recurse then perform a
   * FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the local and foreign cells to recurse on. */
    struct cell *restrict local_cell = cell_pairs[ind].local;
    struct cell *restrict foreign_cell = cell_pairs[ind].foreign;

    rec_fof_search_pair_foreign(props, dim, search_r2, periodic, gparts,
                                nr_gparts, local_cell, foreign_cell,
                                &local_link_count, &local_group_links,
                                &local_group_links_size);
  }

  /* Add links found by this thread to the global link list. */
  /* Lock to prevent race conditions while adding to the global link list.*/
  if (lock_lock(&s->lock) == 0) {

    /* Get pointers to global arrays. */
    int *group_links_size = &props->group_links_size;
    int *group_link_count = &props->group_link_count;
    struct fof_mpi **group_links = &props->group_links;

    /* If the global group_links array is not big enough re-allocate it. */
    if (*group_link_count + local_link_count > *group_links_size) {

      const int old_size = *group_links_size;
      const int new_size =
          max(*group_link_count + local_link_count, 2 * old_size);

      (*group_links) = (struct fof_mpi *)realloc(
          *group_links, new_size * sizeof(struct fof_mpi));

      *group_links_size = new_size;

      message("Re-allocating global group links from %d to %d elements.",
              old_size, new_size);
    }

    /* Copy the local links to the global list. */
    for (int i = 0; i < local_link_count; i++) {

      int found = 0;

      /* Check that the links have not already been added to the list by another
       * thread. */
      for (int l = 0; l < *group_link_count; l++) {
        if ((*group_links)[l].group_i == local_group_links[i].group_i &&
            (*group_links)[l].group_j == local_group_links[i].group_j) {
          found = 1;
          break;
        }
      }

      if (!found) {

        (*group_links)[*group_link_count].group_i =
            local_group_links[i].group_i;
        (*group_links)[*group_link_count].group_i_size =
            local_group_links[i].group_i_size;

        (*group_links)[*group_link_count].group_j =
            local_group_links[i].group_j;
        (*group_links)[*group_link_count].group_j_size =
            local_group_links[i].group_j_size;

        (*group_link_count) = (*group_link_count) + 1;
      }
    }
  }

  /* Release lock. */
  if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space");

  swift_free("fof_local_group_links", local_group_links);
#endif
}

void fof_seed_black_holes(const struct fof_props *props,
                          const struct black_holes_props *bh_props,
                          const struct phys_const *constants,
                          const struct cosmology *cosmo, struct space *s,
                          int num_groups, struct group_length *group_sizes) {

  const long long *max_part_density_index = props->max_part_density_index;

  /* Count the number of black holes to seed */
  int num_seed_black_holes = 0;
  for (int i = 0; i < num_groups + props->extra_bh_seed_count; i++) {
    if (max_part_density_index[i] >= 0) ++num_seed_black_holes;
  }

#ifdef WITH_MPI
  int total_num_seed_black_holes = 0;
  /* Sum the total number of black holes over each MPI rank. */
  MPI_Reduce(&num_seed_black_holes, &total_num_seed_black_holes, 1, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
#else
  int total_num_seed_black_holes = num_seed_black_holes;
#endif

  if (engine_rank == 0)
    message("Seeding %d black hole(s)", total_num_seed_black_holes);

  /* Anything to do this time on this rank? */
  if (num_seed_black_holes == 0) return;

  /* Do we need to reallocate the black hole array for the new particles? */
  if (s->nr_bparts + num_seed_black_holes > s->size_bparts) {
    const size_t nr_bparts_new = s->nr_bparts + num_seed_black_holes;

    s->size_bparts = engine_parts_size_grow * nr_bparts_new;

    struct bpart *bparts_new = NULL;
    if (swift_memalign("bparts", (void **)&bparts_new, bpart_align,
                       sizeof(struct bpart) * s->size_bparts) != 0)
      error("Failed to allocate new bpart data.");
    memcpy(bparts_new, s->bparts, sizeof(struct bpart) * s->nr_bparts);
    swift_free("bparts", s->bparts);

    s->bparts = bparts_new;
  }

  int k = s->nr_bparts;

  /* Loop over the local groups */
  for (int i = 0; i < num_groups + props->extra_bh_seed_count; i++) {

    const long long part_index = max_part_density_index[i];

    /* Should we seed? */
    if (part_index >= 0) {

      /* Handle on the particle to convert */
      struct part *p = &s->parts[part_index];
      struct gpart *gp = p->gpart;

      /* Let's destroy the gas particle */
      p->time_bin = time_bin_inhibited;
      p->gpart = NULL;

      /* Mark the gpart as black hole */
      gp->type = swift_type_black_hole;

      /* Basic properties of the black hole */
      struct bpart *bp = &s->bparts[k];
      bzero(bp, sizeof(struct bpart));
      bp->time_bin = gp->time_bin;

      /* Re-link things */
      bp->gpart = gp;
      gp->id_or_neg_offset = -(bp - s->bparts);

      /* Synchronize masses, positions and velocities */
      bp->mass = gp->mass;
      bp->x[0] = gp->x[0];
      bp->x[1] = gp->x[1];
      bp->x[2] = gp->x[2];
      bp->v[0] = gp->v_full[0];
      bp->v[1] = gp->v_full[1];
      bp->v[2] = gp->v_full[2];

      /* Set a smoothing length */
      bp->h = p->h;

#ifdef SWIFT_DEBUG_CHECKS
      bp->ti_kick = p->ti_kick;
      bp->ti_drift = p->ti_drift;
#endif

      /* Copy over all the gas properties that we want */
      black_holes_create_from_gas(bp, bh_props, constants, cosmo, p);

      /* Move to the next BH slot */
      k++;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if ((int)(s->nr_bparts) + num_seed_black_holes != k) {
    error("Seeded the wrong number of black holes!");
  }
#endif

  /* Update the count of black holes. */
  s->nr_bparts = k;
}

/* Dump FOF group data. */
void fof_dump_group_data(const struct fof_props *props,
                         const char *out_file_name, struct space *s,
                         int num_groups, struct group_length *group_sizes) {

  FILE *file = fopen(out_file_name, "w");

  struct gpart *gparts = s->gparts;
  struct part *parts = s->parts;
  size_t *group_size = props->group_size;
  double *group_mass = props->group_mass;
  const long long *max_part_density_index = props->max_part_density_index;
  const float *max_part_density = props->max_part_density;

  fprintf(file, "# %8s %12s %12s %12s %18s %18s %12s\n", "Group ID",
          "Group Size", "Group Mass", "Max Density", "Max Density Index",
          "Particle ID", "Particle Density");
  fprintf(file,
          "#-------------------------------------------------------------------"
          "-------------\n");

  int bh_seed_count = 0;

  for (int i = 0; i < num_groups; i++) {

    const size_t group_offset = group_sizes[i].index;
    const long long part_id = max_part_density_index[i] >= 0
                                  ? parts[max_part_density_index[i]].id
                                  : -1;
#ifdef WITH_MPI
    fprintf(file, "  %8zu %12zu %12e %12e %18lld %18lld\n",
            gparts[group_offset - node_offset].group_id,
            group_size[group_offset - node_offset], group_mass[i],
            max_part_density[i], max_part_density_index[i], part_id);
#else
    fprintf(file, "  %8zu %12zu %12e %12e %18lld %18lld\n",
            gparts[group_offset].group_id, group_size[group_offset],
            group_mass[i], max_part_density[i], max_part_density_index[i],
            part_id);
#endif

    if (max_part_density_index[i] >= 0) bh_seed_count++;
  }

  /* Dump the extra black hole seeds. */
  for (int i = num_groups; i < num_groups + props->extra_bh_seed_count; i++) {
    const long long part_id = max_part_density_index[i] >= 0
                                  ? parts[max_part_density_index[i]].id
                                  : -1;
    fprintf(file, "  %8zu %12zu %12e %12e %18lld %18lld\n", 0UL, 0UL, 0., 0.,
            0LL, part_id);

    if (max_part_density_index[i] >= 0) bh_seed_count++;
  }

  int total_bh_seed_count = 0;

#ifdef WITH_MPI
  /* Sum the total number of black holes over each MPI rank. */
  MPI_Reduce(&bh_seed_count, &total_bh_seed_count, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
#else
  total_bh_seed_count = bh_seed_count;
#endif

  if (engine_rank == 0)
    message("Seeding %d black hole(s).", total_bh_seed_count);

  fclose(file);
}

/**
 * @brief Search foreign cells for links and communicate any found to the
 * appropriate node.
 *
 * @param props the properties of the FOF scheme.
 * @param s Pointer to a #space.
 */
void fof_search_foreign_cells(struct fof_props *props, const struct space *s) {

#ifdef WITH_MPI

  struct engine *e = s->e;
  int verbose = e->verbose;
  size_t *group_index = props->group_index;
  size_t *group_size = props->group_size;
  const size_t nr_gparts = s->nr_gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = props->l_x2;

  ticks tic = getticks();

  /* Make group IDs globally unique. */
  for (size_t i = 0; i < nr_gparts; i++) group_index[i] += node_offset;

  struct cell_pair_indices *cell_pairs = NULL;
  int group_link_count = 0;
  int cell_pair_count = 0;

  props->group_links_size = fof_props_default_group_link_size;
  props->group_link_count = 0;

  int num_cells_out = 0;
  int num_cells_in = 0;

  /* Find the maximum no. of cell pairs. */
  for (int i = 0; i < e->nr_proxies; i++) {

    for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_out_type[j] & proxy_cell_type_gravity)
        num_cells_out++;
    }

    for (int j = 0; j < e->proxies[i].nr_cells_in; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_in_type[j] & proxy_cell_type_gravity)
        num_cells_in++;
    }
  }

  const int cell_pair_size = num_cells_in * num_cells_out;

  if (swift_memalign("fof_group_links", (void **)&props->group_links,
                     SWIFT_STRUCT_ALIGNMENT,
                     props->group_links_size * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for FOF links over an MPI domain");

  if (swift_memalign("fof_cell_pairs", (void **)&cell_pairs,
                     SWIFT_STRUCT_ALIGNMENT,
                     cell_pair_size * sizeof(struct cell_pair_indices)) != 0)
    error("Error while allocating memory for FOF cell pair indices");

  /* Loop over cells_in and cells_out for each proxy and find which cells are in
   * range of each other to perform the FOF search. Store local cells that are
   * touching foreign cells in a list. */
  for (int i = 0; i < e->nr_proxies; i++) {

    /* Only find links across an MPI domain on one rank. */
    if (engine_rank == min(engine_rank, e->proxies[i].nodeID)) {

      for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

        /* Skip non-gravity cells. */
        if (!(e->proxies[i].cells_out_type[j] & proxy_cell_type_gravity))
          continue;

        struct cell *restrict local_cell = e->proxies[i].cells_out[j];

        /* Skip empty cells. */
        if (local_cell->grav.count == 0) continue;

        for (int k = 0; k < e->proxies[i].nr_cells_in; k++) {

          /* Skip non-gravity cells. */
          if (!(e->proxies[i].cells_in_type[k] & proxy_cell_type_gravity))
            continue;

          struct cell *restrict foreign_cell = e->proxies[i].cells_in[k];

          /* Skip empty cells. */
          if (foreign_cell->grav.count == 0) continue;

          /* Check if local cell has already been added to the local list of
           * cells. */
          const double r2 = cell_min_dist(local_cell, foreign_cell, dim);
          if (r2 < search_r2) {
            cell_pairs[cell_pair_count].local = local_cell;
            cell_pairs[cell_pair_count++].foreign = foreign_cell;
          }
        }
      }
    }
  }

  /* Set the root of outgoing particles. */
  for (int i = 0; i < e->nr_proxies; i++) {

    for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

      struct cell *restrict local_cell = e->proxies[i].cells_out[j];
      struct gpart *gparts = local_cell->grav.parts;

      /* Make a list of particle offsets into the global gparts array. */
      size_t *const offset = group_index + (ptrdiff_t)(gparts - s->gparts);

      /* Set each particle's root and group properties found in the local FOF.*/
      for (int k = 0; k < local_cell->grav.count; k++) {
        const size_t root =
            fof_find_global(offset[k] - node_offset, group_index, nr_gparts);
        gparts[k].group_id = root;
        gparts[k].group_size = group_size[root - node_offset];
      }
    }
  }

  if (verbose)
    message(
        "Finding local/foreign cell pairs and initialising particle roots "
        "took: "
        "%.3f %s.",
        clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  struct scheduler *sched = &e->sched;
  struct task *tasks = sched->tasks;

  /* Activate the send and receive tasks for the gparts. */
  for (int i = 0; i < sched->nr_tasks; i++) {

    struct task *t = &tasks[i];

    if ((t->type == task_type_send && t->subtype == task_subtype_gpart) ||
        (t->type == task_type_recv && t->subtype == task_subtype_gpart)) {
      scheduler_activate(sched, t);
    } else
      t->skip = 1;
  }

  if (verbose)
    message("MPI send/recv task activation took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  ticks local_fof_tic = getticks();

  MPI_Barrier(MPI_COMM_WORLD);

  if (verbose)
    message("Local FOF imbalance: %.3f %s.",
            clocks_from_ticks(getticks() - local_fof_tic), clocks_getunit());

  tic = getticks();

  /* Perform send and receive tasks. */
  engine_launch(e);

  if (verbose)
    message("MPI send/recv comms took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Perform search of group links between local and foreign cells with the
   * threadpool. */
  threadpool_map(&s->e->threadpool, fof_find_foreign_links_mapper, cell_pairs,
                 cell_pair_count, sizeof(struct cell_pair_indices), 1,
                 (struct space *)s);

  group_link_count = props->group_link_count;

  /* Clean up memory. */
  swift_free("fof_cell_pairs", cell_pairs);

  if (verbose)
    message("Searching for foreign links took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  struct fof_mpi *global_group_links = NULL;
  int *displ = NULL, *group_link_counts = NULL;
  int global_group_link_count = 0;

  ticks comms_tic = getticks();

  MPI_Barrier(MPI_COMM_WORLD);

  if (verbose)
    message("Imbalance took: %.3f %s.",
            clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

  comms_tic = getticks();

  /* Sum the total number of links across MPI domains over each MPI rank. */
  MPI_Allreduce(&group_link_count, &global_group_link_count, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);

  /* Unique set of links is half of all group links as each link is found twice
   * by opposing MPI ranks. */
  if (swift_memalign("fof_global_group_links", (void **)&global_group_links,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_link_count * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for the global list of group links");

  if (posix_memalign((void **)&group_link_counts, SWIFT_STRUCT_ALIGNMENT,
                     e->nr_nodes * sizeof(int)) != 0)
    error(
        "Error while allocating memory for the number of group links on each "
        "MPI rank");

  if (posix_memalign((void **)&displ, SWIFT_STRUCT_ALIGNMENT,
                     e->nr_nodes * sizeof(int)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  /* Gather the total number of links on each rank. */
  MPI_Allgather(&group_link_count, 1, MPI_INT, group_link_counts, 1, MPI_INT,
                MPI_COMM_WORLD);

  /* Set the displacements into the global link list using the link counts from
   * each rank */
  displ[0] = 0;
  for (int i = 1; i < e->nr_nodes; i++)
    displ[i] = displ[i - 1] + group_link_counts[i - 1];

  /* Gather the global link list on all ranks. */
  MPI_Allgatherv(props->group_links, group_link_count, fof_mpi_type,
                 global_group_links, group_link_counts, displ, fof_mpi_type,
                 MPI_COMM_WORLD);

  /* Clean up memory. */
  free(displ);
  swift_free("fof_group_links", props->group_links);
  props->group_links = NULL;

  if (verbose) {
    message("Communication took: %.3f %s.",
            clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

    message("Global comms took: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }

  tic = getticks();

  /* Transform the group IDs to a local list going from 0-group_count so a
   * union-find can be performed. */
  size_t *global_group_index = NULL, *global_group_id = NULL,
         *global_group_size = NULL;
  const int global_group_list_size = 2 * global_group_link_count;
  int group_count = 0;

  if (swift_memalign("fof_global_group_index", (void **)&global_group_index,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (swift_memalign("fof_global_group_id", (void **)&global_group_id,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (swift_memalign("fof_global_group_size", (void **)&global_group_size,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  bzero(global_group_size, global_group_list_size * sizeof(size_t));

  /* Create hash table. */
  hashmap_t map;
  hashmap_init(&map);

  /* Store each group ID and its properties. */
  for (int i = 0; i < global_group_link_count; i++) {

    size_t group_i = global_group_links[i].group_i;
    size_t group_j = global_group_links[i].group_j;

    global_group_size[group_count] += global_group_links[i].group_i_size;
    global_group_id[group_count] = group_i;
    hashmap_add_group(group_i, group_count++, &map);

    global_group_size[group_count] += global_group_links[i].group_j_size;
    global_group_id[group_count] = group_j;
    hashmap_add_group(group_j, group_count++, &map);
  }

  if (verbose)
    message("Global list compression took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Create a global_group_index list of groups across MPI domains so that you
   * can perform a union-find locally on each node. */
  /* The value of which is an offset into global_group_id, which is the actual
   * root. */
  for (int i = 0; i < group_count; i++) global_group_index[i] = i;

  /* Store the original group size before incrementing in the Union-Find. */
  size_t *orig_global_group_size = NULL;

  if (swift_memalign("fof_orig_global_group_size",
                     (void **)&orig_global_group_size, SWIFT_STRUCT_ALIGNMENT,
                     group_count * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  for (int i = 0; i < group_count; i++)
    orig_global_group_size[i] = global_group_size[i];

  /* Perform a union-find on the group links. */
  for (int i = 0; i < global_group_link_count; i++) {

    /* Use the hash table to find the group offsets in the index array. */
    size_t find_i =
        hashmap_find_group_offset(global_group_links[i].group_i, &map);
    size_t find_j =
        hashmap_find_group_offset(global_group_links[i].group_j, &map);

    /* Use the offset to find the group's root. */
    size_t root_i = fof_find(find_i, global_group_index);
    size_t root_j = fof_find(find_j, global_group_index);

    size_t group_i = global_group_id[root_i];
    size_t group_j = global_group_id[root_j];

    if (group_i == group_j) continue;

    /* Update roots accordingly. */
    size_t size_i = global_group_size[root_i];
    size_t size_j = global_group_size[root_j];
#ifdef UNION_BY_SIZE_OVER_MPI
    if (size_i < size_j) {
      global_group_index[root_i] = root_j;
      global_group_size[root_j] += size_i;
    } else {
      global_group_index[root_j] = root_i;
      global_group_size[root_i] += size_j;
    }
#else
    if (group_j < group_i) {
      global_group_index[root_i] = root_j;
      global_group_size[root_j] += size_i;
    } else {
      global_group_index[root_j] = root_i;
      global_group_size[root_i] += size_j;
    }
#endif
  }

  hashmap_free(&map);

  if (verbose)
    message("global_group_index construction took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Update each group locally with new root information. */
  for (int i = 0; i < group_count; i++) {

    size_t group_id = global_group_id[i];
    size_t offset = fof_find(global_group_index[i], global_group_index);
    size_t new_root = global_group_id[offset];

    /* If the group is local update its root and size. */
    if (is_local(group_id, nr_gparts) && new_root != group_id) {
      group_index[group_id - node_offset] = new_root;

      group_size[group_id - node_offset] -= orig_global_group_size[i];
    }

    /* If the group linked to a local root update its size. */
    if (is_local(new_root, nr_gparts) && new_root != group_id) {

      /* Use group sizes before Union-Find */
      group_size[new_root - node_offset] += orig_global_group_size[i];
    }
  }

  if (verbose)
    message("Updating groups locally took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Clean up memory. */
  swift_free("fof_global_group_links", global_group_links);
  swift_free("fof_global_group_index", global_group_index);
  swift_free("fof_global_group_size", global_group_size);
  swift_free("fof_global_group_id", global_group_id);
  swift_free("fof_orig_global_group_size", orig_global_group_size);

#endif /* WITH_MPI */
}

/**
 * @brief Perform a FOF search on gravity particles using the cells and applying
 * the Union-Find algorithm.
 *
 * @param props The properties of the FOF scheme.
 * @param bh_props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param s The #space containing the particles.
 * @param dump_results Do we want to write the group catalogue to a file?
 * @param seed_black_holes Do we want to seed black holes in haloes?
 */
void fof_search_tree(struct fof_props *props,
                     const struct black_holes_props *bh_props,
                     const struct phys_const *constants,
                     const struct cosmology *cosmo, struct space *s,
                     const int dump_results, const int seed_black_holes) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t min_group_size = props->min_group_size;
  const size_t group_id_offset = props->group_id_offset;
#ifdef WITH_MPI
  const int nr_nodes = s->e->nr_nodes;
#endif
  struct gpart *gparts = s->gparts;
  size_t *group_index, *group_size;
  int num_groups = 0, num_parts_in_groups = 0, max_group_size = 0;
  int verbose = s->e->verbose;
  ticks tic_total = getticks();

  char output_file_name[PARSER_MAX_LINE_SIZE];
  snprintf(output_file_name, PARSER_MAX_LINE_SIZE, "%s", props->base_name);

  if (verbose)
    message("Searching %zu gravity particles for links with l_x: %lf",
            nr_gparts, sqrt(props->l_x2));

  if (engine_rank == 0 && verbose)
    message("Size of hash table element: %ld", sizeof(hashmap_element_t));

  const size_t group_id_default = props->group_id_default;

#ifdef WITH_MPI

  /* Reset global variable */
  node_offset = 0;

  /* Determine number of gparts on lower numbered MPI ranks */
  long long nr_gparts_cumulative;
  long long nr_gparts_local = s->nr_gparts;
  MPI_Scan(&nr_gparts_local, &nr_gparts_cumulative, 1, MPI_LONG_LONG, MPI_SUM,
           MPI_COMM_WORLD);
  node_offset = nr_gparts_cumulative - nr_gparts_local;

  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
           "_mpi_rank_%d.dat", engine_rank);
#else
  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
           ".dat");
#endif

  /* Local copy of the arrays */
  group_index = props->group_index;
  group_size = props->group_size;

  ticks tic_calc_group_size = getticks();

  threadpool_map(&s->e->threadpool, fof_calc_group_size_mapper, gparts,
                 nr_gparts, sizeof(struct gpart), 0, s);

  if (verbose)
    message("FOF calc group size took (scaling): %.3f %s.",
            clocks_from_ticks(getticks() - tic_calc_group_size),
            clocks_getunit());

#ifdef WITH_MPI
  if (nr_nodes > 1) {

    ticks tic_mpi = getticks();

    /* Search for group links across MPI domains. */
    fof_search_foreign_cells(props, s);

    if (verbose)
      message("fof_search_foreign_cells() took: %.3f %s.",
              clocks_from_ticks(getticks() - tic_mpi), clocks_getunit());
  }
#endif

  size_t num_groups_local = 0, num_parts_in_groups_local = 0,
         max_group_size_local = 0;

  for (size_t i = 0; i < nr_gparts; i++) {

#ifdef WITH_MPI
    /* Find the total number of groups. */
    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size)
      num_groups_local++;
#else
    /* Find the total number of groups. */
    if (group_index[i] == i && group_size[i] >= min_group_size)
      num_groups_local++;
#endif

    /* Find the total number of particles in groups. */
    if (group_size[i] >= min_group_size)
      num_parts_in_groups_local += group_size[i];

    /* Find the largest group. */
    if (group_size[i] > max_group_size_local)
      max_group_size_local = group_size[i];
  }

  /* Sort the groups in descending order based upon size and re-label their IDs
   * 0-num_groups. */
  struct group_length *high_group_sizes = NULL;
  int group_count = 0;

  if (swift_memalign("fof_high_group_sizes", (void **)&high_group_sizes, 32,
                     num_groups_local * sizeof(struct group_length)) != 0)
    error("Failed to allocate list of large groups.");

  /* Store the group_sizes and their offset. */
  for (size_t i = 0; i < nr_gparts; i++) {

#ifdef WITH_MPI
    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size) {
      high_group_sizes[group_count].index = node_offset + i;
      high_group_sizes[group_count++].size = group_size[i];
    }
#else
    if (group_index[i] == i && group_size[i] >= min_group_size) {
      high_group_sizes[group_count].index = i;
      high_group_sizes[group_count++].size = group_size[i];
    }
#endif
  }

  ticks tic = getticks();

  /* Find global properties. */
#ifdef WITH_MPI
  MPI_Allreduce(&num_groups_local, &num_groups, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Reduce(&num_parts_in_groups_local, &num_parts_in_groups, 1, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_group_size_local, &max_group_size, 1, MPI_INT, MPI_MAX, 0,
             MPI_COMM_WORLD);
#else
  num_groups = num_groups_local;
  num_parts_in_groups = num_parts_in_groups_local;
  max_group_size = max_group_size_local;
#endif /* WITH_MPI */
  props->num_groups = num_groups;

  /* Find number of groups on lower numbered MPI ranks */
#ifdef WITH_MPI
  long long nglocal = num_groups_local;
  long long ngsum;
  MPI_Scan(&nglocal, &ngsum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  const size_t num_groups_prev = (size_t)(ngsum - nglocal);
#endif /* WITH_MPI */

  /* Sort local groups into descending order of size */
  qsort(high_group_sizes, num_groups_local, sizeof(struct group_length),
        cmp_func_group_size);

  /* Set default group ID for all particles */
  for (size_t i = 0; i < nr_gparts; i++) gparts[i].group_id = group_id_default;

  /*
    Assign final group IDs to local root particles where the global root is on
    this node and the group is large enough. Within a node IDs are assigned in
    descending order of particle number.
  */
  for (size_t i = 0; i < num_groups_local; i++) {
#ifdef WITH_MPI
    gparts[high_group_sizes[i].index - node_offset].group_id =
        group_id_offset + i + num_groups_prev;
#else
    gparts[high_group_sizes[i].index].group_id = group_id_offset + i;
#endif
  }

#ifdef WITH_MPI
  /*
     Now, for each local root where the global root is on some other node
     AND the total size of the group is >= min_group_size we need to retrieve
     the gparts.group_id we just assigned to the global root.

     Will do that by sending the group_index of these lcoal roots to the node
     where their global root is stored and receiving back the new group_id
     associated with that particle.
  */

  /*
     Identify local roots with global root on another node and large enough
     group_size. Store index of the local and global roots in these cases.

     NOTE: if group_size only contains the total FoF mass for global roots,
     then we have to communicate ALL fragments where the global root is not
     on this node. Hence the commented out extra conditions below.
  */
  size_t nsend = 0;
  for (size_t i = 0; i < nr_gparts; i += 1) {
    if ((!is_local(group_index[i],
                   nr_gparts))) { /* && (group_size[i] >= min_group_size)) { */
      nsend += 1;
    }
  }
  struct fof_final_index *fof_index_send =
      swift_malloc("fof_index_send", sizeof(struct fof_final_index) * nsend);
  nsend = 0;
  for (size_t i = 0; i < nr_gparts; i += 1) {
    if ((!is_local(group_index[i],
                   nr_gparts))) { /* && (group_size[i] >= min_group_size)) { */
      fof_index_send[nsend].local_root = node_offset + i;
      fof_index_send[nsend].global_root = group_index[i];
      nsend += 1;
    }
  }

  /* Sort by global root - this puts the groups in order of which node they're
   * stored on */
  qsort(fof_index_send, nsend, sizeof(struct fof_final_index),
        compare_fof_final_index_global_root);

  /* Determine range of global indexes (i.e. particles) on each node */
  size_t *num_on_node = malloc(nr_nodes * sizeof(size_t));
  MPI_Allgather(&nr_gparts, sizeof(size_t), MPI_BYTE, num_on_node,
                sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);
  size_t *first_on_node = malloc(nr_nodes * sizeof(size_t));
  first_on_node[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1)
    first_on_node[i] = first_on_node[i - 1] + num_on_node[i - 1];

  /* Determine how many entries go to each node */
  int *sendcount = malloc(nr_nodes * sizeof(int));
  for (int i = 0; i < nr_nodes; i += 1) sendcount[i] = 0;
  int dest = 0;
  for (size_t i = 0; i < nsend; i += 1) {
    while ((fof_index_send[i].global_root >=
            first_on_node[dest] + num_on_node[dest]) ||
           (num_on_node[dest] == 0))
      dest += 1;
    if (dest >= nr_nodes) error("Node index out of range!");
    sendcount[dest] += 1;
  }

  int *recvcount = NULL, *sendoffset = NULL, *recvoffset = NULL;
  size_t nrecv = 0;

  fof_compute_send_recv_offsets(nr_nodes, sendcount, &recvcount, &sendoffset,
                                &recvoffset, &nrecv);

  struct fof_final_index *fof_index_recv =
      swift_malloc("fof_index_recv", nrecv * sizeof(struct fof_final_index));

  /* Exchange group indexes */
  MPI_Alltoallv(fof_index_send, sendcount, sendoffset, fof_final_index_type,
                fof_index_recv, recvcount, recvoffset, fof_final_index_type,
                MPI_COMM_WORLD);

  /* For each received global root, look up the group ID we assigned and store
   * it in the struct */
  for (size_t i = 0; i < nrecv; i += 1) {
    if ((fof_index_recv[i].global_root < node_offset) ||
        (fof_index_recv[i].global_root >= node_offset + nr_gparts)) {
      error("Received global root index out of range!");
    }
    fof_index_recv[i].global_root =
        gparts[fof_index_recv[i].global_root - node_offset].group_id;
  }

  /* Send the result back */
  MPI_Alltoallv(fof_index_recv, recvcount, recvoffset, fof_final_index_type,
                fof_index_send, sendcount, sendoffset, fof_final_index_type,
                MPI_COMM_WORLD);

  /* Update local gparts.group_id */
  for (size_t i = 0; i < nsend; i += 1) {
    if ((fof_index_send[i].local_root < node_offset) ||
        (fof_index_send[i].local_root >= node_offset + nr_gparts)) {
      error("Sent local root index out of range!");
    }
    gparts[fof_index_send[i].local_root - node_offset].group_id =
        fof_index_send[i].global_root;
  }

  free(sendcount);
  free(recvcount);
  free(sendoffset);
  free(recvoffset);
  swift_free("fof_index_send", fof_index_send);
  swift_free("fof_index_recv", fof_index_recv);

#endif /* WITH_MPI */

  /* Assign every particle the group_id of its local root. */
  for (size_t i = 0; i < nr_gparts; i++) {
    const size_t root = fof_find_local(i, nr_gparts, group_index);
    gparts[i].group_id = gparts[root].group_id;
  }

  if (verbose)
    message("Group sorting took: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Allocate and initialise a group mass array. */
  if (swift_memalign("group_mass", (void **)&props->group_mass, 32,
                     num_groups_local * sizeof(double)) != 0)
    error("Failed to allocate list of group masses for FOF search.");

  bzero(props->group_mass, num_groups_local * sizeof(double));

  ticks tic_seeding = getticks();

  double *group_mass = props->group_mass;
#ifdef WITH_MPI
  fof_calc_group_mass(props, s, num_groups_local, num_groups_prev, num_on_node,
                      first_on_node, group_mass);
  free(num_on_node);
  free(first_on_node);
#else
  fof_calc_group_mass(props, s, num_groups_local, 0, NULL, NULL, group_mass);
#endif

  if (verbose)
    message("Black hole seeding took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_seeding), clocks_getunit());

  /* Dump group data. */
  if (dump_results) {
    fof_dump_group_data(props, output_file_name, s, num_groups_local,
                        high_group_sizes);
  }

  /* Seed black holes */
  if (seed_black_holes) {
    fof_seed_black_holes(props, bh_props, constants, cosmo, s, num_groups_local,
                         high_group_sizes);
  }

  /* Free the left-overs */
  swift_free("fof_high_group_sizes", high_group_sizes);
  swift_free("fof_group_index", props->group_index);
  swift_free("fof_group_size", props->group_size);
  swift_free("fof_group_mass", props->group_mass);
  swift_free("fof_max_part_density_index", props->max_part_density_index);
  swift_free("fof_max_part_density", props->max_part_density);
  props->group_index = NULL;
  props->group_size = NULL;
  props->group_mass = NULL;
  props->max_part_density_index = NULL;
  props->max_part_density = NULL;

  if (engine_rank == 0) {
    message(
        "No. of groups: %d. No. of particles in groups: %d. No. of particles "
        "not in groups: %lld.",
        num_groups, num_parts_in_groups,
        s->e->total_nr_gparts - num_parts_in_groups);

    message("Largest group by size: %d", max_group_size);
  }
  if (verbose)
    message("FOF search took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_total), clocks_getunit());

#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void fof_struct_dump(const struct fof_props *props, FILE *stream) {

  struct fof_props temp = *props;
  temp.num_groups = 0;
  temp.group_link_count = 0;
  temp.group_links_size = 0;
  temp.group_index = NULL;
  temp.group_size = NULL;
  temp.group_mass = NULL;
  temp.max_part_density_index = NULL;
  temp.max_part_density = NULL;
  temp.group_links = NULL;

  restart_write_blocks((void *)&temp, sizeof(struct fof_props), 1, stream,
                       "fof_props", "fof_props");
}

void fof_struct_restore(struct fof_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct fof_props), 1, stream, NULL,
                      "fof_props");
}
