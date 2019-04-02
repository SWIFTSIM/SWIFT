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
#include "common_io.h"
#include "engine.h"
#include "proxy.h"
#include "threadpool.h"
#include "c_hashmap/hashmap.h"

#ifdef WITH_MPI
MPI_Datatype fof_mpi_type;
MPI_Datatype group_length_mpi_type;
MPI_Datatype fof_final_index_type;
MPI_Datatype fof_final_mass_type;
#endif
size_t node_offset;

#define UNION_BY_SIZE_OVER_MPI (1)
#define FOF_COMPRESS_PATHS_MIN_LENGTH (2)

/* Initialises parameters for the FOF search. */
void fof_init(struct space *s) {

  struct engine *e = s->e;

  /* Check that we can write outputs by testing if the output
   * directory exists and is searchable and writable. */
  parser_get_param_string(e->parameter_file, "FOF:basename",
                          s->fof_data.base_name);

  const char *dirp = dirname(s->fof_data.base_name);
  if (access(dirp, W_OK | X_OK) != 0) {
    error("Cannot write FOF outputs in directory %s (%s)", dirp,
          strerror(errno));
  }

  /* Read the minimum group size. */
  s->fof_data.min_group_size =
      parser_get_opt_param_int(e->parameter_file, "FOF:min_group_size", 20);

  /* Read the default group ID of particles in groups below the minimum group
   * size. */
  const int default_id = parser_get_opt_param_int(
      e->parameter_file, "FOF:group_id_default", 2147483647);

  /* Make sure default group ID is positive. */
  if (default_id < 0)
    error("The default group ID set: %d, has to be positive.", default_id);

  s->fof_data.group_id_default = default_id;

  /* Read the starting group ID. */
  s->fof_data.group_id_offset =
      parser_get_opt_param_int(e->parameter_file, "FOF:group_id_offset", 1);

  /* Read the linking length scale. */
  const double l_x_scale = parser_get_opt_param_double(
      e->parameter_file, "FOF:linking_length_scale", 0.2);

  /* Calculate the particle linking length based upon the mean inter-particle
   * spacing of the DM particles. */
  const long long total_nr_dmparts =
      e->total_nr_gparts - e->total_nr_parts - e->total_nr_sparts;
  double l_x = l_x_scale * (s->dim[0] / cbrt(total_nr_dmparts));

  l_x = parser_get_opt_param_double(e->parameter_file,
                                    "FOF:absolute_linking_length", l_x);

  s->l_x2 = l_x * l_x;

  /* Read the initial group_links array size. */
  s->fof_data.group_links_size_default = parser_get_opt_param_double(
      e->parameter_file, "FOF:group_links_size_default", 20000);

  const size_t nr_local_gparts = s->nr_gparts;

  /* Allocate and initialise a group index array. */
  if (posix_memalign((void **)&s->fof_data.group_index, 32,
                     nr_local_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle group indices for FOF search.");

  /* Allocate and initialise a group size array. */
  if (posix_memalign((void **)&s->fof_data.group_size, 32,
                     nr_local_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of group size for FOF search.");

  /* Set initial group index to particle offset into array and set default group
   * ID. */
  for (size_t i = 0; i < nr_local_gparts; i++) {
    s->fof_data.group_index[i] = i;
    s->gparts[i].group_id = default_id;
    s->fof_data.group_size[i] = 1;
  }

#ifdef WITH_MPI
  /* Check size of linking length against the top-level cell dimensions. */
  if (l_x > s->width[0])
    error(
        "Linking length greater than the width of a top-level cell. Need to "
        "check more than one layer of top-level cells for links.");

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
  if (MPI_Type_contiguous(sizeof(struct fof_final_index), MPI_BYTE, &fof_final_index_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_final_index_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof_final_index.");
  }
  /* Define type for sending fof_final_mass struct */  
  if (MPI_Type_contiguous(sizeof(struct fof_final_mass), MPI_BYTE, &fof_final_mass_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_final_mass_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof_final_mass.");
  }
#endif

#ifdef UNION_BY_SIZE_OVER_MPI
  message(
      "Performing FOF over MPI using union by size and union by rank locally.");
#else
  message("Performing FOF using union by rank.");
#endif
}

/*
 * Sort elements in descending order.
 *
 * Return value meaning
 * <0 The element pointed by b goes before the element pointed by a
 * 0  The element pointed by b is equivalent to the element pointed by a
 * >0 The element pointed by b goes after the element pointed by a
 *
 */
int cmp_func(const void *a, const void *b) {
  struct group_length *a_group_size = (struct group_length *)a;
  struct group_length *b_group_size = (struct group_length *)b;
  if(b_group_size->size > a_group_size->size)
    return 1;
  else if(b_group_size->size < a_group_size->size)
    return -1;
  else
    return 0;
}

#ifdef WITH_MPI
/*
 * Sort elements of struct fof_final_index in ASCENDING order of global_root.
 *
 * Return value meaning
 * <0 The element pointed by b goes before the element pointed by a
 * 0  The element pointed by b is equivalent to the element pointed by a
 * >0 The element pointed by b goes after the element pointed by a
 *
 */
int compare_fof_final_index_global_root(const void *a, const void *b) {
  struct fof_final_index *fof_final_index_a = (struct fof_final_index *)a;
  struct fof_final_index *fof_final_index_b = (struct fof_final_index *)b;
  if(fof_final_index_b->global_root < fof_final_index_a->global_root)
    return 1;
  else if(fof_final_index_b->global_root > fof_final_index_a->global_root)
    return -1;
  else
    return 0;
}

int compare_fof_final_mass_global_root(const void *a, const void *b) {
  struct fof_final_mass *fof_final_mass_a = (struct fof_final_mass *)a;
  struct fof_final_mass *fof_final_mass_b = (struct fof_final_mass *)b;
  if(fof_final_mass_b->global_root < fof_final_mass_a->global_root)
    return 1;
  else if(fof_final_mass_b->global_root > fof_final_mass_a->global_root)
    return -1;
  else
    return 0;
}
#endif

/* Checks whether the group is on the local node. */
__attribute__((always_inline)) INLINE static int is_local(
    const size_t group_id, const size_t nr_gparts) {
  return (group_id >= node_offset && group_id < node_offset + nr_gparts);
}

/* Finds the global root ID of the group a particle exists in. */
__attribute__((always_inline)) INLINE static size_t fof_find_global(
    const size_t i, size_t *group_index, const size_t nr_gparts) {

  size_t root = node_offset + i;
  if (!is_local(root, nr_gparts)) return root;

  while (root != group_index[root - node_offset]) {
    root = group_index[root - node_offset];
    if (!is_local(root, nr_gparts)) break;
  }

  /* Perform path compression. */
  // int index = i;
  // while(index != root) {
  //  int next = group_index[index];
  //  group_index[index] = root;
  //  index = next;
  //}

  return root;
}

/*
  Finds the local root ID of the group a particle exists in
  when group_index contains globally unique identifiers -
  i.e. we stop *before* we advance to a foreign root.

  Here we assume that the input i is a local index and we
  return the local index of the root.
*/
__attribute__((always_inline)) INLINE static size_t fof_find_local(const size_t i, 
                                                                   const size_t nr_gparts, 
                                                                   size_t *group_index) {
  size_t root = node_offset + i;

  while ((group_index[root - node_offset] != root) && 
         (group_index[root - node_offset] >= node_offset) &&
         (group_index[root - node_offset] < node_offset + nr_gparts)) {
    root = group_index[root - node_offset];
  }

  return root - node_offset;
}

/* Finds the local root ID of the group a particle exists in. */
__attribute__((always_inline)) INLINE static size_t fof_find(
    const size_t i, size_t *group_index) {

  size_t root = i;
  size_t tree_depth = 0;
  while (root != group_index[root]) {
#ifdef PATH_HALVING
    atomic_cas(&group_index[root], group_index[root],
               group_index[group_index[root]]);
#endif
    root = group_index[root];
    tree_depth++;
  }

  /* Only perform path compression on trees with a depth of FOF_COMPRESS_PATHS_MIN_LENGTH or higher. */
  if(tree_depth >= FOF_COMPRESS_PATHS_MIN_LENGTH) atomic_cas(&group_index[i], group_index[i], root);

  return root;
}

/* Updates the root and checks that its value has not been changed since being
 * read. */
__attribute__((always_inline)) INLINE static size_t update_root(
    volatile size_t *address, size_t y) {

  size_t *size_t_ptr = (size_t *)address;

  size_t test_val, old_val, new_val;
  old_val = *address;

  test_val = old_val;
  new_val = y;

  /* atomic_cas returns old_val if *size_t_ptr has not changed since being
   * read.*/
  old_val = atomic_cas(size_t_ptr, test_val, new_val);

  if (test_val == old_val)
    return 1;
  else
    return 0;
}

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
      result = update_root(&group_index[root_i_new], root_j_new);

      /* Update root_i on the fly. */
      *root_i = root_j_new;
    } else {

      /* Updates the root and checks that its value has not been changed since
       * being read. */
      result = update_root(&group_index[root_j_new], root_i_new);

      /* Update root_i on the fly. */
      *root_i = root_i_new;
    }
  } while (result != 1);
}

/* Find the shortest distance between cells, remembering to account for boundary
 * conditions. */
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
   * boundary conditions. */
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
  hashmap_value_t *offset = hashmap_get_new(map, group_id, &created_new_element);

  if(offset != NULL) {

    /* If the element is a new entry set its value. */
    if(created_new_element) {
      (*offset).value_st = group_offset;
    }
  }
  else error("Couldn't find key (%zu) or create new one.", group_id);

 }

/* Find a group in the hash table. */
__attribute__((always_inline)) INLINE static size_t hashmap_find_group_offset(
    const size_t group_id, hashmap_t *map) {

  hashmap_value_t *group_offset = hashmap_get(map, group_id);

  if(group_offset == NULL)
    error("Couldn't find key (%zu) or create new one.", group_id);

  return (size_t)(*group_offset).value_st;
}

/* Compute send/recv offsets for MPI communication. */
__attribute__((always_inline)) INLINE static void fof_compute_send_recv_offsets(const int nr_nodes, int *sendcount, int **recvcount, int **sendoffset, int **recvoffset, size_t *nrecv) {

  /* Determine number of entries to receive */
  *recvcount = malloc(nr_nodes * sizeof(int));
  MPI_Alltoall(sendcount, 1, MPI_INT, *recvcount, 1, MPI_INT, MPI_COMM_WORLD);

  /* Compute send/recv offsets */
  *sendoffset = malloc(nr_nodes * sizeof(int));
  
  (*sendoffset)[0] = 0;
  for(int i=1; i<nr_nodes; i++)
    (*sendoffset)[i] = (*sendoffset)[i-1] + sendcount[i-1];
  
  *recvoffset = malloc(nr_nodes * sizeof(int));
  
  (*recvoffset)[0] = 0;
  for(int i=1; i<nr_nodes; i++)
    (*recvoffset)[i] = (*recvoffset)[i-1] + (*recvcount)[i-1];

  /* Allocate receive buffer */
  *nrecv = 0;
  for(int i=0; i<nr_nodes; i++)
    (*nrecv) += (*recvcount)[i];

}
#endif

/* Recurse on a pair of cells and perform a FOF search between cells that are
 * within range. */
void rec_fof_search_pair(struct cell *restrict ci, struct cell *restrict cj,
                         struct space *s, const double dim[3],
                         const double search_r2) {

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

  if (ci == cj) error("Pair FOF called on same cell!!!");

  /* Return if cells are out of range of each other. */
  if (r2 > search_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_pair(ci->progeny[k], cj->progeny[l], s, dim,
                                search_r2);
      }
    }
  } else if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair(ci->progeny[k], cj, s, dim, search_r2);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair(ci, cj->progeny[k], s, dim, search_r2);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells(s, ci, cj);
  }
}

#ifdef WITH_MPI
/* Recurse on a pair of cells (one local, one foreign) and perform a FOF search
 * between cells that are within range. */
static void rec_fof_search_pair_foreign(struct cell *ci, struct cell *cj,
                                        struct space *s, const double *dim,
                                        const double search_r2, int *link_count,
                                        struct fof_mpi **group_links,
                                        int *group_links_size) {

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

  if (ci == cj) error("Pair FOF called on same cell!!!");

  /* Return if cells are out of range of each other. */
  if (r2 > search_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_pair_foreign(ci->progeny[k], cj->progeny[l], s, dim,
                                        search_r2, link_count, group_links,
                                        group_links_size);
      }
    }
  } else if (ci->split) {

    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair_foreign(ci->progeny[k], cj, s, dim, search_r2,
                                    link_count, group_links, group_links_size);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair_foreign(ci, cj->progeny[k], s, dim, search_r2,
                                    link_count, group_links, group_links_size);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells_foreign(s, ci, cj, link_count, group_links,
                                  group_links_size);
  }
}
#endif

/* Recurse on a cell and perform a FOF search between cells that are within
 * range. */
void rec_fof_search_self(struct cell *c, struct space *s, const double dim[3],
                         const double search_r2) {

  /* Recurse? */
  if (c->split) {

    /* Loop over all progeny. Perform pair and self recursion on progenies.*/
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        rec_fof_search_self(c->progeny[k], s, dim, search_r2);

        for (int l = k + 1; l < 8; l++)
          if (c->progeny[l] != NULL)
            rec_fof_search_pair(c->progeny[k], c->progeny[l], s, dim,
                                search_r2);
      }
    }
  }
  /* Otherwise, compute self-interaction. */
  else
    fof_search_cell(s, c);
}

/* Perform a FOF search on a single cell using the Union-Find algorithm.*/
void fof_search_cell(struct space *s, struct cell *c) {

  const size_t count = c->grav.count;
  struct gpart *gparts = c->grav.parts;
  const double l_x2 = s->l_x2;
  size_t *group_index = s->fof_data.group_index;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset = group_index + (ptrdiff_t)(gparts - s->gparts);

  if (c->nodeID != engine_rank)
    error("Performing self FOF search on foreign cell.");

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Find the root of pi. */
    size_t root_i = fof_find(offset[i], group_index);

    for (size_t j = i + 1; j < count; j++) {

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      struct gpart *pj = &gparts[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute the pairwise distance */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) fof_union(&root_i, root_j, group_index);
    }
  }
}

/* Perform a FOF search on a pair of cells using the Union-Find algorithm.*/
void fof_search_pair_cells(struct space *s, struct cell *restrict ci,
                           struct cell *restrict cj) {

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart *gparts_i = ci->grav.parts;
  struct gpart *gparts_j = cj->grav.parts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const float l_x2 = s->l_x2;
  size_t *group_index = s->fof_data.group_index;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset_i = group_index + (ptrdiff_t)(gparts_i - s->gparts);
  size_t *const offset_j = group_index + (ptrdiff_t)(gparts_j - s->gparts);

#ifdef SWIFT_DEBUG_CHECKS
  if (offset_j > offset_i && (offset_j < offset_i + count_i))
    error("Overlapping cells");
  if (offset_i > offset_j && (offset_i < offset_j + count_j))
    error("Overlapping cells");
#endif

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  const int periodic = s->periodic;
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
    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Find the root of pi. */
    size_t root_i = fof_find(offset_i[i], group_index);

    for (size_t j = 0; j < count_j; j++) {

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset_j[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      struct gpart *pj = &gparts_j[j];
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
      if (r2 < l_x2) fof_union(&root_i, root_j, group_index);
    }
  }
}

/* Perform a FOF search between a local and foreign cell using the Union-Find
 * algorithm. Store any links found between particles.*/
void fof_search_pair_cells_foreign(struct space *s, struct cell *ci,
                                   struct cell *cj, int *link_count,
                                   struct fof_mpi **group_links,
                                   int *group_links_size) {

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart *gparts_i = ci->grav.parts;
  struct gpart *gparts_j = cj->grav.parts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;
  size_t *group_index = s->fof_data.group_index;
  size_t *group_size = s->fof_data.group_size;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset_i = group_index + (ptrdiff_t)(gparts_i - s->gparts);

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Check whether cells are local to the node. */
  const int ci_local = (ci->nodeID == engine_rank);
  const int cj_local = (cj->nodeID == engine_rank);

  if ((ci_local && cj_local) || (!ci_local && !cj_local))
    error(
        "FOF search of foreign cells called on two local cells or two foreign "
        "cells.");

  /* Get the relative distance between the pairs, wrapping. */
  const int periodic = s->periodic;
  double diff[3];

  if (ci_local) {

    for (int k = 0; k < 3; k++) {
      diff[k] = cj->loc[k] - ci->loc[k];
      if (periodic && diff[k] < -dim[k] / 2)
        shift[k] = dim[k];
      else if (periodic && diff[k] > dim[k] / 2)
        shift[k] = -dim[k];
      else
        shift[k] = 0.0;
      diff[k] += shift[k];
    }

    /* Loop over particles and find which particles belong in the same group. */
    for (size_t i = 0; i < count_i; i++) {

      struct gpart *pi = &gparts_i[i];
      const double pix = pi->x[0] - shift[0];
      const double piy = pi->x[1] - shift[1];
      const double piz = pi->x[2] - shift[2];

      /* Find the root of pi. */
      const size_t root_i =
          fof_find_global(offset_i[i] - node_offset, group_index, s->nr_gparts);

      for (size_t j = 0; j < count_j; j++) {

        struct gpart *pj = &gparts_j[j];
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
          for (int l = 0; l < *link_count; l++) {
            if ((*group_links)[l].group_i == root_i &&
                (*group_links)[l].group_j == pj->group_id) {
              found = 1;
              break;
            }
          }

          if (!found) {

            /* If the group_links array is not big enough re-allocate it. */
            if (*link_count + 1 > *group_links_size) {

              int new_size = 2 * (*group_links_size);

              *group_links_size = new_size;

              (*group_links) = (struct fof_mpi *)realloc(
                  *group_links, new_size * sizeof(struct fof_mpi));

              message("Re-allocating local group links from %d to %d elements.",
                      *link_count, new_size);
            }

            /* Store the particle group properties for communication. */
            (*group_links)[*link_count].group_i = root_i;
            (*group_links)[*link_count].group_i_size =
                group_size[root_i - node_offset];

            (*group_links)[*link_count].group_j = pj->group_id;
            (*group_links)[*link_count].group_j_size = pj->group_size;

            (*link_count)++;
          }
        }
      }
    }
  } else
    error("Cell ci, is not local.");
}

/* Mapper function to atomically update the group size array. */
void fof_update_group_size_mapper(hashmap_key_t key, hashmap_value_t *value, void *data) {

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
  size_t *group_index = s->fof_data.group_index;
  size_t *group_size = s->fof_data.group_size;

  /* Offset into gparts array. */
  ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);
  size_t *const group_index_offset = group_index + gparts_offset;

  /* Create hash table. */
  hashmap_t map;
  hashmap_init(&map);

  /* Loop over particles and find which cells are in range of each other to perform
   * the FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    hashmap_key_t root = (hashmap_key_t)fof_find(group_index_offset[ind], group_index);
    const size_t gpart_index = gparts_offset + ind;

    /* Only add particles which aren't the root of a group. Stops groups of size 1 being added to the hash table. */ 
    if(root != gpart_index) { 
      hashmap_value_t *size = hashmap_get(&map, root);

      if(size != NULL) (*size).value_st++;
      else error("Couldn't find key (%zu) or create new one.", root);
    }

  }
 
  /* Update the group size array. */
  if (map.size > 0) hashmap_iterate(&map, fof_update_group_size_mapper, group_size);

  hashmap_free(&map);

}

#ifdef WITH_MPI
/* Mapper function to unpack hash table into array. */
void fof_unpack_group_mass_mapper(hashmap_key_t key, hashmap_value_t *value, void *data) {

  struct fof_mass_send_hashmap *fof_mass_send = (struct fof_mass_send_hashmap *)data;
  struct fof_final_mass *mass_send = fof_mass_send->mass_send;
  size_t *nsend = &fof_mass_send->nsend;

  /* Store elements from hash table in array. */
  mass_send[*nsend].global_root = key;
  mass_send[(*nsend)++].group_mass = value->value_dbl;

}
#endif

/* 
 * @brief Calculates the total mass of each group above min_group_size.
 */
void fof_calc_group_mass(struct space *s, const size_t num_groups_local, const size_t num_groups_prev, size_t *num_on_node, size_t *first_on_node, double *group_mass) {

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;
  const size_t group_id_offset = s->fof_data.group_id_offset;
  const size_t group_id_default = s->fof_data.group_id_default;

#ifdef WITH_MPI
  size_t *group_index = s->fof_data.group_index;
  int nr_nodes = s->e->nr_nodes;

  hashmap_t map;
  hashmap_init(&map);

  for(size_t i=0; i<nr_gparts; i++) {

    /* Check if the particle is in a group above the threshold. */
    if(gparts[i].group_id != group_id_default) {

      size_t root = fof_find_global(i, group_index, nr_gparts);

      /* Increment the mass of groups that are local */
      if(is_local(root, nr_gparts)) {

        size_t offset = gparts[i].group_id - group_id_offset - num_groups_prev;

        group_mass[offset] += gparts[i].mass;

      }
      /* Add mass fragments of groups that have a foreign root to a hash table. */
      else {
      
        hashmap_value_t *mass = hashmap_get(&map, (hashmap_key_t)root);

        if(mass != NULL) (*mass).value_dbl += gparts[i].mass;
        else error("Couldn't find key (%zu) or create new one.", root);
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
  if (map.size > 0) hashmap_iterate(&map, fof_unpack_group_mass_mapper, &hashmap_mass_send);
 
  nsend = hashmap_mass_send.nsend;

  if(nsend != map.size) error("No. of mass fragments to send != elements in hash table.");

  hashmap_free(&map);

  /* Sort by global root - this puts the groups in order of which node they're stored on */
  qsort(fof_mass_send, nsend, sizeof(struct fof_final_mass), 
        compare_fof_final_mass_global_root);
 
  /* Determine how many entries go to each node */
  int *sendcount = malloc(nr_nodes*sizeof(int));
  for(int i=0; i<nr_nodes; i+=1)
    sendcount[i] = 0;
  int dest = 0;
  for(size_t i=0;i<nsend;i+=1) {
    while((fof_mass_send[i].global_root >= first_on_node[dest] + num_on_node[dest]) || (num_on_node[dest]==0))
      dest += 1;
    if(dest>=nr_nodes)
      error("Node index out of range!");
    sendcount[dest] += 1;
  } 

  int *recvcount = NULL, *sendoffset = NULL, *recvoffset = NULL;
  size_t nrecv = 0;
  
  fof_compute_send_recv_offsets(nr_nodes, sendcount, &recvcount, &sendoffset, &recvoffset, &nrecv);

  struct fof_final_mass *fof_mass_recv = malloc(nrecv * sizeof(struct fof_final_mass));

  /* Exchange group mass */
  MPI_Alltoallv(fof_mass_send, sendcount, sendoffset, fof_final_mass_type,
                fof_mass_recv, recvcount, recvoffset, fof_final_mass_type,
                MPI_COMM_WORLD);

  /* For each received global root, look up the group ID we assigned and increment the group mass */
  for(size_t i=0; i<nrecv; i+=1) {
    if((fof_mass_recv[i].global_root < node_offset) || (fof_mass_recv[i].global_root >= node_offset+nr_gparts)) {
      error("Received global root index out of range!");
    }
    group_mass[gparts[fof_mass_recv[i].global_root - node_offset].group_id - group_id_offset - num_groups_prev] += fof_mass_recv[i].group_mass;
  }
  
  /* For each received global root, look up the group ID we assigned and send the total group mass back */
  //for(size_t i=0; i<nrecv; i+=1) {
  //  if((fof_mass_recv[i].global_root < node_offset) || (fof_mass_recv[i].global_root >= node_offset+nr_gparts)) {
  //    error("Received global root index out of range!");
  //  }
  //  fof_mass_recv[i].group_mass = group_mass[fof_mass_recv[i].global_root - node_offset];
  //}

  /* Send the result back */
  //MPI_Alltoallv(fof_mass_recv, recvcount, recvoffset, fof_final_mass_type,
  //              fof_mass_send, sendcount, sendoffset, fof_final_mass_type,
  //              MPI_COMM_WORLD);
  //  
  ///* Update local gparts.group_id */
  //for(size_t i=0; i<nsend; i+=1){
  //  if((fof_mass_send[i].local_root < node_offset) || (fof_mass_send[i].local_root >= node_offset + nr_gparts)) {
  //    error("Sent local root index out of range!");
  //  }    
  //  //gparts[fof_mass_send[i].local_root-node_offset].group_id = fof_mass_send[i].global_root;
  //  group_mass[fof_mass_send[i].local_root - node_offset] += fof_mass_send[i].group_mass;
  //}

  MPI_Type_free(&fof_final_mass_type);
  free(sendcount);
  free(recvcount);
  free(sendoffset);
  free(recvoffset);
  free(fof_mass_send);
  free(fof_mass_recv);
#else
  /* Loop over particles and increment the group mass for groups above min_group_size. */
  for(size_t i=0; i<nr_gparts; i++) {
    
    if(gparts[i].group_id != group_id_default) {
      
      group_mass[gparts[i].group_id - group_id_offset] += gparts[i].mass;

    }
  
  }
#endif

}

#ifdef WITH_MPI
/**
 * @brief Mapper function to perform FOF search.
 *
 * @param map_data An array of #cell pair indices.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_find_foreign_links_mapper(void *map_data, int num_elements,
                                   void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct cell_pair_indices *cell_pairs = (struct cell_pair_indices *)map_data;

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;

  /* Store links in an array local to this thread. */
  int local_link_count = 0;
  int local_group_links_size = s->fof_data.group_links_size / s->e->nr_threads;

  /* Init the local group links buffer. */
  struct fof_mpi *local_group_links =
      (struct fof_mpi *)calloc(sizeof(struct fof_mpi), local_group_links_size);
  if (local_group_links == NULL)
    error("Failed to allocate temporary group links buffer.");

  /* Loop over all pairs of local and foreign cells, recurse then perform a
   * FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the local and foreign cells to recurse on. */
    struct cell *restrict local_cell = cell_pairs[ind].local;
    struct cell *restrict foreign_cell = cell_pairs[ind].foreign;

    rec_fof_search_pair_foreign(local_cell, foreign_cell, s, dim, search_r2,
                                &local_link_count, &local_group_links,
                                &local_group_links_size);
  }

  /* Add links found by this thread to the global link list. */
  /* Lock to prevent race conditions while adding to the global link list.*/
  if (lock_lock(&s->lock) == 0) {

    /* Get pointers to global arrays. */
    int *group_links_size = &s->fof_data.group_links_size;
    int *group_link_count = &s->fof_data.group_link_count;
    struct fof_mpi **group_links = &s->fof_data.group_links;

    /* If the global group_links array is not big enough re-allocate it. */
    if (*group_link_count + local_link_count > *group_links_size) {

      const int old_size = *group_links_size; 
      const int new_size = max(*group_link_count + local_link_count, 2 * old_size);

      (*group_links) = (struct fof_mpi *)realloc(
          *group_links, new_size * sizeof(struct fof_mpi));
      
      *group_links_size = new_size;

      message("Re-allocating global group links from %d to %d elements.",
              old_size, new_size);
    }

    /* Copy the local links to the global list. */
    for (int i = 0; i < local_link_count; i++) {

      int found = 0;

      /* Check that the links have not already been added to the list by another thread. */
      for(int l=0; l<*group_link_count; l++) {
        if((*group_links)[l].group_i == local_group_links[i].group_i && (*group_links)[l].group_j == local_group_links[i].group_j) {
          found = 1;
          break;
        }
      }

      if(!found) {

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

  free(local_group_links);
}
#endif

/**
 * @brief Search foreign cells for links and communicate any found to the
 * appropriate node.
 *
 * @param s Pointer to a #space.
 */
void fof_search_foreign_cells(struct space *s) {

#ifdef WITH_MPI

  struct engine *e = s->e;
  size_t *group_index = s->fof_data.group_index;
  size_t *group_size = s->fof_data.group_size;
  const size_t nr_gparts = s->nr_gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;

  message("Searching foreign cells for links.");

  ticks tic = getticks();

  /* Make group IDs globally unique. */
  for (size_t i = 0; i < nr_gparts; i++) group_index[i] += node_offset;

  struct cell_pair_indices *cell_pairs = NULL;
  int group_link_count = 0;
  int cell_pair_count = 0;

  s->fof_data.group_links_size = s->fof_data.group_links_size_default;
  s->fof_data.group_link_count = 0;

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

  if (posix_memalign((void **)&s->fof_data.group_links, SWIFT_STRUCT_ALIGNMENT,
                     s->fof_data.group_links_size * sizeof(struct fof_mpi)) !=
      0)
    error("Error while allocating memory for FOF links over an MPI domain");

  if (posix_memalign((void **)&cell_pairs, SWIFT_STRUCT_ALIGNMENT,
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

  message(
      "Finding local/foreign cell pairs and initialising particle roots took: "
      "%.3f %s.",
      clocks_from_ticks(getticks() - tic), clocks_getunit());

  message("Pairs of touching cells: %d", cell_pair_count);

  tic = getticks();

  struct scheduler *sched = &e->sched;
  struct task *tasks = sched->tasks;

  /* Activate the send and receive tasks for the gparts. */
  for (int i = 0; i < sched->nr_tasks; i++) {

    struct task *t = &tasks[i];

    if (t->type == task_type_send && t->subtype == task_subtype_gpart) {
      scheduler_activate(sched, t);
    }

    if (t->type == task_type_recv && t->subtype == task_subtype_gpart) {
      scheduler_activate(sched, t);
    }
  }

  message("MPI send/recv task activation took: %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Perform send and receive tasks. */
  engine_launch(e);

  message("MPI send/recv comms took: %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Perform search of group links between local and foreign cells with the
   * threadpool. */
  threadpool_map(&s->e->threadpool, fof_find_foreign_links_mapper, cell_pairs,
                 cell_pair_count, sizeof(struct cell_pair_indices), 1, s);

  group_link_count = s->fof_data.group_link_count;

  /* Clean up memory. */
  free(cell_pairs);

  message("Searching for foreign links took: %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());

  message(
      "Rank %d found %d unique group links between local and foreign groups.",
      engine_rank, group_link_count);

  tic = getticks();

  struct fof_mpi *global_group_links = NULL;
  int *displ = NULL, *group_link_counts = NULL;
  int global_group_link_count = 0;

  ticks comms_tic = getticks();

  MPI_Barrier(MPI_COMM_WORLD);

  message("Imbalance took: %.3f %s.", clocks_from_ticks(getticks() - comms_tic),
          clocks_getunit());

  comms_tic = getticks();

  /* Sum the total number of links across MPI domains over each MPI rank. */
  MPI_Allreduce(&group_link_count, &global_group_link_count, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);

  /* Unique set of links is half of all group links as each link is found twice
   * by opposing MPI ranks. */
  if (posix_memalign((void **)&global_group_links, SWIFT_STRUCT_ALIGNMENT,
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
  MPI_Allgatherv(s->fof_data.group_links, group_link_count, fof_mpi_type,
                 global_group_links, group_link_counts, displ, fof_mpi_type,
                 MPI_COMM_WORLD);

  /* Clean up memory. */
  free(displ);
  free(s->fof_data.group_links);

  message("Communication took: %.3f %s.",
          clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

  message("Global comms took: %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());

  tic = getticks();

  /* Transform the group IDs to a local list going from 0-group_count so a
   * union-find can be performed. */
  size_t *global_group_index = NULL, *global_group_id = NULL,
         *global_group_size = NULL;
  const int global_group_list_size = 2 * global_group_link_count;
  int group_count = 0;

  if (posix_memalign((void **)&global_group_index, SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (posix_memalign((void **)&global_group_id, SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (posix_memalign((void **)&global_group_size, SWIFT_STRUCT_ALIGNMENT,
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

  if (posix_memalign((void **)&orig_global_group_size, SWIFT_STRUCT_ALIGNMENT,
                     group_count * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  for (int i = 0; i < group_count; i++)
    orig_global_group_size[i] = global_group_size[i];

  /* Perform a union-find on the group links. */
  for (int i = 0; i < global_group_link_count; i++) {

    /* Use the hash table to find the group offsets in the index array. */
    size_t find_i = hashmap_find_group_offset(global_group_links[i].group_i, &map);
    size_t find_j = hashmap_find_group_offset(global_group_links[i].group_j, &map);

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

  message("Updating groups locally took: %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Clean up memory. */
  free(global_group_links);
  free(global_group_index);
  free(global_group_size);
  free(global_group_id);
  free(orig_global_group_size);

  message("Rank %d finished linking local roots to foreign roots.",
          engine_rank);

#endif /* WITH_MPI */
}

/* Perform a FOF search on gravity particles using the cells and applying the
 * Union-Find algorithm.*/
void fof_search_tree(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t min_group_size = s->fof_data.min_group_size;

  const size_t group_id_offset = s->fof_data.group_id_offset;
#ifdef WITH_MPI
  const int nr_nodes = s->e->nr_nodes;
#endif
  struct gpart *gparts = s->gparts;
  size_t *group_index, *group_size;
  int num_groups = 0, num_parts_in_groups = 0, max_group_size = 0;
  ticks tic_total = getticks();

  char output_file_name[PARSER_MAX_LINE_SIZE];
  snprintf(output_file_name, PARSER_MAX_LINE_SIZE, "%s", s->fof_data.base_name);

  message("Searching %zu gravity particles for links with l_x2: %lf", nr_gparts,
          s->l_x2);
  
  message("Size of hash table element: %ld", sizeof(hashmap_element_t));

  node_offset = 0;

  const size_t group_id_default = s->fof_data.group_id_default;

#ifdef WITH_MPI
  /* Determine number of gparts on lower numbered MPI ranks */
  long long nr_gparts_cumulative;
  long long nr_gparts_local = s->nr_gparts;
  MPI_Scan(&nr_gparts_local, &nr_gparts_cumulative, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); 
  node_offset = nr_gparts_cumulative - nr_gparts_local;

  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
           "_mpi_rank_%d.dat", engine_rank);
#else
  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
           ".dat");
#endif

  group_index = s->fof_data.group_index;
  group_size = s->fof_data.group_size;

  ticks tic_calc_group_size = getticks();

  threadpool_map(&s->e->threadpool, fof_calc_group_size_mapper,
                 gparts, nr_gparts, sizeof(struct gpart), nr_gparts / s->e->nr_threads, s);
  
  message("FOF calc group size took (scaling): %.3f %s.",
      clocks_from_ticks(getticks() - tic_calc_group_size), clocks_getunit());
  
#ifdef WITH_MPI
  if (nr_nodes > 1) {

    ticks tic_mpi = getticks();

    /* Search for group links across MPI domains. */
    fof_search_foreign_cells(s);

    message("fof_search_foreign_cells() took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_mpi), clocks_getunit());
  }
#endif

  message("Calculating group properties...");

  size_t num_groups_local = 0, num_parts_in_groups_local = 0,
         max_group_size_local = 0;

  for (size_t i = 0; i < nr_gparts; i++) {

    /* Find the total number of groups. */
    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size)
      num_groups_local++;

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

  if (posix_memalign((void **)&high_group_sizes, 32,
                     num_groups_local * sizeof(struct group_length)) != 0)
    error("Failed to allocate list of large groups.");

  /* Store the group_sizes and their offset. */
  for (size_t i = 0; i < nr_gparts; i++) {

    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size) {
      high_group_sizes[group_count].index = node_offset + i;
      high_group_sizes[group_count++].size = group_size[i];
    }
  }

  message("Sorting groups...");

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
  s->fof_data.num_groups = num_groups;

  /* Find number of groups on lower numbered MPI ranks */
  size_t num_groups_prev = 0;
#ifdef WITH_MPI
  long long nglocal = num_groups_local;
  long long ngsum;
  MPI_Scan(&nglocal, &ngsum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); /* TODO: define MPI_SIZE_T? */
  num_groups_prev = (size_t) (ngsum-nglocal);
#endif /* WITH_MPI */

  /* Sort local groups into descending order of size */
  qsort(high_group_sizes, num_groups_local, sizeof(struct group_length),
        cmp_func);

  /* Set default group ID for all particles */
  for (size_t i = 0; i < nr_gparts; i++)
    gparts[i].group_id = group_id_default;

  /*
    Assign final group IDs to local root particles where the global root is on this node
    and the group is large enough. Within a node IDs are assigned in descending order of
    particle number.
  */
  for (size_t i = 0; i < num_groups_local; i++)
    gparts[high_group_sizes[i].index-node_offset].group_id = group_id_offset + i + num_groups_prev;

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
     Identify local roots with global root on another node and large enough group_size.
     Store index of the local and global roots in these cases.
  
     NOTE: if group_size only contains the total FoF mass for global roots,
     then we have to communicate ALL fragments where the global root is not
     on this node. Hence the commented out extra conditions below.
  */
  size_t nsend = 0;
  for(size_t i=0; i<nr_gparts;i+=1) {
    if((!is_local(group_index[i], nr_gparts))){ /* && (group_size[i] >= min_group_size)) { */
      nsend += 1;
    }
  }
  struct fof_final_index *fof_index_send = malloc(sizeof(struct fof_final_index)*nsend);
  nsend = 0;
  for(size_t i=0; i<nr_gparts;i+=1) {
    if((!is_local(group_index[i], nr_gparts))){ /* && (group_size[i] >= min_group_size)) { */
      fof_index_send[nsend].local_root  = node_offset + i;
      fof_index_send[nsend].global_root = group_index[i];
      nsend += 1;
    }
  }

  /* Sort by global root - this puts the groups in order of which node they're stored on */
  qsort(fof_index_send, nsend, sizeof(struct fof_final_index), 
        compare_fof_final_index_global_root);

  /* Determine range of global indexes (i.e. particles) on each node */
  size_t *num_on_node = malloc(nr_nodes*sizeof(size_t));
  MPI_Allgather(&nr_gparts,  sizeof(size_t), MPI_BYTE,
                num_on_node, sizeof(size_t), MPI_BYTE,
                MPI_COMM_WORLD);
  size_t *first_on_node = malloc(nr_nodes*sizeof(size_t));
  first_on_node[0] = 0;
  for(int i=1;i<nr_nodes;i+=1)
    first_on_node[i] = first_on_node[i-1] + num_on_node[i-1];

  /* Determine how many entries go to each node */
  int *sendcount = malloc(nr_nodes*sizeof(int));
  for(int i=0; i<nr_nodes; i+=1)
    sendcount[i] = 0;
  int dest = 0;
  for(size_t i=0;i<nsend;i+=1) {
    while((fof_index_send[i].global_root >= first_on_node[dest] + num_on_node[dest]) || (num_on_node[dest]==0))
      dest += 1;
    if(dest>=nr_nodes)
      error("Node index out of range!");
    sendcount[dest] += 1;
  } 
  
  int *recvcount = NULL, *sendoffset = NULL, *recvoffset = NULL;
  size_t nrecv = 0;

  fof_compute_send_recv_offsets(nr_nodes, sendcount, &recvcount, &sendoffset, &recvoffset, &nrecv);

  struct fof_final_index *fof_index_recv = malloc(nrecv*sizeof(struct fof_final_index));

  /* Exchange group indexes */
  MPI_Alltoallv(fof_index_send, sendcount, sendoffset, fof_final_index_type,
                fof_index_recv, recvcount, recvoffset, fof_final_index_type,
                MPI_COMM_WORLD);

  /* For each received global root, look up the group ID we assigned and store it in the struct */
  for(size_t i=0; i<nrecv; i+=1) {
    if((fof_index_recv[i].global_root < node_offset) || (fof_index_recv[i].global_root >= node_offset+nr_gparts)) {
      error("Received global root index out of range!");
    }
    fof_index_recv[i].global_root = gparts[fof_index_recv[i].global_root-node_offset].group_id;
  }

  /* Send the result back */
  MPI_Alltoallv(fof_index_recv, recvcount, recvoffset, fof_final_index_type,
                fof_index_send, sendcount, sendoffset, fof_final_index_type,
                MPI_COMM_WORLD);
    
  /* Update local gparts.group_id */
  for(size_t i=0; i<nsend; i+=1){
    if((fof_index_send[i].local_root < node_offset) || (fof_index_send[i].local_root >= node_offset+nr_gparts)) {
      error("Sent local root index out of range!");
    }    
    gparts[fof_index_send[i].local_root-node_offset].group_id = fof_index_send[i].global_root;
  }

  MPI_Type_free(&fof_final_index_type);
  free(sendcount);
  free(recvcount);
  free(sendoffset);
  free(recvoffset);
  free(fof_index_send);
  free(fof_index_recv);

#endif /* WITH_MPI */

  /* Assign every particle the group_id of its local root. */
  for (size_t i = 0; i < nr_gparts; i++) {
    const size_t root = fof_find_local(i, nr_gparts, group_index);
    gparts[i].group_id = gparts[root].group_id;
  }

  message("Group sorting took: %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());

  double *group_mass;
  /* Allocate and initialise a group mass array. */
  if (posix_memalign((void **)&group_mass, 32,
                     num_groups_local * sizeof(double)) != 0)
    error("Failed to allocate list of group masses for FOF search.");

  bzero(group_mass, num_groups_local * sizeof(double));

#ifdef WITH_MPI
  fof_calc_group_mass(s, num_groups_local, num_groups_prev, num_on_node, first_on_node, group_mass);
  free(num_on_node);
  free(first_on_node);
#else
  fof_calc_group_mass(s, num_groups_local, num_groups_prev, NULL, NULL, group_mass);
#endif

  message("Dumping data...");

  /* Dump group data. */
  fof_dump_group_data(output_file_name, s, num_groups_local, high_group_sizes, group_mass);

  /* Free the left-overs */
  free(high_group_sizes);
  free(group_mass);

  if (engine_rank == 0) {
    message(
        "No. of groups: %d. No. of particles in groups: %d. No. of particles "
        "not in groups: %lld.",
        num_groups, num_parts_in_groups,
        s->e->total_nr_gparts - num_parts_in_groups);

    message("Largest group by size: %d", max_group_size);
  }
  message("FOF search took: %.3f %s.",
          clocks_from_ticks(getticks() - tic_total), clocks_getunit());

#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/* Dump FOF group data. */
void fof_dump_group_data(char *out_file, struct space *s, int num_groups, 
                         struct group_length *group_sizes, double *group_mass) {

  FILE *file = fopen(out_file, "w");

  struct gpart *gparts = s->gparts;
  size_t *group_size = s->fof_data.group_size;

  fprintf(file, "# %8s %10s %10s\n", "Group ID", "Group Size", "Group Mass");
  fprintf(file,
          "#-------------------------------------------------------------------"
          "-------------\n");

  for (int i = 0; i < num_groups; i++) {

    const size_t group_offset = group_sizes[i].index;

    fprintf(file, "  %8zu %10zu %10lf\n",
            gparts[group_offset - node_offset].group_id,
            group_size[group_offset - node_offset], 
            group_mass[i]);
  }

  fclose(file);
}
