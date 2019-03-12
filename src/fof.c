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
#include "c_hashmap/hashmap.c"

#ifdef WITH_MPI
MPI_Datatype fof_mpi_type;
MPI_Datatype group_length_mpi_type;
#endif
size_t node_offset;

typedef struct data_struct_s {
  size_t group_id;
  size_t group_size;
} data_struct_t;

#define UNION_BY_SIZE_OVER_MPI 1

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

  /* Allocate and initialise a group mass array. */
  if (posix_memalign((void **)&s->fof_data.group_mass, 32,
                     nr_local_gparts * sizeof(double)) != 0)
    error("Failed to allocate list of group masses for FOF search.");

  /* Allocate and initialise a group mass array. */
  if (posix_memalign((void **)&s->fof_data.group_CoM, 32,
                     nr_local_gparts * sizeof(struct fof_CoM)) != 0)
    error("Failed to allocate list of group CoM for FOF search.");

  /* Set initial group index to particle offset into array and set default group
   * ID. */
  for (size_t i = 0; i < nr_local_gparts; i++) {
    s->fof_data.group_index[i] = i;
    s->gparts[i].group_id = default_id;
  }

  bzero(s->fof_data.group_size, nr_local_gparts * sizeof(size_t));
  bzero(s->fof_data.group_mass, nr_local_gparts * sizeof(double));
  bzero(s->fof_data.group_CoM, nr_local_gparts * sizeof(struct fof_CoM));

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

/* Finds the global root ID of the group a particle exists in. */
__attribute__((always_inline)) INLINE static size_t fof_find_global(
    const size_t i, size_t *group_index) {

  size_t root = node_offset + i;
  if (root < node_offset) return root;

  while (root != group_index[root - node_offset]) {
    root = group_index[root - node_offset];
    if (root < node_offset) break;
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

/* Finds the local root ID of the group a particle exists in. */
__attribute__((always_inline)) INLINE static size_t fof_find(
    const size_t i, size_t *group_index) {

  size_t root = i;
  while (root != group_index[root]) {
#ifdef PATH_HALVING
    atomic_cas(&group_index[root], group_index[root],
               group_index[group_index[root]]);
#endif
    root = group_index[root];
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

/* Checks whether the group is on the local node. */
__attribute__((always_inline)) INLINE static int is_local(
    const size_t group_id, const size_t nr_gparts) {
  return (group_id >= node_offset && group_id < node_offset + nr_gparts);
}

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
  double *group_mass = s->fof_data.group_mass;
  struct fof_CoM *group_CoM = s->fof_data.group_CoM;

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
          fof_find_global(offset_i[i] - node_offset, group_index);

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

              int new_size =
                  *group_links_size + ceil(0.1 * (double)(*group_links_size));

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
            (*group_links)[*link_count].group_i_mass =
                group_mass[root_i - node_offset];
            (*group_links)[*link_count].group_i_CoM.x =
                group_CoM[root_i - node_offset].x;
            (*group_links)[*link_count].group_i_CoM.y =
                group_CoM[root_i - node_offset].y;
            (*group_links)[*link_count].group_i_CoM.z =
                group_CoM[root_i - node_offset].z;

            (*group_links)[*link_count].group_j = pj->group_id;
            (*group_links)[*link_count].group_j_size = pj->group_size;
            (*group_links)[*link_count].group_j_mass = pj->group_mass;
            (*group_links)[*link_count].group_j_CoM.x = pj->group_CoM.x;
            (*group_links)[*link_count].group_j_CoM.y = pj->group_CoM.y;
            (*group_links)[*link_count].group_j_CoM.z = pj->group_CoM.z;

            (*link_count)++;
          }
        }
      }
    }
  } else
    error("Cell ci, is not local.");
}

/**
 * @brief Mapper function to perform FOF search.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_search_tree_mapper(void *map_data, int num_elements,
                            void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  int *local_cells = (int *)map_data;

  const size_t nr_local_cells = s->nr_local_cells;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;

  /* Make a list of cell offsets into the local top-level cell array. */
  int *const offset =
      s->cell_index + (ptrdiff_t)(local_cells - s->local_cells_top);

  /* Loop over cells and find which cells are in range of each other to perform
   * the FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell. */
    struct cell *restrict ci = &s->cells_top[local_cells[ind]];

    /* Only perform FOF search on local cells. */
    if (ci->nodeID == engine_rank) {

      /* Skip empty cells. */
      if (ci->grav.count == 0) continue;

      /* Perform FOF search on local particles within the cell. */
      rec_fof_search_self(ci, s, dim, search_r2);

      /* Loop over all top-level cells skipping over the cells already searched.
       */
      for (size_t cjd = offset[ind] + 1; cjd < nr_local_cells; cjd++) {

        struct cell *restrict cj = &s->cells_top[s->local_cells_top[cjd]];

        /* Only perform FOF search on local cells. */
        if (cj->nodeID == engine_rank) {

          /* Skip empty cells. */
          if (cj->grav.count == 0) continue;

          rec_fof_search_pair(ci, cj, s, dim, search_r2);
        }
      }
    }
  }
}

/**
 * @brief Mapper function to perform FOF search.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_calc_group_props_mapper(void *map_data, int num_elements,
                            void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct gpart *gparts = (struct gpart *)map_data;
  size_t *group_index = s->fof_data.group_index;
  size_t *group_size = s->fof_data.group_size;
  //double *group_mass = s->fof_data.group_mass;
  //struct fof_CoM *group_CoM = s->fof_data.group_CoM;

  size_t *const offset = s->gpart_index + (ptrdiff_t)(gparts - s->gparts);
  //const ptrdiff_t offset = (ptrdiff_t)(gparts - s->gparts);

  map_t mymap;
  data_struct_t* value;
  
  mymap = hashmap_new();

  /* Loop over particles and find which cells are in range of each other to perform
   * the FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    size_t root = fof_find(offset[ind], group_index);
    //double x = gparts[ind].x[0];
    //double y = gparts[ind].x[1];
    //double z = gparts[ind].x[2];
    //double mass = gparts[ind].mass;
      
    int error = hashmap_get(mymap, (char *)&root, (void**)(&value));
    if(error == MAP_MISSING) {
      value = malloc(sizeof(data_struct_t));
      value->group_id = root;
      value->group_size = 1;
      error = hashmap_put(mymap, (char *)&value->group_id, value);
      if(error != MAP_OK) error("Issue with hash table. Error code: %d", error);
    }
    else if (error == MAP_OK) {
      value->group_size++;
    }
    else {
      error("Issue with hash table. Error code: %d", error);
    }

    /* Use the CoM location set by the first particle added to the group. */
    //if (com_set[root]) {

    //  /* Periodically wrap particle positions if they are located on the other
    //   * side of the domain than the CoM location. */
    //  if (group_bc[root].x > 0.5 * dim[0] && x < 0.5 * dim[0])
    //    x += dim[0];
    //  else if (group_bc[root].x == 0.0 && x > 0.5 * dim[0])
    //    x -= dim[0];

    //  if (group_bc[root].y > 0.5 * dim[1] && y < 0.5 * dim[1])
    //    y += dim[1];
    //  else if (group_bc[root].y == 0.0 && y > 0.5 * dim[1])
    //    y -= dim[1];

    //  if (group_bc[root].z > 0.5 * dim[2] && z < 0.5 * dim[2])
    //    z += dim[2];
    //  else if (group_bc[root].z == 0.0 && z > 0.5 * dim[2])
    //    z -= dim[2];

    //} else {

    //  com_set[root] = 1;

    //  /* Use the first particle to set the CoM location in the domain. */
    //  if (x > 0.5 * dim[0])
    //    group_bc[root].x = dim[0];
    //  else
    //    group_bc[root].x = 0.0;

    //  if (y > 0.5 * dim[1])
    //    group_bc[root].y = dim[1];
    //  else
    //    group_bc[root].y = 0.0;

    //  if (z > 0.5 * dim[2])
    //    group_bc[root].z = dim[2];
    //  else
    //    group_bc[root].z = 0.0;
    //}

    //atomic_add_d(&group_CoM[root].x, mass * x);
    //atomic_add_d(&group_CoM[root].y, mass * y);
    //atomic_add_d(&group_CoM[root].z, mass * z);
    
  }
  
  hashmap_map *m = (hashmap_map *) mymap;

	if (hashmap_length(m) <= 0)
		error("Empty");	

  //for(int i=0; i<hashmap_length(mymap); i++) {
  for(int i=0; i<m->table_size; i++) {

    if(m->data[i].in_use != 0) { 
    //if(m->in_use[i] != 0) { 
      
      data_struct_t* element = NULL;
      element = (data_struct_t *) m->data[i].data;

      atomic_add(&group_size[element->group_id], element->group_size);

      free(element);
    }
  }
    
  hashmap_free(mymap);

}

/**
 * @brief Mapper function to perform FOF search.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_calc_group_props_new_mapper(void *map_data, int num_elements,
                            void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  size_t *local_roots = (size_t *)map_data;
  struct gpart *gparts = s->gparts;
  size_t nr_gparts = s->nr_gparts;
  size_t *group_index = s->fof_data.group_index;
  size_t *group_size = s->fof_data.group_size;
  double *group_mass = s->fof_data.group_mass;
  struct fof_CoM *group_CoM = s->fof_data.group_CoM;
    
  /* Loop over particles and find which cells are in range of each other to perform
   * the FOF search. */
  for (size_t i = 0; i < nr_gparts; i++) {

    size_t root = fof_find(i, group_index);

    for (int ind = 0; ind < num_elements; ind++) {

      if(local_roots[ind] == root) {

        double x = gparts[i].x[0];
        double y = gparts[i].x[1];
        double z = gparts[i].x[2];
        double mass = gparts[i].mass;

        group_size[root]++;
        group_mass[root] += mass;

        /* Use the CoM location set by the first particle added to the group. */
        //if (com_set[root]) {

        //  /* Periodically wrap particle positions if they are located on the other
        //   * side of the domain than the CoM location. */
        //  if (group_bc[root].x > 0.5 * dim[0] && x < 0.5 * dim[0])
        //    x += dim[0];
        //  else if (group_bc[root].x == 0.0 && x > 0.5 * dim[0])
        //    x -= dim[0];

        //  if (group_bc[root].y > 0.5 * dim[1] && y < 0.5 * dim[1])
        //    y += dim[1];
        //  else if (group_bc[root].y == 0.0 && y > 0.5 * dim[1])
        //    y -= dim[1];

        //  if (group_bc[root].z > 0.5 * dim[2] && z < 0.5 * dim[2])
        //    z += dim[2];
        //  else if (group_bc[root].z == 0.0 && z > 0.5 * dim[2])
        //    z -= dim[2];

        //} else {

        //  com_set[root] = 1;

        //  /* Use the first particle to set the CoM location in the domain. */
        //  if (x > 0.5 * dim[0])
        //    group_bc[root].x = dim[0];
        //  else
        //    group_bc[root].x = 0.0;

        //  if (y > 0.5 * dim[1])
        //    group_bc[root].y = dim[1];
        //  else
        //    group_bc[root].y = 0.0;

        //  if (z > 0.5 * dim[2])
        //    group_bc[root].z = dim[2];
        //  else
        //    group_bc[root].z = 0.0;
        //}

        group_CoM[root].x += mass * x;
        group_CoM[root].y += mass * y;
        group_CoM[root].z += mass * z;
      }
    }
  }
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

      int new_size =
          *group_links_size + ceil(0.1 * (double)(*group_links_size));

      *group_links_size = new_size;

      (*group_links) = (struct fof_mpi *)realloc(
          *group_links, new_size * sizeof(struct fof_mpi));

      message("Re-allocating global group links from %d to %d elements.",
              *group_link_count, new_size);
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
        (*group_links)[*group_link_count].group_i_mass =
          local_group_links[i].group_i_mass;
        (*group_links)[*group_link_count].group_i_CoM.x =
          local_group_links[i].group_i_CoM.x;
        (*group_links)[*group_link_count].group_i_CoM.y =
          local_group_links[i].group_i_CoM.y;
        (*group_links)[*group_link_count].group_i_CoM.z =
          local_group_links[i].group_i_CoM.z;

        (*group_links)[*group_link_count].group_j =
          local_group_links[i].group_j;
        (*group_links)[*group_link_count].group_j_size =
          local_group_links[i].group_j_size;
        (*group_links)[*group_link_count].group_j_mass =
          local_group_links[i].group_j_mass;
        (*group_links)[*group_link_count].group_j_CoM.x =
          local_group_links[i].group_j_CoM.x;
        (*group_links)[*group_link_count].group_j_CoM.y =
          local_group_links[i].group_j_CoM.y;
        (*group_links)[*group_link_count].group_j_CoM.z =
          local_group_links[i].group_j_CoM.z;

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
 * @param local_roots List of local roots.
 */
size_t fof_search_foreign_cells(struct space *s, size_t **local_roots) {

  size_t local_root_count = 0;

#ifdef WITH_MPI

  struct engine *e = s->e;
  const size_t min_group_size = s->fof_data.min_group_size;
  size_t *group_index = s->fof_data.group_index;
  size_t *group_size = s->fof_data.group_size;
  double *group_mass = s->fof_data.group_mass;
  struct fof_CoM *group_CoM = s->fof_data.group_CoM;
  const size_t nr_gparts = s->nr_gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;

  message("Searching foreign cells for links.");

  ticks tic = getticks();

  /* Make group IDs globally unique. */
  for (size_t i = 0; i < nr_gparts; i++) group_index[i] += node_offset;

  struct fof_mpi *group_links = s->fof_data.group_links;
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

  group_links = s->fof_data.group_links;

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
            fof_find_global(offset[k] - node_offset, group_index);
        gparts[k].group_id = root;
        gparts[k].group_size = group_size[root - node_offset];
        gparts[k].group_mass = group_mass[root - node_offset];
        gparts[k].group_CoM.x = group_CoM[root - node_offset].x;
        gparts[k].group_CoM.y = group_CoM[root - node_offset].y;
        gparts[k].group_CoM.z = group_CoM[root - node_offset].z;
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
  MPI_Allgatherv(group_links, group_link_count, fof_mpi_type,
                 global_group_links, group_link_counts, displ, fof_mpi_type,
                 MPI_COMM_WORLD);

  /* Clean up memory. */
  free(displ);
  free(group_links);

  message("Communication took: %.3f %s.",
          clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

  message("Global comms took: %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());

  tic = getticks();

  /* Transform the group IDs to a local list going from 0-group_count so a
   * union-find can be performed. */
  size_t *global_group_index = NULL, *global_group_id = NULL,
         *global_group_size = NULL;
  double *global_group_mass = NULL;
  struct fof_CoM *global_group_CoM = NULL;
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

  if (posix_memalign((void **)&global_group_mass, SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(double)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (posix_memalign((void **)&global_group_CoM, SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(struct fof_CoM)) != 0)
    error("Error while allocating memory for the global group CoM.");

  bzero(global_group_size, global_group_list_size * sizeof(size_t));
  bzero(global_group_mass, global_group_list_size * sizeof(double));
  bzero(global_group_CoM, global_group_list_size * sizeof(struct fof_CoM));

  /* Store each group ID and its properties. */
  for (int i = 0; i < global_group_link_count; i++) {

    size_t group_i = global_group_links[i].group_i;
    size_t group_j = global_group_links[i].group_j;

    global_group_size[group_count] += global_group_links[i].group_i_size;
    global_group_mass[group_count] += global_group_links[i].group_i_mass;
    global_group_CoM[group_count].x += global_group_links[i].group_i_CoM.x;
    global_group_CoM[group_count].y += global_group_links[i].group_i_CoM.y;
    global_group_CoM[group_count].z += global_group_links[i].group_i_CoM.z;
    global_group_id[group_count++] = group_i;

    global_group_size[group_count] += global_group_links[i].group_j_size;
    global_group_mass[group_count] += global_group_links[i].group_j_mass;
    global_group_CoM[group_count].x += global_group_links[i].group_j_CoM.x;
    global_group_CoM[group_count].y += global_group_links[i].group_j_CoM.y;
    global_group_CoM[group_count].z += global_group_links[i].group_j_CoM.z;
    global_group_id[group_count++] = group_j;
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

    int find_i = 0, find_j = 0;

    /* TODO: Sort array to get offset or allocate array for initial root
     * positions. */
    /* Find group offset into global_group_id */
    for (int j = 0; j < group_count; j++) {
      if (global_group_id[j] == global_group_links[i].group_i) {
        find_i = j;
        break;
      }
    }

    for (int j = 0; j < group_count; j++) {
      if (global_group_id[j] == global_group_links[i].group_j) {
        find_j = j;
        break;
      }
    }

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

  message("global_group_index construction took: %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  if (posix_memalign((void **)local_roots, SWIFT_STRUCT_ALIGNMENT,
                     group_count * sizeof(size_t)) != 0)
    error("Error while allocating memory for local roots");

  /* Update each group locally with new root information. */
  for (int i = 0; i < group_count; i++) {

    size_t group_id = global_group_id[i];
    size_t offset = fof_find(global_group_index[i], global_group_index);
    size_t new_root = global_group_id[offset];

    /* If the group is local update its root and size. */
    if (is_local(group_id, nr_gparts) && new_root != group_id) {
      group_index[group_id - node_offset] = new_root;

      group_mass[group_id - node_offset] -= global_group_mass[i];
      group_CoM[group_id - node_offset].x -= global_group_CoM[i].x;
      group_CoM[group_id - node_offset].y -= global_group_CoM[i].y;
      group_CoM[group_id - node_offset].z -= global_group_CoM[i].z;
      group_size[group_id - node_offset] -= orig_global_group_size[i];

      /* Store local root if the group is big enough. */
      if (global_group_size[offset] >= min_group_size)
        (*local_roots)[local_root_count++] = group_id - node_offset;
    }

    /* If the group linked to a local root update its size. */
    if (is_local(new_root, nr_gparts) && new_root != group_id) {

      /* Calculate the CoM of each distinct group. */
      double CoM_x_old = group_CoM[new_root - node_offset].x /
                         group_mass[new_root - node_offset];
      double CoM_y_old = group_CoM[new_root - node_offset].y /
                         group_mass[new_root - node_offset];
      double CoM_z_old = group_CoM[new_root - node_offset].z /
                         group_mass[new_root - node_offset];
      double CoM_x = global_group_CoM[i].x / global_group_mass[i];
      double CoM_y = global_group_CoM[i].y / global_group_mass[i];
      double CoM_z = global_group_CoM[i].z / global_group_mass[i];

      /* Periodically wrap the CoM of the group being linked across the domain
       * based upon the location of the other CoM location. */
      if (CoM_x_old > 0.5 * dim[0] && CoM_x < 0.5 * dim[0])
        CoM_x += dim[0];
      else if (CoM_x_old <= 0.5 * dim[0] && CoM_x > 0.5 * dim[0])
        CoM_x -= dim[0];

      if (CoM_y_old > 0.5 * dim[1] && CoM_y < 0.5 * dim[1])
        CoM_y += dim[1];
      else if (CoM_y_old <= 0.5 * dim[1] && CoM_y > 0.5 * dim[1])
        CoM_y -= dim[1];

      if (CoM_z_old > 0.5 * dim[2] && CoM_z < 0.5 * dim[2])
        CoM_z += dim[2];
      else if (CoM_z_old <= 0.5 * dim[2] && CoM_z > 0.5 * dim[2])
        CoM_z -= dim[2];

      /* Update the CoM. */
      group_CoM[new_root - node_offset].x += (CoM_x * global_group_mass[i]);
      group_CoM[new_root - node_offset].y += (CoM_y * global_group_mass[i]);
      group_CoM[new_root - node_offset].z += (CoM_z * global_group_mass[i]);

      /* Use group sizes before Union-Find */
      group_size[new_root - node_offset] += orig_global_group_size[i];
      group_mass[new_root - node_offset] += global_group_mass[i];
    }
  }

  message("Updating groups locally took: %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Clean up memory. */
  free(global_group_links);
  free(global_group_index);
  free(global_group_size);
  free(global_group_mass);
  free(global_group_id);
  free(orig_global_group_size);

  message("Rank %d finished linking local roots to foreign roots.",
          engine_rank);

#endif /* WITH_MPI */

  return local_root_count;
}

/* Perform a FOF search on gravity particles using the cells and applying the
 * Union-Find algorithm.*/
void fof_search_tree(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t min_group_size = s->fof_data.min_group_size;

  const size_t group_id_offset = s->fof_data.group_id_offset;
  const int nr_nodes = s->e->nr_nodes;
  struct gpart *gparts = s->gparts;
  size_t *group_index, *group_size;
  double *group_mass;
  struct fof_CoM *group_CoM;
  int num_groups = 0, num_parts_in_groups = 0, max_group_size = 0;
  double max_group_mass = 0;
  ticks tic_total = getticks();

  char output_file_name[PARSER_MAX_LINE_SIZE];
  snprintf(output_file_name, PARSER_MAX_LINE_SIZE, "%s", s->fof_data.base_name);

  message("Searching %zu gravity particles for links with l_x2: %lf", nr_gparts,
          s->l_x2);
  
  message("Size of hash table element: %ld", sizeof(data_struct_t));

  node_offset = 0;

#ifdef WITH_MPI
  const size_t group_id_default = s->fof_data.group_id_default;
  if (nr_nodes > 1) {

    int *global_nr_gparts = NULL;

    /* Find the total no. of particles on each node and calculate your unique
     * offset. */
    if (posix_memalign((void **)&global_nr_gparts, SWIFT_STRUCT_ALIGNMENT,
                       nr_nodes * sizeof(int)) != 0)
      error("Error while allocating memory for global_nr_gparts array.");

    bzero(global_nr_gparts, nr_nodes * sizeof(int));

    MPI_Allgather(&s->nr_gparts, 1, MPI_INT, global_nr_gparts, 1, MPI_INT,
                  MPI_COMM_WORLD);

    if (engine_rank == 0) {
      printf("No. of particles: [%d", global_nr_gparts[0]);
      for (int i = 1; i < nr_nodes; i++) printf(", %d", global_nr_gparts[i]);
      printf("]\n");
    }

    for (int i = 0; i < engine_rank; i++) node_offset += global_nr_gparts[i];

    /* Clean up memory. */
    free(global_nr_gparts);
  }

  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
           "_mpi_rank_%d.dat", engine_rank);
#else
  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
           ".dat");
#endif

  group_index = s->fof_data.group_index;
  group_size = s->fof_data.group_size;
  group_mass = s->fof_data.group_mass;
  group_CoM = s->fof_data.group_CoM;

  ticks tic = getticks();

  engine_maketasks(s->e);
  engine_marktasks(s->e);
  engine_print_task_counts(s->e);

  message("Making FOF tasks took: %.3f %s.", clocks_from_ticks(getticks() - tic),
      clocks_getunit());

  tic = getticks();

  engine_launch(s->e);

  /* Perform local FOF using the threadpool. */
  //threadpool_map(&s->e->threadpool, fof_search_tree_mapper,
  //s->local_cells_top,
  //              s->nr_local_cells, sizeof(int), 1, s);

  message("Local FOF took: %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());

  //size_t *local_roots = NULL;

  //if (posix_memalign((void **)&local_roots, SWIFT_STRUCT_ALIGNMENT,
  //                   nr_gparts * sizeof(size_t)) != 0)
  //  error("Error while allocating memory for local roots");

  //size_t root_count = 0;

  //for (size_t i = 0; i < nr_gparts; i++) {
  //  if(group_index[i] == i) local_roots[root_count++] = i;
  //}

  threadpool_map(&s->e->threadpool, fof_calc_group_props_mapper,
                 gparts, nr_gparts, sizeof(struct gpart), nr_gparts / s->e->nr_threads, s);

  //threadpool_map(&s->e->threadpool, fof_calc_group_props_new_mapper,
  //               local_roots, root_count, sizeof(size_t), 1, s);

  //struct fof_CoM *group_bc = NULL;
  //int *com_set = NULL;
  //const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  ///* Allocate and initialise a group boundary condition array. */
  //if (posix_memalign((void **)&group_bc, 32,
  //                   nr_gparts * sizeof(struct fof_CoM)) != 0)
  //  error("Failed to allocate list of group CoM for FOF search.");

  //if (posix_memalign((void **)&com_set, 32, nr_gparts * sizeof(int)) != 0)
  //  error("Failed to allocate list of whether the CoM has been initialised.");

  //bzero(group_bc, nr_gparts * sizeof(struct fof_CoM));
  //bzero(com_set, nr_gparts * sizeof(int));

  ///* Calculate the total number of particles in each group, group mass, group ID
  // and group CoM. */
  //for (size_t i = 0; i < nr_gparts; i++) {
  //  size_t root = fof_find(i, group_index);

  //  group_size[root]++;
  //  group_mass[root] += gparts[i].mass;

  //  double x = gparts[i].x[0];
  //  double y = gparts[i].x[1];
  //  double z = gparts[i].x[2];

  //  /* Use the CoM location set by the first particle added to the group. */
  //  if (com_set[root]) {

  //    /* Periodically wrap particle positions if they are located on the other
  //     * side of the domain than the CoM location. */
  //    if (group_bc[root].x > 0.5 * dim[0] && x < 0.5 * dim[0])
  //      x += dim[0];
  //    else if (group_bc[root].x == 0.0 && x > 0.5 * dim[0])
  //      x -= dim[0];

  //    if (group_bc[root].y > 0.5 * dim[1] && y < 0.5 * dim[1])
  //      y += dim[1];
  //    else if (group_bc[root].y == 0.0 && y > 0.5 * dim[1])
  //      y -= dim[1];

  //    if (group_bc[root].z > 0.5 * dim[2] && z < 0.5 * dim[2])
  //      z += dim[2];
  //    else if (group_bc[root].z == 0.0 && z > 0.5 * dim[2])
  //      z -= dim[2];

  //  } else {

  //    com_set[root] = 1;

  //    /* Use the first particle to set the CoM location in the domain. */
  //    if (x > 0.5 * dim[0])
  //      group_bc[root].x = dim[0];
  //    else
  //      group_bc[root].x = 0.0;

  //    if (y > 0.5 * dim[1])
  //      group_bc[root].y = dim[1];
  //    else
  //      group_bc[root].y = 0.0;

  //    if (z > 0.5 * dim[2])
  //      group_bc[root].z = dim[2];
  //    else
  //      group_bc[root].z = 0.0;
  //  }

  //  group_CoM[root].x += gparts[i].mass * x;
  //  group_CoM[root].y += gparts[i].mass * y;
  //  group_CoM[root].z += gparts[i].mass * z;
  //}

#ifdef WITH_MPI
  size_t num_local_roots = 0;
  size_t *local_roots = NULL;
  if (nr_nodes > 1) {

    ticks tic_mpi = getticks();

    /* Search for group links across MPI domains. */
    num_local_roots = fof_search_foreign_cells(s, &local_roots);

    message("fof_search_foreign_cells() took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_mpi), clocks_getunit());
  }
#endif

  message("Calculating group properties...");

  size_t num_groups_local = 0, num_parts_in_groups_local = 0,
         max_group_size_local = 0;
  double max_group_mass_local = 0;

  for (size_t i = 0; i < nr_gparts; i++) {

    /* Find the group's CoM. */
    if (group_size[i] >= min_group_size) {
      group_CoM[i].x = group_CoM[i].x / group_mass[i];
      group_CoM[i].y = group_CoM[i].y / group_mass[i];
      group_CoM[i].z = group_CoM[i].z / group_mass[i];
    }

    /* Find the total number of groups. */
    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size)
      num_groups_local++;

    /* Find the total number of particles in groups. */
    if (group_size[i] >= min_group_size)
      num_parts_in_groups_local += group_size[i];

    /* Find the largest group. */
    if (group_size[i] > max_group_size_local)
      max_group_size_local = group_size[i];

    /* Find the largest group by mass. */
    if (group_mass[i] > max_group_mass_local)
      max_group_mass_local = group_mass[i];
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

  tic = getticks();

  /* Find global properties. */
#ifdef WITH_MPI
  MPI_Allreduce(&num_groups_local, &num_groups, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Reduce(&num_parts_in_groups_local, &num_parts_in_groups, 1, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_group_size_local, &max_group_size, 1, MPI_INT, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&max_group_mass_local, &max_group_mass, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
#else
  num_groups += num_groups_local;
  num_parts_in_groups += num_parts_in_groups_local;
  max_group_size = max_group_size_local;
  max_group_mass = max_group_mass_local;
#endif /* WITH_MPI */

  /* For non-MPI runs and MPI runs with 1 rank sort the groups locally. */
  if (nr_nodes == 1) {
    /* Sort the groups */
    qsort(high_group_sizes, num_groups_local, sizeof(struct group_length),
          cmp_func);

    /* Update group roots with new value. */
    for (size_t i = 0; i < num_groups_local; i++)
      gparts[high_group_sizes[i].index].group_id = group_id_offset + i;

    /* Re-label the groups between group_id_offset - num_groups. */
    for (size_t i = 0; i < nr_gparts; i++) {
      const size_t root = fof_find(i, group_index);

      if (group_size[root] >= min_group_size)
        gparts[i].group_id = gparts[root].group_id;
    }
  }

#ifdef WITH_MPI
  struct group_length *global_group_sizes = NULL;

  if (nr_nodes > 1) {

    ticks sort_tic = getticks();

    int *displ = NULL, *group_counts = NULL;

    if (posix_memalign((void **)&global_group_sizes, 32,
                       num_groups * sizeof(struct group_length)) != 0)
      error("Failed to allocate list of global large groups.");

    if (posix_memalign((void **)&group_counts, SWIFT_STRUCT_ALIGNMENT,
                       nr_nodes * sizeof(int)) != 0)
      error(
          "Error while allocating memory for the number of groups on each MPI "
          "rank");

    if (posix_memalign((void **)&displ, SWIFT_STRUCT_ALIGNMENT,
                       nr_nodes * sizeof(int)) != 0)
      error(
          "Error while allocating memory for the displacement in memory for "
          "the global group size list");

    /* Gather the total number of links on each rank. */
    MPI_Allgather(&num_groups_local, 1, MPI_INT, group_counts, 1, MPI_INT,
                  MPI_COMM_WORLD);

    /* Set the displacements into the global link list using the link counts
     * from each rank */
    displ[0] = 0;
    for (int i = 1; i < nr_nodes; i++)
      displ[i] = displ[i - 1] + group_counts[i - 1];

    /* Gather the global link list on all ranks. */
    MPI_Allgatherv(high_group_sizes, num_groups_local, group_length_mpi_type,
                   global_group_sizes, group_counts, displ,
                   group_length_mpi_type, MPI_COMM_WORLD);

    message("Global sorting comms took: %.3f %s.",
            clocks_from_ticks(getticks() - sort_tic), clocks_getunit());

    sort_tic = getticks();

    /* Sort the groups */
    qsort(global_group_sizes, num_groups, sizeof(struct group_length),
          cmp_func);

    message("qsort took: %.3f %s.", clocks_from_ticks(getticks() - sort_tic),
            clocks_getunit());

    sort_tic = getticks();

    ticks tic1 = getticks();

    /* Set default group ID again. */
    for (size_t i = 0; i < nr_gparts; i++)
      gparts[i].group_id = group_id_default;

    /* Re-label the groups between group_id_offset - num_groups. */
    for (int i = 0; i < num_groups; i++) {
      const size_t root = global_group_sizes[i].index;

      /* Only update local roots. */
      if (is_local(root, nr_gparts))
        gparts[root - node_offset].group_id = group_id_offset + i;
    }

    message("Re-labelling roots took: %.3f %s.",
            clocks_from_ticks(getticks() - tic1), clocks_getunit());

    ticks tic2 = getticks();

    /* Update local roots that have foreign roots. */
    for (size_t i = 0; i < num_local_roots; i++) {
      size_t root = local_roots[i];

      /* For particles that have foreign roots update their group ID. */
      if (!is_local(group_index[root], nr_gparts)) {

        for (int j = 0; j < num_groups; j++) {
          if (group_index[root] == global_group_sizes[j].index) {
            gparts[root].group_id = group_id_offset + j;
            break;
          }
        }
      }
    }

    message("Updating local roots that have foreign took: %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

    ticks tic3 = getticks();

    /* Update group ID of all particles that are not the root of a group. */
    for (size_t i = 0; i < nr_gparts; i++) {
      size_t root = node_offset + i;
      size_t local_root = node_offset + i;

      /* FOF find */
      while (root != group_index[root - node_offset]) {
        root = group_index[root - node_offset];
        if (!is_local(root, nr_gparts)) break;

        /* Store local root before index becomes foreign. */
        local_root = root;
      }

      /* For particles that have foreign roots use local root to set group ID.
       */
      gparts[i].group_id = gparts[local_root - node_offset].group_id;
    }

    message("Setting all part group IDs took: %.3f %s.",
            clocks_from_ticks(getticks() - tic3), clocks_getunit());

    message("Setting group IDs took: %.3f %s.",
            clocks_from_ticks(getticks() - sort_tic), clocks_getunit());

    /* Clean up memory. */
    free(group_counts);
    free(displ);
    free(local_roots);
  }
#endif

  message("Group sorting took: %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());

  s->fof_data.num_groups = num_groups;

  message("Dumping data...");

  /* Dump group data. */
#ifdef WITH_MPI
  if (nr_nodes > 1) {
    fof_dump_group_data(output_file_name, s, global_group_sizes);
    free(global_group_sizes);
  } else
    fof_dump_group_data(output_file_name, s, high_group_sizes);
#else
  fof_dump_group_data(output_file_name, s, high_group_sizes);
#endif

  if (engine_rank == 0) {
    message(
        "No. of groups: %d. No. of particles in groups: %d. No. of particles "
        "not "
        "in groups: %lld.",
        num_groups, num_parts_in_groups,
        s->e->total_nr_gparts - num_parts_in_groups);

    message("Largest group by size: %d", max_group_size);
    message("Largest group by mass: %e", max_group_mass);
  }
  message("FOF search took: %.3f %s.",
          clocks_from_ticks(getticks() - tic_total), clocks_getunit());

#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/* Dump FOF group data. */
void fof_dump_group_data(char *out_file, struct space *s,
                         struct group_length *group_sizes) {

  FILE *file = fopen(out_file, "w");
  const size_t nr_gparts = s->nr_gparts;

  struct gpart *gparts = s->gparts;
  size_t *group_size = s->fof_data.group_size;
  double *group_mass = s->fof_data.group_mass;
  struct fof_CoM *group_CoM = s->fof_data.group_CoM;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  fprintf(file, "# %8s %10s %12s %12s %12s %12s\n", "Group ID", "Group Size",
          "Group Mass", "x CoM", "y CoM", "z CoM");
  fprintf(file,
          "#-------------------------------------------------------------------"
          "-------------\n");

  for (int i = 0; i < s->fof_data.num_groups; i++) {

    const size_t group_offset = group_sizes[i].index;

    if (is_local(group_offset, nr_gparts)) {

      /* Box wrap the CoM. */
      const double CoM_x =
          box_wrap(group_CoM[group_offset - node_offset].x, 0., dim[0]);
      const double CoM_y =
          box_wrap(group_CoM[group_offset - node_offset].y, 0., dim[1]);
      const double CoM_z =
          box_wrap(group_CoM[group_offset - node_offset].z, 0., dim[2]);

      fprintf(file, "  %8zu %10zu %12e %12e %12e %12e\n",
              gparts[group_offset - node_offset].group_id,
              group_size[group_offset - node_offset],
              group_mass[group_offset - node_offset], CoM_x, CoM_y, CoM_z);
    }
  }

  fclose(file);
}
