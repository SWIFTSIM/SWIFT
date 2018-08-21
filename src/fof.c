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
#include "threadpool.h"
#include "engine.h"
#include "proxy.h"
#include "common_io.h"

MPI_Datatype fof_mpi_type;
FILE *fof_file;
int node_offset;

/* Initialises parameters for the FOF search. */
void fof_init(struct space *s, long long Ngas, long long Ngparts) {

  struct engine *e = s->e;

  /* Check that we can write outputs by testing if the output
   * directory exists and is searchable and writable. */
  parser_get_param_string(e->parameter_file, "FOF:basename", s->fof_data.base_name);
  
  const char *dirp = dirname(s->fof_data.base_name);
  if (access(dirp, W_OK | X_OK) != 0) {
    error("Cannot write FOF outputs in directory %s (%s)", dirp, strerror(errno));
  }

  s->fof_data.min_group_size = parser_get_opt_param_int(e->parameter_file, "FOF:min_group_size", 20);
  const double l_x_scale = parser_get_opt_param_double(e->parameter_file, "FOF:linking_length_scale", 0.2);

  /* Calculate the particle linking length based upon the mean inter-particle
   * spacing of the DM particles. */
  const int total_nr_dmparts = Ngparts - Ngas;
  const double l_x = l_x_scale * (s->dim[0] / cbrt(total_nr_dmparts));
  s->l_x2 = l_x * l_x;

#ifdef WITH_MPI
  /* Check size of linking length against the top-level cell dimensions. */
  if(l_x > s->width[0]) error("Linking length greater than the width of a top-level cell. Need to check more than one layer of top-level cells for links.");
  
  if (MPI_Type_contiguous(sizeof(struct fof_mpi) / sizeof(unsigned char), MPI_BYTE,
        &fof_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof.");
  }
#endif

}

/* Finds the global root ID of the group a particle exists in. */
__attribute__((always_inline)) INLINE static int fof_find_global(const int i,
                                                          int *group_index) {

  int root = node_offset + i;
  if (root < node_offset) return root;
  /* TODO: May need this atomic here:
   * while (root != atomic_cas(&group_index[root], group_index[root], group_index[root])) 
   *   root = atomic_cas(&group_index[root], group_index[root], group_index[root]); */
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
__attribute__((always_inline)) INLINE static int fof_find(const int i,
                                                          int *group_index) {

  int root = i;
  /* TODO: May need this atomic here:
   * while (root != atomic_cas(&group_index[root], group_index[root], group_index[root])) 
   *   root = atomic_cas(&group_index[root], group_index[root], group_index[root]); */
  while (root != group_index[root]) {
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

/* Updates the root and checks that its value has not been changed since being read. */
__attribute__((always_inline)) INLINE static int update_root(
    volatile int* address, int y) {

  int* int_ptr = (int*)address;

  int test_val, old_val, new_val;
  old_val = *address;

  test_val = old_val;
  new_val = min(old_val, y);

  /* atomic_cas returns old_val if *int_ptr has not changed since being read.*/
  old_val = atomic_cas(int_ptr, test_val, new_val);

  if(test_val == old_val) return 1;
  else return 0;

}

__attribute__((always_inline)) INLINE static void fof_union(int *root_i, const int root_j,
                                                          int *group_index) {

  int result = 0;

  /* Loop until the root can be set to a new value. */
  do {
    int root_i_new = fof_find(*root_i, group_index);
    const int root_j_new = fof_find(root_j, group_index);

    /* Skip particles in the same group. */
    if(root_i_new == root_j_new) return;

    /* If the root ID of pj is lower than pi's root ID set pi's root to point to pj's. 
     * Otherwise set pj's to root to point to pi's.*/
    if(root_j_new < root_i_new) {
      
      /* Updates the root and checks that its value has not been changed since being read. */
      result = update_root(&group_index[root_i_new], root_j_new);
      
      /* Update root_i on the fly. */
      *root_i = root_j_new;
    }
    else {
      
      /* Updates the root and checks that its value has not been changed since being read. */
      result = update_root(&group_index[root_j_new], root_i_new);
      
      /* Update root_i on the fly. */
      *root_i = root_i_new;
    }
  } while (result != 1);
}

void fof_print_group_list(struct space *s, const size_t nr_gparts, int *group_size, int *group_index, long long *group_id, const int find_group_size, char *out_file, int (*fof_find_func)(int, int *)) {

  FILE *file;

  file = fopen(out_file, "w");

  int this_root = 0;
  for (size_t i = 0; i < nr_gparts; i++) {
    if(group_size[i] == find_group_size) {
      this_root = i;
      break;
    }
  }

  for (size_t i = 0; i < nr_gparts; i++) {
    const int root = (*fof_find_func)(i, group_index);
    if(root == this_root) {
      if(s->e->nr_nodes > 1) fprintf(file, "%7lld\n", group_id[i]);
      else fprintf(file, "%7lld\n", s->gparts[i].id_or_neg_offset);
    }
  }

  fclose(file);
}

/* Find the shortest distance between cells, remembering to account for boundary
 * conditions. */
__attribute__((always_inline)) INLINE static double cell_min_dist(
    const struct cell *ci, const struct cell *cj, const double *dim) {

  /* Get cell locations. */
  const double cix = ci->loc[0];
  const double ciy = ci->loc[1];
  const double ciz = ci->loc[2];
  const double cjx = cj->loc[0];
  const double cjy = cj->loc[1];
  const double cjz = cj->loc[2];

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  double dx[3], r2 = 0.0f;
  dx[0] = min3(fabs(nearest(cix - cjx, dim[0])),
               fabs(nearest(cix - (cjx + cj->width[0]), dim[0])),
               fabs(nearest((cix + ci->width[0]) - cjx, dim[0])));
  dx[0] = min(
      dx[0], fabs(nearest((cix + ci->width[0]) - (cjx + cj->width[0]), dim[0])));

  dx[1] = min3(fabs(nearest(ciy - cjy, dim[1])),
               fabs(nearest(ciy - (cjy + cj->width[1]), dim[1])),
               fabs(nearest((ciy + ci->width[1]) - cjy, dim[1])));
  dx[1] = min(
      dx[1], fabs(nearest((ciy + ci->width[1]) - (cjy + cj->width[1]), dim[1])));

  dx[2] = min3(fabs(nearest(ciz - cjz, dim[2])),
               fabs(nearest(ciz - (cjz + cj->width[2]), dim[2])),
               fabs(nearest((ciz + ci->width[2]) - cjz, dim[2])));
  dx[2] = min(
      dx[2], fabs(nearest((ciz + ci->width[2]) - (cjz + cj->width[2]), dim[2])));

  for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

  return r2;
}

/* Recurse on a pair of cells and perform a FOF search between cells that are within
 * range. */
static void rec_fof_search_pair(struct cell *ci, struct cell *cj, struct space *s,
                           const double *dim,
                           const double search_r2) {

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
            rec_fof_search_pair(ci->progeny[k], cj->progeny[l], s, dim, search_r2);
      }
    }
  }
  /* Perform FOF search between pairs of cells that are within the linking
   * length and not the same cell. */
  else if (ci != cj)
    fof_search_pair_cells(s, ci, cj);
  else if (ci == cj) error("Pair FOF called on same cell!!!");

}

/* Recurse on a pair of cells (one local, one foreign) and perform a FOF search between cells that are within
 * range. */
static void rec_fof_search_pair_foreign(struct cell *ci, struct cell *cj, struct space *s,
                           const double *dim,
                           const double search_r2, size_t *link_count, struct fof_mpi *part_links) {

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
            rec_fof_search_pair_foreign(ci->progeny[k], cj->progeny[l], s, dim, search_r2, link_count, part_links);
      }
    }
  }
  /* Perform FOF search between pairs of cells that are within the linking
   * length and not the same cell. */
  else if (ci != cj)
    fof_search_pair_cells_foreign(s, ci, cj, link_count, part_links);
  else if (ci == cj) error("Pair FOF called on same cell!!!");

}

/* Recurse on a pair of cells (one local, one foreign) and count the total number of possible links between them. */
static void rec_fof_search_pair_foreign_count(struct cell *ci, struct cell *cj, struct space *s,
                           const double *dim,
                           const double search_r2, size_t *nr_links) {

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
            rec_fof_search_pair_foreign_count(ci->progeny[k], cj->progeny[l], s, dim, search_r2, nr_links);
      }
    }
  }
  /* Perform FOF search between pairs of cells that are within the linking
   * length and not the same cell. */
  else if (ci != cj) {
    *nr_links += ci->gcount * cj->gcount;
    return;
  }
  else if (ci == cj) error("Pair FOF called on same cell!!!");

}

/* Recurse on a cell and perform a FOF search between cells that are within
 * range. */
static void rec_fof_search_self(struct cell *ci, struct space *s,
                           const double *dim,
                           const double search_r2) {

  /* Recurse? */
  if (ci->split) {
    
    /* Loop over all progeny. Perform pair and self recursion on progenies.*/
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        rec_fof_search_self(ci->progeny[k], s, dim, search_r2);

        for (int l = k+1; l < 8; l++)
          if (ci->progeny[l] != NULL)
            rec_fof_search_pair(ci->progeny[k], ci->progeny[l], s, dim, search_r2);
      }
    }
  }
  /* Otherwise, compute self-interaction. */
  else 
    fof_search_cell(s, ci);
}

/* Perform a FOF search on a single cell using the Union-Find algorithm.*/
void fof_search_cell(struct space *s, struct cell *c) {

  const size_t count = c->gcount;
  struct gpart *gparts = c->gparts;
  const double l_x2 = s->l_x2;
  int *group_index = s->fof_data.group_index;
  //long long *group_id = s->fof_data.group_id;

  /* Make a list of particle offsets into the global gparts array. */
  int *const offset = group_index + (ptrdiff_t)(gparts - s->gparts);

  if(c->nodeID != engine_rank) error("Performing self FOF search on foreign cell.");

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Find the root of pi. */
    int root_i = fof_find(offset[i], group_index);

    for (size_t j = i + 1; j < count; j++) {

      /* Find the root of pj. */
      const int root_j = fof_find(offset[j], group_index);
      //long long group_j = group_id[root_j];

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
      if (r2 < l_x2) {
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

/* Perform a FOF search on a pair of cells using the Union-Find algorithm.*/
void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj) {

  const size_t count_i = ci->gcount;
  const size_t count_j = cj->gcount;
  struct gpart *gparts_i = ci->gparts;
  struct gpart *gparts_j = cj->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;
  int *group_index = s->fof_data.group_index;
  //long long *group_id = s->fof_data.group_id;
  
  /* Make a list of particle offsets into the global gparts array. */
  int *const offset_i = group_index + (ptrdiff_t)(gparts_i - s->gparts);
  int *const offset_j = group_index + (ptrdiff_t)(gparts_j - s->gparts);

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  const int periodic = s->periodic;
  double diff[3];
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
    int root_i = fof_find(offset_i[i], group_index);
    
    for (size_t j = 0; j < count_j; j++) {

      /* Find the root of pj. */
      const int root_j = fof_find(offset_j[j], group_index);

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
      if (r2 < l_x2) {
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

/* Perform a FOF search between a local and foreign cell using the Union-Find algorithm. 
 * Store any links found between particles.*/
void fof_search_pair_cells_foreign(struct space *s, struct cell *ci, struct cell *cj, size_t *link_count, struct fof_mpi *part_links ) {

  const size_t count_i = ci->gcount;
  const size_t count_j = cj->gcount;
  struct gpart *gparts_i = ci->gparts;
  struct gpart *gparts_j = cj->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;
  int *group_index = s->fof_data.group_index;
  
  /* Make a list of particle offsets into the global gparts array. */
  int *const offset_i = group_index + (ptrdiff_t)(gparts_i - s->gparts);
  int *const offset_j = group_index + (ptrdiff_t)(gparts_j - s->gparts);

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Check whether cells are local to the node. */
  const int ci_local = (ci->nodeID == engine_rank);
  const int cj_local = (cj->nodeID == engine_rank);

  if((ci_local && cj_local) || (!ci_local && !cj_local)) error("FOF search of foreign cells called on two local cells or two foreign cells.");

  /* Get the relative distance between the pairs, wrapping. */
  const int periodic = s->periodic;
  double diff[3];
  
  if(ci_local) {

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
      const int root_i = fof_find_global(offset_i[i] - node_offset, group_index);

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
          
          /* Store the particle group IDs for communication. */
          part_links[*link_count].group_i = root_i;        
          part_links[*link_count].group_j = pj->root;        
          (*link_count)++;
        
        }
      }
    }

  }

  if(cj_local) {

    for (int k = 0; k < 3; k++) {
      diff[k] = ci->loc[k] - cj->loc[k];
      if (periodic && diff[k] < -dim[k] / 2)
        shift[k] = dim[k];
      else if (periodic && diff[k] > dim[k] / 2)
        shift[k] = -dim[k];
      else
        shift[k] = 0.0;
      diff[k] += shift[k];
    }

    /* Loop over particles and find which particles belong in the same group. */
    for (size_t j = 0; j < count_j; j++) {

      struct gpart *pj = &gparts_j[j];
      const double pjx = pj->x[0] - shift[0];
      const double pjy = pj->x[1] - shift[1];
      const double pjz = pj->x[2] - shift[2];

      /* Find the root of pj. */
      int root_j = fof_find_global(offset_j[j] - node_offset, group_index);

      for (size_t i = 0; i < count_i; i++) {

        struct gpart *pi = &gparts_i[i];
        const double pix = pi->x[0];
        const double piy = pi->x[1];
        const double piz = pi->x[2];

        /* Compute pairwise distance, remembering to account for boundary
         * conditions. */
        float dx[3], r2 = 0.0f;
        dx[0] = pjx - pix;
        dx[1] = pjy - piy;
        dx[2] = pjz - piz;

        for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

        /* Hit or miss? */
        if (r2 < l_x2) {

          /* Store the particle group IDs for communication. */
          part_links[*link_count].group_i = root_j;        
          part_links[*link_count].group_j = pi->root; 
          (*link_count)++;
          
        }
      }
    }
  }

}

/* Perform a FOF search on gravity particles using the cells and applying the
 * Union-Find algorithm.*/
void fof_search_tree_serial(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_cells = s->nr_cells;
  struct gpart *gparts = s->gparts;
  long long *group_id;
  int *group_index;
  int *group_size;
  double *group_mass;
  int num_groups = 0;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;
  ticks tic = getticks();

  message("Searching %zu gravity particles for links with l_x2: %lf", nr_gparts,
          s->l_x2);

  /* Allocate and initialise array of particle group IDs. */
  if(s->fof_data.group_index != NULL) free(s->fof_data.group_index);
  if(s->fof_data.group_id != NULL) free(s->fof_data.group_id);

  if (posix_memalign((void **)&s->fof_data.group_index, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group indices for FOF search.");

  if (posix_memalign((void **)&s->fof_data.group_id, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");
  
  /* Initial group ID is particle offset into array. */
  for (size_t i = 0; i < nr_gparts; i++) s->fof_data.group_index[i] = i;
  for (size_t i = 0; i < nr_gparts; i++) s->fof_data.group_id[i] = gparts[i].id_or_neg_offset;
  
  group_index = s->fof_data.group_index;
  group_id = s->fof_data.group_id;
  
  message("Rank: %d, Allocated group_index array of size %zu", engine_rank, s->nr_gparts);

  /* Allocate and initialise a group size array. */
  if (posix_memalign((void **)&group_size, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of group size for FOF search.");
  
  /* Allocate and initialise a group mass array. */
  if (posix_memalign((void **)&group_mass, 32, nr_gparts * sizeof(double)) != 0)
    error("Failed to allocate list of group masses for FOF search.");

  bzero(group_size, nr_gparts * sizeof(int));
  bzero(group_mass, nr_gparts * sizeof(double));

  /* Loop over cells and find which cells are in range of each other to perform
   * the FOF search. */
  for (size_t cid = 0; cid < nr_cells; cid++) {

    struct cell *restrict ci = &s->cells_top[cid];
    
    /* Only perform FOF search on local cells. */
    if(ci->nodeID == engine_rank) {

      /* Skip empty cells. */
      if(ci->gcount == 0) continue;

      /* Perform FOF search on local particles within the cell. */
      rec_fof_search_self(ci, s, dim, search_r2);

      /* Loop over all top-level cells skipping over the cells already searched.
      */
      for (size_t cjd = cid + 1; cjd < nr_cells; cjd++) {

        struct cell *restrict cj = &s->cells_top[cjd];

        /* Only perform FOF search on local cells. */
        if(cj->nodeID == engine_rank) {
          /* Skip empty cells. */
          if(cj->gcount == 0) continue;

          rec_fof_search_pair(ci, cj, s, dim, search_r2);
        }
      }
    }   
  }

  /* Calculate the total number of particles in each group, group mass and the total number of groups. */
  for (size_t i = 0; i < nr_gparts; i++) {
    const int root = fof_find(i, group_index);
    group_size[root]++;
    group_mass[root] += gparts[i].mass;
    if(group_index[i] == i) num_groups++;
    group_id[i] = group_id[root - node_offset];
  }

  fof_dump_group_data("fof_output_tree_serial.dat", nr_gparts, group_index,
                      group_size, group_id, group_mass, 1);

  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_index = 0, max_group_mass_id = 0;
  double max_group_mass = 0;
  for (size_t i = 0; i < nr_gparts; i++) {

    /* Find the total number of particles in groups. */
    if (group_size[i] > 1) num_parts_in_groups += group_size[i];

    /* Find the largest group. */
    if (group_size[i] > max_group_size) {
      max_group_size = group_size[i];
      max_group_index = i;
    }
    
    /* Find the largest group by mass. */
    if (group_mass[i] > max_group_mass) {
      max_group_mass = group_mass[i];
      max_group_mass_id = i;
    }
  }

  message(
      "No. of groups: %d. No. of particles in groups: %d. No. of particles not "
      "in groups: %zu.",
      num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_index);
  message("Biggest group by mass: %e with ID: %d", max_group_mass, max_group_mass_id);

  free(group_size);
  free(group_mass);
  
  message("Serial FOF search took: %.3f %s.",
      clocks_from_ticks(getticks() - tic), clocks_getunit());
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
  struct cell *cells = (struct cell *)map_data;

  const size_t nr_cells = s->nr_cells;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;
 
  /* Make a list of cell offsets into the top-level cell array. */
  int *const offset = s->cell_index + (ptrdiff_t)(cells - s->cells_top);

  /* Loop over cells and find which cells are in range of each other to perform
   * the FOF search. */
  for (size_t ind = 0; ind < num_elements; ind++) {
    
    /* Get the cell. */
    struct cell *restrict ci = &cells[ind];
    
    /* Only perform FOF search on local cells. */
    if(ci->nodeID == engine_rank) {

      /* Skip empty cells. */
      if(ci->gcount == 0) continue;

      /* Perform FOF search on local particles within the cell. */
      rec_fof_search_self(ci, s, dim, search_r2);

      /* Loop over all top-level cells skipping over the cells already searched.
      */
      for (size_t cjd = offset[ind] + 1; cjd < nr_cells; cjd++) {

        struct cell *restrict cj = &s->cells_top[cjd];

        /* Only perform FOF search on local cells. */
        if(cj->nodeID == engine_rank) {
          
          /* Skip empty cells. */
          if(cj->gcount == 0) continue;

          rec_fof_search_pair(ci, cj, s, dim, search_r2);
        }
      }
    }
  }
}

/**
 * @brief Search foreign cells for links and communicate any found to the appropriate node.
 *
 * @param s Pointer to a #space.
 */
void fof_search_foreign_cells(struct space *s) {

#ifdef WITH_MPI

  struct engine *e = s->e;
  struct cell *cells = s->cells_top;
  int *group_index = s->fof_data.group_index;
  int *group_size = s->fof_data.group_size;
  double *group_mass = s->fof_data.group_mass;
  const int nr_gparts = s->nr_gparts;
  const size_t nr_cells = s->nr_cells;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;

  message("Searching foreign cells for links.");

  size_t nr_links = 0;
  int count = 0;

  /* Make group IDs globally unique. */  
  for (size_t i = 0; i < nr_gparts; i++) group_index[i] += node_offset;

  /* Loop over top-level cells and find local cells that touch foreign cells. Calculate the total number of possible links between these cells. */
  for (size_t cid = 0; cid < nr_cells; cid++) {

    struct cell *restrict ci = &cells[cid];

    /* Skip empty cells. */
    if(ci->gcount == 0) continue;

    /* Loop over all top-level cells skipping over the cells already searched. */
    for (size_t cjd = cid + 1; cjd < nr_cells; cjd++) {

      struct cell *restrict cj = &cells[cjd];

      /* Only perform pair FOF search between a local and foreign cell. */
      if((ci->nodeID == engine_rank && cj->nodeID != engine_rank) || (ci->nodeID != engine_rank && cj->nodeID == engine_rank)) {
        /* Skip empty cells. */
        if(cj->gcount == 0) continue;

        const double r2 = cell_min_dist(ci, cj, dim);
        /* Assume there are only links between cells that touch. */
        if(r2 < search_r2) {
          rec_fof_search_pair_foreign_count(ci, cj, s, dim, search_r2, &nr_links);
          count++;
        }
      }
    }
  }

  message("Rank: %d, Total no. of possible links: %zu, cells touching: %d", engine_rank, nr_links, count);

  struct fof_mpi *part_links;
  struct cell **interface_cells;
  int interface_cell_count = 0;
  
  if (posix_memalign((void**)&part_links, SWIFT_STRUCT_ALIGNMENT,
                       nr_links * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for FOF links over an MPI domain");

  if (posix_memalign((void**)&interface_cells, SWIFT_STRUCT_ALIGNMENT,
                       count * sizeof(struct cell *)) != 0)
    error("Error while allocating memory for FOF interface cells");
  
  /* Loop over cells_in and cells_out for each proxy and find which cells are in range of each other to perform
   * the FOF search. Store local cells that are touching foreign cells in a list. */
  for(int i=0; i<e->nr_proxies; i++) {
  
    message("Rank %d, proxy: %d has %d cells_out and %d cells_in.", engine_rank, i, e->proxies[i].nr_cells_out, e->proxies[i].nr_cells_in);

    for(int j=0; j<e->proxies[i].nr_cells_out; j++) {

      struct cell *restrict local_cell = e->proxies[i].cells_out[j];
      int found = 0;

      /* Skip empty cells. */
      if(local_cell->gcount == 0) continue;

      for(int k=0; k<e->proxies[i].nr_cells_in; k++) {
      
        struct cell *restrict foreign_cell = e->proxies[i].cells_in[k];
      
        /* Skip empty cells. */
        if(foreign_cell->gcount == 0) continue;
        
        /* Check if local cell has already been added to the local list of cells. */
        if(!found) {
          const double r2 = cell_min_dist(local_cell, foreign_cell, dim);
          if(r2 < search_r2) {        
            interface_cells[interface_cell_count++] = local_cell;
            found = 1;
          }
        }
      }
    }
  }

  message("No. of interface cells: %d", interface_cell_count);

  /* Set the root of outgoing particles. */
  for(int i=0; i<e->nr_proxies; i++) {
  
    for(int j=0; j<e->proxies[i].nr_cells_out; j++) {

      struct cell *restrict local_cell = e->proxies[i].cells_out[j];
      struct gpart *gparts = local_cell->gparts;
  
      /* Make a list of particle offsets into the global gparts array. */
      int *const offset = group_index + (ptrdiff_t)(gparts - s->gparts);

      /* Set each particle's root found in the local FOF.*/
      for(int k=0; k<local_cell->gcount; k++) {
        const int root = fof_find_global(offset[k] - node_offset, group_index);  
        gparts[k].root = root;
      }
    }
  }

  struct scheduler *sched = &e->sched;
  struct task *tasks = sched->tasks;
  int nr_sends = 0, nr_recvs = 0;

  /* Activate the send and receive tasks for the gparts. */
  for(int i=0; i<sched->nr_tasks; i++) {
  
    struct task *t = &tasks[i];

    if(t->type == task_type_send && t->subtype == task_subtype_gpart) {
      scheduler_activate(sched,t);
      nr_sends++;
    }

    if(t->type == task_type_recv && t->subtype == task_subtype_gpart) {
      scheduler_activate(sched,t);
      nr_recvs++;
    } 
  }

  message("Rank %d. No. of sends activated: %d, no. of receives activated: %d", engine_rank, nr_sends, nr_recvs);

  /* Perform send and receive tasks. */
  engine_launch(e);

  size_t part_link_count = 0;

  /* Loop over each interface cell and find all particle links with foreign cells. */
  for(int i=0; i<interface_cell_count; i++) {

    struct cell *restrict local_cell = interface_cells[i];
      
    /* Skip empty cells. */
    if(local_cell->gcount == 0) continue;

    for(int j=0; j<e->nr_proxies; j++) {

      for(int k=0; k<e->proxies[j].nr_cells_in; k++) {

        struct cell *restrict foreign_cell = e->proxies[j].cells_in[k];

        /* Skip empty cells. */
        if(foreign_cell->gcount == 0) continue;

        rec_fof_search_pair_foreign(local_cell, foreign_cell, s, dim, search_r2, &part_link_count, part_links);

      }
    }
  }

  message("Rank %d found %zu links between local and foreign groups.", engine_rank, part_link_count);

  struct fof_mpi *group_links;
  int group_link_count = 0;
  
  if (posix_memalign((void**)&group_links, SWIFT_STRUCT_ALIGNMENT,
                       part_link_count * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for unique set of group links across MPI domains");

  /* Compress the particle links to a unique set of links between groups over an MPI domain. */
  for(int i=0; i<part_link_count; i++) {
    
    int found = 0;
    int local_root = part_links[i].group_i;
    int foreign_root = part_links[i].group_j;

    if(local_root == foreign_root) error("Particles have same root. local_root: %d, foreign_root: %d", local_root, foreign_root);

    /* Loop over list and check that the link doesn't already exist. */
    for(int j=0; j<group_link_count; j++) {
      if(group_links[j].group_i == local_root && group_links[j].group_j == foreign_root) {
        found = 1;
        break;
      }
    }

    /* If it doesn't already exist in the list add it. */
    if(!found) {
      group_links[group_link_count].group_i = local_root;
      group_links[group_link_count].group_i_size = group_size[local_root - node_offset];
      group_links[group_link_count].group_i_mass = group_mass[local_root - node_offset];
      group_links[group_link_count++].group_j = foreign_root;
    }

  }

  message("Rank %d found %d unique group links.", engine_rank, group_link_count);
 
  struct fof_mpi *global_group_links = NULL;
  int *displ = NULL, *group_link_counts = NULL;
  int global_group_link_count = 0;

  /* Sum the total number of links across MPI domains over each MPI rank. */
  MPI_Allreduce(&group_link_count, &global_group_link_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
  message("Global list count: %d", global_group_link_count);

  /* Unique set of links is half of all group links as each link is found twice by opposing MPI ranks. */
  if (posix_memalign((void**)&global_group_links, SWIFT_STRUCT_ALIGNMENT,
                       global_group_link_count * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for the global list of group links");
  
  if (posix_memalign((void**)&group_link_counts, SWIFT_STRUCT_ALIGNMENT,
                       e->nr_nodes * sizeof(int)) != 0)
    error("Error while allocating memory for the number of group links on each MPI rank");

  if (posix_memalign((void**)&displ, SWIFT_STRUCT_ALIGNMENT,
                       e->nr_nodes * sizeof(int)) != 0)
    error("Error while allocating memory for the displacement in memory for the global group link list");
  
  /* Gather the total number of links on each rank. */
  MPI_Allgather(&group_link_count, 1, MPI_INT, group_link_counts, 1, MPI_INT, MPI_COMM_WORLD);

  /* Set the displacements into the global link list using the link counts from each rank */
  displ[0] = 0;
  for(int i=1; i<e->nr_nodes; i++) displ[i] = displ[i-1] + group_link_counts[i-1];

  /* Gather the global link list on all ranks. */
  MPI_Allgatherv(group_links, group_link_count, fof_mpi_type, global_group_links, group_link_counts, displ, fof_mpi_type, MPI_COMM_WORLD);

  /* Transform the group IDs to a local list going from 0-group_count so a union-find can be performed. */
  int *global_group_index = NULL, *global_group_id = NULL, *global_group_size = NULL;
  double *global_group_mass = NULL;
  int group_count = 0;

  if (posix_memalign((void**)&global_group_index, SWIFT_STRUCT_ALIGNMENT,
                      global_group_link_count  * sizeof(int)) != 0)
    error("Error while allocating memory for the displacement in memory for the global group link list");
  
  if (posix_memalign((void**)&global_group_id, SWIFT_STRUCT_ALIGNMENT,
                      global_group_link_count  * sizeof(int)) != 0)
    error("Error while allocating memory for the displacement in memory for the global group link list");
  
  if (posix_memalign((void**)&global_group_size, SWIFT_STRUCT_ALIGNMENT,
                      global_group_link_count  * sizeof(int)) != 0)
    error("Error while allocating memory for the displacement in memory for the global group link list");
  
  if (posix_memalign((void**)&global_group_mass, SWIFT_STRUCT_ALIGNMENT,
                      global_group_link_count  * sizeof(double)) != 0)
    error("Error while allocating memory for the displacement in memory for the global group link list");

  bzero(global_group_size, global_group_link_count * sizeof(int));
  bzero(global_group_mass, global_group_link_count * sizeof(double));

  /* Compress the list of group links across an MPI domain by removing the symmetric cases. */
  /* Store each group ID once along with its size. */
  for(int i=0; i<global_group_link_count; i++) {
    
    int found_i = 0, found_j = 0;
    int group_i = global_group_links[i].group_i;
    int group_j = global_group_links[i].group_j;

    /* Loop over list and check that the link doesn't already exist. */
    for(int j=0; j<group_count; j++) {
      if(global_group_id[j] == group_i) {
        found_i = 1;
        break;
      }
    }
    
    for(int j=0; j<group_count; j++) {
      if(global_group_id[j] == group_j) {
        found_j = 1;
        break;
      }
    }

    /* If it doesn't already exist in the list add it. */
    if(!found_i) {
      global_group_size[group_count] += global_group_links[i].group_i_size;
      global_group_mass[group_count] += global_group_links[i].group_i_mass;
      global_group_id[group_count++] = group_i;
    }
    
    if(!found_j) {
      /* Need to search for group_j sizes as the group size is only known on the local node. */
      for(int j=0; j<global_group_link_count; j++) {
        if(global_group_links[j].group_i == group_j) {
          global_group_size[group_count] += global_group_links[j].group_i_size;
          global_group_mass[group_count] += global_group_links[j].group_i_mass;
          break;
        }
      }

      global_group_id[group_count++] = group_j;
    }

  }

  /* Create a global_group_index list of groups across MPI domains so that you can perform a union-find locally on each node. */
  /* The value of which is an offset into global_group_id, which is the actual root. */
  for(int i=0; i<group_count; i++) global_group_index[i] = i;

  /* Save group links for each rank to a file. */
  //char fof_map_filename[200] = "group_map";
  //sprintf(fof_map_filename + strlen(fof_map_filename), "_%d.dat",
  //    engine_rank);
  //
  //fof_file = fopen(fof_map_filename, "w");
  //fprintf(fof_file, "# %7s %7s %7s %7s\n", "Index", "Group ID", "Group Size", "Group Mass");
  //fprintf(fof_file, "#-------------------------------\n");

  //for(int i=0; i<group_count; i++) {
  //  fprintf(fof_file, "  %7d %7d %7d %7e\n", global_group_index[i], global_group_id[i], global_group_size[i], global_group_mass[i]);
  //}
  //
  //fclose(fof_file);

  /* Perform a union-find on the group links. */
  for(int i=0; i<global_group_link_count; i++) {
  
    int find_i = 0, find_j = 0;

    /* Find group offset into global_group_id */
    for(int j=0; j<group_count; j++) {  
      if(global_group_id[j] == global_group_links[i].group_i) {
        find_i = j;
        break;
      }
    }
    
    for(int j=0; j<group_count; j++) {  
      if(global_group_id[j] == global_group_links[i].group_j) {
        find_j = j;
        break;
      }
    }

    /* Use the offset to find the group's root. */
    int root_i = fof_find(find_i, global_group_index);
    int root_j = fof_find(find_j, global_group_index);
    
    int group_i = global_group_id[root_i];
    int group_j = global_group_id[root_j];

    /* Update roots accordingly. */
    if(group_j < group_i) global_group_index[root_i] = root_j;
    else global_group_index[root_j] = root_i;

  }

  /* Update each group locally with new root information. */
  for(int i=0; i<group_count; i++) {
    
    int group_id = global_group_id[i];
    int new_root = global_group_id[fof_find(global_group_index[i], global_group_index)];
    
    /* If the group is local update its root and size. */
    if(group_id >= node_offset && group_id < node_offset + nr_gparts &&
       new_root != group_id) {
        group_index[group_id - node_offset] = new_root;
        group_size[group_id - node_offset] -= global_group_size[i];
        group_mass[group_id - node_offset] -= global_group_mass[i];
    }
     
    /* If the group linked to a local root update its size. */
    if(new_root >= node_offset && new_root < node_offset + nr_gparts &&
       new_root != group_id) {
      group_size[new_root - node_offset] += global_group_size[i];
      group_mass[new_root - node_offset] += global_group_mass[i];
    }

  }

  /* Clean up memory. */
  free(part_links);
  free(interface_cells);
  free(group_links);
  free(global_group_links);
  free(displ);
  free(global_group_index);
  free(global_group_size);
  free(global_group_mass);
  free(global_group_id);

  message("Rank %d finished linking local roots to foreign roots.", engine_rank);

#endif /* WITH_MPI */

}

/* Perform a FOF search on gravity particles using the cells and applying the
 * Union-Find algorithm.*/
void fof_search_tree(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_cells = s->nr_cells;
  const int min_group_size = s->fof_data.min_group_size;
  struct gpart *gparts = s->gparts;
  long long *group_id;
  int *group_index, *group_size;
  double *group_mass;
  int num_groups = 0, num_parts_in_groups = 0, max_group_size = 0;
  double max_group_mass = 0;
  ticks tic = getticks();

  char output_file_name[FILENAME_BUFFER_SIZE];
  snprintf(output_file_name, FILENAME_BUFFER_SIZE, s->fof_data.base_name);

  message("Searching %zu gravity particles for links with l_x2: %lf", nr_gparts,
          s->l_x2);

  node_offset = 0;

#ifdef WITH_MPI
  if (s->e->nr_nodes > 1) {
  
    int *global_nr_gparts;

    /* Find the total no. of particles on each node and calculate your unique offset. */
    if (posix_memalign((void**)&global_nr_gparts, SWIFT_STRUCT_ALIGNMENT,
          s->e->nr_nodes * sizeof(int)) != 0)
      error("Error while allocating memory for global_nr_gparts array.");

    bzero(global_nr_gparts, s->e->nr_nodes * sizeof(int));

    MPI_Allgather(&s->nr_gparts, 1, MPI_INT, global_nr_gparts, 1, MPI_INT, MPI_COMM_WORLD);

    if(engine_rank == 0) {
      printf("No. of particles: [%d", global_nr_gparts[0]);
      for(int i = 1; i<s->e->nr_nodes; i++) printf(", %d", global_nr_gparts[i]);
      printf("]\n");
    }

    for(int i = 0; i<engine_rank; i++) node_offset += global_nr_gparts[i];
    
    /* Clean up memory. */
    free(global_nr_gparts);
  }
#endif

  message("Node offset: %d", node_offset);

  /* Allocate and initialise array of particle group IDs. */
  if(s->fof_data.group_index != NULL) free(s->fof_data.group_index);
  if(s->fof_data.group_id != NULL) free(s->fof_data.group_id);

  /* Allocate and initialise a group index array. */
  if (posix_memalign((void **)&s->fof_data.group_index, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group indices for FOF search.");
  
  /* Allocate and initialise a group size array. */
  if (posix_memalign((void **)&s->fof_data.group_size, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of group size for FOF search.");
 
  /* Allocate and initialise a group ID array. */
  if (posix_memalign((void **)&s->fof_data.group_id, 32, nr_gparts * sizeof(long long)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");

  /* Allocate and initialise a group mass array. */
  if (posix_memalign((void **)&s->fof_data.group_mass, 32, nr_gparts * sizeof(double)) != 0)
    error("Failed to allocate list of group masses for FOF search.");

  /* Initial group ID is particle offset into array. */
  for (size_t i = 0; i < nr_gparts; i++) s->fof_data.group_index[i] = i;
  for (size_t i = 0; i < nr_gparts; i++) s->fof_data.group_id[i] = gparts[i].id_or_neg_offset;

  group_index = s->fof_data.group_index;
  group_size = s->fof_data.group_size;
  group_id = s->fof_data.group_id;
  group_mass = s->fof_data.group_mass;
  
  message("Rank: %d, Allocated group_index array of size %zu", engine_rank, s->nr_gparts);

  bzero(group_size, nr_gparts * sizeof(int));
  bzero(group_mass, nr_gparts * sizeof(double));

  /* Activate all the regular tasks */
  threadpool_map(&s->e->threadpool, fof_search_tree_mapper, s->cells_top,
                 nr_cells, sizeof(struct cell), 1, s);

  /* Calculate the total number of particles in each group, group mass and group ID. */
  for (size_t i = 0; i < nr_gparts; i++) {
    int root = fof_find(i, group_index);
    group_size[root]++;
    group_mass[root] += gparts[i].mass;
    group_id[i] = group_id[root];
  }

#ifdef WITH_MPI
  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE, "_mpi_rank_%d.dat", engine_rank);
  
  if (s->e->nr_nodes > 1) {
    /* Search for group links across MPI domains. */
    fof_search_foreign_cells(s);
  }  
#else
  snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE, ".dat");
#endif
 
  /* Dump group data. */ 
  fof_dump_group_data(output_file_name, nr_gparts, group_index,
      group_size, group_id, group_mass, min_group_size);

  int num_groups_local = 0, num_parts_in_groups_local = 0, max_group_size_local = 0;
  double max_group_mass_local = 0;
  
  for (size_t i = 0; i < nr_gparts; i++) {
    
    /* Find the total number of groups. */
    if(group_index[i] == i + node_offset && group_size[i] >= min_group_size) num_groups_local++;

    /* Find the total number of particles in groups. */
    if (group_size[i] >= min_group_size) num_parts_in_groups_local += group_size[i];

    /* Find the largest group. */
    if (group_size[i] > max_group_size_local) max_group_size_local = group_size[i];
    
    /* Find the largest group by mass. */
    if (group_mass[i] > max_group_mass_local) max_group_mass_local = group_mass[i];
  }

  /* Find global properties. */
#ifdef WITH_MPI
  MPI_Reduce(&num_groups_local, &num_groups, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_parts_in_groups_local, &num_parts_in_groups, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_group_size_local, &max_group_size, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_group_mass_local, &max_group_mass, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else
  num_groups += num_groups_local;
  num_parts_in_groups += num_parts_in_groups_local;
  max_group_size = max_group_size_local;
  max_group_mass = max_group_mass_local;
#endif /* WITH_MPI */
  
  if(engine_rank == 0) {
    message(
        "No. of groups: %d. No. of particles in groups: %d. No. of particles not "
        "in groups: %zu.",
        num_groups, num_parts_in_groups, s->e->total_nr_gparts - num_parts_in_groups);

    message("Largest group by size: %d", max_group_size);
    message("Largest group by mass: %e", max_group_mass);
  
    message("FOF search took: %.3f %s.",
        clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

}

/* Dump FOF group data. */
void fof_dump_group_data(char *out_file, const size_t nr_gparts, int *group_index,
                         int *group_size, long long *group_id, double *group_mass, const int min_group_size) {

  FILE *file = fopen(out_file, "w");
  fprintf(file, "# %7s %7s %7s %7s %7s\n", "ID", "Root ID", "Group Size", "Group Mass", "Group ID");
  fprintf(file, "#---------------------------------------\n");

  for (size_t i = 0; i < nr_gparts; i++) {
    if(group_size[i] >= min_group_size) {
      const int root = fof_find_global(i - node_offset, group_index);
      fprintf(file, "  %7zu %7d %7d %7e    %10lld\n", i, root, group_size[i], group_mass[i], group_id[i]);
    }
  }

  fclose(file);
}
