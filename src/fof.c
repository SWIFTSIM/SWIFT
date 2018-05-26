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

/* This object's header. */
#include "fof.h"

/* Local headers. */
//#include "active.h"

/* Finds the root ID of the group a particle exists in. */
__attribute__((always_inline)) INLINE static int fof_find(const int i,
                                                          int *pid) {

  int root = i;
  while (root != pid[root]) root = pid[root];

  /* Perform path compression. */
  // int index = i;
  // while(index != root) {
  //  int next = pid[index];
  //  pid[index] = root;
  //  index = next;
  //}

  return root;
}

/* Find the shortest distance between cells, remembering to account for boundary
 * conditions. */
__attribute__((always_inline)) INLINE static double cell_min_dist(
    const struct cell *ci, const struct cell *cj, const double cix,
    const double ciy, const double ciz, const double *dim) {

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

/* Recurse on a cell and perform a FOF search between cells that are within
 * range. */
static void rec_fof_search_sub(struct cell *ci, struct cell *cj,
                               struct space *s, int *pid, int *num_groups,
                               const double *dim, const double search_r2) {

  /* Recurse on cj. */
  if (cj->split)
    for (int k = 0; k < 8; k++)
      if (cj->progeny[k] != NULL)
        rec_fof_search_sub(ci, cj->progeny[k], s, pid, num_groups, dim,
                           search_r2);

  /* No progeny? */
  if (!cj->split) {
    const double cix = ci->loc[0];
    const double ciy = ci->loc[1];
    const double ciz = ci->loc[2];

    /* Find the shortest distance between cells, remembering to account for
     * boundary conditions. */
    const double r2 = cell_min_dist(ci, cj, cix, ciy, ciz, dim);

    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    if (r2 < search_r2 && ci != cj) {
      fof_search_pair_cells(s, ci, cj, pid, num_groups);
    }
  }
}

/* Recurse on a top-level cell and perform a search of all other top-level
 * cells. Recurse on top-level cells that are within range and perform a FOF
 * search between cells. */
static void rec_fof_search(struct cell *ci, const int cid, struct space *s,
                           int *pid, int *num_groups, const double *dim,
                           const double search_r2) {

  /* Recurse on ci. */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        rec_fof_search(ci->progeny[k], cid, s, pid, num_groups, dim, search_r2);

  /* No progeny? */
  if (!ci->split) {
    const double cix = ci->loc[0];
    const double ciy = ci->loc[1];
    const double ciz = ci->loc[2];

    /* Perform FOF search on local particles within the cell. */
    fof_search_cell(s, ci, pid, num_groups);

    /* Loop over all top-level cells skipping over the cells already searched.
     */
    for (int cjd = cid; cjd < s->nr_cells; cjd++) {

      struct cell *restrict cj = &s->cells_top[cjd];

      /* Find the shortest distance between cells, remembering to account for
       * boundary conditions. */
      const double r2 = cell_min_dist(ci, cj, cix, ciy, ciz, dim);

      /* If the leaf cell of ci is within the linking length of the top-level
       * cell cj recurse on cj. Otherwise skip over cj as all of its leaf cells
       * will also be out of range. */
      if (r2 > search_r2)
        continue;
      else if (cj->split)
        rec_fof_search_sub(ci, cj, s, pid, num_groups, dim, search_r2);
      /* Perform FOF search between pairs of cells that are within the linking
       * length and not the same cell. */
      else if (ci != cj)
        fof_search_pair_cells(s, ci, cj, pid, num_groups);
    }
  }
}

static void rec_fof_search_ci_cj(struct cell *ci, struct cell *cj, struct space *s,
                           int *pid, int *num_groups, const double *dim,
                           const double search_r2) {

  const double cix = ci->loc[0];
  const double ciy = ci->loc[1];
  const double ciz = ci->loc[2];

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, cix, ciy, ciz, dim);

  if (r2 > search_r2)
    return;
  /* Recurse on both cells if they are both split. */
  else if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_ci_cj(ci->progeny[k], cj->progeny[l], s, pid, num_groups, dim, search_r2);
      }
    }
  }
  /* Perform FOF search between pairs of cells that are within the linking
   * length and not the same cell. */
  else if (ci != cj)
    fof_search_pair_cells(s, ci, cj, pid, num_groups);

}

static void rec_fof_search_self_ci(struct cell *ci, struct space *s,
                           int *pid, int *num_groups, const double *dim,
                           const double search_r2) {

  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        rec_fof_search_self_ci(ci->progeny[k], s, pid, num_groups, dim, search_r2);

        for (int l = k+1; l < 8; l++)
          if (ci->progeny[l] != NULL)
            rec_fof_search_ci_cj(ci->progeny[k], ci->progeny[l], s, pid, num_groups, dim, search_r2);
      }
    }
  }
  else 
    fof_search_cell(s, ci, pid, num_groups);
}

/* Perform naive N^2 FOF search on gravity particles using the Union-Find
 * algorithm.*/
void fof_search_serial(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;
  int *pid;
  int *group_sizes;
  int num_groups = nr_gparts;

  message("Searching %ld gravity particles for links with l_x2: %lf", nr_gparts,
          l_x2);

  /* Allocate and populate array of particle group IDs. */
  if (posix_memalign((void **)&pid, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");

  for (size_t i = 0; i < nr_gparts; i++) pid[i] = i;

  /* Allocate and initialise a group size array. */
  if (posix_memalign((void **)&group_sizes, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of number in groups for FOF search.");

  bzero(group_sizes, nr_gparts * sizeof(int));

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < nr_gparts; i++) {

    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    for (size_t j = i + 1; j < nr_gparts; j++) {

      /* Find the roots of pi and pj. */
      const int root_i = fof_find(i, pid);
      const int root_j = fof_find(j, pid);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      struct gpart *pj = &gparts[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance, remembering to account for boundary
       * conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) {
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < l_x2) {

        if (root_j < root_i)
          pid[root_i] = root_j;
        else
          pid[root_j] = root_i;

        num_groups--;
      }
    }
  }

  /* Calculate the total number of particles in each group. */
  for (size_t i = 0; i < nr_gparts; i++) group_sizes[fof_find(i, pid)]++;

  fof_dump_group_data("fof_output_serial.dat", nr_gparts, pid, group_sizes);

  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_id = 0;
  for (size_t i = 0; i < nr_gparts; i++) {

    if (group_sizes[i] > 1) num_parts_in_groups += group_sizes[i];
    if (group_sizes[i] > max_group_size) {
      max_group_size = group_sizes[i];
      max_group_id = i;
    }
  }

  message(
      "No. of groups: %d. No. of particles in groups: %d. No. of particles not "
      "in groups: %ld.",
      num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_id);

  free(pid);
  free(group_sizes);
}

/* Perform a FOF search on a single cell using the Union-Find algorithm.*/
void fof_search_cell(struct space *s, struct cell *c, int *pid,
                     int *num_groups) {

  const size_t count = c->gcount;
  struct gpart *gparts = c->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];
    const size_t offset_i = pi->offset;

    for (size_t j = i + 1; j < count; j++) {

      struct gpart *pj = &gparts[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const size_t offset_j = pj->offset;

      /* Compute pairwise distance, remembering to account for boundary
       * conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) {
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Find the roots of pi and pj. */
      const int root_i = fof_find(offset_i, pid);
      const int root_j = fof_find(offset_j, pid);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Hit or miss? */
      if (r2 < l_x2) {

        if (root_j < root_i)
          pid[root_i] = root_j;
        else
          pid[root_j] = root_i;

        (*num_groups)--;
      }
    }
  }
}

/* Perform a FOF search on a pair of cells using the Union-Find algorithm.*/
void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj,
                           int *pid, int *num_groups) {

  const size_t count_i = ci->gcount;
  const size_t count_j = cj->gcount;
  struct gpart *gparts_i = ci->gparts;
  struct gpart *gparts_j = cj->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;

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
    const size_t offset_i = pi->offset;

    for (size_t j = 0; j < count_j; j++) {

      struct gpart *pj = &gparts_j[j];
      const size_t offset_j = pj->offset;

      /* Find the roots of pi and pj. */
      const int root_i = fof_find(offset_i, pid);
      const int root_j = fof_find(offset_j, pid);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Compute pairwise distance, remembering to account for boundary
       * conditions. */
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        if (root_j < root_i)
          pid[root_i] = root_j;
        else
          pid[root_j] = root_i;

        (*num_groups)--;
      }
    }
  }
}

/* Perform a FOF search on gravity particles using the cells and applying the
 * Union-Find algorithm.*/
void fof_search_tree_serial(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_cells = s->nr_cells;
  int *pid, *group_sizes;
  int num_groups = nr_gparts;
  struct gpart *gparts = s->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->l_x2;

  message("Searching %ld gravity particles for links with l_x2: %lf", nr_gparts,
          s->l_x2);

  /* Allocate and populate array of particle group IDs. */
  if (posix_memalign((void **)&pid, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");

  for (size_t i = 0; i < nr_gparts; i++) {
    gparts[i].offset = i;
    pid[i] = i;
  }

  /* Allocate and initialise a group size array. */
  if (posix_memalign((void **)&group_sizes, 32, nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of group size for FOF search.");

  bzero(group_sizes, nr_gparts * sizeof(int));

  /* Loop over cells and find which cells are in range of each other to perform
   * the FOF search. */
  for (size_t cid = 0; cid < nr_cells; cid++) {

    struct cell *restrict ci = &s->cells_top[cid];

    /* Perform FOF search on local particles within the cell. */
    rec_fof_search_self_ci(ci, s, pid, &num_groups, dim, search_r2);

    /* Loop over all top-level cells skipping over the cells already searched.
    */
    for (int cjd = cid + 1; cjd < nr_cells; cjd++) {

      struct cell *restrict cj = &s->cells_top[cjd];

      rec_fof_search_ci_cj(ci, cj, s, pid, &num_groups, dim, search_r2);
    }   
  }

  /* Calculate the total number of particles in each group. */
  for (size_t i = 0; i < nr_gparts; i++) group_sizes[fof_find(i, pid)]++;

  fof_dump_group_data("fof_output_tree_serial.dat", nr_gparts, pid,
                      group_sizes);

  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_id = 0;
  for (size_t i = 0; i < nr_gparts; i++) {

    if (group_sizes[i] > 1) num_parts_in_groups += group_sizes[i];
    if (group_sizes[i] > max_group_size) {
      max_group_size = group_sizes[i];
      max_group_id = i;
    }
  }

  message(
      "No. of groups: %d. No. of particles in groups: %d. No. of particles not "
      "in groups: %ld.",
      num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_id);

  free(pid);
  free(group_sizes);
}

/* Dump FOF group data. */
void fof_dump_group_data(char *out_file, const size_t nr_gparts, int *pid,
                         int *group_sizes) {

  FILE *file = fopen(out_file, "w");
  fprintf(file, "# %7s %7s %7s\n", "ID", "Root ID", "Group Size");
  fprintf(file, "#-------------------------------\n");

  for (size_t i = 0; i < nr_gparts; i++) {
    fprintf(file, "  %7ld %7d %7d\n", i, pid[i], group_sizes[i]);
  }

  fclose(file);
}
