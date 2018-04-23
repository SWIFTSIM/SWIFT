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

__attribute__((always_inline)) INLINE int fof_find(const int i, const int *pid) {
    
  int root = i;
  while(root != pid[root])
    root = pid[root];

  return root;
}

__attribute__((always_inline)) INLINE int fof_find_path_comp(const int i, int *pid) {
    
  int root = i;
  while(root != pid[root])
    root = pid[root];

  /* Perform path compression. */
  int index = i;
  while(index != root) {
    int next = pid[index];
    pid[index] = root;
    index = next;
  }
  
  return root;
}

void fof_search_naive(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;
  int *pid;
  int *num_in_groups;
  int num_groups = nr_gparts;

  message("Searching %ld gravity particles for links with l_x2: %lf", nr_gparts, l_x2);

  /* Allocate and populate array of particle group IDs. */
  if (posix_memalign((void **)&pid, 32,
        nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");

  for(size_t i=0; i<nr_gparts; i++) pid[i] = i;    

  if (posix_memalign((void **)&num_in_groups, 32,
        nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of number in groups for FOF search.");

  for(size_t i=0; i<nr_gparts; i++) num_in_groups[i] = 1;    

  for(size_t runs=0; runs<2; runs++) 
  /* Loop over particles and find which particles belong in the same group. */
  for(size_t i=0; i<nr_gparts; i++) {
  
    //message("Searching for particle: %ld groups.", i);

    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    for(size_t j=0; j<nr_gparts; j++) { 

      /* Skip yourself. */
      if(i == j) continue;

      struct gpart *pj = &gparts[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance, remembering to account for boundary conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) {
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < l_x2 && pid[j] < pid[i]) {

        num_in_groups[pid[i]]--;
        num_in_groups[pid[j]]++;
        
        pid[i] = pid[j]; 

        num_groups--;
      }
    }
  }

  fof_dump_group_data("fof_output_naive.dat", nr_gparts, pid, num_in_groups);
  
  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_id = 0;
  for(size_t i=0; i<nr_gparts; i++) {

      if(num_in_groups[i] > 1) num_parts_in_groups += num_in_groups[i];
      if( num_in_groups[i] > max_group_size) {
        max_group_size = num_in_groups[i];
        max_group_id = i;
      }
  }
  
  message("No. of groups: %d. No. of particles in groups: %d. No. of particles not in groups: %ld.", num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_id);

  free(pid);
  free(num_in_groups);
}

void fof_search_serial(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;
  int *pid;
  int *num_in_groups;
  int num_groups = nr_gparts;

  message("Searching %ld gravity particles for links with l_x2: %lf", nr_gparts, l_x2);

  /* Allocate and populate array of particle group IDs. */
  if (posix_memalign((void **)&pid, 32,
        nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");

  for(size_t i=0; i<nr_gparts; i++) pid[i] = i;    

  if (posix_memalign((void **)&num_in_groups, 32,
        nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of number in groups for FOF search.");

  for(size_t i=0; i<nr_gparts; i++) num_in_groups[i] = 1;    

  /* Loop over particles and find which particles belong in the same group. */
  for(size_t i=0; i<nr_gparts; i++) {
  
    //message("Searching for particle: %ld groups.", i);

    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    for(size_t j=0; j<nr_gparts; j++) { 

      /* Skip yourself. */
      if(i == j) continue;

      struct gpart *pj = &gparts[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance, remembering to account for boundary conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) {
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Find the roots of pi and pj. */
      const int root_i = fof_find(i, pid);
      const int root_j = fof_find(j, pid);
      
      //const int root_i = fof_find_path_comp(i, pid);
      //const int root_j = fof_find_path_comp(j, pid);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Hit or miss? */
      if (r2 < l_x2)  {

        if(root_j < root_i) {
          pid[root_i] = root_j; 
          num_in_groups[root_i]--;
          num_in_groups[root_j]++;
        }
        else {
          pid[root_j] = root_i;
          num_in_groups[root_j]--;
          num_in_groups[root_i]++;
        }
        
        num_groups--;
      }
    }
  }

  fof_dump_group_data("fof_output_serial.dat", nr_gparts, pid, num_in_groups);

  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_id = 0;
  for(size_t i=0; i<nr_gparts; i++) {

      if(num_in_groups[i] > 1) num_parts_in_groups += num_in_groups[i];
      if( num_in_groups[i] > max_group_size) {
        max_group_size = num_in_groups[i];
        max_group_id = i;
      }
  }


  message("No. of groups: %d. No. of particles in groups: %d. No. of particles not in groups: %ld.", num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_id);

  free(pid);
  free(num_in_groups);
}

void fof_search_cell(struct space *s, struct cell *c, int *pid, int *num_in_groups, int *num_groups) {

  const size_t count = c->gcount;
  struct gpart *gparts = c->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;

  /* Loop over particles and find which particles belong in the same group. */
  for(size_t i=0; i<count; i++) {
  
    struct gpart *pi = &gparts[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];
    const size_t offset_i = pi->offset;

    for(size_t j=0; j<count; j++) { 

      /* Skip yourself. */
      if(i == j) continue;
      
      struct gpart *pj = &gparts[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const size_t offset_j = pj->offset;

      /* Compute pairwise distance, remembering to account for boundary conditions. */
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
      
      //const int root_i = fof_find_path_comp(offset_i, pid);
      //const int root_j = fof_find_path_comp(offset_j, pid);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Hit or miss? */
      if (r2 < l_x2)  {

        if(root_j < root_i) {
          pid[root_i] = root_j; 
          num_in_groups[root_i]--;
          num_in_groups[root_j]++;
        }
        else {
          pid[root_j] = root_i;
          num_in_groups[root_j]--;
          num_in_groups[root_i]++;
        }
        
        (*num_groups)--;
      }
    }
  }

}

void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj, int *pid, int *num_in_groups, int *num_groups) {

  const size_t count_i = ci->gcount;
  const size_t count_j = cj->gcount;
  struct gpart *gparts_i = ci->gparts;
  struct gpart *gparts_j = cj->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double l_x2 = s->l_x2;

  /* Loop over particles and find which particles belong in the same group. */
  for(size_t i=0; i<count_i; i++) {
  
    struct gpart *pi = &gparts_i[i];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];
    const size_t offset_i = pi->offset;

    for(size_t j=0; j<count_j; j++) { 

      struct gpart *pj = &gparts_j[j];
      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const size_t offset_j = pj->offset;

      /* Compute pairwise distance, remembering to account for boundary conditions. */
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
      
      //const int root_i = fof_find_path_comp(offset_i, pid);
      //const int root_j = fof_find_path_comp(offset_j, pid);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      /* Hit or miss? */
      if (r2 < l_x2)  {

        if(root_j < root_i) {
          pid[root_i] = root_j; 
          num_in_groups[root_i]--;
          num_in_groups[root_j]++;
        }
        else {
          pid[root_j] = root_i;
          num_in_groups[root_j]--;
          num_in_groups[root_i]++;
        }
        
        (*num_groups)--;
      }
    }
  }

  /* Loop over particles and find which particles belong in the same group. */
  for(size_t j=0; j<count_j; j++) {
  
    //message("Searching for particle: %ld groups.", i);

    struct gpart *pj = &gparts_j[j];
    const double pjx = pj->x[0];
    const double pjy = pj->x[1];
    const double pjz = pj->x[2];
    const size_t offset_j = pj->offset;

    for(size_t i=0; i<count_i; i++) { 

      struct gpart *pi = &gparts_i[i];
      const double pix = pi->x[0];
      const double piy = pi->x[1];
      const double piz = pi->x[2];
      const size_t offset_i = pi->offset;

      /* Compute pairwise distance, remembering to account for boundary conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = pjx - pix;
      dx[1] = pjy - piy;
      dx[2] = pjz - piz;

      for (int k = 0; k < 3; k++) {
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Find the roots of pi and pj. */
      const int root_j = fof_find(offset_j, pid);
      const int root_i = fof_find(offset_i, pid);
      
      //const int root_i = fof_find_path_comp(i, pid);
      //const int root_j = fof_find_path_comp(j, pid);

      /* Skip particles in the same group. */
      if (root_j == root_i) continue;

      /* Hit or miss? */
      if (r2 < l_x2)  {

        if(root_i < root_j) {
          pid[root_j] = root_i; 
          num_in_groups[root_j]--;
          num_in_groups[root_i]++;
        }
        else {
          pid[root_i] = root_j;
          num_in_groups[root_i]--;
          num_in_groups[root_j]++;
        }
        
        (*num_groups)--;
      }
    }
  }
}

void fof_search_tree_serial(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_cells = s->nr_cells;
  int *pid, *num_in_groups;
  int num_groups = nr_gparts;
  struct gpart *gparts = s->gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->cell_search_r2;

  message("Searching %ld gravity particles for links with l_x2: %lf", nr_gparts, s->l_x2);

  /* Allocate and populate array of particle group IDs. */
  if (posix_memalign((void **)&pid, 32,
        nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of particle group IDs for FOF search.");

  for(size_t i=0; i<nr_gparts; i++) {
    gparts[i].offset = i;
    pid[i] = i;
  }

  if (posix_memalign((void **)&num_in_groups, 32,
        nr_gparts * sizeof(int)) != 0)
    error("Failed to allocate list of number in groups for FOF search.");

  for(size_t i=0; i<nr_gparts; i++) num_in_groups[i] = 1;    

  /* Loop over particles and find which particles belong in the same group. */
  for(int cid=0; cid<nr_cells; cid++) {
    
    struct cell *restrict ci = &s->cells_top[cid];
    const double cix = ci->loc[0];
    const double ciy = ci->loc[1];
    const double ciz = ci->loc[2];
 
    fof_search_cell(s, ci, pid, num_in_groups, &num_groups);
    
    for(int cjd=cid + 1; cjd<nr_cells; cjd++) {      
    
      struct cell *restrict cj = &s->cells_top[cjd];
      const double cjx = cj->loc[0];
      const double cjy = cj->loc[1];
      const double cjz = cj->loc[2];

      /* Compute distance between cells, remembering to account for boundary conditions. */
      float dx[3], r2 = 0.0f;
      dx[0] = cix - cjx;
      dx[1] = ciy - cjy;
      dx[2] = ciz - cjz;

      for (int k = 0; k < 3; k++) {
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Perform FOF search between pairs of cells that are within the linking length. */
      if (r2 < search_r2)  {
        fof_search_pair_cells(s, ci, cj, pid, num_in_groups, &num_groups);
      }
    }
  }
  
  fof_dump_group_data("fof_output_tree_serial.dat", nr_gparts, pid, num_in_groups);

  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_id = 0;
  for(size_t i=0; i<nr_gparts; i++) {

      if(num_in_groups[i] > 1) num_parts_in_groups += num_in_groups[i];
      if( num_in_groups[i] > max_group_size) {
        max_group_size = num_in_groups[i];
        max_group_id = i;
      }
  }

  message("No. of groups: %d. No. of particles in groups: %d. No. of particles not in groups: %ld.", num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_id);

  free(pid);
  free(num_in_groups);
}

/* Dump FOF group data. */
void fof_dump_group_data(char *out_file, const int nr_gparts, int *pid, int *num_in_groups) {

  FILE *file = fopen(out_file,"w");
  fprintf(file, "# %7s %7s %7s\n", "ID", "Root ID","Group Size"); 
  fprintf(file, "#-------------------------------\n"); 

  for(size_t i=0; i<nr_gparts; i++) {

      fprintf(file, "  %7ld %7d %7d\n", i, pid[i], num_in_groups[i]); 

  }
  
  fclose(file);
  
}
