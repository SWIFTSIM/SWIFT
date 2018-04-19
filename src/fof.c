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

void fof_search(struct engine *e) {

  const size_t nr_gparts = e->s->nr_gparts;
  struct gpart *gparts = e->s->gparts;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const double l_x = 0.2 * (dim[0] / pow(nr_gparts, 1./3.));
  const double l_x2 = l_x * l_x;
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

      /* Hit or miss? */
      if (r2 < l_x2 && pid[j] < pid[i]) {
          pid[i] = pid[j]; 
          num_in_groups[i]--;
          num_in_groups[j]++;
          num_groups--;
      }
    }
  }

  int num_parts_in_groups = 0;
  int max_group_size = 0, max_group_id = 0;
  for(size_t i=0; i<nr_gparts; i++) {

      if(num_in_groups[i] > 1) num_parts_in_groups += num_in_groups[i];
      if( num_in_groups[i] > max_group_size) {
        max_group_size = num_in_groups[i];
        max_group_id = i;
      }
  }
  
  message("No. of groups: %d. No. of particles in groups: %d. No. of particles not in groups: %d.", num_groups, num_parts_in_groups, nr_gparts - num_parts_in_groups);
  message("Biggest group size: %d with ID: %d", max_group_size, max_group_id);

  free(pid);
  free(num_in_groups);
}
