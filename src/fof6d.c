/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 James Willis (james.s.willis@durham.ac.uk)
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

#ifdef WITH_FOF

/* Some standard headers. */
#include <errno.h>
#include <libgen.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "fof6d.h"

/* Local headers. */
#include "fof.h"
#include "black_holes.h"
#include "common_io.h"
#include "engine.h"
#include "hashmap.h"
#include "memuse.h"
#include "proxy.h"
#include "threadpool.h"

void fof6d_calc_vel_disp(struct fof_props *props, struct space *s, const size_t num_parts_in_groups) {

  const int num_groups = props->num_groups;
  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;
  //size_t *group_index = props->group_index;
  //size_t *group_size = props->group_size;
  double *group_mass = props->group_mass;
  double *v_disp = NULL;
  double *v_mean[3] = {NULL, NULL, NULL};
  size_t *part_index = NULL;

  /* Allocate and initialise a velocity dispersion array. */
  if (swift_memalign("fof6d_v_disp", (void **)&v_disp, SWIFT_STRUCT_ALIGNMENT,
                     num_groups * sizeof(double)) != 0)
    error("Failed to allocate list of group velocity dispersions for 6DFOF search.");

  if (swift_memalign("fof6d_v_mean[0]", (void **)&v_mean[0], SWIFT_STRUCT_ALIGNMENT,
                     num_groups * sizeof(double)) != 0)
    error("Failed to allocate list of mean velocity for 6DFOF search.");

  if (swift_memalign("fof6d_v_mean[1]", (void **)&v_mean[1], SWIFT_STRUCT_ALIGNMENT,
                     num_groups * sizeof(double)) != 0)
    error("Failed to allocate list of mean velocity for 6DFOF search.");

  if (swift_memalign("fof6d_v_mean[2]", (void **)&v_mean[2], SWIFT_STRUCT_ALIGNMENT,
                     num_groups * sizeof(double)) != 0)
    error("Failed to allocate list of mean velocity for 6DFOF search.");

  if (swift_memalign("fof6d_part_index", (void **)&part_index, SWIFT_STRUCT_ALIGNMENT,
                     num_parts_in_groups * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle indices in groups for 6DFOF search.");

  bzero(v_disp, num_groups * sizeof(double));
  bzero(v_mean[0], num_groups * sizeof(double));
  bzero(v_mean[1], num_groups * sizeof(double));
  bzero(v_mean[2], num_groups * sizeof(double));
  bzero(part_index, num_parts_in_groups * sizeof(size_t));
  
  size_t part_ctr = 0;

  /* Calculate the mean velocity for each group. */
  for (size_t i = 0; i < nr_gparts; i++) {
    
    const size_t group_id = gparts[i].fof_data.group_id;

    if(group_id != fof_props_default_group_id) {
       v_mean[0][group_id - 1] += gparts[i].mass * gparts[i].v_full[0]; 
       v_mean[1][group_id - 1] += gparts[i].mass * gparts[i].v_full[1]; 
       v_mean[2][group_id - 1] += gparts[i].mass * gparts[i].v_full[2];

       /* JSW TODO: Could be calculated in fof_search_tree */
       part_index[part_ctr++] = i;
    }
  }
  
  for (int i = 0; i < num_groups; i++) {
    const double one_over_mass = 1.0 / group_mass[i];

    v_mean[0][i] *= one_over_mass; 
    v_mean[1][i] *= one_over_mass; 
    v_mean[2][i] *= one_over_mass;
  }

  /* Calculate the velocity dispersion for each group. */
  for (size_t i = 0; i < num_parts_in_groups; i++) {
    
    const size_t index = part_index[i];
    const size_t group_id = gparts[index].fof_data.group_id - 1;

    const double v_diff[3] = {gparts[index].v_full[0] - v_mean[0][group_id],
                              gparts[index].v_full[1] - v_mean[1][group_id],
                              gparts[index].v_full[2] - v_mean[2][group_id]};

    v_disp[group_id] += (v_diff[0] * v_diff[0] + 
                         v_diff[1] * v_diff[1] + 
                         v_diff[2] * v_diff[2]) * gparts[index].mass;
  }

  for (int i = 0; i < num_groups; i++) {
    v_disp[i] /= group_mass[i];
    //message("v_disp[%d]: %lf, v_mean: [%lf, %lf, %lf], group_mass: %lf", i, v_disp[i], v_mean[0][i], v_mean[1][i], v_mean[2][i], group_mass[i]);
  }

  fof6d_split_groups(props, s, num_parts_in_groups, v_disp, part_index);

  swift_free("fof6d_v_disp", v_disp);
  swift_free("fof6d_v_mean[0]", v_mean[0]);
  swift_free("fof6d_v_mean[1]", v_mean[1]);
  swift_free("fof6d_v_mean[2]", v_mean[2]);
  swift_free("fof6d_part_index", part_index);

}


void fof6d_n2_search(struct fof_6d *groups, struct space *s, const int num_groups, const size_t num_parts_in_groups, const double *v_disp, size_t *group_index, const size_t *group_size, const double l_x2) {

  double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Perform a neighbour search over each group. */ 
  for (int i = 0; i < num_groups; i++) {
  
    const double l_v2 = v_disp[i];
    const double l_xv2 = l_x2 * l_v2;

    message("Performing N^2 neighbour search over group %d which has %zu particles.", i, group_size[i]);

    for (size_t j = 0; j < group_size[i]; j++) {
    
      struct gpart *pi = groups[i].gparts[j];
      const double pix = pi->x[0], piy = pi->x[1], piz = pi->x[2];
      const double vix = pi->v_full[0], viy = pi->v_full[1], viz = pi->v_full[2];
    
      /* Find the root of pi. */
      size_t root_i = fof_find(group_index[i], group_index);
      
      for (size_t k = j + 1; k < group_size[i]; k++) {
      
        struct gpart *pj = groups[i].gparts[k];
        const double pjx = pj->x[0], pjy = pj->x[1], pjz = pj->x[2];
        const double vjx = pj->v_full[0], vjy = pj->v_full[1], vjz = pj->v_full[2];

        /* Find the root of pj. */
        const size_t root_j = fof_find(group_index[j], group_index);

        /* Compute pairwise distance, remembering to account for boundary
         * conditions. */
        float dx[3] = {nearest(pix - pjx, dim[0]), 
                       nearest(piy - pjy, dim[1]), 
                       nearest(piz - pjz, dim[2])};

        /* Compute pairwise velocity difference. */
        float dv[3] = {vix - vjx, viy - vjy, viz - vjz};
        
        float dx2 = 0.f;
        for (int l = 0; l < 3; l++) dx2 += dx[l] * dx[l];

        float dv2 = 0.f;
        for (int l = 0; l < 3; l++) dv2 += dv[l] * dv[l];

        /* Hit or miss? */
        if ((dx2 * l_v2 + dv2 * l_x2) < l_xv2) {
        
          //message("dx2: %f, dv2: %f, l_x2: %lf, l_v2: %lf, l_x2 * l_v2: %lf", dx2, dv2, l_x2, l_v2, l_xv2);

          /* Merge the groups */
          fof_union(&root_i, root_j, group_index);
        }
      }
    }
  }

}

void fof6d_split_groups(struct fof_props *props, struct space *s, const size_t num_parts_in_groups, const double *v_disp, const size_t *part_index) {

  const int num_groups = props->num_groups;
  struct gpart *gparts = s->gparts;
  //const size_t nr_gparts = s->nr_gparts;
  size_t *group_index = NULL;
  //size_t *group_size = props->group_size;
  size_t *group_size = NULL;
  //double *group_mass = props->group_mass;
  struct fof_6d groups[num_groups];
  size_t *group_part_ctr = NULL;
  const double l_x2 = props->l_x2;

  if (swift_memalign("fof6d_group_part_ctr", (void **)&group_part_ctr, SWIFT_STRUCT_ALIGNMENT,
                     num_groups * sizeof(size_t)) != 0)
    error("Failed to allocate list of group particle counters for 6DFOF search.");

  if (swift_memalign("fof6d_group_size", (void **)&group_size, SWIFT_STRUCT_ALIGNMENT,
                     num_groups * sizeof(size_t)) != 0)
    error("Failed to allocate list of group sizes for 6DFOF search.");


  bzero(group_size, num_groups * sizeof(size_t));

  /* Get ptrs to particles in groups */
  for (size_t i = 0; i < num_parts_in_groups; i++) {
    
    const size_t index = part_index[i];
    const size_t group_id = gparts[index].fof_data.group_id - 1;
    group_size[group_id]++;
  }
 
  /* Allocate gpart* arrays */
  for (int i = 0; i < num_groups; i++) {
    if (posix_memalign((void **)&groups[i].gparts, 32,
          group_size[i] * sizeof(struct gpart*)) != 0)
      error("Failed to allocate list of group masses for FOF search.");
  }

  bzero(group_part_ctr, num_groups * sizeof(size_t));

  /* Allocate and initialise a new group index array. */
  if (swift_memalign("fof6d_group_index", (void **)&group_index, SWIFT_STRUCT_ALIGNMENT,
                     num_parts_in_groups * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle group indices for 6DFOF search.");

  /* Set initial group index */
  threadpool_map(&s->e->threadpool, fof_set_initial_group_index_mapper,
                 group_index, num_parts_in_groups, sizeof(size_t),
                 threadpool_auto_chunk_size, group_index);

  /* Get ptrs to particles in groups */
  for (size_t i = 0; i < num_parts_in_groups; i++) {
    
    const size_t index = part_index[i];
    const size_t group_id = gparts[index].fof_data.group_id - 1;
    const size_t part_ctr = group_part_ctr[group_id];
    groups[group_id].gparts[part_ctr] = &gparts[index];
    group_part_ctr[group_id] = group_part_ctr[group_id] + 1;
  }
 
  fof6d_n2_search(groups, s, num_groups, num_parts_in_groups, v_disp, group_index, group_size, l_x2);

  int num_6d_groups = 0;
  for (size_t i = 0; i < num_parts_in_groups; i++) {
    if(group_index[i] == i) num_6d_groups++;
  }

  size_t *group_6d_size = NULL;

  if (swift_memalign("fof6d_group_6d_size", (void **)&group_6d_size, SWIFT_STRUCT_ALIGNMENT,
                     num_parts_in_groups * sizeof(size_t)) != 0)
    error("Failed to allocate list of group sizes for 6DFOF search.");

  bzero(group_6d_size, num_parts_in_groups * sizeof(size_t));

  for (size_t i = 0; i < num_parts_in_groups; i++) {
    size_t root = fof_find(group_index[i], group_index);
    group_6d_size[root]++;
  }

  struct group_length *high_group_sizes = NULL;
  int group_count = 0;

  if (swift_memalign("fof6d_high_group_sizes", (void **)&high_group_sizes, SWIFT_STRUCT_ALIGNMENT,
                     num_6d_groups * sizeof(struct group_length)) != 0)
    error("Failed to allocate list of large groups.");

  for (size_t i = 0; i < num_parts_in_groups; i++) {
    gparts[part_index[i]].fof_data.group_id = fof_props_default_group_id;
    if(group_6d_size[i] > 1) {
      gparts[part_index[i]].fof_data.group_id = group_count;
      high_group_sizes[group_count].index = i;
      high_group_sizes[group_count++].size = group_6d_size[i];
    }
  }

  for (size_t i = 0; i < num_parts_in_groups; i++) {
    size_t root = fof_find(group_index[i], group_index);
    gparts[part_index[i]].fof_data.group_id = gparts[part_index[root]].fof_data.group_id;
  }
  
  FILE *file = fopen("fof6d_output.dat", "w");

  fprintf(file, "# %8s %12s %12s %18s\n", "Group ID",
          "Group Size", "Group Mass",
          "Particle ID");
  fprintf(file,
          "#-------------------------------------------------------------------"
          "-------------\n");

  //for (int i = 0; i < num_6d_groups; i++) {
  for (int i = 0; i < group_count; i++) {

    //const size_t group_offset = group_sizes[i].index;
    fprintf(file, "  %8d %12zu %12e %18d\n",
            i, high_group_sizes[i].size, 0.0, 0);
  }
  
  fclose(file);

  message("No. of 3D FoF groups: %d", num_groups);
  message("No. of 6D FoF groups: %d", num_6d_groups);

  swift_free("fof6d_group_size", group_size);
  swift_free("fof6d_group_6d_size", group_6d_size);

  for (int i = 0; i < num_groups; i++) {
    free(groups[i].gparts);
  }
}
#endif /* WITH_FOF6D */
