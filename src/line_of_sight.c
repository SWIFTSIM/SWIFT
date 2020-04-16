/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "engine.h"
#include "kernel_hydro.h"
#include "line_of_sight.h"
#include "periodic.h"

/**
 * @brief Generates random sightline positions.
 *
 * Independent sightlines are made for the XY, YZ and XZ planes.
 *
 * @param LOS Structure to store sightlines.
 * @param params Sightline parameters.
 */
void generate_line_of_sights(struct line_of_sight *Los,
                             struct line_of_sight_params *params) {

  /* Keep track of number of sightlines. */
  int count = 0;

  /* Sightlines in XY plane. */
  for (int i = 0; i < params->num_along_xy; i++) {
    Los[count].Xpos =
        ((float)rand() / (float)(RAND_MAX) * (params->xmax - params->xmin)) +
        params->xmin;
    Los[count].Ypos =
        ((float)rand() / (float)(RAND_MAX) * (params->ymax - params->ymin)) +
        params->ymin;
    Los[count].particles_in_los_local = 0;
    Los[count].particles_in_los_total = 0;
    Los[count].xaxis = 0;
    Los[count].yaxis = 1;
    Los[count].zaxis = 2;
    count += 1;
  }

  /* Sightlines in YZ plane. */
  for (int i = 0; i < params->num_along_yz; i++) {
    Los[count].Xpos =
        ((float)rand() / (float)(RAND_MAX) * (params->ymax - params->ymin)) +
        params->ymin;
    Los[count].Ypos =
        ((float)rand() / (float)(RAND_MAX) * (params->zmax - params->zmin)) +
        params->zmin;
    Los[count].particles_in_los_local = 0;
    Los[count].particles_in_los_total = 0;
    Los[count].xaxis = 1;
    Los[count].yaxis = 2;
    Los[count].zaxis = 0;
    count += 1;
  }

  /* Sightlines in XZ plane. */
  for (int i = 0; i < params->num_along_xz; i++) {
    Los[count].Xpos =
        ((float)rand() / (float)(RAND_MAX) * (params->xmax - params->xmin)) +
        params->xmin;
    Los[count].Ypos =
        ((float)rand() / (float)(RAND_MAX) * (params->zmax - params->zmin)) +
        params->zmin;
    Los[count].particles_in_los_local = 0;
    Los[count].particles_in_los_total = 0;
    Los[count].xaxis = 0;
    Los[count].yaxis = 2;
    Los[count].zaxis = 1;
    count += 1;
  }
}

/**
 * @brief Print line_of_sight structure (for debugging).
 *
 * @param Los Structure of sightlines to print.
 * @param params Sightline parameters.
 */
void print_los_info(struct line_of_sight *Los,
                    struct line_of_sight_params *params) {

  printf("\nPrinting LOS information...\n");
  for (int i = 0; i < params->num_tot; i++) {
    printf("[LOS %i] Xpos:%g Ypos:%g particles_in_los_total:%i\n", i,
           Los[i].Xpos, Los[i].Ypos, Los[i].particles_in_los_total);
  }
}

void do_line_of_sight(struct engine *e) {

  const struct space *s = e->s;
  const size_t nr_parts = s->nr_parts;

  /* LOS params, dummy until included in parameter file. */
  struct line_of_sight_params LOS_params = {.num_along_xy = 1,
                                            .num_along_yz = 1,
                                            .num_along_xz = 1,
                                            .xmin = 0,
                                            .xmax = 10,
                                            .ymin = 0,
                                            .ymax = 10,
                                            .zmin = 0,
                                            .zmax = 10};
  LOS_params.num_tot = LOS_params.num_along_xy + LOS_params.num_along_yz +
                       LOS_params.num_along_xz;

  /* Start by generating the random sightline positions. */
  struct line_of_sight *LOS_list = (struct line_of_sight *)malloc(
      LOS_params.num_tot * sizeof(struct line_of_sight));
  if (e->nodeID == 0) generate_line_of_sights(LOS_list, &LOS_params);
#ifdef WITH_MPI
  MPI_Bcast(LOS_list, LOS_params.num_tot * sizeof(struct line_of_sight),
            MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  double dx, dy, r2, hsml;

  FILE *f;

  if (e->nodeID == 0) {
    char filename[200];
    sprintf(filename, "los_%04i.csv", engine_current_step);
    f = fopen(filename, "w");
    if (f == NULL) error("Error opening los file.");

    /* Write header */
    fprintf(f, "los_id,x,y,z,hsml\n");
  }

  /* Loop over each random LOS. */
  for (int j = 0; j < LOS_params.num_tot; j++) {

    /* Loop over each gas particle to find those in LOS. */
    for (size_t i = 0; i < nr_parts; i++) {
      if (s->parts[i].gpart->type == swift_type_gas) {
        /* Distance from this particle to LOS along x dim. */
        dx = s->parts[i].x[LOS_list[j].xaxis] - LOS_list[j].Xpos;

        /* Periodic wrap. */
        if (e->s->periodic) dx = nearest(dx, s->dim[LOS_list[j].xaxis]);

        /* Smoothing length of this part. */
        hsml = s->parts[i].h * kernel_gamma;

        /* Does this particle fall into our LOS? */
        if (dx <= hsml) {
          /* Distance from this particle to LOS along y dim. */
          dy = s->parts[i].x[LOS_list[j].yaxis] - LOS_list[j].Ypos;
          
          /* Periodic wrap. */
          if (e->s->periodic) dy = nearest(dy, s->dim[LOS_list[j].yaxis]);

          /* Does this particle still fall into our LOS? */
          if (dy <= hsml) {
            /* 2D distance to LOS. */
            r2 = dx * dx + dy * dy;

            if (r2 <= hsml * hsml) {

              /* We've found one. */
              LOS_list[j].particles_in_los_local += 1;
            }
          }
        }
      }
    } /* End of loop over all gas particles */

#ifdef WITH_MPI
    int LOS_counts[e->nr_nodes];
    int LOS_disps[e->nr_nodes];

    /* How many particles does each rank have for this LOS? */
    MPI_Allgather(&LOS_list[j].particles_in_los_local, 1, MPI_INT,
                  &LOS_counts, 1, MPI_INT, MPI_COMM_WORLD);

    for (int k = 0, disp_count = 0; k < e->nr_nodes; k++) {
      /* Total particles in this LOS. */
      LOS_list[j].particles_in_los_total += LOS_counts[k];

      /* Counts and disps for Gatherv. */
      LOS_disps[k] = disp_count;
      disp_count += LOS_counts[k];
    }
#else
    LOS_list[j].particles_in_los_total = LOS_list[j].particles_in_los_local;
#endif

    /* Setup los particles structure. */
    struct line_of_sight_particles *LOS_particles;
    if (e->nodeID == 0) {
      LOS_particles = (struct line_of_sight_particles *)malloc(
          LOS_list[j].particles_in_los_total *
          sizeof(struct line_of_sight_particles));
    } else { 
      LOS_particles = (struct line_of_sight_particles *)malloc(
          LOS_list[j].particles_in_los_local *
          sizeof(struct line_of_sight_particles));
    }

    /* Loop over each gas particle again to store properties. */
    int count = 0;

    for (size_t i = 0; i < nr_parts; i++) {

      if (s->parts[i].gpart->type == swift_type_gas) {
        /* Distance from this particle to LOS along x dim. */
        dx = s->parts[i].x[LOS_list[j].xaxis] - LOS_list[j].Xpos;

        /* Periodic wrap. */
        if (e->s->periodic) dx = nearest(dx, s->dim[LOS_list[j].xaxis]);

        /* Smoothing length of this part. */
        hsml = s->parts[i].h * kernel_gamma;

        /* Does this particle fall into our LOS? */
        if (dx <= hsml) {
          /* Distance from this particle to LOS along y dim. */
          dy = s->parts[i].x[LOS_list[j].yaxis] - LOS_list[j].Ypos;
        
          /* Periodic wrap. */
          if (e->s->periodic) dy = nearest(dy, s->dim[LOS_list[j].yaxis]);

          /* Does this particle still fall into our LOS? */
          if (dy <= hsml) {
            /* 2D distance to LOS. */
            r2 = dx * dx + dy * dy;

            if (r2 <= hsml * hsml) {

              /* Store particle properties. */
              LOS_particles[count].pos[0] = s->parts[i].x[0];
              LOS_particles[count].pos[1] = s->parts[i].x[1];
              LOS_particles[count].pos[2] = s->parts[i].x[2];

              LOS_particles[count].h = s->parts[i].h;

              count++;
            }
          }
        }
      }
    } /* End of loop over all gas particles. */

#ifdef WITH_MPI
    /* Collect all particles in LOS to rank 0. */
    if (e->nodeID == 0) {
      MPI_Gatherv(MPI_IN_PLACE, 0,
                  lospart_mpi_type, LOS_particles, LOS_counts, LOS_disps,
                  lospart_mpi_type, 0, MPI_COMM_WORLD);
    } else {
      MPI_Gatherv(LOS_particles, LOS_list[j].particles_in_los_local,
                  lospart_mpi_type, LOS_particles, LOS_counts, LOS_disps,
                  lospart_mpi_type, 0, MPI_COMM_WORLD);
    }
#endif
    
    /* Write particles to file. */
    if (e->nodeID == 0) {
      for (int kk = 0; kk < LOS_list[j].particles_in_los_total; kk++) {
        fprintf(f, "%i,%g,%g,%g,%g\n", j, LOS_particles[kk].pos[0],
                LOS_particles[kk].pos[1], LOS_particles[kk].pos[2],
                LOS_particles[kk].h);
      }
    }
  } /* End of loop over each LOS */

  if (e->nodeID == 0) fclose(f);

  //#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) print_los_info(LOS_list, &LOS_params);
  //#endif
}
