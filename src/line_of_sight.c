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
#include "io_properties.h"
#include "hydro_io.h"

void los_init(double dim[3], struct los_props *los_params,
        struct swift_params *params) {
  /* How many line of sights in each plane. */
  los_params->num_along_xy =
      parser_get_opt_param_int(params, "LineOfSight:num_along_xy", 0);
  los_params->num_along_yz =
      parser_get_opt_param_int(params, "LineOfSight:num_along_yz", 0);
  los_params->num_along_xz =
      parser_get_opt_param_int(params, "LineOfSight:num_along_xz", 0);

  /* Min/max range across x,y and z where random LOS's are allowed. */
  los_params->xmin =
            parser_get_opt_param_double(params, "LineOfSight:xmin", 0.);
  los_params->xmax =
            parser_get_opt_param_double(params, "LineOfSight:xmax", dim[0]);
  los_params->ymin =
            parser_get_opt_param_double(params, "LineOfSight:ymin", 0.);
  los_params->ymax =
            parser_get_opt_param_double(params, "LineOfSight:ymax", dim[1]);
  los_params->zmin =
            parser_get_opt_param_double(params, "LineOfSight:zmin", 0.);
  los_params->zmax =
            parser_get_opt_param_double(params, "LineOfSight:zmax", dim[2]);

  /* Compute total number of sightlines. */
  los_params->num_tot = los_params->num_along_xy +
                        los_params->num_along_yz +
                        los_params->num_along_xz;
} 

/**
 * @brief Generates random sightline positions.
 *
 * Independent sightlines are made for the XY, YZ and XZ planes.
 *
 * @param LOS Structure to store sightlines.
 * @param params Sightline parameters.
 */
void generate_line_of_sights(struct line_of_sight *Los,
                             const struct los_props *params) {

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
void print_los_info(const struct line_of_sight *Los,
                    const struct los_props *params) {

  printf("\nPrinting LOS information...\n");
  for (int i = 0; i < params->num_tot; i++) {
    printf("[LOS %i] Xpos:%g Ypos:%g particles_in_los_total:%i\n", i,
           Los[i].Xpos, Los[i].Ypos, Los[i].particles_in_los_total);
  }
}

void do_line_of_sight(struct engine *e) {

  /* Start counting. */
  const ticks tic = getticks();

  const struct space *s = e->s;
  const struct part *p = s->parts;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const size_t nr_parts = s->nr_parts;
  const struct los_props *LOS_params = e->los_properties;

  double dx, dy, r2, hsml;

  /* HDF5 vars. */
  hid_t h_file, h_grp;
  char fileName[200], groupName[200];

  /* Start by generating the random sightline positions. */
  struct line_of_sight *LOS_list = (struct line_of_sight *)malloc(
      LOS_params->num_tot * sizeof(struct line_of_sight));
  if (e->nodeID == 0) generate_line_of_sights(LOS_list, LOS_params);
#ifdef WITH_MPI
  MPI_Bcast(LOS_list, LOS_params->num_tot * sizeof(struct line_of_sight),
            MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Create HDF5 file. */
  if (e->nodeID == 0) {
    sprintf(fileName, "los_%04i.hdf5", e->los_output_count);
    h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Error while opening file '%s'.", fileName);
  }

  /* Loop over each random LOS. */
  for (int j = 0; j < LOS_params->num_tot; j++) {

    /* Loop over each gas particle to find those in LOS. */
    for (size_t i = 0; i < nr_parts; i++) {
      if (p[i].gpart->type == swift_type_gas) {
        /* Distance from this particle to LOS along x dim. */
        dx = p[i].x[LOS_list[j].xaxis] - LOS_list[j].Xpos;

        /* Periodic wrap. */
        if (periodic) dx = nearest(dx, dim[LOS_list[j].xaxis]);

        /* Smoothing length of this part. */
        hsml = p[i].h * kernel_gamma;

        /* Does this particle fall into our LOS? */
        if (dx <= hsml) {
          /* Distance from this particle to LOS along y dim. */
          dy = p[i].x[LOS_list[j].yaxis] - LOS_list[j].Ypos;
          
          /* Periodic wrap. */
          if (periodic) dy = nearest(dy, dim[LOS_list[j].yaxis]);

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
    struct part *LOS_particles;
    if (e->nodeID == 0) {
      LOS_particles = (struct part *)malloc(
          LOS_list[j].particles_in_los_total *
          sizeof(struct part));
    } else {
      LOS_particles = (struct part *)malloc(
          LOS_list[j].particles_in_los_local *
          sizeof(struct part));
    }

    /* Loop over each gas particle again to store properties. */
    int count = 0;

    for (size_t i = 0; i < nr_parts; i++) {

      if (p[i].gpart->type == swift_type_gas) {
        /* Distance from this particle to LOS along x dim. */
        dx = p[i].x[LOS_list[j].xaxis] - LOS_list[j].Xpos;

        /* Periodic wrap. */
        if (periodic) dx = nearest(dx, dim[LOS_list[j].xaxis]);

        /* Smoothing length of this part. */
        hsml = p[i].h * kernel_gamma;

        /* Does this particle fall into our LOS? */
        if (dx <= hsml) {
          /* Distance from this particle to LOS along y dim. */
          dy = p[i].x[LOS_list[j].yaxis] - LOS_list[j].Ypos;
        
          /* Periodic wrap. */
          if (periodic) dy = nearest(dy, dim[LOS_list[j].yaxis]);

          /* Does this particle still fall into our LOS? */
          if (dy <= hsml) {
            /* 2D distance to LOS. */
            r2 = dx * dx + dy * dy;

            if (r2 <= hsml * hsml) {

              /* Store particle properties. */
              LOS_particles[count].x[0] = p[i].x[0];
              LOS_particles[count].x[1] = p[i].x[1];
              LOS_particles[count].x[2] = p[i].x[2];

              LOS_particles[count].h = hsml;

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
                  part_mpi_type, LOS_particles, LOS_counts, LOS_disps,
                  part_mpi_type, 0, MPI_COMM_WORLD);
    } else {
      MPI_Gatherv(LOS_particles, LOS_list[j].particles_in_los_local,
                  part_mpi_type, LOS_particles, LOS_counts, LOS_disps,
                  part_mpi_type, 0, MPI_COMM_WORLD);
    }
#endif
    
    /* Write particles to file. */
    if (e->nodeID == 0) {
      sprintf(groupName, "/LOS_%04i", j);
      h_grp = H5Gcreate(h_file, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating LOS HDF5 group\n");
      write_los_hdf5_datasets(h_grp, j, LOS_list[j].particles_in_los_total, LOS_particles, e);
      H5Gclose(h_grp);
    }
  } /* End of loop over each LOS */

  if (e->nodeID == 0) {
    /* Close file */
    H5Fclose(h_file);
  }

  //#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) print_los_info(LOS_list, LOS_params);
  //#endif
 
  /* Up the count. */
  e->los_output_count++;

  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void write_los_hdf5_datasets(hid_t grp, int j, int N, struct part* parts,
        struct engine* e) {

  struct io_props p;

  /* Coordinates. */
  p = io_make_output_field_convert_part(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, parts, NULL,
      convert_part_pos, "Co-moving positions of the particles");
  write_los_hdf5_dataset(p, N, j, e, grp);
  bzero(&p, sizeof(struct io_props));

  /* Smoothing Lengths. */
  //p = io_make_output_field(
  //    "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, parts, h,
  //    "Co-moving smoothing lengths (FWHM of the kernel) of the particles");
  //message("%s", p.description);
  ////write_los_hdf5_dataset(p, N, j, e, grp);
  //bzero(&p, sizeof(struct io_props));
}

void write_los_hdf5_dataset(const struct io_props p, int N, int j, struct engine* e,
        hid_t grp) {

  hsize_t dims[2];
  hid_t dataset_id, dataspace_id;
  herr_t status;
  char att_name[200];

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     N * io_sizeof_type(p.type)) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy particle data to temp buffer */
  io_copy_temp_buffer(temp, e, p, N, e->internal_units, e->snapshot_units);

  /* Prepare dataset */
  dims[0] = N;
  dims[1] = p.dimension;
  dataspace_id = H5Screate_simple(2, dims, NULL);

  /* Write dataset */
  sprintf(att_name, "/LOS_%04i/%s", j, p.name);
  dataset_id = H5Dcreate(grp, att_name, io_hdf5_type(p.type), dataspace_id,    
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite(dataset_id, io_hdf5_type(p.type), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp);
  if (status < 0) error("Error while writing data array '%s'.", p.name);
 
  /* Free and close everything */
  swift_free("writebuff", temp);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

}
