/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>

//#include "cooling_io.h"
#include "engine.h"
#include "kernel_hydro.h"
#include "line_of_sight.h"
#include "periodic.h"
#include "io_properties.h"
#include "hydro_io.h"
#include "chemistry_io.h"
#include "fof_io.h"
#include "star_formation_io.h"
#include "tracers_io.h"
#include "velociraptor_io.h"

/**
 * @brief Reads the LOS properties from the param file.
 *
 * @param dim Space dimensions.
 * @param los_params Sightline parameters to save into.
 * @param params Swift params to read from.
 */
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

  /* Where are we saving them? */
  parser_get_param_string(params, "LineOfSight:basename", los_params->basename);
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

/**
 * @brief Main work function for computing line of sights.
 * 
 * 1) Construct random line of sight positions.
 * 2) Loop over each line of sight.
 * -  a) Loop over each part to see who falls within this LOS.
 * -  b) Use this count to construct a LOS parts array.
 * -  c) Loop over each part and extract those in LOS.
 * -  d) Save LOS particles to HDF5 file.
 * 
 * @param e The engine.
 */
void do_line_of_sight(struct engine *e) {

  /* Start counting. */
  const ticks tic = getticks();

  const struct space *s = e->s;
  const struct part *p = s->parts;
  const struct xpart *xp = s->xparts;
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
    sprintf(fileName, "%s_%04i.hdf5", LOS_params->basename, e->los_output_count);
    h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Error while opening file '%s'.", fileName);

    /* Write header */
    write_hdf5_header(h_file, e, LOS_params);
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
    struct xpart *LOS_particles_xparts;
    if (e->nodeID == 0) {
      LOS_particles = (struct part *)malloc(
          LOS_list[j].particles_in_los_total *
          sizeof(struct part));
      LOS_particles_xparts = (struct xpart *)malloc(
          LOS_list[j].particles_in_los_total *
          sizeof(struct xpart));
    } else {
      LOS_particles = (struct part *)malloc(
          LOS_list[j].particles_in_los_local *
          sizeof(struct part));
      LOS_particles_xparts = (struct xpart *)malloc(
          LOS_list[j].particles_in_los_local *
          sizeof(struct xpart));
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
              LOS_particles[count] = p[i];
              LOS_particles_xparts[count] = xp[i];

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
      MPI_Gatherv(MPI_IN_PLACE, 0,
                  xpart_mpi_type, LOS_particles_xparts, LOS_counts, LOS_disps,
                  xpart_mpi_type, 0, MPI_COMM_WORLD);
    } else {
      MPI_Gatherv(LOS_particles, LOS_list[j].particles_in_los_local,
                  part_mpi_type, LOS_particles, LOS_counts, LOS_disps,
                  part_mpi_type, 0, MPI_COMM_WORLD);
      MPI_Gatherv(LOS_particles_xparts, LOS_list[j].particles_in_los_local,
                  xpart_mpi_type, LOS_particles_xparts, LOS_counts, LOS_disps,
                  xpart_mpi_type, 0, MPI_COMM_WORLD);
    }
#endif
    
    /* Write particles to file. */
    if (e->nodeID == 0) {
      /* Create HDF5 group for this LOS */
      sprintf(groupName, "/LOS_%04i", j);
      h_grp = H5Gcreate(h_file, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating LOS HDF5 group\n");

      /* Record this LOS attributes */
      io_write_attribute(h_grp, "NumParts", INT, &LOS_list[j].particles_in_los_total, 1);
      io_write_attribute(h_grp, "Xaxis", INT, &LOS_list[j].xaxis, 1);
      io_write_attribute(h_grp, "Yaxis", INT, &LOS_list[j].yaxis, 1);
      io_write_attribute(h_grp, "Zaxis", INT, &LOS_list[j].zaxis, 1);
      io_write_attribute(h_grp, "Xpos", FLOAT, &LOS_list[j].Xpos, 1);
      io_write_attribute(h_grp, "Ypos", FLOAT, &LOS_list[j].Ypos, 1);

      /* Write the data for this LOS */
      write_los_hdf5_datasets(h_grp, j, LOS_list[j].particles_in_los_total, LOS_particles, e,
        LOS_particles_xparts);

      /* Close HDF5 group */
      H5Gclose(h_grp);
    }
  } /* End of loop over each LOS */

  if (e->nodeID == 0) {
    /* Close HDF5 file */
    H5Fclose(h_file);
  }

  //#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) print_los_info(LOS_list, LOS_params);
  //#endif
 
  /* Up the count. */
  e->los_output_count++;

  /* How long did we take? */
  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void write_los_hdf5_datasets(hid_t grp, int j, int N, const struct part* parts,
        struct engine* e, const struct xpart* xparts) {

  /* What kind of run are we working with? */
  const int with_cosmology = e->policy & engine_policy_cosmology;
  //const int with_cooling = e->policy & engine_policy_cooling;
  //const int with_temperature = e->policy & engine_policy_temperature;
  const int with_fof = e->policy & engine_policy_fof;
#ifdef HAVE_VELOCIRAPTOR
  const int with_stf = (e->policy & engine_policy_structure_finding) &&
                       (e->s->gpart_group_data != NULL);
#else
  const int with_stf = 0;
#endif

  int num_fields = 0;
  struct io_props list[100];

  /* Find the output fields */
  hydro_write_particles(parts, xparts, list, &num_fields);
  num_fields += chemistry_write_particles(parts, list + num_fields);
  //if (with_cooling || with_temperature) {
  //  num_fields += cooling_write_particles(
  //      parts, xparts, list + num_fields, e->cooling_func);
  //}
  if (with_fof) {
    num_fields += fof_write_parts(parts, xparts, list + num_fields);
  }
  if (with_stf) {
    num_fields +=
        velociraptor_write_parts(parts, xparts, list + num_fields);
  }
  num_fields += tracers_write_particles(
      parts, xparts, list + num_fields, with_cosmology);
  num_fields += star_formation_write_particles(parts, xparts,
                                                list + num_fields);

  /* Loop over and write each output field */
  for (int i = 0; i < num_fields; i++) {
    write_los_hdf5_dataset(list[i], N, j, e, grp);
  }
}

/**
 * @brief Writes dataset for a given part attribute.
 *
 * @param p io_props dataset for this attribute.
 * @param N number of parts in this line of sight.
 * @param j Line of sight ID.
 * @param e The engine.
 * @param grp HDF5 group to write to.
 */
void write_los_hdf5_dataset(const struct io_props p, int N, int j, struct engine* e,
        hid_t grp) {

  hid_t dataset_id, dataspace_id;
  herr_t status;
  char att_name[200];

  int rank;
  hsize_t shape[2];
  if (p.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = p.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
  }

  /* Allocate temporary buffer */
  const size_t num_elements = N * p.dimension;
  const size_t typeSize = io_sizeof_type(p.type);
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy particle data to temp buffer */
  io_copy_temp_buffer(temp, e, p, N, e->internal_units, e->snapshot_units);

  /* Prepare dataset */
  dataspace_id = H5Screate_simple(rank, shape, NULL);

  /* Write dataset */
  sprintf(att_name, "/LOS_%04i/%s", j, p.name);
  dataset_id = H5Dcreate(grp, att_name, io_hdf5_type(p.type), dataspace_id,    
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite(dataset_id, io_hdf5_type(p.type), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp);
  if (status < 0) error("Error while writing data array '%s'.", p.name);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, e->snapshot_units, p.units,
                              p.scale_factor_exponent);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, p.units);
  io_write_attribute_f(dataset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(dataset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(dataset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(dataset_id, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(dataset_id, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(dataset_id, "h-scale exponent", 0.f);
  io_write_attribute_f(dataset_id, "a-scale exponent", p.scale_factor_exponent);
  io_write_attribute_s(dataset_id, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double factor =
      units_cgs_conversion_factor(e->snapshot_units, p.units);
  io_write_attribute_d(
      dataset_id,
      "Conversion factor to CGS (not including cosmological corrections)",
      factor);
  io_write_attribute_d(
      dataset_id,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor * pow(e->cosmology->a, p.scale_factor_exponent));

#ifdef SWIFT_DEBUG_CHECKS
  if (strlen(p.description) == 0)
    error("Invalid (empty) description of the field '%s'", p.name);
#endif

  /* Write the full description */
  io_write_attribute_s(dataset_id, "Description", p.description);

  /* Free and close everything */
  swift_free("writebuff", temp);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

/**
 * @brief Writes HDF5 headers and information groups for this line of sight.
 *
 * @param h_file HDF5 file reference.
 * @param e The engine.
 * @param LOS_params The line of sight params.
 */
void write_hdf5_header(hid_t h_file, const struct engine *e, const struct los_props *LOS_params) {
  /* Open header to write simulation properties */
  hid_t h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time =
      units_conversion_factor(e->internal_units, e->snapshot_units, UNIT_CONV_TIME);
  const double factor_length = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
  const int dimension = (int)hydro_dimension;
  io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);
  io_write_attribute(h_grp, "Redshift", DOUBLE, &e->cosmology->z, 1);
  io_write_attribute(h_grp, "Scale-factor", DOUBLE, &e->cosmology->a, 1);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);

  /* Store the time at which the snapshot was written */
  time_t tm = time(NULL);
  struct tm* timeinfo = localtime(&tm);
  char snapshot_date[64];
  strftime(snapshot_date, 64, "%T %F %Z", timeinfo);
  io_write_attribute_s(h_grp, "Snapshot date", snapshot_date);

  /* Close group */
  H5Gclose(h_grp);

  /* Print the code version */
  io_write_code_description(h_file);

  /* Print the run's policy */
  io_write_engine_policy(h_file, e);

  /* Print the physical constants */
  phys_const_print_snapshot(h_file, e->physical_constants);

  /* Print the SPH parameters */
  if (e->policy & engine_policy_hydro) {
    h_grp = H5Gcreate(h_file, "/HydroScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating SPH group");
    hydro_props_print_snapshot(h_grp, e->hydro_properties);
    hydro_write_flavour(h_grp);
    H5Gclose(h_grp);
  }

  /* Print the subgrid parameters */
  h_grp = H5Gcreate(h_file, "/SubgridScheme", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating subgrid group");
  hid_t h_grp_columns =
      H5Gcreate(h_grp, "NamedColumns", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp_columns < 0) error("Error while creating named columns group");
  entropy_floor_write_flavour(h_grp);
  //cooling_write_flavour(h_grp, h_grp_columns, e->cooling_func);
  H5Gclose(h_grp_columns);
  H5Gclose(h_grp);

  /* Print the gravity parameters */
  if (e->policy & engine_policy_self_gravity) {
    h_grp = H5Gcreate(h_file, "/GravityScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating gravity group");
    gravity_props_print_snapshot(h_grp, e->gravity_properties);
    H5Gclose(h_grp);
  }

  /* Print the cosmological model */
  h_grp =
      H5Gcreate(h_file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cosmology group");
  if (e->policy & engine_policy_cosmology)
    io_write_attribute_i(h_grp, "Cosmological run", 1);
  else
    io_write_attribute_i(h_grp, "Cosmological run", 0);
  cosmology_write_model(h_grp, e->cosmology);
  H5Gclose(h_grp);

  /* Print the runtime parameters */
  h_grp =
      H5Gcreate(h_file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, 1);
  H5Gclose(h_grp);

  /* Print the runtime unused parameters */
  h_grp = H5Gcreate(h_file, "/UnusedParameters", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, 0);
  H5Gclose(h_grp);

  /* Print the system of Units used in the spashot */
  io_write_unit_system(h_file, e->snapshot_units, "Units");

  /* Print the system of Units used internally */
  io_write_unit_system(h_file, e->internal_units, "InternalCodeUnits");

  /* Print the LOS properties */
  h_grp = H5Gcreate(h_file, "/LineOfSightParameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating LOS group");

  /* Record this LOS attributes */
  io_write_attribute(h_grp, "NumAlongXY", INT, &LOS_params->num_along_xy, 1);
  io_write_attribute(h_grp, "NumAlongYZ", INT, &LOS_params->num_along_yz, 1);
  io_write_attribute(h_grp, "NumAlongXZ", INT, &LOS_params->num_along_xz, 1);
  io_write_attribute(h_grp, "NumLineOfSight", INT, &LOS_params->num_tot, 1);
  io_write_attribute(h_grp, "Minx", FLOAT, &LOS_params->xmin, 1);
  io_write_attribute(h_grp, "Maxx", FLOAT, &LOS_params->xmax, 1);
  io_write_attribute(h_grp, "Miny", FLOAT, &LOS_params->ymin, 1);
  io_write_attribute(h_grp, "Maxy", FLOAT, &LOS_params->ymax, 1);
  io_write_attribute(h_grp, "Minz", FLOAT, &LOS_params->zmin, 1);
  io_write_attribute(h_grp, "Maxz", FLOAT, &LOS_params->zmax, 1);
  H5Gclose(h_grp);
}
