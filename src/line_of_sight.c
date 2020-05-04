/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "cooling_io.h"
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
                             const struct los_props *params,
                             const int periodic, const double dim[3]) {

  /* Keep track of number of sightlines. */
  int count = 0;

  double Xpos, Ypos;

  /* Sightlines in XY plane. */
  for (int i = 0; i < params->num_along_xy; i++) {
    Xpos = ((float)rand() / (float)(RAND_MAX) * (params->xmax - params->xmin)) +
        params->xmin;
    Ypos = ((float)rand() / (float)(RAND_MAX) * (params->ymax - params->ymin)) +
        params->ymin;
    create_line_of_sight(Xpos, Ypos, 0, 1, 2, periodic, dim, &Los[count]);
    count += 1;
  }

  /* Sightlines in YZ plane. */
  for (int i = 0; i < params->num_along_yz; i++) {
    Xpos = ((float)rand() / (float)(RAND_MAX) * (params->ymax - params->ymin)) +
        params->ymin;
    Ypos = ((float)rand() / (float)(RAND_MAX) * (params->zmax - params->zmin)) +
        params->zmin;
    create_line_of_sight(Xpos, Ypos, 1, 2, 0, periodic, dim, &Los[count]);
    count += 1;
  }

  /* Sightlines in XZ plane. */
  for (int i = 0; i < params->num_along_xz; i++) {
    Xpos = ((float)rand() / (float)(RAND_MAX) * (params->xmax - params->xmin)) +
        params->xmin;
    Ypos = ((float)rand() / (float)(RAND_MAX) * (params->zmax - params->zmin)) +
        params->zmin;
    create_line_of_sight(Xpos, Ypos, 0, 2, 1, periodic, dim, &Los[count]);
    count += 1;
  }

  /* Make sure we made the correct ammount */
  if (count != params->num_tot) error("Could not make the right number of sightlines");
}

void create_line_of_sight(const double Xpos, const double Ypos, 
        const int xaxis, const int yaxis, const int zaxis,
        const int periodic, const double dim[3], struct line_of_sight *los) {
  los->Xpos = Xpos;
  los->Ypos = Ypos;
  los->particles_in_los_local = 0;
  los->particles_in_los_total = 0;
  los->xaxis = xaxis;
  los->yaxis = yaxis;
  los->zaxis = zaxis;
  los->periodic = periodic;
  los->dim[0] = dim[0];
  los->dim[1] = dim[1];
  los->dim[2] = dim[2];
}

/**
 * @brief Print line_of_sight information.
 *
 * @param Los Structure to print.
 */
void print_los_info(const struct line_of_sight *Los, const int i) {

  printf("[LOS %i] Xpos:%g Ypos:%g particles_in_los_total:%i\n", i,
        Los[i].Xpos, Los[i].Ypos, Los[i].particles_in_los_total);
  fflush(stdout);
}

/**
 * @brief Loop over each part to see which ones intersect the LOS.
 * 
 * @param map_data The parts.
 * @param count The number of parts.
 * @param extra_data The line_of_sight structure for this LOS.
 */
void los_first_loop_mapper(void *restrict map_data, int count,
                           void *restrict extra_data) {

  int los_particle_count = 0;
  double dx, dy, r2, hsml;
  struct line_of_sight *restrict LOS_list = (struct line_of_sight *)extra_data;
  struct part *restrict p = (struct part *)map_data;

  /* Loop over each part to find those in LOS. */
  for (size_t i = 0; i < count; i++) {
    if (p[i].gpart->type == swift_type_gas) {
      /* Distance from this part to LOS along x dim. */
      dx = p[i].x[LOS_list->xaxis] - LOS_list->Xpos;

      /* Periodic wrap. */
      if (LOS_list->periodic) dx = nearest(dx, LOS_list->dim[LOS_list->xaxis]);

      /* Smoothing length of this part. */
      hsml = p[i].h * kernel_gamma;

      /* Does this particle fall into our LOS? */
      if (dx <= hsml) {
        /* Distance from this part to LOS along y dim. */
        dy = p[i].x[LOS_list->yaxis] - LOS_list->Ypos;
      
        /* Periodic wrap. */
        if (LOS_list->periodic) dy = nearest(dy, LOS_list->dim[LOS_list->yaxis]);

        /* Does this part still fall into our LOS? */
        if (dy <= hsml) {
          /* 2D distance to LOS. */
          r2 = dx * dx + dy * dy;

          if (r2 <= hsml * hsml) {

            /* We've found one. */
            los_particle_count++;
          }
        }
      }
    }
  } /* End of loop over all gas particles */

  atomic_add(&LOS_list->particles_in_los_local, los_particle_count);
}

/**
 * @brief Main work function for computing line of sights.
 * 
 * 1) Construct N random line of sight positions.
 * 2) Loop over each line of sight.
 * -  2.1) Loop over each part to see which fall within this LOS.
 * -  2.2) Use this count to construct a LOS parts array.
 * -  2.3) Loop over each part and extract those in LOS to new array.
 * -  2.4) Save LOS parts to HDF5 file.
 * 
 * @param e The engine.
 */
void do_line_of_sight(struct engine *e) {

  /* Start counting. */
  const ticks tic = getticks();

  const struct space *s = e->s;
  struct part *p = s->parts;
  const struct xpart *xp = s->xparts;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const size_t nr_parts = s->nr_parts;
  const struct los_props *LOS_params = e->los_properties;
  const int verbose = e->verbose;

  double dx, dy, r2, hsml;

  /* HDF5 vars. */
  hid_t h_file = -1, h_grp = -1;
  char fileName[256], groupName[200];

  /* Start by generating the random sightline positions. */
  struct line_of_sight *LOS_list = (struct line_of_sight *)malloc(
      LOS_params->num_tot * sizeof(struct line_of_sight));
  if (e->nodeID == 0) {
      if (verbose) message("Generating random sightlines...");
      generate_line_of_sights(LOS_list, LOS_params, periodic, dim);
      if (verbose) message("Generated %i random sightlines.", LOS_params->num_tot);
  }
#ifdef WITH_MPI
  MPI_Bcast(LOS_list, LOS_params->num_tot * sizeof(struct line_of_sight),
            MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Node 0 creates the HDF5 file. */
  if (e->nodeID == 0) {
    sprintf(fileName, "%s_%04i.hdf5", LOS_params->basename, e->los_output_count);
    if (verbose) message("Creating LOS file: %s", fileName);
    h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Error while opening file '%s'.", fileName);

    /* Write header */
    write_hdf5_header(h_file, e, LOS_params);
  }
#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Main loop over each random LOS. */
  for (int j = 0; j < LOS_params->num_tot; j++) {

    /* First count all parts that intersect with this line of sight */
    threadpool_map(&s->e->threadpool, los_first_loop_mapper, p,
              nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
              &LOS_list[j]);

#ifdef WITH_MPI
    /* Make sure all nodes know how many parts are in this LOS */
    int LOS_counts[e->nr_nodes];
    int LOS_disps[e->nr_nodes];

    /* How many parts does each rank have for this LOS? */
    MPI_Allgather(&LOS_list[j].particles_in_los_local, 1, MPI_INT,
                  &LOS_counts, 1, MPI_INT, MPI_COMM_WORLD);

    for (int k = 0, disp_count = 0; k < e->nr_nodes; k++) {
      /* Total parts in this LOS. */
      LOS_list[j].particles_in_los_total += LOS_counts[k];

      /* Counts and disps for Gatherv. */
      LOS_disps[k] = disp_count;
      disp_count += LOS_counts[k];
    }
#else
    LOS_list[j].particles_in_los_total = LOS_list[j].particles_in_los_local;
#endif

    //#ifdef SWIFT_DEBUG_CHECKS
    if (e->nodeID == 0) print_los_info(LOS_list, j);
    //#endif

    /* Don't work with empty LOS */
    if (LOS_list[j].particles_in_los_total == 0) {
        if (e->nodeID == 0) {
          message("*WARNING* LOS %i is empty", j);
          print_los_info(LOS_list, j);
        }
        continue;
    }

    /* Setup LOS part and xpart structures. */
    struct part *LOS_parts = NULL;
    struct xpart *LOS_xparts = NULL;

    if (e->nodeID == 0) {
      if ((LOS_parts = (struct part *)swift_malloc(
            "los_parts_array", sizeof(struct part) * LOS_list[j].particles_in_los_total)) == NULL)
        error("Failed to allocate LOS part memory.");
      if ((LOS_xparts = (struct xpart *)swift_malloc(
            "los_xparts_array", sizeof(struct xpart) * LOS_list[j].particles_in_los_total)) == NULL)
        error("Failed to allocate LOS xpart memory.");
    } else {
      if ((LOS_parts = (struct part *)swift_malloc(
            "los_parts_array", sizeof(struct part) * LOS_list[j].particles_in_los_local)) == NULL)
        error("Failed to allocate LOS part memory.");
      if ((LOS_xparts = (struct xpart *)swift_malloc(
            "los_xparts_array", sizeof(struct xpart) * LOS_list[j].particles_in_los_local)) == NULL)
        error("Failed to allocate LOS xpart memory.");
    }

    /* Loop over each part again, pulling out those in LOS. */
    int count = 0;
    
    for (size_t i = 0; i < nr_parts; i++) {

      /* Distance from this part to LOS along x dim. */
      dx = p[i].x[LOS_list[j].xaxis] - LOS_list[j].Xpos;

      /* Periodic wrap. */
      if (LOS_list[j].periodic) dx = nearest(dx, LOS_list[j].dim[LOS_list[j].xaxis]);

      /* Smoothing length of this part. */
      hsml = p[i].h * kernel_gamma;

      /* Does this part fall into our LOS? */
      if (dx <= hsml) {
        /* Distance from this part to LOS along y dim. */
        dy = p[i].x[LOS_list[j].yaxis] - LOS_list[j].Ypos;

        /* Periodic wrap. */
        if (LOS_list[j].periodic) dy = nearest(dy, LOS_list[j].dim[LOS_list[j].yaxis]);

        /* Does this part still fall into our LOS? */
        if (dy <= hsml) {
          /* 2D distance to LOS. */
          r2 = dx * dx + dy * dy;

          if (r2 <= hsml * hsml) {

            /* Store part and xpart properties. */
            memcpy(&LOS_parts[count], &p[i], sizeof(struct part));
            memcpy(&LOS_xparts[count], &xp[i], sizeof(struct xpart));

            count++;
          }
        }
      }
    }

//#ifdef SWIFT_DEBUG_CHECKS
    if (count != LOS_list[j].particles_in_los_local) error("LOS counts don't add up");
//#endif

#ifdef WITH_MPI
    /* Collect all parts in this LOS to rank 0. */
    if (e->nodeID == 0) {
      MPI_Gatherv(MPI_IN_PLACE, 0,
                  part_mpi_type, LOS_parts, LOS_counts, LOS_disps,
                  part_mpi_type, 0, MPI_COMM_WORLD);
      MPI_Gatherv(MPI_IN_PLACE, 0,
                  xpart_mpi_type, LOS_xparts, LOS_counts, LOS_disps,
                  xpart_mpi_type, 0, MPI_COMM_WORLD);
    } else {
      MPI_Gatherv(LOS_parts, LOS_list[j].particles_in_los_local,
                  part_mpi_type, LOS_parts, LOS_counts, LOS_disps,
                  part_mpi_type, 0, MPI_COMM_WORLD);
      MPI_Gatherv(LOS_xparts, LOS_list[j].particles_in_los_local,
                  xpart_mpi_type, LOS_xparts, LOS_counts, LOS_disps,
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
      io_write_attribute(h_grp, "Xpos", DOUBLE, &LOS_list[j].Xpos, 1);
      io_write_attribute(h_grp, "Ypos", DOUBLE, &LOS_list[j].Ypos, 1);

      /* Write the data for this LOS */
      write_los_hdf5_datasets(h_grp, j, LOS_list[j].particles_in_los_total, LOS_parts, e,
        LOS_xparts);

      /* Close HDF5 group */
      H5Gclose(h_grp);
    }

    /* Free up some memory */
    swift_free("los_parts_array", LOS_parts);
    swift_free("los_xparts_array", LOS_xparts);

  } /* End of loop over each LOS */

  /* Close HDF5 file */
  if (e->nodeID == 0) H5Fclose(h_file);
 
  /* Up the count. */
  e->los_output_count++;

  /* How long did we take? */
  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Write parts in LOS to HDF5 file.
 *
 * @param grp HDF5 group of this LOS.
 * @param j LOS ID.
 * @param N number of parts in this line of sight.
 * @param parts the list of parts in this LOS.
 * @param e The engine.
 * @param xparts the list of xparts in this LOS.
 */
void write_los_hdf5_datasets(hid_t grp, int j, int N, const struct part* parts,
        struct engine* e, const struct xpart* xparts) {

  /* What kind of run are we working with? */
  struct swift_params* params = e->parameter_file;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int with_fof = e->policy & engine_policy_fof;
#ifdef HAVE_VELOCIRAPTOR
  const int with_stf = (e->policy & engine_policy_structure_finding) &&
                       (e->s->gpart_group_data != NULL);
#else
  const int with_stf = 0;
#endif

  int num_fields = 0;
  struct io_props list[100];

  /* Find all the gas output fields */
  hydro_write_particles(parts, xparts, list, &num_fields);
  num_fields += chemistry_write_particles(parts, list + num_fields);
  if (with_cooling || with_temperature) {
    num_fields += cooling_write_particles(
        parts, xparts, list + num_fields, e->cooling_func);
  }
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

  /* Loop over each output field */
  for (int i = 0; i < num_fields; i++) {

    /* Did the user cancel this field? */
    char field[PARSER_MAX_LINE_SIZE];
    sprintf(field, "SelectOutputLOS:%.*s", FIELD_BUFFER_SIZE,
            list[i].name);
    int should_write = parser_get_opt_param_int(params, field, 1);  
    
    /* Write (if selected) */
    if (should_write) write_los_hdf5_dataset(list[i], N, j, e, grp);
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
void write_los_hdf5_dataset(const struct io_props props, int N, int j, struct engine* e,
        hid_t grp) {

  /* Create data space */
  const hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  int rank = 0;
  hsize_t shape[2];
  hsize_t chunk_shape[2];
  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = props.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = 0;
  }

  /* Make sure the chunks are not larger than the dataset */
  if (chunk_shape[0] > N) chunk_shape[0] = N;

  /* Change shape of data space */
  hid_t h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
  if (h_err < 0)
    error("Error while changing data space shape for field '%s'.", props.name);

  /* Dataset properties */
  const hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Set chunk size */
  h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
  if (h_err < 0)
    error("Error while setting chunk size (%llu, %llu) for field '%s'.",
          chunk_shape[0], chunk_shape[1], props.name);

  /* Impose check-sum to verify data corruption */
  h_err = H5Pset_fletcher32(h_prop);
  if (h_err < 0)
    error("Error while setting checksum options for field '%s'.", props.name);

  /* Impose data compression */
  if (e->snapshot_compression > 0) {
    h_err = H5Pset_shuffle(h_prop);
    if (h_err < 0)
      error("Error while setting shuffling options for field '%s'.",
            props.name);

    h_err = H5Pset_deflate(h_prop, e->snapshot_compression);
    if (h_err < 0)
      error("Error while setting compression options for field '%s'.",
            props.name);
  }

  /* Allocate temporary buffer */
  const size_t num_elements = N * props.dimension;
  const size_t typeSize = io_sizeof_type(props.type);
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy particle data to temp buffer */
  io_copy_temp_buffer(temp, e, props, N, e->internal_units, e->snapshot_units);

  /* Create dataset */
  char att_name[200];
  sprintf(att_name, "/LOS_%04i/%s", j, props.name);
  const hid_t h_data = H5Dcreate(grp, att_name, io_hdf5_type(props.type),
                                 h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", props.name);

  /* Write dataset */
  herr_t status = H5Dwrite(h_data, io_hdf5_type(props.type), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp);
  if (status < 0) error("Error while writing data array '%s'.", props.name);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, e->snapshot_units, props.units,
                              props.scale_factor_exponent);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, props.units);
  io_write_attribute_f(h_data, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(h_data, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(h_data, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(h_data, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(h_data, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(h_data, "h-scale exponent", 0.f);
  io_write_attribute_f(h_data, "a-scale exponent", props.scale_factor_exponent);
  io_write_attribute_s(h_data, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double factor =
      units_cgs_conversion_factor(e->snapshot_units, props.units);
  io_write_attribute_d(
      h_data,
      "Conversion factor to CGS (not including cosmological corrections)",
      factor);
  io_write_attribute_d(
      h_data,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor * pow(e->cosmology->a, props.scale_factor_exponent));

#ifdef SWIFT_DEBUG_CHECKS
  if (strlen(props.description) == 0)
    error("Invalid (empty) description of the field '%s'", props.name);
#endif

  /* Write the full description */
  io_write_attribute_s(h_data, "Description", props.description);

  /* Free and close everything */
  swift_free("writebuff", temp);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);
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
  cooling_write_flavour(h_grp, h_grp_columns, e->cooling_func);
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
  io_write_attribute(h_grp, "Minx", DOUBLE, &LOS_params->xmin, 1);
  io_write_attribute(h_grp, "Maxx", DOUBLE, &LOS_params->xmax, 1);
  io_write_attribute(h_grp, "Miny", DOUBLE, &LOS_params->ymin, 1);
  io_write_attribute(h_grp, "Maxy", DOUBLE, &LOS_params->ymax, 1);
  io_write_attribute(h_grp, "Minz", DOUBLE, &LOS_params->zmin, 1);
  io_write_attribute(h_grp, "Maxz", DOUBLE, &LOS_params->zmax, 1);
  H5Gclose(h_grp);
}

/**
 * @brief Write a los_props struct to the given FILE as a stream of bytes.
 *
 * @param internal_los the struct
 * @param stream the file stream
 */
void los_struct_dump(const struct los_props *internal_los,
                            FILE *stream) {
  restart_write_blocks((void *)internal_los, sizeof(struct los_props), 1,
                       stream, "losparams", "los params");
}

/**
 * @brief Restore a los_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param internal_los the struct
 * @param stream the file stream
 */
void los_struct_restore(const struct los_props *internal_los,
                               FILE *stream) {
  restart_read_blocks((void *)internal_los, sizeof(struct los_props), 1,
                      stream, NULL, "los params");
}

