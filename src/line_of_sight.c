/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "atomic.h"
#include "engine.h"
#include "hydro_io.h"
#include "io_properties.h"
#include "kernel_hydro.h"
#include "line_of_sight.h"
#include "periodic.h"
#include "version.h"

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Will the line of sight intersect a given cell?
 *
 * Also return 0 if the cell is empty.
 *
 * @param c The top level cell.
 * @param los The line of sight structure.
 */
static INLINE int does_los_intersect(const struct cell *c,
                                     const struct line_of_sight *los) {

  /* Empty cell? */
  if (c->hydro.count == 0) return 0;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.h_max <= 0.) error("Invalid h_max for does_los_intersect");
#endif

  /* Distance from LOS to left and bottom cell edge. */
  const double cx_min = c->loc[los->xaxis];
  const double cy_min = c->loc[los->yaxis];

  /* Distance from LOS to right and top cell edge. */
  const double cx_max = c->loc[los->xaxis] + c->width[los->xaxis];
  const double cy_max = c->loc[los->yaxis] + c->width[los->yaxis];

  /* Maximum smoothing length of a part in this top level cell. */
  const double hsml = c->hydro.h_max * kernel_gamma;
  const double hsml2 = hsml * hsml;

  double dx, dy;

  if (los->periodic) {
    dx = min(fabs(nearest(cx_min - los->Xpos, los->dim[los->xaxis])),
             fabs(nearest(cx_max - los->Xpos, los->dim[los->xaxis])));
    dy = min(fabs(nearest(cy_min - los->Ypos, los->dim[los->yaxis])),
             fabs(nearest(cy_max - los->Ypos, los->dim[los->yaxis])));
  } else {
    dx = fabs(cx_max - los->Xpos);
    dy = fabs(cy_max - los->Ypos);
  }

  /* Is sightline directly within this top level cell? */
  if (dx < (1.01 * c->width[los->xaxis]) / 2. &&
      dy < (1.01 * c->width[los->yaxis]) / 2.) {
    return 1;
    /* Could a part from this top level cell smooth into the sightline? */
  } else if (dx * dx + dy * dy < hsml2) {
    return 1;
    /* Don't need to work with this top level cell. */
  } else {
    return 0;
  }
}

/**
 * @brief Reads the LOS properties from the param file.
 *
 * @param dim Space dimensions.
 * @param los_params Sightline parameters to save into.
 * @param params Swift params to read from.
 */
void los_init(const double dim[3], struct los_props *los_params,
              struct swift_params *params) {

  /* How many line of sights in each plane. */
  los_params->num_along_x =
      parser_get_param_int(params, "LineOfSight:num_along_x");
  los_params->num_along_y =
      parser_get_param_int(params, "LineOfSight:num_along_y");
  los_params->num_along_z =
      parser_get_param_int(params, "LineOfSight:num_along_z");

  /* Min/max range across x,y and z (simulation axes) where random
   * LOS's are allowed. */
  los_params->allowed_losrange_x[0] = 0.;
  los_params->allowed_losrange_x[1] = dim[0];
  parser_get_opt_param_double_array(params, "LineOfSight:allowed_los_range_x",
                                    2, los_params->allowed_losrange_x);
  los_params->allowed_losrange_y[0] = 0.;
  los_params->allowed_losrange_y[1] = dim[1];
  parser_get_opt_param_double_array(params, "LineOfSight:allowed_los_range_y",
                                    2, los_params->allowed_losrange_y);
  los_params->allowed_losrange_z[0] = 0.;
  los_params->allowed_losrange_z[1] = dim[2];
  parser_get_opt_param_double_array(params, "LineOfSight:allowed_los_range_z",
                                    2, los_params->allowed_losrange_z);

  /* Compute total number of sightlines. */
  los_params->num_tot = los_params->num_along_z + los_params->num_along_x +
                        los_params->num_along_y;

  /* Where are we saving them? */
  parser_get_param_string(params, "LineOfSight:basename", los_params->basename);

  /* Min/max range allowed when sightline is shooting down x,y and z
   * (simulation axes). */
  los_params->range_when_shooting_down_axis[0][0] = 0.;
  los_params->range_when_shooting_down_axis[0][1] = dim[0];
  parser_get_opt_param_double_array(
      params, "LineOfSight:range_when_shooting_down_x", 2,
      los_params->range_when_shooting_down_axis[0]);
  los_params->range_when_shooting_down_axis[1][0] = 0.;
  los_params->range_when_shooting_down_axis[1][1] = dim[1];
  parser_get_opt_param_double_array(
      params, "LineOfSight:range_when_shooting_down_y", 2,
      los_params->range_when_shooting_down_axis[1]);
  los_params->range_when_shooting_down_axis[2][0] = 0.;
  los_params->range_when_shooting_down_axis[2][1] = dim[2];
  parser_get_opt_param_double_array(
      params, "LineOfSight:range_when_shooting_down_z", 2,
      los_params->range_when_shooting_down_axis[2]);
}

/**
 *  @brief Create a #line_of_sight object from its attributes
 */
void create_sightline(const double Xpos, const double Ypos,
                      enum los_direction xaxis, enum los_direction yaxis,
                      enum los_direction zaxis, const int periodic,
                      const double dim[3], struct line_of_sight *los,
                      const double range_when_shooting_down_axis[2]) {
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
  los->num_intersecting_top_level_cells = 0;
  los->range_when_shooting_down_axis[0] = range_when_shooting_down_axis[0];
  los->range_when_shooting_down_axis[1] = range_when_shooting_down_axis[1];
}

/**
 * @brief Generates random sightline positions.
 *
 * Independent sightlines are made for the XY, YZ and XZ planes.
 *
 * @param Los Structure to store sightlines.
 * @param params Sightline parameters.
 * @param periodic Is this calculation using periodic BCs.
 * @param dim The dimension of the volume along the three axis.
 */
void generate_sightlines(struct line_of_sight *Los,
                         const struct los_props *params, const int periodic,
                         const double dim[3]) {

  /* Keep track of number of sightlines. */
  int count = 0;

  /* Sightlines in XY plane, shoots down Z. */
  for (int i = 0; i < params->num_along_z; i++) {
    double Xpos =
        ((float)rand() / (float)(RAND_MAX) *
         (params->allowed_losrange_x[1] - params->allowed_losrange_x[0])) +
        params->allowed_losrange_x[0];
    double Ypos =
        ((float)rand() / (float)(RAND_MAX) *
         (params->allowed_losrange_y[1] - params->allowed_losrange_y[0])) +
        params->allowed_losrange_y[0];

    create_sightline(Xpos, Ypos, simulation_x_axis, simulation_y_axis,
                     simulation_z_axis, periodic, dim, &Los[count],
                     params->range_when_shooting_down_axis[simulation_z_axis]);
    count++;
  }

  /* Sightlines in YZ plane, shoots down X. */
  for (int i = 0; i < params->num_along_x; i++) {
    double Xpos =
        ((float)rand() / (float)(RAND_MAX) *
         (params->allowed_losrange_y[1] - params->allowed_losrange_y[0])) +
        params->allowed_losrange_y[0];
    double Ypos =
        ((float)rand() / (float)(RAND_MAX) *
         (params->allowed_losrange_z[1] - params->allowed_losrange_z[0])) +
        params->allowed_losrange_z[0];

    create_sightline(Xpos, Ypos, simulation_y_axis, simulation_z_axis,
                     simulation_x_axis, periodic, dim, &Los[count],
                     params->range_when_shooting_down_axis[simulation_x_axis]);
    count++;
  }

  /* Sightlines in XZ plane, shoots down Y. */
  for (int i = 0; i < params->num_along_y; i++) {
    double Xpos =
        ((float)rand() / (float)(RAND_MAX) *
         (params->allowed_losrange_x[1] - params->allowed_losrange_x[0])) +
        params->allowed_losrange_x[0];
    double Ypos =
        ((float)rand() / (float)(RAND_MAX) *
         (params->allowed_losrange_z[1] - params->allowed_losrange_z[0])) +
        params->allowed_losrange_z[0];

    create_sightline(Xpos, Ypos, simulation_x_axis, simulation_z_axis,
                     simulation_y_axis, periodic, dim, &Los[count],
                     params->range_when_shooting_down_axis[simulation_y_axis]);
    count++;
  }

  /* Make sure we made the correct ammount */
  if (count != params->num_tot)
    error("Could not make the right number of sightlines");
}

/**
 * @brief Print #line_of_sight information.
 *
 * @param Los Structure to print.
 * @param i The index of the #line_of_sight to dump.
 */
void print_los_info(const struct line_of_sight *Los, const int i) {

  message(
      "[LOS %i] Xpos:%g Ypos:%g parts_in_los:%i "
      "num_intersecting_top_level_cells:%i",
      i, Los[i].Xpos, Los[i].Ypos, Los[i].particles_in_los_total,
      Los[i].num_intersecting_top_level_cells);
}

/**
 * @brief Writes dataset for a given part attribute.
 *
 * @param props dataset for this attribute.
 * @param N number of parts in this line of sight.
 * @param j Line of sight ID.
 * @param e The engine.
 * @param grp HDF5 group to write to.
 */
void write_los_hdf5_dataset(const struct io_props props, const size_t N,
                            const int j, const struct engine *e, hid_t grp) {

  /* Create data space */
  const hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  /* Decide what chunk size to use based on compression */
  int log2_chunk_size = e->snapshot_compression > 0 ? 12 : 18;

  int rank = 0;
  hsize_t shape[2];
  hsize_t chunk_shape[2];
  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    chunk_shape[0] = 1 << log2_chunk_size;
    chunk_shape[1] = props.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << log2_chunk_size;
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
          (unsigned long long)chunk_shape[0],
          (unsigned long long)chunk_shape[1], props.name);

  /* Impose check-sum to verify data corruption */
  h_err = H5Pset_fletcher32(h_prop);
  if (h_err < 0)
    error("Error while setting checksum options for field '%s'.", props.name);

  /* Impose data compression */
  char comp_buffer[32] = "None";
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
  void *temp = NULL;
  if (swift_memalign("writebuff", (void **)&temp, IO_BUFFER_ALIGNMENT,
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
  herr_t status = H5Dwrite(h_data, io_hdf5_type(props.type), H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, temp);
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
  io_write_attribute_s(h_data, "Lossy compression filter", comp_buffer);

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
 * @brief Write parts in LOS to HDF5 file.
 *
 * @param grp HDF5 group of this LOS.
 * @param j LOS ID.
 * @param N number of parts in this line of sight.
 * @param parts the list of parts in this LOS.
 * @param e The engine.
 * @param xparts the list of xparts in this LOS.
 */
void write_los_hdf5_datasets(hid_t grp, const int j, const size_t N,
                             const struct part *parts, const struct engine *e,
                             const struct xpart *xparts) {

  /* What kind of run are we working with? */
  struct swift_params *params = e->parameter_file;
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
  const int with_rt = e->policy & engine_policy_rt;

  int num_fields = 0;
  struct io_props list[100];

  /* Find all the gas output fields */
  io_select_hydro_fields(parts, xparts, with_cosmology, with_cooling,
                         with_temperature, with_fof, with_stf, with_rt, e,
                         &num_fields, list);

  /* Loop over each output field */
  for (int i = 0; i < num_fields; i++) {

    /* Did the user cancel this field? */
    char field[PARSER_MAX_LINE_SIZE];
    sprintf(field, "SelectOutputLOS:%.*s", FIELD_BUFFER_SIZE, list[i].name);
    int should_write = parser_get_opt_param_int(params, field, 1);

    /* Write (if selected) */
    if (should_write) write_los_hdf5_dataset(list[i], N, j, e, grp);
  }
}

/**
 * @brief Writes HDF5 headers and information groups for this line of sight.
 *
 * @param h_file HDF5 file reference.
 * @param e The engine.
 * @param LOS_params The line of sight params.
 * @param total_num_parts_in_los The total number of particles in all the LoS.
 */
void write_hdf5_header(hid_t h_file, const struct engine *e,
                       const struct los_props *LOS_params,
                       const size_t total_num_parts_in_los) {

  /* Open header to write simulation properties */
  hid_t h_grp =
      H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_TIME);
  const double factor_length = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute_d(h_grp, "Time", dblTime);
  io_write_attribute_d(h_grp, "Dimension", (int)hydro_dimension);
  io_write_attribute_d(h_grp, "Redshift", e->cosmology->z);
  io_write_attribute_d(h_grp, "Scale-factor", e->cosmology->a);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);
  io_write_attribute_s(h_grp, "System", hostname());
  io_write_attribute(h_grp, "Shift", DOUBLE, e->s->initial_shift, 3);

  /* Write out the particle types */
  io_write_part_type_names(h_grp);

  /* Write out the time-base */
  if (e->policy & engine_policy_cosmology) {
    io_write_attribute_d(h_grp, "TimeBase_dloga", e->time_base);
    const double delta_t = cosmology_get_timebase(e->cosmology, e->ti_current);
    io_write_attribute_d(h_grp, "TimeBase_dt", delta_t);
  } else {
    io_write_attribute_d(h_grp, "TimeBase_dloga", 0);
    io_write_attribute_d(h_grp, "TimeBase_dt", e->time_base);
  }

  /* Store the time at which the snapshot was written */
  time_t tm = time(NULL);
  struct tm *timeinfo = localtime(&tm);
  char snapshot_date[64];
  strftime(snapshot_date, 64, "%T %F %Z", timeinfo);
  io_write_attribute_s(h_grp, "SnapshotDate", snapshot_date);

  /* GADGET-2 legacy values */
  /* Number of particles of each type */
  long long N_total[swift_type_count] = {0};
  N_total[0] = total_num_parts_in_los;
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};
  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);
  }
  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, N_total,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total", UINT, numParticles,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total_HighWord", UINT,
                     numParticlesHighWord, swift_type_count);
  io_write_attribute(h_grp, "TotalNumberOfParticles", LONGLONG, N_total,
                     swift_type_count);
  double MassTable[swift_type_count] = {0};
  io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
  io_write_attribute(h_grp, "InitialMassTable", DOUBLE,
                     e->s->initial_mean_mass_particles, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                     swift_type_count);
  io_write_attribute_i(h_grp, "NumFilesPerSnapshot", 1);
  io_write_attribute_i(h_grp, "ThisFile", 0);
  io_write_attribute_s(h_grp, "SelectOutput", "Default");
  io_write_attribute_i(h_grp, "Virtual", 0);
  const int to_write[swift_type_count] = {1}; /* We can only have gas */
  io_write_attribute(h_grp, "CanHaveTypes", INT, to_write, swift_type_count);
  io_write_attribute_s(h_grp, "OutputType", "LineOfSight");

  /* Close group */
  H5Gclose(h_grp);

  /* Copy metadata from ICs to the file */
  ic_info_write_hdf5(e->ics_metadata, h_file);

  /* Write all the meta-data */
  io_write_meta_data(h_file, e, e->internal_units, e->snapshot_units,
                     /*fof=*/0);

  /* Print the LOS properties */
  h_grp = H5Gcreate(h_file, "/LineOfSightParameters", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating LOS group");

  /* Record this LOS attributes */
  const int num_los_per_axis[3] = {LOS_params->num_along_x,
                                   LOS_params->num_along_y,
                                   LOS_params->num_along_z};
  io_write_attribute(h_grp, "NumLineOfSight_PerAxis", INT, num_los_per_axis, 3);
  io_write_attribute(h_grp, "NumLineOfSight_Total", INT, &LOS_params->num_tot,
                     1);
  io_write_attribute(h_grp, "AllowedLOSRangeX", DOUBLE,
                     LOS_params->allowed_losrange_x, 2);
  io_write_attribute(h_grp, "AllowedLOSRangeY", DOUBLE,
                     LOS_params->allowed_losrange_y, 2);
  io_write_attribute(h_grp, "AllowedLOSRangeZ", DOUBLE,
                     LOS_params->allowed_losrange_z, 2);
  H5Gclose(h_grp);
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

  struct line_of_sight *LOS_list = (struct line_of_sight *)extra_data;
  const struct part *parts = (struct part *)map_data;

  size_t los_particle_count = 0;

  /* Loop over each part to find those in LOS. */
  for (int i = 0; i < count; i++) {

    /* Don't consider inhibited parts. */
    if (parts[i].time_bin == time_bin_inhibited) continue;

    /* Don't consider part if outwith allowed z-range. */
    if (parts[i].x[LOS_list->zaxis] <
            LOS_list->range_when_shooting_down_axis[0] ||
        parts[i].x[LOS_list->zaxis] >
            LOS_list->range_when_shooting_down_axis[1])
      continue;

    /* Distance from this part to LOS along x dim. */
    double dx = parts[i].x[LOS_list->xaxis] - LOS_list->Xpos;

    /* Periodic wrap. */
    if (LOS_list->periodic) dx = nearest(dx, LOS_list->dim[LOS_list->xaxis]);

    /* Square. */
    const double dx2 = dx * dx;

    /* Smoothing length of this part. */
    const double hsml = parts[i].h * kernel_gamma;
    const double hsml2 = hsml * hsml;

    /* Does this particle fall into our LOS? */
    if (dx2 < hsml2) {

      /* Distance from this part to LOS along y dim. */
      double dy = parts[i].x[LOS_list->yaxis] - LOS_list->Ypos;

      /* Periodic wrap. */
      if (LOS_list->periodic) dy = nearest(dy, LOS_list->dim[LOS_list->yaxis]);

      /* Square. */
      const double dy2 = dy * dy;

      /* Does this part still fall into our LOS? */
      if (dy2 < hsml2) {

        /* 2D distance to LOS. */
        if (dx2 + dy2 <= hsml2) {

          /* We've found one. */
          los_particle_count++;
        }
      }
    }
  } /* End of loop over all parts */

  atomic_add(&LOS_list->particles_in_los_local, los_particle_count);
}

/**
 * @brief Find all top level cells that a LOS will intersect.
 *
 * This includes the top level cells the LOS directly passes through
 * and the neighbouring top level cells whose parts could smooth into the LOS.
 *
 * @param e The engine.
 * @param los The line_of_sight structure.
 * @param los_cells_top (return) Array indicating whether this cell is
 * intersected.
 * @param cells The array of top-level cells.
 * @param local_cells_with_particles The list of local non-empty top-level
 * cells.
 * @param nr_local_cells_with_particles The number of local non-empty top-level
 * cells.
 */
void find_intersecting_top_level_cells(
    const struct engine *e, struct line_of_sight *los,
    int *restrict los_cells_top, const struct cell *cells,
    const int *restrict local_cells_with_particles,
    const int nr_local_cells_with_particles) {

  /* Keep track of how many top level cells we intersect. */
  int num_intersecting_top_level_cells = 0;

  /* Loop over each top level cell */
  for (int n = 0; n < nr_local_cells_with_particles; n++) {

    /* This top level cell. */
    const struct cell *c = &cells[local_cells_with_particles[n]];

    if (does_los_intersect(c, los)) {
      num_intersecting_top_level_cells++;
      los_cells_top[n] = 1;
    }
  }

#ifdef WITH_MPI
  if (MPI_Allreduce(MPI_IN_PLACE, &num_intersecting_top_level_cells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce num_intersecting_top_level_cells.");
#endif

  /* Store how many top level cells this LOS intersects. */
  los->num_intersecting_top_level_cells = num_intersecting_top_level_cells;
}

/**
 * @brief Main work function for computing line of sights.
 *
 * 1) Construct N random line of sight positions.
 * 2) Loop over each line of sight.
 *  - 2.1) Find which top level cells sightline intersects.
 *  - 2.2) Loop over each part in these top level cells to see which intersect
 * sightline.
 *  - 2.3) Use this count to construct a LOS parts/xparts array.
 *  - 2.4) Loop over each part and extract those in sightline to new array.
 *  - 2.5) Save sightline parts to HDF5 file.
 *
 * @param e The engine.
 */
void do_line_of_sight(struct engine *e) {

  /* Start counting. */
  const ticks tic = getticks();

  const struct space *s = e->s;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const struct los_props *LOS_params = e->los_properties;
  const int verbose = e->verbose;

  /* Start by generating the random sightline positions. */
  struct line_of_sight *LOS_list = (struct line_of_sight *)malloc(
      LOS_params->num_tot * sizeof(struct line_of_sight));

  if (e->nodeID == 0) {
    generate_sightlines(LOS_list, LOS_params, periodic, dim);
    if (verbose)
      message("Generated %i random sightlines.", LOS_params->num_tot);
  }

#ifdef WITH_MPI
  /* Share the list of LoS with all the MPI ranks */
  MPI_Bcast(LOS_list, LOS_params->num_tot * sizeof(struct line_of_sight),
            MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Node 0 creates the HDF5 file. */
  hid_t h_file = -1, h_grp = -1;
  char fileName[256], groupName[200];

  if (e->nodeID == 0) {
    sprintf(fileName, "%s_%04i.hdf5", LOS_params->basename,
            e->los_output_count);
    if (verbose) message("Creating LOS file: %s", fileName);
    h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Error while opening file '%s'.", fileName);
  }

#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Keep track of the total number of parts in all sightlines. */
  size_t total_num_parts_in_los = 0;

  /* ------------------------------- */
  /* Main loop over each random LOS. */
  /* ------------------------------- */

  /* Get list of local non-empty top level cells. */
  const struct cell *cells = e->s->cells_top;
  const int *local_cells_with_particles = e->s->local_cells_with_particles_top;
  const int nr_local_cells_with_particles = s->nr_local_cells_with_particles;

  /* Loop over each random LOS. */
  for (int j = 0; j < LOS_params->num_tot; j++) {

    /* Create empty top level cell list for this LOS */
    int *los_cells_top = (int *)swift_malloc(
        "tl_cells_los", nr_local_cells_with_particles * sizeof(int));
    bzero(los_cells_top, nr_local_cells_with_particles * sizeof(int));

    /* Find all top level cells this LOS will intersect. */
    find_intersecting_top_level_cells(e, &LOS_list[j], los_cells_top, cells,
                                      local_cells_with_particles,
                                      nr_local_cells_with_particles);

    /* Next count all the parts that intersect with this line of sight */
    for (int n = 0; n < nr_local_cells_with_particles; n++) {
      if (los_cells_top[n] == 1) {
        const struct cell *c = &cells[local_cells_with_particles[n]];
        threadpool_map(&s->e->threadpool, los_first_loop_mapper, c->hydro.parts,
                       c->hydro.count, sizeof(struct part),
                       threadpool_auto_chunk_size, &LOS_list[j]);
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    /* Confirm we are capturing all the parts that intersect the LOS by redoing
     * the count looping over all parts in the space (not just those in the
     * subset of top level cells). */

    struct part *parts = s->parts;
    const size_t nr_parts = s->nr_parts;
    const int old_particles_in_los_local = LOS_list[j].particles_in_los_local;
    LOS_list[j].particles_in_los_local = 0;

    /* Count all parts that intersect with this line of sight. */
    threadpool_map(&s->e->threadpool, los_first_loop_mapper, parts, nr_parts,
                   sizeof(struct part), threadpool_auto_chunk_size,
                   &LOS_list[j]);

    /* Make sure we get the same answer as above. */
    if (old_particles_in_los_local != LOS_list[j].particles_in_los_local)
      error("Space vs cells don't match s:%d != c:%d",
            LOS_list[j].particles_in_los_local, old_particles_in_los_local);
#endif

#ifdef WITH_MPI
    /* Make sure all nodes know how many parts are in this LOS */
    int *counts = (int *)malloc(sizeof(int) * e->nr_nodes);
    int *offsets = (int *)malloc(sizeof(int) * e->nr_nodes);

    /* How many parts does each rank have for this LOS? */
    MPI_Allgather(&LOS_list[j].particles_in_los_local, 1, MPI_INT, counts, 1,
                  MPI_INT, MPI_COMM_WORLD);

    int offset_count = 0;
    for (int k = 0; k < e->nr_nodes; k++) {

      /* Total parts in this LOS. */
      LOS_list[j].particles_in_los_total += counts[k];

      /* Counts and offsets for Gatherv. */
      offsets[k] = offset_count;
      offset_count += counts[k];
    }
#else
    LOS_list[j].particles_in_los_total = LOS_list[j].particles_in_los_local;
#endif
    total_num_parts_in_los += LOS_list[j].particles_in_los_total;

    /* Print information about this LOS */
    if (e->nodeID == 0) print_los_info(LOS_list, j);

    /* Don't work with empty LOS */
    if (LOS_list[j].particles_in_los_total == 0) {
      if (e->nodeID == 0) {
        message("*WARNING* LOS %i is empty", j);
        print_los_info(LOS_list, j);
      }
#ifdef WITH_MPI
      free(counts);
      counts = NULL;
      free(offsets);
      offsets = NULL;
#endif
      swift_free("tl_cells_los", los_cells_top);
      continue;
    }

    /* Setup LOS part and xpart structures. */
    struct part *LOS_parts = NULL;
    struct xpart *LOS_xparts = NULL;
    struct gpart *LOS_gparts = NULL;

    /* Rank 0 allocates more space as it will gather all the data for writing */
    if (e->nodeID == 0) {
      if ((LOS_parts = (struct part *)swift_malloc(
               "los_parts_array",
               sizeof(struct part) * LOS_list[j].particles_in_los_total)) ==
          NULL)
        error("Failed to allocate LOS part memory.");
      if ((LOS_xparts = (struct xpart *)swift_malloc(
               "los_xparts_array",
               sizeof(struct xpart) * LOS_list[j].particles_in_los_total)) ==
          NULL)
        error("Failed to allocate LOS xpart memory.");
      if ((LOS_gparts = (struct gpart *)swift_malloc(
               "los_gparts_array",
               sizeof(struct gpart) * LOS_list[j].particles_in_los_total)) ==
          NULL)
        error("Failed to allocate LOS gpart memory.");
    } else {
      if ((LOS_parts = (struct part *)swift_malloc(
               "los_parts_array",
               sizeof(struct part) * LOS_list[j].particles_in_los_local)) ==
          NULL)
        error("Failed to allocate LOS part memory.");
      if ((LOS_xparts = (struct xpart *)swift_malloc(
               "los_xparts_array",
               sizeof(struct xpart) * LOS_list[j].particles_in_los_local)) ==
          NULL)
        error("Failed to allocate LOS xpart memory.");
      if ((LOS_gparts = (struct gpart *)swift_malloc(
               "los_gparts_array",
               sizeof(struct gpart) * LOS_list[j].particles_in_los_local)) ==
          NULL)
        error("Failed to allocate LOS gpart memory.");
    }

    /* Loop over each part again, pulling out those in LOS. */
    int count = 0;

    for (int n = 0; n < e->s->nr_local_cells_with_particles; ++n) {

      if (los_cells_top[n] == 0) continue;

      const struct cell *c = &cells[local_cells_with_particles[n]];
      const struct part *cell_parts = c->hydro.parts;
      const struct xpart *cell_xparts = c->hydro.xparts;
      const size_t num_parts_in_cell = c->hydro.count;

      for (size_t i = 0; i < num_parts_in_cell; i++) {

        /* Don't consider inhibited parts. */
        if (cell_parts[i].time_bin == time_bin_inhibited) continue;
        if (cell_parts[i].time_bin == time_bin_not_created) continue;

        /* Don't consider part if outwith allowed z-range. */
        if (cell_parts[i].x[LOS_list[j].zaxis] <
                LOS_list[j].range_when_shooting_down_axis[0] ||
            cell_parts[i].x[LOS_list[j].zaxis] >
                LOS_list[j].range_when_shooting_down_axis[1])
          continue;

        /* Distance from this part to LOS along x dim. */
        double dx = cell_parts[i].x[LOS_list[j].xaxis] - LOS_list[j].Xpos;

        /* Periodic wrap. */
        if (LOS_list[j].periodic)
          dx = nearest(dx, LOS_list[j].dim[LOS_list[j].xaxis]);

        /* Square */
        const double dx2 = dx * dx;

        /* Smoothing length of this part. */
        const double hsml = cell_parts[i].h * kernel_gamma;
        const double hsml2 = hsml * hsml;

        /* Does this part fall into our LOS? */
        if (dx2 < hsml2) {

          /* Distance from this part to LOS along y dim. */
          double dy = cell_parts[i].x[LOS_list[j].yaxis] - LOS_list[j].Ypos;

          /* Periodic wrap. */
          if (LOS_list[j].periodic)
            dy = nearest(dy, LOS_list[j].dim[LOS_list[j].yaxis]);

          /* Square */
          const double dy2 = dy * dy;

          /* Does this part still fall into our LOS? */
          if (dy2 < hsml2) {
            /* 2D distance to LOS. */

            if (dx2 + dy2 <= hsml2) {

              /* Store part and xpart properties. */
              memcpy(&LOS_parts[count], &cell_parts[i], sizeof(struct part));
              memcpy(&LOS_xparts[count], &cell_xparts[i], sizeof(struct xpart));
              memcpy(&LOS_gparts[count], cell_parts[i].gpart,
                     sizeof(struct gpart));

              count++;
            }
          }
        }
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count != LOS_list[j].particles_in_los_local)
      error("LOS counts don't add up");
#endif

#ifdef WITH_MPI
    /* Collect all parts in this LOS to rank 0. */
    if (e->nodeID == 0) {
      MPI_Gatherv(MPI_IN_PLACE, 0, part_mpi_type, LOS_parts, counts, offsets,
                  part_mpi_type, 0, MPI_COMM_WORLD);
      MPI_Gatherv(MPI_IN_PLACE, 0, xpart_mpi_type, LOS_xparts, counts, offsets,
                  xpart_mpi_type, 0, MPI_COMM_WORLD);
      MPI_Gatherv(MPI_IN_PLACE, 0, gpart_mpi_type, LOS_gparts, counts, offsets,
                  gpart_mpi_type, 0, MPI_COMM_WORLD);

    } else {
      MPI_Gatherv(LOS_parts, LOS_list[j].particles_in_los_local, part_mpi_type,
                  LOS_parts, counts, offsets, part_mpi_type, 0, MPI_COMM_WORLD);
      MPI_Gatherv(LOS_xparts, LOS_list[j].particles_in_los_local,
                  xpart_mpi_type, LOS_xparts, counts, offsets, xpart_mpi_type,
                  0, MPI_COMM_WORLD);
      MPI_Gatherv(LOS_gparts, LOS_list[j].particles_in_los_local,
                  gpart_mpi_type, LOS_gparts, counts, offsets, gpart_mpi_type,
                  0, MPI_COMM_WORLD);
    }
#endif

    /* Re-instate part->gpart pointer on teh receiving side */
    if (e->nodeID == 0) {
#ifdef WITH_MPI
      for (int i = 0; i < LOS_list[j].particles_in_los_total; ++i) {
        LOS_parts[i].gpart = &LOS_gparts[i];
        LOS_gparts[i].id_or_neg_offset = -i;
      }
#endif
    }

    /* Rank 0 writes particles to file. */
    if (e->nodeID == 0) {

      /* Create HDF5 group for this LOS */
      sprintf(groupName, "/LOS_%04i", j);
      h_grp =
          H5Gcreate(h_file, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating LOS HDF5 group\n");

      /* Record this LOS attributes */
      io_write_attribute(h_grp, "NumParts", INT,
                         &LOS_list[j].particles_in_los_total, 1);
      io_write_attribute(h_grp, "Xaxis", INT, &LOS_list[j].xaxis, 1);
      io_write_attribute(h_grp, "Yaxis", INT, &LOS_list[j].yaxis, 1);
      io_write_attribute(h_grp, "Zaxis", INT, &LOS_list[j].zaxis, 1);
      io_write_attribute(h_grp, "Xpos", DOUBLE, &LOS_list[j].Xpos, 1);
      io_write_attribute(h_grp, "Ypos", DOUBLE, &LOS_list[j].Ypos, 1);

#ifdef SWIFT_DEBUG_CHECKS
      for (int i = 0; i < LOS_list[j].particles_in_los_total; ++i) {
        if (LOS_parts[i].gpart != &LOS_gparts[i]) error("Incorrect pointers!");
      }
#endif

      /* Write the data for this LOS */
      write_los_hdf5_datasets(h_grp, j, LOS_list[j].particles_in_los_total,
                              LOS_parts, e, LOS_xparts);

      /* Close HDF5 group */
      H5Gclose(h_grp);
    }

    /* Free up some memory */
#ifdef WITH_MPI
    free(counts);
    free(offsets);
#endif
    swift_free("tl_cells_los", los_cells_top);
    swift_free("los_parts_array", LOS_parts);
    swift_free("los_xparts_array", LOS_xparts);
    swift_free("los_gparts_array", LOS_gparts);

  } /* End of loop over each LOS */

  if (e->nodeID == 0) {
    /* Write header */
    write_hdf5_header(h_file, e, LOS_params, total_num_parts_in_los);

    /* Close HDF5 file */
    H5Fclose(h_file);
  }

  /* Up the LOS counter. */
  e->los_output_count++;

  /* How long did we take? */
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Write a los_props struct to the given FILE as a stream of bytes.
 *
 * @param internal_los the struct
 * @param stream the file stream
 */
void los_struct_dump(const struct los_props *internal_los, FILE *stream) {
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
void los_struct_restore(const struct los_props *internal_los, FILE *stream) {
  restart_read_blocks((void *)internal_los, sizeof(struct los_props), 1, stream,
                      NULL, "los params");
}
