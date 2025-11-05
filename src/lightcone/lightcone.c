/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Some standard headers. */
#include <hdf5.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* HEALPix C API */
#ifdef HAVE_CHEALPIX
#include <chealpix.h>
#endif

/* This object's header. */
#include "lightcone/lightcone.h"

/* Local headers */
#include "common_io.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "extra_io.h"
#include "gravity_io.h"
#include "hydro.h"
#include "lightcone/lightcone_particle_io.h"
#include "lightcone/lightcone_replications.h"
#include "lock.h"
#include "neutrino_io.h"
#include "parser.h"
#include "part_type.h"
#include "particle_buffer.h"
#include "periodic.h"
#include "restart.h"
#include "space.h"
#include "timeline.h"
#include "tools.h"
#include "units.h"

/* Whether to dump the replication list */
// #define DUMP_REPLICATIONS
#ifdef DUMP_REPLICATIONS
static int output_nr = 0;
#endif

/* MPI rank for diagnostic messages */
extern int engine_rank;

#ifdef HAVE_CHEALPIX
/**
 * @brief Read in map types and compression info from a text file
 *
 * The first column in the file is the map type name and the second
 * column is the compression method. Columns are separated by
 * whitespace.
 *
 * @param map_types_file file with the map type names and compression methods
 * @param nr_map_types returns the number of map types read
 * @param map_types returns an array of nr_map_types map type structs
 */
static void read_map_types_file(const char *map_types_file, int *nr_map_types,
                                struct lightcone_map_type **map_types) {

  int map_type_nr = 0;
  if (engine_rank == 0) {

    FILE *fd = fopen(map_types_file, "r");
    if (!fd)
      error("Failed to open lightcone map types file %s", map_types_file);

    /* Count number of non-zero length lines */
    size_t len = 0;
    char *line = NULL;
    int nr_lines = 0;
    while (getline(&line, &len, fd) != -1) nr_lines += 1;
    rewind(fd);

    /* Allocate output arrays */
    *map_types = (struct lightcone_map_type *)calloc(
        nr_lines, sizeof(struct lightcone_map_type));

    /* Read lines */
    for (int i = 0; i < nr_lines; i += 1) {

      /* Get name and compression type from this line */
      char compression[PARSER_MAX_LINE_SIZE];
      if (fscanf(fd, "%s %s", (*map_types)[map_type_nr].name, compression) != 2)
        error("Failed to read line from map types file");

      /* Look up compression scheme */
      (*map_types)[map_type_nr].compression =
          compression_scheme_from_name(compression);

      /* Only keep maps which have not been disabled */
      if ((*map_types)[map_type_nr].compression != compression_do_not_write)
        map_type_nr += 1;
    }
    fclose(fd);
    free(line);
  }

#ifdef WITH_MPI
  MPI_Bcast(&map_type_nr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (engine_rank != 0)
    *map_types = (struct lightcone_map_type *)calloc(
        map_type_nr, sizeof(struct lightcone_map_type));
  MPI_Bcast(*map_types, sizeof(struct lightcone_map_type) * map_type_nr,
            MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Return number of enabled map types */
  *nr_map_types = map_type_nr;
}
#endif

/**
 * @brief Identify which healpix map types we're making
 *
 * @param props the #lightcone_props structure
 *
 * For each requested map type find the update functions by matching names.
 * Map types are defined in lightcone_map_types.h and there may be extra
 * types defined by various physics modules. The array map_type_array
 * below determines where we look for extra map types.
 *
 * This function assumes that props->map_type is already allocated and
 * props->map_type[:].name has been set to the list of map names from the
 * .yml file. It sets the update_map, ptype_contributes and units fields
 * in the props->map_type array.
 *
 */
static void lightcone_identify_map_types(struct lightcone_props *props) {

  /* Loop over requested map types */
  for (int map_nr = 0; map_nr < props->nr_maps; map_nr += 1) {

    /* Use null function pointer to indicate not found yet */
    props->map_type[map_nr].update_map = NULL;

    /* Places to search for lightcone map types:
       extra map types are provided by various physics modules. */
    const int num_places = 3;
    const struct lightcone_map_type *map_type_array[] = {
        lightcone_map_types, extra_lightcone_map_types,
        neutrino_lightcone_map_types};

    /* Loop over places to search for map types */
    for (int i = 0; i < num_places; i += 1) {

      int type_nr = 0;
      const struct lightcone_map_type *map_types_to_search = map_type_array[i];
      while (map_types_to_search[type_nr].update_map) {
        if (strcmp(map_types_to_search[type_nr].name,
                   props->map_type[map_nr].name) == 0) {
          props->map_type[map_nr] = map_types_to_search[type_nr];
          if (engine_rank == 0)
            message("lightcone %d: lightcone map %d is of type %s",
                    props->index, map_nr, map_types_to_search[type_nr].name);
        }
        type_nr += 1;
      }

    } /* Next place to search */

    if (!props->map_type[map_nr].update_map)
      error("Unable to locate lightcone map type %s",
            props->map_type[map_nr].name);
  }
}

/**
 * @brief Allocate particle I/O buffers for a lightcone
 *
 * @param props the #lightcone_props structure
 */
static void lightcone_allocate_buffers(struct lightcone_props *props) {

  /* Initialize particle output buffers */
  const size_t elements_per_block = (size_t)props->buffer_chunk_size;

  if (props->use_type[swift_type_gas]) {
    particle_buffer_init(&props->buffer[swift_type_gas],
                         sizeof(struct lightcone_gas_data), elements_per_block,
                         "lightcone_gas");
  }

  if (props->use_type[swift_type_dark_matter]) {
    particle_buffer_init(&props->buffer[swift_type_dark_matter],
                         sizeof(struct lightcone_dark_matter_data),
                         elements_per_block, "lightcone_dm");
  }

  if (props->use_type[swift_type_dark_matter_background]) {
    particle_buffer_init(&props->buffer[swift_type_dark_matter_background],
                         sizeof(struct lightcone_dark_matter_data),
                         elements_per_block, "lightcone_dm_bg");
  }

  if (props->use_type[swift_type_stars]) {
    particle_buffer_init(&props->buffer[swift_type_stars],
                         sizeof(struct lightcone_stars_data),
                         elements_per_block, "lightcone_stars");
  }

  if (props->use_type[swift_type_black_hole]) {
    particle_buffer_init(&props->buffer[swift_type_black_hole],
                         sizeof(struct lightcone_black_hole_data),
                         elements_per_block, "lightcone_bh");
  }

  if (props->use_type[swift_type_neutrino]) {
    particle_buffer_init(&props->buffer[swift_type_neutrino],
                         sizeof(struct lightcone_neutrino_data),
                         elements_per_block, "lightcone_neutrino");
  }
}

/**
 * @brief Dump lightcone_props struct to the output stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to write to.
 */
void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  /* Don't dump the replication list - will regenerate it as needed */
  struct lightcone_props tmp = *props;
  tmp.replication_list.nrep = 0;
  tmp.replication_list.replication = NULL;
  tmp.have_replication_list = 0;

  /* Don't write out particle buffers - must flush before dumping restart. */
  memset(tmp.buffer, 0, sizeof(struct particle_buffer) * swift_type_count);

  /* Don't write array pointers */
  tmp.shell = NULL;
  tmp.map_type = NULL;
  for (int ptype = 0; ptype < swift_type_count; ptype += 1)
    tmp.part_type[ptype].map_index = NULL;

  /* Dump the lightcone struct */
  restart_write_blocks((void *)&tmp, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");

  /* Dump the array of map types */
  restart_write_blocks((void *)props->map_type,
                       sizeof(struct lightcone_map_type), props->nr_maps,
                       stream, "lightcone_props", "lightcone_props");

  /* Dump the array of shells */
  lightcone_shell_array_dump(props->shell, props->nr_shells, stream);

  /* For each particle type we have an array of lightcone map indexes to update.
   * Dump these. */
  for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
    const struct lightcone_particle_type *this_type =
        &(props->part_type[ptype]);
    restart_write_blocks((void *)this_type->map_index, sizeof(int),
                         this_type->nr_maps, stream, "lightcone_props",
                         "lightcone_props");
  }
}

/**
 * @brief Initialise the particle output fields for each particle type.
 *
 * @param props the #lightcone_props structure
 */
void lightcone_define_output_fields(struct lightcone_props *props) {

  for (int ptype = 0; ptype < swift_type_count; ptype += 1)
    lightcone_io_field_list_init(&props->particle_fields[ptype]);

  /* Add the default set of fields for all models, from lightcone_particle_io.c
   */
  lightcone_io_append_gas_output_fields(
      &props->particle_fields[swift_type_gas]);
  lightcone_io_append_dark_matter_output_fields(
      &props->particle_fields[swift_type_dark_matter]);
  lightcone_io_append_dark_matter_background_output_fields(
      &props->particle_fields[swift_type_dark_matter_background]);
  lightcone_io_append_stars_output_fields(
      &props->particle_fields[swift_type_stars]);
  lightcone_io_append_black_hole_output_fields(
      &props->particle_fields[swift_type_black_hole]);
  lightcone_io_append_neutrino_output_fields(
      &props->particle_fields[swift_type_neutrino]);
}

/**
 * @brief Restore lightcone_props struct from the input stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to read from.
 */
void lightcone_struct_restore(struct lightcone_props *props, FILE *stream) {

  /* Restore lightcone struct */
  restart_read_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                      NULL, "lightcone_props");

  /* Read in the map types */
  props->map_type = (struct lightcone_map_type *)malloc(
      sizeof(struct lightcone_map_type) * props->nr_maps);
  restart_read_blocks((void *)props->map_type,
                      sizeof(struct lightcone_map_type), props->nr_maps, stream,
                      NULL, "lightcone_props");

  /* Read in the shells */
  props->shell = lightcone_shell_array_restore(
      stream, props->nr_shells, props->part_type, props->buffer_chunk_size);

  /* For each particle type we have an array of lightcone map indexes to update.
   * Restore these. */
  for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
    struct lightcone_particle_type *this_type = &(props->part_type[ptype]);
    this_type->map_index = (int *)malloc(sizeof(int) * this_type->nr_maps);
    restart_read_blocks((void *)this_type->map_index, sizeof(int),
                        this_type->nr_maps, stream, NULL, "lightcone_props");
  }

  /* Restore pointers to functions for updating healpix maps */
  lightcone_identify_map_types(props);

  /* Update function pointers for each map */
  for (int shell_nr = 0; shell_nr < props->nr_shells; shell_nr += 1) {
    for (int map_nr = 0; map_nr < props->nr_maps; map_nr += 1) {
      props->shell[shell_nr].map[map_nr].type = props->map_type[map_nr];
    }
  }

  /* Re-allocate particle data buffers */
  lightcone_allocate_buffers(props);

  /* Define output quantities */
  lightcone_define_output_fields(props);

  /* Tabulate the projected kernel */
  projected_kernel_init(&props->kernel_table);
}

#ifdef HAVE_CHEALPIX
/**
 * @brief Locate a lightcone parameter in the .yml file
 *
 * First check the section specific to this lightcone then
 * fall back to LightconeCommon if not found.
 *
 * @param params the swift parameters struct
 * @param index index of the lightcone
 * @param name name of the parameter to find
 * @param outbuf returns the parameter value
 *
 */
static char *find_parameter(struct swift_params *params, const int index,
                            const char *name, char *outbuf) {

  char full_name[PARSER_MAX_LINE_SIZE];

  /* Check section specific to this lightcone */
  check_snprintf(full_name, PARSER_MAX_LINE_SIZE, "Lightcone%d:%s", index,
                 name);
  if (parser_does_param_exist(params, full_name)) {
    strcpy(outbuf, full_name);
    return outbuf;
  }

  /* Will look in LightconeCommon section if parameter was not found */
  check_snprintf(full_name, PARSER_MAX_LINE_SIZE, "LightconeCommon:%s", name);
  strcpy(outbuf, full_name);
  return outbuf;
}
#endif

/**
 * @brief Initialise the properties of the lightcone code.
 *
 * @param props the #lightcone_props structure to fill.
 * @param index index of the lightcone to initialize
 * @param s the #space structure.
 * @param cosmo the #cosmology structure.
 * @param params the parameter file parser.
 * @param internal_units swift internal unit system
 * @param physical_constants swift physical constant values
 * @param verbose the verbosity flag
 */
void lightcone_init(struct lightcone_props *props, const int index,
                    const struct space *s, const struct cosmology *cosmo,
                    struct swift_params *params,
                    const struct unit_system *internal_units,
                    const struct phys_const *physical_constants,
                    const int verbose) {

#ifdef HAVE_CHEALPIX

  /* Macro to generate parameter names given section name */
  char buf[PARSER_MAX_LINE_SIZE];
#define YML_NAME(x) find_parameter(params, index, x, buf)

  /* Store index of this lightcone in the .yml file */
  props->index = index;

  /* Verbose lightcone output - use passed in value of --verbose flag */
  props->verbose = verbose;

  /* Define output quantities */
  lightcone_define_output_fields(props);

  /* For each particle type, get redshift range for lightcone particle output */
  for (int i = 0; i < swift_type_count; i += 1) {
    const int len = PARSER_MAX_LINE_SIZE;
    char param_name[len];
    double zrange[2] = {0.0, -1.0}; /* default max < min means do not output */
    check_snprintf(param_name, len, "z_range_for_%s", part_type_names[i]);
    parser_get_opt_param_double_array(params, YML_NAME(param_name), 2, zrange);
    props->z_min_for_type[i] = zrange[0];
    props->z_max_for_type[i] = zrange[1];
    /* Will only output types with z_max > z_min */
    props->use_type[i] = props->z_max_for_type[i] > props->z_min_for_type[i];
    if (engine_rank == 0 && verbose) {
      if (props->use_type[i]) {
        message("lightcone %d: %s particles will be output from z=%f to z=%f",
                props->index, part_type_names[i], zrange[0], zrange[1]);
      } else {
        message("lightcone %d: %s particle output is disabled", props->index,
                part_type_names[i]);
      }
    }
  }

  /* For each type, find range in comoving distance squared in which we output
   * particles */
  for (int i = 0; i < swift_type_count; i += 1) {
    if (props->use_type[i]) {
      const double a_min = 1.0 / (1.0 + props->z_max_for_type[i]);
      props->r2_max_for_type[i] =
          pow(cosmology_get_comoving_distance(cosmo, a_min), 2.0);
      const double a_max = 1.0 / (1.0 + props->z_min_for_type[i]);
      props->r2_min_for_type[i] =
          pow(cosmology_get_comoving_distance(cosmo, a_max), 2.0);
    } else {
      props->r2_min_for_type[i] = 0.0;
      props->r2_max_for_type[i] = 0.0;
    }
  }

  /*
    Allow selective output of gas particles at high redshift.
    Will output gas particles if redshift < min_z_for_gas_filtering OR
    (temperature > min_temp_for_filtered_gas AND nh >
    min_nh_for_filtered_gas*(1+z)^4)
  */
  props->gas_filtering_enabled =
      parser_get_opt_param_int(params, YML_NAME("gas_filtering_enabled"), 0);
  if (props->gas_filtering_enabled) {
    props->min_z_for_gas_filtering =
        parser_get_param_double(params, YML_NAME("min_z_for_gas_filtering"));
    props->min_temp_for_filtered_gas =
        parser_get_param_double(params, YML_NAME("min_temp_for_filtered_gas"));
    props->min_nh_for_filtered_gas =
        parser_get_param_double(params, YML_NAME("min_nh_for_filtered_gas"));
    props->max_a_for_gas_filtering =
        1.0 / (1.0 + props->min_z_for_gas_filtering);
    /* Convert temperature and density thresholds to internal units, assuming
     * they're input in CGS */
    props->min_temp_for_filtered_gas /=
        units_cgs_conversion_factor(internal_units, UNIT_CONV_TEMPERATURE);
    props->min_nh_for_filtered_gas /=
        units_cgs_conversion_factor(internal_units, UNIT_CONV_NUMBER_DENSITY);
  }

  /* Exclude particles from xray and sz maps if they have been recently AGN
   * heated */
  props->xray_maps_recent_AGN_injection_exclusion_time =
      parser_get_opt_param_double(
          params, YML_NAME("xray_maps_recent_AGN_injection_exclusion_time_myr"),
          -1.0);
  /* Assume supplied value is in megayears and physical constants are in
   * internal units */
  props->xray_maps_recent_AGN_injection_exclusion_time *=
      1.0e6 * physical_constants->const_year;

  /*
    Temperature limits for recently AGN heated gas to be excluded from xray and
    sz maps.

    Gas is excluded if it has log(temperature) which is greater than
    log(AGN_Delta_T)+xray_maps_recent_AGN_logdT_min and less than
    log(AGN_Delta_T)+xray_maps_recent_AGN_logdT_max.

    Only takes effect if xray_maps_recent_AGN_injection_exclusion_time_myr is
    set.
  */
  if (props->xray_maps_recent_AGN_injection_exclusion_time > 0.0) {
    double delta_logt_min = parser_get_param_double(
        params, YML_NAME("xray_maps_recent_AGN_injection_delta_logT_min"));
    if (delta_logt_min > 0.0)
      error("xray_maps_recent_AGN_injection_delta_logT_min should be negative");
    props->xray_maps_recent_AGN_min_temp_factor = pow(10.0, delta_logt_min);
    double delta_logt_max = parser_get_param_double(
        params, YML_NAME("xray_maps_recent_AGN_injection_delta_logT_max"));
    if (delta_logt_max < 0.0)
      error("xray_maps_recent_AGN_injection_delta_logT_max should be positive");
    props->xray_maps_recent_AGN_max_temp_factor = pow(10.0, delta_logt_max);
    if (delta_logt_max < delta_logt_min)
      error(
          "xray_maps_recent_AGN_injection_delta_logT_max should be greater "
          "than _min!");
  }

  /* Directory in which to write this lightcone */
  parser_get_opt_param_string(params, YML_NAME("subdir"), props->subdir, ".");

  /* Base name for output files */
  parser_get_param_string(params, YML_NAME("basename"), props->basename);

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, YML_NAME("observer_position"), 3,
                                props->observer_position);

  /* Write particles to disk if this many or more are in the buffer */
  props->max_particles_buffered = parser_get_opt_param_int(
      params, YML_NAME("max_particles_buffered"), 100000);

  /* Chunk size for particles buffered in memory */
  props->buffer_chunk_size =
      parser_get_opt_param_int(params, YML_NAME("buffer_chunk_size"), 20000);

  /* Chunk size for particles in the HDF5 output files */
  props->hdf5_chunk_size =
      parser_get_opt_param_int(params, YML_NAME("hdf5_chunk_size"), 16384);

  /* Maximum amount of data (in megabytes) to send from any one rank when
   * updating healpix maps */
  props->max_map_update_send_size_mb = parser_get_opt_param_double(
      params, YML_NAME("max_map_update_send_size_mb"), 512.0);

  /* Compression options */
  props->particles_lossy_compression = parser_get_opt_param_int(
      params, YML_NAME("particles_lossy_compression"), 0);
  props->particles_gzip_level =
      parser_get_opt_param_int(params, YML_NAME("particles_gzip_level"), 0);
  props->maps_gzip_level =
      parser_get_opt_param_int(params, YML_NAME("maps_gzip_level"), 0);

  /* Get the size of the simulation box */
  props->boxsize = s->dim[0];
  if (s->dim[1] != s->dim[0] || s->dim[2] != s->dim[0])
    error("Lightcones require a cubic simulation box.");

  /* Get top level cell size */
  props->cell_width = s->width[0];
  if (s->width[1] != s->width[0] || s->width[2] != s->width[0])
    error("Lightcones require cubic top level cells.");

  /* Initially have no replication list */
  props->have_replication_list = 0;
  props->ti_old = 0;
  props->ti_current = 0;

  /* Initialize various counters */
  for (int i = 0; i < swift_type_count; i += 1) {
    props->num_particles_written_this_rank[i] = 0;
    props->num_particles_written_to_file[i] = 0;
  }
  props->current_file = -1;
  props->file_needs_finalizing = 0;

  /* Always start a new file initially */
  props->start_new_file = 1;

  /*
     Healpix map parameters for this lightcone
  */

  /* Healpix nside parameter */
  props->nside = parser_get_param_int(params, YML_NAME("nside"));

  /* Update lightcone pixel data if more than this number of updates are
   * buffered */
  props->max_updates_buffered = parser_get_opt_param_int(
      params, YML_NAME("max_updates_buffered"), 1000000);

  /*! Whether to write distributed maps in MPI mode */
  props->distributed_maps =
      parser_get_opt_param_int(params, YML_NAME("distributed_maps"), 1);

  /* Name of the file with radii of spherical shells */
  parser_get_param_string(params, YML_NAME("radius_file"), props->radius_file);

  /* Get names of the healpix maps to make for this lightcone */
  char map_types_file[FILENAME_BUFFER_SIZE];
  parser_get_param_string(params, YML_NAME("map_names_file"), map_types_file);
  read_map_types_file(map_types_file, &props->nr_maps, &props->map_type);

  /* For each requested map type find the update function by matching names */
  lightcone_identify_map_types(props);

  /* For each particle type, determine which healpix maps will be updated */
  const int nr_maps = props->nr_maps;
  for (int ptype = 0; ptype < swift_type_count; ptype += 1) {

    struct lightcone_particle_type *this_type = &(props->part_type[ptype]);

    /* Count maps updated by this particle type */
    this_type->nr_maps = 0;
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      if (props->map_type[map_nr].ptype_contributes(ptype))
        this_type->nr_maps += 1;
    }

    /* Store indexes of maps to update for this particle type */
    this_type->map_index = (int *)malloc(sizeof(int) * this_type->nr_maps);
    this_type->nr_maps = 0;
    this_type->nr_smoothed_maps = 0;
    this_type->nr_unsmoothed_maps = 0;

    /* First the smoothed maps */
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      const struct lightcone_map_type *map_type = &(props->map_type[map_nr]);
      if (map_type->ptype_contributes(ptype) &&
          map_type->smoothing == map_smoothed) {
        this_type->map_index[this_type->nr_maps] = map_nr;
        this_type->nr_maps += 1;
        this_type->nr_smoothed_maps += 1;
      }
    }

    /* Then the un-smoothed maps */
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      const struct lightcone_map_type *map_type = &(props->map_type[map_nr]);
      if (map_type->ptype_contributes(ptype) &&
          map_type->smoothing == map_unsmoothed) {
        this_type->map_index[this_type->nr_maps] = map_nr;
        this_type->nr_maps += 1;
        this_type->nr_unsmoothed_maps += 1;
      }
    }

    /* Determine how much data we need to store per particle of this type.
       We need theta and phi angular coordinates, angular size of the particle,
       and the values to be added to the healpix maps */
    this_type->buffer_element_size =
        (3 + this_type->nr_maps) * sizeof(union lightcone_map_buffer_entry);
  }

  /* Check the number of healpix pixels doesn't overflow pixel_index_t */
  const unsigned long long nside_ull = props->nside;
  const unsigned long long npix_ull = 12ull * nside_ull * nside_ull;
  if (npix_ull > MAX_PIXEL_INDEX)
    error(
        "Number of HEALPix pixels is to large for pixel_index_t (see "
        "lightcone/pixel_index.h)");

  /* Set up the array of lightcone shells for this lightcone */
  const pixel_index_t total_nr_pix = nside2npix64(props->nside);
  props->shell = lightcone_shell_array_init(
      cosmo, props->radius_file, props->nr_maps, props->map_type, props->nside,
      total_nr_pix, props->part_type, props->buffer_chunk_size,
      &props->nr_shells);

  /* Compute area of a healpix pixel */
  props->pixel_area_steradians = 4 * M_PI / total_nr_pix;

  /* Report shell radii */
  const int nr_shells = props->nr_shells;
  if (engine_rank == 0) {
    for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
      message("lightcone %d: shell %d has inner radius %e and outer radius %e",
              index, shell_nr, props->shell[shell_nr].rmin,
              props->shell[shell_nr].rmax);
    }
  }
  if (engine_rank == 0)
    message("lightcone %d: there are %d lightcone shells and %d maps per shell",
            index, nr_shells, nr_maps);

  /* For each particle type, find the full redshift range to search for
   * lightcone crossings */
  int have_particle_output = 0;
  for (int i = 0; i < swift_type_count; i += 1) {

    /* Initially set range to search to range used for particle output, if any
     */
    if (props->use_type[i]) {
      props->a_min_search_for_type[i] = 1.0 / (1.0 + props->z_max_for_type[i]);
      props->a_max_search_for_type[i] = 1.0 / (1.0 + props->z_min_for_type[i]);
      have_particle_output = 1;
    } else {
      props->a_min_search_for_type[i] = DBL_MAX;
      props->a_max_search_for_type[i] = 0.0;
    }

    /* Then expand the range to include any healpix maps this type contributes
     * to */
    for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
      const double shell_a_min = props->shell[shell_nr].amin;
      const double shell_a_max = props->shell[shell_nr].amax;
      for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
        const struct lightcone_map_type *map_type = &(props->map_type[map_nr]);
        if (map_type->ptype_contributes(i)) {
          if (shell_a_min < props->a_min_search_for_type[i])
            props->a_min_search_for_type[i] = shell_a_min;
          if (shell_a_max > props->a_max_search_for_type[i])
            props->a_max_search_for_type[i] = shell_a_max;
        }
      }
    }
    /* Next particle type */
  }

  /* Determine the full redshift range to search for all particle types */
  double a_min = DBL_MAX;
  double a_max = 0.0;
  for (int i = 0; i < swift_type_count; i += 1) {
    if (props->a_max_search_for_type[i] > props->a_min_search_for_type[i]) {
      if (props->a_min_search_for_type[i] < a_min)
        a_min = props->a_min_search_for_type[i];
      if (props->a_max_search_for_type[i] > a_max)
        a_max = props->a_max_search_for_type[i];
    }
  }

  /* Check we have a valid range in expansion factor for the lightcone */
  if (a_min > a_max)
    error(
        "Code was run with --lightcone but no particle outputs or healpix maps "
        "are enabled");
  props->a_min = a_min;
  props->a_max = a_max;
  if (engine_rank == 0) {
    for (int i = 0; i < swift_type_count; i += 1) {
      if (props->a_max_search_for_type[i] > props->a_min_search_for_type[i]) {
        message("lightcone %d: range in expansion factor for %s: %e to %e",
                index, part_type_names[i], props->a_min_search_for_type[i],
                props->a_max_search_for_type[i]);
      } else {
        message("lightcone %d: no lightcone output for %s", index,
                part_type_names[i]);
      }
    }
    message("lightcone %d: range in expansion factor overall: %e to %e", index,
            a_min, a_max);
  }

  /* Store the corresponding comoving distance squared */
  props->r2_max = pow(cosmology_get_comoving_distance(cosmo, a_min), 2.0);
  props->r2_min = pow(cosmology_get_comoving_distance(cosmo, a_max), 2.0);

  /* Allocate lightcone output buffers */
  lightcone_allocate_buffers(props);

  /* Tabulate the projected kernel */
  projected_kernel_init(&props->kernel_table);

  /* Ensure that the output directories exist */
  if (engine_rank == 0) {
    const int len = FILENAME_BUFFER_SIZE;
    char dirname[len];
    safe_checkdir(props->subdir, 1);
    /* Directory for particle outputs */
    if (have_particle_output) {
      check_snprintf(dirname, len, "%s/%s_particles", props->subdir,
                     props->basename);
      safe_checkdir(dirname, 1);
    }
    /* Directory for shell outputs */
    if ((props->nr_shells > 0) && (props->nr_maps > 0)) {
      check_snprintf(dirname, len, "%s/%s_shells", props->subdir,
                     props->basename);
      safe_checkdir(dirname, 1);
    }
  }
#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#else
  error("Need HEALPix C API to make lightcones");
#endif
}

/**
 * @brief Return the name of a lightcone particle output file
 *
 * @param buf returns the filename
 * @param len length of the buffer buf
 * @param subdir subdirectory in which to write output
 * @param basename base name of this lightcone
 * @param current_file lightcone particle file index, which is
 *        incremented after each restart dump
 * @param comm_rank rank of this MPI communicator
 */
static void particle_file_name(char *buf, int len, char *subdir, char *basename,
                               int current_file, int comm_rank) {

  check_snprintf(buf, len, "%s/%s_particles/%s_%04d.%d.hdf5", subdir, basename,
                 basename, current_file, comm_rank);
}

/**
 * @brief Flush any buffers which exceed the specified size.
 *
 * Also used to flush buffers before dumping restart files, in
 * which case we should have flush_all=1 and end_file=1 so that
 * buffers are flushed regardless of size and we will start a
 * new set of lightcone files after the restart dump.
 *
 * @param props the #lightcone_props structure.
 * @param a the current expansion factor
 * @param internal_units swift internal unit system
 * @param snapshot_units swift snapshot unit system
 * @param flush_all flag to force flush of all buffers
 * @param end_file if true, subsequent calls write to a new file
 *
 */
void lightcone_flush_particle_buffers(struct lightcone_props *props, double a,
                                      const struct unit_system *internal_units,
                                      const struct unit_system *snapshot_units,
                                      int flush_all, int end_file) {

  ticks tic = getticks();

  /* Should never be called with end_file=1 and flush_all=0 */
  if (end_file && (!flush_all))
    error("Finalizing file without flushing buffers!");

  /* Will flush any buffers with more particles than this */
  size_t max_to_buffer = (size_t)props->max_particles_buffered;
  if (flush_all) max_to_buffer = 0;

  /* Count how many types have data to write out */
  int types_to_flush = 0;
  for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
    if (props->use_type[ptype]) {
      const size_t num_to_write =
          particle_buffer_num_elements(&props->buffer[ptype]);
      if (num_to_write >= max_to_buffer && num_to_write > 0)
        types_to_flush += 1;
    }
  }

  /* Check if there's anything to do */
  if ((types_to_flush > 0) || (end_file && props->file_needs_finalizing)) {

    /* We have data to flush, so open or create the output file */
    hid_t file_id, h_props;
    char fname[FILENAME_BUFFER_SIZE];
    if (props->start_new_file) {

      /* Get the name of the next file */
      props->current_file += 1;
      particle_file_name(fname, FILENAME_BUFFER_SIZE, props->subdir,
                         props->basename, props->current_file, engine_rank);

      h_props = H5Pcreate(H5P_FILE_ACCESS);
      herr_t err =
          H5Pset_libver_bounds(h_props, HDF5_LOWEST_FILE_FORMAT_VERSION,
                               HDF5_HIGHEST_FILE_FORMAT_VERSION);
      if (err < 0) error("Error setting the hdf5 API version");

      /* Create the file */
      file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, h_props);
      if (file_id < 0) error("Unable to create new lightcone file: %s", fname);

      /* This new file has not been finalized yet */
      props->file_needs_finalizing = 1;

      /* We have now written no particles to the current file */
      for (int ptype = 0; ptype < swift_type_count; ptype += 1)
        props->num_particles_written_to_file[ptype] = 0;

      /* Write the system of Units used in the snapshot */
      io_write_unit_system(file_id, snapshot_units, "Units");

      /* Write the system of Units used internally */
      io_write_unit_system(file_id, internal_units, "InternalCodeUnits");

      /* Write the observer position and redshift limits */
      hid_t group_id = H5Gcreate(file_id, "Lightcone", H5P_DEFAULT, H5P_DEFAULT,
                                 H5P_DEFAULT);
      io_write_attribute(group_id, "observer_position", DOUBLE,
                         props->observer_position, 3);

      for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
        char name[PARSER_MAX_LINE_SIZE];
        check_snprintf(name, PARSER_MAX_LINE_SIZE, "minimum_redshift_%s",
                       part_type_names[ptype]);
        io_write_attribute_d(group_id, name, props->z_min_for_type[ptype]);
        check_snprintf(name, PARSER_MAX_LINE_SIZE, "maximum_redshift_%s",
                       part_type_names[ptype]);
        io_write_attribute_d(group_id, name, props->z_max_for_type[ptype]);
      }

      /* Record number of MPI ranks so we know how many files there are */
      int comm_rank = 0;
      int comm_size = 1;
#ifdef WITH_MPI
      MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
      MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif
      io_write_attribute_i(group_id, "mpi_rank", comm_rank);
      io_write_attribute_i(group_id, "nr_mpi_ranks", comm_size);
      io_write_attribute_i(group_id, "file_index", props->current_file);

      H5Gclose(group_id);

      /* We no longer need to create a new file */
      props->start_new_file = 0;

    } else {

      h_props = H5Pcreate(H5P_FILE_ACCESS);
      herr_t err =
          H5Pset_libver_bounds(h_props, HDF5_LOWEST_FILE_FORMAT_VERSION,
                               HDF5_HIGHEST_FILE_FORMAT_VERSION);
      if (err < 0) error("Error setting the hdf5 API version");

      /* Re-open an existing file */
      particle_file_name(fname, FILENAME_BUFFER_SIZE, props->subdir,
                         props->basename, props->current_file, engine_rank);
      file_id = H5Fopen(fname, H5F_ACC_RDWR, h_props);
      if (file_id < 0)
        error("Unable to open current lightcone file: %s", fname);
    }

    /* Loop over particle types */
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      if (props->use_type[ptype]) {
        const size_t num_to_write =
            particle_buffer_num_elements(&props->buffer[ptype]);
        if (num_to_write >= max_to_buffer && num_to_write > 0) {
          lightcone_write_particles(props, internal_units, snapshot_units,
                                    ptype, file_id);
          particle_buffer_empty(&props->buffer[ptype]);
          props->num_particles_written_to_file[ptype] += num_to_write;
          props->num_particles_written_this_rank[ptype] += num_to_write;
        }
      }
    }

    /* Check if this is the last write to this file */
    if (end_file) {
      hid_t group_id = H5Gopen(file_id, "Lightcone", H5P_DEFAULT);
      /* Flag the file as complete */
      io_write_attribute_i(group_id, "file_complete", 1);
      /* Write the expected number of particles in all files written by this
         rank up to and including this one. */
      for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
        char name[PARSER_MAX_LINE_SIZE];
        check_snprintf(name, PARSER_MAX_LINE_SIZE, "cumulative_count_%s",
                       part_type_names[ptype]);
        io_write_attribute_ll(group_id, name,
                              props->num_particles_written_this_rank[ptype]);
      }
      /* Write the expansion factor at which we closed this file */
      io_write_attribute_d(group_id, "expansion_factor", a);
      H5Gclose(group_id);
      props->file_needs_finalizing = 0;
    }

    /* We're done updating the output file */
    H5Fclose(file_id);
    H5Pclose(h_props);
  }

  /* If we need to start a new file next time, record this */
  if (end_file) props->start_new_file = 1;

  if (props->verbose && engine_rank == 0 && types_to_flush > 0)
    message("lightcone %d: Flushing particle buffers took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Flush lightcone map update buffers for all shells
 *
 * @param props the #lightcone_props structure.
 * @param tp the swift #threadpool struct to use
 *
 */
void lightcone_flush_map_updates(struct lightcone_props *props,
                                 struct threadpool *tp) {

  ticks tic = getticks();

  /* Apply updates to all current shells */
  for (int shell_nr = 0; shell_nr < props->nr_shells; shell_nr += 1) {
    if (props->shell[shell_nr].state == shell_current) {
      lightcone_shell_flush_map_updates(&props->shell[shell_nr], tp,
                                        props->part_type,
                                        props->max_map_update_send_size_mb,
                                        &props->kernel_table, props->verbose);
    }
  }

  /* Report runtime */
  if (props->verbose && engine_rank == 0)
    message("lightcone %d: Applying lightcone map updates took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Write and deallocate any completed lightcone shells
 *
 * @param props the #lightcone_props structure.
 * @param tp the swift #threadpool struct to use
 * @param c the #cosmology structure
 * @param internal_units swift internal unit system
 * @param snapshot_units swift snapshot unit system
 * @param dump_all flag to indicate that all shells should be dumped
 * @param need_flush whether there might be buffered updates to apply
 *
 */
void lightcone_dump_completed_shells(struct lightcone_props *props,
                                     struct threadpool *tp,
                                     const struct cosmology *c,
                                     const struct unit_system *internal_units,
                                     const struct unit_system *snapshot_units,
                                     const int dump_all, const int need_flush) {
#ifdef HAVE_HDF5

  ticks tic = getticks();

  int comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif

  /* Get number of shells and maps per shell */
  const int nr_shells = props->nr_shells;
  const int nr_maps = props->nr_maps;

  /* Get conversion factor for shell radii */
  const double length_conversion_factor =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);

  /* Compute expansion factor corresponding to time props->ti_old,
     which is the earliest time any particle might have been drifted
     from on this step. Here we assume that no particle remains to
     be drifted from any time earlier than this so that any shell
     whose redshift range is entirely before ti_old can be now be
     written out and deallocated. */
  const double a_complete = c->a_begin * exp(props->ti_old * c->time_base);

  int num_shells_written = 0;
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {

    /* Will write out this shell if it has been updated but not written
       out yet and either we advanced past its redshift range or we're
       dumping all remaining shells at the end of the simulation */
    if (props->shell[shell_nr].state == shell_current) {
      if (props->shell[shell_nr].amax < a_complete || dump_all) {

        if (props->verbose && engine_rank == 0)
          message("lightcone %d: writing out completed shell %d at a=%f",
                  props->index, shell_nr, c->a);

        num_shells_written += 1;

        /* Apply any buffered updates for this shell, if we didn't already */
        if (need_flush) {
          lightcone_shell_flush_map_updates(
              &props->shell[shell_nr], tp, props->part_type,
              props->max_map_update_send_size_mb, &props->kernel_table,
              props->verbose);
        }

        /* Set the baseline value for the maps */
        for (int map_nr = 0; map_nr < nr_maps; map_nr += 1)
          lightcone_map_set_baseline(c, props,
                                     &(props->shell[shell_nr].map[map_nr]));

        /* Ensure output directory exists */
        char fname[FILENAME_BUFFER_SIZE];
        check_snprintf(fname, FILENAME_BUFFER_SIZE, "%s/%s_shells/shell_%d",
                       props->subdir, props->basename, shell_nr);
        if (engine_rank == 0) safe_checkdir(fname, 1);
#ifdef WITH_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        /* Get the name of the file to write:
           In collective mode all ranks get the same file name.
           In distributed mode we include engine_rank in the file name. */
        int file_num = props->distributed_maps ? engine_rank : 0;
        check_snprintf(fname, FILENAME_BUFFER_SIZE,
                       "%s/%s_shells/shell_%d/%s.shell_%d.%d.hdf5",
                       props->subdir, props->basename, shell_nr,
                       props->basename, shell_nr, file_num);

        /* Create the output file for this shell */
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);

        /* Set the minimal API version to avoid issues with advanced features */
        herr_t err =
            H5Pset_libver_bounds(fapl_id, HDF5_LOWEST_FILE_FORMAT_VERSION,
                                 HDF5_HIGHEST_FILE_FORMAT_VERSION);
        if (err < 0) error("Error setting the hdf5 API version");

        /* Set MPI collective mode, if necessary */
        int collective = 0;
#ifdef WITH_MPI
#ifdef HAVE_PARALLEL_HDF5
        if (!props->distributed_maps) {
          if (H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
            error("Unable to set HDF5 MPI-IO file access mode");
          collective = 1;
        }
#else
        if (!props->distributed_maps)
          error(
              "Writing lightcone maps in MPI collective mode requires parallel "
              "HDF5");
#endif
#endif
        /* Number of files to write */
        int nr_files_per_shell = collective ? 1 : comm_size;

        /* Create the output file(s) */
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        if (file_id < 0) error("Unable to create file %s", fname);

        /* Write header with metadata */
        hid_t header =
            H5Gcreate(file_id, "Shell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        io_write_attribute_i(header, "nr_files_per_shell", nr_files_per_shell);
        io_write_attribute_d(
            header, "comoving_inner_radius",
            props->shell[shell_nr].rmin * length_conversion_factor);
        io_write_attribute_d(
            header, "comoving_outer_radius",
            props->shell[shell_nr].rmax * length_conversion_factor);
        H5Gclose(header);

        /* Write the system of Units used in the snapshot */
        io_write_unit_system(file_id, snapshot_units, "Units");

        /* Write the system of Units used internally */
        io_write_unit_system(file_id, internal_units, "InternalCodeUnits");

        /* Write the lightcone maps for this shell */
        for (int map_nr = 0; map_nr < nr_maps; map_nr += 1)
          lightcone_map_write(&(props->shell[shell_nr].map[map_nr]), file_id,
                              props->map_type[map_nr].name, internal_units,
                              snapshot_units, collective,
                              props->maps_gzip_level, props->hdf5_chunk_size,
                              props->map_type[map_nr].compression);

        /* Close the file */
        H5Pclose(fapl_id);
        H5Fclose(file_id);

        /* Free the pixel data associated with this shell */
        for (int map_nr = 0; map_nr < nr_maps; map_nr += 1)
          lightcone_map_free_pixels(&(props->shell[shell_nr].map[map_nr]));

        /* Update status of this shell */
        props->shell[shell_nr].state = shell_complete;
      }
    }
  }

  if (props->verbose && engine_rank == 0 && num_shells_written > 0)
    message("lightcone %d: Writing completed lightcone shells took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("Need HDF5 to write out lightcone maps");
#endif
}

/**
 * @brief Deallocate lightcone data.
 *
 * @param props the #lightcone_props structure.
 *
 */
void lightcone_clean(struct lightcone_props *props) {

  /* Deallocate particle buffers */
  for (int i = 0; i < swift_type_count; i += 1) {
    if (props->use_type[i]) particle_buffer_free(&props->buffer[i]);
  }

  /* Free replication list, if we have one */
  if (props->have_replication_list)
    replication_list_clean(&props->replication_list);

  /* Clean lightcone maps and free the structs */
  const int nr_shells = props->nr_shells;
  const int nr_maps = props->nr_maps;
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      lightcone_map_clean(&(props->shell[shell_nr].map[map_nr]));
    }
    free(props->shell[shell_nr].map);
  }

  /* Free buffers associated with each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      particle_buffer_free(&(props->shell[shell_nr].buffer[ptype]));
    }
  }

  /* Free array of shells */
  free(props->shell);

  /* Free array of lightcone map types */
  free(props->map_type);

  /* Free data associated with particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
    struct lightcone_particle_type *this_type = &(props->part_type[ptype]);
    free(this_type->map_index);
  }

  /* Free lists of output quantities */
  for (int ptype = 0; ptype < swift_type_count; ptype += 1)
    lightcone_io_field_list_clean(&props->particle_fields[ptype]);

  /* Free the projected kernel */
  projected_kernel_clean(&props->kernel_table);
}

/**
 * @brief Determine periodic copies of the simulation box which could
 * contribute to the lightcone.
 *
 *                     \
 *           \          \
 *            |         |
 * Obs      A |    B    | C
 *            |         |
 *           /          /
 *          R1         /
 *                    R0
 *
 * Consider a single particle being drifted. Here R0 is the comoving
 * distance to the time the particle is drifted FROM. R1 is the comoving
 * distance to the time the particle is drifted TO on this step.
 *
 * Particles which are beyond the lightcone surface at the start of
 * their drift (C) cannot cross the lightcone on this step if v < c.
 * Particles between the lightcone surfaces at the start and end of
 * their drift (B) may cross the lightcone (and certainly will if they
 * have zero velocity).
 *
 * Particles just within the lightcone surface at the start of their
 * drift (A) may be able to cross the lightcone due to their velocity so
 * we need to allow a boundary layer on the inside edge of the shell.
 * If we assume v < c, then we can use a layer of thickness R0-R1.
 *
 * Here we compute the earliest and latest times particles may be drifted
 * between, find the corresponding comoving distances R0 and R1, reduce
 * the inner distance by R0-R1, and find all periodic copies of the
 * simulation box which overlap this spherical shell.
 *
 * Later we use this list to know which periodic copies to check when
 * particles are drifted.
 *
 * This routine also determines which lightcone healpix maps might be
 * updated on this time step and allocates the pixel data if necessary.
 *
 * @param props The #lightcone_props structure
 * @param cosmo The #cosmology structure
 * @param ti_earliest_undrifted earliest integer time any particle might
 *        be drifted from on this step
 * @param ti_current End of the timestep
 *
 */
void lightcone_prepare_for_step(struct lightcone_props *props,
                                const struct cosmology *cosmo,
                                const integertime_t ti_earliest_undrifted,
                                const integertime_t ti_current) {
  ticks tic = getticks();

  /* Deallocate the old list, if there is one */
  if (props->have_replication_list)
    replication_list_clean(&props->replication_list);

  /* Get the size of the simulation box */
  const double boxsize = props->boxsize;

  /* Get a lower limit on earliest time particle may be drifted from */
  const integertime_t ti_lim = ti_earliest_undrifted;

  /* Get expansion factor at earliest and latest times particles might be
   * drifted between */
  double a_current = cosmo->a_begin * exp(ti_current * cosmo->time_base);
  double a_old = cosmo->a_begin * exp(ti_lim * cosmo->time_base);
  if (a_old < cosmo->a_begin) a_old = cosmo->a_begin;

  /* Convert redshift range to a distance range */
  double lightcone_rmin = cosmology_get_comoving_distance(cosmo, a_current);
  double lightcone_rmax = cosmology_get_comoving_distance(cosmo, a_old);
  if (lightcone_rmin > lightcone_rmax) error("Lightcone has rmin > rmax");

  /* Allow inner boundary layer, assuming all particles have v < c.
     This is to account for particles moving during the time step. */
  double boundary = lightcone_rmax - lightcone_rmin;
  lightcone_rmin -= boundary;
  if (lightcone_rmin < 0) lightcone_rmin = 0;

  if (a_current < props->a_min || a_old > props->a_max) {
    /* Timestep does not overlap the lightcone redshift range */
    replication_list_init_empty(&props->replication_list);
  } else {
    /* Timestep may contribute particles to the lightcone */
    replication_list_init(&props->replication_list, boxsize, props->cell_width,
                          props->observer_position, lightcone_rmin,
                          lightcone_rmax);
  }

  /* Record that we made the list */
  props->have_replication_list = 1;

  /* Store times we used to make the list, for consistency check later */
  props->ti_old = ti_lim;
  props->ti_current = ti_current;

  /* Report the size of the list */
#ifdef DUMP_REPLICATIONS
  if (engine_rank == 0) {
    message("lightcone %d: no. of replications to check: %d", props->index,
            props->replication_list.nrep);
    message("lightcone %d: shell to search inner radius=%e, outer radius=%e",
            props->index, lightcone_rmin, lightcone_rmax);
  }
#endif

  /* Write out the list, if required */
#ifdef DUMP_REPLICATIONS
  if (engine_rank == 0) {
    char fname[500];
    sprintf(fname, "replication_list.%d.txt", output_nr);
    FILE *fd_rep = fopen(fname, "w");
    fprintf(fd_rep, "# Observer x, y, z\n");
    fprintf(fd_rep, "%e, %e, %e\n", props->observer_position[0],
            props->observer_position[1], props->observer_position[2]);
    fprintf(fd_rep, "# Box size, inner radius, outer radius\n");
    fprintf(fd_rep, "%e, %e, %e\n", boxsize, lightcone_rmin - boundary,
            lightcone_rmax);
    fprintf(fd_rep, "# x, y, z, rmin2, rmax2\n");
    replication_list_write(&props->replication_list, fd_rep);
    fclose(fd_rep);
    output_nr += 1;
  }
#endif

  /* Number of shells and maps per shell */
  const int nr_maps = props->nr_maps;
  const int nr_shells = props->nr_shells;

  /* Range of shells that might be updated this step */
  int shell_nr_min = nr_shells;
  int shell_nr_max = -1;

  /* Loop over healpix map shells */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {

    const double shell_amin = props->shell[shell_nr].amin;
    const double shell_amax = props->shell[shell_nr].amax;
    const double step_amin = a_old;
    const double step_amax = a_current;

    /* Check if this shell might be updated */
    if (step_amin <= shell_amax && step_amax >= shell_amin) {

      switch (props->shell[shell_nr].state) {
        case shell_uninitialized:
          /* This shell has not been allocated yet, so allocate it */
          if (props->verbose && engine_rank == 0)
            message("lightcone %d: allocating pixels for shell %d at a=%f",
                    props->index, shell_nr, cosmo->a);
          for (int map_nr = 0; map_nr < nr_maps; map_nr += 1)
            lightcone_map_allocate_pixels(&(props->shell[shell_nr].map[map_nr]),
                                          /* zero_pixels = */ 1);
          props->shell[shell_nr].state = shell_current;
          break;
        case shell_complete:
          /* Shell has already been written out and freed - should never happen
           */
          error(
              "Lightcone shell has been written out while particles could "
              "still contribute");
          break;
        case shell_current:
          /* Already initialized, nothing to do */
          break;
      }

      /* Record range of shells that might be updated this step */
      if (shell_nr < shell_nr_min) shell_nr_min = shell_nr;
      if (shell_nr > shell_nr_max) shell_nr_max = shell_nr;
    }
  }
  props->shell_nr_min = shell_nr_min;
  props->shell_nr_max = shell_nr_max;

  /* Determine which particle types might contribute to lightcone outputs at
   * this step */
  for (int i = 0; i < swift_type_count; i += 1) {
    props->check_type_for_crossing[i] = 0;
    if (props->a_max_search_for_type[i] >= props->a_min_search_for_type[i]) {
      if (a_current >= props->a_min_search_for_type[i] &&
          a_old <= props->a_max_search_for_type[i])
        props->check_type_for_crossing[i] = 1;
    }
  }

  if (props->verbose && engine_rank == 0)
    message("lightcone %d: Lightcone timestep preparations took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Determine whether lightcone map buffers should be flushed this step.
 *
 * @param props The #lightcone_props structure
 *
 */
int lightcone_trigger_map_update(struct lightcone_props *props) {

  size_t total_updates = 0;
  const int nr_shells = props->nr_shells;
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    if (props->shell[shell_nr].state == shell_current) {
      for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
        total_updates += particle_buffer_num_elements(
            &(props->shell[shell_nr].buffer[ptype]));
      }
    }
  }
  return total_updates >= ((size_t)props->max_updates_buffered);
}

/**
 * @brief Add a particle to the output buffer
 *
 * @param props The #lightcone_props structure
 * @param e The #engine structure
 * @param gp The #gpart to buffer
 * @param a_cross Expansion factor of lightcone crossing
 * @param x_cross Position of the gpart at lightcone crossing
 */
void lightcone_buffer_particle(struct lightcone_props *props,
                               const struct engine *e, const struct gpart *gp,
                               const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  switch (gp->type) {
    case swift_type_gas: {

      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];
      struct lightcone_gas_data data;
      if (lightcone_store_gas(e, props, gp, p, xp, a_cross, x_cross, &data))
        particle_buffer_append(props->buffer + swift_type_gas, &data);

    } break;

    case swift_type_stars: {

      const struct spart *sp = &sparts[-gp->id_or_neg_offset];
      struct lightcone_stars_data data;
      if (lightcone_store_stars(e, props, gp, sp, a_cross, x_cross, &data))
        particle_buffer_append(props->buffer + swift_type_stars, &data);

    } break;

    case swift_type_black_hole: {

      const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
      struct lightcone_black_hole_data data;
      if (lightcone_store_black_hole(e, props, gp, bp, a_cross, x_cross, &data))
        particle_buffer_append(props->buffer + swift_type_black_hole, &data);

    } break;

    case swift_type_dark_matter: {

      struct lightcone_dark_matter_data data;
      if (lightcone_store_dark_matter(e, props, gp, a_cross, x_cross, &data))
        particle_buffer_append(props->buffer + swift_type_dark_matter, &data);

    } break;

    case swift_type_dark_matter_background: {

      /* Assumed to have same properties as DM particles */
      struct lightcone_dark_matter_data data;
      if (lightcone_store_dark_matter(e, props, gp, a_cross, x_cross, &data))
        particle_buffer_append(
            props->buffer + swift_type_dark_matter_background, &data);

    } break;

    case swift_type_neutrino: {

      struct lightcone_neutrino_data data;
      if (lightcone_store_neutrino(e, props, gp, a_cross, x_cross, &data))
        particle_buffer_append(props->buffer + swift_type_neutrino, &data);

    } break;

    default:
      error("Particle type not supported in lightcones");
  }
}

#ifdef HAVE_CHEALPIX
/**
 * @brief Compute the angular smoothing length of a particle
 *
 * @param pos particle position vector relative to observer
 * @param hsml physical smoothing length of the particle
 *
 */
static double angular_smoothing_scale(const double *pos, const double hsml) {

  /* Compute distance to particle */
  double dist = 0;
  for (int i = 0; i < 3; i += 1) dist += pos[i] * pos[i];
  dist = sqrt(dist);

  /* Avoid trig call for small angles (accurate to about 0.3%) */
  if (dist > 10.0 * hsml)
    return hsml / dist;
  else
    return atan(hsml / dist);
}
#endif

/**
 * @brief Buffer a particle's contribution to the healpix map(s)
 *
 * @param props The #lightcone_props structure
 * @param e The #engine structure
 * @param gp The #gpart to buffer
 * @param a_cross Expansion factor of lightcone crossing
 * @param x_cross Position of the gpart at lightcone crossing
 *
 */
void lightcone_buffer_map_update(struct lightcone_props *props,
                                 const struct engine *e, const struct gpart *gp,
                                 const double a_cross,
                                 const double x_cross[3]) {
#ifdef HAVE_CHEALPIX

  /* Find information on healpix maps this particle type contributes to */
  const struct lightcone_particle_type *part_type_info =
      &(props->part_type[gp->type]);

  /* If this particle type contributes to no healpix maps, do nothing */
  if (part_type_info->nr_maps == 0) return;

  /* Get angular coordinates of the particle */
  double theta, phi;
  vec2ang(x_cross, &theta, &phi);

  /* Get angular size of the particle */
  double radius;
  if (gp->type == swift_type_gas) {
    const struct part *parts = e->s->parts;
    const struct part *p = &parts[-gp->id_or_neg_offset];
    radius = angular_smoothing_scale(x_cross, p->h);
  } else {
    radius = 0.0;
  }

  /* Loop over shells to update */
  for (int shell_nr = props->shell_nr_min; shell_nr <= props->shell_nr_max;
       shell_nr += 1) {
    if (a_cross > props->shell[shell_nr].amin &&
        a_cross <= props->shell[shell_nr].amax) {

      /* Make sure this shell is available for updating */
      if (props->shell[shell_nr].state == shell_uninitialized)
        error("Attempt to update shell which has not been allocated");
      if (props->shell[shell_nr].state == shell_complete)
        error("Attempt to update shell which has been written out");

      /* Allocate storage for updates and set particle coordinates and radius */
      union lightcone_map_buffer_entry *data =
          (union lightcone_map_buffer_entry *)malloc(
              part_type_info->buffer_element_size);
      data[0].i = angle_to_int(theta);
      data[1].i = angle_to_int(phi);
      data[2].f = radius;

      /* Loop over healpix maps which this particle type contributes to and find
       * values to add */
      for (int i = 0; i < part_type_info->nr_maps; i += 1) {
        int map_nr = part_type_info->map_index[i];
        /* The value to add to the map may need to be scaled to fit in a float
         */
        const double fac = props->map_type[map_nr].buffer_scale_factor;
        /* Fetch the value to add to the map */
        const double val =
            props->map_type[map_nr].update_map(e, props, gp, a_cross, x_cross);
        /* Store the scaled value */
        data[3 + i].f = fac * val;
#ifdef LIGHTCONE_MAP_CHECK_TOTAL
        /* Accumulate total quantity added to each map for consistency check */
        atomic_add_d(&props->shell[shell_nr].map[map_nr].total, val);
#endif
      }

      /* Buffer the updates */
      particle_buffer_append(&(props->shell[shell_nr].buffer[gp->type]), data);

      /* Free update info */
      free(data);
    }
  } /* Next shell */
#else
  error("Need HEALPix C API to make lightcones");
#endif
}

/**
 * @brief Compute memory used by lightcones on this rank
 *
 * @param props The #lightcone_props structure
 * @param particle_buffer_bytes returns bytes used to buffer particles
 * @param map_buffer bytes returns bytes used to buffer map updates
 * @param pixel_data_bytes returns bytes used to store map pixels
 *
 */
void lightcone_memory_use(struct lightcone_props *props,
                          size_t *particle_buffer_bytes,
                          size_t *map_buffer_bytes, size_t *pixel_data_bytes) {

  *particle_buffer_bytes = 0;
  *map_buffer_bytes = 0;
  *pixel_data_bytes = 0;

  /* Accumulate memory used by particle buffers - one buffer per particle type
   */
  for (int i = 0; i < swift_type_count; i += 1) {
    if (props->use_type[i])
      *particle_buffer_bytes += particle_buffer_memory_use(props->buffer + i);
  }

  /* Accumulate memory used by map update buffers and pixel data */
  const int nr_maps = props->nr_maps;
  const int nr_shells = props->nr_shells;
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {

    /* Healpix map updates - one buffer per particle type per shell */
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      *map_buffer_bytes +=
          particle_buffer_memory_use(&(props->shell[shell_nr].buffer[ptype]));
    }

    /* Pixel data - one buffer per map per shell */
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      struct lightcone_map *map = &(props->shell[shell_nr].map[map_nr]);
      if (map->data) *pixel_data_bytes += map->local_nr_pix * sizeof(double);
    }
  }
}

/**
 * @brief Write out number of files per rank for this lightcone
 *
 * @param props The #lightcone_props structure
 * @param internal_units swift internal unit system
 * @param snapshot_units swift snapshot unit system
 *
 */
void lightcone_write_index(struct lightcone_props *props,
                           const struct unit_system *internal_units,
                           const struct unit_system *snapshot_units) {
  int comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif

  /* Collect current file index on each rank */
  int *current_file_on_rank = (int *)malloc(sizeof(int) * comm_size);
#ifdef WITH_MPI
  MPI_Gather(&props->current_file, 1, MPI_INT, current_file_on_rank, 1, MPI_INT,
             0, MPI_COMM_WORLD);
#else
  current_file_on_rank[0] = props->current_file;
#endif

  if (engine_rank == 0) {

    /* Get conversion factor for shell radii */
    const double length_conversion_factor = units_conversion_factor(
        internal_units, snapshot_units, UNIT_CONV_LENGTH);

    /* Get the name of the index file */
    char fname[FILENAME_BUFFER_SIZE];
    check_snprintf(fname, FILENAME_BUFFER_SIZE, "%s/%s_index.hdf5",
                   props->subdir, props->basename);

    hid_t h_props = H5Pcreate(H5P_FILE_ACCESS);
    herr_t err = H5Pset_libver_bounds(h_props, HDF5_LOWEST_FILE_FORMAT_VERSION,
                                      HDF5_HIGHEST_FILE_FORMAT_VERSION);
    if (err < 0) error("Error setting the hdf5 API version");

    /* Create the file */
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, h_props);

    /* Write number of MPI ranks and number of files */
    hid_t group_id =
        H5Gcreate(file_id, "Lightcone", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    io_write_attribute_i(group_id, "nr_mpi_ranks", comm_size);
    io_write_attribute(group_id, "final_particle_file_on_rank", INT,
                       current_file_on_rank, comm_size);

    /* Write number of files the lightcone maps are distributed over */
    int nr_files_per_shell = props->distributed_maps ? comm_size : 1;
    io_write_attribute_i(group_id, "nr_files_per_shell", nr_files_per_shell);

    /* Write observer position and redshift limits */
    io_write_attribute(group_id, "observer_position", DOUBLE,
                       props->observer_position, 3);
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      char name[PARSER_MAX_LINE_SIZE];
      check_snprintf(name, PARSER_MAX_LINE_SIZE, "minimum_redshift_%s",
                     part_type_names[ptype]);
      io_write_attribute_d(group_id, name, props->z_min_for_type[ptype]);
      check_snprintf(name, PARSER_MAX_LINE_SIZE, "maximum_redshift_%s",
                     part_type_names[ptype]);
      io_write_attribute_d(group_id, name, props->z_max_for_type[ptype]);
    }

    /* Write the number of shells and their radii */
    const int nr_shells = props->nr_shells;
    io_write_attribute_i(group_id, "nr_shells", nr_shells);
    double *shell_inner_radii = (double *)malloc(sizeof(double) * nr_shells);
    double *shell_outer_radii = (double *)malloc(sizeof(double) * nr_shells);
    for (int i = 0; i < nr_shells; i += 1) {
      shell_inner_radii[i] = props->shell[i].rmin * length_conversion_factor;
      shell_outer_radii[i] = props->shell[i].rmax * length_conversion_factor;
    }
    io_write_attribute(group_id, "shell_inner_radii", DOUBLE, shell_inner_radii,
                       nr_shells);
    io_write_attribute(group_id, "shell_outer_radii", DOUBLE, shell_outer_radii,
                       nr_shells);
    free(shell_outer_radii);
    free(shell_inner_radii);

    H5Gclose(group_id);
    H5Fclose(file_id);
    H5Pclose(h_props);
  }

  free(current_file_on_rank);
}

/**
 * @brief Add the baseline value to a lightcone map
 *
 * @param c the #cosmology struct
 * @param props the properties of this lightcone
 * @param map the #lightcone_map structure
 */
void lightcone_map_set_baseline(const struct cosmology *c,
                                struct lightcone_props *props,
                                struct lightcone_map *map) {

  /* Nothing to do if there is no baseline function */
  if (map->type.baseline_func == NULL) return;

  /* Fetch the baseline value */
  double baseline_value = map->type.baseline_func(c, props, map);

  /* Add it to the map if necessary */
  if (baseline_value != 0.0) {
    for (pixel_index_t i = 0; i < map->local_nr_pix; i += 1) {
#ifdef LIGHTCONE_MAP_CHECK_TOTAL
      map->total += baseline_value;
#endif
      map->data[i] += baseline_value;
    }
  }
}
