/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "common_io.h"

/* Local includes. */
#include "error.h"
#include "io_properties.h"
#include "output_options.h"
#include "units.h"

/* Some standard headers. */
#include <string.h>

/**
 * @brief Return the particle type code of a select_output parameter
 *
 * @param name The name of the parameter under consideration.
 *
 * @return The (integer) particle type of the parameter.
 */
int io_get_param_ptype(const char* name) {

  const int name_len = strlen(name);

  for (int ptype = 0; ptype < swift_type_count; ptype++) {
    const int ptype_name_len = strlen(part_type_names[ptype]);
    if (name_len >= ptype_name_len &&
        strcmp(&name[name_len - ptype_name_len], part_type_names[ptype]) == 0)
      return ptype;
  }

  /* If we get here, we could not match the name, so something's gone wrong. */
  error("Could not determine the particle type for parameter '%s'.", name);

  /* We can never get here, but the compiler may complain if we don't return
   * an int after promising to do so... */
  return -1;
}

/**
 * @brief Return the number and names of all output fields of a given ptype.
 *
 * @param ptype The index of the particle type under consideration.
 * @param list An io_props list that will hold the individual fields.
 * @param with_cosmology Use cosmological name variant?
 * @param with_fof Include FoF related fields?
 * @param with_stf Include STF related fields?
 *
 * @return The total number of fields that can be written for the ptype.
 */
int io_get_ptype_fields(const int ptype, struct io_props* list,
                        const int with_cosmology, const int with_fof,
                        const int with_stf) {

  int num_fields = 0;

  switch (ptype) {

    case swift_type_gas:
      io_select_hydro_fields(NULL, NULL, with_cosmology, /*with_cooling=*/1,
                             /*with_temperature=*/1, with_fof, with_stf,
                             /*with_rt=*/1, /*e=*/NULL, &num_fields, list);
      break;

    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
    case swift_type_neutrino:
      io_select_dm_fields(NULL, NULL, with_fof, with_stf, /*e=*/NULL,
                          &num_fields, list);
      break;

    case swift_type_stars:
      io_select_star_fields(NULL, with_cosmology, with_fof, with_stf,
                            /*with_rt=*/1,
                            /*e=*/NULL, &num_fields, list);
      break;

    case swift_type_sink:
      io_select_sink_fields(NULL, with_cosmology, with_fof, with_stf,
                            /*e=*/NULL, &num_fields, list);
      break;

    case swift_type_black_hole:
      io_select_bh_fields(NULL, with_cosmology, with_fof, with_stf, /*e=*/NULL,
                          &num_fields, list);
      break;

    default:
      error("Particle Type %d not yet supported. Aborting", ptype);
  }

  return num_fields;
}

/**
 * @brief Prepare the output option fields according to the user's choices and
 * verify that they are valid.
 *
 * @param output_options The #output_options for this run
 * @param with_cosmology Ran with cosmology?
 * @param with_fof Are we running with on-the-fly Fof?
 * @param with_stf Are we running with on-the-fly structure finder?
 * @param verbose The verbose level
 */
void io_prepare_output_fields(struct output_options* output_options,
                              const int with_cosmology, const int with_fof,
                              const int with_stf, int verbose) {

  const int MAX_NUM_PTYPE_FIELDS = 100;

  /* Parameter struct for the output options */
  struct swift_params* params = output_options->select_output;

  /* Get all possible outputs per particle type */
  int ptype_num_fields_total[swift_type_count] = {0};
  struct io_props field_list[swift_type_count][MAX_NUM_PTYPE_FIELDS];

  for (int ptype = 0; ptype < swift_type_count; ptype++)
    ptype_num_fields_total[ptype] = io_get_ptype_fields(
        ptype, field_list[ptype], with_cosmology, with_fof, with_stf);

  /* Check for whether we have a `Default` section */
  int have_default = 0;

  /* Loop over each section, i.e. different class of output */
  for (int section_id = 0; section_id < params->sectionCount; section_id++) {

    /* Get the name of current (selection) section, without a trailing colon */
    char section_name[FIELD_BUFFER_SIZE];
    strcpy(section_name, params->section[section_id].name);
    section_name[strlen(section_name) - 1] = 0;

    /* Is this the `Default` section? */
    if (strcmp(section_name, select_output_header_default_name) == 0)
      have_default = 1;

    /* How many fields should each ptype write by default? */
    int ptype_num_fields_to_write[swift_type_count];

    /* What is the default writing status for each ptype (on/off)? */
    int ptype_default_write_status[swift_type_count];

    /* Initialise section-specific writing counters for each particle type.
     * If default is 'write', then we start from the total to deduct any fields
     * that are switched off. If the default is 'off', we have to start from
     * zero and then count upwards for each field that is switched back on. */
    for (int ptype = 0; ptype < swift_type_count; ptype++) {

      /* Internally also verifies that the default level is allowed */
      const enum lossy_compression_schemes compression_level_current_default =
          output_options_get_ptype_default_compression(
              params, section_name, (enum part_type)ptype, verbose);

      if (compression_level_current_default == compression_do_not_write) {
        ptype_default_write_status[ptype] = 0;
        ptype_num_fields_to_write[ptype] = 0;
      } else {
        ptype_default_write_status[ptype] = 1;
        ptype_num_fields_to_write[ptype] = ptype_num_fields_total[ptype];
      }

    } /* ends loop over particle types */

    /* Loop over each parameter */
    for (int param_id = 0; param_id < params->paramCount; param_id++) {

      /* Full name of the parameter to check */
      const char* param_name = params->data[param_id].name;

      /* Check whether the file still contains the old, now inappropriate
       * 'SelectOutput' section */
      if (strstr(param_name, "SelectOutput:") != NULL) {
        error(
            "Output selection files no longer require the use of top level "
            "SelectOutput; see the documentation for changes.");
      }

      /* Skip if the parameter belongs to another output class or is a
       * 'Standard' parameter */
      if (strstr(param_name, section_name) == NULL) continue;
      if (strstr(param_name, ":Standard_") != NULL) continue;
      if (strstr(param_name, ":basename") != NULL) continue;
      if (strstr(param_name, ":subdir") != NULL) continue;

      /* Get the particle type for current parameter
       * (raises an error if it could not determine it) */
      const int param_ptype = io_get_param_ptype(param_name);

      /* Issue a warning if this parameter does not pertain to any of the
       * known fields from this ptype. */
      int field_id = 0;
      char field_name[PARSER_MAX_LINE_SIZE];
      for (field_id = 0; field_id < ptype_num_fields_total[param_ptype];
           field_id++) {

        sprintf(field_name, "%s:%.*s_%s", section_name, FIELD_BUFFER_SIZE,
                field_list[param_ptype][field_id].name,
                part_type_names[param_ptype]);

        if (strcmp(param_name, field_name) == 0) break;
      }

      int param_is_known = 0; /* Update below if it is a known one */
      if (field_id < ptype_num_fields_total[param_ptype])
        param_is_known = 1;
      else
        message(
            "WARNING: Trying to change behaviour of field '%s' (read from "
            "'%s') that does not exist. This may be because you are not "
            "running with all of the physics that you compiled the code with.",
            param_name, params->fileName);

      /* Perform a correctness check on the _value_ of the parameter */
      char param_value[FIELD_BUFFER_SIZE];
      parser_get_param_string(params, param_name, param_value);

      int value_id = 0;
      for (value_id = 0; value_id < compression_level_count; value_id++)
        if (strcmp(param_value, lossy_compression_schemes_names[value_id]) == 0)
          break;

      if (value_id == compression_level_count)
        error("Choice of output selection parameter %s ('%s') is invalid.",
              param_name, param_value);

      /* Adjust number of fields to be written for param_ptype, if this field's
       * status is different from default and it is a known one. */
      if (param_is_known) {
        const int is_on =
            strcmp(param_value,
                   lossy_compression_schemes_names[compression_do_not_write]) !=
            0;

        if (is_on && !ptype_default_write_status[param_ptype]) {
          /* Particle should be written even though default is off:
           * increase field count */
          ptype_num_fields_to_write[param_ptype] += 1;
        }
        if (!is_on && ptype_default_write_status[param_ptype]) {
          /* Particle should not be written, even though default is on:
           * decrease field count */
          ptype_num_fields_to_write[param_ptype] -= 1;
        }
      }
    } /* ends loop over parameters */

    /* Second loop over ptypes, to write out total number of fields to write */
    for (int ptype = 0; ptype < swift_type_count; ptype++) {

#ifdef SWIFT_DEBUG_CHECKS
      /* Sanity check: is the number of fields to write non-negative? */
      if (ptype_num_fields_to_write[ptype] < 0)
        error(
            "We seem to have subtracted too many fields for particle "
            "type %d in output class %s (total to write is %d)",
            ptype, section_name, ptype_num_fields_to_write[ptype]);
#endif
      output_options->num_fields_to_write[section_id][ptype] =
          ptype_num_fields_to_write[ptype];
    }
  } /* Ends loop over sections, for different output classes */

  /* Add field numbers for (possible) implicit `Default` output class */
  if (!have_default) {
    const int default_id = output_options->select_output->sectionCount;
    for (int ptype = 0; ptype < swift_type_count; ptype++)
      output_options->num_fields_to_write[default_id][ptype] =
          ptype_num_fields_total[ptype];
  }
}

/**
 * @brief Write the output field parameters file
 *
 * @param filename The file to write.
 * @param with_cosmology Use cosmological name variant?
 * @param with_fof Use fof?
 * @param with_stf Using Velociraptor STF?
 */
void io_write_output_field_parameter(const char* filename, int with_cosmology,
                                     int with_fof, int with_stf) {

  FILE* file = fopen(filename, "w");
  if (file == NULL) error("Error opening file '%s'", filename);

  /* Create a fake unit system for the snapshots */
  struct unit_system snapshot_units;
  units_init_cgs(&snapshot_units);

  /* Loop over all particle types */
  fprintf(file, "Default:\n");
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    struct io_props list[100];
    int num_fields =
        io_get_ptype_fields(ptype, list, with_cosmology, with_fof, with_stf);

    if (num_fields == 0) continue;

    /* Output a header for that particle type */
    fprintf(file, "  # Particle Type %s\n", part_type_names[ptype]);

    /* Write all the fields of this particle type */
    for (int i = 0; i < num_fields; ++i) {

      char unit_buffer[FIELD_BUFFER_SIZE] = {0};
      units_cgs_conversion_string(unit_buffer, &snapshot_units, list[i].units,
                                  list[i].scale_factor_exponent);

      /* Need to buffer with a maximal size - otherwise we can't read in again
       * because comments are too long */
      char comment_write_buffer[PARSER_MAX_LINE_SIZE / 2];

      sprintf(comment_write_buffer, "%.*s", PARSER_MAX_LINE_SIZE / 2 - 1,
              list[i].description);

      /* If our string is too long, replace the last few characters (before
       * \0) with ... for 'fancy printing' */
      if (strlen(comment_write_buffer) > PARSER_MAX_LINE_SIZE / 2 - 3) {
        strcpy(&comment_write_buffer[PARSER_MAX_LINE_SIZE / 2 - 4], "...");
      }

      fprintf(file, "  %s_%s: %s  # %s : %s\n", list[i].name,
              part_type_names[ptype], "on", comment_write_buffer, unit_buffer);
    }

    fprintf(file, "\n");
  }

  fclose(file);

  printf(
      "List of valid ouput fields for the particle in snapshots dumped in "
      "'%s'.\n",
      filename);
}
