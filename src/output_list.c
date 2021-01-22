/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausamman (loic.hausammann@epfl.ch)
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
#include "output_list.h"

/* Local includes. */
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "restart.h"
#include "tools.h"

/* Some standard headers. */
#include <string.h>

/**
 * @brief Read a file containing a list of time
 *
 * @param output_list The #output_list to fill.
 * @param filename The file to read.
 * @param cosmo The #cosmology model.
 */
void output_list_read_file(struct output_list *output_list,
                           const char *filename, struct cosmology *cosmo) {

  /* Open file */
  FILE *file = fopen(filename, "r");
  if (file == NULL) error("Error opening file '%s'", filename);

  /* Count number of lines */
  size_t len = 0;
  char *line = NULL;
  size_t nber_line = 0;
  while (getline(&line, &len, file) != -1) nber_line++;

  output_list->size = nber_line - 1; /* Do not count header */

  /* Return to start of file and initialize time array */
  fseek(file, 0, SEEK_SET);
  output_list->times = (double *)malloc(sizeof(double) * output_list->size);
  output_list->snapshot_labels = (int *)malloc(sizeof(int) * output_list->size);
  output_list->select_output_indices =
      (int *)malloc(sizeof(int) * output_list->size);

  if ((!output_list->times) || (!output_list->select_output_indices)) {
    error(
        "Unable to malloc output_list. "
        "Try reducing the number of lines in %s",
        filename);
  }

  /* Read header */
  if (getline(&line, &len, file) == -1)
    error("Unable to read header in file '%s'", filename);

  /* Remove end of line character */
  line[strcspn(line, "\n")] = 0;

  /* Find type of data in file */
  int type = -1;
  output_list->select_output_on = 0;
  output_list->select_output_number_of_names = 0;
  output_list->alternative_labels_on = 0;

  trim_trailing(line);

  if (strcasecmp(line, "# Redshift") == 0) {
    type = OUTPUT_LIST_REDSHIFT;
  } else if (strcasecmp(line, "# Time") == 0) {
    type = OUTPUT_LIST_AGE;
  } else if (strcasecmp(line, "# Scale Factor") == 0) {
    type = OUTPUT_LIST_SCALE_FACTOR;
  } else if (strcasecmp(line, "# Redshift, Select Output") == 0) {
    type = OUTPUT_LIST_REDSHIFT;
    output_list->select_output_on = 1;
  } else if (strcasecmp(line, "# Time, Select Output") == 0) {
    type = OUTPUT_LIST_AGE;
    output_list->select_output_on = 1;
  } else if (strcasecmp(line, "# Scale Factor, Select Output") == 0) {
    type = OUTPUT_LIST_SCALE_FACTOR;
    output_list->select_output_on = 1;
  } else if (strcasecmp(line, "# Redshift, Select Output, Label") == 0) {
    type = OUTPUT_LIST_REDSHIFT;
    output_list->select_output_on = 1;
    output_list->alternative_labels_on = 1;
  } else if (strcasecmp(line, "# Time, Select Output, Label") == 0) {
    type = OUTPUT_LIST_AGE;
    output_list->select_output_on = 1;
    output_list->alternative_labels_on = 1;
  } else if (strcasecmp(line, "# Scale Factor, Select Output, Label") == 0) {
    type = OUTPUT_LIST_SCALE_FACTOR;
    output_list->select_output_on = 1;
    output_list->alternative_labels_on = 1;
  } else {
    error("Unable to interpret the header (%s) in file '%s'", line, filename);
  }

  if (!cosmo &&
      (type == OUTPUT_LIST_SCALE_FACTOR || type == OUTPUT_LIST_REDSHIFT))
    error(
        "Unable to compute a redshift or a scale factor without cosmology. "
        "Please change the header in '%s'",
        filename);

  if (!output_list->select_output_on && output_list->alternative_labels_on)
    error(
        "Found an output list with alternative labels but not individual "
        "output selections");

  /* Read file */
  size_t ind = 0;
  int read_successfully = 0;
  int found_select_output = 0;
  char select_output_buffer[FIELD_BUFFER_SIZE] =
      select_output_header_default_name;
  while (getline(&line, &len, file) != -1) {

    double *time = &output_list->times[ind];
    int *label = &output_list->snapshot_labels[ind];

    /* Write data to output_list */
    if (output_list->select_output_on && output_list->alternative_labels_on) {
      read_successfully = sscanf(line, "%lf, %[^,], %d", time,
                                 select_output_buffer, label) == 3;
    } else if (output_list->select_output_on) {
      read_successfully =
          sscanf(line, "%lf, %s", time, select_output_buffer) == 2;
    } else {
      read_successfully = sscanf(line, "%lf", time) == 1;
    }

    if (!read_successfully) {
      error(
          "Tried parsing output_list but found '%s' with illegal "
          "characters in file '%s'.",
          line, filename);
    }

    /* Transform input into correct time (e.g. ages or scale factor) */
    if (type == OUTPUT_LIST_REDSHIFT) *time = 1. / (1. + *time);

    if (cosmo && type == OUTPUT_LIST_AGE)
      *time = cosmology_get_scale_factor(cosmo, *time);

    /* Search to find index for select output - select_output_index is the index
     * in the select_output_names array that corresponds to this select output
     * name. */
    found_select_output = 0;
    for (int i = 0; i < output_list->select_output_number_of_names; i++) {

      if (!strcmp(select_output_buffer, output_list->select_output_names[i])) {
        /* We already have this select output list string in the buffer! */
        output_list->select_output_indices[ind] = i;
        found_select_output = 1;
      }
    }

    /* If we did not assign it above, we haven't encountered this name before
     * and we need to create this name in the array */
    if (!found_select_output) {
      strcpy(
          output_list
              ->select_output_names[output_list->select_output_number_of_names],
          select_output_buffer);
      output_list->select_output_indices[ind] =
          output_list->select_output_number_of_names;
      output_list->select_output_number_of_names += 1;
    }

    /* Update size */
    ind += 1;
  }

  /* Cleanup */
  free(line);

  if (ind != output_list->size)
    error("Did not read the correct number of output times.");

  /* Check that the list is in monotonic order */
  for (size_t i = 1; i < output_list->size; ++i) {

    if ((type == OUTPUT_LIST_REDSHIFT) &&
        (output_list->times[i] <= output_list->times[i - 1]))
      error("Output list not having monotonically decreasing redshifts.");

    if ((type == OUTPUT_LIST_AGE) &&
        (output_list->times[i] <= output_list->times[i - 1]))
      error("Output list not having monotonically increasing ages.");

    if ((type == OUTPUT_LIST_SCALE_FACTOR) &&
        (output_list->times[i] <= output_list->times[i - 1]))
      error("Output list not having monotonically increasing scale-factors.");
  }

  /* set current indice to 0 */
  output_list->cur_ind = 0;

  /* We check if this is true later */
  output_list->final_step_dump = 0;

  fclose(file);
}

/**
 * @brief Read the next time for an output
 *
 * @param t The #output_list
 * @param e The #engine.
 * @param name The name of the output (e.g. 'stats', 'snapshots', 'stf')
 * @param ti_next updated to the next output time
 */
void output_list_read_next_time(struct output_list *t, const struct engine *e,
                                const char *name, integertime_t *ti_next) {
  int is_cosmo = e->policy & engine_policy_cosmology;

  /* Find upper-bound on last output */
  double time_end;
  if (is_cosmo)
    time_end = e->cosmology->a_end;
  else
    time_end = e->time_end;

  /* Find next snapshot above current time */
  double time = t->times[t->cur_ind];
  size_t ind = t->cur_ind;
  while (time <= time_end) {

    /* Output time on the integer timeline */
    if (is_cosmo)
      *ti_next = log(time / e->cosmology->a_begin) / e->time_base;
    else
      *ti_next = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (*ti_next > e->ti_current) break;

    ind += 1;
    if (ind == t->size) break;

    time = t->times[ind];
    t->cur_ind = ind;
  }

  /* Do we need to do a dump at the end of the last timestep?
   * Note that what this really does is given that a t=t_max, z=0,
   * or a=1 is found in output_list.txt set the flag `final_step_dump`
   * to 1 - this is not special behaviour that is controlled by a
   * parameter file flag. */
  if (time == time_end ||
      (time > time_end && time - time_end < OUTPUT_LIST_EPS_TIME_END)) {
    t->final_step_dump = 1;
    if (e->verbose) {
      if (is_cosmo) {
        message("Next output time for %s set to a=%e.", name, time_end);
      } else {
        message("Next output time for %s set to t=%e.", name, time_end);
      }
    }
  }

  /* Deal with last statistics */
  if (*ti_next >= max_nr_timesteps || ind == t->size || time >= time_end) {
    *ti_next = -1;
    if (e->verbose && t->final_step_dump != 1) {
      message("No further output time for %s.", name);

      /* Do not print anything about the output style (below) */
      return;
    }
  } else {

    /* Be nice, talk... */
    if (is_cosmo) {
      const double next_time =
          exp(*ti_next * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next output time for %s set to a=%e.", name, next_time);
    } else {
      const double next_time = *ti_next * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next output time for %s set to t=%e.", name, next_time);
    }
  }

  /* Finally, talk if we are a snapshot and we are using SelectOutput */
  if (e->verbose) {
    if (t->select_output_number_of_names > 1) {
      char select_output_style[FIELD_BUFFER_SIZE];
      output_list_get_current_select_output(t, select_output_style);
      message("Next output style for %s set to %s.", name, select_output_style);
    }
  }
}

/**
 * @brief Copys the string describing the current output name into the
 *        buffer described by select_output_name.
 *
 * @param t The #output_list
 * @param select_output_name updated to the current Select Output choice.
 **/
void output_list_get_current_select_output(struct output_list *t,
                                           char *select_output_name) {
  strcpy(select_output_name,
         t->select_output_names[t->select_output_indices[t->cur_ind]]);
}

/**
 * @brief initialize an output list
 *
 * @param list The output list to initialize.
 * @param e The #engine.
 * @param name The name of the section in the param file.
 * @param delta_time (return) The delta between the first two outputs
 */
void output_list_init(struct output_list **list, const struct engine *e,
                      const char *name, double *const delta_time) {

  struct swift_params *params = e->parameter_file;

  if (*list != NULL) error("Output list already allocated!");

  /* Get cosmo */
  struct cosmology *cosmo = NULL;
  if (e->policy & engine_policy_cosmology) cosmo = e->cosmology;

  /* Read output on/off */
  char param_name[PARSER_MAX_LINE_SIZE];
  sprintf(param_name, "%s:output_list_on", name);
  const int output_list_on = parser_get_opt_param_int(params, param_name, 0);

  /* Check if read output_list */
  if (!output_list_on) return;

  /* Read output_list for snapshots */
  *list = (struct output_list *)malloc(sizeof(struct output_list));
  bzero(*list, sizeof(struct output_list));
  (*list)->output_list_on = output_list_on;

  /* Read filename */
  char filename[PARSER_MAX_LINE_SIZE];
  sprintf(param_name, "%s:output_list", name);
  parser_get_param_string(params, param_name, filename);

  if (e->verbose) message("Reading %s output file.", name);
  output_list_read_file(*list, filename, cosmo);

  if ((*list)->size < 2)
    error("You need to provide more snapshots in '%s'", filename);

  /* Set data for later checks */
  if (cosmo) {
    *delta_time = (*list)->times[1] / (*list)->times[0];
  } else {
    *delta_time = (*list)->times[1] - (*list)->times[0];
  }
}

/**
 * @brief Print an #output_list
 */
void output_list_print(const struct output_list *output_list) {

  printf("/*\t Total number of Select Output options: %d \t */\n",
         output_list->select_output_number_of_names);
  printf("/*\t Select Output options found in output list: \t */\n");
  for (int i = 0; i < output_list->select_output_number_of_names; i++) {
    printf("\t Index: %d, Name: %s\n", i, output_list->select_output_names[i]);
  }
  printf("\n/*\t Time Array, Select Output \t */\n");
  printf("Number of Lines: %zu\n", output_list->size);
  for (size_t ind = 0; ind < output_list->size; ind++) {
    if (ind == output_list->cur_ind)
      printf(
          "\t %lf, %s <-- Current\n", output_list->times[ind],
          output_list
              ->select_output_names[output_list->select_output_indices[ind]]);
    else
      printf(
          "\t %lf, %s\n", output_list->times[ind],
          output_list
              ->select_output_names[output_list->select_output_indices[ind]]);
  }
}

/**
 * @brief Clean an #output_list
 */
void output_list_clean(struct output_list **output_list) {
  if (*output_list) {
    free((*output_list)->times);
    free((*output_list)->snapshot_labels);
    free((*output_list)->select_output_indices);
    free(*output_list);
    *output_list = NULL;
  }
}

/**
 * @brief Dump an #output_list in a restart file
 */
void output_list_struct_dump(struct output_list *list, FILE *stream) {
  restart_write_blocks(list, sizeof(struct output_list), 1, stream,
                       "output_list", "output_list struct");

  restart_write_blocks(list->times, list->size, sizeof(double), stream,
                       "output_list", "times");

  restart_write_blocks(list->select_output_indices, list->size, sizeof(int),
                       stream, "output_list", "select_output_indices");
}

/**
 * @brief Restore an #output_list from a restart file
 */
void output_list_struct_restore(struct output_list *list, FILE *stream) {
  restart_read_blocks(list, sizeof(struct output_list), 1, stream, NULL,
                      "output_list struct");

  list->times = (double *)malloc(sizeof(double) * list->size);
  restart_read_blocks(list->times, list->size, sizeof(double), stream, NULL,
                      "times");

  list->select_output_indices = (int *)malloc(sizeof(int) * list->size);
  restart_read_blocks(list->select_output_indices, list->size, sizeof(int),
                      stream, NULL, "select_output_indices");
}
