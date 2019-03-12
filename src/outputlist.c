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
#include "outputlist.h"

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
 * @param outputlist The #output_list to fill.
 * @param filename The file to read.
 * @param cosmo The #cosmology model.
 */
void output_list_read_file(struct output_list *outputlist, const char *filename,
                           struct cosmology *cosmo) {

  /* Open file */
  FILE *file = fopen(filename, "r");
  if (file == NULL) error("Error opening file '%s'", filename);

  /* Count number of lines */
  size_t len = 0;
  char *line = NULL;
  size_t nber_line = 0;
  while (getline(&line, &len, file) != -1) nber_line++;

  outputlist->size = nber_line - 1; /* Do not count header */

  /* Return to start of file and initialize time array */
  fseek(file, 0, SEEK_SET);
  outputlist->times = (double *)malloc(sizeof(double) * outputlist->size);
  if (!outputlist->times)
    error(
        "Unable to malloc output_list. "
        "Try reducing the number of lines in %s",
        filename);

  /* Read header */
  if (getline(&line, &len, file) == -1)
    error("Unable to read header in file '%s'", filename);

  /* Remove end of line character */
  line[strcspn(line, "\n")] = 0;

  /* Find type of data in file */
  int type = -1;

  trim_trailing(line);
  if (strcasecmp(line, "# Redshift") == 0)
    type = OUTPUT_LIST_REDSHIFT;
  else if (strcasecmp(line, "# Time") == 0)
    type = OUTPUT_LIST_AGE;
  else if (strcasecmp(line, "# Scale Factor") == 0)
    type = OUTPUT_LIST_SCALE_FACTOR;
  else
    error("Unable to interpret the header (%s) in file '%s'", line, filename);

  if (!cosmo &&
      (type == OUTPUT_LIST_SCALE_FACTOR || type == OUTPUT_LIST_REDSHIFT))
    error(
        "Unable to compute a redshift or a scale factor without cosmology. "
        "Please change the header in '%s'",
        filename);

  /* Read file */
  size_t ind = 0;
  while (getline(&line, &len, file) != -1) {
    double *time = &outputlist->times[ind];
    /* Write data to outputlist */
    if (sscanf(line, "%lf", time) != 1)
      error(
          "Tried parsing double but found '%s' with illegal double "
          "characters in file '%s'.",
          line, filename);

    /* Transform input into correct time (e.g. ages or scale factor) */
    if (type == OUTPUT_LIST_REDSHIFT) *time = 1. / (1. + *time);

    if (cosmo && type == OUTPUT_LIST_AGE)
      *time = cosmology_get_scale_factor(cosmo, *time);

    /* Update size */
    ind += 1;
  }

  /* Cleanup */
  free(line);

  if (ind != outputlist->size)
    error("Did not read the correct number of output times.");

  /* Check that the list is in monotonic order */
  for (size_t i = 1; i < outputlist->size; ++i) {

    if ((type == OUTPUT_LIST_REDSHIFT) &&
        (outputlist->times[i] <= outputlist->times[i - 1]))
      error("Output list not having monotonically decreasing redshifts.");

    if ((type == OUTPUT_LIST_AGE) &&
        (outputlist->times[i] <= outputlist->times[i - 1]))
      error("Output list not having monotonically increasing ages.");

    if ((type == OUTPUT_LIST_SCALE_FACTOR) &&
        (outputlist->times[i] <= outputlist->times[i - 1]))
      error("Output list not having monotonically increasing scale-factors.");
  }

  /* set current indice to 0 */
  outputlist->cur_ind = 0;

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

  /* Find next snasphot above current time */
  double time = t->times[t->cur_ind];
  size_t ind = t->cur_ind;
  while (time < time_end) {

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

  /* Deal with last statistics */
  if (*ti_next >= max_nr_timesteps || ind == t->size || time >= time_end) {
    *ti_next = -1;
    if (e->verbose) message("No further output time for %s.", name);
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
}

/**
 * @brief initialize an output list
 *
 * @param list The output list to initialize
 * @param e The #engine
 * @param name The name of the section in params
 * @param delta_time updated to the initial delta time
 * @param time_first updated to the time of first output (scale factor or cosmic
 * time)
 */
void output_list_init(struct output_list **list, const struct engine *e,
                      const char *name, double *delta_time,
                      double *time_first) {
  struct swift_params *params = e->parameter_file;

  /* get cosmo */
  struct cosmology *cosmo = NULL;
  if (e->policy & engine_policy_cosmology) cosmo = e->cosmology;

  /* Read output on/off */
  char param_name[PARSER_MAX_LINE_SIZE];
  sprintf(param_name, "%s:output_list_on", name);
  int outputlist_on = parser_get_opt_param_int(params, param_name, 0);

  /* Check if read outputlist */
  if (!outputlist_on) return;

  /* Read outputlist for snapshots */
  *list = (struct output_list *)malloc(sizeof(struct output_list));

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
    *time_first = (*list)->times[0];
  } else {
    *delta_time = (*list)->times[1] - (*list)->times[0];
    *time_first = (*list)->times[0];
  }
}

/**
 * @brief Print an #output_list
 */
void output_list_print(const struct output_list *outputlist) {

  printf("/*\t Time Array\t */\n");
  printf("Number of Line: %lu\n", outputlist->size);
  for (size_t ind = 0; ind < outputlist->size; ind++) {
    if (ind == outputlist->cur_ind)
      printf("\t%lf <-- Current\n", outputlist->times[ind]);
    else
      printf("\t%lf\n", outputlist->times[ind]);
  }
}

/**
 * @brief Clean an #output_list
 */
void output_list_clean(struct output_list **outputlist) {
  if (*outputlist) {
    free((*outputlist)->times);
    free(*outputlist);
    *outputlist = NULL;
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
}
