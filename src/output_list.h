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
#ifndef SWIFT_OUTPUT_LIST_H
#define SWIFT_OUTPUT_LIST_H

/* Config parameters. */
#include <config.h>

/* Local headers */
#include "common_io.h"
#include "timeline.h"

/* Pre-declarations */
struct cosmology;
struct engine;

/*! Maximal number of output lists */
#define OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES 8

/*! How far past time_end do we check for an output time? */
#define OUTPUT_LIST_EPS_TIME_END 1e-4

/**
 * @brief the different output_list type
 */
enum output_list_type {
  OUTPUT_LIST_AGE,
  OUTPUT_LIST_REDSHIFT,
  OUTPUT_LIST_SCALE_FACTOR,
  OUTPUT_LIST_COMOVING_DISTANCE,
};

/**
 * @brief the array containing the output times
 */
struct output_list {

  /* Select output names. */
  char select_output_names[OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES]
                          [FIELD_BUFFER_SIZE];

  /* Time array */
  double *times;

  /* Select output indices - each index corresponds to a string
   * in select_output. Chosen to be this instead of an array of
   * pointers because of restarts. */
  int *select_output_indices;

  /* List of snapshot labels if not using the defaults */
  int *snapshot_labels;

  /* Total number of currently used select output names */
  int select_output_number_of_names;

  /* Size of the time array (i.e. number of outputs) */
  size_t size;

  /* Current index */
  size_t cur_ind;

  /* Was the Select Output option used? */
  int select_output_on;

  /* Are we using individual labels for the runs? */
  int alternative_labels_on;

  /* Is this output list activated? */
  int output_list_on;

  /* Dump on final timestep? */
  int final_step_dump;
};

void output_list_read_file(struct output_list *output_list,
                           const char *filename, struct cosmology *cosmo);
void output_list_read_next_time(struct output_list *t, const struct engine *e,
                                const char *name, integertime_t *ti_next);
void output_list_get_current_select_output(struct output_list *t,
                                           char *select_output_name);
void output_list_init(struct output_list **list, const struct engine *e,
                      const char *name, double *const delta_time);
void output_list_print(const struct output_list *output_list);
void output_list_clean(struct output_list **output_list);
void output_list_struct_dump(struct output_list *list, FILE *stream);
void output_list_struct_restore(struct output_list *list, FILE *stream);
void output_list_check_selection(const struct output_list *list,
                                 const struct output_options *output_options);

#endif /* SWIFT_OUTPUT_LIST_H */
