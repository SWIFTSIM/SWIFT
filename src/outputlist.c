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
#include "error.h"
#include "restart.h"

/* Some standard headers. */
#include <string.h>

/**
 * @brief Read a file containing a list of time
 *
 * @param outputlist The #outputlist to fill.
 * @param filename The file to read.
 * @param cosmo The #cosmology model.
 */
void outputlist_read_file(struct outputlist *outputlist, const char *filename,
                          struct cosmology *cosmo) {

  /* Open file */
  FILE *file = fopen(filename, "r");
  if (file == NULL) error("Error opening file '%s'", filename);

  /* Declare reading variables */
  ssize_t read;
  size_t len = 0;
  char *line = NULL;

  /* Count number of lines */
  size_t nber_line = -1; /* Do not count header */
  while (getline(&line, &len, file) != -1) {
    nber_line++;
  }
  outputlist->size = nber_line;

  /* Return to start of file and initialize time array */
  fseek(file, 0, SEEK_SET);
  outputlist->times = (double *)malloc(sizeof(double) * nber_line);
  if (!outputlist->times)
    error(
        "Unable to malloc outputlist time array. "
        "Try reducing the number of lines in %s",
        filename);

  /* Read header */
  if ((read = getline(&line, &len, file)) == -1)
    error("Unable to read header in file '%s'", filename);

  /* Remove end of line character */
  line[strcspn(line, "\n")] = 0;

  /* Find type of data in file */
  int type = -1;
  if (!strcmp(line, "# Redshift"))
    type = OUTPUTLIST_REDSHIFT;
  else if (!strcmp(line, "# Time"))
    type = OUTPUTLIST_AGE;
  else if (!strcmp(line, "# Scale Factor"))
    type = OUTPUTLIST_SCALE_FACTOR;
  else
    error("Unable to interpret the header (%s) in file '%s'", line, filename);

  if (!cosmo &&
      (type == OUTPUTLIST_SCALE_FACTOR || type == OUTPUTLIST_REDSHIFT))
    error(
        "Unable to compute a redshift or a scale factor without cosmology. "
        "Please change the header in '%s'",
        filename);

  /* Read file */
  size_t ind = 0;
  while ((read = getline(&line, &len, file)) != -1) {
    double *time = &outputlist->times[ind];
    /* Write data to outputlist */
    if (sscanf(line, "%lf", time) != 1) {
      error(
          "Tried parsing double but found '%s' with illegal double "
          "characters in file '%s'.",
          line, filename);
    }

    /* Transform input into correct time (e.g. ages or scale factor) */
    if (type == OUTPUTLIST_REDSHIFT) *time = 1. / (1. + *time);

    if (cosmo && type == OUTPUTLIST_AGE)
      *time = cosmology_get_scale_factor(cosmo, *time);

    /* Update size */
    ind += 1;
  }

  fclose(file);
}

/**
 * @brief Print an #outputlist
 */
void outputlist_print(const struct outputlist *outputlist) {

  printf("/*\t Time Array\t */\n");
  printf("Number of Line: %lu\n", outputlist->size);
  for (size_t ind = 0; ind < outputlist->size; ind++) {
    printf("\t%lf\n", outputlist->times[ind]);
  }
}

/**
 * @brief Clean an #outputlist
 */
void outputlist_clean(struct outputlist *outputlist) {
  free(outputlist->times);
}

/**
 * @brief Dump an #outputlist in a restart file
 */
void outputlist_struct_dump(struct outputlist *list, FILE *stream) {
  restart_write_blocks(list, sizeof(struct outputlist), 1, stream, "outputlist",
                       "outputlist struct");

  restart_write_blocks(list->times, list->size, sizeof(double), stream,
                       "outputlist", "times");
}

/**
 * @brief Restore an #outputlist from a restart file
 */
void outputlist_struct_restore(struct outputlist *list, FILE *stream) {
  restart_read_blocks(list, sizeof(struct outputlist), 1, stream, NULL,
                      "outputlist struct");

  list->times = (double *)malloc(sizeof(double) * list->size);
  restart_read_blocks(list->times, list->size, sizeof(double), stream, NULL,
                      "times");
}
