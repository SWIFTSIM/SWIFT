/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "snaplist.h"

/* Local includes. */
#include "cosmology.h"
#include "error.h"
#include "restart.h"

/* Some standard headers. */
#include <string.h>

/**
 * @brief Read a file containing a list of time
 *
 * @param snaplist The #snaplist to fill.
 * @param filename The file to read.
 * @param cosmo The #cosmology model.
 * @param max_size The maximal size for the time array.
 */
void snaplist_read_file(struct snaplist *snaplist, const char* filename,
			struct cosmology *cosmo, size_t max_size) {
  /* initialize snaplist */
  snaplist->size = 0;
  snaplist->times = (double *) malloc(sizeof(double) * max_size);
  snaplist->max_size = max_size;
  
  /* Open file */
  FILE* file = fopen(filename, "r");
  if (file == NULL) error("Error opening file '%s'", filename);

  /* Read header */
  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  if((read = getline(&line, &len, file)) == -1)
    error("Unable to read header in file '%s'", filename);

  /* Remove end of line character */
  line[strcspn(line, "\n")] = 0;

  /* Find type of data in file */
  int type;
  if (!strcmp(line, "# Redshift"))
    type = SNAPLIST_REDSHIFT;
  else if (!strcmp(line, "# Time"))
    type = SNAPLIST_AGE;
  else if (!strcmp(line, "# Scale Factor"))
    type = SNAPLIST_SCALE_FACTOR;
  else
    error("Unable to interpret the header (%s) in file '%s'",
	  line, filename);

  if (!cosmo &&
      (type == SNAPLIST_SCALE_FACTOR || type == SNAPLIST_REDSHIFT))
    error("Unable to compute a redshift or a scale factor without cosmology. "
	  "Please change the header in '%s'", filename);
    

  /* Read file */
  while ((read = getline(&line, &len, file)) != -1) {

    /* Check data size */
    if (snaplist->size == max_size)
      error("Not enough memory to write the snaplist buffer. "
	    "Please decrease the total number of output required in '%s'.",
	    filename);

    double *time = &snaplist->times[snaplist->size];
    /* Write data to snaplist */
    if (sscanf(line, "%lf", time) != 1) {
      error(
            "Tried parsing double but found '%s' with illegal double "
            "characters in file '%s'.",
	    line, filename);
    }

    /* Transform input into correct time (e.g. ages or scale factor) */
    if (type == SNAPLIST_REDSHIFT)
      *time = cosmology_get_a_from_z(cosmo, *time);

    if (cosmo && type == SNAPLIST_AGE)
      *time = cosmology_get_scale_factor(cosmo, *time);

    /* Update size */
    snaplist->size += 1;
  }

  fclose(file);
}

/**
 * @brief Print a snaplist
 */
void snaplist_print(const struct snaplist *snaplist) {

  printf("/*\t Time Array\t */\n");
  for(size_t ind = 0; ind < snaplist->size; ind++) {
    printf("\t%lf\n", snaplist->times[ind]);
  }
}

/**
 * @brief Clean a snaplist
 */
void snaplist_clean(struct snaplist *snaplist) {
  free(snaplist->times);
}

/**
 * @brief Dump a snaplist in a restart file
 */
void snaplist_struct_dump(struct snaplist *list, FILE *stream) {
  restart_write_blocks(list, sizeof(struct snaplist), 1, stream, "snaplist", "snaplist struct");

  restart_write_blocks(list, list->max_size, sizeof(double), stream, "snaplist", "times");
}


/**
 * @brief Restore a snaplist from a restart file
 */
void snaplist_struct_restore(struct snaplist * list, FILE *stream) {
  restart_read_blocks(list, sizeof(struct snaplist), 1, stream, NULL, "snaplist struct");
  
  list->times = (double *) malloc(sizeof(double) * list->max_size);
  restart_read_blocks(list->times, list->max_size, sizeof(double), stream, NULL, "times");
}
