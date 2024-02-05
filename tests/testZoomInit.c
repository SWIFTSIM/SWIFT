/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

/* Standard headers. */
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "parser.h"
#include "space.h"
#include "zoom_region/zoom_init.h"

void make_mock_space(struct space *s) {

  /* Define the members we need for the test. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;
  s->nr_gparts = 18;

  /* Allocate memory for the gparts. */
  struct gpart *gparts =
      (struct gpart *)malloc(s->nr_gparts * sizeof(struct gpart));

  /* Define the corners of the region */
  double cube_corners[8][3] = {
      {560, 560, 560}, {560, 640, 560}, {640, 560, 560}, {640, 640, 560},
      {560, 560, 640}, {560, 640, 640}, {640, 560, 640}, {640, 640, 640}};

  /* Loop over the gparts and set up baxckground and zoom particles. */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    gparts[i].mass = 1.0;

    /* Handle background and zoom region particles differently. */
    if (i < 10) {
      /* Set background particles to be evenly spaced. */
      gparts[i].x[0] = s->dim[0] / s->nr_gparts * i;
      gparts[i].x[1] = s->dim[1] / s->nr_gparts * i;
      gparts[i].x[2] = s->dim[2] / s->nr_gparts * i;
      gparts[i].type = swift_type_dark_matter_background;

    } else {
      /* Set zoom region particles to be at the corners of the region. */
      gparts[i].x[0] = cube_corners[i - 10][0];
      gparts[i].x[1] = cube_corners[i - 10][1];
      gparts[i].x[2] = cube_corners[i - 10][2];
      gparts[i].type = swift_type_dark_matter;
    }

    s->gparts = gparts;
  }
}

int main(int argc, char *argv[]) {

  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Create a space structure. */
  struct space *s = malloc(sizeof(struct space));
  make_mock_space(s);

  /* Run the zoom_init function. */
  zoom_region_init(&param_file, s, 0);

  /* Test what we've calculated and ensure the centre is in the centre of the
   * box. This ensures the dimensions, bounds and cdims have all been
   * calculated correctly. */
  assert(s->zoom_props->region_lower_bounds[0] + (s->zoom_props->dim[0] / 2) ==
         500);
  assert(s->zoom_props->region_lower_bounds[1] + (s->zoom_props->dim[1] / 2) ==
         500);
  assert(s->zoom_props->region_lower_bounds[2] + (s->zoom_props->dim[2] / 2) ==
         500);

  free(s);

  return 0;
}
