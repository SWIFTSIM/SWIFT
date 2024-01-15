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

/* Function to initialize a dummy gravity_props structure */
void create_dummy_gravity_props(struct gravity_props *props) {

  /* Set up the dummy props */
  props->r_s = 1.25;
  props->r_cut_max_ratio = 4.5;
  props->mesh_size = 512;
}

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
  for (int i = 0; i < s->nr_gparts; i++) {
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

    printf("gpart %d: (%f, %f, %f)\n", i, gparts[i].x[0], gparts[i].x[1],
           gparts[i].x[2]);

    s->gparts = gparts;
  }
}

int main(int argc, char *argv[]) {

  const char *input_file = argv[1];
  const int test_type = atoi(argv[2]);

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Create a space structure. */
  struct space *s = malloc(sizeof(struct space));
  make_mock_space(s);

  /* Create a dummy gravity_props structure. */
  struct gravity_props *grav_props = malloc(sizeof(struct gravity_props));
  create_dummy_gravity_props(grav_props);

  printf("Running zoom_init test.\n");

  /* Run the zoom_init function. */
  zoom_region_init(&param_file, s, (const struct gravity_props *)grav_props, 1);

  /* Test what we've calculated and ensure it's as expected */
  assert(s->nr_gparts == 18);
  assert(s->with_zoom_region == 1);
  assert(s->zoom_props->cdim[0] == 16);
  assert(s->zoom_props->cdim[1] == 16);
  assert(s->zoom_props->cdim[2] == 16);
  assert(s->zoom_props->region_bounds[0] + (s->zoom_props->dim[0] / 2) == 500);
  assert(s->zoom_props->region_bounds[1] + (s->zoom_props->dim[1] / 2) == 500);
  assert(s->zoom_props->region_bounds[2] + (s->zoom_props->dim[2] / 2) == 500);
  assert(s->zoom_props->width[0] ==
         s->zoom_props->dim[0] / s->zoom_props->cdim[0]);
  assert(s->zoom_props->width[1] ==
         s->zoom_props->dim[1] / s->zoom_props->cdim[1]);
  assert(s->zoom_props->width[2] ==
         s->zoom_props->dim[2] / s->zoom_props->cdim[2]);

  if (test_type == 0) {
    assert(s->zoom_props->dim[0] == 100);
    assert(s->zoom_props->dim[1] == 100);
    assert(s->zoom_props->dim[2] == 100);
  } else if (test_type == 1) {
    assert(s->zoom_props->dim[0] == 100);
    assert(s->zoom_props->dim[1] == 100);
    assert(s->zoom_props->dim[2] == 100);
  } else if (test_type == 2) {
    assert(s->zoom_props->dim[0] == 111.111111);
    assert(s->zoom_props->dim[1] == 111.111111);
    assert(s->zoom_props->dim[2] == 111.111111);
  }

  free(s);
  free(grav_props);

  return 0;
}
