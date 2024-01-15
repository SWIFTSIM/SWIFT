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
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "parser.h"
#include "space.h"
#include "zoom_region/zoom_init.h"

struct gravity_props {
  double r_s;              // Placeholder for scale radius
  double r_cut_max_ratio;  // Placeholder for maximum cutoff ratio
  int mesh_size;           // Placeholder for mesh size
};

// Function to initialize a dummy gravity_props structure
struct gravity_props create_dummy_gravity_props() {
  struct gravity_props dummy_props;

  /* Set up the dummy props */
  dummy_props.r_s = 1.25;
  dummy_props.r_cut_max_ratio = 4.5;
  dummy_props.mesh_size = 512;

  return dummy_props;
}

struct space *make_mock_space(struct space *s) {

  /* Define the members we need for the test. */
  s->dim = 1000;
  s->nr_gparts = 18;

  /* Allocate memory for the gparts. */
  s->gparts = (struct gpart *)malloc(s->nr_gparts * sizeof(struct gpart));

  /* Define the corners of the region */
  double cube_corners[8][3] = {
      {560, 560, 560}, {560, 640, 560}, {640, 560, 560}, {640, 640, 560},
      {560, 560, 640}, {560, 640, 640}, {640, 560, 640}, {640, 640, 640}};

  /* Loop over the gparts and set up baxckground and zoom particles. */
  for (int i = 0; i < s->nr_gparts; i++) {
    s->gparts[i].mass = 1.0;

    /* Handle background and zoom region particles differently. */
    if (i < 10) {
      /* Set background particles to be evenly spaced. */
      s->gparts[i].x[0] = s->dim / s->nr_gparts * i;
      s->gparts[i].x[1] = s->dim / s->nr_gparts * i;
      s->gparts[i].x[2] = s->dim / s->nr_gparts * i;
      s->gparts[i].type = 2;

    } else {
      /* Set zoom region particles to be at the corners of the region. */
      memcpy(s->gparts[i].x, cube_corners[i], 3 * sizeof(double));
      s->gparts[i].type = 1;
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
    s = make_mock_space(s);

    /* Create a dummy gravity_props structure. */
    struct gravity_props grav_props = create_dummy_gravity_props();

    /* Run the zoom_init function. */
    zoom_region_init(&param_file, s, &grav_props, 1);
  }
