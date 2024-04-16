/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2024 Will J. Roper (w.roper@sussex.ac.uk).
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
#include <fenv.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "parser.h"
#include "space.h"
#include "swift.h"
#include "zoom_region/zoom.h"

void make_mock_space(struct space *s, const double zoom_width) {

  /* Define the members we need for the test. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;
  s->nr_gparts = 18;

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Allocate memory for the gparts. */
  struct gpart *gparts =
      (struct gpart *)malloc(s->nr_gparts * sizeof(struct gpart));
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* Define the corners of the region */
  double cube_corners[8][3] = {
      {550 - (zoom_width / 2), 550 - (zoom_width / 2), 550 - (zoom_width / 2)},
      {550 - (zoom_width / 2), 550 + (zoom_width / 2), 550 - (zoom_width / 2)},
      {550 + (zoom_width / 2), 550 - (zoom_width / 2), 550 - (zoom_width / 2)},
      {550 + (zoom_width / 2), 550 + (zoom_width / 2), 550 - (zoom_width / 2)},
      {550 - (zoom_width / 2), 550 - (zoom_width / 2), 550 + (zoom_width / 2)},
      {550 - (zoom_width / 2), 550 + (zoom_width / 2), 550 + (zoom_width / 2)},
      {550 + (zoom_width / 2), 550 - (zoom_width / 2), 550 + (zoom_width / 2)},
      {550 + (zoom_width / 2), 550 + (zoom_width / 2), 550 + (zoom_width / 2)}};

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
  }

  s->gparts = gparts;
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const char *input_file = argv[1];
  const double zoom_width = atof(argv[2]);

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Create a space structure. */
  struct space *s = malloc(sizeof(struct space));
  bzero(s, sizeof(struct space));
  make_mock_space(s, zoom_width);

  /* Flag that we are running a zoom. */
  s->with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_props_init(&param_file, s, 0);
  zoom_region_init(s, 0);

  /* Test what we've calculated and ensure the centre is in the centre of the
   * box. This ensures the dimensions, bounds and cdims have all been
   * calculated correctly. */
  assert(s->zoom_props->region_lower_bounds[0] + (s->zoom_props->dim[0] / 2) ==
         500);
  assert(s->zoom_props->region_lower_bounds[1] + (s->zoom_props->dim[1] / 2) ==
         500);
  assert(s->zoom_props->region_lower_bounds[2] + (s->zoom_props->dim[2] / 2) ==
         500);

  /* Ensure the cell boundaries line up. */
  if (s->zoom_props->with_buffer_cells) {
    int found_bkg_bufferi_low = 0;
    int found_bkg_bufferj_low = 0;
    int found_bkg_bufferk_low = 0;
    int found_bkg_bufferi_up = 0;
    int found_bkg_bufferj_up = 0;
    int found_bkg_bufferk_up = 0;
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {

          if (s->width[0] * i == s->zoom_props->buffer_lower_bounds[0])
            found_bkg_bufferi_low = 1;
          if (s->width[1] * j == s->zoom_props->buffer_lower_bounds[1])
            found_bkg_bufferj_low = 1;
          if (s->width[2] * k == s->zoom_props->buffer_lower_bounds[2])
            found_bkg_bufferk_low = 1;

          if (s->width[0] * i == s->zoom_props->buffer_upper_bounds[0])
            found_bkg_bufferi_up = 1;
          if (s->width[1] * j == s->zoom_props->buffer_upper_bounds[1])
            found_bkg_bufferj_up = 1;
          if (s->width[2] * k == s->zoom_props->buffer_upper_bounds[2])
            found_bkg_bufferk_up = 1;
        }
      }
    }
    /* Did we find the boundaries?. */
    assert(found_bkg_bufferi_low && found_bkg_bufferj_low &&
           found_bkg_bufferk_low && found_bkg_bufferi_up &&
           found_bkg_bufferj_up && found_bkg_bufferk_up);

    /* And for the zoom cells. */
    int found_buffer_zoomi_low = 0;
    int found_buffer_zoomj_low = 0;
    int found_buffer_zoomk_low = 0;
    int found_buffer_zoomi_up = 0;
    int found_buffer_zoomj_up = 0;
    int found_buffer_zoomk_up = 0;
    for (int i = 0; i < s->zoom_props->buffer_cdim[0]; i++) {
      for (int j = 0; j < s->zoom_props->buffer_cdim[1]; j++) {
        for (int k = 0; k < s->zoom_props->buffer_cdim[2]; k++) {
          if (s->zoom_props->buffer_lower_bounds[0] +
                  s->zoom_props->buffer_width[0] * i ==
              s->zoom_props->region_lower_bounds[0])
            found_buffer_zoomi_low = 1;
          if (s->zoom_props->buffer_lower_bounds[0] +
                  s->zoom_props->buffer_width[1] * j ==
              s->zoom_props->region_lower_bounds[1])
            found_buffer_zoomj_low = 1;
          if (s->zoom_props->buffer_lower_bounds[0] +
                  s->zoom_props->buffer_width[2] * k ==
              s->zoom_props->region_lower_bounds[2])
            found_buffer_zoomk_low = 1;

          if (s->zoom_props->buffer_lower_bounds[0] +
                  s->zoom_props->buffer_width[0] * i ==
              s->zoom_props->region_upper_bounds[0])
            found_buffer_zoomi_up = 1;
          if (s->zoom_props->buffer_lower_bounds[0] +
                  s->zoom_props->buffer_width[1] * j ==
              s->zoom_props->region_upper_bounds[1])
            found_buffer_zoomj_up = 1;
          if (s->zoom_props->buffer_lower_bounds[0] +
                  s->zoom_props->buffer_width[2] * k ==
              s->zoom_props->region_upper_bounds[2])
            found_buffer_zoomk_up = 1;
        }
      }
    }
    /* Did we find the boundaries?. */
    assert(found_buffer_zoomi_low && found_buffer_zoomj_low &&
           found_buffer_zoomk_low && found_buffer_zoomi_up &&
           found_buffer_zoomj_up && found_buffer_zoomk_up);
  } else {

    /* Check the background and zoom cells align. */
    int found_bkg_zoomi_low = 0;
    int found_bkg_zoomj_low = 0;
    int found_bkg_zoomk_low = 0;
    int found_bkg_zoomi_up = 0;
    int found_bkg_zoomj_up = 0;
    int found_bkg_zoomk_up = 0;
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {
          if (s->width[0] * i == s->zoom_props->region_lower_bounds[0])
            found_bkg_zoomi_low = 1;
          if (s->width[1] * j == s->zoom_props->region_lower_bounds[1])
            found_bkg_zoomj_low = 1;
          if (s->width[2] * k == s->zoom_props->region_lower_bounds[2])
            found_bkg_zoomk_low = 1;

          if (s->width[0] * i == s->zoom_props->region_upper_bounds[0])
            found_bkg_zoomi_up = 1;
          if (s->width[1] * j == s->zoom_props->region_upper_bounds[1])
            found_bkg_zoomj_up = 1;
          if (s->width[2] * k == s->zoom_props->region_upper_bounds[2])
            found_bkg_zoomk_up = 1;
        }
      }
    }
    /* Did we find the boundaries?. */
    assert(found_bkg_zoomi_low && found_bkg_zoomj_low && found_bkg_zoomk_low &&
           found_bkg_zoomi_up && found_bkg_zoomj_up && found_bkg_zoomk_up);
  }

  free(s->gparts);
  free(s->zoom_props);
  free(s);

  return 0;
}
