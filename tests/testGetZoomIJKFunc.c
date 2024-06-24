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
#include "cell.h"
#include "parser.h"
#include "space.h"
#include "swift.h"
#include "zoom_region/zoom.h"

void make_mock_space(struct space *s) {

  /* Define the members we need for the test. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;
  s->nr_gparts = 18;

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Allocate memory for the gparts. */
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0) {
    error("Failed to allocate memory for gparts");
  }
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Define the corners of the region */
  double cube_corners[8][3] = {{487.5, 487.5, 487.5}, {487.5, 512.5, 487.5},
                               {512.5, 487.5, 487.5}, {512.5, 512.5, 487.5},
                               {487.5, 487.5, 512.5}, {487.5, 512.5, 512.5},
                               {512.5, 487.5, 512.5}, {512.5, 512.5, 512.5}};

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

void make_mock_cells(struct space *s) {
  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Calculate the number of cells. */
  s->nr_cells =
      s->cdim[0] * s->cdim[1] * s->cdim[2] +
      zoom_props->cdim[0] * zoom_props->cdim[1] * zoom_props->cdim[2] +
      zoom_props->buffer_cdim[0] * zoom_props->buffer_cdim[1] *
          zoom_props->buffer_cdim[2];

  /* Allocate cells. */
  if (posix_memalign((void **)&s->cells_top, cell_align,
                     s->nr_cells * sizeof(struct cell)) != 0) {
    error("Couldn't allocate the cell");
  }
  bzero(s->cells_top, s->nr_cells * sizeof(struct cell));

  /* Get some zoom region properties */
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const int buffer_offset = zoom_props->buffer_cell_offset;
  const double zoom_bounds[3] = {zoom_props->region_lower_bounds[0],
                                 zoom_props->region_lower_bounds[1],
                                 zoom_props->region_lower_bounds[2]};
  const double buffer_bounds[3] = {zoom_props->buffer_lower_bounds[0],
                                   zoom_props->buffer_lower_bounds[1],
                                   zoom_props->buffer_lower_bounds[2]};

  /* Loop over zoom cells and set locations and initial values */
  for (int i = 0; i < zoom_props->cdim[0]; i++) {
    for (int j = 0; j < zoom_props->cdim[1]; j++) {
      for (int k = 0; k < zoom_props->cdim[2]; k++) {
        const size_t cid = cell_getid(zoom_props->cdim, i, j, k);
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->width[0]) + zoom_bounds[0];
        c->loc[1] = (j * zoom_props->width[1]) + zoom_bounds[1];
        c->loc[2] = (k * zoom_props->width[2]) + zoom_bounds[2];
        c->width[0] = zoom_props->width[0];
        c->width[1] = zoom_props->width[1];
        c->width[2] = zoom_props->width[2];

        c->type = cell_type_zoom;
        c->subtype = cell_subtype_regular;
      }
    }
  }

  /* Loop over natural cells and set locations and initial values */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        const size_t cid = cell_getid_offset(s->cdim, bkg_cell_offset, i, j, k);
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = i * s->width[0];
        c->loc[1] = j * s->width[1];
        c->loc[2] = k * s->width[2];
        c->width[0] = s->width[0];
        c->width[1] = s->width[1];
        c->width[2] = s->width[2];
        c->type = cell_type_bkg;
        c->subtype = cell_subtype_regular;
      }
    }
  }

  /* Loop over buffer cells and set locations and initial values */
  for (int i = 0; i < zoom_props->buffer_cdim[0]; i++) {
    for (int j = 0; j < zoom_props->buffer_cdim[1]; j++) {
      for (int k = 0; k < zoom_props->buffer_cdim[2]; k++) {
        const size_t cid =
            cell_getid_offset(zoom_props->buffer_cdim, buffer_offset, i, j, k);
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->buffer_width[0]) + buffer_bounds[0];
        c->loc[1] = (j * zoom_props->buffer_width[1]) + buffer_bounds[1];
        c->loc[2] = (k * zoom_props->buffer_width[2]) + buffer_bounds[2];
        c->width[0] = zoom_props->buffer_width[0];
        c->width[1] = zoom_props->buffer_width[1];
        c->width[2] = zoom_props->buffer_width[2];
        c->type = cell_type_buffer;
        c->subtype = cell_subtype_regular;
      }
    }
  }

  /* Label void cells. */
  for (int cid = zoom_props->buffer_cell_offset;
       cid < zoom_props->buffer_cell_offset + zoom_props->nr_buffer_cells;
       cid++) {

    /* Get the cell */
    struct cell *c = &s->cells_top[cid];

    /* Get the middle of the cell. */
    double mid[3] = {c->loc[0] + 0.5 * c->width[0],
                     c->loc[1] + 0.5 * c->width[1],
                     c->loc[2] + 0.5 * c->width[2]};

    /* Label this cell if it contains the zoom region. */
    if ((mid[0] > s->zoom_props->region_lower_bounds[0]) &&
        (mid[0] < s->zoom_props->region_upper_bounds[0]) &&
        (mid[1] > s->zoom_props->region_lower_bounds[1]) &&
        (mid[1] < s->zoom_props->region_upper_bounds[1]) &&
        (mid[2] > s->zoom_props->region_lower_bounds[2]) &&
        (mid[2] < s->zoom_props->region_upper_bounds[2])) {
      c->subtype = cell_subtype_void;
    }
  }

  /* Label the empty cells. */
  for (int cid = zoom_props->bkg_cell_offset;
       cid < zoom_props->bkg_cell_offset + zoom_props->nr_bkg_cells; cid++) {

    /* Get this cell. */
    struct cell *c = &s->cells_top[cid];

    /* Get the middle of the cell. */
    double mid[3] = {c->loc[0] + 0.5 * c->width[0],
                     c->loc[1] + 0.5 * c->width[1],
                     c->loc[2] + 0.5 * c->width[2]};

    /* Assign the cell type. */
    if ((mid[0] > s->zoom_props->buffer_lower_bounds[0]) &&
        (mid[0] < s->zoom_props->buffer_upper_bounds[0]) &&
        (mid[1] > s->zoom_props->buffer_lower_bounds[1]) &&
        (mid[1] < s->zoom_props->buffer_upper_bounds[1]) &&
        (mid[2] > s->zoom_props->buffer_lower_bounds[2]) &&
        (mid[2] < s->zoom_props->buffer_upper_bounds[2])) {
      c->subtype = cell_subtype_empty;
    }
  }
}

int main(int argc, char *argv[]) {

  /* Initialise CPU frequency, this also starts time. */
  unsigned long cpu_freq = 0;
  clocks_set_cpufreq(cpu_freq);

  /* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Create a space structure. */
  struct space s;
  bzero(&s, sizeof(struct space));
  make_mock_space(&s);

  /* Flag that we are running a zoom. */
  s.with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_props_init(&param_file, &s, 0);
  zoom_region_init(&s, 0);

  /* Make the cells. */
  make_mock_cells(&s);

  /* Define coordinates to test. */
  const double zoom_coords[3] = {500.0, 500.0, 500.0};
  const double buffer_coords[3] = {450.0, 550.0, 450.0};
  const double bkg_coords[3] = {10.0, 990.0, 800.0};

  /* Get a zoom, buffer and background cell using cell_getid functions. */
  const size_t zoom_cell_id =
      cell_getid_from_pos(&s, zoom_coords[0], zoom_coords[1], zoom_coords[2]);
  const size_t buffer_cell_id = cell_getid_from_pos(
      &s, buffer_coords[0], buffer_coords[1], buffer_coords[2]);
  const size_t bkg_cell_id =
      cell_getid_from_pos(&s, bkg_coords[0], bkg_coords[1], bkg_coords[2]);

  /* Get the cells. */
  const struct cell *zoom_cell = &s.cells_top[zoom_cell_id];
  const struct cell *buffer_cell = &s.cells_top[buffer_cell_id];
  const struct cell *bkg_cell = &s.cells_top[bkg_cell_id];

  /* Test that the right type of cell was returned. */
  assert(zoom_cell->type == cell_type_zoom);
  assert(buffer_cell->type == cell_type_buffer);
  assert(bkg_cell->type == cell_type_bkg);

  /* Test that the coordinate is actually inside the returned cell. */
  assert(zoom_coords[0] >= zoom_cell->loc[0] &&
         zoom_coords[0] < zoom_cell->loc[0] + zoom_cell->width[0]);
  assert(zoom_coords[1] >= zoom_cell->loc[1] &&
         zoom_coords[1] < zoom_cell->loc[1] + zoom_cell->width[1]);
  assert(zoom_coords[2] >= zoom_cell->loc[2] &&
         zoom_coords[2] < zoom_cell->loc[2] + zoom_cell->width[2]);
  assert(buffer_coords[0] >= buffer_cell->loc[0] &&
         buffer_coords[0] < buffer_cell->loc[0] + buffer_cell->width[0]);
  assert(buffer_coords[1] >= buffer_cell->loc[1] &&
         buffer_coords[1] < buffer_cell->loc[1] + buffer_cell->width[1]);
  assert(buffer_coords[2] >= buffer_cell->loc[2] &&
         buffer_coords[2] < buffer_cell->loc[2] + buffer_cell->width[2]);
  assert(bkg_coords[0] >= bkg_cell->loc[0] &&
         bkg_coords[0] < bkg_cell->loc[0] + bkg_cell->width[0]);
  assert(bkg_coords[1] >= bkg_cell->loc[1] &&
         bkg_coords[1] < bkg_cell->loc[1] + bkg_cell->width[1]);
  assert(bkg_coords[2] >= bkg_cell->loc[2] &&
         bkg_coords[2] < bkg_cell->loc[2] + bkg_cell->width[2]);

  free(s.cells_top);
  free(s.zoom_props);

  return 0;
}
