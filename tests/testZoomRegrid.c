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

/* Local headers. */
#include "swift.h"
#include "zoom_region/zoom.h"

/* This script tests the process of regridding a zoom simulation. It'll first
 * calculate the zoom geometry, construct the cell structures as would be done
 * in a run during a regrid, populates 1 particle per cell and then tests that
 * all cells (and different cell types) containing the expected number
 * particles. This ensures the regrid produces the expected cell structure. */

void make_mock_space(struct space *s) {

  /* Define the boxsize. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;

  /* Define the gpart count (this is a background particle per cell + the 8
   * corner high res particles used to define the region.) */
  s->nr_gparts = (2 * 10 * 10 * 10) + (16 * 16 * 16) + 8;

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Allocate memory for the gparts. */
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0) {
    error("Failed to allocate memory for gparts");
  }
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* Define the corners of the region */
  const double mid = 500;
  const double offset = 12.5;
  double cube_corners[8][3] = {{mid - offset, mid - offset, mid - offset},
                               {mid - offset, mid - offset, mid + offset},
                               {mid - offset, mid + offset, mid - offset},
                               {mid - offset, mid + offset, mid + offset},
                               {mid + offset, mid - offset, mid - offset},
                               {mid + offset, mid - offset, mid + offset},
                               {mid + offset, mid + offset, mid - offset},
                               {mid + offset, mid + offset, mid + offset}};

  /* Loop over the cell grids making a background particle per cell. */
  int gpart_count = 0;
  const double bkg_width = 100;
  const double buffer_width = 20;
  const double zoom_width = 2.5;
  const double buffer_lower[3] = {400, 400, 400};
  const double zoom_lower[3] = {480, 480, 480};
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      for (int k = 0; k < 10; k++) {
        gparts[gpart_count].mass = 1.0;

        /* Set background particles to be at the center of the cell. */
        gparts[gpart_count].x[0] = bkg_width * (i + 0.5);
        gparts[gpart_count].x[1] = bkg_width * (j + 0.5);
        gparts[gpart_count].x[2] = bkg_width * (k + 0.5);
        gparts[gpart_count++].type = swift_type_dark_matter_background;
      }
    }
  }
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      for (int k = 0; k < 10; k++) {
        gparts[gpart_count].mass = 1.0;

        /* Set buffer particles to be at the center of the cell. */
        gparts[gpart_count].x[0] = (buffer_width * (i + 0.5)) + buffer_lower[0];
        gparts[gpart_count].x[1] = (buffer_width * (j + 0.5)) + buffer_lower[1];
        gparts[gpart_count].x[2] = (buffer_width * (k + 0.5)) + buffer_lower[2];
        gparts[gpart_count++].type = swift_type_dark_matter_background;
      }
    }
  }
  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 16; j++) {
      for (int k = 0; k < 16; k++) {
        gparts[gpart_count].mass = 1.0;

        /* Set zoom particles to be at the center of the cell. */
        gparts[gpart_count].x[0] = (zoom_width * (i + 0.5)) + zoom_lower[0];
        gparts[gpart_count].x[1] = (zoom_width * (j + 0.5)) + zoom_lower[1];
        gparts[gpart_count].x[2] = (zoom_width * (k + 0.5)) + zoom_lower[2];
        gparts[gpart_count++].type = swift_type_dark_matter_background;
      }
    }
  }

  /* Make a zoom particle in each corner to define the zoom region. */
  for (size_t i = 0; i < 8; i++) {
    gparts[gpart_count].mass = 1.0;

    /* Set zoom region particles to be at the corners of the region. */
    gparts[gpart_count].x[0] = cube_corners[i][0];
    gparts[gpart_count].x[1] = cube_corners[i][1];
    gparts[gpart_count].x[2] = cube_corners[i][2];
    gparts[gpart_count++].type = swift_type_dark_matter;
  }

  s->gparts = gparts;
}

void associate_gparts_to_cells(struct space *s) {
  for (size_t i = 0; i < s->nr_gparts; i++) {
    struct gpart *gpart = &s->gparts[i];

    int cid = cell_getid_from_pos(s, gpart->x[0], gpart->x[1], gpart->x[2]);

    struct cell *c = &s->cells_top[cid];
    if (c == NULL) {
      error("Failed to find cell for gpart.");
    }
    c->grav.count++;
  }
}

void make_gparts_grid(struct space *s) {}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(argv[1], &param_file);

  /* Create a space structure. (this creates fake particles to get the layout
   * of cells right but we will modify this shortly to get a particle per
   * cell) */
  struct space *s = malloc(sizeof(struct space));
  bzero(s, sizeof(struct space));
  make_mock_space(s);

  /* Flag that we are running a zoom. */
  s->with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_props_init(&param_file, s, /*verbose*/ 0);

  /* Run the regridding. */
  space_regrid(s, /*verbose*/ 0);

  /* Associate gparts. */
  associate_gparts_to_cells(s);

  /* Test all cells that aren't void or empty have at least one particle. */
  for (int i = 0; i < s->nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    if (c->subtype != cell_subtype_void && c->subtype != cell_subtype_empty) {
      if (c->grav.count == 0) {
        error(
            "Cell %d has no particles (c->type = %s, c->subtype = %s, c->loc = "
            "[%f, %f, %f])",
            i, cellID_names[c->type], subcellID_names[c->subtype], c->loc[0],
            c->loc[1], c->loc[2]);
      }
    }

    /* If we have a background cell ensure the count is 1 (because of nesting
     * zoom and buffer cells can have more). */
    if (c->type == cell_type_bkg && c->subtype != cell_subtype_empty &&
        c->grav.count != 1) {
      error(
          "Cell %d has %d particles (c->type = %s, c->subtype = %s, c->loc = "
          "[%f, %f, %f])",
          i, c->grav.count, cellID_names[c->type], subcellID_names[c->subtype],
          c->loc[0], c->loc[1], c->loc[2]);
    }

    /* Ensure empty and void cells have 0 counts. */
    if (c->subtype == cell_subtype_empty || c->subtype == cell_subtype_void) {
      if (c->grav.count != 0) {
        error(
            "Cell %d has %d particles (c->type = %s, c->subtype = %s, c->loc = "
            "[%f, %f, %f])",
            i, c->grav.count, cellID_names[c->type],
            subcellID_names[c->subtype], c->loc[0], c->loc[1], c->loc[2]);
      }
    }
  }

  /* Test each cell type contains the expected counts. (This ensures the zoom
   * and void cells make sense too.)*/
  int zoom_count = 0;
  int bkg_count = 0;
  int buffer_count = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    if (c->type == cell_type_zoom) {
      zoom_count += c->grav.count;
    } else if (c->type == cell_type_bkg) {
      bkg_count += c->grav.count;
    } else if (c->type == cell_type_buffer) {
      buffer_count += c->grav.count;
    } else {
      error("Cell %d has unexpected type %s", i, cellID_names[c->type]);
    }
  }

  /* Zoom region should have zoomcdim^3 + 8 high res corners + 8 void cells. */
  assert(zoom_count == 16 * 16 * 16 + 8 + 8);

  /* Background cell should have bkgcdim^3 - 8 empty */
  assert(bkg_count == 10 * 10 * 10 - 8);

  /* Buffer cell should have buffercdim^3 - 8 void + 8 empty */
  assert(buffer_count == 10 * 10 * 10);

  /* Free the space. */
  free(s->local_cells_top);
  free(s->multipoles_top);
  free(s->local_cells_with_tasks_top);
  free(s->cells_with_particles_top);
  free(s->local_cells_with_particles_top);
  free(s->zoom_props->local_zoom_cells_top);
  free(s->zoom_props->void_cell_indices);
  free(s->zoom_props->neighbour_cells_top);
  free(s->zoom_props->local_bkg_cells_top);
  free(s->zoom_props->local_buffer_cells_top);
  free(s->zoom_props->local_zoom_cells_with_particles_top);
  free(s->zoom_props->local_bkg_cells_with_particles_top);
  free(s->zoom_props->local_buffer_cells_with_particles_top);
  free(s->cells_top);
  free(s->gparts);
  free(s->zoom_props);
  free(s);

  return 0;
}
