/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Will J. Roper (w.roper@sussex.ac.uk).
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
#include "gravity_properties.h"
#include "parser.h"
#include "space.h"
#include "swift.h"
#include "zoom_region/zoom.h"

struct truncation_result {
  double box_dim;
  double high_res_dim;
  size_t nr_inhibited_background;
  size_t nr_active_background;
};

static void make_mock_space(struct space *s, const double zoom_width) {

  s->dim[0] = 1000.0;
  s->dim[1] = 1000.0;
  s->dim[2] = 1000.0;
  s->periodic = 1;
  s->nr_gparts = 18;

  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate memory for gparts");
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  double cube_corners[8][3] = {
      {500 - (zoom_width / 2), 500 - (zoom_width / 2), 500 - (zoom_width / 2)},
      {500 - (zoom_width / 2), 500 + (zoom_width / 2), 500 - (zoom_width / 2)},
      {500 + (zoom_width / 2), 500 - (zoom_width / 2), 500 - (zoom_width / 2)},
      {500 + (zoom_width / 2), 500 + (zoom_width / 2), 500 - (zoom_width / 2)},
      {500 - (zoom_width / 2), 500 - (zoom_width / 2), 500 + (zoom_width / 2)},
      {500 - (zoom_width / 2), 500 + (zoom_width / 2), 500 + (zoom_width / 2)},
      {500 + (zoom_width / 2), 500 - (zoom_width / 2), 500 + (zoom_width / 2)},
      {500 + (zoom_width / 2), 500 + (zoom_width / 2), 500 + (zoom_width / 2)}};

  for (size_t i = 0; i < s->nr_gparts; i++) {
    gparts[i].mass = 1.0;
    gparts[i].time_bin = 1;

    if (i < 10) {
      gparts[i].x[0] = s->dim[0] / s->nr_gparts * i;
      gparts[i].x[1] = s->dim[1] / s->nr_gparts * i;
      gparts[i].x[2] = s->dim[2] / s->nr_gparts * i;
      gparts[i].type = swift_type_dark_matter_background;
    } else {
      gparts[i].x[0] = cube_corners[i - 10][0];
      gparts[i].x[1] = cube_corners[i - 10][1];
      gparts[i].x[2] = cube_corners[i - 10][2];
      gparts[i].type = swift_type_dark_matter;
    }
  }

  s->gparts = gparts;
}

static void free_mock_space(struct space *s) {
  free(s->zoom_props);
  free(s->gparts);
  free(s);
}

static struct truncation_result run_truncation_case(const char *param_file_name,
                                                    const double zoom_width) {

  struct swift_params params;
  parser_read_file(param_file_name, &params);

  struct space *s = malloc(sizeof(struct space));
  struct engine *e = malloc(sizeof(struct engine));
  struct gravity_props *gravity = malloc(sizeof(struct gravity_props));
  if (s == NULL || e == NULL || gravity == NULL)
    error("Failed to allocate truncation test structures.");

  bzero(s, sizeof(struct space));
  bzero(e, sizeof(struct engine));
  bzero(gravity, sizeof(struct gravity_props));

  make_mock_space(s, zoom_width);
  s->with_zoom_region = 1;
  s->e = e;

  gravity->r_s = 1.0;
  gravity->r_cut_max_ratio = 1.0;
  e->gravity_properties = gravity;
  e->restarting = 0;
  threadpool_init(&e->threadpool, 1);

  zoom_props_init(&params, s, /*verbose=*/0);
  zoom_region_init(s, /*regridding=*/0, /*verbose=*/0);

  struct truncation_result result = {0.0, 0.0, 0, 0};
  result.box_dim = s->dim[0];
  result.high_res_dim =
      max3(s->zoom_props->part_dim[0], s->zoom_props->part_dim[1],
           s->zoom_props->part_dim[2]);

  for (size_t i = 0; i < s->nr_gparts; i++) {
    const struct gpart *gp = &s->gparts[i];
    if (gp->type != swift_type_dark_matter_background) continue;

    if (gp->time_bin == time_bin_inhibited) {
      result.nr_inhibited_background++;
    } else {
      result.nr_active_background++;
      for (int j = 0; j < 3; j++) {
        if (gp->x[j] < 0.0 || gp->x[j] >= s->dim[j]) {
          error("Active background particle lies outside the retained box.");
        }
      }
    }
  }

  for (size_t i = 0; i < s->nr_gparts; i++) {
    const struct gpart *gp = &s->gparts[i];
    if (gp->type == swift_type_dark_matter_background) continue;
    if (gp->time_bin == time_bin_inhibited)
      error(
          "High-resolution particle was incorrectly inhibited by truncation.");
  }

  threadpool_clean(&e->threadpool);
  free(gravity);
  free(e);
  free_mock_space(s);

  return result;
}

int main(int argc, char *argv[]) {

  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  if (argc != 5)
    error("Usage: %s <base.yml> <tight.yml> <strong.yml> <zoom_width>",
          argv[0]);

  const double zoom_width = atof(argv[4]);

  const struct truncation_result base =
      run_truncation_case(argv[1], zoom_width);
  const struct truncation_result tight =
      run_truncation_case(argv[2], zoom_width);
  const struct truncation_result strong =
      run_truncation_case(argv[3], zoom_width);
  const struct truncation_result wide =
      run_truncation_case(argv[1], 2.0 * zoom_width);

  assert(base.high_res_dim == zoom_width);
  assert(wide.high_res_dim == 2.0 * zoom_width);

  assert(tight.box_dim > base.box_dim);
  assert(tight.nr_active_background >= base.nr_active_background);

  assert(strong.box_dim > base.box_dim);
  assert(strong.nr_active_background >= base.nr_active_background);

  assert(wide.box_dim > base.box_dim);
  assert(wide.nr_active_background >= base.nr_active_background);

  return 0;
}
