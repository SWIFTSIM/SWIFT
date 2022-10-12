/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#include <config.h>

/* Some standard headers. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

const int L = 8;

struct cell *make_cell(const size_t n, const double offset[3],
                       const double size, long long *partId) {

  const size_t count = n * n * n;
  struct cell *cell = NULL;
  if (posix_memalign((void **)&cell, cell_align, sizeof(struct cell)) != 0) {
    error("Couldn't allocate the cell");
  }
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->grav.parts, gpart_align,
                     count * sizeof(struct gpart)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->grav.parts, count * sizeof(struct gpart));

  if (posix_memalign((void **)&cell->grav.multipole, multipole_align,
                     sizeof(struct gravity_tensors)) != 0) {
    error("Couldn't allocate the mulipole");
  }

  /* Construct the parts */
  struct gpart *gpart = cell->grav.parts;
  for (size_t x = 0; x < n; ++x) {
    for (size_t y = 0; y < n; ++y) {
      for (size_t z = 0; z < n; ++z) {
        gpart->x[0] = offset[0] + random_uniform(0., size);
        gpart->x[1] = offset[1] + random_uniform(0., size);
        gpart->x[2] = offset[2] + random_uniform(0., size);

        gpart->v_full[0] = 0.f;
        gpart->v_full[1] = 0.f;
        gpart->v_full[2] = 0.f;

        gpart->a_grav[0] = 0.f;
        gpart->a_grav[1] = 0.f;
        gpart->a_grav[2] = 0.f;

        gpart->id_or_neg_offset = -(++(*partId));
        gpart->mass = 1.f;
        gpart->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        gpart->ti_drift = 8;
        gpart->ti_kick = 8;
#endif
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->grav.count = count;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->grav.super = cell;
  cell->grav.ti_old_part = 8;
  cell->grav.ti_end_min = 8;

  return cell;
}

void clean_up(struct cell *c) {
  free(c->grav.parts);
  free(c->grav.multipole);
  free(c);
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Initialise a few things to get us going */
  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.max_active_bin = num_time_bins;
  e.time = 0.1f;
  e.ti_current = 8;
  e.time_base = 1e-10;

  struct pm_mesh mesh;
  mesh.periodic = 0;
  mesh.dim[0] = 10.;
  mesh.dim[1] = 10.;
  mesh.dim[2] = 10.;
  mesh.r_s_inv = 0.;
  mesh.r_cut_min = 0.;
  e.mesh = &mesh;

  struct gravity_props props;
  props.a_smooth = 1.25;
  props.epsilon_DM_cur = 1e-10;
  props.epsilon_baryon_cur = 1e-10;
  e.gravity_properties = &props;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &e;

  /* Init the cache for gravity interaction */
  gravity_cache_init(&runner.ci_gravity_cache, L * L * L);
  gravity_cache_init(&runner.cj_gravity_cache, L * L * L);

  long long partID = 0LL;

  const double offset_i[3] = {0., 0., 0.};
  const double offset_j[3] = {2., 0., 0.};

  struct cell *ci = make_cell(L, offset_i, 1.0, &partID);
  struct cell *cj = make_cell(L, offset_j, 1.0, &partID);

  gravity_P2M(ci->grav.multipole, ci->grav.parts, ci->grav.count, &props);
  gravity_P2M(cj->grav.multipole, cj->grav.parts, cj->grav.count, &props);

  runner_dopair_grav_pp(&runner, ci, cj, /*symmetric=*/0, /*allow_mpole=*/1);

  
  
  /* Clean things to make the sanitizer happy ... */
  clean_up(ci);
  clean_up(cj);

  return 0;
}
