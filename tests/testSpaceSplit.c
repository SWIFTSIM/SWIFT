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
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "engine.h"
#include "parser.h"
#include "space.h"
#include "swift.h"
#include "zoom_region/zoom_init.h"

double generate_gaussian_coordinate(const double mean, const double std,
                                    const double cell_width) {
  double u1 = (double)rand() / RAND_MAX;
  double u2 = (double)rand() / RAND_MAX;
  double z0 = (sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2)) * std + mean;

  /* We only want to go out at most by the size of a cell. If we've got a
   * coordinate out too far we should try again. */
  if (z0 < mean - cell_width / 2 || z0 > mean + cell_width / 2) {
    return generate_gaussian_coordinate(mean, std, cell_width);
  }

  return z0;
}

void make_mock_space(struct space *s, struct engine *e, const double std) {

  /* Define the members we need for the test. */
  s->with_self_gravity = 1;
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;
  s->cdim[0] = 16;
  s->cdim[1] = 16;
  s->cdim[2] = 16;
  s->nr_gparts = 100000;
  s->width[0] = s->dim[0] / s->cdim[0];
  s->width[1] = s->dim[1] / s->cdim[1];
  s->width[2] = s->dim[2] / s->cdim[2];

  /* Attach the engine (all members have been zeroed which should be sufficient
   * for the test). */
  e->s = s;
  s->e = e;

  /* Allocate memory for the gparts. */
  struct gpart *gparts =
      (struct gpart *)malloc(s->nr_gparts * sizeof(struct gpart));
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* Create gparts randomly sampled from a normal distribution centred on
   * the middle of the box with a width. */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    gparts[i].x[0] =
        generate_gaussian_coordinate(s->dim[0] / 2, std, s->width[0]);
    gparts[i].x[1] =
        generate_gaussian_coordinate(s->dim[1] / 2, std, s->width[1]);
    gparts[i].x[2] =
        generate_gaussian_coordinate(s->dim[2] / 2, std, s->width[2]);
    gparts[i].mass = 1.0;
    gparts[i].type = swift_type_dark_matter;
  }

  s->gparts = gparts;

  /* Allocate sub cells and multipoles. */
  s->cells_sub = (struct cell **)calloc(2, sizeof(struct cell *));
  s->multipoles_sub =
      (struct gravity_tensors **)calloc(2, sizeof(struct gravity_tensors *));
}

void make_cell(struct cell *c, struct space *s) {

  /* Allocate a multipole for this cell. */
  struct gravity_tensors *m = malloc(sizeof(struct gravity_tensors));
  bzero(m, sizeof(struct gravity_tensors));

  c->loc[0] = s->dim[0] / 2 - s->width[0] / 2;
  c->loc[1] = s->dim[1] / 2 - s->width[1] / 2;
  c->loc[2] = s->dim[2] / 2 - s->width[2] / 2;
  c->width[0] = s->width[0];
  c->width[1] = s->width[1];
  c->width[2] = s->width[2];
  c->dmin = s->width[0];
  if (s->with_self_gravity) c->grav.multipole = m;
  c->type = cell_type_regular;
  c->subtype = cell_subtype_regular;
  c->depth = 0;
  c->split = 0;
  c->hydro.count = 0;
  c->grav.count = s->nr_gparts;
  c->stars.count = 0;
  c->sinks.count = 0;
  c->top = c;
  c->super = c;
  c->hydro.super = c;
  c->grav.super = c;
  c->hydro.ti_old_part = s->e->ti_current;
  c->grav.ti_old_part = s->e->ti_current;
  c->stars.ti_old_part = s->e->ti_current;
  c->sinks.ti_old_part = s->e->ti_current;
  c->black_holes.ti_old_part = s->e->ti_current;
  c->grav.ti_old_multipole = s->e->ti_current;
#ifdef WITH_MPI
  c->mpi.tag = -1;
  c->mpi.recv = NULL;
  c->mpi.send = NULL;
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
  c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  cell_assign_top_level_cell_index(c, s);
#endif

  /* Attach the gparts. */
  c->grav.parts = s->gparts;
}

void test_cell_tree(struct cell *c, struct space *s) {

  /* Recurse if the the cell is split. */
  if (c->split) {
    for (int i = 0; i < 8; i++) {
      if (c->progeny[i] == NULL) {
        continue;
      }
      test_cell_tree(c->progeny[i], s);
    }

    /* Ensure the particle counts agree. */
    int count = 0;
    for (int i = 0; i < 8; i++) {
      if (c->progeny[i] == NULL) {
        continue;
      }
      count += c->progeny[i]->grav.count;
    }
    assert(count == c->grav.count);

/* Ensure the multipole counts agree (the member to test this only exists when
 * configured with debugging checks). */
#ifdef SWIFT_DEBUG_CHECKS
    int mpole_count = 0;
    for (int i = 0; i < 8; i++) {
      if (c->progeny[i] == NULL) {
        continue;
      }
      mpole_count += c->progeny[i]->grav.multipole->m_pole.num_gpart;
    }
    assert(mpole_count == c->grav.multipole->m_pole.num_gpart);
    assert(mpole_count == c->grav.count);
    assert(count == mpole_count);
#endif
  }

  /* Nothing to check if we're in a leaf */
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* Define the standard deviation of the particle distribution with some
   * randomness to test different tree depths (1-33 for default RAND_MAX). The
   * smaller stds will result in a more concentrated particle distribution. */
  const double std = 49 * ((double)rand() / RAND_MAX) + 1;

  /* Create a space and engine structure. */
  struct space s;
  struct engine e;
  bzero(&s, sizeof(struct space));
  bzero(&e, sizeof(struct engine));
  make_mock_space(&s, &e, std);

  /* Create a cell to hold all the particles. */
  struct cell *c = NULL;  // malloc(sizeof(struct cell));
  if (posix_memalign((void **)&c, cell_align, sizeof(struct cell)) != 0) {
    error("Couldn't allocate the cell");
  }
  // message("%p",c);
  bzero(c, sizeof(struct cell));
  make_cell(c, &s);

  /* Allocate the gpart buffer and set the others to NULL. */
  struct cell_buff *buff = NULL, *sbuff = NULL, *bbuff = NULL, *gbuff = NULL,
                   *sink_buff = NULL;
  space_allocate_and_fill_buffers(c, &buff, &sbuff, &bbuff, &gbuff, &sink_buff);

  /* Recursively split the cell. */
  space_split_recursive(&s, c, buff, sbuff, bbuff, gbuff, sink_buff,
                        /*tpid*/ 0);

  /* Free the particle buffers. */
  swift_free("tempgbuff", gbuff);

  /* Test the cell tree. */
  test_cell_tree(c, &s);

  /* Clean up. */
  free(s.gparts);
  free(s.cells_sub);
  free(s.multipoles_sub);
  free(c->grav.multipole);
  free(c);
}
