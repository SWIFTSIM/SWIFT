/*******************************************************************************
 * This file is part of SWIFT.
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

/* Local headers. */
#include "swift.h"

/* Some standard headers. */
#include <math.h>
#include <stdio.h>
#include <string.h>

/**
 * Check that a position wrapped back into the periodic box always gets a
 * valid top-level cell index, also when it lies a rounding error away from
 * the lower or the upper face of the box.
 *
 * The positions are given to space_parts_get_cell_index(). A particle at
 * x = -8e-19 is wrapped to dim - 1 ulp. That value, multiplied by the inverse
 * cell width, rounds to exactly cdim and gives a cell index one past the end
 * of the grid.
 */

#define MAX_PARTS 64

/* Number of failed checks. */
static int num_fail = 0;

/* Run the cell index on a set of positions along one axis and check them. */
static void check_axis(struct space *s, const int axis, const double *x_in,
                       const int n) {

  struct part parts[MAX_PARTS];
  int ind[MAX_PARTS];
  int cell_counts[512];
  size_t inhibited = 0, extra = 0;

  bzero(parts, sizeof(parts));
  bzero(cell_counts, sizeof(cell_counts));
  for (int k = 0; k < n; k++) {
    parts[k].x[axis] = x_in[k];
    parts[k].time_bin = 1;
    hydro_set_mass(&parts[k], 1.f);
  }
  s->parts = parts;
  s->nr_parts = n;

  space_parts_get_cell_index(s, ind, cell_counts, &inhibited, &extra, 0);

  for (int k = 0; k < n; k++) {
    const double x = parts[k].x[axis];
    const double dim = s->dim[axis];
    const double scaled = x * s->iwidth[axis];

    /* Periodic distance between the input and the stored position. */
    double d = fabs(x - x_in[k]);
    d = fmin(d, fabs(d - dim));

    const int expected =
        cell_getid(s->cdim, parts[k].x[0] * s->iwidth[0],
                   parts[k].x[1] * s->iwidth[1], parts[k].x[2] * s->iwidth[2]);

    int bad = 0;
    if (ind[k] < 0 || ind[k] >= s->nr_cells) bad = 1;
    if (!(x >= 0. && x < dim)) bad = 1;
    if (!(scaled < s->cdim[axis])) bad = 1;
    if (ind[k] != expected) bad = 1;
    if (!(d < 1e-12 * dim)) bad = 1;

    if (bad) {
      ++num_fail;
      printf(
          "FAIL axis=%d cdim=%d dim=%.17g x_in=%.17g x_out=%.17g "
          "x_out*ih=%.17g index=%d (nr_cells=%d)\n",
          axis, s->cdim[axis], dim, x_in[k], x, scaled, ind[k], s->nr_cells);
    }
  }
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  struct engine e;
  bzero(&e, sizeof(e));
  threadpool_init(&e.threadpool, 1);

  static struct space s;
  bzero(&s, sizeof(s));
  s.e = &e;
  s.periodic = 1;

  /* The first box size is the one of the crash that revealed the problem. */
  const double dims[4] = {0.0050992011205797626, 1., 100., 142.5};

  for (int id = 0; id < 4; ++id) {
    for (int cdim = 3; cdim <= 8; ++cdim) {
      for (int axis = 0; axis < 3; ++axis) {

        const double dim = dims[id];
        for (int i = 0; i < 3; ++i) {
          s.dim[i] = dim;
          s.cdim[i] = cdim;
          s.width[i] = dim / cdim;
          s.iwidth[i] = 1.0 / s.width[i];
        }
        s.nr_cells = cdim * cdim * cdim;

        double x[MAX_PARTS];
        int n = 0;

        /* On the lower face */
        x[n++] = 0.;

        /* Just below the lower face: wraps to dim or dim - 1 ulp */
        for (int k = 1; k <= 8; ++k) x[n++] = -ldexp(dim, -52) * 0.125 * k;
        x[n++] = -8.0547470915348642e-19 * dim / dims[0];

        /* Just below the upper face, already inside the box */
        double xu = dim;
        for (int k = 0; k < 8; ++k) {
          xu = nextafter(xu, 0.);
          x[n++] = xu;
        }

        /* Just above the upper face, wraps to just above 0 */
        x[n++] = dim;
        x[n++] = nextafter(dim, 2. * dim);

        check_axis(&s, axis, x, n);
      }
    }
  }

  threadpool_clean(&e.threadpool);

  if (num_fail) {
    printf("%d checks failed.\n", num_fail);
    return 1;
  }
  printf("All checks passed.\n");
  return 0;
}
