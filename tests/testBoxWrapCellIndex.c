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
#include <sys/wait.h>
#include <unistd.h>

/**
 * Check the top-level cell index of positions that lie a rounding error away
 * from a face of the box, for the five mappers (parts, gparts, sparts, bparts,
 * sinks), in a periodic box and in a box with walls.
 *
 * A particle at x = -8e-19 is wrapped to dim - 1 ulp. That value, multiplied
 * by the inverse cell width, rounds to exactly cdim and gives a cell index one
 * past the end of the grid. In a periodic box the particle goes to x = 0. With
 * walls it moves down by one ulp, onto the last cell. A particle more than one
 * box length outside stays outside, so that the debug checks still stop the
 * run.
 */

#define MAX_PARTS 64
#define NUM_TYPES 5

static const char *type_name[NUM_TYPES] = {"part", "gpart", "spart", "bpart",
                                           "sink"};

/* Number of failed checks. */
static int num_fail = 0;

static struct part P[MAX_PARTS];
static struct gpart G[MAX_PARTS];
static struct spart S[MAX_PARTS];
static struct bpart B[MAX_PARTS];
static struct sink K[MAX_PARTS];

/* Position array of particle k of a given type. */
static double *pos(const int type, const int k) {
  switch (type) {
    case 0:
      return P[k].x;
    case 1:
      return G[k].x;
    case 2:
      return S[k].x;
    case 3:
      return B[k].x;
    default:
      return K[k].x;
  }
}

/* Reset n particles of a type and put x_in on one axis. */
static void fill(const int type, const int axis, const double *x_in,
                 const int n) {
  for (int k = 0; k < n; k++) {
    switch (type) {
      case 0:
        bzero(&P[k], sizeof(P[k]));
        P[k].time_bin = 1;
        hydro_set_mass(&P[k], 1.f);
        break;
      case 1:
        bzero(&G[k], sizeof(G[k]));
        G[k].time_bin = 1;
        G[k].mass = 1.f;
        G[k].type = swift_type_dark_matter;
        break;
      case 2:
        bzero(&S[k], sizeof(S[k]));
        S[k].time_bin = 1;
        S[k].mass = 1.f;
        break;
      case 3:
        bzero(&B[k], sizeof(B[k]));
        B[k].time_bin = 1;
        B[k].mass = 1.f;
        break;
      default:
        bzero(&K[k], sizeof(K[k]));
        K[k].time_bin = 1;
        K[k].mass = 1.f;
        break;
    }
    pos(type, k)[axis] = x_in[k];
  }
}

/* Call the cell index mapper of a type on the n filled particles. */
static void get_index(struct space *s, const int type, int *ind, const int n) {
  int cell_counts[512];
  size_t inhibited = 0, extra = 0;
  bzero(cell_counts, sizeof(cell_counts));

  switch (type) {
    case 0:
      s->parts = P;
      s->nr_parts = n;
      space_parts_get_cell_index(s, ind, cell_counts, &inhibited, &extra, 0);
      break;
    case 1:
      s->gparts = G;
      s->nr_gparts = n;
      space_gparts_get_cell_index(s, ind, cell_counts, &inhibited, &extra, 0);
      break;
    case 2:
      s->sparts = S;
      s->nr_sparts = n;
      space_sparts_get_cell_index(s, ind, cell_counts, &inhibited, &extra, 0);
      break;
    case 3:
      s->bparts = B;
      s->nr_bparts = n;
      space_bparts_get_cell_index(s, ind, cell_counts, &inhibited, &extra, 0);
      break;
    default:
      s->sinks = K;
      s->nr_sinks = n;
      space_sinks_get_cell_index(s, ind, cell_counts, &inhibited, &extra, 0);
      break;
  }
}

/* Set the box and the top-level grid like space_regrid() does. */
static void set_grid(struct space *s, const double dim, const int cdim,
                     const int periodic) {
  for (int i = 0; i < 3; ++i) {
    s->dim[i] = dim;
    s->cdim[i] = cdim;
    s->width[i] = dim / cdim;
    s->iwidth[i] = 1.0 / s->width[i];
  }
  s->nr_cells = cdim * cdim * cdim;
  s->periodic = periodic;
}

/* Run the cell index on positions along one axis and check them. */
static void check_axis(struct space *s, const int type, const int axis,
                       const double *x_in, const int n) {

  int ind[MAX_PARTS];
  fill(type, axis, x_in, n);
  get_index(s, type, ind, n);

  for (int k = 0; k < n; k++) {
    const double *p = pos(type, k);
    const double x = p[axis];
    const double dim = s->dim[axis];
    const double scaled = x * s->iwidth[axis];

    /* Distance between the input and the stored position: periodic distance
     * in a periodic box, plain distance with walls. */
    double d = fabs(x - x_in[k]);
    if (s->periodic) d = fmin(d, fabs(d - dim));

    const int expected = cell_getid(s->cdim, p[0] * s->iwidth[0],
                                    p[1] * s->iwidth[1], p[2] * s->iwidth[2]);

    int bad = 0;
    if (ind[k] < 0 || ind[k] >= s->nr_cells) bad = 1;
    if (!(x >= 0. && x < dim)) bad = 1;
    if (!(scaled < s->cdim[axis])) bad = 1;
    if (ind[k] != expected) bad = 1;
    if (!(d < 1e-12 * dim)) bad = 1;

    if (bad) {
      ++num_fail;
      printf(
          "FAIL %s periodic=%d axis=%d cdim=%d dim=%.17g x_in=%.17g "
          "x_out=%.17g x_out*ih=%.17g index=%d (nr_cells=%d)\n",
          type_name[type], s->periodic, axis, s->cdim[axis], dim, x_in[k], x,
          scaled, ind[k], s->nr_cells);
    }
  }
}

/* A particle a rounding error below the upper wall must stay in the last cell
 * and move by one ulp, not jump across the box. */
static void check_upper_wall(struct space *s, const int type) {

  const double dim = 0.0050992011205797626;
  set_grid(s, dim, 3, /*periodic=*/0);

  const double x_in = nextafter(dim, 0.);
  int ind[1];
  fill(type, 0, &x_in, 1);
  get_index(s, type, ind, 1);
  const double x_out = pos(type, 0)[0];

  /* The case only means something if the product rounds to cdim. */
  const int precondition = (x_in * s->iwidth[0] >= 3.);
  const int bad = !precondition || x_out != nextafter(x_in, 0.) ||
                  ind[0] != 2 * s->cdim[1] * s->cdim[2];
  if (bad) {
    ++num_fail;
    printf(
        "FAIL %s upper wall: x_in*ih=%.17g x_in=%.17g x_out=%.17g index=%d "
        "(expected %d)\n",
        type_name[type], x_in * s->iwidth[0], x_in, x_out, ind[0],
        2 * s->cdim[1] * s->cdim[2]);
  }
}

/* A particle more than one box length outside is not wrapped by box_wrap()
 * into the box. The guard must leave it there. */
static void check_escaped(struct space *s, const int type) {

  set_grid(s, 1., 4, /*periodic=*/1);

  /* Along z, so that the (invalid) index stays inside the cell buffer. */
  const double x_in = 2.5;
  int ind[1];

#ifdef SWIFT_DEBUG_CHECKS
  /* The debug check must stop the run. */
  int fd[2];
  if (pipe(fd) != 0) error("pipe failed");
  fflush(stdout);
  const pid_t pid = fork();
  if (pid == 0) {
    dup2(fd[1], 2);
    close(fd[0]);
    close(fd[1]);
    alarm(60);
    fill(type, 2, &x_in, 1);
    get_index(s, type, ind, 1);
    _exit(0);
  }
  close(fd[1]);
  char msg[8192];
  ssize_t len = 0, r;
  while (len < (ssize_t)sizeof(msg) - 1 &&
         (r = read(fd[0], msg + len, sizeof(msg) - 1 - len)) > 0)
    len += r;
  msg[len] = '\0';
  close(fd[0]);
  int status = 0;
  waitpid(pid, &status, 0);
  const int died = !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
  const int said = strstr(msg, "Particle outside of simulation box") != NULL;
  if (!died || !said) {
    ++num_fail;
    printf(
        "FAIL %s escaped: debug check did not stop the run (died=%d said=%d)\n",
        type_name[type], died, said);
  }
#else
  /* Without the debug checks, the position must be left outside the box. */
  fill(type, 2, &x_in, 1);
  get_index(s, type, ind, 1);
  const double x_out = pos(type, 0)[2];
  if (x_out != 1.5) {
    ++num_fail;
    printf("FAIL %s escaped: x_in=%.17g x_out=%.17g (expected 1.5)\n",
           type_name[type], x_in, x_out);
  }
#endif
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

  /* The first box size is the one of the crash that revealed the problem. */
  const double dims[4] = {0.0050992011205797626, 1., 100., 142.5};

  for (int periodic = 0; periodic <= 1; ++periodic) {
    for (int type = 0; type < NUM_TYPES; ++type) {
      for (int id = 0; id < 4; ++id) {
        for (int cdim = 3; cdim <= 8; ++cdim) {
          for (int axis = 0; axis < 3; ++axis) {

            const double dim = dims[id];
            set_grid(&s, dim, cdim, periodic);

            double x[MAX_PARTS];
            int n = 0;

            /* On the lower face */
            x[n++] = 0.;

            /* Just below the upper face, already inside the box */
            double xu = dim;
            for (int k = 0; k < 8; ++k) {
              xu = nextafter(xu, 0.);
              x[n++] = xu;
            }

            if (periodic) {
              /* Just below the lower face: wraps to dim or dim - 1 ulp */
              for (int k = 1; k <= 8; ++k)
                x[n++] = -ldexp(dim, -52) * 0.125 * k;
              x[n++] = -8.0547470915348642e-19 * dim / dims[0];

              /* On and just above the upper face: wraps to just above 0 */
              x[n++] = dim;
              x[n++] = nextafter(dim, 2. * dim);
            }

            check_axis(&s, type, axis, x, n);
          }
        }
      }
      if (!periodic) check_upper_wall(&s, type);
      if (periodic) check_escaped(&s, type);
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
