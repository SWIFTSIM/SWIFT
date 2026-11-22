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
#include <stdio.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/**
 * Check the debug checks of the position of gas, stars and black holes in
 * their leaf cell.
 *
 * The cell edge is the sum (0.04 + 0.01) + 0.01 = 0.060000000000000005. A
 * particle placed at 0.06 lies 6.9e-18 outside the cell although it sits on
 * its face, and the check must accept it. A particle 1e-3 of a cell width
 * outside must still stop the run. The black hole case needs a black hole
 * model.
 */

enum { GAS, STAR, BH };
static const char *name[3] = {"gas", "star", "black hole"};

static int num_fail = 0;

/* Run the check in a child, since error() ends the process. Return 1 if the
 * check stopped the run. */
static int check_stops_run(const int type, const double offset_in_widths) {

#ifdef SWIFT_DEBUG_CHECKS
  static struct cell c;
  static struct part p;
  static struct spart sp;
  static struct bpart bp;

  /* Not constant-folded: the rounding must happen at run time. */
  volatile double a = 0.04, b = 0.01;
  const double width = b;
  const double edge = (a + b) + b;

  bzero(&c, sizeof(c));
  c.split = 0;
  c.nodeID = 0;
  c.width[0] = c.width[1] = c.width[2] = width;
  for (int k = 0; k < 3; k++) c.loc[k] = edge;

  double x[3];
  for (int k = 0; k < 3; k++) x[k] = 0.06 - offset_in_widths * width;
  if (offset_in_widths == 0. && !(c.loc[0] - x[0] > 0.))
    error("The edge sum does not round outside the face: test is void.");

  bzero(&p, sizeof(p));
  bzero(&sp, sizeof(sp));
  bzero(&bp, sizeof(bp));
  p.time_bin = sp.time_bin = bp.time_bin = 1;
  p.id = sp.id = bp.id = 1;
  memcpy(p.x, x, sizeof(x));
  memcpy(sp.x, x, sizeof(x));
  memcpy(bp.x, x, sizeof(x));
  c.hydro.count = c.stars.count = c.black_holes.count = 1;
  c.hydro.parts = &p;
  c.stars.parts = &sp;
  c.black_holes.parts = &bp;

  fflush(stdout);
  const pid_t pid = fork();
  if (pid == 0) {
    alarm(60);
    if (!freopen("/dev/null", "w", stderr)) _exit(2);
    if (type == GAS) cell_check_hydro_dx_outside_max(&c);
    if (type == STAR) cell_check_stars_dx_outside_max(&c);
    if (type == BH) cell_check_black_holes_dx_outside_max(&c);
    _exit(0);
  }
  int status = 0;
  waitpid(pid, &status, 0);
  return !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
#else
  return 0;
#endif
}

int main(int argc, char *argv[]) {

#ifndef SWIFT_DEBUG_CHECKS
  printf("Skipped: the checks exist only with --enable-debugging-checks.\n");
  return 0;
#endif

  for (int type = GAS; type <= BH; type++) {

#ifdef BLACK_HOLES_NONE
    /* All black hole members of a cell share their bytes in this model. */
    if (type == BH) continue;
#endif

    const int stops_on_face = check_stops_run(type, 0.);
    printf("%-10s on a face, edge rounded out by 1 ulp: %s\n", name[type],
           stops_on_face ? "STOPS (wrong)" : "accepted");
    if (stops_on_face) ++num_fail;

    const int stops_outside = check_stops_run(type, 1e-3);
    printf("%-10s 1e-3 cell widths outside:             %s\n", name[type],
           stops_outside ? "stops" : "ACCEPTED (wrong)");
    if (!stops_outside) ++num_fail;
  }

  return num_fail != 0;
}
