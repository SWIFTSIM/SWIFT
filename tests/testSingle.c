/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

/* Local headers. */
#include "swift.h"

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif

#ifdef DEFAULT_SPH

/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */

int main(int argc, char *argv[]) {

  int k, N = 100;
  struct part p1, p2;
  float x, w, dwdx, r2, dx[3] = {0.0f, 0.0f, 0.0f}, gradw[3];

  /* Greeting message */
  printf("This is %s\n", package_description());

  /* Init the particles. */
  for (k = 0; k < 3; k++) {
    p1.a_hydro[k] = 0.0f;
    p1.v[k] = 0.0f;
    p1.x[k] = 0.0;
    p2.a_hydro[k] = 0.0f;
    p2.v[k] = 0.0f;
    p2.x[k] = 0.0;
  }
  p1.v[0] = 100.0f;
  p1.id = 0;
  p2.id = 1;
  p1.density.wcount = 48.0f;
  p2.density.wcount = 48.0f;
  p1.rho = 1.0f;
  p1.mass = 9.7059e-4;
  p1.h = 0.222871287 / 2;
  p2.rho = 1.0f;
  p2.mass = 9.7059e-4;
  p2.h = 0.222871287 / 2;
  p1.force.soundspeed = 0.0040824829f;
  p1.force.balsara = 0.0f;
  p2.force.soundspeed = 58.8972740361f;
  p2.force.balsara = 0.0f;
  p1.u = 1.e-5 / (hydro_gamma_minus_one * p1.rho);
  p2.u = 1.e-5 / (hydro_gamma_minus_one * p2.rho) + 100.0f / (33 * p2.mass);
  p1.force.P_over_rho2 = p1.u * hydro_gamma_minus_one / p1.rho;
  p2.force.P_over_rho2 = p2.u * hydro_gamma_minus_one / p2.rho;

  /* Dump a header. */
  // printParticle_single(&p1, NULL);
  // printParticle_single(&p2, NULL);
  printf("# r a_1 udt_1 a_2 udt_2\n");

  /* Loop over the different radii. */
  for (k = 1; k <= N; k++) {

    /* Set the distance/radius. */
    dx[0] = -((float)k) / N * fmaxf(p1.h, p2.h) * kernel_gamma;
    r2 = dx[0] * dx[0];

    /* Clear the particle fields. */
    p1.a_hydro[0] = 0.0f;
    p1.force.u_dt = 0.0f;
    p2.a_hydro[0] = 0.0f;
    p2.force.u_dt = 0.0f;

    /* Interact the particles. */
    runner_iact_force(r2, dx, p1.h, p2.h, &p1, &p2);

    /* Clear the particle fields. */
    /* p1.rho = 0.0f; p1.density.wcount = 0.0f;
    p2.rho = 0.0f; p2.density.wcount = 0.0f; */

    /* Interact the particles. */
    // runner_iact_density( r2 , dx , p1.h , p2.h , &p1 , &p2 );

    /* Evaluate just the kernel. */
    x = fabsf(dx[0]) / p1.h;
    kernel_deval(x, &w, &dwdx);
    gradw[0] = dwdx / (p1.h * p1.h * p1.h * p1.h) * dx[0] /
               sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    gradw[1] = dwdx / (p1.h * p1.h * p1.h * p1.h) * dx[1] /
               sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    gradw[2] = dwdx / (p1.h * p1.h * p1.h * p1.h) * dx[2] /
               sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

    /* Output the results. */
    printf(
        "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", -dx[0],
        p1.a_hydro[0], p1.a_hydro[1], p1.a_hydro[2], p1.force.u_dt,
        /// -dx[0] , p1.rho , p1.density.wcount , p2.rho , p2.density.wcount ,
        w, dwdx, gradw[0], gradw[1], gradw[2]);

  } /* loop over radii. */

  /* All is calm, all is bright. */
  return 0;
}
#else

int main(int argc, char *argv[]) { return 0; }

#endif
