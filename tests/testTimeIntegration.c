/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include "swift.h"

#include <stdlib.h>
#include <string.h>

/**
 * @brief Test the kick-drift-kick leapfrog integration
 * via a Sun-Earth simulation
 */
int main(int argc, char *argv[]) {

  struct cell c;
  int i;

  /* Orbital parameters */
  int N_orbits = 8;            /* Number of orbits */
  float G = 6.67384e-11;       /* Newton's constant */
  float M_sun = 1.9885e30;     /* Sun mass [kg] */
  float M_earth = 5.97219e24;  /* Earth mass [kg] */
  float r_max = 152097701000.; /* [m] */
  float r_min = 147098074000.; /* [m] */
  float v_max = 30287.;        /* [m/s] */
  // float v_min = 29291.;           /* [m/s] */

  /* Derived quantities */
  float e = (r_max - r_min) / (r_max + r_min); /* Eccentricity */
  float b = sqrtf(r_max * r_min);              /* Semi-minor axis */
  float p = b * sqrtf(1 - e * e);              /* Semi-lactus rectum */
  float a = p / (1 - e * e);                   /* Semi-major axis */
  float T = sqrtf(4 * M_PI * M_PI * a * a * a /
                  (G * (M_sun + M_earth))); /* Period [s] */

  /* Print some info */
  message("Semi-major axis: a=%e [m]", a);
  message("Semi-minor axis: b=%e [m]", b);
  message("Eccentricity e=%f", e);
  message("Period T=%f [s] = %f days", T, T / (60 * 60 * 24));

  /* Time-step size */
  float dt = 0.001 * T;
  int N = N_orbits * T / dt;

  message("Running for %d steps with dt=%e", N, dt);

  /* Create a particle */
  struct part *parts = NULL;
  parts = (struct part *)malloc(sizeof(struct part));
  bzero(parts, sizeof(struct part));
  struct xpart *xparts = NULL;
  xparts = (struct xpart *)malloc(sizeof(struct xpart));
  bzero(xparts, sizeof(struct xpart));

  /* Put the particle on the orbit */
  parts[0].x[0] = r_max;
  parts[0].x[1] = 0.;
  parts[0].x[2] = 0.;

  parts[0].v[0] = 0.;
  parts[0].v[1] = v_max;
  parts[0].v[2] = 0.;

  xparts[0].v_full[0] = 0.;
  xparts[0].v_full[1] = v_max;
  xparts[0].v_full[2] = 0.;

  /* Set the particle in the cell */
  c.hydro.parts = parts;
  c.hydro.xparts = xparts;
  c.hydro.count = 1;
  c.split = 0;

  /* Create an engine and a fake runner */
  struct runner run;
  struct engine eng;

  run.e = &eng;

  eng.time = 0.;
  eng.time_begin = 0.;
  eng.time_end = N_orbits * T;
  eng.dt_min = dt; /* This forces the time-step to be dt        */
  eng.dt_max = dt; /* irrespective of the state of the particle */

  /* Simulate ! */
  for (i = 0; i < N; i++) {

    /* Move forward in time */
    eng.time_old = eng.time;
    eng.time += dt;

    /* Compute gravitational acceleration */
    float r2 = c.hydro.parts[0].x[0] * c.hydro.parts[0].x[0] +
               c.hydro.parts[0].x[1] * c.hydro.parts[0].x[1];
    float r = sqrtf(r2);
    c.hydro.parts[0].a_hydro[0] =
        -(G * M_sun * c.hydro.parts[0].x[0] / r * r * r);
    c.hydro.parts[0].a_hydro[1] =
        -(G * M_sun * c.hydro.parts[0].x[1] / r * r * r);

    /* Kick... */
    runner_do_kick2(&run, &c, 0);
  }

  /* Clean-up */
  free(parts);
  free(xparts);

  return 0;
}
