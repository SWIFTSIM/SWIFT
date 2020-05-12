/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

const int num_M2L_runs = 1 << 23;
const int num_M2P_runs = 1 << 23;
const int num_PP_runs = 1;  // << 8;

void make_cell(struct cell *c, int N, const double loc[3], double width,
               int id_base, const struct gravity_props *grav_props) {

  bzero(c, sizeof(struct cell));

  /* Start by setting the basics */
  c->loc[0] = loc[0];
  c->loc[1] = loc[1];
  c->loc[2] = loc[2];
  c->width[0] = width;
  c->width[1] = width;
  c->width[2] = width;

  /* Initialise the locks */
  lock_init(&c->grav.plock);
  lock_init(&c->grav.mlock);

  /* Set the time bins */
  c->grav.ti_end_min = 1;
  c->grav.ti_end_max = 1;
  c->grav.ti_beg_max = 1;
  c->grav.ti_old_part = 1;
  c->grav.ti_old_multipole = 1;

  /* Create the particles */
  c->grav.count = N;
  c->grav.count_total = N;
  c->grav.parts = malloc(N * sizeof(struct gpart));
  bzero(c->grav.parts, N * sizeof(struct gpart));
  for (int i = 0.; i < N; ++i) {

    c->grav.parts[i].id_or_neg_offset = id_base + i;
    c->grav.parts[i].x[0] = loc[0] + width * rand() / ((double)RAND_MAX);
    c->grav.parts[i].x[1] = loc[1] + width * rand() / ((double)RAND_MAX);
    c->grav.parts[i].x[2] = loc[2] + width * rand() / ((double)RAND_MAX);
    c->grav.parts[i].mass = 1.;
    c->grav.parts[i].type = swift_type_dark_matter;
    c->grav.parts[i].time_bin = 1;
  }

  /* Create the multipoles */
  c->grav.multipole = malloc(sizeof(struct gravity_tensors));
  gravity_reset(c->grav.multipole);
  gravity_P2M(c->grav.multipole, c->grav.parts, N, grav_props);
  gravity_multipole_compute_power(&c->grav.multipole->m_pole);
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

  /* Construct gravity properties */
  struct gravity_props grav_props;
  bzero(&grav_props, sizeof(struct gravity_props));
  grav_props.use_advanced_MAC = 1;
  grav_props.use_adaptive_tolerance = 1;
  grav_props.adaptive_tolerance = 1e-4;
  grav_props.theta_crit = 0.5;
  grav_props.G_Newton = 1.;
  grav_props.mesh_size = 64;
  grav_props.a_smooth = 1.25;

  /* Space properites */
  const double dim[3] = {100., 100., 100.};
  const double r_s = grav_props.a_smooth * dim[0] / grav_props.mesh_size;
  const double r_s_inv = 1. / r_s;

  /* Mesh structure */
  struct pm_mesh mesh;
  mesh.periodic = 0;
  mesh.dim[0] = dim[0];
  mesh.dim[1] = dim[1];
  mesh.dim[2] = dim[2];

  /* Construct an engine */
  struct engine e;
  e.mesh = &mesh;
  e.max_active_bin = 56;

  /* Construct a runner */
  struct runner r;
  r.e = &e;

  /* Construct two cells */
  struct cell ci;
  struct cell cj;
  const double loc_i[3] = {0., 0., 0.};
  const double loc_j[3] = {1., 1., 1.};
  const int num_particles = 8;
  make_cell(&ci, num_particles, loc_i, 1., 0, &grav_props);
  make_cell(&cj, num_particles, loc_j, 1., num_particles, &grav_props);

  message("Number of runs: %d", num_M2L_runs);

  /* Construct arrays of multipoles to prevent too much optimization */
  struct gravity_tensors *tensors_i =
      malloc(num_M2L_runs * sizeof(struct gravity_tensors));
  struct gravity_tensors *tensors_j =
      malloc(num_M2L_runs * sizeof(struct gravity_tensors));
  for (int n = 0; n < num_M2L_runs; ++n) {

    memcpy(&tensors_i[n], ci.grav.multipole, sizeof(struct gravity_tensors));
    memcpy(&tensors_j[n], cj.grav.multipole, sizeof(struct gravity_tensors));

    /* Move the values a bit to prevent optimization in the actual loops */
    tensors_i[n].CoM[0] += rand() / ((double)RAND_MAX);
    tensors_i[n].CoM[1] += rand() / ((double)RAND_MAX);
    tensors_i[n].CoM[1] += rand() / ((double)RAND_MAX);

    tensors_j[n].CoM[0] += rand() / ((double)RAND_MAX);
    tensors_j[n].CoM[1] += rand() / ((double)RAND_MAX);
    tensors_j[n].CoM[1] += rand() / ((double)RAND_MAX);

    tensors_i[n].m_pole.M_000 += rand() / ((double)RAND_MAX);
    tensors_j[n].m_pole.M_000 += rand() / ((double)RAND_MAX);

#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    tensors_i[n].m_pole.M_200 += rand() / ((double)RAND_MAX);
    tensors_i[n].m_pole.M_020 += rand() / ((double)RAND_MAX);
    tensors_i[n].m_pole.M_002 += rand() / ((double)RAND_MAX);

    tensors_j[n].m_pole.M_200 += rand() / ((double)RAND_MAX);
    tensors_j[n].m_pole.M_020 += rand() / ((double)RAND_MAX);
    tensors_j[n].m_pole.M_002 += rand() / ((double)RAND_MAX);
#endif
  }

  /* Now run a series of M2L kernels */

  /********
   * Symmetric non-periodic M2L
   ********/
  ticks tic = getticks();
  for (int n = 0; n < num_M2L_runs; ++n) {

    gravity_M2L_symmetric(&tensors_i[n].pot,     //
                          &tensors_j[n].pot,     //
                          &tensors_i[n].m_pole,  //
                          &tensors_j[n].m_pole,  //
                          tensors_i[n].CoM,      //
                          tensors_j[n].CoM,      //
                          &grav_props, /* periodic=*/0, dim, r_s_inv);
  }
  ticks toc = getticks();
  message("%30s at order %d took %4d %s.", "Symmetric non-periodic M2L",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_M2L_runs), "ns");

  /********
   * Symmetric periodic M2L
   ********/
  tic = getticks();
  for (int n = 0; n < num_M2L_runs; ++n) {

    gravity_M2L_symmetric(&tensors_i[n].pot,     //
                          &tensors_j[n].pot,     //
                          &tensors_i[n].m_pole,  //
                          &tensors_j[n].m_pole,  //
                          tensors_i[n].CoM,      //
                          tensors_j[n].CoM,      //
                          &grav_props, /* periodic=*/1, dim, r_s_inv);
  }
  toc = getticks();
  message("%30s at order %d took %4d %s.", "Symmetric periodic M2L",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_M2L_runs), "ns");

  /********
   * Non-symmetric non-periodic M2L
   ********/
  tic = getticks();
  for (int n = 0; n < num_M2L_runs; ++n) {

    gravity_M2L_nonsym(&tensors_i[n].pot,     //
                       &tensors_j[n].m_pole,  //
                       tensors_i[n].CoM,      //
                       tensors_j[n].CoM,      //
                       &grav_props, /* periodic=*/0, dim, r_s_inv);
  }
  toc = getticks();
  message("%30s at order %d took %4d %s.", "Non-symmetric non-periodic M2L",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_M2L_runs), "ns");

  /********
   * Non-symmetric periodic M2L
   ********/
  tic = getticks();
  for (int n = 0; n < num_M2L_runs; ++n) {

    gravity_M2L_nonsym(&tensors_i[n].pot,     //
                       &tensors_j[n].m_pole,  //
                       tensors_i[n].CoM,      //
                       tensors_j[n].CoM,      //
                       &grav_props, /* periodic=*/1, dim, r_s_inv);
  }
  toc = getticks();
  message("%30s at order %d took %4d %s.", "Non-symmetric periodic M2L",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_M2L_runs), "ns");

  /* Now run a series of M2L kernels */

  /********
   * Non-periodic M2P
   ********/
  tic = getticks();
  for (int n = 0; n < num_M2P_runs; ++n) {

    const int index = n % num_particles;

    const float r_x = tensors_j[n].CoM[0] - ci.grav.parts[index].x[0];
    const float r_y = tensors_j[n].CoM[1] - ci.grav.parts[index].x[1];
    const float r_z = tensors_j[n].CoM[2] - ci.grav.parts[index].x[2];
    const float r2 = r_x * r_x + r_y * r_y + r_z * r_z;
    const float eps = gravity_get_softening(&ci.grav.parts[index], &grav_props);

    struct reduced_grav_tensor l = {0.f, 0.f, 0.f, 0.f};
    gravity_M2P(&tensors_j[n].m_pole, r_x, r_y, r_z, r2, eps,
                /*periodic=*/0, r_s_inv, &l);

    ci.grav.parts[index].a_grav[0] += l.F_100;
    ci.grav.parts[index].a_grav[1] += l.F_010;
    ci.grav.parts[index].a_grav[2] += l.F_001;
  }
  toc = getticks();
  message("%30s at order %d took %4d %s.", "Non-periodic M2P",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_M2P_runs), "ns");

  /********
   * Periodic M2P
   ********/
  tic = getticks();
  for (int n = 0; n < num_M2P_runs; ++n) {

    const int index = n % num_particles;

    const float r_x = tensors_j[n].CoM[0] - ci.grav.parts[index].x[0];
    const float r_y = tensors_j[n].CoM[1] - ci.grav.parts[index].x[1];
    const float r_z = tensors_j[n].CoM[2] - ci.grav.parts[index].x[2];
    const float r2 = r_x * r_x + r_y * r_y + r_z * r_z;
    const float eps = gravity_get_softening(&ci.grav.parts[index], &grav_props);

    struct reduced_grav_tensor l = {0.f, 0.f, 0.f, 0.f};
    gravity_M2P(&tensors_j[n].m_pole, r_x, r_y, r_z, r2, eps,
                /*periodic=*/1, r_s_inv, &l);

    ci.grav.parts[index].a_grav[0] += l.F_100;
    ci.grav.parts[index].a_grav[1] += l.F_010;
    ci.grav.parts[index].a_grav[2] += l.F_001;
  }
  toc = getticks();
  message("%30s at order %d took %4d %s.", "Periodic M2P",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_M2P_runs), "ns");

  /* Print out to avoid optimization */
  // gravity_field_tensors_print(&ci.grav.multipole->pot);
  // gravity_field_tensors_print(&cj.grav.multipole->pot);

  tic = getticks();
  for (int n = 0; n < num_PP_runs; ++n) {
    runner_dopair_grav_pp(&r, &ci, &cj, 1, 0);
  }
  toc = getticks();
  message("%30s at order %d took %4d %s.", "dopair_grav (no mpole)",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_PP_runs), "ns");

  tic = getticks();
  runner_dopair_grav_pp(&r, &ci, &cj, 1, 1);
  toc = getticks();
  message("%30s at order %d took %4d %s.", "dopair_grav (mpole)",
          SELF_GRAVITY_MULTIPOLE_ORDER,
          (int)(1e6 * clocks_from_ticks(toc - tic) / num_PP_runs), "ns");

  return 0;
}
