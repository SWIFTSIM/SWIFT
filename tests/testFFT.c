/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Some standard headers. */
#include "../config.h"

// MATTHIEU fix this test
#if 1 || !defined(HAVE_FFTW)

int main(int argc, char *argv[]) { return 0; }

#else

/* Some standard headers. */
#include <fenv.h>
#include <stdlib.h>
#include <string.h>

/* Includes. */
#include "runner_doiact_fft.h"
#include "swift.h"

__attribute__((always_inline)) INLINE static int row_major_id(int i, int j,
                                                              int k, int N) {
  return (((i + N) % N) * N * N + ((j + N) % N) * N + ((k + N) % N));
}

int is_close(double x, double y, double abs_err) {
  return (abs(x - y) < abs_err);
}

int main(int argc, char *argv[]) {
  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Make one particle */
  int nr_gparts = 1;
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, gpart_align,
                     nr_gparts * sizeof(struct gpart)) != 0)
    error("Impossible to allocate memory for gparts.");
  bzero(gparts, nr_gparts * sizeof(struct gpart));

  gparts[0].x[0] = 0.3;
  gparts[0].x[1] = 0.8;
  gparts[0].x[2] = 0.2;
  gparts[0].mass = 1.f;

  /* Read the parameter file */
  struct swift_params *params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  parser_read_file("fft_params.yml", params);

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  /* Initialise the gravity properties */
  struct gravity_props gravity_properties;
  gravity_props_init(&gravity_properties, params, &cosmo);

  /* Build the infrastructure */
  struct space space;
  double dim[3] = {1., 1., 1.};
  space_init(&space, params, &cosmo, dim, NULL, gparts, NULL, 0, nr_gparts, 0,
             1, 1, 0, 0, 1, 1, 0);

  struct engine engine;
  engine.s = &space;
  space.e = &engine;
  engine.time = 0.1f;
  engine.ti_current = 0;
  engine.ti_old = 0;
  engine.max_active_bin = num_time_bins;
  engine.gravity_properties = &gravity_properties;
  engine.nr_threads = 1;
  engine.nodeID = 0;
  engine_rank = 0;
  engine.verbose = 1;

  struct runner runner;
  runner.e = &engine;

  /* Initialize the threadpool. */
  threadpool_init(&engine.threadpool, engine.nr_threads);

  /* Construct the space and all the multipoles. */
  space_rebuild(&space, 1);

/* Initialise the Ewald correction table */
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  gravity_exact_force_ewald_init(dim[0]);
#endif

  /* Run the FFT task */
  runner_do_grav_fft(&runner, 1);

  /* Now check that we got the right answer */
  int nr_cells = space.nr_cells;
  double *r = (double *)malloc(nr_cells * sizeof(double));
  double *m = (double *)malloc(nr_cells * sizeof(double));
  double *pot = (double *)malloc(nr_cells * sizeof(double));
  double *pot_exact = (double *)malloc(nr_cells * sizeof(double));

  FILE *file = fopen("potential.dat", "w");
  for (int i = 0; i < nr_cells; ++i) {
    pot[i] = space.multipoles_top[i].pot.F_000;
    m[i] = space.multipoles_top[i].m_pole.M_000;
    double dx =
        nearest(space.multipoles_top[i].CoM[0] - gparts[0].x[0], dim[0]);
    double dy =
        nearest(space.multipoles_top[i].CoM[1] - gparts[0].x[1], dim[1]);
    double dz =
        nearest(space.multipoles_top[i].CoM[2] - gparts[0].x[2], dim[2]);

    /* Distance */
    r[i] = sqrt(dx * dx + dy * dy + dz * dz);

    /* Potential with correction */
    if (r[i] > 0) pot_exact[i] = 1. / r[i];

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
    /* Get Ewald periodic correction */
    double f_corr[3], pot_corr;
    gravity_exact_force_ewald_evaluate(dx, dy, dz, f_corr, &pot_corr);
    pot_exact[i] -= pot_corr;
#endif

    fprintf(file, "%e %e %e %e\n", r[i], m[i], pot[i], pot_exact[i]);
  }
  fclose(file);

  /* Let's now check the interpolation functions */
  int cdim[3] = {space.cdim[0], space.cdim[1], space.cdim[2]};

  /* Constant function --> Derivatives must be 0 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] = 1.;
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;

        if (!is_close(f->F_000, -1., 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, 0., 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Linear function in x --> Derivatives must be 1 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] = i / ((double)cdim[0]);
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 2; i < cdim[0] - 3; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;

        if (!is_close(f->F_000, -i / ((double)cdim[0]), 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, -1., 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Linear function in y --> Derivatives must be 1 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] = j / ((double)cdim[0]);
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 2; j < cdim[1] - 3; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;

        if (!is_close(f->F_000, -j / ((double)cdim[0]), 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, -1., 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Linear function in z --> Derivatives must be 1 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] = k / ((double)cdim[0]);
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 2; k < cdim[2] - 3; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;

        if (!is_close(f->F_000, -k / ((double)cdim[0]), 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, -1., 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Quadratic function in x --> Derivatives must be 2 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] =
            i * i / ((double)cdim[0] * cdim[0]);
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 2; i < cdim[0] - 3; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;
        const double val = i / ((double)cdim[0]);
        const double val2 = val * val;

        if (!is_close(f->F_000, -val2, 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, -2 * val, 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 2.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Quadratic function in y --> Derivatives must be 2 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] =
            j * j / ((double)cdim[0] * cdim[0]);
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 2; j < cdim[1] - 3; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;
        const double val = j / ((double)cdim[0]);
        const double val2 = val * val;

        if (!is_close(f->F_000, -val2, 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, -2 * val, 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 2.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Quadratic function in z --> Derivatives must be 2 */
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 0; k < cdim[2]; ++k) {
        pot[row_major_id(i, j, k, cdim[0])] =
            k * k / ((double)cdim[0] * cdim[0]);
      }
    }
  }
  for (int i = 0; i < nr_cells; ++i)
    gravity_field_tensors_init(&space.multipoles_top[i].pot, engine.ti_current);
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&space.multipoles_top[i], pot, cdim[0], cdim[0], dim);
  for (int i = 0; i < cdim[0]; ++i) {
    for (int j = 0; j < cdim[1]; ++j) {
      for (int k = 2; k < cdim[2] - 3; ++k) {
        const struct grav_tensor *f =
            &space.multipoles_top[row_major_id(i, j, k, cdim[0])].pot;
        const double val = k / ((double)cdim[0]);
        const double val2 = val * val;

        if (!is_close(f->F_000, -val2, 1e-10))
          error("Invalid value for (%d %d %d) F_000 (%e)", i, j, k, f->F_000);
        if (!is_close(f->F_100, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_100 (%e)", i, j, k, f->F_100);
        if (!is_close(f->F_010, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_010 (%e)", i, j, k, f->F_010);
        if (!is_close(f->F_001, -2 * val, 1e-10))
          error("Invalid value for (%d %d %d) F_001 (%e)", i, j, k, f->F_001);
        if (!is_close(f->F_200, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_200 (%e)", i, j, k, f->F_200);
        if (!is_close(f->F_020, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_020 (%e)", i, j, k, f->F_020);
        if (!is_close(f->F_002, 2.f, 1e-10))
          error("Invalid value for (%d %d %d) F_002 (%e)", i, j, k, f->F_002);
        if (!is_close(f->F_011, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_011 (%e)", i, j, k, f->F_011);
        if (!is_close(f->F_101, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_101 (%e)", i, j, k, f->F_101);
        if (!is_close(f->F_110, 0.f, 1e-10))
          error("Invalid value for (%d %d %d) F_110 (%e)", i, j, k, f->F_110);
      }
    }
  }

  /* Clean up */
  free(r);
  free(m);
  free(pot);
  free(pot_exact);
  free(params);
  free(gparts);
  return 0;
}

#endif
