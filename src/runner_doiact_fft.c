/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* This object's header. */
#include "runner_doiact_fft.h"

/* Local includes. */
#include "engine.h"
#include "error.h"
#include "kernel_long_gravity.h"
#include "runner.h"
#include "space.h"
#include "timers.h"

#ifdef HAVE_FFTW

/**
 * @brief Returns 1D index of a 3D NxNxN array using row-major style.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline)) INLINE static int row_major_id(int i, int j,
                                                              int k, int N) {
  return ((i % N) * N * N + (j % N) * N + (k % N));
}

/**
 * @brief Assigns a given multipole to a density mesh using the CIC method.
 *
 * @param m The #multipole.
 * @param rho The density mesh.
 * @param N the size of the mesh along one axis.
 * @param fac The width of a mesh cell.
 */
__attribute__((always_inline)) INLINE static void multipole_to_mesh_CIC(
    const struct gravity_tensors* m, double* rho, int N, double fac) {

  int i = (int)(fac * m->CoM[0]);
  if (i >= N) i = N - 1;
  const double dx = fac * m->CoM[0] - i;
  const double tx = 1. - dx;

  int j = (int)(fac * m->CoM[1]);
  if (j >= N) j = N - 1;
  const double dy = fac * m->CoM[1] - j;
  const double ty = 1. - dy;

  int k = (int)(fac * m->CoM[2]);
  if (k >= N) k = N - 1;
  const double dz = fac * m->CoM[2] - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i >= N) error("Invalid multipole position in x");
  if (j < 0 || j >= N) error("Invalid multipole position in y");
  if (k < 0 || k >= N) error("Invalid multipole position in z");
#endif

  /* CIC ! */
  rho[row_major_id(i + 0, j + 0, k + 0, N)] += m->m_pole.M_000 * tx * ty * tz;
  rho[row_major_id(i + 0, j + 0, k + 1, N)] += m->m_pole.M_000 * tx * ty * dz;
  rho[row_major_id(i + 0, j + 1, k + 0, N)] += m->m_pole.M_000 * tx * dy * tz;
  rho[row_major_id(i + 0, j + 1, k + 1, N)] += m->m_pole.M_000 * tx * dy * dz;
  rho[row_major_id(i + 1, j + 0, k + 0, N)] += m->m_pole.M_000 * dx * ty * tz;
  rho[row_major_id(i + 1, j + 0, k + 1, N)] += m->m_pole.M_000 * dx * ty * dz;
  rho[row_major_id(i + 1, j + 1, k + 0, N)] += m->m_pole.M_000 * dx * dy * tz;
  rho[row_major_id(i + 1, j + 1, k + 1, N)] += m->m_pole.M_000 * dx * dy * dz;
}

/**
 * @brief Computes the potential on a multipole from a given mesh using the CIC
 * method.
 *
 * @param m The #multipole.
 * @param pot The potential mesh.
 * @param N the size of the mesh along one axis.
 * @param fac width of a mesh cell.
 */
__attribute__((always_inline)) INLINE static void mesh_to_multipole_CIC(
    struct gravity_tensors* m, double* pot, int N, double fac) {

  int i = (int)(fac * m->CoM[0]);
  if (i >= N) i = N - 1;
  const double dx = fac * m->CoM[0] - i;
  const double tx = 1. - dx;

  int j = (int)(fac * m->CoM[1]);
  if (j >= N) j = N - 1;
  const double dy = fac * m->CoM[1] - j;
  const double ty = 1. - dy;

  int k = (int)(fac * m->CoM[2]);
  if (k >= N) k = N - 1;
  const double dz = fac * m->CoM[2] - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i >= N) error("Invalid multipole position in x");
  if (j < 0 || j >= N) error("Invalid multipole position in y");
  if (k < 0 || k >= N) error("Invalid multipole position in z");
#endif

  /* CIC ! */
  m->pot.F_000 += pot[row_major_id(i + 0, j + 0, k + 0, N)] * tx * ty * tz;
  m->pot.F_000 += pot[row_major_id(i + 0, j + 0, k + 1, N)] * tx * ty * dz;
  m->pot.F_000 += pot[row_major_id(i + 0, j + 1, k + 0, N)] * tx * dy * tz;
  m->pot.F_000 += pot[row_major_id(i + 0, j + 1, k + 1, N)] * tx * dy * dz;
  m->pot.F_000 += pot[row_major_id(i + 1, j + 0, k + 0, N)] * dx * ty * tz;
  m->pot.F_000 += pot[row_major_id(i + 1, j + 0, k + 1, N)] * dx * ty * dz;
  m->pot.F_000 += pot[row_major_id(i + 1, j + 1, k + 0, N)] * dx * dy * tz;
  m->pot.F_000 += pot[row_major_id(i + 1, j + 1, k + 1, N)] * dx * dy * dz;
}

#endif

/**
 * @brief Computes the potential on the top multipoles using a Fourier transform
 *
 * @param r The #runner task
 * @param timer Are we timing this ?
 */
void runner_do_grav_fft(struct runner* r, int timer) {

#ifdef HAVE_FFTW

  const struct engine* e = r->e;
  const struct space* s = e->s;
  const integertime_t ti_current = e->ti_current;
  const double a_smooth = e->gravity_properties->a_smooth;
  const double box_size = s->dim[0];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};

  TIMER_TIC;

  if (cdim[0] != cdim[1] || cdim[0] != cdim[2]) error("Non-square mesh");

  /* Some useful constants */
  const int N = cdim[0];
  const int N_half = N / 2;
  const double cell_fac = N / box_size;

  /* Recover the list of top-level multipoles */
  const int nr_cells = s->nr_cells;
  struct gravity_tensors* restrict multipoles = s->multipoles_top;
  struct cell* cells = s->cells_top;

  /* Make sure everything has been drifted to the current point */
  for (int i = 0; i < nr_cells; ++i)
    if (cells[i].ti_old_multipole != ti_current)
      cell_drift_multipole(&cells[i], e);
  // error("Top-level multipole %d not drifted", i);

  /* Allocates some memory for the density mesh */
  double* restrict rho = fftw_malloc(sizeof(double) * N * N * N);
  if (rho == NULL) error("Error allocating memory for density mesh");

  /* Allocates some memory for the mesh in Fourier space */
  fftw_complex* restrict frho =
      fftw_malloc(sizeof(fftw_complex) * N * N * (N_half + 1));
  if (frho == NULL)
    error("Error allocating memory for transform of density mesh");

  /* Prepare the FFT library */
  fftw_plan forward_plan = fftw_plan_dft_r2c_3d(
      N, N, N, rho, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  fftw_plan inverse_plan = fftw_plan_dft_c2r_3d(
      N, N, N, frho, rho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  /* Do a CIC mesh assignment of the multipoles */
  bzero(rho, N * N * N * sizeof(double));
  for (int i = 0; i < nr_cells; ++i)
    multipole_to_mesh_CIC(&multipoles[i], rho, N, cell_fac);

  /* Fourier transform to go to magic-land */
  fftw_execute(forward_plan);

  /* frho now contains the Fourier transform of the density field */
  /* frho contains NxNx(N/2+1) complex numbers */

  /* Some common factors */
  const double green_fac = -1. / (M_PI * box_size);
  const double a_smooth2 = 4. * M_PI * a_smooth * a_smooth / ((double)(N * N));
  const double k_fac = M_PI / (double)N;

  /* Now de-convolve the CIC kernel and apply the Green function */
  for (int i = 0; i < N; ++i) {

    /* kx component of vector in Fourier space and 1/sinc(kx) */
    const int kx = (i > N_half ? i - N : i);
    const double kx_d = (double)kx;
    const double fx = k_fac * kx_d;
    const double sinc_kx_inv = (kx != 0) ? fx / sin(fx) : 1.;

    for (int j = 0; j < N; ++j) {

      /* ky component of vector in Fourier space and 1/sinc(ky) */
      const int ky = (j > N_half ? j - N : j);
      const double ky_d = (double)ky;
      const double fy = k_fac * ky_d;
      const double sinc_ky_inv = (ky != 0) ? fy / sin(fy) : 1.;

      for (int k = 0; k < N_half + 1; ++k) {

        /* kz component of vector in Fourier space and 1/sinc(kz) */
        const int kz = (k > N_half ? k - N : k);
        const double kz_d = (double)kz;
        const double fz = k_fac * kz_d;
        const double sinc_kz_inv = (kz != 0) ? fz / sin(fz) : 1.;

        /* Norm of vector in Fourier space */
        const double k2 = (kx_d * kx_d + ky_d * ky_d + kz_d * kz_d);

        /* Avoid FPEs... */
        if (k2 == 0.) continue;

        /* Green function */
        double W;
        fourier_kernel_long_grav_eval(k2 * a_smooth2, &W);
        const double green_cor = green_fac * W / k2;

        /* Deconvolution of CIC */
        const double CIC_cor = sinc_kx_inv * sinc_ky_inv * sinc_kz_inv;
        const double CIC_cor2 = CIC_cor * CIC_cor;
        const double CIC_cor4 = CIC_cor2 * CIC_cor2;

        /* Combined correction */
        const double total_cor = green_cor * CIC_cor4;

        /* Apply to the mesh */
        const int index = N * (N_half + 1) * i + (N_half + 1) * j + k;
        frho[index][0] *= total_cor;
        frho[index][1] *= total_cor;
      }
    }
  }

  /* Correct singularity at (0,0,0) */
  frho[0][0] = 0.;
  frho[0][1] = 0.;

  /* Fourier transform to come back from magic-land */
  fftw_execute(inverse_plan);

  /* rho now contains the potential */
  /* This array is now again NxNxN real numbers */

  /* Get the potential from the mesh using CIC */
  for (int i = 0; i < nr_cells; ++i)
    mesh_to_multipole_CIC(&multipoles[i], rho, N, cell_fac);

  /* Clean-up the mess */
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);
  fftw_free(rho);
  fftw_free(frho);

  /* Time the whole thing */
  if (timer) TIMER_TOC(timer_dograv_top_level);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

#ifdef HAVE_FFTW
void print_array(double* array, int N) {

  for (int k = N - 1; k >= 0; --k) {
    printf("--- z = %d ---------\n", k);
    for (int j = N - 1; j >= 0; --j) {
      for (int i = 0; i < N; ++i) {
        printf("%f ", array[i * N * N + j * N + k]);
      }
      printf("\n");
    }
  }
}

void print_carray(fftw_complex* array, int N) {

  for (int k = N - 1; k >= 0; --k) {
    printf("--- z = %d ---------\n", k);
    for (int j = N - 1; j >= 0; --j) {
      for (int i = 0; i < N; ++i) {
        printf("(%f %f) ", array[i * N * N + j * N + k][0],
               array[i * N * N + j * N + k][1]);
      }
      printf("\n");
    }
  }
}
#endif /* HAVE_FFTW */
