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

/* Some standard headers. */
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

/* Includes. */
#include "swift.h"

const double G = 1.;

const size_t N = 16;
const size_t PMGRID = 8;

// const double asmth = 2. * M_PI * const_gravity_a_smooth / boxSize;
// const double asmth2 = asmth * asmth;
// const double fact = G / (M_PI * boxSize) * (1. / (2. * boxSize / PMGRID));

int main() {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Simulation properties */
  const size_t count = N * N * N;
  const double boxSize = 1.;

  /* Create some particles */
  struct gpart* gparts = malloc(count * sizeof(struct gpart));
  bzero(gparts, count * sizeof(struct gpart));
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      for (size_t k = 0; k < N; ++k) {

        struct gpart* gp = &gparts[i * N * N + j * N + k];

        gp->x[0] = i * boxSize / N + boxSize / (2 * N);
        gp->x[1] = j * boxSize / N + boxSize / (2 * N);
        gp->x[2] = k * boxSize / N + boxSize / (2 * N);

        gp->mass = 1. / count;

        gp->id_or_neg_offset = i * N * N + j * N + k;
      }
    }
  }

  /* Properties of the mesh */
  const size_t meshmin[3] = {0, 0, 0};
  const size_t meshmax[3] = {PMGRID - 1, PMGRID - 1, PMGRID - 1};

  const size_t dimx = meshmax[0] - meshmin[0] + 2;
  const size_t dimy = meshmax[1] - meshmin[1] + 2;
  const size_t dimz = meshmax[2] - meshmin[2] + 2;

  const double fac = PMGRID / boxSize;
  const size_t PMGRID2 = 2 * (PMGRID / 2 + 1);

  /* message("dimx=%zd dimy=%zd dimz=%zd", dimx, dimy, dimz); */

  /* Allocate and empty the workspace mesh */
  const size_t workspace_size = (dimx + 4) * (dimy + 4) * (dimz + 4);
  double* workspace = fftw_malloc(workspace_size * sizeof(double));
  bzero(workspace, workspace_size * sizeof(double));

  /* Do CIC with the particles */
  for (size_t pid = 0; pid < count; ++pid) {

    const struct gpart* const gp = &gparts[pid];

    const size_t slab_x =
        (fac * gp->x[0] >= PMGRID) ? PMGRID - 1 : fac * gp->x[0];
    const size_t slab_y =
        (fac * gp->x[1] >= PMGRID) ? PMGRID - 1 : fac * gp->x[1];
    const size_t slab_z =
        (fac * gp->x[2] >= PMGRID) ? PMGRID - 1 : fac * gp->x[2];

    const double dx = fac * gp->x[0] - (double)slab_x;
    const double dy = fac * gp->x[1] - (double)slab_y;
    const double dz = fac * gp->x[2] - (double)slab_z;

    const size_t slab_xx = slab_x + 1;
    const size_t slab_yy = slab_y + 1;
    const size_t slab_zz = slab_z + 1;

    workspace[(slab_x * dimy + slab_y) * dimz + slab_z] +=
        gp->mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
    workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] +=
        gp->mass * (1.0 - dx) * dy * (1.0 - dz);
    workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] +=
        gp->mass * (1.0 - dx) * (1.0 - dy) * dz;
    workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] +=
        gp->mass * (1.0 - dx) * dy * dz;
    workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] +=
        gp->mass * (dx) * (1.0 - dy) * (1.0 - dz);
    workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] +=
        gp->mass * (dx)*dy * (1.0 - dz);
    workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] +=
        gp->mass * (dx) * (1.0 - dy) * dz;
    workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] +=
        gp->mass * (dx)*dy * dz;
  }

  /* for(size_t i = 0 ; i < dimx*dimy*dimz; ++i) */
  /*   message("workspace[%zd] = %f", i, workspace[i]); */

  /* Prepare the force grid */
  const size_t fft_size = workspace_size;
  double* forcegrid = fftw_malloc(fft_size * sizeof(double));
  bzero(forcegrid, fft_size * sizeof(double));

  const size_t sendmin = 0, recvmin = 0;
  const size_t sendmax = PMGRID, recvmax = PMGRID;

  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
         (sendmax - sendmin + 1) * dimy * dimz * sizeof(double));
  
  /* for (size_t i = 0; i < fft_size; ++i) */
  /*   if (forcegrid[i] != workspace[i]) error("wrong"); */

  /* Prepare the density grid */
  double* rhogrid = fftw_malloc(fft_size * sizeof(double));
  bzero(rhogrid, fft_size * sizeof(double));
  
  /* Now get the density */
  for (size_t slab_x = recvmin; slab_x <= recvmax; slab_x++) {

    const size_t slab_xx = slab_x % PMGRID;

    for (size_t slab_y = recvmin; slab_y <= recvmax; slab_y++) {

      const size_t slab_yy = slab_y % PMGRID;

      for (size_t slab_z = recvmin; slab_z <= recvmax; slab_z++) {

        const size_t slab_zz = slab_z % PMGRID;

        rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz] +=
            forcegrid[((slab_x - recvmin) * dimy + (slab_y - recvmin)) * dimz +
                      (slab_z - recvmin)];
      }
    }
  }

  /* for (size_t i = 0; i < 640; i++) { */
  /*   if (rhogrid[i] != workspace[i]) { */
  /*     message("rhogrid[%zd]= %f workspace[%zd]= %f forcegrid[%zd]= %f", i, */
  /*             rhogrid[i], i, workspace[i], i, forcegrid[i]); */
  /*   } */
  /* } */


  /* FFT of the density field */
  fftw_complex* fftgrid = fftw_malloc(fft_size * sizeof(fftw_complex));
  fftw_plan plan_forward = fftw_plan_dft_r2c_3d(PMGRID, PMGRID, PMGRID, rhogrid, fftgrid, FFTW_ESTIMATE);
  fftw_execute(plan_forward);


  for(size_t i = 0; i < 640; i++) {
    message("workspace[%zd]= %f", i, fftgrid[i][0]);
  }
  
  
  /* Clean-up */
  fftw_destroy_plan(plan_forward);
  fftw_free(forcegrid);
  fftw_free(rhogrid);
  fftw_free(workspace);
  free(gparts);
  return 0;
}
