/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifdef HAVE_FFTW
#include <fftw3.h>
#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
#include <fftw3-mpi.h>
#endif
#endif

/* This object's header. */
#include "mesh_gravity.h"

/* Local includes. */
#include "active.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "gravity_properties.h"
#include "kernel_long_gravity.h"
#include "mesh_gravity.h"
#include "mesh_gravity_mpi.h"
#include "mesh_gravity_patch.h"
#include "neutrino.h"
#include "part.h"
#include "restart.h"
#include "row_major_id.h"
#include "runner.h"
#include "space.h"
#include "threadpool.h"
#include "zoom_region.h"

/* Standard includes */
#include <math.h>

/**
 * @brief Compute the mesh forces and potential, including periodic correction.
 *
 * Interpolates the top-level multipoles on-to a mesh, move to Fourier space,
 * compute the potential including short-range correction and move back
 * to real space. We use CIC for the interpolation.
 *
 * This version stores the full N*N*N mesh on each MPI rank and uses the
 * non-MPI version of FFTW.
 *
 * The particles mesh accelerations and potentials are also updated.
 *
 * @param mesh The #pm_mesh used to store the potential.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
void compute_potential_zoom(struct pm_mesh* mesh, const struct space* s,
                            struct threadpool* tp, const int verbose) {

#ifdef HAVE_FFTW

  const double r_s = mesh->r_s;
  const double region_size = mesh->dim[0];
  const double dim[3] = {mesh->dim[0], mesh->dim[1], mesh->dim[2]};
  const int* local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  if (r_s <= 0.) error("Invalid value of a_smooth");
  if (mesh->dim[0] != dim[0] || mesh->dim[1] != dim[1] ||
      mesh->dim[2] != dim[2])
    error("Domain size does not match the value stored in the space.");

  /* Some useful constants */
  const int N = mesh->N;
  const int N_half = N / 2;
  const double cell_fac = N / region_size;

  /* Use the memory allocated for the potential to temporarily store rho */
  double* restrict rho = mesh->potential_global;
  if (rho == NULL) error("Error allocating memory for density mesh");

  /* Allocates some memory for the mesh in Fourier space */
  fftw_complex* restrict frho =
      (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N_half + 1));
  if (frho == NULL)
    error("Error allocating memory for transform of density mesh");
  memuse_log_allocation("fftw_frho", frho, 1,
                        sizeof(fftw_complex) * N * N * (N_half + 1));

  /* Prepare the FFT library */
  fftw_plan forward_plan = fftw_plan_dft_r2c_3d(
      N, N, N, rho, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  fftw_plan inverse_plan = fftw_plan_dft_c2r_3d(
      N, N, N, frho, rho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  ticks tic = getticks();

  /* Zero everything */
  bzero(rho, N * N * N * sizeof(double));

  /* Gather some neutrino constants if using delta-f weighting on the mesh */
  struct neutrino_model nu_model;
  bzero(&nu_model, sizeof(struct neutrino_model));
  if (s->e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(s, &nu_model);

  /* Gather the mesh shared information to be used by the threads */
  struct cic_mapper_data data;
  data.cells = s->cells_top;
  data.rho = rho;
  data.potential = NULL;
  data.N = N;
  data.use_local_patches = mesh->use_local_patches;
  data.fac = cell_fac;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];
  data.const_G = 0.f;
  data.nu_model = &nu_model;

  if (nr_local_cells == 0) {

    /* We don't have a cell infrastructure in place so we need to
     * directly loop over the particles */
    threadpool_map(tp, gpart_to_mesh_CIC_mapper, s->gparts, s->nr_gparts,
                   sizeof(struct gpart), threadpool_auto_chunk_size,
                   (void*)&data);

  } else { /* Normal case */

    /* Do a parallel CIC mesh assignment of the gparts but only using
     * the local top-level cells */
    threadpool_map(tp, cell_gpart_to_mesh_CIC_mapper, (void*)local_cells,
                   nr_local_cells, sizeof(int), threadpool_auto_chunk_size,
                   (void*)&data);
  }

  if (verbose)
    message("Gpart assignment took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#ifdef WITH_MPI

  MPI_Barrier(MPI_COMM_WORLD);
  tic = getticks();

  /* Merge everybody's share of the density mesh */
  MPI_Allreduce(MPI_IN_PLACE, rho, N * N * N, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  if (verbose)
    message("Mesh MPI-reduction took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#endif

  // message("\n\n\n DENSITY");
  // print_array(rho, N);

  tic = getticks();

  /* Fourier transform to go to magic-land */
  fftw_execute(forward_plan);

  if (verbose)
    message("Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* frho now contains the Fourier transform of the density field */
  /* frho contains NxNx(N/2+1) complex numbers */

  tic = getticks();

  /* Now de-convolve the CIC kernel and apply the Green function */
  mesh_apply_Green_function(tp, frho, /*slice_offset=*/0, /*slice_width=*/N,
                            /* mesh_size=*/N, r_s, region_size);

  if (verbose)
    message("Applying Green function took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* If using linear response neutrinos, apply the response to the mesh */
  if (s->e->neutrino_properties->use_linear_response) {
    neutrino_response_compute(s, mesh, tp, frho, /*slice_offset=*/0,
                              /*slice_width=*/N, verbose);

    if (verbose)
      message("Applying neutrino response took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());

    tic = getticks();
  }

  /* Fourier transform to come back from magic-land */
  fftw_execute(inverse_plan);

  if (verbose)
    message("Reverse Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* rho now contains the potential */
  /* This array is now again NxNxN real numbers */

  /* Let's store it in the structure */
  mesh->potential_global = rho;

  /* message("\n\n\n POTENTIAL"); */
  /* print_array(mesh->potential_global, N); */

  tic = getticks();

  /* Gather the mesh shared information to be used by the threads */
  data.cells = s->cells_top;
  data.rho = NULL;
  data.potential = mesh->potential_global;
  data.N = N;
  data.fac = cell_fac;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];
  data.const_G = s->e->physical_constants->const_newton_G;

  if (nr_local_cells == 0) {

    /* We don't have a cell infrastructure in place so we need to
     * directly loop over the particles */
    threadpool_map(tp, mesh_to_gpart_CIC_mapper, s->gparts, s->nr_gparts,
                   sizeof(struct gpart), threadpool_auto_chunk_size,
                   (void*)&data);

  } else { /* Normal case */

    /* Do a parallel CIC mesh interpolation onto the gparts but only using
       the local top-level cells */
    threadpool_map(tp, cell_mesh_to_gpart_CIC_mapper, (void*)local_cells,
                   nr_local_cells, sizeof(int), threadpool_auto_chunk_size,
                   (void*)&data);
  }

  if (verbose)
    message("Computing mesh accelerations took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Clean-up the mess */
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);
  memuse_log_allocation("fftw_frho", frho, 0, 0);
  fftw_free(frho);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/**
 * @brief Initialises the mesh used for the long-range periodic forces
 *
 * @param mesh The #pm_mesh to initialise.
 * @param props The propoerties of the gravity scheme.
 * @param dim The (comoving) side-lengths of the simulation volume.
 * @param nr_threads The number of threads on this MPI rank.
 */
void pm_zoom_mesh_init(struct pm_mesh* mesh, const struct gravity_props* props,
                       int nr_threads) {

#ifdef HAVE_FFTW

  const int N = props->zoom_props->mesh_size;
  const double region_size = props->zoom_props->dim;

  mesh->nr_threads = nr_threads;
  mesh->periodic = 1;
  mesh->N = N;
  mesh->distributed_mesh = props->distributed_mesh;
  mesh->use_local_patches = props->mesh_uses_local_patches;
  mesh->dim[0] = region_size;
  mesh->dim[1] = region_size;
  mesh->dim[2] = region_size;
  mesh->cell_fac = N / region_size;
  mesh->r_s = props->a_smooth * region_size / N;
  mesh->r_s_inv = 1. / mesh->r_s;
  mesh->r_cut_max = mesh->r_s * props->r_cut_max_ratio;
  mesh->r_cut_min = mesh->r_s * props->r_cut_min_ratio;
  mesh->potential_global = NULL;
  mesh->ti_beg_mesh_last = -1;
  mesh->ti_end_mesh_last = -1;
  mesh->ti_beg_mesh_next = -1;
  mesh->ti_end_mesh_next = -1;

  if (!mesh->distributed_mesh && mesh->N > 1290)
    error(
        "Mesh too big. The number of cells is larger than 2^31. "
        "A mesh side-length > 1290 currently not supported in the high res "
        "region.");

  if (2. * mesh->r_cut_max > region_size)
    error("Mesh too small or r_cut_max too big for this box size");

  pm_mesh_allocate(mesh);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}
