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
#include "mesh_gravity.h"

/* Local includes. */
#include "active.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "gravity_properties.h"
#include "kernel_long_gravity.h"
#include "part.h"
#include "runner.h"
#include "space.h"

#ifdef HAVE_FFTW

/**
 * @brief Returns 1D index of a 3D NxNxN array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices is >= N
 * or < 0.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline)) INLINE static int row_major_id_periodic(int i,
                                                                       int j,
                                                                       int k,
                                                                       int N) {
  return (((i + N) % N) * N * N + ((j + N) % N) * N + ((k + N) % N));
}

/**
 * @brief Interpolate values from a the mesh using CIC.
 *
 * @param mesh The mesh to read from.
 * @param i The index of the cell along x
 * @param j The index of the cell along y
 * @param k The index of the cell along z
 * @param tx First CIC coefficient along x
 * @param ty First CIC coefficient along y
 * @param tz First CIC coefficient along z
 * @param dx Second CIC coefficient along x
 * @param dy Second CIC coefficient along y
 * @param dz Second CIC coefficient along z
 */
__attribute__((always_inline)) INLINE static double CIC_get(
    double mesh[6][6][6], int i, int j, int k, double tx, double ty, double tz,
    double dx, double dy, double dz) {

  double temp;
  temp = mesh[i + 0][j + 0][k + 0] * tx * ty * tz;
  temp += mesh[i + 0][j + 0][k + 1] * tx * ty * dz;
  temp += mesh[i + 0][j + 1][k + 0] * tx * dy * tz;
  temp += mesh[i + 0][j + 1][k + 1] * tx * dy * dz;
  temp += mesh[i + 1][j + 0][k + 0] * dx * ty * tz;
  temp += mesh[i + 1][j + 0][k + 1] * dx * ty * dz;
  temp += mesh[i + 1][j + 1][k + 0] * dx * dy * tz;
  temp += mesh[i + 1][j + 1][k + 1] * dx * dy * dz;

  return temp;
}

/**
 * @brief Interpolate a value to a mesh using CIC.
 *
 * @param mesh The mesh to write to
 * @param N The side-length of the mesh
 * @param i The index of the cell along x
 * @param j The index of the cell along y
 * @param k The index of the cell along z
 * @param tx First CIC coefficient along x
 * @param ty First CIC coefficient along y
 * @param tz First CIC coefficient along z
 * @param dx Second CIC coefficient along x
 * @param dy Second CIC coefficient along y
 * @param dz Second CIC coefficient along z
 * @param value The value to interpolate.
 */
__attribute__((always_inline)) INLINE static void CIC_set(
    double* mesh, int N, int i, int j, int k, double tx, double ty, double tz,
    double dx, double dy, double dz, double value) {

  /* Classic CIC interpolation */
  atomic_add_d(&mesh[row_major_id_periodic(i + 0, j + 0, k + 0, N)],
               value * tx * ty * tz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 0, j + 0, k + 1, N)],
               value * tx * ty * dz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 0, j + 1, k + 0, N)],
               value * tx * dy * tz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 0, j + 1, k + 1, N)],
               value * tx * dy * dz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 1, j + 0, k + 0, N)],
               value * dx * ty * tz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 1, j + 0, k + 1, N)],
               value * dx * ty * dz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 1, j + 1, k + 0, N)],
               value * dx * dy * tz);
  atomic_add_d(&mesh[row_major_id_periodic(i + 1, j + 1, k + 1, N)],
               value * dx * dy * dz);
}

/**
 * @brief Assigns a given #gpart to a density mesh using the CIC method.
 *
 * @param gp The #gpart.
 * @param rho The density mesh.
 * @param N the size of the mesh along one axis.
 * @param fac The width of a mesh cell.
 * @param dim The dimensions of the simulation box.
 */
INLINE static void gpart_to_mesh_CIC(const struct gpart* gp, double* rho, int N,
                                     double fac, const double dim[3]) {

  /* Box wrap the multipole's position */
  const double pos_x = box_wrap(gp->x[0], 0., dim[0]);
  const double pos_y = box_wrap(gp->x[1], 0., dim[1]);
  const double pos_z = box_wrap(gp->x[2], 0., dim[2]);

  /* Workout the CIC coefficients */
  int i = (int)(fac * pos_x);
  if (i >= N) i = N - 1;
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int)(fac * pos_y);
  if (j >= N) j = N - 1;
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;

  int k = (int)(fac * pos_z);
  if (k >= N) k = N - 1;
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i >= N) error("Invalid gpart position in x");
  if (j < 0 || j >= N) error("Invalid gpart position in y");
  if (k < 0 || k >= N) error("Invalid gpart position in z");
#endif

  const double mass = gp->mass;

  /* CIC ! */
  CIC_set(rho, N, i, j, k, tx, ty, tz, dx, dy, dz, mass);
}

/**
 * @brief Assigns all the #gpart of a #cell to a density mesh using the CIC
 * method.
 *
 * @param c The #cell.
 * @param rho The density mesh.
 * @param N the size of the mesh along one axis.
 * @param fac The width of a mesh cell.
 * @param dim The dimensions of the simulation box.
 */
void cell_gpart_to_mesh_CIC(const struct cell* c, double* rho, int N,
                            double fac, const double dim[3]) {
  const int gcount = c->grav.count;
  const struct gpart* gparts = c->grav.parts;

  /* Assign all the gpart of that cell to the mesh */
  for (int i = 0; i < gcount; ++i)
    gpart_to_mesh_CIC(&gparts[i], rho, N, fac, dim);
}

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct cic_mapper_data {
  const struct cell* cells;
  double* rho;
  int N;
  double fac;
  double dim[3];
};

/**
 * @brief Threadpool mapper function for the mesh CIC assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
void cell_gpart_to_mesh_CIC_mapper(void* map_data, int num, void* extra) {

  /* Unpack the shared information */
  const struct cic_mapper_data* data = (struct cic_mapper_data*)extra;
  const struct cell* cells = data->cells;
  double* rho = data->rho;
  const int N = data->N;
  const double fac = data->fac;
  const double dim[3] = {data->dim[0], data->dim[1], data->dim[2]};

  /* Pointer to the chunk to be processed */
  int* local_cells = (int*)map_data;

  // MATTHIEU: This could in principle be improved by creating a local mesh
  //           with just the extent required for the cell. Assignment can
  //           then be done without atomics. That local mesh is then added
  //           atomically to the global one.

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell* c = &cells[local_cells[i]];

    /* Assign this cell's content to the mesh */
    cell_gpart_to_mesh_CIC(c, rho, N, fac, dim);
  }
}

/**
 * @brief Computes the potential on a gpart from a given mesh using the CIC
 * method.
 *
 * Debugging routine.
 *
 * @param gp The #gpart.
 * @param pot The potential mesh.
 * @param N the size of the mesh along one axis.
 * @param fac width of a mesh cell.
 * @param dim The dimensions of the simulation box.
 */
void mesh_to_gparts_CIC(struct gpart* gp, const double* pot, int N, double fac,
                        const double dim[3]) {

  /* Box wrap the gpart's position */
  const double pos_x = box_wrap(gp->x[0], 0., dim[0]);
  const double pos_y = box_wrap(gp->x[1], 0., dim[1]);
  const double pos_z = box_wrap(gp->x[2], 0., dim[2]);

  int i = (int)(fac * pos_x);
  if (i >= N) i = N - 1;
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int)(fac * pos_y);
  if (j >= N) j = N - 1;
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;

  int k = (int)(fac * pos_z);
  if (k >= N) k = N - 1;
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i >= N) error("Invalid gpart position in x");
  if (j < 0 || j >= N) error("Invalid gpart position in y");
  if (k < 0 || k >= N) error("Invalid gpart position in z");
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (gp->a_grav_PM[0] != 0. || gp->potential_PM != 0.)
    error("Particle with non-initalised stuff");
#endif

  /* First, copy the necessary part of the mesh for stencil operations */
  /* This includes box-wrapping in all 3 dimensions. */
  double phi[6][6][6];
  for (int iii = -2; iii <= 3; ++iii) {
    for (int jjj = -2; jjj <= 3; ++jjj) {
      for (int kkk = -2; kkk <= 3; ++kkk) {
        phi[iii + 2][jjj + 2][kkk + 2] =
            pot[row_major_id_periodic(i + iii, j + jjj, k + kkk, N)];
      }
    }
  }

  /* Some local accumulators */
  double p = 0.;
  double a[3] = {0.};

  /* Indices of (i,j,k) in the local copy of the mesh */
  const int ii = 2, jj = 2, kk = 2;

  /* Simple CIC for the potential itself */
  p += CIC_get(phi, ii, jj, kk, tx, ty, tz, dx, dy, dz);

  /* ---- */

  /* 5-point stencil along each axis for the accelerations */
  a[0] += (1. / 12.) * CIC_get(phi, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (2. / 3.) * CIC_get(phi, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] += (2. / 3.) * CIC_get(phi, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (1. / 12.) * CIC_get(phi, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);

  a[1] += (1. / 12.) * CIC_get(phi, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (2. / 3.) * CIC_get(phi, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] += (2. / 3.) * CIC_get(phi, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (1. / 12.) * CIC_get(phi, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);

  a[2] += (1. / 12.) * CIC_get(phi, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  a[2] -= (2. / 3.) * CIC_get(phi, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  a[2] += (2. / 3.) * CIC_get(phi, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  a[2] -= (1. / 12.) * CIC_get(phi, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);

  /* ---- */

  /* Store things back */
  gravity_add_comoving_potential(gp, p);
  gp->a_grav[0] += fac * a[0];
  gp->a_grav[1] += fac * a[1];
  gp->a_grav[2] += fac * a[2];
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  gp->potential_PM = p;
  gp->a_grav_PM[0] = fac * a[0];
  gp->a_grav_PM[1] = fac * a[1];
  gp->a_grav_PM[2] = fac * a[2];
#endif
}

#endif

/**
 * @brief Compute the potential, including periodic correction on the mesh.
 *
 * Interpolates the top-level multipoles on-to a mesh, move to Fourier space,
 * compute the potential including short-range correction and move back
 * to real space. We use CIC for the interpolation.
 *
 * Note that there is no multiplication by G_newton at this stage.
 *
 * @param mesh The #pm_mesh used to store the potential.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
void pm_mesh_compute_potential(struct pm_mesh* mesh, const struct space* s,
                               struct threadpool* tp, int verbose) {

#ifdef HAVE_FFTW

  const double r_s = mesh->r_s;
  const double box_size = s->dim[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int* local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  if (r_s <= 0.) error("Invalid value of a_smooth");
  if (mesh->dim[0] != dim[0] || mesh->dim[1] != dim[1] ||
      mesh->dim[2] != dim[2])
    error("Domain size does not match the value stored in the space.");

  /* Some useful constants */
  const int N = mesh->N;
  const int N_half = N / 2;
  const double cell_fac = N / box_size;

  /* Use the memory allocated for the potential to temporarily store rho */
  double* restrict rho = mesh->potential;
  if (rho == NULL) error("Error allocating memory for density mesh");
  bzero(rho, N * N * N * sizeof(double));

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

  /* Gather the mesh shared information to be used by the threads */
  struct cic_mapper_data data;
  data.cells = s->cells_top;
  data.rho = rho;
  data.N = N;
  data.fac = cell_fac;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];

  /* Do a parallel CIC mesh assignment of the gparts but only using
     the local top-level cells */
  threadpool_map(tp, cell_gpart_to_mesh_CIC_mapper, (void*)local_cells,
                 nr_local_cells, sizeof(int), 0, (void*)&data);

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
    message("Mesh comunication took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#endif

  /* message("\n\n\n DENSITY"); */
  /* print_array(rho, N); */

  tic = getticks();

  /* Fourier transform to go to magic-land */
  fftw_execute(forward_plan);

  if (verbose)
    message("Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* frho now contains the Fourier transform of the density field */
  /* frho contains NxNx(N/2+1) complex numbers */

  tic = getticks();

  /* Some common factors */
  const double green_fac = -1. / (M_PI * box_size);
  const double a_smooth2 = 4. * M_PI * M_PI * r_s * r_s / (box_size * box_size);
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
        const double sinc_kz_inv = (kz != 0) ? fz / (sin(fz) + FLT_MIN) : 1.;

        /* Norm of vector in Fourier space */
        const double k2 = (kx_d * kx_d + ky_d * ky_d + kz_d * kz_d);

        /* Avoid FPEs... */
        if (k2 == 0.) continue;

        /* Green function */
        double W = 1.;
        fourier_kernel_long_grav_eval(k2 * a_smooth2, &W);
        const double green_cor = green_fac * W / (k2 + FLT_MIN);

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

  if (verbose)
    message("Applying Green function took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Fourier transform to come back from magic-land */
  fftw_execute(inverse_plan);

  if (verbose)
    message("Backwards Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* rho now contains the potential */
  /* This array is now again NxNxN real numbers */

  /* Let's store it in the structure */
  mesh->potential = rho;

  /* message("\n\n\n POTENTIAL"); */
  /* print_array(potential, N); */

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
 * @brief Interpolate the forces and potential from the mesh to the #gpart.
 *
 * We use CIC interpolation. The resulting accelerations and potential must
 * be multiplied by G_newton.
 *
 * @param mesh The #pm_mesh (containing the potential) to interpolate from.
 * @param e The #engine (to check active status).
 * @param gparts The #gpart to interpolate to.
 * @param gcount The number of #gpart.
 */
void pm_mesh_interpolate_forces(const struct pm_mesh* mesh,
                                const struct engine* e, struct gpart* gparts,
                                int gcount) {

#ifdef HAVE_FFTW

  const int N = mesh->N;
  const double cell_fac = mesh->cell_fac;
  const double* potential = mesh->potential;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};

  /* Get the potential from the mesh to the active gparts using CIC */
  for (int i = 0; i < gcount; ++i) {
    struct gpart* gp = &gparts[i];

    if (gpart_is_active(gp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (gp->ti_drift != e->ti_current)
        error("gpart not drifted to current time");

      /* Check that the particle was initialised */
      if (gp->initialised == 0)
        error("Adding forces to an un-initialised gpart.");
#endif

      mesh_to_gparts_CIC(gp, potential, N, cell_fac, dim);
    }
  }
#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/**
 * @brief Initialisses the mesh used for the long-range periodic forces
 *
 * @param mesh The #pm_mesh to initialise.
 * @param props The propoerties of the gravity scheme.
 * @param dim The (comoving) side-lengths of the simulation volume.
 * @param nr_threads The number of threads on this MPI rank.
 */
void pm_mesh_init(struct pm_mesh* mesh, const struct gravity_props* props,
                  double dim[3], int nr_threads) {

#ifdef HAVE_FFTW

  if (dim[0] != dim[1] || dim[0] != dim[2])
    error("Doing mesh-gravity on a non-cubic domain");

  const int N = props->mesh_size;
  const double box_size = dim[0];

  mesh->nr_threads = nr_threads;
  mesh->periodic = 1;
  mesh->N = N;
  mesh->dim[0] = dim[0];
  mesh->dim[1] = dim[1];
  mesh->dim[2] = dim[2];
  mesh->cell_fac = N / box_size;
  mesh->r_s = props->a_smooth * box_size / N;
  mesh->r_s_inv = 1. / mesh->r_s;
  mesh->r_cut_max = mesh->r_s * props->r_cut_max_ratio;
  mesh->r_cut_min = mesh->r_s * props->r_cut_min_ratio;

  if (2. * mesh->r_cut_max > box_size)
    error("Mesh too small or r_cut_max too big for this box size");

#ifdef HAVE_THREADED_FFTW
  /* Initialise the thread-parallel FFTW version */
  if (N >= 64) {
    fftw_init_threads();
    fftw_plan_with_nthreads(nr_threads);
  }
#endif

  /* Allocate the memory for the combined density and potential array */
  mesh->potential = (double*)fftw_malloc(sizeof(double) * N * N * N);
  if (mesh->potential == NULL)
    error("Error allocating memory for the long-range gravity mesh.");
  memuse_log_allocation("fftw_mesh.potential", mesh->potential, 1,
                        sizeof(double) * N * N * N);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/**
 * @brief Initialises the mesh for the case where we don't do mesh gravity
 * calculations
 *
 * Crucially this set the 'periodic' propoerty to 0 and all the relevant values
 * to a
 * state where all calculations will default to pure non-periodic Newtonian.
 *
 * @param mesh The #pm_mesh to initialise.
 * @param dim The (comoving) side-lengths of the simulation volume.
 */
void pm_mesh_init_no_mesh(struct pm_mesh* mesh, double dim[3]) {

  bzero(mesh, sizeof(struct pm_mesh));

  /* Fill in non-zero properties */
  mesh->dim[0] = dim[0];
  mesh->dim[1] = dim[1];
  mesh->dim[2] = dim[2];
  mesh->r_s = FLT_MAX;
  mesh->r_cut_min = FLT_MAX;
  mesh->r_cut_max = FLT_MAX;
}

/**
 * @brief Frees the memory allocated for the long-range mesh.
 */
void pm_mesh_clean(struct pm_mesh* mesh) {

#ifdef HAVE_THREADED_FFTW
  fftw_cleanup_threads();
#endif

  if (mesh->potential) {
    memuse_log_allocation("fftw_mesh.potential", mesh->potential, 0, 0);
    free(mesh->potential);
  }
  mesh->potential = 0;
}

/**
 * @brief Write a #pm_mesh struct to the given FILE as a stream of bytes.
 *
 * @param mesh the struct
 * @param stream the file stream
 */
void pm_mesh_struct_dump(const struct pm_mesh* mesh, FILE* stream) {
  restart_write_blocks((void*)mesh, sizeof(struct pm_mesh), 1, stream,
                       "gravity", "gravity props");
}

/**
 * @brief Restore a #pm_mesh struct from the given FILE as a stream of
 * bytes.
 *
 * @param mesh the struct
 * @param stream the file stream
 */
void pm_mesh_struct_restore(struct pm_mesh* mesh, FILE* stream) {

  restart_read_blocks((void*)mesh, sizeof(struct pm_mesh), 1, stream, NULL,
                      "gravity props");

  if (mesh->periodic) {

#ifdef HAVE_FFTW
    const int N = mesh->N;

#ifdef HAVE_THREADED_FFTW
    /* Initialise the thread-parallel FFTW version */
    if (N >= 64) {
      fftw_init_threads();
      fftw_plan_with_nthreads(mesh->nr_threads);
    }
#endif

    /* Allocate the memory for the combined density and potential array */
    mesh->potential = (double*)fftw_malloc(sizeof(double) * N * N * N);
    if (mesh->potential == NULL)
      error("Error allocating memory for the long-range gravity mesh.");
    memuse_log_allocation("fftw_mesh.potential", mesh->potential, 1,
                          sizeof(double) * N * N * N);
#else
    error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
  }
}
