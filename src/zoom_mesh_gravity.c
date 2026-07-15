/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will J. Roper (w.roper@sussex.ac.uk).
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
#endif
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "zoom_mesh_gravity.h"

/* Local headers. */
#include "atomic.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "gravity_properties.h"
#include "inline.h"
#include "kernel_long_gravity.h"
#include "memuse.h"
#include "parser.h"
#include "part.h"
#include "restart.h"
#include "row_major_id.h"
#include "space.h"
#include "threadpool.h"

/* Standard headers. */
#include <float.h>
#include <math.h>
#include <string.h>

/* Local function prototypes. */
static void zoom_mesh_allocate(struct zoom_pm_mesh *mesh);
static void zoom_mesh_free(struct zoom_pm_mesh *mesh);
static void zoom_mesh_init_no_mesh(struct zoom_pm_mesh *mesh);
static int zoom_mesh_cell_is_inside(const struct zoom_pm_mesh *mesh,
                                    const struct cell *c,
                                    const int use_max_dx);

/**
 * @brief Initialize a zoom mesh structure for the case where it is disabled.
 *
 * @param mesh The #zoom_pm_mesh to initialise.
 */
static void zoom_mesh_init_no_mesh(struct zoom_pm_mesh *mesh) {

  /* Zero out everything. */
  bzero(mesh, sizeof(struct zoom_pm_mesh));
}

/**
 * @brief Initialize the high-resolution zoom gravity mesh geometry.
 *
 * The mesh is disabled unless the parameter Gravity:zoom_mesh is non-zero. When
 * enabled, the mesh covers the current zoom region plus a buffer measured in
 * background top-level cell shells.
 *
 * @param mesh The #zoom_pm_mesh to initialise.
 * @param params The parsed parameter file.
 * @param props The properties of the gravity scheme.
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void zoom_mesh_init(struct zoom_pm_mesh *mesh, struct swift_params *params,
                    const struct gravity_props *props, const struct space *s,
                    const int verbose) {

  /* Start with a clean mesh which is identical to there being no mesh. */
  zoom_mesh_init_no_mesh(mesh);

  /* Simulations without self-gravity never create a zoom gravity mesh. */
  if (params == NULL) return;

  /* Check if the zoom mesh is enabled, if not we are done. */
  const int enabled = parser_get_opt_param_int(params, "Gravity:zoom_mesh", 0);
  if (!enabled) return;

  /* Getting here without a zoom region is nonsensical. */
  if (!s->with_zoom_region) error("Gravity:zoom_mesh requires a zoom region.");

  /* Ensure that the box is periodic */
  if (!s->periodic)
    error(
        "Can't run with a high resolution mesh without periodic boundary "
        "conditions to support the low resolution mesh.");

  /* Ensure that the zoom region has been initialized. */
  if (s->zoom_props == NULL)
    error("Gravity:zoom_mesh requested before zoom region initialization.");

  /* Unpack the active mesh parameters and sanity check them. The user-provided
   * side length refers to the inner N^3 zoom mesh, not the zero-padded FFT
   * domain used by the zero padded Hockney-Eastwood method. */
  const int side_length =
      parser_get_param_int(params, "Gravity:zoom_mesh_side_length");
  if (side_length < 8)
    error("Gravity:zoom_mesh_side_length must be at least 8.");

  /* This implementation stores the full zero-padded (2N)^3 FFT mesh on every
   * rank. Keep the doubled mesh side length within the same practical serial
   * FFT limit used by the regular non-distributed PM mesh. */
  if (2 * side_length > 1290)
    error(
        "Gravity:zoom_mesh_side_length is too large for the current "
        "zero-padded zoom mesh implementation. The active side length N is "
        "%d, but the zero padded Hockney-Eastwood FFT uses side length 2N=%d. "
        "Use N <= 645 until a distributed zoom mesh is implemented.",
        side_length, 2 * side_length);

  const double r_cut_max_ratio =
      parser_get_param_double(params, "Gravity:zoom_r_cut_max");
  if (r_cut_max_ratio <= 0.)
    error("Gravity:zoom_r_cut_max must be positive.");

  const double r_cut_min_ratio =
      parser_get_opt_param_double(params, "Gravity:zoom_r_cut_min", 0.);
  if (r_cut_min_ratio < 0.) error("Gravity:zoom_r_cut_min must be >= 0.");
  if (r_cut_min_ratio > r_cut_max_ratio)
    error("Gravity:zoom_r_cut_min must be <= Gravity:zoom_r_cut_max.");

  const int buffer_cells =
      parser_get_opt_param_int(params, "Gravity:zoom_mesh_bkg_buffer_cells", 1);
  if (buffer_cells < 1)
    error("Gravity:zoom_mesh_bkg_buffer_cells must be >= 1.");

  if (s->zoom_props->nr_bkg_cells <= 0 || s->zoom_props->bkg_cells_top == NULL)
    error("Gravity:zoom_mesh requires background top-level cells.");
  if (s->zoom_props->nr_zoom_cells <= 0 ||
      s->zoom_props->zoom_cells_top == NULL)
    error("Gravity:zoom_mesh requires zoom top-level cells.");

  const double zoom_cell_width[3] = {
      s->zoom_props->zoom_cells_top[0].width[0],
      s->zoom_props->zoom_cells_top[0].width[1],
      s->zoom_props->zoom_cells_top[0].width[2]};
  if (zoom_cell_width[0] != zoom_cell_width[1] ||
      zoom_cell_width[0] != zoom_cell_width[2])
    error("Gravity:zoom_mesh requires cubic zoom top-level cells.");
  const double zoom_r_s = props->a_smooth * zoom_cell_width[0];

  /* Compute the buffer width in each dimension. */
  const double buffer[3] = {
      buffer_cells * s->zoom_props->bkg_cells_top[0].width[0],
      buffer_cells * s->zoom_props->bkg_cells_top[0].width[1],
      buffer_cells * s->zoom_props->bkg_cells_top[0].width[2]};

  /* Fill in the mesh structure. */
  mesh->enabled = 1;
  mesh->r_s = zoom_r_s;
  mesh->r_s_inv = 1. / mesh->r_s;
  mesh->r_cut_max = mesh->r_s * r_cut_max_ratio;
  mesh->r_cut_min = mesh->r_s * r_cut_min_ratio;
  mesh->potential_global = NULL;

  for (int i = 0; i < 3; ++i) {
    mesh->N[i] = side_length;
    mesh->loc[i] = s->zoom_props->region_lower_bounds[i] - buffer[i];
    mesh->dim[i] = s->zoom_props->dim[i] + 2. * buffer[i];
    mesh->cell_fac[i] = side_length / mesh->dim[i];
  }

  /* Report the mesh geometry to the user if requested. */
  if (verbose)
    message(
        "Zoom gravity mesh: enabled with N=[%d %d %d], loc=[%.6e %.6e %.6e], "
        "dim=[%.6e %.6e %.6e], buffer_cells=%d, r_cut_max=%.6e.",
        mesh->N[0], mesh->N[1], mesh->N[2], mesh->loc[0], mesh->loc[1],
        mesh->loc[2], mesh->dim[0], mesh->dim[1], mesh->dim[2], buffer_cells,
        mesh->r_cut_max);

  /* Allocate the zero-padded FFT mesh. */
  zoom_mesh_allocate(mesh);
}

/**
 * @brief Clean the zoom mesh structure.
 *
 * @param mesh The #zoom_pm_mesh to clean.
 */
void zoom_mesh_clean(struct zoom_pm_mesh *mesh) {

  /* Nothing to do if the mesh was never created. */
  if (mesh == NULL) return;

  /* Release any allocated FFT storage and reset to a no-mesh state. */
  zoom_mesh_free(mesh);
  zoom_mesh_init_no_mesh(mesh);
}

/**
 * @brief Allocates the zero-padded zoom mesh grid.
 *
 * @param mesh The #zoom_pm_mesh structure.
 */
static void zoom_mesh_allocate(struct zoom_pm_mesh *mesh) {

  /* Nothing to allocate for a disabled or absent mesh. */
  if (mesh == NULL || !mesh->enabled) return;

#ifdef HAVE_FFTW

  /* The zero padded Hockney-Eastwood method embeds the active N^3 mesh in a
   * (2N)^3 FFT domain. */
  const int N = mesh->N[0];
  const int M = 2 * N;
  const size_t nr_mesh_cells = (size_t)M * M * M;

  /* Allocate storage used for both the density and potential fields. */
  mesh->potential_global =
      (double *)fftw_malloc(sizeof(double) * nr_mesh_cells);
  if (mesh->potential_global == NULL)
    error("Error allocating memory for the zoom gravity mesh.");

  /* Record the FFTW allocation in the memory-use logger. */
  memuse_log_allocation("fftw_zoom_mesh.potential", mesh->potential_global, 1,
                        sizeof(double) * nr_mesh_cells);
#else
  error("No FFTW library found. Cannot compute zoom mesh gravity.");
#endif
}

/**
 * @brief Frees the zero-padded zoom mesh grid.
 *
 * @param mesh The #zoom_pm_mesh structure.
 */
static void zoom_mesh_free(struct zoom_pm_mesh *mesh) {

  /* Nothing to release if the mesh was never allocated. */
  if (mesh == NULL || mesh->potential_global == NULL) return;

#ifdef HAVE_FFTW

  /* Update memory accounting before releasing the FFTW allocation. */
  memuse_log_allocation("fftw_zoom_mesh.potential", mesh->potential_global, 0,
                        0);
  fftw_free(mesh->potential_global);

  /* Avoid leaving a dangling pointer in restart/output code paths. */
  mesh->potential_global = NULL;
#else
  error("No FFTW library found. Cannot compute zoom mesh gravity.");
#endif
}

/**
 * @brief Write a #zoom_pm_mesh struct to the given FILE as a stream of bytes.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param stream The file stream.
 */
void zoom_mesh_struct_dump(const struct zoom_pm_mesh *mesh, FILE *stream) {

  /* The FFT work array is not serialized; it is reallocated on restore. */
  restart_write_blocks((void *)mesh, sizeof(struct zoom_pm_mesh), 1, stream,
                       "gravity", "zoom gravity mesh");
}

/**
 * @brief Restore a #zoom_pm_mesh struct from the given FILE as a stream of
 * bytes.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param stream The file stream.
 */
void zoom_mesh_struct_restore(struct zoom_pm_mesh *mesh, FILE *stream) {

  /* Restore the POD state first. */
  restart_read_blocks((void *)mesh, sizeof(struct zoom_pm_mesh), 1, stream,
                      NULL, "zoom gravity mesh");

  /* Recreate the transient FFT storage from the restored geometry. */
  mesh->potential_global = NULL;
  if (mesh->enabled) zoom_mesh_allocate(mesh);
}

#ifdef HAVE_FFTW

/**
 * @brief Returns the 1D index of a 3D MxMxM zoom mesh array.
 *
 * Unlike the regular PM mesh helpers, this does not wrap periodically. The
 * zero padded Hockney-Eastwood embedding stores negative Green-function offsets
 * explicitly in the upper half of the zero-padded domain.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param M Size of the zero-padded array along one axis.
 * @return The row-major 1D index.
 */
__attribute__((always_inline, const)) INLINE static size_t zoom_mesh_index(
    const int i, const int j, const int k, const int M) {

  /* Standard row-major indexing without periodic wrapping. */
  return ((size_t)i * M + (size_t)j) * M + (size_t)k;
}

/**
 * @brief Interpolate a value to the zero-padded zoom mesh using CIC.
 *
 * @param mesh The mesh to write to.
 * @param M The size of the zero-padded mesh along one axis.
 * @param i The index of the cell along x.
 * @param j The index of the cell along y.
 * @param k The index of the cell along z.
 * @param tx First CIC coefficient along x.
 * @param ty First CIC coefficient along y.
 * @param tz First CIC coefficient along z.
 * @param dx Second CIC coefficient along x.
 * @param dy Second CIC coefficient along y.
 * @param dz Second CIC coefficient along z.
 * @param value The value to interpolate.
 */
__attribute__((always_inline)) INLINE static void zoom_mesh_CIC_set(
    double *mesh, const int M, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz, const double value) {

  /* Classic CIC assignment to the eight surrounding mesh vertices. */
  atomic_add_d(&mesh[zoom_mesh_index(i + 0, j + 0, k + 0, M)],
               value * tx * ty * tz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 0, j + 0, k + 1, M)],
               value * tx * ty * dz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 0, j + 1, k + 0, M)],
               value * tx * dy * tz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 0, j + 1, k + 1, M)],
               value * tx * dy * dz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 1, j + 0, k + 0, M)],
               value * dx * ty * tz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 1, j + 0, k + 1, M)],
               value * dx * ty * dz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 1, j + 1, k + 0, M)],
               value * dx * dy * tz);
  atomic_add_d(&mesh[zoom_mesh_index(i + 1, j + 1, k + 1, M)],
               value * dx * dy * dz);
}

/**
 * @brief Interpolate values from the zero-padded zoom mesh using CIC.
 *
 * @param mesh The mesh to read from.
 * @param M The size of the zero-padded mesh along one axis.
 * @param i The index of the cell along x.
 * @param j The index of the cell along y.
 * @param k The index of the cell along z.
 * @param tx First CIC coefficient along x.
 * @param ty First CIC coefficient along y.
 * @param tz First CIC coefficient along z.
 * @param dx Second CIC coefficient along x.
 * @param dy Second CIC coefficient along y.
 * @param dz Second CIC coefficient along z.
 * @return The interpolated value.
 */
__attribute__((always_inline, const)) INLINE static double zoom_mesh_CIC_get(
    const double *mesh, const int M, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz) {

  /* Classic CIC interpolation from the eight surrounding mesh vertices. */
  double value = mesh[zoom_mesh_index(i + 0, j + 0, k + 0, M)] * tx * ty * tz;
  value += mesh[zoom_mesh_index(i + 0, j + 0, k + 1, M)] * tx * ty * dz;
  value += mesh[zoom_mesh_index(i + 0, j + 1, k + 0, M)] * tx * dy * tz;
  value += mesh[zoom_mesh_index(i + 0, j + 1, k + 1, M)] * tx * dy * dz;
  value += mesh[zoom_mesh_index(i + 1, j + 0, k + 0, M)] * dx * ty * tz;
  value += mesh[zoom_mesh_index(i + 1, j + 0, k + 1, M)] * dx * ty * dz;
  value += mesh[zoom_mesh_index(i + 1, j + 1, k + 0, M)] * dx * dy * tz;
  value += mesh[zoom_mesh_index(i + 1, j + 1, k + 1, M)] * dx * dy * dz;

  return value;
}

/**
 * @brief Shared information about the zoom mesh to be used by all threads in
 * the pool.
 */
struct zoom_cic_mapper_data {

  /*! The zoom mesh geometry. */
  const struct zoom_pm_mesh *mesh;

  /*! The density mesh to write to. */
  double *rho;

  /*! The potential mesh to read from. */
  const double *potential;

  /*! The size of the zero-padded mesh along one axis. */
  int M;

  /*! Newton's constant in internal units. */
  float const_G;
};

/**
 * @brief Compute the zoom mesh CIC index and coefficients for a #gpart.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param gp The #gpart.
 * @param ind The integer mesh index along each axis (return).
 * @param d The second CIC coefficient along each axis (return).
 * @param t The first CIC coefficient along each axis (return).
 * @param need_stencil Should the point be far enough from the active mesh
 * edges for the five-point force stencil?
 * @return 1 if the particle can be represented on the mesh, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int zoom_mesh_gpart_index(
    const struct zoom_pm_mesh *mesh, const struct gpart *gp, int ind[3],
    double d[3], double t[3], const int need_stencil) {

  /* Convert the particle position from box coordinates to mesh coordinates. */
  for (int axis = 0; axis < 3; ++axis) {
    const double x = (gp->x[axis] - mesh->loc[axis]) * mesh->cell_fac[axis];
    ind[axis] = (int)floor(x);
    d[axis] = x - ind[axis];
    t[axis] = 1. - d[axis];

    /* Force interpolation needs two mesh cells on either side for the
     * five-point stencil. Density assignment only needs the CIC cell. */
    if (need_stencil) {
      if (ind[axis] < 2 || ind[axis] >= mesh->N[axis] - 3) return 0;
    } else {
      if (ind[axis] < 0 || ind[axis] >= mesh->N[axis] - 1) return 0;
    }
  }

  return 1;
}

/**
 * @brief Is a cell fully covered by the high-resolution zoom mesh?
 *
 * @param mesh The #zoom_pm_mesh.
 * @param c The #cell.
 * @param use_max_dx Should we include the maximal particle displacement since
 * rebuild?
 * @return 1 if the cell is fully covered by the zoom mesh, 0 otherwise.
 */
int zoom_mesh_cell_is_covered(const struct zoom_pm_mesh *mesh,
                              const struct cell *c, const int use_max_dx) {

  /* Disabled or absent meshes cannot cover any cell. */
  if (mesh == NULL) return 0;
  if (!mesh->enabled) return 0;

  /* Defer to the geometric containment check. */
  return zoom_mesh_cell_is_inside(mesh, c, use_max_dx);
}

/**
 * @brief Are two cells fully covered by the high-resolution zoom mesh?
 *
 * @param mesh The #zoom_pm_mesh.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param use_max_dx Should we include the maximal particle displacement since
 * rebuild?
 * @return 1 if both cells are fully covered by the zoom mesh, 0 otherwise.
 */
int zoom_mesh_cells_are_covered(const struct zoom_pm_mesh *mesh,
                                const struct cell *ci, const struct cell *cj,
                                const int use_max_dx) {

  /* Both cells must be fully contained before the zoom mesh is applicable. */
  return zoom_mesh_cell_is_covered(mesh, ci, use_max_dx) &&
         zoom_mesh_cell_is_covered(mesh, cj, use_max_dx);
}

/**
 * @brief Recursively deposit covered local cells on the zoom mesh.
 *
 * @param c The #cell to deposit.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_deposit_cell(const struct cell *c,
                                   const struct zoom_cic_mapper_data *data) {

  /* Recurse until we find covered leaves. */
  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL) zoom_mesh_deposit_cell(c->progeny[k], data);
    return;
  }

  /* The zoom mesh correction follows the same cell coverage predicate as the
   * tree split. Particles in partially covered cells keep the global PM split. */
  if (!zoom_mesh_cell_is_covered(data->mesh, c, /*use_max_dx=*/0)) return;

  struct gpart *gparts = c->grav.parts;
  for (int n = 0; n < c->grav.count; ++n) {
    const struct gpart *gp = &gparts[n];
    if (gp->time_bin == time_bin_inhibited) continue;

    int ind[3];
    double d[3], t[3];
    if (!zoom_mesh_gpart_index(data->mesh, gp, ind, d, t, /*need_stencil=*/0))
      error("Covered cell contains a particle outside the zoom mesh.");

    zoom_mesh_CIC_set(data->rho, data->M, ind[0], ind[1], ind[2], t[0], t[1],
                      t[2], d[0], d[1], d[2], gp->mass);
  }
}

/**
 * @brief Threadpool mapper function for zoom mesh force interpolation.
 *
 * Interpolates the isolated zoom mesh correction to particles and adds it to
 * the existing mesh acceleration and potential fields.
 *
 * @param map_data A chunk of #gpart objects.
 * @param num The number of #gpart objects in the chunk.
 * @param extra The shared #zoom_cic_mapper_data.
 */
static void zoom_mesh_to_gpart_CIC_mapper(void *map_data, int num,
                                          void *extra) {

  /* Unpack the shared mapper data. */
  const struct zoom_cic_mapper_data *data =
      (const struct zoom_cic_mapper_data *)extra;
  const struct zoom_pm_mesh *mesh = data->mesh;
  struct gpart *gparts = (struct gpart *)map_data;
  const double *potential = data->potential;
  const int M = data->M;
  const float const_G = data->const_G;

  /* Interpolate the force correction onto all valid particles in this chunk. */
  for (int n = 0; n < num; ++n) {
    struct gpart *gp = &gparts[n];
    if (gp->time_bin == time_bin_inhibited) continue;

    /* Skip particles outside the active mesh or too close to the stencil edge. */
    int ind[3];
    double d[3], t[3];
    if (!zoom_mesh_gpart_index(mesh, gp, ind, d, t, /*need_stencil=*/1))
      continue;

    double p = 0.;
    double a[3] = {0.};

    /* Interpolate the potential using CIC. */
    p = zoom_mesh_CIC_get(potential, M, ind[0], ind[1], ind[2], t[0], t[1],
                          t[2], d[0], d[1], d[2]);

    /* Compute accelerations using the same five-point stencil as the regular
     * PM mesh, with CIC interpolation applied at each stencil point. */
    a[0] +=
        (1. / 12.) * zoom_mesh_CIC_get(potential, M, ind[0] + 2, ind[1], ind[2],
                                       t[0], t[1], t[2], d[0], d[1], d[2]);
    a[0] -=
        (2. / 3.) * zoom_mesh_CIC_get(potential, M, ind[0] + 1, ind[1], ind[2],
                                      t[0], t[1], t[2], d[0], d[1], d[2]);
    a[0] +=
        (2. / 3.) * zoom_mesh_CIC_get(potential, M, ind[0] - 1, ind[1], ind[2],
                                      t[0], t[1], t[2], d[0], d[1], d[2]);
    a[0] -=
        (1. / 12.) * zoom_mesh_CIC_get(potential, M, ind[0] - 2, ind[1], ind[2],
                                       t[0], t[1], t[2], d[0], d[1], d[2]);

    a[1] +=
        (1. / 12.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1] + 2, ind[2],
                                       t[0], t[1], t[2], d[0], d[1], d[2]);
    a[1] -=
        (2. / 3.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1] + 1, ind[2],
                                      t[0], t[1], t[2], d[0], d[1], d[2]);
    a[1] +=
        (2. / 3.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1] - 1, ind[2],
                                      t[0], t[1], t[2], d[0], d[1], d[2]);
    a[1] -=
        (1. / 12.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1] - 2, ind[2],
                                       t[0], t[1], t[2], d[0], d[1], d[2]);

    a[2] +=
        (1. / 12.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1], ind[2] + 2,
                                       t[0], t[1], t[2], d[0], d[1], d[2]);
    a[2] -=
        (2. / 3.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1], ind[2] + 1,
                                      t[0], t[1], t[2], d[0], d[1], d[2]);
    a[2] +=
        (2. / 3.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1], ind[2] - 1,
                                      t[0], t[1], t[2], d[0], d[1], d[2]);
    a[2] -=
        (1. / 12.) * zoom_mesh_CIC_get(potential, M, ind[0], ind[1], ind[2] - 2,
                                       t[0], t[1], t[2], d[0], d[1], d[2]);

    /* Add the correction to the existing mesh contribution. */
    gp->a_grav_mesh[0] += const_G * mesh->cell_fac[0] * a[0];
    gp->a_grav_mesh[1] += const_G * mesh->cell_fac[1] * a[1];
    gp->a_grav_mesh[2] += const_G * mesh->cell_fac[2] * a[2];
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
    gp->potential_mesh += const_G * p;
#endif
  }
}

/**
 * @brief Recursively interpolate the zoom mesh correction to covered cells.
 *
 * @param c The #cell to update.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_interpolate_cell(struct cell *c,
                                       const struct zoom_cic_mapper_data *data) {

  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL) zoom_mesh_interpolate_cell(c->progeny[k], data);
    return;
  }

  if (!zoom_mesh_cell_is_covered(data->mesh, c, /*use_max_dx=*/0)) return;

  zoom_mesh_to_gpart_CIC_mapper(c->grav.parts, c->grav.count, (void *)data);
}

/**
 * @brief Deposit all local covered cells on the zoom mesh.
 *
 * @param s The #space.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_deposit_cells(const struct space *s,
                                    const struct zoom_cic_mapper_data *data) {

  if (s->with_zoom_region && s->zoom_props != NULL) {
    for (int n = 0; n < s->zoom_props->nr_local_zoom_cells_with_particles; ++n) {
      const int cid = s->zoom_props->local_zoom_cells_with_particles_top[n];
      zoom_mesh_deposit_cell(&s->cells_top[cid], data);
    }
    for (int n = 0; n < s->zoom_props->nr_local_bkg_cells_with_particles; ++n) {
      const int cid = s->zoom_props->local_bkg_cells_with_particles_top[n];
      zoom_mesh_deposit_cell(&s->cells_top[cid], data);
    }
  } else if (s->cells_top != NULL) {
    for (int n = 0; n < s->nr_local_cells_with_particles; ++n) {
      const int cid = s->local_cells_with_particles_top[n];
      zoom_mesh_deposit_cell(&s->cells_top[cid], data);
    }
  }
}

/**
 * @brief Interpolate the zoom mesh correction to all local covered cells.
 *
 * @param s The #space.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_interpolate_cells(const struct space *s,
                                        const struct zoom_cic_mapper_data *data) {

  if (s->with_zoom_region && s->zoom_props != NULL) {
    for (int n = 0; n < s->zoom_props->nr_local_zoom_cells_with_particles; ++n) {
      const int cid = s->zoom_props->local_zoom_cells_with_particles_top[n];
      zoom_mesh_interpolate_cell(&s->cells_top[cid], data);
    }
    for (int n = 0; n < s->zoom_props->nr_local_bkg_cells_with_particles; ++n) {
      const int cid = s->zoom_props->local_bkg_cells_with_particles_top[n];
      zoom_mesh_interpolate_cell(&s->cells_top[cid], data);
    }
  } else if (s->cells_top != NULL) {
    for (int n = 0; n < s->nr_local_cells_with_particles; ++n) {
      const int cid = s->local_cells_with_particles_top[n];
      zoom_mesh_interpolate_cell(&s->cells_top[cid], data);
    }
  }
}

/**
 * @brief Construct the zero padded Hockney-Eastwood isolated-convolution kernel.
 *
 * The active zoom mesh has side length N but the FFT domain has side length
 * 2N. Positive Green-function offsets are stored in [0, N-1] and negative
 * offsets are stored in [2N-(N-1), 2N-1]. This makes the circular convolution
 * performed by the FFT equivalent to the isolated linear convolution on the
 * active N^3 mesh.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param kernel The zero-padded real-space kernel to fill.
 * @param s The #space.
 */
static void zoom_mesh_apply_kernel(struct zoom_pm_mesh *mesh, double *kernel,
                                   const struct space *s) {

  /* Unpack the active and zero-padded mesh sizes. */
  const int N = mesh->N[0];
  const int M = 2 * N;

  /* Unpack the split scales for the zoom and global PM meshes. */
  const double r_s_inv = mesh->r_s_inv;
  const double global_r_s_inv = s->periodic ? s->e->mesh->r_s_inv : 0.;
  const double global_r_s = global_r_s_inv > 0. ? 1. / global_r_s_inv : 0.;

  /* Fill the embedded real-space correction Green's function. */
  for (int i = 0; i < M; ++i) {
    /* Indices in the upper half of the padded domain are negative offsets. */
    const int ii = i < N ? i : i - M;
    const double dx = ii / mesh->cell_fac[0];
    for (int j = 0; j < M; ++j) {
      const int jj = j < N ? j : j - M;
      const double dy = jj / mesh->cell_fac[1];
      for (int k = 0; k < M; ++k) {
        const int kk = k < N ? k : k - M;
        const double dz = kk / mesh->cell_fac[2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);

        /* Long-range potential represented by the zoom split. */
        double zoom_long;
        if (r == 0.)
          zoom_long = -M_2_SQRTPI * 0.5 * r_s_inv;
        else
          zoom_long = -erf(0.5 * r * r_s_inv) / r;

        /* Long-range potential already represented by the global PM split. */
        double global_long = 0.;
        if (global_r_s_inv > 0.) {
          if (r == 0.)
            global_long = -M_2_SQRTPI * 0.5 / global_r_s;
          else
            global_long = -erf(0.5 * r * global_r_s_inv) / r;
        }

        /* Store only the correction that must be added to the global PM field. */
        kernel[zoom_mesh_index(i, j, k, M)] = zoom_long - global_long;
      }
    }
  }
}

/**
 * @brief Compute the zoom mesh correction to the gravity forces.
 *
 * This uses the zero padded Hockney-Eastwood method for isolated boundary
 * conditions: particle mass is deposited on the active N^3 zoom mesh, both the
 * density and Green's function are embedded in a zero-padded (2N)^3 FFT domain,
 * and the resulting circular convolution is therefore the isolated linear
 * convolution on the active mesh. The resulting correction is then interpolated
 * back to particles and added to their existing mesh acceleration and potential.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
void zoom_mesh_compute_potential(struct zoom_pm_mesh *mesh,
                                 const struct space *s, struct threadpool *tp,
                                 int verbose) {

  (void)tp;

  /* Nothing to compute if the zoom mesh is disabled. */
  if (mesh == NULL || !mesh->enabled) return;

  /* Unpack sizes and memory requirements for the zero-padded FFT domain. */
  const int N = mesh->N[0];
  const int M = 2 * N;
  const int M_half = M / 2;
  const size_t real_size = (size_t)M * M * M;
  const size_t complex_size = (size_t)M * M * (M_half + 1);
  if (mesh->potential_global == NULL) zoom_mesh_allocate(mesh);
  double *rho = mesh->potential_global;

  /* Clear the mesh before assigning the current particle distribution. */
  ticks tic = getticks();
  bzero(rho, real_size * sizeof(double));

  /* Fill in the data shared by the CIC density-assignment mapper. */
  struct zoom_cic_mapper_data data;
  data.mesh = mesh;
  data.rho = rho;
  data.potential = NULL;
  data.M = M;
  data.const_G = 0.f;

  /* Assign the local covered cells to the active zoom mesh. This mirrors the
   * cell predicate used to select the zoom split in the tree. */
  zoom_mesh_deposit_cells(s, &data);

#ifdef WITH_MPI
  /* This implementation stores the full mesh on every rank, so sum the mass
   * assignment from all ranks before doing the serial FFT. */
  MPI_Allreduce(MPI_IN_PLACE, rho, real_size, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  if (verbose)
    message("Zoom mesh gpart assignment took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Allocate the temporary real-space kernel and Fourier-space work arrays. */
  double *kernel = (double *)fftw_malloc(real_size * sizeof(double));
  fftw_complex *frho =
      (fftw_complex *)fftw_malloc(complex_size * sizeof(fftw_complex));
  fftw_complex *fkernel =
      (fftw_complex *)fftw_malloc(complex_size * sizeof(fftw_complex));
  if (kernel == NULL || frho == NULL || fkernel == NULL)
    error("Error allocating zoom mesh FFT work arrays.");

  /* Record the temporary allocations in the memory-use logger. */
  memuse_log_allocation("fftw_zoom_mesh.kernel", kernel, 1,
                        real_size * sizeof(double));
  memuse_log_allocation("fftw_zoom_mesh.frho", frho, 1,
                        complex_size * sizeof(fftw_complex));
  memuse_log_allocation("fftw_zoom_mesh.fkernel", fkernel, 1,
                        complex_size * sizeof(fftw_complex));

  /* Construct the isolated-convolution correction kernel. */
  zoom_mesh_apply_kernel(mesh, kernel, s);

  /* Prepare the forward transforms of the density and kernel and the inverse
   * transform of their product. */
  fftw_plan rho_forward =
      fftw_plan_dft_r2c_3d(M, M, M, rho, frho, FFTW_ESTIMATE);
  fftw_plan kernel_forward =
      fftw_plan_dft_r2c_3d(M, M, M, kernel, fkernel, FFTW_ESTIMATE);
  fftw_plan inverse = fftw_plan_dft_c2r_3d(M, M, M, frho, rho, FFTW_ESTIMATE);

  /* Transform density and kernel to Fourier space. */
  fftw_execute(rho_forward);
  fftw_execute(kernel_forward);

  /* Multiply the Fourier-space density and kernel. The inverse FFTW transform
   * is unnormalised, hence the explicit 1 / M^3 factor here. */
  const double norm = 1. / (double)real_size;
  for (size_t i = 0; i < complex_size; ++i) {
    const double a = frho[i][0];
    const double b = frho[i][1];
    const double c = fkernel[i][0];
    const double d = fkernel[i][1];
    frho[i][0] = norm * (a * c - b * d);
    frho[i][1] = norm * (a * d + b * c);
  }

  /* Transform back to real space; rho now stores the potential correction. */
  fftw_execute(inverse);

  /* Destroy the plans before releasing the temporary arrays. */
  fftw_destroy_plan(rho_forward);
  fftw_destroy_plan(kernel_forward);
  fftw_destroy_plan(inverse);

  /* Release all temporary FFT storage. */
  memuse_log_allocation("fftw_zoom_mesh.kernel", kernel, 0, 0);
  memuse_log_allocation("fftw_zoom_mesh.frho", frho, 0, 0);
  memuse_log_allocation("fftw_zoom_mesh.fkernel", fkernel, 0, 0);
  fftw_free(kernel);
  fftw_free(frho);
  fftw_free(fkernel);

  if (verbose)
    message("Zoom mesh convolution took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Reuse the mapper data to interpolate the potential correction back to the
   * particles as accelerations and potentials. */
  data.rho = NULL;
  data.potential = rho;
  data.const_G = s->e->physical_constants->const_newton_G;

  zoom_mesh_interpolate_cells(s, &data);

  if (verbose)
    message("Zoom mesh acceleration interpolation took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

#endif /* HAVE_FFTW */

#ifndef HAVE_FFTW
/**
 * @brief Compute the zoom mesh correction to the gravity forces.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
void zoom_mesh_compute_potential(struct zoom_pm_mesh *mesh,
                                 const struct space *s, struct threadpool *tp,
                                 int verbose) {

  /* Disabled zoom meshes must remain harmless in non-FFTW builds. */
  if (mesh == NULL || !mesh->enabled) return;

  error("No FFTW library found. Cannot compute zoom mesh gravity.");
}
#endif

/**
 * @brief Test whether a distance is safely beyond a mesh cutoff.
 *
 * Matches the tolerance convention used by the regular gravity mesh
 * predicate.
 *
 * @param min_radius2 The squared minimal pair separation.
 * @param max_distance2 The squared PM truncation cutoff distance.
 * @return 1 if the pair is safely beyond the cutoff, 0 otherwise.
 */
__attribute__((const)) INLINE static int zoom_mesh_distance_is_above_cutoff(
    const double min_radius2, const double max_distance2) {

  /* Derive a scale-aware tolerance. (128 is arbitrarily "big enough") */
  const double scale = fmax(1.0, fmax(fabs(min_radius2), fabs(max_distance2)));
  const double tol = 128.0 * DBL_EPSILON * scale;

  /* Only accept distances safely above the cutoff. */
  return min_radius2 > max_distance2 + tol;
}

/**
 * @brief Is a cell fully inside the zoom mesh bounds?
 *
 * @param mesh The #zoom_pm_mesh.
 * @param c The #cell.
 * @param use_max_dx Should we include the maximal particle displacement since
 * rebuild?
 * @return 1 if the full cell is inside the mesh, 0 otherwise.
 */
__attribute__((nonnull)) static int zoom_mesh_cell_is_inside(
    const struct zoom_pm_mesh *mesh, const struct cell *c,
    const int use_max_dx) {

  /* Check containment independently along each axis. */
  for (int i = 0; i < 3; ++i) {
    /* Mesh bounds in box coordinates, restricted to the active-mesh region
     * where the five-point force stencil is valid. The negative side needs two
     * cells for i-2. The positive side needs three cells because the i+2
     * stencil point is CIC-interpolated and therefore also reads i+3. */
    const double cell_width = 1. / mesh->cell_fac[i];
    const double mesh_min = mesh->loc[i] + 2. * cell_width;
    const double mesh_max = mesh->loc[i] + mesh->dim[i] - 3. * cell_width;

    /* Inflate the cell by the maximal particle displacement if requested. */
    const double dx_max = use_max_dx ? c->grav.multipole->dx_max[i] : 0.;
    const double cell_min = c->loc[i] - dx_max;
    const double cell_max = c->loc[i] + c->width[i] + dx_max;

    /* Allow a small geometry tolerance for roundoff in cell construction. */
    const double tol = 1e-6 * mesh->dim[i];

    /* Any excursion outside the mesh bounds disqualifies the cell. */
    if (cell_min < mesh_min - tol) return 0;
    if (cell_max > mesh_max + tol) return 0;
  }

  return 1;
}

/**
 * @brief Can a pair of cells use the zoom mesh for their interaction?
 *
 * This is the rebuild-time predicate used when creating the task graph. It is
 * deliberately conservative: only pairs fully contained in the zoom mesh and
 * safely beyond the zoom mesh cutoff can be skipped.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @return 1 if the zoom mesh can be used, 0 otherwise.
 */
int zoom_mesh_can_use_mesh(const struct zoom_pm_mesh *mesh,
                           const struct space *s, const struct cell *ci,
                           const struct cell *cj) {

  /* A disabled mesh cannot replace any tree interaction. */
  if (mesh == NULL) return 0;
  if (!mesh->enabled) return 0;

  /* Both cells must be fully covered by the high-resolution mesh. */
  if (!zoom_mesh_cells_are_covered(mesh, ci, cj, /*use_max_dx=*/0)) return 0;

  /* The gravity task recursion must only compare equal-sized cell pairs here. */
  if (ci->width[0] != cj->width[0] || ci->width[1] != cj->width[1] ||
      ci->width[2] != cj->width[2])
    error("Cells of different size in zoom mesh pair predicate!");

  /* Check whether the cells are safely beyond the zoom mesh cutoff. */
  const double max_distance2 = mesh->r_cut_max * mesh->r_cut_max;
  const double min_radius2 = cell_min_dist2(ci, cj, /*periodic=*/0, s->dim);

  return zoom_mesh_distance_is_above_cutoff(min_radius2, max_distance2);
}

/**
 * @brief Can a pair of cells still use the zoom mesh between rebuilds?
 *
 * This mirrors #zoom_mesh_can_use_mesh but includes the maximal particle
 * displacement since the last rebuild in the distance estimate.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @return 1 if the zoom mesh can still be used, 0 otherwise.
 */
int zoom_mesh_can_use_mesh_between_rebuilds(const struct zoom_pm_mesh *mesh,
                                            const struct space *s,
                                            const struct cell *ci,
                                            const struct cell *cj) {

  /* A disabled mesh cannot replace any tree interaction. */
  if (mesh == NULL) return 0;
  if (!mesh->enabled) return 0;

  /* Include particle motion since rebuild when checking mesh coverage. */
  if (!zoom_mesh_cells_are_covered(mesh, ci, cj, /*use_max_dx=*/1)) return 0;

  /* The gravity task recursion must only compare equal-sized cell pairs here. */
  if (ci->width[0] != cj->width[0] || ci->width[1] != cj->width[1] ||
      ci->width[2] != cj->width[2])
    error("Cells of different size in zoom mesh pair predicate!");

  /* Check whether the cells remain safely beyond the cutoff after drift. */
  const double max_distance2 = mesh->r_cut_max * mesh->r_cut_max;
  const double min_radius2 =
      cell_min_dist2_with_max_dx(ci, cj, /*periodic=*/0, s->dim);

  return zoom_mesh_distance_is_above_cutoff(min_radius2, max_distance2);
}
