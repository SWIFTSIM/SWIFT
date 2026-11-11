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
#include "zoom_mesh_gravity_mpi.h"

/* Local headers. */
#include "atomic.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "gravity_properties.h"
#include "inline.h"
#include "kernel_long_gravity.h"
#include "memuse.h"
#include "mesh_gravity.h"
#include "mesh_gravity_patch.h"
#include "neutrino.h"
#include "parser.h"
#include "part.h"
#include "restart.h"
#include "row_major_id.h"
#include "space.h"
#include "threadpool.h"

/* Standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/* Local function prototypes. */
static void zoom_mesh_allocate(struct zoom_pm_mesh *mesh);
static void zoom_mesh_free(struct zoom_pm_mesh *mesh);
static void zoom_mesh_init_no_mesh(struct zoom_pm_mesh *mesh);
static int zoom_mesh_cell_is_inside(const struct zoom_pm_mesh *mesh,
                                    const struct cell *c, const int use_max_dx);

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
 * @brief Set the zoom mesh geometry from the current zoom region.
 *
 * The mesh covers the zoom region plus a buffer of background top-level
 * cells. The side length and buffer size must already be set.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space.
 * @param props The properties of the gravity scheme.
 * @param r_cut_max_ratio The maximal truncation distance in units of r_s.
 * @param r_cut_min_ratio The minimal truncation distance in units of r_s.
 */
static void zoom_mesh_set_geometry(struct zoom_pm_mesh *mesh,
                                   const struct space *s,
                                   const struct gravity_props *props,
                                   const double r_cut_max_ratio,
                                   const double r_cut_min_ratio) {

  const struct zoom_region_properties *z = s->zoom_props;

  /* The physical domain actually covered by the active zoom mesh. */
  for (int i = 0; i < 3; ++i) {
    const double buffer = mesh->buffer_cells * z->bkg_cells_top[0].width[i];
    mesh->loc[i] = z->region_lower_bounds[i] - buffer;
    mesh->dim[i] = z->dim[i] + 2. * buffer;
    mesh->cell_fac[i] = mesh->N[i] / mesh->dim[i];
  }
  if (mesh->dim[0] != mesh->dim[1] || mesh->dim[0] != mesh->dim[2])
    error("Gravity:zoom_mesh requires a cubic zoom mesh domain.");

  /* The smoothing scale must track the zoom mesh's own grid resolution
   * (mirroring how the background mesh derives r_s from dim/mesh_size in
   * gravity_props_init), NOT the zoom region's top-level tree cell width.
   * The tree cdim is set independently for domain-decomposition/load
   * balancing and can be far coarser than the requested mesh resolution. */
  mesh->r_s = props->a_smooth * mesh->dim[0] / mesh->N[0];
  mesh->r_s_inv = 1. / mesh->r_s;
  mesh->r_cut_max = mesh->r_s * r_cut_max_ratio;
  mesh->r_cut_min = mesh->r_s * r_cut_min_ratio;

  /* The zoom mesh only makes sense if it is strictly higher resolution than
   * the background mesh it is correcting. If it isn't, the mesh-tree
   * splitting criterion becomes *more* conservative inside the zoom region
   * (not less), silently exploding the gravity task count instead of
   * reducing it. */
  if (mesh->r_s >= props->r_s)
    error(
        "The zoom gravity mesh must be higher resolution than the "
        "background mesh, but its smoothing scale (r_s=%.6e, from "
        "Gravity:zoom_mesh_side_length=%d over a domain of %.6e) is not "
        "smaller than the background mesh's (r_s=%.6e, from "
        "Gravity:mesh_side_length). Increase Gravity:zoom_mesh_side_length "
        "or shrink the zoom region.",
        mesh->r_s, mesh->N[0], mesh->dim[0], props->r_s);
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
 * @param nr_threads Number of threads used by FFTW plans.
 * @param verbose Are we talkative?
 */
void zoom_mesh_init(struct zoom_pm_mesh *mesh, struct swift_params *params,
                    const struct gravity_props *props, const struct space *s,
                    const int nr_threads, const int verbose) {

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

  /* Replicated transforms retain the global mesh's practical serial limit. */
  if (!props->distributed_mesh && 2 * side_length > 1290)
    error(
        "Gravity:zoom_mesh_side_length is too large for the current "
        "zero-padded zoom mesh implementation. The active side length N is "
        "%d, but the zero padded Hockney-Eastwood FFT uses side length 2N=%d. "
        "Use N <= 645 or enable Gravity:distributed_mesh.",
        side_length, 2 * side_length);

  const double r_cut_max_ratio =
      parser_get_param_double(params, "Gravity:zoom_r_cut_max");
  if (r_cut_max_ratio <= 0.) error("Gravity:zoom_r_cut_max must be positive.");

  const double r_cut_min_ratio =
      parser_get_opt_param_double(params, "Gravity:zoom_r_cut_min", 0.);
  if (r_cut_min_ratio < 0.) error("Gravity:zoom_r_cut_min must be >= 0.");
  if (r_cut_min_ratio > r_cut_max_ratio)
    error("Gravity:zoom_r_cut_min must be <= Gravity:zoom_r_cut_max.");

  /* The truncated tree force must fit inside the mesh domain, mirroring the
   * background mesh's own r_cut_max check in gravity_props_init. */
  if (2. * props->a_smooth * r_cut_max_ratio > side_length)
    error(
        "Zoom mesh too small given Gravity:zoom_r_cut_max. "
        "Gravity:zoom_mesh_side_length should be at least %d cells wide "
        "(currently: %d).",
        (int)(2. * props->a_smooth * r_cut_max_ratio) + 1, side_length);

  const int buffer_cells =
      parser_get_opt_param_int(params, "Gravity:zoom_mesh_bkg_buffer_cells", 1);
  if (buffer_cells < 1)
    error("Gravity:zoom_mesh_bkg_buffer_cells must be >= 1.");

  if (s->zoom_props->nr_bkg_cells <= 0 || s->zoom_props->bkg_cells_top == NULL)
    error("Gravity:zoom_mesh requires background top-level cells.");
  if (s->zoom_props->nr_zoom_cells <= 0 ||
      s->zoom_props->zoom_cells_top == NULL)
    error("Gravity:zoom_mesh requires zoom top-level cells.");

  /* Fill in the mesh structure. */
  mesh->enabled = 1;
  mesh->nr_threads = nr_threads;
  mesh->distributed_mesh = props->distributed_mesh;
  mesh->use_local_patches = props->mesh_uses_local_patches;
  mesh->buffer_cells = buffer_cells;
  mesh->potential_global = NULL;
  for (int i = 0; i < 3; ++i) mesh->N[i] = side_length;
  zoom_mesh_set_geometry(mesh, s, props, r_cut_max_ratio, r_cut_min_ratio);

  /* Report the mesh geometry to the user if requested. */
  if (verbose)
    message(
        "Zoom gravity mesh: enabled with N=[%d %d %d], loc=[%.6e %.6e %.6e], "
        "dim=[%.6e %.6e %.6e], buffer_cells=%d, r_cut_max=%.6e.",
        mesh->N[0], mesh->N[1], mesh->N[2], mesh->loc[0], mesh->loc[1],
        mesh->loc[2], mesh->dim[0], mesh->dim[1], mesh->dim[2], buffer_cells,
        mesh->r_cut_max);

  /* The FFTs are performed on the zero-padded mesh. */
  initialise_fftw(2 * side_length, mesh->nr_threads);

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
 * @brief Rebuild the zoom mesh after rebuilding the zoom top-level cells.
 *
 * The Fourier kernel depends on the mesh spacing, so invalidate it along with
 * the old potential before updating the mesh position and smoothing scale.
 *
 * @param mesh The #zoom_pm_mesh to rebuild.
 * @param s The #space containing the new zoom region geometry.
 */
void zoom_mesh_regrid(struct zoom_pm_mesh *mesh, const struct space *s) {

  if (mesh == NULL || !mesh->enabled) return;

  /* Keep the configured cut-offs in units of the smoothing scale. */
  const double r_cut_max_ratio = mesh->r_cut_max / mesh->r_s;
  const double r_cut_min_ratio = mesh->r_cut_min / mesh->r_s;

  zoom_mesh_free(mesh);
  zoom_mesh_set_geometry(mesh, s, s->e->gravity_properties, r_cut_max_ratio,
                         r_cut_min_ratio);
  zoom_mesh_allocate(mesh);
}

/**
 * @brief Allocates the zero-padded zoom mesh grid.
 *
 * @param mesh The #zoom_pm_mesh structure.
 */
static void zoom_mesh_allocate(struct zoom_pm_mesh *mesh) {

  /* Nothing to allocate for a disabled or absent mesh. */
  if (mesh == NULL || !mesh->enabled || mesh->distributed_mesh) return;

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
  if (mesh == NULL ||
      (mesh->potential_global == NULL && mesh->kernel_fft == NULL))
    return;

#ifdef HAVE_FFTW

  if (mesh->potential_global != NULL) {
    memuse_log_allocation("fftw_zoom_mesh.potential", mesh->potential_global, 0,
                          0);
    fftw_free(mesh->potential_global);
    mesh->potential_global = NULL;
  }

  if (mesh->kernel_fft != NULL) {
    memuse_log_allocation("fftw_zoom_mesh.kernel_fft", mesh->kernel_fft, 0, 0);
    fftw_free(mesh->kernel_fft);
    mesh->kernel_fft = NULL;
  }
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

  /* Transient FFT arrays are not serialized; recreate them on restore. */
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

  /* The FFT arrays are not part of the dump. */
  mesh->potential_global = NULL;
  mesh->kernel_fft = NULL;

  if (mesh->enabled) {

    initialise_fftw(2 * mesh->N[0], mesh->nr_threads);
    zoom_mesh_allocate(mesh);
  }
}

/**
 * @brief Return the Fourier-space correction for CIC assignment and sampling.
 *
 * Same correction as applied to the periodic mesh in
 * mesh_apply_Green_function_mapper().
 *
 * @param kx Fourier mode along x.
 * @param ky Fourier mode along y.
 * @param kz Fourier mode along z.
 * @param N Mesh side length.
 * @return The CIC deconvolution factor.
 */
double zoom_mesh_CIC_deconvolution(const int kx, const int ky, const int kz,
                                   const int N) {

  const double k_fac = M_PI / (double)N;
  const double fx = k_fac * (double)kx;
  const double fy = k_fac * (double)ky;
  const double fz = k_fac * (double)kz;
  const double sinc_kx_inv = (kx != 0) ? fx / sin(fx) : 1.;
  const double sinc_ky_inv = (ky != 0) ? fy / sin(fy) : 1.;
  const double sinc_kz_inv = (kz != 0) ? fz / (sin(fz) + FLT_MIN) : 1.;
  const double CIC_cor = sinc_kx_inv * sinc_ky_inv * sinc_kz_inv;
  const double CIC_cor2 = CIC_cor * CIC_cor;
  return CIC_cor2 * CIC_cor2;
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

  /*! The top-level cell array. */
  struct cell *cells;

  /*! The density mesh to write to. */
  double *rho;

  /*! The potential mesh to read from. */
  const double *potential;

  /*! The size of the zero-padded mesh along one axis. */
  int M;

  /*! Newton's constant in internal units. */
  float const_G;

  /*! Constants used for neutrino delta-f mesh weighting. */
  const struct neutrino_model *nu_model;
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

  /* Long-range gravity can interact at the top-cell level. Do not apply a
   * fine-mesh correction to a covered leaf if its parent top cell still uses
   * the coarse tree split. */
  return zoom_mesh_cell_is_inside(mesh, c->top, use_max_dx) &&
         zoom_mesh_cell_is_inside(mesh, c, use_max_dx);
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
 * @brief Assigns all the #gpart of a #cell to the zoom density mesh using the
 * CIC method.
 *
 * @param c The #cell.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_cell_gpart_to_mesh_CIC(
    const struct cell *c, const struct zoom_cic_mapper_data *data) {

  const int gcount = c->grav.count;
  const struct gpart *gparts = c->grav.parts;

  /* Assign all the gpart of that cell to the mesh */
  for (int n = 0; n < gcount; ++n) {

    const struct gpart *gp = &gparts[n];

    if (gp->time_bin == time_bin_inhibited) continue;

    /* Workout the CIC coefficients */
    int ind[3];
    double d[3], t[3];
    if (!zoom_mesh_gpart_index(data->mesh, gp, ind, d, t, /*need_stencil=*/0))
      error(
          "Covered cell contains a particle outside the zoom mesh "
          "(x=[%g %g %g], index=[%d %d %d]).",
          gp->x[0], gp->x[1], gp->x[2], ind[0], ind[1], ind[2]);

    /* Compute weight (for neutrino delta-f weighting) */
    double weight = 1.0;
    if (gp->type == swift_type_neutrino)
      gpart_neutrino_weight_mesh_only(gp, data->nu_model, &weight);

    const double mass = gp->mass;
    const double value = mass * weight;

    /* CIC ! */
    zoom_mesh_CIC_set(data->rho, data->M, ind[0], ind[1], ind[2], t[0], t[1],
                      t[2], d[0], d[1], d[2], value);
  }
}

/**
 * @brief Initialize a mesh patch to cover the particles of a cell on the zoom
 * mesh.
 *
 * Zoom mesh equivalent of pm_mesh_patch_init(). The patch coordinates are
 * indices into the zero-padded zoom mesh, which is never wrapped.
 *
 * @param patch A pointer to the mesh patch.
 * @param mesh The #zoom_pm_mesh.
 * @param c The #cell which the patch should cover.
 * @return 1 if a patch was allocated, 0 if the cell has no particle to assign.
 */
static int zoom_mesh_patch_init(struct pm_mesh_patch *patch,
                                const struct zoom_pm_mesh *mesh,
                                const struct cell *c) {

  const int gcount = c->grav.count;
  const struct gpart *gparts = c->grav.parts;

  /* The patch indexes the zero-padded mesh */
  patch->N = 2 * mesh->N[0];
  patch->fac = mesh->cell_fac[0];

  /* No wrapping on the isolated mesh */
  for (int i = 0; i < 3; i++) {
    patch->wrap_min[i] = mesh->loc[i];
    patch->wrap_max[i] = mesh->loc[i] + mesh->dim[i];
  }

  /* Find the extent of the particle distribution in the cell */
  int found = 0;
  for (int i = 0; i < 3; i++) {
    patch->mesh_min[i] = INT_MAX;
    patch->mesh_max[i] = INT_MIN;
  }
  for (int ipart = 0; ipart < gcount; ipart++) {

    const struct gpart *gp = &gparts[ipart];

    if (gp->time_bin == time_bin_inhibited) continue;

    int ind[3];
    double d[3], t[3];
    if (!zoom_mesh_gpart_index(mesh, gp, ind, d, t, /*need_stencil=*/0))
      error(
          "Covered cell contains a particle outside the zoom mesh "
          "(x=[%g %g %g], index=[%d %d %d]).",
          gp->x[0], gp->x[1], gp->x[2], ind[0], ind[1], ind[2]);

    for (int i = 0; i < 3; i++) {
      if (ind[i] < patch->mesh_min[i]) patch->mesh_min[i] = ind[i];
      /* CIC interpolation requires one extra element in the positive
       * direction */
      if (ind[i] + 1 > patch->mesh_max[i]) patch->mesh_max[i] = ind[i] + 1;
    }
    found = 1;
  }
  if (!found) return 0;

  /* Determine the integer size of the mesh */
  int num_cells = 1;
  for (int i = 0; i < 3; i++) {
    patch->mesh_size[i] = patch->mesh_max[i] - patch->mesh_min[i] + 1;
    num_cells *= patch->mesh_size[i];
  }

  /* Allocate the mesh */
  if (swift_memalign("mesh_patch", (void **)&patch->mesh, SWIFT_CACHE_ALIGNMENT,
                     num_cells * sizeof(double)) != 0)
    error("Failed to allocate array for mesh patch!");

  return 1;
}

/**
 * @brief Accumulate contributions from a cell to a local zoom mesh patch.
 *
 * Zoom mesh equivalent of accumulate_cell_to_local_patch().
 *
 * @param c The #cell containing the particles.
 * @param data The shared CIC mapper data.
 * @param patch The local mesh patch.
 * @return 1 if a patch was allocated, 0 if the cell has no particle to assign.
 */
static int zoom_mesh_accumulate_cell_to_local_patch(
    const struct cell *c, const struct zoom_cic_mapper_data *data,
    struct pm_mesh_patch *patch) {

  /* Initialise the local mesh patch */
  if (!zoom_mesh_patch_init(patch, data->mesh, c)) return 0;
  pm_mesh_patch_zero(patch);

  /* Loop over particles in this cell */
  for (int ipart = 0; ipart < c->grav.count; ipart++) {

    const struct gpart *gp = &(c->grav.parts[ipart]);

    if (gp->time_bin == time_bin_inhibited) continue;

    /* Workout the CIC coefficients */
    int ind[3];
    double d[3], t[3];
    zoom_mesh_gpart_index(data->mesh, gp, ind, d, t, /*need_stencil=*/0);

    /* Get coordinates within the mesh patch */
    const int ii = ind[0] - patch->mesh_min[0];
    const int jj = ind[1] - patch->mesh_min[1];
    const int kk = ind[2] - patch->mesh_min[2];

    /* Compute weight (for neutrino delta-f weighting) */
    double weight = 1.0;
    if (gp->type == swift_type_neutrino)
      gpart_neutrino_weight_mesh_only(gp, data->nu_model, &weight);

    /* Accumulate contributions to the local mesh patch */
    const double mass = gp->mass;
    const double value = mass * weight;
    pm_mesh_patch_CIC_set(patch, ii, jj, kk, t[0], t[1], t[2], d[0], d[1], d[2],
                          value);
  }

  return 1;
}

/**
 * @brief Threadpool mapper function for the zoom mesh CIC assignment of a
 * cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
static void zoom_mesh_deposit_cells_mapper(void *map_data, int num,
                                           void *extra) {

  /* Unpack the shared information */
  const struct zoom_cic_mapper_data *data =
      (const struct zoom_cic_mapper_data *)extra;
  const struct cell *cells = data->cells;

  /* Pointer to the chunk to be processed */
  int *local_cells = (int *)map_data;

  /* A temporary patch of the zoom mesh */
  struct pm_mesh_patch patch;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell *c = &cells[local_cells[i]];

    /* Skip empty cells */
    if (c->grav.count == 0) continue;

    /* Particles in cells the zoom mesh does not cover keep the global PM
     * split. This is the same predicate as used for the tree split. */
    if (!zoom_mesh_cell_is_covered(data->mesh, c, /*use_max_dx=*/0)) continue;

    if (data->mesh->use_local_patches) {

      /* Do a CIC interpolation of all the particles in this cell onto
         the local patch (allocates memory in the patch) */
      if (!zoom_mesh_accumulate_cell_to_local_patch(c, data, &patch)) continue;

      /* Copy the local patch values back onto the zoom mesh */
      pm_add_patch_to_global_mesh(data->rho, &patch);

      /* Free the allocated memory */
      pm_mesh_patch_clean(&patch);

    } else {

      /* Assign this cell's content directly atomically to the mesh */
      zoom_mesh_cell_gpart_to_mesh_CIC(c, data);
    }
  }
}

/**
 * @brief Label zoom-cell density assignment in threadpool diagnostics.
 *
 * @param map_data The zoom-cell indices.
 * @param num The number of indices.
 * @param extra The shared #zoom_cic_mapper_data.
 */
static void zoom_mesh_deposit_zoom_cells_mapper(void *map_data, int num,
                                                 void *extra) {
  zoom_mesh_deposit_cells_mapper(map_data, num, extra);
}

/**
 * @brief Label background-cell density assignment in threadpool diagnostics.
 *
 * @param map_data The background-cell indices.
 * @param num The number of indices.
 * @param extra The shared #zoom_cic_mapper_data.
 */
static void zoom_mesh_deposit_bkg_cells_mapper(void *map_data, int num,
                                                void *extra) {
  zoom_mesh_deposit_cells_mapper(map_data, num, extra);
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

    /* Skip particles outside the active mesh or too close to the stencil edge.
     */
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
    gravity_add_comoving_mesh_potential(gp, const_G * p);
  }
}

/**
 * @brief Threadpool mapper function for the zoom mesh CIC interpolation of a
 * cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
static void zoom_mesh_interpolate_cells_mapper(void *map_data, int num,
                                               void *extra) {

  /* Unpack the shared information */
  const struct zoom_cic_mapper_data *data =
      (const struct zoom_cic_mapper_data *)extra;
  const struct cell *cells = data->cells;

  /* Pointer to the chunk to be processed */
  int *local_cells = (int *)map_data;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell *c = &cells[local_cells[i]];

    /* Only cells covered by the zoom mesh receive the correction */
    if (!zoom_mesh_cell_is_covered(data->mesh, c, /*use_max_dx=*/0)) continue;

    /* Interpolate the correction onto this cell's content */
    zoom_mesh_to_gpart_CIC_mapper(c->grav.parts, c->grav.count, (void *)data);
  }
}

/**
 * @brief Label zoom-cell force interpolation in threadpool diagnostics.
 *
 * @param map_data The zoom-cell indices.
 * @param num The number of indices.
 * @param extra The shared #zoom_cic_mapper_data.
 */
static void zoom_mesh_interpolate_zoom_cells_mapper(void *map_data, int num,
                                                     void *extra) {
  zoom_mesh_interpolate_cells_mapper(map_data, num, extra);
}

/**
 * @brief Label background-cell force interpolation in diagnostics.
 *
 * @param map_data The background-cell indices.
 * @param num The number of indices.
 * @param extra The shared #zoom_cic_mapper_data.
 */
static void zoom_mesh_interpolate_bkg_cells_mapper(void *map_data, int num,
                                                    void *extra) {
  zoom_mesh_interpolate_cells_mapper(map_data, num, extra);
}

/**
 * @brief Deposit all local covered cells on the zoom mesh.
 *
 * @param s The #space.
 * @param tp The #threadpool used to process the cells.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_deposit_cells(const struct space *s, struct threadpool *tp,
                                     const struct zoom_cic_mapper_data *data) {

  threadpool_map(tp, zoom_mesh_deposit_zoom_cells_mapper,
                 s->zoom_props->local_zoom_cells_top,
                 s->zoom_props->nr_local_zoom_cells, sizeof(int),
                 threadpool_auto_chunk_size, (void *)data);
  threadpool_map(tp, zoom_mesh_deposit_bkg_cells_mapper,
                 s->zoom_props->local_bkg_cells_top,
                 s->zoom_props->nr_local_bkg_cells, sizeof(int),
                 threadpool_auto_chunk_size, (void *)data);
}

/**
 * @brief Interpolate the zoom mesh correction to all local covered cells.
 *
 * @param s The #space.
 * @param tp The #threadpool used to process the cells.
 * @param data The shared CIC mapper data.
 */
static void zoom_mesh_interpolate_cells(
    const struct space *s, struct threadpool *tp,
    const struct zoom_cic_mapper_data *data) {

  threadpool_map(tp, zoom_mesh_interpolate_zoom_cells_mapper,
                 s->zoom_props->local_zoom_cells_top,
                 s->zoom_props->nr_local_zoom_cells, sizeof(int),
                 threadpool_auto_chunk_size, (void *)data);
  threadpool_map(tp, zoom_mesh_interpolate_bkg_cells_mapper,
                 s->zoom_props->local_bkg_cells_top,
                 s->zoom_props->nr_local_bkg_cells, sizeof(int),
                 threadpool_auto_chunk_size, (void *)data);
}

/** @brief Shared data for constructing the isolated kernel. */
struct zoom_kernel_mapper_data {
  const struct zoom_pm_mesh *mesh;
  double *kernel;
  double global_r_s_inv;
};

/**
 * @brief Threadpool mapper constructing slices of the isolated kernel.
 *
 * @param map_data The first real-space kernel slab in this chunk.
 * @param num The number of slabs in the chunk.
 * @param extra The mesh geometry and global smoothing scale.
 */
static void zoom_mesh_apply_kernel_mapper(void *map_data, int num,
                                          void *extra) {

  const struct zoom_kernel_mapper_data *data =
      (const struct zoom_kernel_mapper_data *)extra;
  const struct zoom_pm_mesh *mesh = data->mesh;
  double *kernel = (double *)map_data;

  /* Unpack the active and zero-padded mesh sizes. */
  const int N = mesh->N[0];
  const int M = 2 * N;

  /* Unpack the split scales for the zoom and global PM meshes. */
  const double r_s_inv = mesh->r_s_inv;
  const double global_r_s_inv = data->global_r_s_inv;
  const double global_r_s = global_r_s_inv > 0. ? 1. / global_r_s_inv : 0.;
  const int i_start = (int)((kernel - data->kernel) / ((size_t)M * M));

  /* Fill the embedded real-space correction Green's function. */
  for (int i = i_start; i < i_start + num; ++i) {
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

        /* Store only the correction that must be added to the global PM field.
         */
        data->kernel[zoom_mesh_index(i, j, k, M)] = zoom_long - global_long;
      }
    }
  }
}

/**
 * @brief Construct the zero padded Hockney-Eastwood isolated-convolution
 * kernel using the threadpool.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param kernel The zero-padded real-space kernel to fill.
 * @param s The #space.
 * @param tp The #threadpool object used for parallelisation.
 */
static void zoom_mesh_apply_kernel(struct zoom_pm_mesh *mesh, double *kernel,
                                   const struct space *s,
                                   struct threadpool *tp) {

  const int M = 2 * mesh->N[0];
  struct zoom_kernel_mapper_data data = {
      .mesh = mesh,
      .kernel = kernel,
      .global_r_s_inv = s->periodic ? s->e->mesh->r_s_inv : 0.};

  threadpool_map(tp, zoom_mesh_apply_kernel_mapper, kernel, M,
                 M * M * (int)sizeof(double), threadpool_auto_chunk_size,
                 &data);
}

struct zoom_multiply_mapper_data {
  fftw_complex *density;
  const fftw_complex *kernel;
  double norm;
};

struct zoom_deconvolution_mapper_data {
  fftw_complex *kernel;
  int M;
};

/**
 * @brief Apply the global PM CIC deconvolution to the cached zoom kernel.
 *
 * @param map_data The Fourier coefficients in this chunk.
 * @param num The number of coefficients.
 * @param extra The full kernel and its mesh side length.
 */
static void zoom_mesh_deconvolve_mapper(void *map_data, int num, void *extra) {

  fftw_complex *kernel = (fftw_complex *)map_data;
  const struct zoom_deconvolution_mapper_data *data =
      (const struct zoom_deconvolution_mapper_data *)extra;
  const size_t offset = kernel - data->kernel;
  const int M = data->M;
  const int complex_width = M / 2 + 1;

  for (int i = 0; i < num; ++i) {
    const size_t index = offset + i;
    const int kx_index = index / ((size_t)M * complex_width);
    const int remainder = index % ((size_t)M * complex_width);
    const int ky_index = remainder / complex_width;
    const int kz = remainder % complex_width;
    const int kx = kx_index > M / 2 ? kx_index - M : kx_index;
    const int ky = ky_index > M / 2 ? ky_index - M : ky_index;
    const double CIC_cor4 = zoom_mesh_CIC_deconvolution(kx, ky, kz, M);
    kernel[i][0] *= CIC_cor4;
    kernel[i][1] *= CIC_cor4;
  }
}

/**
 * @brief Threadpool mapper multiplying density by the Fourier-space kernel.
 *
 * @param map_data The density coefficients in this chunk.
 * @param num The number of coefficients.
 * @param extra The kernel and FFT normalisation.
 */
static void zoom_mesh_multiply_mapper(void *map_data, int num, void *extra) {

  fftw_complex *density = (fftw_complex *)map_data;
  const struct zoom_multiply_mapper_data *data =
      (const struct zoom_multiply_mapper_data *)extra;
  const size_t offset = density - data->density;

  for (int i = 0; i < num; ++i) {
    const size_t k = offset + i;
    const double a = density[i][0];
    const double b = density[i][1];
    const double c = data->kernel[k][0];
    const double d = data->kernel[k][1];
    density[i][0] = data->norm * (a * c - b * d);
    density[i][1] = data->norm * (a * d + b * c);
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
 * back to particles and added to their existing mesh acceleration and
 * potential.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
static void zoom_mesh_compute_potential_global(struct zoom_pm_mesh *mesh,
                                               const struct space *s,
                                               struct threadpool *tp,
                                               int verbose) {

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

  ticks tic = getticks();
  bzero(rho, real_size * sizeof(double));

  struct zoom_cic_mapper_data data;
  data.mesh = mesh;
  data.cells = s->cells_top;
  data.rho = rho;
  data.potential = NULL;
  data.M = M;
  data.const_G = 0.f;

  /* Match the global PM mesh's delta-f neutrino mass weighting. */
  struct neutrino_model nu_model;
  bzero(&nu_model, sizeof(struct neutrino_model));
  if (s->e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(s, &nu_model);
  data.nu_model = &nu_model;

  /* Assign the local covered cells to the active zoom mesh. This mirrors the
   * cell predicate used to select the zoom split in the tree. */
  zoom_mesh_deposit_cells(s, tp, &data);

  if (verbose)
    message("Zoom mesh gpart assignment took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#ifdef WITH_MPI

  MPI_Barrier(MPI_COMM_WORLD);
  tic = getticks();

  /* Merge everybody's share of the density mesh */
  MPI_Allreduce(MPI_IN_PLACE, rho, real_size, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  if (verbose)
    message("Zoom mesh MPI-reduction took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#endif

  tic = getticks();

  fftw_complex *frho =
      (fftw_complex *)fftw_malloc(complex_size * sizeof(fftw_complex));
  if (frho == NULL)
    error("Error allocating zoom mesh FFT work arrays.");

  /* Record the temporary allocations in the memory-use logger. */
  memuse_log_allocation("fftw_zoom_mesh.frho", frho, 1,
                        complex_size * sizeof(fftw_complex));

  /* Prepare transforms of the density and resulting potential. */
  fftw_plan rho_forward =
      fftw_plan_dft_r2c_3d(M, M, M, rho, frho,
                           FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  fftw_plan inverse = fftw_plan_dft_c2r_3d(
      M, M, M, frho, rho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  /* Build and transform the static correction kernel only once. */
  if (mesh->kernel_fft == NULL) {
    double *kernel = (double *)fftw_malloc(real_size * sizeof(double));
    fftw_complex *fkernel =
        (fftw_complex *)fftw_malloc(complex_size * sizeof(fftw_complex));
    if (kernel == NULL || fkernel == NULL)
      error("Error allocating zoom mesh kernel arrays.");
    memuse_log_allocation("fftw_zoom_mesh.kernel", kernel, 1,
                          real_size * sizeof(double));
    memuse_log_allocation("fftw_zoom_mesh.kernel_fft", fkernel, 1,
                          complex_size * sizeof(fftw_complex));

    zoom_mesh_apply_kernel(mesh, kernel, s, tp);
    fftw_plan kernel_forward =
        fftw_plan_dft_r2c_3d(M, M, M, kernel, fkernel,
                             FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftw_execute(kernel_forward);
    fftw_destroy_plan(kernel_forward);
    struct zoom_deconvolution_mapper_data deconvolution_data = {
        .kernel = fkernel, .M = M};
    threadpool_map(tp, zoom_mesh_deconvolve_mapper, fkernel, complex_size,
                   sizeof(fftw_complex), threadpool_auto_chunk_size,
                   &deconvolution_data);
    memuse_log_allocation("fftw_zoom_mesh.kernel", kernel, 0, 0);
    fftw_free(kernel);
    mesh->kernel_fft = fkernel;
  }

  /* Transform density to Fourier space. */
  fftw_execute(rho_forward);

  /* Multiply the Fourier-space density and kernel. The inverse FFTW transform
   * is unnormalised, hence the explicit 1 / M^3 factor here. */
  struct zoom_multiply_mapper_data multiply_data = {
      .density = frho,
      .kernel = (const fftw_complex *)mesh->kernel_fft,
      .norm = 1. / (double)real_size};
  threadpool_map(tp, zoom_mesh_multiply_mapper, frho, complex_size,
                 sizeof(fftw_complex), threadpool_auto_chunk_size,
                 &multiply_data);

  /* Transform back to real space; rho now stores the potential correction. */
  fftw_execute(inverse);

  /* Destroy the plans before releasing the temporary arrays. */
  fftw_destroy_plan(rho_forward);
  fftw_destroy_plan(inverse);

  /* Release all temporary FFT storage. */
  memuse_log_allocation("fftw_zoom_mesh.frho", frho, 0, 0);
  fftw_free(frho);

  if (verbose)
    message("Zoom mesh convolution took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Reuse the mapper data to interpolate the potential correction back to the
   * particles as accelerations and potentials. */
  data.rho = NULL;
  data.potential = rho;
  data.const_G = s->e->physical_constants->const_newton_G;

  zoom_mesh_interpolate_cells(s, tp, &data);

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
static void zoom_mesh_compute_potential_global(struct zoom_pm_mesh *mesh,
                                               const struct space *s,
                                               struct threadpool *tp,
                                               int verbose) {

  /* Disabled zoom meshes must remain harmless in non-FFTW builds. */
  if (mesh == NULL || !mesh->enabled) return;

  error("No FFTW library found. Cannot compute zoom mesh gravity.");
}
#endif

/**
 * @brief Add the isolated zoom mesh correction to covered particles.
 *
 * Check the softening of particles that actually receive the correction, then
 * select the replicated or distributed FFT solve.
 *
 * @param mesh The #zoom_pm_mesh, or NULL when disabled.
 * @param s The #space containing the particles.
 * @param tp The #threadpool used for deposition and interpolation.
 * @param verbose Are we talkative?
 */
void zoom_mesh_compute_potential(struct zoom_pm_mesh *mesh,
                                 const struct space *s, struct threadpool *tp,
                                 int verbose) {

  if (mesh == NULL || !mesh->enabled) return;

  /* Background particles can have mass-dependent softenings larger than the
   * configured high-resolution values. Check only cells receiving the zoom
   * correction; uncovered cells keep the global tree-PM split. */
  double max_softening = 0.;
  if (s->cells_top != NULL && s->zoom_props != NULL) {
    const struct zoom_region_properties *z = s->zoom_props;
    const int *lists[2] = {z->local_zoom_cells_top,
                           z->local_bkg_cells_top};
    const int counts[2] = {z->nr_local_zoom_cells,
                           z->nr_local_bkg_cells};
    for (int type = 0; type < 2; ++type) {
      for (int i = 0; i < counts[type]; ++i) {
        const struct cell *c = &s->cells_top[lists[type][i]];
        if (!zoom_mesh_cell_is_covered(mesh, c, /*use_max_dx=*/0)) continue;
        for (int j = 0; j < c->grav.count; ++j) {
          const struct gpart *gp = &c->grav.parts[j];
          if (gp->time_bin == time_bin_inhibited) continue;
          max_softening = max(max_softening,
                              (double)gravity_get_softening(
                                  gp, s->e->gravity_properties));
        }
      }
    }
  }
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &max_softening, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
#endif
  if (max_softening > mesh->r_s)
    error("Zoom mesh smoothing scale r_s=%g is smaller than the softening "
          "radius h=%g of a covered gravity particle.",
          mesh->r_s, max_softening);

  if (mesh->distributed_mesh)
    zoom_mesh_compute_potential_distributed(mesh, s, tp, verbose);
  else
    zoom_mesh_compute_potential_global(mesh, s, tp, verbose);
}

/**
 * @brief Bounds of the zoom mesh region in which forces can be interpolated.
 *
 * The five-point force stencil needs two cells on the low side for i-2 and
 * three on the high side because the i+2 stencil point is CIC-interpolated
 * and therefore also reads i+3. A small tolerance allows for roundoff in the
 * cell construction.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param i The axis.
 * @param active_min The lower bound in box coordinates (return).
 * @param active_max The upper bound in box coordinates (return).
 */
static void zoom_mesh_get_active_bounds(const struct zoom_pm_mesh *mesh,
                                        const int i, double *active_min,
                                        double *active_max) {

  const double cell_width = 1. / mesh->cell_fac[i];
  const double tol = 1e-6 * mesh->dim[i];
  *active_min = mesh->loc[i] + 2. * cell_width - tol;
  *active_max = mesh->loc[i] + mesh->dim[i] - 3. * cell_width + tol;
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
     * where the five-point force stencil is valid. */
    double mesh_min, mesh_max;
    zoom_mesh_get_active_bounds(mesh, i, &mesh_min, &mesh_max);

    /* Inflate the cell by the maximal particle displacement if requested. */
    const double dx_max = use_max_dx ? c->grav.multipole->dx_max[i] : 0.;
    const double cell_min = c->loc[i] - dx_max;
    const double cell_max = c->loc[i] + c->width[i] + dx_max;

    /* Any excursion outside the mesh bounds disqualifies the cell. */
    if (cell_min < mesh_min) return 0;
    if (cell_max > mesh_max) return 0;
  }

  return 1;
}

/**
 * @brief Check the zoom mesh geometry against the zoom region it must cover.
 *
 * When the zoom mesh is enabled it, rather than the background mesh, supplies
 * the long-range forces inside the zoom region. It therefore inherits the two
 * requirements the background mesh has to satisfy for the background cells
 * (see gravity_props_init):
 *
 *   1. It must resolve the top-level cells it is responsible for, i.e. its
 *      cells must be no wider than the zoom top-level cells.
 *   2. Its active region must actually contain those cells. Only cells fully
 *      inside the active region are covered (zoom_mesh_cell_is_covered); any
 *      zoom cell left outside would fall back on the background mesh, which is
 *      not guaranteed to resolve it.
 *
 * Both depend on the zoom region geometry, so this is called from space_regrid
 * and re-validated whenever the region is rebuilt.
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space.
 */
void zoom_mesh_check_geometry(const struct zoom_pm_mesh *mesh,
                              const struct space *s) {

  /* A disabled or absent mesh covers nothing and constrains nothing. */
  if (mesh == NULL || !mesh->enabled) return;

  const struct zoom_region_properties *z = s->zoom_props;

  for (int i = 0; i < 3; ++i) {

    const double cell_width = 1. / mesh->cell_fac[i];

    /* 1. The zoom mesh must resolve the top-level zoom cells. */
    if (cell_width > z->width[i])
      error(
          "Zoom mesh too coarse given the size of top-level zoom cells "
          "(zoom cell width= %.4f, zoom mesh cell width= %.4f along axis %d). "
          "Gravity:zoom_mesh_side_length should be at least %d (currently: "
          "%d).",
          z->width[i], cell_width, i, (int)ceil(mesh->dim[i] / z->width[i]),
          mesh->N[i]);

    /* 2. The active region must contain the whole zoom region. */
    double active_min, active_max;
    zoom_mesh_get_active_bounds(mesh, i, &active_min, &active_max);

    if (z->region_lower_bounds[i] < active_min ||
        z->region_upper_bounds[i] > active_max)
      error(
          "The active zoom mesh region ([%.4f, %.4f] along axis %d) does not "
          "contain the zoom region ([%.4f, %.4f]), so some zoom cells would "
          "not be covered by the zoom mesh. Increase "
          "Gravity:zoom_mesh_side_length or "
          "Gravity:zoom_mesh_bkg_buffer_cells.",
          active_min, active_max, i, z->region_lower_bounds[i],
          z->region_upper_bounds[i]);
  }
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

  /* The gravity task recursion must only compare equal-sized cell pairs here.
   */
  if (ci->width[0] != cj->width[0] || ci->width[1] != cj->width[1] ||
      ci->width[2] != cj->width[2])
    error("Cells of different size in zoom mesh pair predicate!");

  /* Check whether the cells are safely beyond the zoom mesh cutoff. */
  const double max_distance2 = mesh->r_cut_max * mesh->r_cut_max;
  const double min_radius2 = cell_min_dist2(ci, cj, /*periodic=*/0, s->dim);

  return cell_mesh_distance_is_above_cutoff(min_radius2, max_distance2);
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

  /* The gravity task recursion must only compare equal-sized cell pairs here.
   */
  if (ci->width[0] != cj->width[0] || ci->width[1] != cj->width[1] ||
      ci->width[2] != cj->width[2])
    error("Cells of different size in zoom mesh pair predicate!");

  /* Check whether the cells remain safely beyond the cutoff after drift. */
  const double max_distance2 = mesh->r_cut_max * mesh->r_cut_max;
  const double min_radius2 =
      cell_min_dist2_with_max_dx(ci, cj, /*periodic=*/0, s->dim);

  return cell_mesh_distance_is_above_cutoff(min_radius2, max_distance2);
}
