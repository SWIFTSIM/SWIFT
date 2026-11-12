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
#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
#include <fftw3-mpi.h>
#endif
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
#include "mesh_gravity_mpi.h"
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
static void zoom_mesh_free_kernel(struct zoom_pm_mesh *mesh);
static int zoom_mesh_cell_is_inside(const struct zoom_pm_mesh *mesh,
                                    const struct cell *c, const int use_max_dx);

/**
 * @brief Initialize a zoom mesh structure for the case where it is disabled.
 *
 * @param mesh The #zoom_pm_mesh to initialise.
 */
void zoom_mesh_init_no_mesh(struct zoom_pm_mesh *mesh) {

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
  }
  if (mesh->dim[0] != mesh->dim[1] || mesh->dim[0] != mesh->dim[2])
    error("Gravity:zoom_mesh requires a cubic zoom mesh domain.");
  mesh->cell_fac = mesh->N / mesh->dim[0];

  for (int i = 0; i < 3; ++i) {

    /* The region in which forces can be interpolated. The five-point force
     * stencil needs two cells on the low side for i-2 and three on the high
     * side because the i+2 stencil point is CIC-interpolated and therefore
     * also reads i+3. A small tolerance allows for roundoff in the cell
     * construction. */
    const double cell_width = 1. / mesh->cell_fac;
    const double tol = 1e-6 * mesh->dim[i];
    mesh->active_min[i] = mesh->loc[i] + 2. * cell_width - tol;
    mesh->active_max[i] = mesh->loc[i] + mesh->dim[i] - 3. * cell_width + tol;
  }

  /* The smoothing scale must track the zoom mesh's own grid resolution
   * (mirroring how the background mesh derives r_s from dim/mesh_size in
   * gravity_props_init), NOT the zoom region's top-level tree cell width.
   * The tree cdim is set independently for domain-decomposition/load
   * balancing and can be far coarser than the requested mesh resolution. */
  mesh->r_s = props->a_smooth * mesh->dim[0] / mesh->N;
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
        mesh->r_s, mesh->N, mesh->dim[0], props->r_s);

  /* The search ranges used when constructing the gravity tasks are set by the
   * cut-off of the background mesh, so the zoom mesh cannot truncate the tree
   * forces further out than the background mesh does. */
  if (mesh->r_cut_max > props->r_s * props->r_cut_max_ratio)
    error(
        "The zoom gravity mesh truncates the tree forces at a larger distance "
        "(%.6e) than the background mesh (%.6e). Decrease "
        "Gravity:zoom_r_cut_max or increase Gravity:zoom_mesh_side_length.",
        mesh->r_cut_max, props->r_s * props->r_cut_max_ratio);

  /* Positions are not wrapped on the zoom mesh, so it has to lie inside the
   * box and no pair of covered particles can be closer through the periodic
   * boundary than the cut-off of the background mesh. */
  for (int i = 0; i < 3; ++i) {
    if (mesh->loc[i] < 0. || mesh->loc[i] + mesh->dim[i] > s->dim[i] ||
        mesh->dim[i] > s->dim[i] - props->r_s * props->r_cut_max_ratio)
      error(
          "The zoom gravity mesh ([%.6e, %.6e] along axis %d) must lie inside "
          "the box (size %.6e) and be smaller than the box by at least the "
          "cut-off of the background mesh (%.6e). Decrease "
          "Gravity:zoom_mesh_bkg_buffer_cells or shrink the zoom region.",
          mesh->loc[i], mesh->loc[i] + mesh->dim[i], i, s->dim[i],
          props->r_s * props->r_cut_max_ratio);
  }
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
  mesh->N = side_length;
  zoom_mesh_set_geometry(mesh, s, props, r_cut_max_ratio, r_cut_min_ratio);

  /* Report the mesh geometry to the user if requested. */
  if (verbose)
    message(
        "Zoom gravity mesh: enabled with N=%d, loc=[%.6e %.6e %.6e], "
        "dim=[%.6e %.6e %.6e], buffer_cells=%d, r_cut_max=%.6e.",
        mesh->N, mesh->loc[0], mesh->loc[1], mesh->loc[2], mesh->dim[0],
        mesh->dim[1], mesh->dim[2], buffer_cells, mesh->r_cut_max);

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
  if (mesh == NULL || !mesh->enabled) return;

  /* Release any allocated FFT storage and reset to a no-mesh state. */
  zoom_mesh_free_kernel(mesh);
  zoom_mesh_free(mesh);
  zoom_mesh_init_no_mesh(mesh);
}

/**
 * @brief Rebuild the zoom mesh after rebuilding the zoom top-level cells.
 *
 * The Fourier kernel depends on the mesh spacing, so invalidate it before
 * updating the mesh position and smoothing scale. The side-length of the mesh
 * does not change so the potential grid is kept.
 *
 * @param mesh The #zoom_pm_mesh to rebuild.
 * @param s The #space containing the new zoom region geometry.
 */
void zoom_mesh_regrid(struct zoom_pm_mesh *mesh, const struct space *s) {

  if (mesh == NULL || !mesh->enabled) return;

  /* Keep the configured cut-offs in units of the smoothing scale. */
  const double r_cut_max_ratio = mesh->r_cut_max / mesh->r_s;
  const double r_cut_min_ratio = mesh->r_cut_min / mesh->r_s;

  zoom_mesh_free_kernel(mesh);
  zoom_mesh_set_geometry(mesh, s, s->e->gravity_properties, r_cut_max_ratio,
                         r_cut_min_ratio);
}

/**
 * @brief Allocates the zero-padded potential grid to be ready for an FFT
 * calculation
 *
 * @param mesh The #zoom_pm_mesh structure.
 */
void zoom_mesh_allocate(struct zoom_pm_mesh *mesh) {

#ifdef HAVE_FFTW

  if (mesh->distributed_mesh) {

  } else {

    /* The active mesh is embedded in a zero-padded mesh twice as large */
    const size_t N = 2 * mesh->N;

    /* Allocate the memory for the combined density and potential array */
    mesh->potential_global = (double *)fftw_malloc(sizeof(double) * N * N * N);
    if (mesh->potential_global == NULL)
      error("Error allocating memory for the zoom gravity mesh.");
    memuse_log_allocation("fftw_zoom_mesh.potential", mesh->potential_global, 1,
                          sizeof(double) * N * N * N);
  }
#else
  error("No FFTW library found. Cannot compute zoom mesh gravity.");
#endif
}

/**
 * @brief Frees the potential grid.
 *
 * @param mesh The #zoom_pm_mesh structure.
 */
void zoom_mesh_free(struct zoom_pm_mesh *mesh) {

#ifdef HAVE_FFTW

  if (!mesh->distributed_mesh && mesh->potential_global) {
    memuse_log_allocation("fftw_zoom_mesh.potential", mesh->potential_global, 0,
                          0);
    fftw_free(mesh->potential_global);
    mesh->potential_global = NULL;
  }

#else
  error("No FFTW library found. Cannot compute zoom mesh gravity.");
#endif
}

/**
 * @brief Frees the cached Fourier transform of the kernel.
 *
 * @param mesh The #zoom_pm_mesh structure.
 */
static void zoom_mesh_free_kernel(struct zoom_pm_mesh *mesh) {

#ifdef HAVE_FFTW

  if (mesh->kernel_fft) {
    memuse_log_allocation("fftw_zoom_mesh.kernel_fft", mesh->kernel_fft, 0, 0);
    fftw_free(mesh->kernel_fft);
    mesh->kernel_fft = NULL;
  }
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

    initialise_fftw(2 * mesh->N, mesh->nr_threads);
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
__attribute__((always_inline, const)) INLINE static double zoom_mesh_CIC_get(
    double mesh[6][6][6], const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz) {

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
__attribute__((always_inline)) INLINE static void zoom_mesh_CIC_set(
    double *mesh, const int N, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz, const double value) {

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

  /*! Newton's constant in internal units. */
  float const_G;

  /*! Constants used for neutrino delta-f mesh weighting. */
  const struct neutrino_model *nu_model;
};

/**
 * @brief Assigns a given #gpart to the zoom density mesh using the CIC method.
 *
 * @param gp The #gpart.
 * @param rho The density mesh.
 * @param mesh The #zoom_pm_mesh.
 * @param nu_model Struct with neutrino constants
 */
INLINE static void zoom_mesh_gpart_to_mesh_CIC(
    const struct gpart *gp, double *rho, const struct zoom_pm_mesh *mesh,
    const struct neutrino_model *nu_model) {

  const int N = mesh->N;
  const double fac = mesh->cell_fac;

  /* Position relative to the corner of the zoom mesh */
  const double pos_x = gp->x[0] - mesh->loc[0];
  const double pos_y = gp->x[1] - mesh->loc[1];
  const double pos_z = gp->x[2] - mesh->loc[2];

  /* Workout the CIC coefficients */
  int i = (int)floor(fac * pos_x);
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int)floor(fac * pos_y);
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;

  int k = (int)floor(fac * pos_z);
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->time_bin == time_bin_not_created)
    error("Found an extra particle in mesh CIC.");

  if (i < 0 || i >= N) error("Invalid gpart position in x");
  if (j < 0 || j >= N) error("Invalid gpart position in y");
  if (k < 0 || k >= N) error("Invalid gpart position in z");
#endif

  /* Compute weight (for neutrino delta-f weighting) */
  double weight = 1.0;
  if (gp->type == swift_type_neutrino)
    gpart_neutrino_weight_mesh_only(gp, nu_model, &weight);

  const double mass = gp->mass;
  const double value = mass * weight;

  /* CIC ! (on the zero-padded mesh) */
  zoom_mesh_CIC_set(rho, 2 * N, i, j, k, tx, ty, tz, dx, dy, dz, value);
}

/**
 * @brief Assigns all the #gpart of a #cell to the zoom density mesh using the
 * CIC method.
 *
 * @param c The #cell.
 * @param rho The density mesh.
 * @param mesh The #zoom_pm_mesh.
 * @param nu_model Struct with neutrino constants
 */
static void zoom_mesh_cell_gpart_to_mesh_CIC(
    const struct cell *c, double *rho, const struct zoom_pm_mesh *mesh,
    const struct neutrino_model *nu_model) {

  const int gcount = c->grav.count;
  const struct gpart *gparts = c->grav.parts;

  /* Assign all the gpart of that cell to the mesh */
  for (int i = 0; i < gcount; ++i) {
    if (gparts[i].time_bin == time_bin_inhibited) continue;
    zoom_mesh_gpart_to_mesh_CIC(&gparts[i], rho, mesh, nu_model);
  }
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
      zoom_mesh_accumulate_cell_to_local_patch(data->mesh, c, &patch,
                                               data->nu_model);

      /* Copy the local patch values back onto the zoom mesh */
      pm_add_patch_to_global_mesh(data->rho, &patch);

      /* Free the allocated memory */
      pm_mesh_patch_clean(&patch);

    } else {

      /* Assign this cell's content directly atomically to the mesh */
      zoom_mesh_cell_gpart_to_mesh_CIC(c, data->rho, data->mesh,
                                       data->nu_model);
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
 * @brief Adds the zoom mesh correction to a gpart from a given mesh using the
 * CIC method.
 *
 * Zoom mesh equivalent of mesh_to_gpart_CIC(). The correction is added to the
 * forces and potential of the periodic mesh.
 *
 * @param gp The #gpart.
 * @param pot The potential mesh.
 * @param mesh The #zoom_pm_mesh.
 * @param const_G Gravitional constant
 */
static void zoom_mesh_to_gpart_CIC(struct gpart *gp, const double *pot,
                                   const struct zoom_pm_mesh *mesh,
                                   const float const_G) {

  /* The potential is stored on the zero-padded mesh */
  const int N = 2 * mesh->N;
  const double fac = mesh->cell_fac;

  /* Position relative to the corner of the zoom mesh */
  const double pos_x = gp->x[0] - mesh->loc[0];
  const double pos_y = gp->x[1] - mesh->loc[1];
  const double pos_z = gp->x[2] - mesh->loc[2];

  int i = (int)floor(fac * pos_x);
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int)floor(fac * pos_y);
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;

  int k = (int)floor(fac * pos_z);
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->time_bin == time_bin_not_created)
    error("Found an extra particle when computing gravity from mesh.");

  if (i < 0 || i >= mesh->N) error("Invalid gpart position in x");
  if (j < 0 || j >= mesh->N) error("Invalid gpart position in y");
  if (k < 0 || k >= mesh->N) error("Invalid gpart position in z");
#endif

  /* First, copy the necessary part of the mesh for stencil operations */
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
  p += zoom_mesh_CIC_get(phi, ii, jj, kk, tx, ty, tz, dx, dy, dz);

  /* ---- */

  /* 5-point stencil along each axis for the accelerations */
  a[0] += (1. / 12.) *
          zoom_mesh_CIC_get(phi, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (2. / 3.) *
          zoom_mesh_CIC_get(phi, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] += (2. / 3.) *
          zoom_mesh_CIC_get(phi, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (1. / 12.) *
          zoom_mesh_CIC_get(phi, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);

  a[1] += (1. / 12.) *
          zoom_mesh_CIC_get(phi, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (2. / 3.) *
          zoom_mesh_CIC_get(phi, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] += (2. / 3.) *
          zoom_mesh_CIC_get(phi, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (1. / 12.) *
          zoom_mesh_CIC_get(phi, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);

  a[2] += (1. / 12.) *
          zoom_mesh_CIC_get(phi, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  a[2] -= (2. / 3.) *
          zoom_mesh_CIC_get(phi, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  a[2] += (2. / 3.) *
          zoom_mesh_CIC_get(phi, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  a[2] -= (1. / 12.) *
          zoom_mesh_CIC_get(phi, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);

  /* ---- */

  /* Add the correction to the existing mesh contribution */
  gp->a_grav_mesh[0] += const_G * fac * a[0];
  gp->a_grav_mesh[1] += const_G * fac * a[1];
  gp->a_grav_mesh[2] += const_G * fac * a[2];
  gravity_add_comoving_mesh_potential(gp, const_G * p);
}

/**
 * @brief Adds the zoom mesh correction to all the #gpart of a #cell.
 *
 * @param c The #cell.
 * @param potential The potential mesh.
 * @param mesh The #zoom_pm_mesh.
 * @param const_G Gravitional constant
 */
static void zoom_mesh_cell_mesh_to_gpart_CIC(const struct cell *c,
                                             const double *potential,
                                             const struct zoom_pm_mesh *mesh,
                                             const float const_G) {

  const int gcount = c->grav.count;
  struct gpart *gparts = c->grav.parts;

  /* Interpolate the mesh onto all the gpart of that cell */
  for (int i = 0; i < gcount; ++i) {

    struct gpart *gp = &gparts[i];

    if (gp->time_bin == time_bin_inhibited) continue;

    zoom_mesh_to_gpart_CIC(gp, potential, mesh, const_G);
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
    zoom_mesh_cell_mesh_to_gpart_CIC(c, data->potential, data->mesh,
                                     data->const_G);
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
static void zoom_mesh_deposit_cells(const struct space *s,
                                    struct threadpool *tp,
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

/**
 * @brief Shared information about the zoom mesh kernel to be used by all the
 * threads in the pool.
 */
struct zoom_kernel_data {

  /*! The zoom mesh geometry and smoothing scale. */
  const struct zoom_pm_mesh *mesh;

  /*! The real-space kernel. */
  double *kernel;

  /*! The Fourier-space kernel. */
  fftw_complex *fkernel;

  /*! Inverse of the smoothing scale of the periodic mesh. */
  double global_r_s_inv;

  /*! Padding of the last dimension of the real-space array. */
  int pad;

  /*! The x coordinate of the start of the slice on this MPI rank. */
  int slice_offset;
};

/**
 * @brief Mapper function for the construction of the real-space kernel.
 *
 * The kernel is the long-range potential of the zoom mesh split minus the
 * long-range potential already provided by the periodic mesh. Indices in the
 * upper half of the zero-padded domain are negative offsets.
 *
 * @param map_data The array of the real-space kernel.
 * @param num The number of elements to iterate on (along the x-axis).
 * @param extra The properties of the kernel.
 */
static void zoom_mesh_kernel_mapper(void *map_data, const int num,
                                    void *extra) {

  const struct zoom_kernel_data *data = (const struct zoom_kernel_data *)extra;

  /* Unpack the array */
  double *const kernel = data->kernel;
  const struct zoom_pm_mesh *mesh = data->mesh;
  const int N = mesh->N;
  const int M = 2 * N;
  const size_t M_pad = M + data->pad;

  /* Unpack the split scales of the zoom and periodic meshes */
  const double r_s_inv = mesh->r_s_inv;
  const double global_r_s_inv = data->global_r_s_inv;

  /* Find what slice of the full mesh is stored on this MPI rank */
  const int slice_offset = data->slice_offset;

  /* Range of x coordinates in the full mesh handled by this call */
  const int i_start =
      (int)(((double *)map_data - kernel) / ((size_t)M * M_pad)) + slice_offset;
  const int i_end = i_start + num;

  /* Loop over the x range corresponding to this thread */
  for (int i = i_start; i < i_end; ++i) {

    const int ii = (i < N ? i : i - M);
    const double dx = ii / mesh->cell_fac;

    for (int j = 0; j < M; ++j) {

      const int jj = (j < N ? j : j - M);
      const double dy = jj / mesh->cell_fac;

      for (int k = 0; k < M; ++k) {

        const int kk = (k < N ? k : k - M);
        const double dz = kk / mesh->cell_fac;

        const double r = sqrt(dx * dx + dy * dy + dz * dz);

        /* Long-range potential of the zoom mesh split */
        const double zoom_long = (r == 0.) ? -M_2_SQRTPI * 0.5 * r_s_inv
                                           : -erf(0.5 * r * r_s_inv) / r;

        /* Long-range potential already provided by the periodic mesh */
        const double global_long = (r == 0.)
                                       ? -M_2_SQRTPI * 0.5 * global_r_s_inv
                                       : -erf(0.5 * r * global_r_s_inv) / r;

        /* Store the correction to the periodic mesh */
        const size_t index =
            ((size_t)(i - slice_offset) * M + j) * M_pad + (size_t)k;
        kernel[index] = zoom_long - global_long;
      }
    }
  }
}

/**
 * @brief Construct the real-space kernel of the zoom mesh correction.
 *
 * @param tp The threadpool.
 * @param mesh The #zoom_pm_mesh.
 * @param kernel The MxMx(M+pad) real array of the kernel.
 * @param global_r_s_inv Inverse of the smoothing scale of the periodic mesh.
 * @param slice_offset The x coordinate of the start of the slice on this MPI
 * rank
 * @param slice_width The width of the local slice on this MPI rank
 * @param pad Padding of the last dimension of the array.
 */
static void zoom_mesh_kernel(struct threadpool *tp,
                             const struct zoom_pm_mesh *mesh, double *kernel,
                             const double global_r_s_inv,
                             const int slice_offset, const int slice_width,
                             const int pad) {

  const int M = 2 * mesh->N;

  struct zoom_kernel_data data;
  data.mesh = mesh;
  data.kernel = kernel;
  data.fkernel = NULL;
  data.global_r_s_inv = global_r_s_inv;
  data.pad = pad;
  data.slice_offset = slice_offset;

  /* Parallelize the kernel construction using the threadpool
     to split the x-axis loop over the threads. */
  threadpool_map(tp, zoom_mesh_kernel_mapper, kernel, slice_width,
                 M * (M + pad) * sizeof(double), threadpool_auto_chunk_size,
                 &data);
}

/**
 * @brief Mapper function for the CIC deconvolution of the kernel.
 *
 * @param map_data The array of the kernel Fourier transform.
 * @param num The number of elements to iterate on (along the x-axis).
 * @param extra The properties of the kernel.
 */
static void zoom_mesh_deconvolve_kernel_mapper(void *map_data, const int num,
                                               void *extra) {

  const struct zoom_kernel_data *data = (const struct zoom_kernel_data *)extra;

  /* Unpack the array */
  fftw_complex *const fkernel = data->fkernel;
  const int M = 2 * data->mesh->N;
  const int M_half = M / 2;

  /* Find what slice of the full mesh is stored on this MPI rank */
  const int slice_offset = data->slice_offset;

  /* Range of x coordinates in the full mesh handled by this call */
  const int i_start = ((fftw_complex *)map_data - fkernel) + slice_offset;
  const int i_end = i_start + num;

  /* Loop over the x range corresponding to this thread */
  for (int i = i_start; i < i_end; ++i) {

    /* kx component of vector in Fourier space */
    const int kx = (i > M_half ? i - M : i);

    for (int j = 0; j < M; ++j) {

      /* ky component of vector in Fourier space */
      const int ky = (j > M_half ? j - M : j);

      for (int k = 0; k < M_half + 1; ++k) {

        /* kz component of vector in Fourier space */
        const int kz = (k > M_half ? k - M : k);

        /* Deconvolution of CIC */
        const double CIC_cor4 = zoom_mesh_CIC_deconvolution(kx, ky, kz, M);

        /* Apply to the mesh */
        const size_t index = (size_t)M * (M_half + 1) * (i - slice_offset) +
                             (size_t)(M_half + 1) * j + k;
        fkernel[index][0] *= CIC_cor4;
        fkernel[index][1] *= CIC_cor4;
      }
    }
  }
}

/**
 * @brief Deconvolve the CIC kernel from the Fourier transform of the zoom
 * mesh kernel.
 *
 * @param tp The threadpool.
 * @param mesh The #zoom_pm_mesh.
 * @param fkernel The MxMx(M/2+1) complex array of the Fourier transform of the
 * kernel.
 * @param slice_offset The x coordinate of the start of the slice on this MPI
 * rank
 * @param slice_width The width of the local slice on this MPI rank
 */
static void zoom_mesh_deconvolve_kernel(struct threadpool *tp,
                                        const struct zoom_pm_mesh *mesh,
                                        fftw_complex *fkernel,
                                        const int slice_offset,
                                        const int slice_width) {

  struct zoom_kernel_data data;
  data.mesh = mesh;
  data.kernel = NULL;
  data.fkernel = fkernel;
  data.global_r_s_inv = 0.;
  data.pad = 0;
  data.slice_offset = slice_offset;

  /* Parallelize the deconvolution using the threadpool
     to split the x-axis loop over the threads.
     The array is M x M x (M/2). We use the thread to each deal with
     a range [i_min, i_max[ x M x (M/2) */
  threadpool_map(tp, zoom_mesh_deconvolve_kernel_mapper, fkernel, slice_width,
                 sizeof(fftw_complex), threadpool_auto_chunk_size, &data);
}

/**
 * @brief Shared information about the kernel multiplication to be used by all
 * the threads in the pool.
 */
struct zoom_multiply_mapper_data {

  /*! The Fourier transform of the density field. */
  fftw_complex *density;

  /*! The Fourier transform of the kernel. */
  const fftw_complex *kernel;

  /*! The normalisation of the inverse transform. */
  double norm;
};

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

  /* Unpack sizes and memory requirements for the zero-padded FFT domain. */
  const int N = mesh->N;
  const int M = 2 * N;
  const int M_half = M / 2;
  const size_t real_size = (size_t)M * M * M;
  const size_t complex_size = (size_t)M * M * (M_half + 1);
  double *rho = mesh->potential_global;

  ticks tic = getticks();
  bzero(rho, real_size * sizeof(double));

  struct zoom_cic_mapper_data data;
  data.mesh = mesh;
  data.cells = s->cells_top;
  data.rho = rho;
  data.potential = NULL;
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
  if (frho == NULL) error("Error allocating zoom mesh FFT work arrays.");

  /* Record the temporary allocations in the memory-use logger. */
  memuse_log_allocation("fftw_zoom_mesh.frho", frho, 1,
                        complex_size * sizeof(fftw_complex));

  /* Prepare transforms of the density and resulting potential. */
  fftw_plan rho_forward = fftw_plan_dft_r2c_3d(
      M, M, M, rho, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  fftw_plan inverse = fftw_plan_dft_c2r_3d(M, M, M, frho, rho,
                                           FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

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

    zoom_mesh_kernel(tp, mesh, kernel, s->e->mesh->r_s_inv, /*slice_offset=*/0,
                     /*slice_width=*/M, /*pad=*/0);
    fftw_plan kernel_forward = fftw_plan_dft_r2c_3d(
        M, M, M, kernel, fkernel, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftw_execute(kernel_forward);
    fftw_destroy_plan(kernel_forward);
    zoom_mesh_deconvolve_kernel(tp, mesh, fkernel, /*slice_offset=*/0,
                                /*slice_width=*/M);
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

  error("No FFTW library found. Cannot compute zoom mesh gravity.");
}
#endif

/**
 * @brief Compute the zoom mesh correction to the gravity forces.
 *
 * This version uses the MPI version of FFTW with the zero-padded mesh
 * distributed in slices over the MPI ranks, as done for the periodic mesh in
 * compute_potential_distributed().
 *
 * @param mesh The #zoom_pm_mesh.
 * @param s The #space containing the particles.
 * @param tp The #threadpool object used for parallelisation.
 * @param verbose Are we talkative?
 */
static void zoom_mesh_compute_potential_distributed(struct zoom_pm_mesh *mesh,
                                                    const struct space *s,
                                                    struct threadpool *tp,
                                                    const int verbose) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  const int nr_local_cells = s->nr_local_cells;

  /* Some useful constants. The FFTs are performed on the zero-padded mesh. */
  const int M = 2 * mesh->N;

  ticks tic = getticks();

  /* Create an array of mesh patches. One per local top-level cell. */
  struct pm_mesh_patch *local_patches = (struct pm_mesh_patch *)malloc(
      nr_local_cells * sizeof(struct pm_mesh_patch));
  if (local_patches == NULL)
    error("Could not allocate array of local mesh patches!");
  memset(local_patches, 0, nr_local_cells * sizeof(struct pm_mesh_patch));

  /* Calculate contributions to density field on this MPI rank */
  zoom_mpi_mesh_accumulate_gparts_to_local_patches(tp, mesh, s, local_patches);
  if (verbose)
    message("Zoom mesh: Accumulating mass to local patches took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Ask FFTW what slice of the density field we need to store on this task.
     Note that fftw_mpi_local_size_3d works in terms of the size of the complex
     output. The last dimension of the real input is padded to 2*(M/2+1). */
  ptrdiff_t local_n0, local_0_start;
  ptrdiff_t nalloc =
      fftw_mpi_local_size_3d((ptrdiff_t)M, (ptrdiff_t)M, (ptrdiff_t)(M / 2 + 1),
                             MPI_COMM_WORLD, &local_n0, &local_0_start);
  if (verbose)
    message("Zoom mesh: Local density field slice has thickness %d.",
            (int)local_n0);
  if (verbose)
    message("Zoom mesh: local patch size = %d, local mesh cells = %lld",
            nr_local_cells, (long long)(local_n0 * M * M));

  /* Allocate storage for mesh slices.
   *
   * Note: nalloc is the number of *complex* values.
   */
  double *rho_slice = (double *)fftw_malloc(2 * nalloc * sizeof(double));
  memset(rho_slice, 0, 2 * nalloc * sizeof(double));

  tic = getticks();

  /* Construct density field slices from contributions stored in the local
   * patches.
   * Note: This cleans up the local_patches entries. */
  mpi_mesh_local_patches_to_slices(M, (int)local_n0, local_patches,
                                   nr_local_cells, rho_slice, tp, verbose);
  if (verbose)
    message("Zoom mesh: Assembling mesh slices took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Build and transform the static correction kernel only once. */
  if (mesh->kernel_fft == NULL) {

    tic = getticks();

    double *kernel_slice = (double *)fftw_malloc(2 * nalloc * sizeof(double));
    fftw_complex *fkernel_slice =
        (fftw_complex *)fftw_malloc(nalloc * sizeof(fftw_complex));
    if (kernel_slice == NULL || fkernel_slice == NULL)
      error("Error allocating zoom mesh kernel arrays.");
    memset(kernel_slice, 0, 2 * nalloc * sizeof(double));
    memuse_log_allocation("fftw_zoom_mesh.kernel_fft", fkernel_slice, 1,
                          nalloc * sizeof(fftw_complex));

    /* The last dimension of the real input is padded to 2*(M/2+1). */
    zoom_mesh_kernel(tp, mesh, kernel_slice, s->e->mesh->r_s_inv,
                     (int)local_0_start, (int)local_n0, /*pad=*/2);

    fftw_plan mpi_kernel_plan = fftw_mpi_plan_dft_r2c_3d(
        M, M, M, kernel_slice, fkernel_slice, MPI_COMM_WORLD,
        FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT | FFTW_DESTROY_INPUT);
    fftw_execute(mpi_kernel_plan);
    fftw_destroy_plan(mpi_kernel_plan);

    zoom_mesh_deconvolve_kernel(tp, mesh, fkernel_slice, (int)local_0_start,
                                (int)local_n0);

    fftw_free(kernel_slice);
    mesh->kernel_fft = fkernel_slice;

    if (verbose)
      message("Constructing the zoom mesh kernel took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

  tic = getticks();

  /* Allocate storage for the slices of the FFT of the density mesh */
  fftw_complex *frho_slice =
      (fftw_complex *)fftw_malloc(nalloc * sizeof(fftw_complex));

  /* Carry out the MPI Fourier transform. We can save a bit of time
   * if we allow FFTW to transpose the first two dimensions of the output.
   *
   * Layout of the MPI FFTW input and output:
   *
   * Input mesh contains M*M*M reals, padded to M*M*(2*(M/2+1)).
   * Output Fourier transform is M*M*(M/2+1) complex values.
   *
   * The first two dimensions of the transform are transposed in
   * the output. Each MPI rank has slice of thickness local_n0
   * starting at local_0_start in the first dimension.
   */
  fftw_plan mpi_plan = fftw_mpi_plan_dft_r2c_3d(
      M, M, M, rho_slice, frho_slice, MPI_COMM_WORLD,
      FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT | FFTW_DESTROY_INPUT);
  fftw_execute(mpi_plan);
  fftw_destroy_plan(mpi_plan);
  if (verbose)
    message("Zoom mesh: MPI Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Multiply the local slice of the MPI mesh by the kernel. The inverse FFTW
   * transform is unnormalised, hence the explicit 1 / M^3 factor here. */
  struct zoom_multiply_mapper_data multiply_data = {
      .density = frho_slice,
      .kernel = (const fftw_complex *)mesh->kernel_fft,
      .norm = 1. / ((double)M * M * M)};
  threadpool_map(tp, zoom_mesh_multiply_mapper, frho_slice,
                 (size_t)local_n0 * M * (M / 2 + 1), sizeof(fftw_complex),
                 threadpool_auto_chunk_size, &multiply_data);
  if (verbose)
    message("Applying zoom mesh kernel took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Carry out the reverse MPI Fourier transform */
  fftw_plan mpi_inverse_plan = fftw_mpi_plan_dft_c2r_3d(
      M, M, M, frho_slice, rho_slice, MPI_COMM_WORLD,
      FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN | FFTW_DESTROY_INPUT);
  fftw_execute(mpi_inverse_plan);
  fftw_destroy_plan(mpi_inverse_plan);

  if (verbose)
    message("Zoom mesh: MPI Reverse Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* We can now free the Fourier-space data */
  fftw_free(frho_slice);

  tic = getticks();

  /* Fetch MPI mesh entries we need on this rank from other ranks */
  zoom_mpi_mesh_fetch_potential(mesh, s, local_0_start, local_n0, rho_slice,
                                local_patches, tp, verbose);

  if (verbose)
    message("Zoom mesh: Fetching local potential took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Free the local slice of the potential */
  fftw_free(rho_slice);

  tic = getticks();

  /* Compute accelerations and potentials for the gparts */
  zoom_mpi_mesh_update_gparts(local_patches, mesh, s, tp);

  /* Clean the local patches array */
  for (int i = 0; i < nr_local_cells; ++i)
    pm_mesh_patch_clean(&local_patches[i]);
  free(local_patches);

  if (verbose)
    message("Zoom mesh: Computing mesh accelerations took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("No FFTW MPI library available. Cannot compute distributed mesh.");
#endif
}

/**
 * @brief Add the isolated zoom mesh correction to covered particles.
 *
 * Check the softening of the cells that actually receive the correction, using
 * the maximal softening recorded in their multipoles, then select the
 * replicated or distributed FFT solve.
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
    const int *lists[2] = {z->local_zoom_cells_top, z->local_bkg_cells_top};
    const int counts[2] = {z->nr_local_zoom_cells, z->nr_local_bkg_cells};
    for (int type = 0; type < 2; ++type) {
      for (int i = 0; i < counts[type]; ++i) {
        const struct cell *c = &s->cells_top[lists[type][i]];
        if (!zoom_mesh_cell_is_covered(mesh, c, /*use_max_dx=*/0)) continue;
        if (c->grav.count == 0) continue;
        max_softening =
            max(max_softening, (double)c->grav.multipole->m_pole.max_softening);
      }
    }
  }
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &max_softening, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
#endif
  if (max_softening > mesh->r_s)
    error(
        "Zoom mesh smoothing scale r_s=%g is smaller than the softening "
        "radius h=%g of a covered gravity particle.",
        mesh->r_s, max_softening);

  if (mesh->distributed_mesh)
    zoom_mesh_compute_potential_distributed(mesh, s, tp, verbose);
  else
    zoom_mesh_compute_potential_global(mesh, s, tp, verbose);
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
    /* Inflate the cell by the maximal particle displacement if requested. */
    const double dx_max = use_max_dx ? c->grav.multipole->dx_max[i] : 0.;
    const double cell_min = c->loc[i] - dx_max;
    const double cell_max = c->loc[i] + c->width[i] + dx_max;

    /* Any excursion outside the mesh bounds disqualifies the cell. */
    if (cell_min < mesh->active_min[i]) return 0;
    if (cell_max > mesh->active_max[i]) return 0;
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

    const double cell_width = 1. / mesh->cell_fac;

    /* 1. The zoom mesh must resolve the top-level zoom cells. */
    if (cell_width > z->width[i])
      error(
          "Zoom mesh too coarse given the size of top-level zoom cells "
          "(zoom cell width= %.4f, zoom mesh cell width= %.4f along axis %d). "
          "Gravity:zoom_mesh_side_length should be at least %d (currently: "
          "%d).",
          z->width[i], cell_width, i, (int)ceil(mesh->dim[i] / z->width[i]),
          mesh->N);

    /* 2. The active region must contain the whole zoom region. */
    const double active_min = mesh->active_min[i];
    const double active_max = mesh->active_max[i];

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

  /* Check whether the cells remain safely beyond the cutoff after drift. */
  const double max_distance2 = mesh->r_cut_max * mesh->r_cut_max;
  const double min_radius2 =
      cell_min_dist2_with_max_dx(ci, cj, /*periodic=*/0, s->dim);

  return cell_mesh_distance_is_above_cutoff(min_radius2, max_distance2);
}
