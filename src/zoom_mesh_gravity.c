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
  if (mesh->dim[0] != mesh->dim[1] || mesh->dim[0] != mesh->dim[2]) {
    error("Gravity:zoom_mesh requires a cubic zoom mesh domain.");
  }
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
  if (mesh->r_s >= props->r_s) {
    error(
        "The zoom gravity mesh must be higher resolution than the "
        "background mesh, but its smoothing scale (r_s=%.6e, from "
        "Gravity:zoom_mesh_side_length=%d over a domain of %.6e) is not "
        "smaller than the background mesh's (r_s=%.6e, from "
        "Gravity:mesh_side_length). Increase Gravity:zoom_mesh_side_length "
        "or shrink the zoom region.",
        mesh->r_s, mesh->N, mesh->dim[0], props->r_s);
  }

  /* The search ranges used when constructing the gravity tasks are set by the
   * cut-off of the background mesh, so the zoom mesh cannot truncate the tree
   * forces further out than the background mesh does. */
  if (mesh->r_cut_max > props->r_s * props->r_cut_max_ratio) {
    error(
        "The zoom gravity mesh truncates the tree forces at a larger distance "
        "(%.6e) than the background mesh (%.6e). Decrease "
        "Gravity:zoom_r_cut_max or increase Gravity:zoom_mesh_side_length.",
        mesh->r_cut_max, props->r_s * props->r_cut_max_ratio);
  }

  /* Positions are not wrapped on the zoom mesh, so it has to lie inside the
   * box and no pair of covered particles can be closer through the periodic
   * boundary than the cut-off of the background mesh. */
  for (int i = 0; i < 3; ++i) {
    if (mesh->loc[i] < 0. || mesh->loc[i] + mesh->dim[i] > s->dim[i] ||
        mesh->dim[i] > s->dim[i] - props->r_s * props->r_cut_max_ratio) {
      error(
          "The zoom gravity mesh ([%.6e, %.6e] along axis %d) must lie inside "
          "the box (size %.6e) and be smaller than the box by at least the "
          "cut-off of the background mesh (%.6e). Decrease "
          "Gravity:zoom_mesh_bkg_buffer_cells or shrink the zoom region.",
          mesh->loc[i], mesh->loc[i] + mesh->dim[i], i, s->dim[i],
          props->r_s * props->r_cut_max_ratio);
    }
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
  if (!s->periodic) {
    error(
        "Can't run with a high resolution mesh without periodic boundary "
        "conditions to support the low resolution mesh.");
  }

  /* Ensure that the zoom region has been initialized. */
  if (s->zoom_props == NULL) {
    error("Gravity:zoom_mesh requested before zoom region initialization.");
  }

  /* Unpack the active mesh parameters and sanity check them. The user-provided
   * side length refers to the inner N^3 zoom mesh, not the zero-padded FFT
   * domain used by the zero padded Hockney-Eastwood method. */
  const int side_length =
      parser_get_param_int(params, "Gravity:zoom_mesh_side_length");
  if (side_length < 8) {
    error("Gravity:zoom_mesh_side_length must be at least 8.");
  }

  /* The zoom mesh is not distributed over MPI ranks */
  if (props->distributed_mesh) {
    error("The zoom gravity mesh does not support Gravity:distributed_mesh.");
  }

  /* Replicated transforms retain the global mesh's practical serial limit. */
  if (2 * side_length > 1290) {
    error(
        "Gravity:zoom_mesh_side_length is too large for the current "
        "zero-padded zoom mesh implementation. The active side length N is "
        "%d, but the zero padded Hockney-Eastwood FFT uses side length 2N=%d. "
        "Use N <= 645.",
        side_length, 2 * side_length);
  }

  const double r_cut_max_ratio =
      parser_get_param_double(params, "Gravity:zoom_r_cut_max");
  if (r_cut_max_ratio <= 0.) error("Gravity:zoom_r_cut_max must be positive.");

  const double r_cut_min_ratio =
      parser_get_opt_param_double(params, "Gravity:zoom_r_cut_min", 0.);
  if (r_cut_min_ratio < 0.) error("Gravity:zoom_r_cut_min must be >= 0.");
  if (r_cut_min_ratio > r_cut_max_ratio) {
    error("Gravity:zoom_r_cut_min must be <= Gravity:zoom_r_cut_max.");
  }

  /* The truncated tree force must fit inside the mesh domain, mirroring the
   * background mesh's own r_cut_max check in gravity_props_init. */
  if (2. * props->a_smooth * r_cut_max_ratio > side_length) {
    error(
        "Zoom mesh too small given Gravity:zoom_r_cut_max. "
        "Gravity:zoom_mesh_side_length should be at least %d cells wide "
        "(currently: %d).",
        (int)(2. * props->a_smooth * r_cut_max_ratio) + 1, side_length);
  }

  const int buffer_cells =
      parser_get_opt_param_int(params, "Gravity:zoom_mesh_bkg_buffer_cells", 1);
  if (buffer_cells < 1) {
    error("Gravity:zoom_mesh_bkg_buffer_cells must be >= 1.");
  }

  if (s->zoom_props->nr_bkg_cells <= 0 ||
      s->zoom_props->bkg_cells_top == NULL) {
    error("Gravity:zoom_mesh requires background top-level cells.");
  }
  if (s->zoom_props->nr_zoom_cells <= 0 ||
      s->zoom_props->zoom_cells_top == NULL) {
    error("Gravity:zoom_mesh requires zoom top-level cells.");
  }

  /* Fill in the mesh structure. */
  mesh->enabled = 1;
  mesh->nr_threads = nr_threads;
  mesh->use_local_patches = props->mesh_uses_local_patches;
  mesh->buffer_cells = buffer_cells;
  mesh->N = side_length;
  zoom_mesh_set_geometry(mesh, s, props, r_cut_max_ratio, r_cut_min_ratio);

  /* Report the mesh geometry to the user if requested. */
  if (verbose) {
    message(
        "Zoom gravity mesh: enabled with N=%d, loc=[%.6e %.6e %.6e], "
        "dim=[%.6e %.6e %.6e], buffer_cells=%d, r_cut_max=%.6e.",
        mesh->N, mesh->loc[0], mesh->loc[1], mesh->loc[2], mesh->dim[0],
        mesh->dim[1], mesh->dim[2], buffer_cells, mesh->r_cut_max);
  }

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

  /* The active mesh is embedded in a zero-padded mesh twice as large */
  const size_t N = 2 * mesh->N;

  /* Allocate the memory for the combined density and potential array */
  mesh->potential_global = (double *)fftw_malloc(sizeof(double) * N * N * N);
  if (mesh->potential_global == NULL) {
    error("Error allocating memory for the zoom gravity mesh.");
  }
  memuse_log_allocation("fftw_zoom_mesh.potential", mesh->potential_global, 1,
                        sizeof(double) * N * N * N);
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

  if (mesh->potential_global) {
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
