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
 *******************************************************************************/

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "zoom_mesh_gravity_mpi.h"

/* Local includes. */
#include "active.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "exchange_structs.h"
#include "lock.h"
#include "mesh_gravity_patch.h"
#include "mesh_gravity_sort.h"
#include "neutrino.h"
#include "part.h"
#include "periodic.h"
#include "row_major_id.h"
#include "space.h"
#include "threadpool.h"
#include "zoom_mesh_gravity.h"

/**
 * @brief Initialize a mesh patch to cover a cell on the zoom mesh
 *
 * Zoom mesh equivalent of pm_mesh_patch_init(). The patch coordinates are
 * indices into the zero-padded zoom mesh and the particles are not wrapped.
 *
 * @param patch A pointer to the mesh patch
 * @param mesh The #zoom_pm_mesh
 * @param cell The cell which the mesh should cover
 * @param boundary_size Size of the boundary layer to include
 */
static void zoom_mesh_patch_init(struct pm_mesh_patch *patch,
                                 const struct zoom_pm_mesh *mesh,
                                 const struct cell *cell,
                                 const int boundary_size) {

  const int gcount = cell->grav.count;
  const struct gpart *gparts = cell->grav.parts;
  const double fac = mesh->cell_fac[0];

  /* The FFTs are performed on the zero-padded mesh */
  patch->N = 2 * mesh->N[0];
  patch->fac = fac;

  /* No wrapping on the isolated mesh, positions are relative to its corner */
  for (int i = 0; i < 3; i++) {
    patch->wrap_min[i] = 0.;
    patch->wrap_max[i] = mesh->dim[i];
  }

  /* Find the extent of the particle distribution in the cell */
  double pos_min[3];
  double pos_max[3];
  for (int i = 0; i < 3; i++) {
    pos_min[i] = patch->wrap_max[i];
    pos_max[i] = patch->wrap_min[i];
  }
  for (int ipart = 0; ipart < gcount; ipart++) {

    const struct gpart *gp = &gparts[ipart];

    if (gp->time_bin == time_bin_inhibited) continue;

    for (int i = 0; i < 3; i++) {
      const double pos = gp->x[i] - mesh->loc[i];
      if (pos < pos_min[i]) pos_min[i] = pos;
      if (pos > pos_max[i]) pos_max[i] = pos;
    }
  }

  /* Determine the integer size and coordinates of the mesh */
  int num_cells = 1;
  for (int i = 0; i < 3; i++) {
    patch->mesh_min[i] = floor(pos_min[i] * fac) - boundary_size;
    /* CIC interpolation requires one extra element in the positive direction */
    patch->mesh_max[i] = floor(pos_max[i] * fac) + boundary_size + 1;
    patch->mesh_size[i] = patch->mesh_max[i] - patch->mesh_min[i] + 1;
    num_cells *= patch->mesh_size[i];
  }

  /* Allocate the mesh */
  if (swift_memalign("mesh_patch", (void **)&patch->mesh, SWIFT_CACHE_ALIGNMENT,
                     num_cells * sizeof(double)) != 0)
    error("Failed to allocate array for mesh patch!");
}

/**
 * @brief Accumulate contributions from cell to the zoom mesh density field
 *
 * Zoom mesh equivalent of accumulate_cell_to_local_patch().
 *
 * Allocates a temporary mesh which covers the top level cell and
 * accumulates mass contributions to this mesh.
 *
 * @param mesh The #zoom_pm_mesh
 * @param cell The #cell containing the particles.
 * @param patch The local mesh patch
 * @param nu_model Struct with neutrino constants
 *
 */
void zoom_mesh_accumulate_cell_to_local_patch(
    const struct zoom_pm_mesh *mesh, const struct cell *cell,
    struct pm_mesh_patch *patch, const struct neutrino_model *nu_model) {

  /* If the cell is empty, then there's nothing to do
     (and the code to find the extent of the cell would fail) */
  if (cell->grav.count == 0) return;

  const double fac = mesh->cell_fac[0];

  /* Initialise the local mesh patch */
  zoom_mesh_patch_init(patch, mesh, cell, /*boundary_size=*/1);
  pm_mesh_patch_zero(patch);

  /* Loop over particles in this cell */
  for (int ipart = 0; ipart < cell->grav.count; ipart++) {

    const struct gpart *gp = &(cell->grav.parts[ipart]);

    if (gp->time_bin == time_bin_inhibited) continue;

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
    if (i < 0 || i >= mesh->N[0]) error("Invalid gpart position in x");
    if (j < 0 || j >= mesh->N[1]) error("Invalid gpart position in y");
    if (k < 0 || k >= mesh->N[2]) error("Invalid gpart position in z");
#endif

    /* Get coordinates within the mesh patch */
    const int ii = i - patch->mesh_min[0];
    const int jj = j - patch->mesh_min[1];
    const int kk = k - patch->mesh_min[2];

    /* Compute weight (for neutrino delta-f weighting) */
    double weight = 1.0;
    if (gp->type == swift_type_neutrino)
      gpart_neutrino_weight_mesh_only(gp, nu_model, &weight);

    /* Accumulate contributions to the local mesh patch */
    const double mass = gp->mass;
    const double value = mass * weight;
    pm_mesh_patch_CIC_set(patch, ii, jj, kk, tx, ty, tz, dx, dy, dz, value);
  }
}

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct zoom_accumulate_mapper_data {
  const struct cell *cells;
  const int *local_cells;
  struct pm_mesh_patch *local_patches;
  const struct zoom_pm_mesh *mesh;
  struct neutrino_model *nu_model;
};

/**
 * @brief Threadpool mapper function for the zoom mesh CIC assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
static void zoom_accumulate_cell_to_local_patches_mapper(void *map_data,
                                                         int num, void *extra) {

  /* Unpack the shared information */
  const struct zoom_accumulate_mapper_data *data =
      (struct zoom_accumulate_mapper_data *)extra;
  const struct cell *cells = data->cells;
  const struct zoom_pm_mesh *mesh = data->mesh;
  const struct neutrino_model *nu_model = data->nu_model;

  /* Pointer to the chunk to be processed */
  int *local_cells = (int *)map_data;

  /* Start at the same position in the list of local patches */
  const size_t offset = local_cells - data->local_cells;
  struct pm_mesh_patch *local_patches = data->local_patches + offset;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell *c = &cells[local_cells[i]];

    /* Skip empty cells */
    if (c->grav.count == 0) continue;

    /* Skip cells not covered by the zoom mesh */
    if (!zoom_mesh_cell_is_covered(mesh, c, /*use_max_dx=*/0)) continue;

    /* Assign this cell's content to the mesh */
    zoom_mesh_accumulate_cell_to_local_patch(mesh, c, &local_patches[i],
                                             nu_model);
  }
}

#endif

/**
 * @brief Accumulate local contributions to the zoom mesh density field
 *
 * Zoom mesh equivalent of mpi_mesh_accumulate_gparts_to_local_patches().
 *
 * Fill the array of local patches with the data corresponding
 * to the local top-level cells covered by the zoom mesh.
 * The patches are stored in the order of the space->local_cells_top list.
 *
 * @param tp The #threadpool object.
 * @param mesh The #zoom_pm_mesh
 * @param s The #space containing the particles.
 * @param local_patches The array of *local* mesh patches.
 *
 */
void zoom_mpi_mesh_accumulate_gparts_to_local_patches(
    struct threadpool *tp, const struct zoom_pm_mesh *mesh,
    const struct space *s, struct pm_mesh_patch *local_patches) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Gather some neutrino constants if using delta-f weighting on the mesh */
  struct neutrino_model nu_model;
  bzero(&nu_model, sizeof(struct neutrino_model));
  if (s->e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(s, &nu_model);

  /* Use the threadpool to parallelize over cells */
  struct zoom_accumulate_mapper_data data;
  data.cells = s->cells_top;
  data.local_cells = local_cells;
  data.local_patches = local_patches;
  data.mesh = mesh;
  data.nu_model = &nu_model;
  threadpool_map(tp, zoom_accumulate_cell_to_local_patches_mapper,
                 (void *)local_cells, nr_local_cells, sizeof(int),
                 threadpool_auto_chunk_size, (void *)&data);

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

/**
 * @brief Count the number of zoom mesh cells we will need to request from
 * other nodes
 *
 * Zoom mesh equivalent of count_required_mesh_cells().
 *
 * @param mesh The #zoom_pm_mesh
 * @param s The #space containing the particles.
 */
static size_t zoom_count_required_mesh_cells(const struct zoom_pm_mesh *mesh,
                                             const struct space *s) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const double fac = mesh->cell_fac[0];

  size_t count = 0;

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell++) {

    struct cell *cell = &(s->cells_top[local_cells[icell]]);

    /* Skip empty cells */
    if (cell->grav.count == 0) continue;

    /* Skip cells not covered by the zoom mesh */
    if (!zoom_mesh_cell_is_covered(mesh, cell, /*use_max_dx=*/0)) continue;

    /* Determine range of FFT mesh cells we need for particles in this top
     * level cell. The 5 point stencil used for accelerations requires
     * 2 neighbouring FFT mesh cells in each direction and for CIC
     * evaluation of the accelerations we need one extra FFT mesh cell
     * in the +ve direction.
     *
     * We also have to add a small buffer to avoid problems with rounding
     *
     * TODO: can we calculate exactly how big the rounding error can be?
     * Will just assume that 1% of a mesh cell is enough for now.*/
    int ixmin[3];
    int ixmax[3];
    for (int idim = 0; idim < 3; idim++) {
      const double xmin = cell->loc[idim] - mesh->loc[idim] - 2.01 / fac;
      const double xmax =
          cell->loc[idim] - mesh->loc[idim] + cell->width[idim] + 3.01 / fac;
      ixmin[idim] = (int)floor(xmin * fac);
      ixmax[idim] = (int)floor(xmax * fac);
    }

    const int delta_i = (ixmax[0] - ixmin[0]) + 1;
    const int delta_j = (ixmax[1] - ixmin[1]) + 1;
    const int delta_k = (ixmax[2] - ixmin[2]) + 1;

#ifdef SWIFT_DEBUG_CHECKS
    if (delta_i > (1 << 12))
      error("Not enough bits to store local x-axis index");
    if (delta_j > (1 << 12))
      error("Not enough bits to store local y-axis index");
    if (delta_k > (1 << 12))
      error("Not enough bits to store local z-axis index");
#endif

    count += delta_i * delta_j * delta_k;
  }
  return count;
}

/**
 * @brief Initialise the zoom mesh cells we will need to request from other
 * nodes
 *
 * Zoom mesh equivalent of init_required_mesh_cells().
 *
 * @param mesh The #zoom_pm_mesh
 * @param s The #space containing the particles.
 * @param send_cells The array of mesh cells to request (to be filled).
 */
static size_t zoom_init_required_mesh_cells(
    const struct zoom_pm_mesh *mesh, const struct space *s,
    struct mesh_key_value_pot *send_cells) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const double fac = mesh->cell_fac[0];

  /* The FFTs are performed on the zero-padded mesh */
  const int N = 2 * mesh->N[0];

  size_t count = 0;

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell++) {

    struct cell *cell = &(s->cells_top[local_cells[icell]]);

    /* Skip empty cells */
    if (cell->grav.count == 0) continue;

    /* Skip cells not covered by the zoom mesh */
    if (!zoom_mesh_cell_is_covered(mesh, cell, /*use_max_dx=*/0)) continue;

    /* Determine range of FFT mesh cells we need for particles in this top
     * level cell. The 5 point stencil used for accelerations requires
     * 2 neighbouring FFT mesh cells in each direction and for CIC
     * evaluation of the accelerations we need one extra FFT mesh cell
     * in the +ve direction.
     *
     * We also have to add a small buffer to avoid problems with rounding
     *
     * TODO: can we calculate exactly how big the rounding error can be?
     * Will just assume that 1% of a mesh cell is enough for now.*/
    int ixmin[3];
    int ixmax[3];
    for (int idim = 0; idim < 3; idim++) {
      const double xmin = cell->loc[idim] - mesh->loc[idim] - 2.01 / fac;
      const double xmax =
          cell->loc[idim] - mesh->loc[idim] + cell->width[idim] + 3.01 / fac;
      ixmin[idim] = (int)floor(xmin * fac);
      ixmax[idim] = (int)floor(xmax * fac);
    }

#ifdef SWIFT_DEBUG_CHECKS
    const int delta_i = (ixmax[0] - ixmin[0]) + 1;
    const int delta_j = (ixmax[1] - ixmin[1]) + 1;
    const int delta_k = (ixmax[2] - ixmin[2]) + 1;

    if (delta_i > (1 << 12))
      error("Not enough bits to store local x-axis index");
    if (delta_j > (1 << 12))
      error("Not enough bits to store local y-axis index");
    if (delta_k > (1 << 12))
      error("Not enough bits to store local z-axis index");
#endif

    /* Add the required cells to the map */
    for (int i = ixmin[0]; i <= ixmax[0]; i++) {
      for (int j = ixmin[1]; j <= ixmax[1]; j++) {
        for (int k = ixmin[2]; k <= ixmax[2]; k++) {
          const size_t index = row_major_id_periodic_size_t_padded(i, j, k, N);

          /* Indices relative to the patch */
          const int ii = i - ixmin[0];
          const int jj = j - ixmin[1];
          const int kk = k - ixmin[2];

          /* Generate a combined index */
          const size_t cell_index =
              cell_index_from_patch_index(icell, ii, jj, kk);

          send_cells[count].cell_index = cell_index;
          send_cells[count].value = 0.;
          send_cells[count].key = index;

#ifdef SWIFT_DEBUG_CHECKS
          int test_ci, test_ii, test_jj, test_kk;
          patch_index_from_cell_index(cell_index, &test_ci, &test_ii, &test_jj,
                                      &test_kk);
          if (icell != test_ci) error("Invalid cell_index!");
          if (ii != test_ii) error("Invalid ii!");
          if (jj != test_jj) error("Invalid jj!");
          if (kk != test_kk) error("Invalid kk!");
#endif

          ++count;
        }
      }
    }
  }

  return count;
}

/**
 * @brief Fill the local patches with the zoom mesh cells we received
 *
 * Zoom mesh equivalent of fill_local_patches_from_mesh_cells().
 *
 * @param mesh The #zoom_pm_mesh
 * @param s The #space containing the particles.
 * @param mesh_cells The mesh cells sorted by top-level cell.
 * @param local_patches The array of *local* mesh patches.
 * @param nr_send_tot The number of mesh cells.
 */
static void zoom_fill_local_patches_from_mesh_cells(
    const struct zoom_pm_mesh *mesh, const struct space *s,
    const struct mesh_key_value_pot *mesh_cells,
    struct pm_mesh_patch *local_patches, const size_t nr_send_tot) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const double fac = mesh->cell_fac[0];

  size_t offset = 0;

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell++) {

    const struct cell *cell = &(s->cells_top[local_cells[icell]]);
    struct pm_mesh_patch *patch = &local_patches[icell];

    /* Skip empty cells */
    if (cell->grav.count == 0) continue;

    /* Skip cells not covered by the zoom mesh */
    if (!zoom_mesh_cell_is_covered(mesh, cell, /*use_max_dx=*/0)) continue;

    /* The FFTs are performed on the zero-padded mesh */
    patch->N = 2 * mesh->N[0];
    patch->fac = fac;

    /* No wrapping on the isolated mesh, positions are relative to its corner */
    for (int i = 0; i < 3; i++) {
      patch->wrap_min[i] = 0.;
      patch->wrap_max[i] = mesh->dim[i];
    }

    int num_cells = 1;
    for (int i = 0; i < 3; i++) {
      const double xmin = cell->loc[i] - mesh->loc[i] - 2.01 / fac;
      const double xmax =
          cell->loc[i] - mesh->loc[i] + cell->width[i] + 3.01 / fac;
      patch->mesh_min[i] = (int)floor(xmin * fac);
      patch->mesh_max[i] = (int)floor(xmax * fac);
      patch->mesh_size[i] = patch->mesh_max[i] - patch->mesh_min[i] + 1;
      num_cells *= patch->mesh_size[i];
    }

    /* Allocate the mesh */
    if (swift_memalign("mesh_patch", (void **)&patch->mesh,
                       SWIFT_CACHE_ALIGNMENT, num_cells * sizeof(double)) != 0)
      error("Failed to allocate array for mesh patch!");

#ifdef SWIFT_DEBUG_CHECKS
    int count = 0;
    for (size_t imesh = 0; imesh < nr_send_tot; ++imesh) {

      const size_t cell_index = mesh_cells[imesh].cell_index;
      const int temp = cell_index_extract_patch_index(cell_index);
      if (temp == icell) ++count;
    }
    if (count != num_cells)
      error(
          "Invalide number of cells to fill the patch! icell=%d count=%d "
          "num_cells=%d",
          icell, count, num_cells);
#endif

    /* Now, we can start filling the mesh patch cells from the array of
     * key-index-value tuples */
    for (size_t imesh = offset; imesh < offset + num_cells; ++imesh) {

      const size_t cell_index = mesh_cells[imesh].cell_index;
      const double pot = mesh_cells[imesh].value;

      /* Recover the patch index (should be icell) and the i,j,k indices */
      int patch_index, i, j, k;
      patch_index_from_cell_index(cell_index, &patch_index, &i, &j, &k);

#ifdef SWIFT_DEBUG_CHECKS
      const int temp = cell_index_extract_patch_index(cell_index);

      if (patch_index != icell)
        error(
            "mesh cell not sorted in the right order icell=%d patch_index=%d "
            "cell_index=%zd i=%d j=%d k=%d imesh=%zd temp=%d",
            icell, patch_index, cell_index, i, j, k, imesh, temp);
#endif

      /* Flatten out the i,j,k index */
      const int local_index = i * (patch->mesh_size[1] * patch->mesh_size[2]) +
                              j * patch->mesh_size[2] + k;

      /* Store the potential */
      patch->mesh[local_index] = pot;
    }

    /* Move to the next cell */
    offset += num_cells;
  }
}

#endif

/**
 * @brief Retrieve the potential in the zoom mesh cells we need to
 * compute the force on particles on this MPI rank. Result is
 * returned in the supplied local patches.
 *
 * Zoom mesh equivalent of mpi_mesh_fetch_potential().
 *
 * We need all cells containing points -2 and +3 mesh cell widths
 * away from each particle along each axis to compute the
 * potential gradient.
 *
 * @param mesh The #zoom_pm_mesh
 * @param s The #space containing the particles.
 * @param local_0_start Offset to the first mesh x coordinate on this rank
 * @param local_n0 Width of the mesh slab on this rank
 * @param potential_slice Array with the potential on the local slice of the
 * mesh
 * @param local_patches The array of *local* mesh patches.
 * @param tp The #threadpool object.
 * @param verbose Are we talkative?
 */
void zoom_mpi_mesh_fetch_potential(const struct zoom_pm_mesh *mesh,
                                   const struct space *s,
                                   const int local_0_start, const int local_n0,
                                   double *potential_slice,
                                   struct pm_mesh_patch *local_patches,
                                   struct threadpool *tp, const int verbose) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* The FFTs are performed on the zero-padded mesh */
  const int N = 2 * mesh->N[0];

  /* Determine rank, number of MPI ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  ticks tic = getticks();

  /* Determine how many mesh cells we will need to request */
  const size_t nr_send_tot = zoom_count_required_mesh_cells(mesh, s);

  if (verbose)
    message(" - Counting required mesh patches took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  struct mesh_key_value_pot *send_cells_unsorted;
  if (swift_memalign("send_cells_unsorted", (void **)&send_cells_unsorted,
                     SWIFT_CACHE_ALIGNMENT,
                     nr_send_tot * sizeof(struct mesh_key_value_pot)) != 0)
    error("Failed to allocate array for cells to request!");

  tic = getticks();

  /* Initialise the mesh cells we will request */
  const size_t check_count =
      zoom_init_required_mesh_cells(mesh, s, send_cells_unsorted);

  if (nr_send_tot != check_count)
    error("Count and initialisation incompatible!");

  if (verbose)
    message(" - Init required mesh patches took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  struct mesh_key_value_pot *send_cells;
  if (swift_memalign("send_cells", (void **)&send_cells, SWIFT_CACHE_ALIGNMENT,
                     nr_send_tot * sizeof(struct mesh_key_value_pot)) != 0)
    error("Failed to allocate array for cells to request!");

  size_t *sorted_offsets = (size_t *)malloc(N * sizeof(size_t));

  /* Do a bucket sort of the mesh elements to have them sorted
   * by global x-coordinate (note we don't care about y,z at this stage) */
  bucket_sort_mesh_key_value_pot(send_cells_unsorted, nr_send_tot, N, tp,
                                 send_cells, sorted_offsets);

  swift_free("send_cells_unsorted", send_cells_unsorted);
  send_cells_unsorted = NULL;

  if (verbose)
    message(" - 1st mesh patches sort took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Get width of the mesh slice on each rank */
  int *slice_width = (int *)malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&local_n0, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine first mesh x coordinate stored on each rank */
  int *slice_offset = (int *)malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i++) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Count how many mesh cells we need to request from each MPI rank */
  size_t *nr_send = (size_t *)calloc(nr_nodes, sizeof(size_t));

  /* Loop over the offsets */
  int dest_node = 0;
  for (int i = 0; i < N; ++i) {

    /* Find the first mesh cell in that bucket */
    const size_t j = sorted_offsets[i];

    /* Get the x coordinate of this mesh cell in the global mesh */
    const int mesh_x =
        get_xcoord_from_padded_row_major_id((size_t)send_cells[j].key, N);

    /* Advance to the destination node that is to contain this x coordinate */
    while ((mesh_x >= slice_offset[dest_node] + slice_width[dest_node]) ||
           (slice_width[dest_node] == 0)) {
      dest_node++;
    }

    /* Add all the mesh cells in this bucket */
    if (i < N - 1)
      nr_send[dest_node] += sorted_offsets[i + 1] - sorted_offsets[i];
    else
      nr_send[dest_node] += nr_send_tot - sorted_offsets[i];
  }

#ifdef SWIFT_DEBUG_CHECKS
  size_t *nr_send_check = (size_t *)calloc(nr_nodes, sizeof(size_t));

  /* Brute-force list without using the offsets */
  int dest_node_check = 0;
  for (size_t i = 0; i < nr_send_tot; i++) {
    while (get_xcoord_from_padded_row_major_id(send_cells[i].key, N) >=
               (slice_offset[dest_node_check] + slice_width[dest_node_check]) ||
           slice_width[dest_node_check] == 0) {
      dest_node_check++;
    }
    if (dest_node_check >= nr_nodes || dest_node_check < 0)
      error("Destination node out of range");
    nr_send_check[dest_node_check]++;
  }

  /* Verify the "smart" list is as good as the brute-force one */
  for (int i = 0; i < nr_nodes; ++i) {
    if (nr_send[i] != nr_send_check[i]) error("Invalid send list!");
  }
  free(nr_send_check);
#endif

  /* We don't need the sorted offsets any more from here onwards */
  free(sorted_offsets);

  /* Determine how many requests we'll receive from each MPI rank */
  size_t *nr_recv = (size_t *)malloc(sizeof(size_t) * nr_nodes);
  MPI_Alltoall(nr_send, sizeof(size_t), MPI_BYTE, nr_recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t nr_recv_tot = 0;
  for (int i = 0; i < nr_nodes; i++) {
    nr_recv_tot += nr_recv[i];
  }

  if (verbose)
    message(" - Preparing comms buffers took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Allocate buffer to receive requests */
  struct mesh_key_value_pot *recv_cells;
  if (swift_memalign("recv_cells", (void **)&recv_cells, SWIFT_CACHE_ALIGNMENT,
                     nr_recv_tot * sizeof(struct mesh_key_value_pot)) != 0)
    error("Failed to allocate array for mesh receive buffer!");

  /* Send requests for cells to other ranks */
  exchange_structs(nr_send, (char *)send_cells, nr_recv, (char *)recv_cells,
                   sizeof(struct mesh_key_value_pot));

  if (verbose)
    message(" - 1st exchange took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Look up potential in the requested cells */
  for (size_t i = 0; i < nr_recv_tot; i++) {
#ifdef SWIFT_DEBUG_CHECKS
    const size_t cells_in_slab = ((size_t)N) * (2 * (N / 2 + 1));
    const size_t first_local_id = local_0_start * cells_in_slab;
    const size_t num_local_ids = local_n0 * cells_in_slab;
    if (recv_cells[i].key < first_local_id ||
        recv_cells[i].key >= first_local_id + num_local_ids) {
      error("Requested potential mesh cell ID is out of range");
    }
#endif
    const size_t local_id =
        get_index_in_local_slice(recv_cells[i].key, N, local_0_start);
#ifdef SWIFT_DEBUG_CHECKS
    const size_t Ns = N;
    if (local_id >= Ns * (2 * (Ns / 2 + 1)) * local_n0)
      error("Local potential mesh cell ID is out of range");
#endif
    recv_cells[i].value = potential_slice[local_id];
  }

  if (verbose)
    message(" - Filling of the potential values took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Return the results */
  exchange_structs(nr_recv, (char *)recv_cells, nr_send, (char *)send_cells,
                   sizeof(struct mesh_key_value_pot));

  if (verbose)
    message(" - 2nd exchange took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Tidy up */
  swift_free("recv_cells", recv_cells);
  free(slice_width);
  free(slice_offset);
  free(nr_send);
  free(nr_recv);

  struct mesh_key_value_pot *send_cells_sorted;
  if (swift_memalign("send_cells_sorted", (void **)&send_cells_sorted,
                     SWIFT_CACHE_ALIGNMENT,
                     nr_send_tot * sizeof(struct mesh_key_value_pot)) != 0)
    error("Failed to allocate array for cells to request!");

  /* Now sort the mesh cells by top-level cell index (i.e. by the patch they
   * belong to) */
  bucket_sort_mesh_key_value_pot_index(
      send_cells, nr_send_tot, s->nr_local_cells, tp, send_cells_sorted);

  swift_free("send_cells", send_cells);
  send_cells = NULL;

  if (verbose)
    message(" - 2nd mesh patches sort took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Initialise the local patches with the data we just received */
  zoom_fill_local_patches_from_mesh_cells(mesh, s, send_cells_sorted,
                                          local_patches, nr_send_tot);

  if (verbose)
    message(" - Filling the local patches took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  swift_free("send_cells_sorted", send_cells_sorted);

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

/**
 * @brief Adds the zoom mesh correction to a gpart from a given mesh patch
 * using the CIC method.
 *
 * Zoom mesh equivalent of mesh_patch_to_gparts_CIC(). The correction is added
 * to the forces and potential of the periodic mesh.
 *
 * @param gp The #gpart.
 * @param patch The local mesh patch
 * @param mesh The #zoom_pm_mesh
 * @param const_G Gravitional constant
 */
static void zoom_mesh_patch_to_gparts_CIC(struct gpart *gp,
                                          const struct pm_mesh_patch *patch,
                                          const struct zoom_pm_mesh *mesh,
                                          const float const_G) {

  const double fac = patch->fac;

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

  /* Some local accumulators */
  double p = 0.;
  double a[3] = {0.};

  /* Get coordinates within the mesh patch */
  const int ii = i - patch->mesh_min[0];
  const int jj = j - patch->mesh_min[1];
  const int kk = k - patch->mesh_min[2];

  /* Simple CIC for the potential itself */
  p += pm_mesh_patch_CIC_get(patch, ii, jj, kk, tx, ty, tz, dx, dy, dz);

  /* 5-point stencil along each axis for the accelerations */
  a[0] += (1. / 12.) *
          pm_mesh_patch_CIC_get(patch, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (2. / 3.) *
          pm_mesh_patch_CIC_get(patch, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] += (2. / 3.) *
          pm_mesh_patch_CIC_get(patch, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (1. / 12.) *
          pm_mesh_patch_CIC_get(patch, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);

  a[1] += (1. / 12.) *
          pm_mesh_patch_CIC_get(patch, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (2. / 3.) *
          pm_mesh_patch_CIC_get(patch, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] += (2. / 3.) *
          pm_mesh_patch_CIC_get(patch, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (1. / 12.) *
          pm_mesh_patch_CIC_get(patch, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);

  a[2] += (1. / 12.) *
          pm_mesh_patch_CIC_get(patch, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  a[2] -= (2. / 3.) *
          pm_mesh_patch_CIC_get(patch, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  a[2] += (2. / 3.) *
          pm_mesh_patch_CIC_get(patch, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  a[2] -= (1. / 12.) *
          pm_mesh_patch_CIC_get(patch, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);

  /* Add the correction to the existing mesh contribution */
  gp->a_grav_mesh[0] += const_G * fac * a[0];
  gp->a_grav_mesh[1] += const_G * fac * a[1];
  gp->a_grav_mesh[2] += const_G * fac * a[2];
  gravity_add_comoving_mesh_potential(gp, const_G * p);
}

/**
 * @brief Interpolate the zoom mesh correction from the mesh to the #gpart.
 *
 * Zoom mesh equivalent of cell_distributed_mesh_to_gpart_CIC(). This function
 * updates the particles in one #cell.
 *
 * @param c The #cell containing the #gpart to update
 * @param patch The local mesh patch to interpolate from.
 * @param mesh The #zoom_pm_mesh
 * @param const_G Gravitional constant
 */
static void zoom_cell_distributed_mesh_to_gpart_CIC(
    const struct cell *c, const struct pm_mesh_patch *patch,
    const struct zoom_pm_mesh *mesh, const float const_G) {

  const int gcount = c->grav.count;
  struct gpart *gparts = c->grav.parts;

  /* Check for empty cell as this would cause problems finding the extent */
  if (gcount == 0) return;

  /* Only cells covered by the zoom mesh receive the correction */
  if (!zoom_mesh_cell_is_covered(mesh, c, /*use_max_dx=*/0)) return;

  /* Get the potential from the mesh patch to the active gparts using CIC */
  for (int i = 0; i < gcount; ++i) {
    struct gpart *gp = &gparts[i];

    if (gp->time_bin == time_bin_inhibited) continue;

    zoom_mesh_patch_to_gparts_CIC(gp, patch, mesh, const_G);
  }
}

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct zoom_distributed_cic_mapper_data {
  const struct cell *cells;
  const int *local_cells;
  const struct pm_mesh_patch *local_patches;
  const struct zoom_pm_mesh *mesh;
  float const_G;
};

/**
 * @brief Threadpool mapper function for the zoom mesh CIC interpolation of a
 * cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
static void zoom_cell_distributed_mesh_to_gpart_CIC_mapper(void *map_data,
                                                           int num,
                                                           void *extra) {

  /* Unpack the shared information */
  const struct zoom_distributed_cic_mapper_data *data =
      (struct zoom_distributed_cic_mapper_data *)extra;
  const struct cell *cells = data->cells;
  const struct zoom_pm_mesh *mesh = data->mesh;
  const float const_G = data->const_G;

  /* Pointer to the chunk to be processed */
  int *local_cells = (int *)map_data;

  /* Start at the same position in the list of local patches */
  const size_t offset = local_cells - data->local_cells;
  const struct pm_mesh_patch *local_patches = data->local_patches + offset;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell *c = &cells[local_cells[i]];

    /* Update acceleration and potential for gparts in this cell */
    zoom_cell_distributed_mesh_to_gpart_CIC(c, &local_patches[i], mesh,
                                            const_G);
  }
}

#endif

/**
 * @brief Interpolate the zoom mesh correction from the local patches to the
 * #gpart.
 *
 * Zoom mesh equivalent of mpi_mesh_update_gparts().
 *
 * @param local_patches The array of *local* mesh patches.
 * @param mesh The #zoom_pm_mesh
 * @param s The #space containing the particles.
 * @param tp The #threadpool object.
 */
void zoom_mpi_mesh_update_gparts(struct pm_mesh_patch *local_patches,
                                 const struct zoom_pm_mesh *mesh,
                                 const struct space *s, struct threadpool *tp) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Gather the mesh shared information to be used by the threads */
  struct zoom_distributed_cic_mapper_data data;
  data.cells = s->cells_top;
  data.local_cells = local_cells;
  data.local_patches = local_patches;
  data.mesh = mesh;
  data.const_G = s->e->physical_constants->const_newton_G;

  if (nr_local_cells == 0) {
    error("Distributed mesh not implemented without cells");
  } else {
    /* Evaluate acceleration and potential for each gpart */
    threadpool_map(tp, zoom_cell_distributed_mesh_to_gpart_CIC_mapper,
                   (void *)local_cells, nr_local_cells, sizeof(int),
                   threadpool_auto_chunk_size, (void *)&data);
  }
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}
