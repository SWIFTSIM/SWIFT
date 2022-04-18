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
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "mesh_gravity_mpi.h"

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

/**
 * @brief Accumulate contributions from cell to density field
 *
 * Allocates a temporary mesh which covers the top level cell,
 * accumulates mass contributions to this mesh, and then
 * adds these contributions to the supplied hashmap.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param dim The dimensions of the simulation box.
 * @param cell The #cell containing the particles.
 * @param patch The local mesh patch
 * @param nu_model Struct with neutrino constants
 *
 */
void accumulate_cell_to_local_patch(const int N, const double fac,
                                    const double *dim, const struct cell *cell,
                                    struct pm_mesh_patch *patch,
                                    const struct neutrino_model *nu_model) {

  /* If the cell is empty, then there's nothing to do
     (and the code to find the extent of the cell would fail) */
  if (cell->grav.count == 0) return;

  /* Initialise the local mesh patch */
  pm_mesh_patch_init(patch, cell, N, fac, dim, /*boundary_size=*/1);
  pm_mesh_patch_zero(patch);

  /* Loop over particles in this cell */
  for (int ipart = 0; ipart < cell->grav.count; ipart++) {

    const struct gpart *gp = &(cell->grav.parts[ipart]);

    /* Box wrap the particle's position to the copy nearest the cell centre */
    const double pos_x =
        box_wrap(gp->x[0], patch->wrap_min[0], patch->wrap_max[0]);
    const double pos_y =
        box_wrap(gp->x[1], patch->wrap_min[1], patch->wrap_max[1]);
    const double pos_z =
        box_wrap(gp->x[2], patch->wrap_min[2], patch->wrap_max[2]);

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

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct accumulate_mapper_data {
  const struct cell *cells;
  const int *local_cells;
  struct pm_mesh_patch *local_patches;
  int N;
  double fac;
  double dim[3];
  struct neutrino_model *nu_model;
};

/**
 * @brief Threadpool mapper function for the mesh CIC assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
void accumulate_cell_to_local_patches_mapper(void *map_data, int num,
                                             void *extra) {

  /* Unpack the shared information */
  const struct accumulate_mapper_data *data =
      (struct accumulate_mapper_data *)extra;
  const struct cell *cells = data->cells;
  const int N = data->N;
  const double fac = data->fac;
  const double dim[3] = {data->dim[0], data->dim[1], data->dim[2]};
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

    /* Assign this cell's content to the mesh */
    accumulate_cell_to_local_patch(N, fac, dim, c, &local_patches[i], nu_model);
  }
}

/**
 * @brief Accumulate local contributions to the density field
 *
 * Fill the array of local patches with the data corresponding
 * to the local top-level cells.
 * The patches are stored in the order of the space->local_cells_top list.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param s The #space containing the particles.
 * @param local_patches The array of *local* mesh patches.
 *
 */
void mpi_mesh_accumulate_gparts_to_local_patches(
    struct threadpool *tp, const int N, const double fac, const struct space *s,
    struct pm_mesh_patch *local_patches) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Gather some neutrino constants if using delta-f weighting on the mesh */
  struct neutrino_model nu_model;
  bzero(&nu_model, sizeof(struct neutrino_model));
  if (s->e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(s, &nu_model);

  /* Use the threadpool to parallelize over cells */
  struct accumulate_mapper_data data;
  data.cells = s->cells_top;
  data.local_cells = local_cells;
  data.local_patches = local_patches;
  data.N = N;
  data.fac = fac;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];
  data.nu_model = &nu_model;
  threadpool_map(tp, accumulate_cell_to_local_patches_mapper,
                 (void *)local_cells, nr_local_cells, sizeof(int),
                 threadpool_auto_chunk_size, (void *)&data);
  if (lock_destroy(&lock) != 0) error("Impossible to destroy lock!");

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

void mesh_patches_to_sorted_array(const struct pm_mesh_patch *local_patches,
                                  const int nr_patches,
                                  struct mesh_key_value_rho *array,
                                  const size_t size) {

  size_t count = 0;
  for (int p = 0; p < nr_patches; ++p) {
    const struct pm_mesh_patch *patch = &local_patches[p];

    /* Loop over all cells in the patch */
    for (int i = 0; i < patch->mesh_size[0]; i++) {
      for (int j = 0; j < patch->mesh_size[1]; j++) {
        for (int k = 0; k < patch->mesh_size[2]; k++) {

          /* Find array index in the mesh patch */
          const int local_index = pm_mesh_patch_index(patch, i, j, k);

          /* Find index in the full mesh */
          const size_t global_index = row_major_id_periodic_size_t_padded(
              i + patch->mesh_min[0], j + patch->mesh_min[1],
              k + patch->mesh_min[2], patch->N);

          /* Get the value */
          const double value = patch->mesh[local_index];

          /* Store everything in the flattened array using
           * the global index as a key */
          array[count].value = value;
          array[count].key = global_index;

          ++count;
        }
      }
    }
  }

  /* quick check... */
  if (count != size) error("Error flattening the mesh patches!");
}

/**
 * @brief Convert the array of local patches to a slab-distributed 3D mesh
 *
 * For FFTW each rank needs to hold a slice of the full mesh.
 * This routine does the necessary communication to convert
 * the per-rank local patches into a slab-distributed mesh.
 *
 * This function will clean the memory allocated by each of the entry
 * in the local_patches array.
 *
 * @param N The size of the mesh
 * @param local_n0 The thickness of the slice to store on this rank
 * @param local_patches The array of local patches.
 * @param nr_patches The number of local patches.
 * @param mesh Pointer to the output data buffer.
 * @param tp The #threadpool object.
 * @param verbose Are we talkative?
 */
void mpi_mesh_local_patches_to_slices(const int N, const int local_n0,
                                      struct pm_mesh_patch *local_patches,
                                      const int nr_patches, double *mesh,
                                      struct threadpool *tp,
                                      const int verbose) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  ticks tic = getticks();

  /* Count the total number of mesh cells we have.
   *
   * Note: There might be duplicates. We don't care at this point. */
  size_t count = 0;
  for (int i = 0; i < nr_patches; ++i) {
    const struct pm_mesh_patch *p = &local_patches[i];
    count += p->mesh_size[0] * p->mesh_size[1] * p->mesh_size[2];
  }

  /* Create an array to contain all the individual mesh cells we have
   * on this node. For now, this is in random order */
  struct mesh_key_value_rho *mesh_sendbuf_unsorted;
  if (swift_memalign("mesh_sendbuf_unsorted", (void **)&mesh_sendbuf_unsorted,
                     SWIFT_CACHE_ALIGNMENT,
                     count * sizeof(struct mesh_key_value_rho)) != 0)
    error("Failed to allocate array for unsorted mesh send buffer!");

  /* Make an array with the (key, value) pairs from the mesh patches.
   *
   * We're going to distribute them between ranks according to their
   * x coordinate, so we later need to put them in order of destination rank. */
  mesh_patches_to_sorted_array(local_patches, nr_patches, mesh_sendbuf_unsorted,
                               count);

  if (verbose)
    message(" - Converting mesh patches to array took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Clean the local patches array */
  for (int i = 0; i < nr_patches; ++i) pm_mesh_patch_clean(&local_patches[i]);

  tic = getticks();

  /* Now, create space for a sorted version of the array of mesh cells */
  struct mesh_key_value_rho *mesh_sendbuf;
  if (swift_memalign("mesh_sendbuf", (void **)&mesh_sendbuf,
                     SWIFT_CACHE_ALIGNMENT,
                     count * sizeof(struct mesh_key_value_rho)) != 0)
    error("Failed to allocate array for unsorted mesh send buffer!");

  size_t *sorted_offsets = (size_t *)malloc(N * sizeof(size_t));

  /* Do a bucket sort of the mesh elements to have them sorted
   * by global x-coordinate (note we don't care about y,z at this stage)
   * Also reover the offsets where we switch from one bin to the next */
  bucket_sort_mesh_key_value_rho(mesh_sendbuf_unsorted, count, N, tp,
                                 mesh_sendbuf, sorted_offsets);

  /* Let's free the unsorted array to keep things lean */
  swift_free("mesh_sendbuf_unsorted", mesh_sendbuf_unsorted);
  mesh_sendbuf_unsorted = NULL;

  if (verbose)
    message(" - Sorting of mesh cells took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Get width of the slice on each rank */
  int *slice_width = (int *)malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&local_n0, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine offset to the slice on each rank */
  int *slice_offset = (int *)malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i++) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Compute how many elements are to be sent to each rank */
  size_t *nr_send = (size_t *)calloc(nr_nodes, sizeof(size_t));

  /* Loop over the offsets */
  int dest_node = 0;
  for (int i = 0; i < N; ++i) {

    /* Find the first mesh cell in that bucket */
    const size_t j = sorted_offsets[i];

    /* Get the x coordinate of this mesh cell in the global mesh */
    const int mesh_x =
        get_xcoord_from_padded_row_major_id((size_t)mesh_sendbuf[j].key, N);

    /* Advance to the destination node that is to contain this x coordinate */
    while ((mesh_x >= slice_offset[dest_node] + slice_width[dest_node]) ||
           (slice_width[dest_node] == 0)) {
      dest_node++;
    }

    /* Add all the mesh cells in this bucket */
    if (i < N - 1)
      nr_send[dest_node] += sorted_offsets[i + 1] - sorted_offsets[i];
    else
      nr_send[dest_node] += count - sorted_offsets[i];
  }

#ifdef SWIFT_DEBUG_CHECKS
  size_t *nr_send_check = (size_t *)calloc(nr_nodes, sizeof(size_t));

  /* Brute-force list without using the offsets */
  int dest_node_check = 0;
  for (size_t i = 0; i < count; i++) {
    /* Get the x coordinate of this mesh cell in the global mesh */
    const int mesh_x =
        get_xcoord_from_padded_row_major_id((size_t)mesh_sendbuf[i].key, N);
    /* Advance to the destination node that is to contain this x coordinate */
    while ((mesh_x >=
            slice_offset[dest_node_check] + slice_width[dest_node_check]) ||
           (slice_width[dest_node_check] == 0)) {
      dest_node_check++;
    }
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

  /* Allocate the receive buffer */
  struct mesh_key_value_rho *mesh_recvbuf;
  if (swift_memalign("mesh_recvbuf", (void **)&mesh_recvbuf,
                     SWIFT_CACHE_ALIGNMENT,
                     nr_recv_tot * sizeof(struct mesh_key_value_rho)) != 0)
    error("Failed to allocate receive buffer for constructing MPI FFT mesh");

  if (verbose)
    message(" - Preparing comms buffers took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Carry out the communication */
  exchange_structs(nr_send, (char *)mesh_sendbuf, nr_recv, (char *)mesh_recvbuf,
                   sizeof(struct mesh_key_value_rho));

  if (verbose)
    message(" - MPI exchange took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Copy received data to the output buffer.
   * This is now a local slice of the global mesh. */
  for (size_t i = 0; i < nr_recv_tot; i++) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Verify that we indeed got a cell that should be in the local mesh slice
     */
    const int xcoord =
        get_xcoord_from_padded_row_major_id(mesh_recvbuf[i].key, N);
    if (xcoord < slice_offset[nodeID])
      error("Received mesh cell is not in the local slice (xcoord too small)");
    if (xcoord >= slice_offset[nodeID] + slice_width[nodeID])
      error("Received mesh cell is not in the local slice (xcoord too large)");
#endif

    /* What cell are we looking at? */
    const size_t local_index = get_index_in_local_slice(
        (size_t)mesh_recvbuf[i].key, N, slice_offset[nodeID]);

    /* Add to the cell*/
    mesh[local_index] += mesh_recvbuf[i].value;
  }

  if (verbose)
    message(" - Filling of the density values took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Tidy up */
  free(slice_width);
  free(slice_offset);
  free(nr_send);
  free(nr_recv);
  swift_free("mesh_recvbuf", mesh_recvbuf);
  swift_free("mesh_sendbuf", mesh_sendbuf);
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

/**
 * @brief Count the number of mesh cells we will need to request from other
 * nodes
 *
 * @param N the mesh size.
 * @param fac Inverse of the FFT mesh cell size
 * @param s The #space containing the particles.
 */
size_t count_required_mesh_cells(const int N, const double fac,
                                 const struct space *s) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  size_t count = 0;

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell++) {

    struct cell *cell = &(s->cells_top[local_cells[icell]]);

    /* Skip empty cells */
    if (cell->grav.count == 0) continue;

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
      const double xmin = cell->loc[idim] - 2.01 / fac;
      const double xmax = cell->loc[idim] + cell->width[idim] + 3.01 / fac;
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

size_t init_required_mesh_cells(const int N, const double fac,
                                const struct space *s,
                                struct mesh_key_value_pot *send_cells) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  size_t count = 0;

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell++) {

    struct cell *cell = &(s->cells_top[local_cells[icell]]);

    /* Skip empty cells */
    if (cell->grav.count == 0) continue;

    /* Determine range of FFT mesh cells we need for particles in this top
       level cell. The 5 point stencil used for accelerations requires
       2 neighbouring FFT mesh cells in each direction and for CIC
       evaluation of the accelerations we need one extra FFT mesh cell
       in the +ve direction.

       We also have to add a small buffer to avoid problems with rounding

       TODO: can we calculate exactly how big the rounding error can be?
       Will just assume that 1% of a mesh cell is enough for now.
    */
    int ixmin[3];
    int ixmax[3];
    for (int idim = 0; idim < 3; idim++) {
      const double xmin = cell->loc[idim] - 2.01 / fac;
      const double xmax = cell->loc[idim] + cell->width[idim] + 3.01 / fac;
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

void fill_local_patches_from_mesh_cells(
    const int N, const double fac, const struct space *s,
    const struct mesh_key_value_pot *mesh_cells,
    struct pm_mesh_patch *local_patches, const size_t nr_send_tot) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  size_t offset = 0;

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell++) {

    const struct cell *cell = &(s->cells_top[local_cells[icell]]);
    struct pm_mesh_patch *patch = &local_patches[icell];

    /* Skip empty cells */
    if (cell->grav.count == 0) continue;

    patch->N = N;
    patch->fac = fac;

    /* Will need to wrap particles to position nearest the cell centre */
    for (int i = 0; i < 3; i++) {
      patch->wrap_min[i] = cell->loc[i] + 0.5 * cell->width[i] - 0.5 * dim[i];
      patch->wrap_max[i] = cell->loc[i] + 0.5 * cell->width[i] + 0.5 * dim[i];
    }

    int num_cells = 1;
    for (int i = 0; i < 3; i++) {
      const double xmin = cell->loc[i] - 2.01 / fac;
      const double xmax = cell->loc[i] + cell->width[i] + 3.01 / fac;
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

/**
 * @brief Retrieve the potential in the mesh cells we need to
 * compute the force on particles on this MPI rank. Result is
 * returned in the supplied hashmap, which should be initially
 * empty.
 *
 * We need all cells containing points -2 and +3 mesh cell widths
 * away from each particle along each axis to compute the
 * potential gradient.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the FFT mesh cell size
 * @param s The #space containing the particles.
 * @param local_0_start Offset to the first mesh x coordinate on this rank
 * @param local_n0 Width of the mesh slab on this rank
 * @param potential_slice Array with the potential on the local slice of the
 * mesh
 * @param tp The #threadpool object.
 * @param verbose Are we talkative?
 */
void mpi_mesh_fetch_potential(const int N, const double fac,
                              const struct space *s, const int local_0_start,
                              const int local_n0, double *potential_slice,
                              struct pm_mesh_patch *local_patches,
                              struct threadpool *tp, const int verbose) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* Determine rank, number of MPI ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  ticks tic = getticks();

  /* Determine how many mesh cells we will need to request */
  const size_t nr_send_tot = count_required_mesh_cells(N, fac, s);

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
      init_required_mesh_cells(N, fac, s, send_cells_unsorted);

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
  fill_local_patches_from_mesh_cells(N, fac, s, send_cells_sorted,
                                     local_patches, nr_send_tot);

  if (verbose)
    message(" - Filling the local patches took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  swift_free("send_cells_sorted", send_cells_sorted);

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

/**
 * @brief Computes the potential on a gpart from a given mesh using the CIC
 * method.
 *
 * @param gp The #gpart.
 * @param patch The local mesh patch
 */
#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
void mesh_patch_to_gparts_CIC(struct gpart *gp,
                              const struct pm_mesh_patch *patch) {

  const double fac = patch->fac;

  /* Box wrap the gpart's position to the copy nearest the cell centre */
  const double pos_x =
      box_wrap(gp->x[0], patch->wrap_min[0], patch->wrap_max[0]);
  const double pos_y =
      box_wrap(gp->x[1], patch->wrap_min[1], patch->wrap_max[1]);
  const double pos_z =
      box_wrap(gp->x[2], patch->wrap_min[2], patch->wrap_max[2]);

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

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (gp->a_grav_mesh[0] != 0.) error("Particle with non-initalised stuff");
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  if (gp->potential_mesh != 0.) error("Particle with non-initalised stuff");
#endif
#endif

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

  /* Store things back */
  gp->a_grav_mesh[0] = fac * a[0];
  gp->a_grav_mesh[1] = fac * a[1];
  gp->a_grav_mesh[2] = fac * a[2];
  gravity_add_comoving_mesh_potential(gp, p);
}
#endif

/**
 * @brief Interpolate the forces and potential from the mesh to the #gpart.
 *
 * This is for the case where the mesh is distributed between MPI ranks
 * and stored in the form of a hashmap. This function updates the particles
 * in one #cell.
 *
 * @param c The #cell containing the #gpart to update
 * @param potential Hashmap containing the potential to interpolate from.
 * @param N Size of the full mesh
 * @param fac Inverse of the FFT mesh cell size
 * @param const_G Gravitional constant
 * @param dim Dimensions of the #space
 */
void cell_distributed_mesh_to_gpart_CIC(const struct cell *c,
                                        const struct pm_mesh_patch *patch,
                                        const int N, const double fac,
                                        const float const_G,
                                        const double dim[3]) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  const int gcount = c->grav.count;
  struct gpart *gparts = c->grav.parts;

  /* Check for empty cell as this would cause problems finding the extent */
  if (gcount == 0) return;

  /* Get the potential from the mesh patch to the active gparts using CIC */
  for (int i = 0; i < gcount; ++i) {
    struct gpart *gp = &gparts[i];

    if (gp->time_bin == time_bin_inhibited) continue;

    gp->a_grav_mesh[0] = 0.f;
    gp->a_grav_mesh[1] = 0.f;
    gp->a_grav_mesh[2] = 0.f;
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
    gp->potential_mesh = 0.f;
#endif

    mesh_patch_to_gparts_CIC(gp, patch);

    gp->a_grav_mesh[0] *= const_G;
    gp->a_grav_mesh[1] *= const_G;
    gp->a_grav_mesh[2] *= const_G;
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
    gp->potential_mesh *= const_G;
#endif
  }

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct distributed_cic_mapper_data {
  const struct cell *cells;
  const int *local_cells;
  const struct pm_mesh_patch *local_patches;
  int N;
  double fac;
  double dim[3];
  float const_G;
};

/**
 * @brief Threadpool mapper function for the mesh CIC assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
void cell_distributed_mesh_to_gpart_CIC_mapper(void *map_data, int num,
                                               void *extra) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* Unpack the shared information */
  const struct distributed_cic_mapper_data *data =
      (struct distributed_cic_mapper_data *)extra;
  const struct cell *cells = data->cells;
  const int N = data->N;
  const double fac = data->fac;
  const double dim[3] = {data->dim[0], data->dim[1], data->dim[2]};
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
    cell_distributed_mesh_to_gpart_CIC(c, &local_patches[i], N, fac, const_G,
                                       dim);
  }

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

void mpi_mesh_update_gparts(struct pm_mesh_patch *local_patches,
                            const struct space *s, struct threadpool *tp,
                            const int N, const double cell_fac) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Gather the mesh shared information to be used by the threads */
  struct distributed_cic_mapper_data data;
  data.cells = s->cells_top;
  data.local_cells = local_cells;
  data.local_patches = local_patches;
  data.N = N;
  data.fac = cell_fac;
  data.dim[0] = s->dim[0];
  data.dim[1] = s->dim[1];
  data.dim[2] = s->dim[2];
  data.const_G = s->e->physical_constants->const_newton_G;

  if (nr_local_cells == 0) {
    error("Distributed mesh not implemented without cells");
  } else {
    /* Evaluate acceleration and potential for each gpart */
    threadpool_map(tp, cell_distributed_mesh_to_gpart_CIC_mapper,
                   (void *)local_cells, nr_local_cells, sizeof(int),
                   threadpool_auto_chunk_size, (void *)&data);
  }
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}
