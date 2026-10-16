/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will J. Roper (w.roper@sussex.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *******************************************************************************/

#include "zoom_mesh_gravity_mpi.h"

#include "error.h"
#include "zoom_mesh_gravity.h"

#include <config.h>

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

#include "cell.h"
#include "clocks.h"
#include "engine.h"
#include "exchange_structs.h"
#include "memuse.h"
#include "mesh_gravity.h"
#include "neutrino.h"
#include "part.h"
#include "space.h"
#include "threadpool.h"

#include <fftw3-mpi.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/** A private, non-periodic mesh patch belonging to one local top-level cell. */
struct zoom_mpi_patch {
  struct cell *cell;
  int min[3];
  int size[3];
  double *values;
};

/** Density contribution routed to the owner of its real-space x slab. */
struct zoom_mpi_density_value {
  size_t key;
  double value;
};

/** Potential request, retaining the destination within the requester's patch.
 */
struct zoom_mpi_potential_value {
  size_t key;
  size_t patch_offset;
  int patch;
  double value;
};

struct zoom_mpi_layout {
  int M;
  int complex_width;
  ptrdiff_t alloc_local;
  ptrdiff_t local_n0;
  ptrdiff_t local_0_start;
  ptrdiff_t local_n1;
  ptrdiff_t local_1_start;
  int nr_nodes;
  int node_id;
  int *width;
  int *offset;
};

struct zoom_mpi_patch_mapper_data {
  const struct zoom_pm_mesh *mesh;
  const struct neutrino_model *nu_model;
};

struct zoom_mpi_interpolate_mapper_data {
  const struct zoom_pm_mesh *mesh;
  float const_G;
};

struct zoom_mpi_fourier_mapper_data {
  fftw_complex *field;
  const fftw_complex *kernel;
  ptrdiff_t local_1_start;
  int M;
  int complex_width;
  double norm;
};

static void *zoom_mpi_malloc(const char *label, size_t count, size_t size) {
  void *p = swift_malloc(label, (count != 0 ? count : 1) * size);
  if (p == NULL) error("Failed to allocate %s.", label);
  return p;
}

static size_t zoom_mpi_index(const int i, const int j, const int k,
                             const int M) {
  return ((size_t)i * M + j) * M + k;
}

static size_t zoom_mpi_real_local_index(const int i, const int j, const int k,
                                        const struct zoom_mpi_layout *layout) {
  return ((size_t)(i - layout->local_0_start) * layout->M + j) *
             (2 * layout->complex_width) +
         k;
}

static size_t zoom_mpi_patch_index(const struct zoom_mpi_patch *patch,
                                   const int i, const int j, const int k) {
  return ((size_t)i * patch->size[1] + j) * patch->size[2] + k;
}

static int zoom_mpi_owner(const struct zoom_mpi_layout *layout, const int x) {
  for (int rank = 0; rank < layout->nr_nodes; ++rank)
    if (layout->width[rank] > 0 && x >= layout->offset[rank] &&
        x < layout->offset[rank] + layout->width[rank])
      return rank;
  error("No MPI slab owns zoom mesh x index %d.", x);
  return -1;
}

static void zoom_mpi_layout_init(struct zoom_mpi_layout *layout, const int M) {
  layout->M = M;
  layout->complex_width = M / 2 + 1;
  MPI_Comm_size(MPI_COMM_WORLD, &layout->nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &layout->node_id);
  layout->alloc_local = fftw_mpi_local_size_3d_transposed(
      M, M, layout->complex_width, MPI_COMM_WORLD, &layout->local_n0,
      &layout->local_0_start, &layout->local_n1, &layout->local_1_start);
  layout->width =
      zoom_mpi_malloc("zoom_mpi_slab_width", layout->nr_nodes, sizeof(int));
  layout->offset =
      zoom_mpi_malloc("zoom_mpi_slab_offset", layout->nr_nodes, sizeof(int));
  const int width = (int)layout->local_n0;
  const int offset = (int)layout->local_0_start;
  MPI_Allgather(&width, 1, MPI_INT, layout->width, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&offset, 1, MPI_INT, layout->offset, 1, MPI_INT,
                MPI_COMM_WORLD);
}

static void zoom_mpi_layout_clean(struct zoom_mpi_layout *layout) {
  swift_free("zoom_mpi_slab_width", layout->width);
  swift_free("zoom_mpi_slab_offset", layout->offset);
}

static int zoom_mpi_gpart_index(const struct zoom_pm_mesh *mesh,
                                const struct gpart *gp, int ind[3], double d[3],
                                double t[3], const int need_stencil) {
  for (int axis = 0; axis < 3; ++axis) {
    const double x = (gp->x[axis] - mesh->loc[axis]) * mesh->cell_fac[axis];
    ind[axis] = (int)floor(x);
    d[axis] = x - ind[axis];
    t[axis] = 1. - d[axis];
    if (need_stencil) {
      if (ind[axis] < 2 || ind[axis] >= mesh->N[axis] - 3) return 0;
    } else if (ind[axis] < 0 || ind[axis] >= mesh->N[axis] - 1) {
      return 0;
    }
  }
  return 1;
}

static void zoom_mpi_density_bounds(const struct cell *c,
                                    const struct zoom_pm_mesh *mesh, int min[3],
                                    int max[3], int *found) {
  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL)
        zoom_mpi_density_bounds(c->progeny[k], mesh, min, max, found);
    return;
  }
  if (!zoom_mesh_cell_is_covered(mesh, c, 0)) return;
  for (int n = 0; n < c->grav.count; ++n) {
    const struct gpart *gp = &c->grav.parts[n];
    if (gp->time_bin == time_bin_inhibited) continue;
    int ind[3];
    double d[3], t[3];
    if (!zoom_mpi_gpart_index(mesh, gp, ind, d, t, 0))
      error("Covered particle lies outside distributed zoom mesh.");
    for (int axis = 0; axis < 3; ++axis) {
      if (ind[axis] < min[axis]) min[axis] = ind[axis];
      if (ind[axis] + 1 > max[axis]) max[axis] = ind[axis] + 1;
    }
    *found = 1;
  }
}

static void zoom_mpi_deposit_leaf(const struct cell *c,
                                  const struct zoom_mpi_patch_mapper_data *data,
                                  struct zoom_mpi_patch *patch) {
  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL)
        zoom_mpi_deposit_leaf(c->progeny[k], data, patch);
    return;
  }
  if (!zoom_mesh_cell_is_covered(data->mesh, c, 0)) return;
  for (int n = 0; n < c->grav.count; ++n) {
    const struct gpart *gp = &c->grav.parts[n];
    if (gp->time_bin == time_bin_inhibited) continue;
    int ind[3];
    double d[3], t[3];
    if (!zoom_mpi_gpart_index(data->mesh, gp, ind, d, t, 0))
      error("Covered particle lies outside distributed zoom mesh.");
    double weight = 1.;
    if (gp->type == swift_type_neutrino)
      gpart_neutrino_weight_mesh_only(gp, data->nu_model, &weight);
    const double value = gp->mass * weight;
    const int i = ind[0] - patch->min[0];
    const int j = ind[1] - patch->min[1];
    const int k = ind[2] - patch->min[2];
#define ZOOM_MPI_ADD(ii, jj, kk, w) \
  patch->values[zoom_mpi_patch_index(patch, ii, jj, kk)] += value * (w)
    ZOOM_MPI_ADD(i, j, k, t[0] * t[1] * t[2]);
    ZOOM_MPI_ADD(i, j, k + 1, t[0] * t[1] * d[2]);
    ZOOM_MPI_ADD(i, j + 1, k, t[0] * d[1] * t[2]);
    ZOOM_MPI_ADD(i, j + 1, k + 1, t[0] * d[1] * d[2]);
    ZOOM_MPI_ADD(i + 1, j, k, d[0] * t[1] * t[2]);
    ZOOM_MPI_ADD(i + 1, j, k + 1, d[0] * t[1] * d[2]);
    ZOOM_MPI_ADD(i + 1, j + 1, k, d[0] * d[1] * t[2]);
    ZOOM_MPI_ADD(i + 1, j + 1, k + 1, d[0] * d[1] * d[2]);
#undef ZOOM_MPI_ADD
  }
}

static void zoom_mpi_build_density_patches_mapper(void *map_data, int num,
                                                  void *extra) {
  struct zoom_mpi_patch *patches = map_data;
  const struct zoom_mpi_patch_mapper_data *data = extra;
  for (int p = 0; p < num; ++p) {
    struct zoom_mpi_patch *patch = &patches[p];
    int max[3] = {INT_MIN, INT_MIN, INT_MIN};
    patch->min[0] = patch->min[1] = patch->min[2] = INT_MAX;
    int found = 0;
    zoom_mpi_density_bounds(patch->cell, data->mesh, patch->min, max, &found);
    if (!found) continue;
    size_t count = 1;
    for (int axis = 0; axis < 3; ++axis) {
      patch->size[axis] = max[axis] - patch->min[axis] + 1;
      count *= patch->size[axis];
    }
    patch->values =
        zoom_mpi_malloc("zoom_mpi_density_patch", count, sizeof(double));
    memset(patch->values, 0, count * sizeof(double));
    zoom_mpi_deposit_leaf(patch->cell, data, patch);
  }
}

static int zoom_mpi_collect_cells(const struct space *s,
                                  struct zoom_mpi_patch **patches_out) {
  if (!s->with_zoom_region || s->zoom_props == NULL)
    error("Distributed zoom mesh requires production zoom-region cell lists.");
  const int nz = s->zoom_props->nr_local_zoom_cells_with_particles;
  const int nb = s->zoom_props->nr_local_bkg_cells_with_particles;
  const int count = nz + nb;
  struct zoom_mpi_patch *patches =
      zoom_mpi_malloc("zoom_mpi_patches", count, sizeof(struct zoom_mpi_patch));
  memset(patches, 0, (count != 0 ? count : 1) * sizeof(struct zoom_mpi_patch));
  for (int i = 0; i < nz; ++i)
    patches[i].cell =
        &s->cells_top[s->zoom_props->local_zoom_cells_with_particles_top[i]];
  for (int i = 0; i < nb; ++i)
    patches[nz + i].cell =
        &s->cells_top[s->zoom_props->local_bkg_cells_with_particles_top[i]];
  *patches_out = patches;
  return count;
}

static void zoom_mpi_exchange_counts(const size_t *send, size_t *recv,
                                     const int nr_nodes) {
  MPI_Alltoall((void *)send, sizeof(size_t), MPI_BYTE, recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
}

static void zoom_mpi_patches_to_density(struct zoom_mpi_patch *patches,
                                        const int nr_patches,
                                        const struct zoom_mpi_layout *layout,
                                        double *rho) {
  size_t *send_count =
      zoom_mpi_malloc("zoom_mpi_send_count", layout->nr_nodes, sizeof(size_t));
  size_t *recv_count =
      zoom_mpi_malloc("zoom_mpi_recv_count", layout->nr_nodes, sizeof(size_t));
  memset(send_count, 0, layout->nr_nodes * sizeof(size_t));
  size_t total = 0;
  for (int p = 0; p < nr_patches; ++p) {
    const struct zoom_mpi_patch *patch = &patches[p];
    if (patch->values == NULL) continue;
    for (int i = 0; i < patch->size[0]; ++i)
      for (int j = 0; j < patch->size[1]; ++j)
        for (int k = 0; k < patch->size[2]; ++k) {
          const double value =
              patch->values[zoom_mpi_patch_index(patch, i, j, k)];
          if (value == 0.) continue;
          ++send_count[zoom_mpi_owner(layout, i + patch->min[0])];
          ++total;
        }
  }
  size_t *cursor =
      zoom_mpi_malloc("zoom_mpi_cursor", layout->nr_nodes, sizeof(size_t));
  cursor[0] = 0;
  for (int rank = 1; rank < layout->nr_nodes; ++rank)
    cursor[rank] = cursor[rank - 1] + send_count[rank - 1];
  struct zoom_mpi_density_value *send_values = zoom_mpi_malloc(
      "zoom_mpi_density_send", total, sizeof(struct zoom_mpi_density_value));
  for (int p = 0; p < nr_patches; ++p) {
    struct zoom_mpi_patch *patch = &patches[p];
    if (patch->values == NULL) continue;
    for (int i = 0; i < patch->size[0]; ++i)
      for (int j = 0; j < patch->size[1]; ++j)
        for (int k = 0; k < patch->size[2]; ++k) {
          const size_t pi = zoom_mpi_patch_index(patch, i, j, k);
          const double value = patch->values[pi];
          if (value == 0.) continue;
          const int x = i + patch->min[0];
          const int y = j + patch->min[1];
          const int z = k + patch->min[2];
          const int owner = zoom_mpi_owner(layout, x);
          send_values[cursor[owner]].key = zoom_mpi_index(x, y, z, layout->M);
          send_values[cursor[owner]].value = value;
          ++cursor[owner];
        }
    swift_free("zoom_mpi_density_patch", patch->values);
    patch->values = NULL;
  }
  zoom_mpi_exchange_counts(send_count, recv_count, layout->nr_nodes);
  size_t recv_total = 0;
  for (int rank = 0; rank < layout->nr_nodes; ++rank)
    recv_total += recv_count[rank];
  struct zoom_mpi_density_value *recv_values =
      zoom_mpi_malloc("zoom_mpi_density_recv", recv_total,
                      sizeof(struct zoom_mpi_density_value));
  exchange_structs(send_count, send_values, recv_count, recv_values,
                   sizeof(struct zoom_mpi_density_value));
  for (size_t q = 0; q < recv_total; ++q) {
    const size_t key = recv_values[q].key;
    const int z = key % layout->M;
    const size_t xy = key / layout->M;
    const int y = xy % layout->M;
    const int x = xy / layout->M;
#ifdef SWIFT_DEBUG_CHECKS
    if (x < layout->local_0_start ||
        x >= layout->local_0_start + layout->local_n0)
      error("Received density outside local distributed zoom slab.");
#endif
    rho[zoom_mpi_real_local_index(x, y, z, layout)] += recv_values[q].value;
  }
  swift_free("zoom_mpi_density_send", send_values);
  swift_free("zoom_mpi_density_recv", recv_values);
  swift_free("zoom_mpi_cursor", cursor);
  swift_free("zoom_mpi_send_count", send_count);
  swift_free("zoom_mpi_recv_count", recv_count);
}

static void zoom_mpi_kernel_real(double *kernel,
                                 const struct zoom_pm_mesh *mesh,
                                 const struct space *s,
                                 const struct zoom_mpi_layout *layout) {
  const int N = mesh->N[0];
  const int M = layout->M;
  const double global_r_s_inv = s->periodic ? s->e->mesh->r_s_inv : 0.;
  for (ptrdiff_t il = 0; il < layout->local_n0; ++il) {
    const int i = (int)(il + layout->local_0_start);
    const int ii = i < N ? i : i - M;
    const double dx = ii / mesh->cell_fac[0];
    for (int j = 0; j < M; ++j) {
      const int jj = j < N ? j : j - M;
      const double dy = jj / mesh->cell_fac[1];
      for (int k = 0; k < M; ++k) {
        const int kk = k < N ? k : k - M;
        const double dz = kk / mesh->cell_fac[2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);
        const double zoom_long = r == 0. ? -M_2_SQRTPI * 0.5 * mesh->r_s_inv
                                         : -erf(0.5 * r * mesh->r_s_inv) / r;
        const double global_long =
            global_r_s_inv <= 0.
                ? 0.
                : (r == 0. ? -M_2_SQRTPI * 0.5 * global_r_s_inv
                           : -erf(0.5 * r * global_r_s_inv) / r);
        kernel[zoom_mpi_real_local_index(i, j, k, layout)] =
            zoom_long - global_long;
      }
    }
  }
}

static void zoom_mpi_deconvolve_mapper(void *map_data, int num, void *extra) {
  fftw_complex *values = map_data;
  const struct zoom_mpi_fourier_mapper_data *data = extra;
  const ptrdiff_t start =
      (values - data->field) / ((ptrdiff_t)data->M * data->complex_width);
  for (int local_y = (int)start; local_y < (int)start + num; ++local_y) {
    const int iy = (int)data->local_1_start + local_y;
    const int ky = iy > data->M / 2 ? iy - data->M : iy;
    for (int ix = 0; ix < data->M; ++ix) {
      const int kx = ix > data->M / 2 ? ix - data->M : ix;
      for (int kz = 0; kz < data->complex_width; ++kz) {
        const size_t q =
            ((size_t)local_y * data->M + ix) * data->complex_width + kz;
        const double correction = mesh_CIC_deconvolution(kx, ky, kz, data->M);
        data->field[q][0] *= correction;
        data->field[q][1] *= correction;
      }
    }
  }
}

static void zoom_mpi_multiply_mapper(void *map_data, int num, void *extra) {
  fftw_complex *field = map_data;
  const struct zoom_mpi_fourier_mapper_data *data = extra;
  const size_t offset = field - data->field;
  for (int i = 0; i < num; ++i) {
    const double a = field[i][0];
    const double b = field[i][1];
    const double c = data->kernel[offset + i][0];
    const double d = data->kernel[offset + i][1];
    field[i][0] = data->norm * (a * c - b * d);
    field[i][1] = data->norm * (a * d + b * c);
  }
}

static void zoom_mpi_build_kernel(struct zoom_pm_mesh *mesh,
                                  const struct space *s, struct threadpool *tp,
                                  const struct zoom_mpi_layout *layout) {
  const size_t alloc_count =
      layout->alloc_local > 0 ? (size_t)layout->alloc_local : 1;
  double *kernel = fftw_malloc(2 * alloc_count * sizeof(double));
  fftw_complex *kernel_fft = fftw_malloc(alloc_count * sizeof(fftw_complex));
  if (kernel == NULL || kernel_fft == NULL)
    error("Failed to allocate distributed zoom kernel.");
  memset(kernel, 0, 2 * alloc_count * sizeof(double));
  zoom_mpi_kernel_real(kernel, mesh, s, layout);
  fftw_plan plan = fftw_mpi_plan_dft_r2c_3d(
      layout->M, layout->M, layout->M, kernel, kernel_fft, MPI_COMM_WORLD,
      FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT | FFTW_DESTROY_INPUT);
  if (plan == NULL) error("Failed to plan distributed zoom kernel FFT.");
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(kernel);
  struct zoom_mpi_fourier_mapper_data data = {
      .field = kernel_fft,
      .local_1_start = layout->local_1_start,
      .M = layout->M,
      .complex_width = layout->complex_width};
  if (layout->local_n1 > 0)
    threadpool_map(tp, zoom_mpi_deconvolve_mapper, kernel_fft,
                   (int)layout->local_n1,
                   layout->M * layout->complex_width * sizeof(fftw_complex),
                   threadpool_auto_chunk_size, &data);
  mesh->kernel_fft = kernel_fft;
  memuse_log_allocation("fftw_zoom_mesh.kernel_fft", kernel_fft, 1,
                        alloc_count * sizeof(fftw_complex));
}

static void zoom_mpi_potential_bounds(const struct cell *c,
                                      const struct zoom_pm_mesh *mesh,
                                      int min[3], int max[3], int *found) {
  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL)
        zoom_mpi_potential_bounds(c->progeny[k], mesh, min, max, found);
    return;
  }
  if (!zoom_mesh_cell_is_covered(mesh, c, 0)) return;
  for (int n = 0; n < c->grav.count; ++n) {
    const struct gpart *gp = &c->grav.parts[n];
    if (gp->time_bin == time_bin_inhibited) continue;
    int ind[3];
    double d[3], t[3];
    if (!zoom_mpi_gpart_index(mesh, gp, ind, d, t, 1))
      error("Covered particle lacks distributed zoom force stencil.");
    for (int axis = 0; axis < 3; ++axis) {
      if (ind[axis] - 2 < min[axis]) min[axis] = ind[axis] - 2;
      if (ind[axis] + 3 > max[axis]) max[axis] = ind[axis] + 3;
    }
    *found = 1;
  }
}

static size_t zoom_mpi_prepare_potential_patches(
    struct zoom_mpi_patch *patches, const int nr_patches,
    const struct zoom_pm_mesh *mesh) {
  size_t total = 0;
  for (int p = 0; p < nr_patches; ++p) {
    struct zoom_mpi_patch *patch = &patches[p];
    int max[3] = {INT_MIN, INT_MIN, INT_MIN};
    patch->min[0] = patch->min[1] = patch->min[2] = INT_MAX;
    int found = 0;
    zoom_mpi_potential_bounds(patch->cell, mesh, patch->min, max, &found);
    if (!found) continue;
    size_t count = 1;
    for (int axis = 0; axis < 3; ++axis) {
      patch->size[axis] = max[axis] - patch->min[axis] + 1;
      if (patch->min[axis] < 0 || max[axis] >= mesh->N[axis])
        error("Distributed zoom potential patch exceeds active mesh.");
      count *= patch->size[axis];
    }
    patch->values =
        zoom_mpi_malloc("zoom_mpi_potential_patch", count, sizeof(double));
    total += count;
  }
  return total;
}

static void zoom_mpi_request_potential(struct zoom_mpi_patch *patches,
                                       const int nr_patches,
                                       const struct zoom_pm_mesh *mesh,
                                       const struct zoom_mpi_layout *layout,
                                       const double *potential) {
  const size_t total =
      zoom_mpi_prepare_potential_patches(patches, nr_patches, mesh);
  size_t *send_count = zoom_mpi_malloc("zoom_mpi_pot_send_count",
                                       layout->nr_nodes, sizeof(size_t));
  size_t *recv_count = zoom_mpi_malloc("zoom_mpi_pot_recv_count",
                                       layout->nr_nodes, sizeof(size_t));
  size_t *cursor =
      zoom_mpi_malloc("zoom_mpi_pot_cursor", layout->nr_nodes, sizeof(size_t));
  memset(send_count, 0, layout->nr_nodes * sizeof(size_t));
  for (int p = 0; p < nr_patches; ++p)
    if (patches[p].values != NULL)
      for (int i = 0; i < patches[p].size[0]; ++i) {
        const size_t yz = (size_t)patches[p].size[1] * patches[p].size[2];
        send_count[zoom_mpi_owner(layout, i + patches[p].min[0])] += yz;
      }
  cursor[0] = 0;
  for (int rank = 1; rank < layout->nr_nodes; ++rank)
    cursor[rank] = cursor[rank - 1] + send_count[rank - 1];
  struct zoom_mpi_potential_value *requests = zoom_mpi_malloc(
      "zoom_mpi_pot_requests", total, sizeof(struct zoom_mpi_potential_value));
  for (int p = 0; p < nr_patches; ++p) {
    const struct zoom_mpi_patch *patch = &patches[p];
    if (patch->values == NULL) continue;
    for (int i = 0; i < patch->size[0]; ++i)
      for (int j = 0; j < patch->size[1]; ++j)
        for (int k = 0; k < patch->size[2]; ++k) {
          const int x = i + patch->min[0];
          const int y = j + patch->min[1];
          const int z = k + patch->min[2];
          const int owner = zoom_mpi_owner(layout, x);
          struct zoom_mpi_potential_value *request = &requests[cursor[owner]++];
          request->key = zoom_mpi_index(x, y, z, layout->M);
          request->patch = p;
          request->patch_offset = zoom_mpi_patch_index(patch, i, j, k);
          request->value = 0.;
        }
  }
  zoom_mpi_exchange_counts(send_count, recv_count, layout->nr_nodes);
  size_t recv_total = 0;
  for (int rank = 0; rank < layout->nr_nodes; ++rank)
    recv_total += recv_count[rank];
  struct zoom_mpi_potential_value *incoming =
      zoom_mpi_malloc("zoom_mpi_pot_incoming", recv_total,
                      sizeof(struct zoom_mpi_potential_value));
  exchange_structs(send_count, requests, recv_count, incoming,
                   sizeof(struct zoom_mpi_potential_value));
  for (size_t q = 0; q < recv_total; ++q) {
    const size_t key = incoming[q].key;
    const int z = key % layout->M;
    const size_t xy = key / layout->M;
    const int y = xy % layout->M;
    const int x = xy / layout->M;
#ifdef SWIFT_DEBUG_CHECKS
    if (x < layout->local_0_start ||
        x >= layout->local_0_start + layout->local_n0)
      error("Potential request outside local distributed zoom slab.");
#endif
    incoming[q].value = potential[zoom_mpi_real_local_index(x, y, z, layout)];
  }
  exchange_structs(recv_count, incoming, send_count, requests,
                   sizeof(struct zoom_mpi_potential_value));
  for (size_t q = 0; q < total; ++q) {
    const int p = requests[q].patch;
    if (p < 0 || p >= nr_patches)
      error("Invalid distributed zoom potential patch response.");
    patches[p].values[requests[q].patch_offset] = requests[q].value;
  }
  swift_free("zoom_mpi_pot_requests", requests);
  swift_free("zoom_mpi_pot_incoming", incoming);
  swift_free("zoom_mpi_pot_cursor", cursor);
  swift_free("zoom_mpi_pot_send_count", send_count);
  swift_free("zoom_mpi_pot_recv_count", recv_count);
}

static double zoom_mpi_patch_CIC_get(const struct zoom_mpi_patch *patch,
                                     const int i, const int j, const int k,
                                     const double t[3], const double d[3]) {
#define ZOOM_MPI_GET(ii, jj, kk, w) \
  patch->values[zoom_mpi_patch_index(patch, ii, jj, kk)] * (w)
  double value = ZOOM_MPI_GET(i, j, k, t[0] * t[1] * t[2]);
  value += ZOOM_MPI_GET(i, j, k + 1, t[0] * t[1] * d[2]);
  value += ZOOM_MPI_GET(i, j + 1, k, t[0] * d[1] * t[2]);
  value += ZOOM_MPI_GET(i, j + 1, k + 1, t[0] * d[1] * d[2]);
  value += ZOOM_MPI_GET(i + 1, j, k, d[0] * t[1] * t[2]);
  value += ZOOM_MPI_GET(i + 1, j, k + 1, d[0] * t[1] * d[2]);
  value += ZOOM_MPI_GET(i + 1, j + 1, k, d[0] * d[1] * t[2]);
  value += ZOOM_MPI_GET(i + 1, j + 1, k + 1, d[0] * d[1] * d[2]);
#undef ZOOM_MPI_GET
  return value;
}

static void zoom_mpi_interpolate_leaf(
    struct cell *c, const struct zoom_mpi_patch *patch,
    const struct zoom_mpi_interpolate_mapper_data *data) {
  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL)
        zoom_mpi_interpolate_leaf(c->progeny[k], patch, data);
    return;
  }
  if (!zoom_mesh_cell_is_covered(data->mesh, c, 0)) return;
  for (int n = 0; n < c->grav.count; ++n) {
    struct gpart *gp = &c->grav.parts[n];
    if (gp->time_bin == time_bin_inhibited) continue;
    int ind[3];
    double d[3], t[3];
    if (!zoom_mpi_gpart_index(data->mesh, gp, ind, d, t, 1))
      error("Covered particle lacks distributed zoom force stencil.");
    const int i = ind[0] - patch->min[0];
    const int j = ind[1] - patch->min[1];
    const int k = ind[2] - patch->min[2];
    const double p = zoom_mpi_patch_CIC_get(patch, i, j, k, t, d);
    double a[3];
    a[0] = (1. / 12.) * zoom_mpi_patch_CIC_get(patch, i + 2, j, k, t, d) -
           (2. / 3.) * zoom_mpi_patch_CIC_get(patch, i + 1, j, k, t, d) +
           (2. / 3.) * zoom_mpi_patch_CIC_get(patch, i - 1, j, k, t, d) -
           (1. / 12.) * zoom_mpi_patch_CIC_get(patch, i - 2, j, k, t, d);
    a[1] = (1. / 12.) * zoom_mpi_patch_CIC_get(patch, i, j + 2, k, t, d) -
           (2. / 3.) * zoom_mpi_patch_CIC_get(patch, i, j + 1, k, t, d) +
           (2. / 3.) * zoom_mpi_patch_CIC_get(patch, i, j - 1, k, t, d) -
           (1. / 12.) * zoom_mpi_patch_CIC_get(patch, i, j - 2, k, t, d);
    a[2] = (1. / 12.) * zoom_mpi_patch_CIC_get(patch, i, j, k + 2, t, d) -
           (2. / 3.) * zoom_mpi_patch_CIC_get(patch, i, j, k + 1, t, d) +
           (2. / 3.) * zoom_mpi_patch_CIC_get(patch, i, j, k - 1, t, d) -
           (1. / 12.) * zoom_mpi_patch_CIC_get(patch, i, j, k - 2, t, d);
    gp->a_grav_mesh[0] += data->const_G * data->mesh->cell_fac[0] * a[0];
    gp->a_grav_mesh[1] += data->const_G * data->mesh->cell_fac[1] * a[1];
    gp->a_grav_mesh[2] += data->const_G * data->mesh->cell_fac[2] * a[2];
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
    gp->potential_mesh += data->const_G * p;
#endif
  }
}

static void zoom_mpi_interpolate_mapper(void *map_data, int num, void *extra) {
  struct zoom_mpi_patch *patches = map_data;
  const struct zoom_mpi_interpolate_mapper_data *data = extra;
  for (int p = 0; p < num; ++p)
    if (patches[p].values != NULL)
      zoom_mpi_interpolate_leaf(patches[p].cell, &patches[p], data);
}

void zoom_mesh_compute_potential_distributed(struct zoom_pm_mesh *mesh,
                                             const struct space *s,
                                             struct threadpool *tp,
                                             int verbose) {
  if (mesh == NULL || !mesh->enabled) return;
  if (mesh->N[0] != mesh->N[1] || mesh->N[0] != mesh->N[2])
    error("Distributed zoom mesh must be cubic.");

  const int M = 2 * mesh->N[0];
  struct zoom_mpi_layout layout;
  zoom_mpi_layout_init(&layout, M);
  if (verbose)
    message("Distributed zoom FFT slab: real x=[%td,%td), Fourier y=[%td,%td).",
            layout.local_0_start, layout.local_0_start + layout.local_n0,
            layout.local_1_start, layout.local_1_start + layout.local_n1);

  ticks tic = getticks();
  struct zoom_mpi_patch *patches;
  const int nr_patches = zoom_mpi_collect_cells(s, &patches);
  struct neutrino_model nu_model;
  memset(&nu_model, 0, sizeof(nu_model));
  if (s->e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(s, &nu_model);
  const struct zoom_mpi_patch_mapper_data patch_data = {mesh, &nu_model};
  if (nr_patches > 0)
    threadpool_map(tp, zoom_mpi_build_density_patches_mapper, patches,
                   nr_patches, sizeof(struct zoom_mpi_patch),
                   threadpool_auto_chunk_size, (void *)&patch_data);
  if (verbose)
    message("Distributed zoom private patch deposition took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  const size_t alloc_count =
      layout.alloc_local > 0 ? (size_t)layout.alloc_local : 1;
  double *rho = fftw_malloc(2 * alloc_count * sizeof(double));
  fftw_complex *frho = fftw_malloc(alloc_count * sizeof(fftw_complex));
  if (rho == NULL || frho == NULL)
    error("Failed to allocate distributed zoom FFT arrays.");
  memset(rho, 0, 2 * alloc_count * sizeof(double));

  tic = getticks();
  zoom_mpi_patches_to_density(patches, nr_patches, &layout, rho);
  if (verbose)
    message("Distributed zoom density exchange took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (mesh->kernel_fft == NULL) {
    tic = getticks();
    zoom_mpi_build_kernel(mesh, s, tp, &layout);
    if (verbose)
      message("Distributed zoom kernel build took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

  tic = getticks();
  fftw_plan forward = fftw_mpi_plan_dft_r2c_3d(
      M, M, M, rho, frho, MPI_COMM_WORLD,
      FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT | FFTW_DESTROY_INPUT);
  if (forward == NULL) error("Failed to plan distributed zoom forward FFT.");
  fftw_execute(forward);
  fftw_destroy_plan(forward);
  if (verbose)
    message("Distributed zoom forward FFT took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
  const size_t fourier_count =
      (size_t)layout.local_n1 * M * layout.complex_width;
  struct zoom_mpi_fourier_mapper_data multiply_data = {
      .field = frho,
      .kernel = (const fftw_complex *)mesh->kernel_fft,
      .M = M,
      .complex_width = layout.complex_width,
      .norm = 1. / ((double)M * M * M)};
  if (fourier_count > 0)
    threadpool_map(tp, zoom_mpi_multiply_mapper, frho, fourier_count,
                   sizeof(fftw_complex), threadpool_auto_chunk_size,
                   &multiply_data);
  if (s->e->neutrino_properties->use_linear_response)
    neutrino_response_compute_mesh(s, tp, frho, layout.local_1_start,
                                   layout.local_n1, M, 2. * mesh->dim[0],
                                   /*clear_dc=*/0, verbose);
  if (verbose)
    message("Distributed zoom Fourier response took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
  fftw_plan inverse = fftw_mpi_plan_dft_c2r_3d(
      M, M, M, frho, rho, MPI_COMM_WORLD,
      FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN | FFTW_DESTROY_INPUT);
  if (inverse == NULL) error("Failed to plan distributed zoom inverse FFT.");
  fftw_execute(inverse);
  fftw_destroy_plan(inverse);
  fftw_free(frho);
  if (verbose)
    message("Distributed zoom inverse FFT took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
  zoom_mpi_request_potential(patches, nr_patches, mesh, &layout, rho);
  fftw_free(rho);
  if (verbose)
    message("Distributed zoom potential patch exchange took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
  const struct zoom_mpi_interpolate_mapper_data interpolate_data = {
      mesh, s->e->physical_constants->const_newton_G};
  if (nr_patches > 0)
    threadpool_map(tp, zoom_mpi_interpolate_mapper, patches, nr_patches,
                   sizeof(struct zoom_mpi_patch), threadpool_auto_chunk_size,
                   (void *)&interpolate_data);
  for (int p = 0; p < nr_patches; ++p)
    if (patches[p].values != NULL)
      swift_free("zoom_mpi_potential_patch", patches[p].values);
  swift_free("zoom_mpi_patches", patches);
  zoom_mpi_layout_clean(&layout);
  if (verbose)
    message("Distributed zoom additive interpolation took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

#else

void zoom_mesh_compute_potential_distributed(struct zoom_pm_mesh *mesh,
                                             const struct space *s,
                                             struct threadpool *tp,
                                             int verbose) {
  (void)s;
  (void)tp;
  (void)verbose;
  if (mesh == NULL || !mesh->enabled) return;
  error("FFTW MPI not found - unable to use distributed zoom mesh");
}

#endif
