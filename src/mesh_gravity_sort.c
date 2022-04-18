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

/* This object's header. */
#include "mesh_gravity_sort.h"

/* Local includes. */
#include "align.h"
#include "atomic.h"
#include "error.h"
#include "row_major_id.h"
#include "threadpool.h"

struct mapper_extra_data {

  /* Mesh size */
  int N;

  /* Buckets */
  size_t *bucket_counts;
};

/**
 * @param Count how may mesh cells will end up in each x-coord bucket.
 */
void bucket_sort_mesh_key_value_rho_count_mapper(void *map_data, int nr_parts,
                                                 void *extra_data) {

  /* Unpack the data */
  const struct mesh_key_value_rho *array_in =
      (const struct mesh_key_value_rho *)map_data;
  struct mapper_extra_data *data = (struct mapper_extra_data *)extra_data;
  const int N = data->N;
  size_t *global_bucket_counts = data->bucket_counts;

  /* Local buckets */
  size_t *local_bucket_counts = (size_t *)calloc(N, sizeof(size_t));

  /* Count how many items will land in each bucket. */
  for (int i = 0; i < nr_parts; ++i) {

    const size_t key = array_in[i].key;

    /* Get the x coordinate of this mesh cell in the global mesh
     * Note: we don't need to sort more precisely than just
     * the x coordinate */
    const int mesh_x = get_xcoord_from_padded_row_major_id(key, N);

#ifdef SWIFT_DEBUG_CHECKS
    if (mesh_x < 0) error("Invalid mesh cell x-coordinate (too small)");
    if (mesh_x >= N) error("Invalid mesh cell x-coordinate (too large)");
#endif

    /* Add a contribution to the bucket count */
    local_bucket_counts[mesh_x]++;
  }

  /* Now write back to memory */
  for (int i = 0; i < N; ++i) {
    atomic_add(&global_bucket_counts[i], local_bucket_counts[i]);
  }

  /* Clean up */
  free(local_bucket_counts);
}

/**
 * @param Count how may mesh cells will end up in each x-coord bucket.
 */
void bucket_sort_mesh_key_value_pot_count_mapper(void *map_data, int nr_parts,
                                                 void *extra_data) {

  /* Unpack the data */
  const struct mesh_key_value_pot *array_in =
      (const struct mesh_key_value_pot *)map_data;
  struct mapper_extra_data *data = (struct mapper_extra_data *)extra_data;
  const int N = data->N;
  size_t *global_bucket_counts = data->bucket_counts;

  /* Local buckets */
  size_t *local_bucket_counts = (size_t *)calloc(N, sizeof(size_t));

  /* Count how many items will land in each bucket. */
  for (int i = 0; i < nr_parts; ++i) {

    const size_t key = array_in[i].key;

    /* Get the x coordinate of this mesh cell in the global mesh
     * Note: we don't need to sort more precisely than just
     * the x coordinate */
    const int mesh_x = get_xcoord_from_padded_row_major_id(key, N);

#ifdef SWIFT_DEBUG_CHECKS
    if (mesh_x < 0) error("Invalid mesh cell x-coordinate (too small)");
    if (mesh_x >= N) error("Invalid mesh cell x-coordinate (too large)");
#endif

    /* Add a contribution to the bucket count */
    local_bucket_counts[mesh_x]++;
  }

  /* Now write back to memory */
  for (int i = 0; i < N; ++i) {
    atomic_add(&global_bucket_counts[i], local_bucket_counts[i]);
  }

  /* Clean up */
  free(local_bucket_counts);
}

/**
 * @param Count how may mesh cells will end up in each cell index
 */
void bucket_sort_mesh_key_value_pot_index_count_mapper(void *map_data,
                                                       int nr_parts,
                                                       void *extra_data) {

  /* Unpack the data */
  const struct mesh_key_value_pot *array_in =
      (const struct mesh_key_value_pot *)map_data;
  struct mapper_extra_data *data = (struct mapper_extra_data *)extra_data;
  const int N = data->N;
  size_t *global_bucket_counts = data->bucket_counts;

  /* Local buckets */
  size_t *local_bucket_counts = (size_t *)calloc(N, sizeof(size_t));

  /* Count how many items will land in each bucket. */
  for (int i = 0; i < nr_parts; ++i) {

    const size_t key = array_in[i].cell_index;

    /* Get the cell_index */
    const int cell_id = cell_index_extract_patch_index(key);

#ifdef SWIFT_DEBUG_CHECKS
    if (cell_id < 0) error("Invalid mesh cell x-coordinate (too small)");
    if (cell_id >= N) error("Invalid mesh cell x-coordinate (too large)");
#endif

    /* Add a contribution to the bucket count */
    local_bucket_counts[cell_id]++;
  }

  /* Now write back to memory */
  for (int i = 0; i < N; ++i) {
    atomic_add(&global_bucket_counts[i], local_bucket_counts[i]);
  }

  /* Clean up */
  free(local_bucket_counts);
}

/**
 * @brief Bucket sort of the array of mesh cells based on their x-coord.
 *
 * Note the two mesh_key_value_rho arrays must be aligned on
 * SWIFT_CACHE_ALIGNMENT.
 *
 * @param array_in The unsorted array of mesh-key value pairs.
 * @param count The number of elements in the mesh-key value pair arrays.
 * @param N The size of the mesh. (i.e. the number of possible x-axis values)
 * @param tp The #threadpool object.
 * @param array_out The sorted array of mesh-key value pairs (to be filled).
 * @param bucket_offsets The offsets in the sorted array where we change x-coord
 * (to be filled).
 */
void bucket_sort_mesh_key_value_rho(const struct mesh_key_value_rho *array_in,
                                    const size_t count, const int N,
                                    struct threadpool *tp,
                                    struct mesh_key_value_rho *array_out,
                                    size_t *bucket_offsets) {

  /* Create an array of bucket counts and one of offsets */
  size_t *bucket_counts = (size_t *)malloc(N * sizeof(size_t));
  memset(bucket_counts, 0, N * sizeof(size_t));

  struct mapper_extra_data extra_data;
  extra_data.N = N;
  extra_data.bucket_counts = bucket_counts;

  /* Collect the number of items that will end up in each bucket */
  threadpool_map(tp, bucket_sort_mesh_key_value_rho_count_mapper,
                 (void *)array_in, count, sizeof(struct mesh_key_value_rho),
                 threadpool_auto_chunk_size, &extra_data);

#ifdef SWIFT_DEBUG_CHECKS
  size_t count_check = 0;
  for (int i = 0; i < N; ++i) count_check += bucket_counts[i];
  if (count_check != count) error("Bucket count is not matching");
#endif

  /* Now we can build the array of offsets (cumsum of the counts) */
  bucket_offsets[0] = 0;
  for (int i = 1; i < N; ++i) {
    bucket_offsets[i] = bucket_offsets[i - 1] + bucket_counts[i - 1];
  }

  /* Remind the compiler that the array is nicely aligned */
  swift_declare_aligned_ptr(struct mesh_key_value_rho, array_out_aligned,
                            array_out, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(const struct mesh_key_value_rho, array_in_aligned,
                            array_in, SWIFT_CACHE_ALIGNMENT);

  /* Now, we can do the actual sorting */
  for (size_t i = 0; i < count; ++i) {

    const size_t key = array_in_aligned[i].key;
    const int mesh_x = get_xcoord_from_padded_row_major_id(key, N);

    /* Where does this element land? */
    const size_t index = bucket_offsets[mesh_x];

    /* Copy the element to its correct position */
    memcpy(&array_out_aligned[index], &array_in_aligned[i],
           sizeof(struct mesh_key_value_rho));

    /* Move the start of this bucket by one */
    bucket_offsets[mesh_x]++;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that things have indeed been sorted */
  int last_mesh_x =
      get_xcoord_from_padded_row_major_id(array_out_aligned[0].key, N);
  ;
  for (size_t i = 1; i < count; ++i) {
    const size_t key = array_out_aligned[i].key;
    const int mesh_x = get_xcoord_from_padded_row_major_id(key, N);
    if (mesh_x < last_mesh_x) error("Unsorted array!");
    last_mesh_x = mesh_x;
  }
#endif

  /* Restore the bucket offsets to their former glory */
  for (int i = 0; i < N; ++i) {
    bucket_offsets[i] -= bucket_counts[i];
  }

  /* Clean up! */
  free(bucket_counts);
}

/**
 * @brief Bucket sort of the array of mesh cells based on their x-coord.
 *
 * Note the two mesh_key_value_pot arrays must be aligned on
 * SWIFT_CACHE_ALIGNMENT.
 *
 * @param array_in The unsorted array of mesh-key value pairs.
 * @param count The number of elements in the mesh-key value pair arrays.
 * @param N The size of the mesh. (i.e. the number of possible x-axis values)
 * @param tp The #threadpool object.
 * @param array_out The sorted array of mesh-key value pairs (to be filled).
 * @param bucket_offsets The offsets in the sorted array where we change x-coord
 * (to be filled).
 */
void bucket_sort_mesh_key_value_pot(const struct mesh_key_value_pot *array_in,
                                    const size_t count, const int N,
                                    struct threadpool *tp,
                                    struct mesh_key_value_pot *array_out,
                                    size_t *bucket_offsets) {

  /* Create an array of bucket counts and one of offsets */
  size_t *bucket_counts = (size_t *)malloc(N * sizeof(size_t));
  memset(bucket_counts, 0, N * sizeof(size_t));

  struct mapper_extra_data extra_data;
  extra_data.N = N;
  extra_data.bucket_counts = bucket_counts;

  /* Collect the number of items that will end up in each bucket */
  threadpool_map(tp, bucket_sort_mesh_key_value_pot_count_mapper,
                 (void *)array_in, count, sizeof(struct mesh_key_value_pot),
                 threadpool_auto_chunk_size, &extra_data);

#ifdef SWIFT_DEBUG_CHECKS
  size_t count_check = 0;
  for (int i = 0; i < N; ++i) count_check += bucket_counts[i];
  if (count_check != count) error("Bucket count is not matching");
#endif

  /* Now we can build the array of offsets (cumsum of the counts) */
  bucket_offsets[0] = 0;
  for (int i = 1; i < N; ++i) {
    bucket_offsets[i] = bucket_offsets[i - 1] + bucket_counts[i - 1];
  }

  /* Remind the compiler that the array is nicely aligned */
  swift_declare_aligned_ptr(struct mesh_key_value_pot, array_out_aligned,
                            array_out, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(const struct mesh_key_value_pot, array_in_aligned,
                            array_in, SWIFT_CACHE_ALIGNMENT);

  /* Now, we can do the actual sorting */
  for (size_t i = 0; i < count; ++i) {

    const size_t key = array_in_aligned[i].key;
    const int mesh_x = get_xcoord_from_padded_row_major_id(key, N);

    /* Where does this element land? */
    const size_t index = bucket_offsets[mesh_x];

    /* Copy the element to its correct position */
    memcpy(&array_out_aligned[index], &array_in_aligned[i],
           sizeof(struct mesh_key_value_pot));

    /* Move the start of this bucket by one */
    bucket_offsets[mesh_x]++;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that things have indeed been sorted */
  int last_mesh_x =
      get_xcoord_from_padded_row_major_id(array_out_aligned[0].key, N);

  for (size_t i = 1; i < count; ++i) {
    const size_t key = array_out_aligned[i].key;
    const int mesh_x = get_xcoord_from_padded_row_major_id(key, N);
    if (mesh_x < last_mesh_x) error("Unsorted array!");
    last_mesh_x = mesh_x;
  }
#endif

  /* Restore the bucket offsets to their former glory */
  for (int i = 0; i < N; ++i) {
    bucket_offsets[i] -= bucket_counts[i];
  }

  /* Clean up! */
  free(bucket_counts);
}

/**
 * @brief Bucket sort of the array of mesh cells based on their cell index.
 *
 * Note the two mesh_key_value_pot arrays must be aligned on
 * SWIFT_CACHE_ALIGNMENT.
 *
 * @param array_in The unsorted array of mesh-key value pairs.
 * @param count The number of elements in the mesh-key value pair arrays.
 * @param N The number of buckets (i.e. the nr. of local top-level cells).
 * @param tp The #threadpool object.
 * @param array_out The sorted array of mesh-key value pairs (to be filled).
 */
void bucket_sort_mesh_key_value_pot_index(
    const struct mesh_key_value_pot *array_in, const size_t count, const int N,
    struct threadpool *tp, struct mesh_key_value_pot *array_out) {

  /* Create an array of bucket counts and one of offsets */
  size_t *bucket_counts = (size_t *)malloc(N * sizeof(size_t));
  size_t *bucket_offsets = (size_t *)malloc(N * sizeof(size_t));
  memset(bucket_counts, 0, N * sizeof(size_t));

  struct mapper_extra_data extra_data;
  extra_data.N = N;
  extra_data.bucket_counts = bucket_counts;

  /* Collect the number of items that will end up in each bucket */
  threadpool_map(tp, bucket_sort_mesh_key_value_pot_index_count_mapper,
                 (void *)array_in, count, sizeof(struct mesh_key_value_pot),
                 threadpool_auto_chunk_size, &extra_data);

#ifdef SWIFT_DEBUG_CHECKS
  size_t count_check = 0;
  for (int i = 0; i < N; ++i) count_check += bucket_counts[i];
  if (count_check != count) error("Bucket count is not matching");
#endif

  /* Now we can build the array of offsets (cumsum of the counts) */
  bucket_offsets[0] = 0;
  for (int i = 1; i < N; ++i) {
    bucket_offsets[i] = bucket_offsets[i - 1] + bucket_counts[i - 1];
  }

  /* Remind the compiler that the array is nicely aligned */
  swift_declare_aligned_ptr(const struct mesh_key_value_pot, array_in_aligned,
                            array_in, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(struct mesh_key_value_pot, array_out_aligned,
                            array_out, SWIFT_CACHE_ALIGNMENT);

  /* Now, we can do the actual sorting */
  for (size_t i = 0; i < count; ++i) {

    const size_t key = array_in_aligned[i].cell_index;
    const int cell_id = cell_index_extract_patch_index(key);

    /* Where does this element land? */
    const size_t index = bucket_offsets[cell_id];

    /* Copy the element to its correct position */
    memcpy(&array_out_aligned[index], &array_in_aligned[i],
           sizeof(struct mesh_key_value_pot));

    /* Move the start of this bucket by one */
    bucket_offsets[cell_id]++;
  }

  /* Clean up! */
  free(bucket_offsets);
  free(bucket_counts);
}
