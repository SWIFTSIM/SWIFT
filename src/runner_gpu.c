/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Sarah Johnston (sarah.c.johnston@durham.ac.uk)
 *                    Will Roper (w.roper@sussex.ac.uk)
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

#include "active.h"
#include "engine.h"
#include "error.h"
#include "gpu_functions.h"
#include "gpu_mapping.h"
#include "runner.h"
#include "runner_doiact_grav.h"
#include "scheduler.h"
#include "timers.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

static int runner_gpu_local_ranks_on_device_for_budget = 1;

/* ------------------------------------------------------------------------- */
/* GPU debugging helpers: functions relating to the debugging checks for GPU */
/* ------------------------------------------------------------------------- */

/**
 * @brief Check the counters of the GPU queues to ensure no counting errors
 *
 * @param r The #runner
 * @param sched The scheduler
 * @param where	The string for the function that causes the error
 */
static void runner_gpu_check_queue_counters(struct runner *r,
                                            struct scheduler *sched,
                                            const char *where) {

  const int self_left = sched->queues[r->qid].gpu_self_tasks_left;
  const int pair_left = sched->queues[r->qid].gpu_pair_tasks_left;

  if (self_left < 0 || pair_left < 0 ||
      self_left > 1000000 || pair_left > 1000000) {
    error("%s: GPU queue counter corrupted: qid=%d "
          "self_left=%d pair_left=%d",
          where,
          r->qid,
          self_left,
          pair_left);
  }
}

/**
 * @brief Record debugging info when a GPU task is completed
 *
 * @param r The #runner
 * @param sched The scheduler
 * @param t The #task
 */
static void runner_gpu_count_self_task(struct runner *r,
                                       struct scheduler *sched,
                                       struct task *t) {

  if (t == NULL)
    error("runner_gpu_count_self_task got NULL task.");

  if (t->gpu_completed) {
    error("Packing already GPU-completed self task: task=%p type=%s subtype=%s "
          "gpu_counted=%d qid=%d self_left=%d",
          (void *)t,
          taskID_names[t->type],
          subtaskID_names[t->subtype],
          t->gpu_counted,
          r->qid,
          sched->queues[r->qid].gpu_self_tasks_left);
  }

  if (t->gpu_counted)
    return;

  t->gpu_counted = 1;

  lock_lock(&sched->queues[r->qid].lock);
  sched->queues[r->qid].gpu_self_tasks_left++;
  runner_gpu_check_queue_counters(r, sched, "runner_gpu_count_self_task");
  (void)lock_unlock(&sched->queues[r->qid].lock);
}

/**
 * @brief Count the number of GPU pair tasks
 *
 * @param r The #runner
 * @param sched The scheduler
 * @param t The #task
 */
static void runner_gpu_count_pair_task(struct runner *r,
                                       struct scheduler *sched,
                                       struct task *t) {

  if (t == NULL)
    error("runner_gpu_count_pair_task got NULL task.");

  if (t->gpu_counted)
    return;

  t->gpu_counted = 1;

  lock_lock(&sched->queues[r->qid].lock);
  sched->queues[r->qid].gpu_pair_tasks_left++;
  runner_gpu_check_queue_counters(r, sched, "runner_gpu_count_pair_task");
  (void)lock_unlock(&sched->queues[r->qid].lock);
}

/**
 * @brief Bind the GPU to the MPI rank
 *
 * @param r The #runner
 */
void runner_gpu_bind_device(struct runner *r) {

#if defined(WITH_CUDA) || defined(WITH_HIP)
  const int device_id = r->gpu.device_id;

  if (device_id < 0)
    error("runner_gpu_bind_device: invalid GPU device_id=%d", device_id);

  const GPUError err = GPUSetDevice(device_id);

  if (err != GPU_SUCCESS)
    error("runner_gpu_bind_device: GPUSetDevice(%d) failed: %s",
          device_id, GPUGetErrorString(err));
#else
  (void)r;
#endif
}

/**
 * @brief GPU error check
 *
 * @param where	The string for the function that causes the error
 */
static inline void runner_gpu_check_error(const char *where) {
  const GPUError err = GPUGetLastError();

  if (err != GPU_SUCCESS)
    error("%s: %s", where, GPUGetErrorString(err));
}


/* ------------------------------------------------------------------------- */
/* GPU timing helpers: thread-local state only, no struct layout changes.     */
/* This avoids changing runner_gpu.h and avoids shared locks in runner threads.*/
/* ------------------------------------------------------------------------- */

#ifdef SWIFT_DEBUG_TASKS
/**
 * @brief Return the elapsed time between two GPU events in seconds.
 *
 * Computes the elapsed time between the supplied GPU start and stop events
 * using GPUEventElapsedTime(). The GPU runtime reports the elapsed time in
 * milliseconds, which is converted to seconds before being returned.
 *
 * @param start The GPU event marking the start of the timed interval.
 * @param stop The GPU event marking the end of the timed interval.
 *
 * @return The elapsed time between @p start and @p stop in seconds.
 */
static double runner_gpu_event_offset_s(
    const GPUEvent start,
    const GPUEvent stop) {

  float elapsed_ms = 0.0f;

  const GPUError err =
      GPUEventElapsedTime(&elapsed_ms, start, stop);

  if (err != GPU_SUCCESS) {
    error("GPUEventElapsedTime failed: %s",
          GPUGetErrorString(err));
  }

  return 1.0e-3 * (double)elapsed_ms;
}

/**
 * @brief Append a GPU timeline entry to the timeline CSV file.
 *
 * Writes timing information for a GPU task to the file specified by the
 * SWIFT_GPU_TIMELINE_FILE environment variable. If this variable is not set,
 * the default file "gpu_timeline.csv" is used.
 *
 * If the file is empty, a CSV header is written before the first data row.
 * The function silently returns if the output file cannot be opened.
 *
 * @param kind String identifying the type of GPU task, for example "self" or "pair".
 * @param step Simulation timestep associated with the task.
 * @param runner_id ID of the runner executing the task.
 * @param substream_id ID of the GPU substream used for the task.
 * @param anchor_tic CPU tick value used as the reference point for the timeline measurements.
 * @param sync_toc CPU tick value recorded after synchronization.
 * @param h2d_end_s Time, in seconds relative to the timeline anchor, at which the host-to-device transfer completed.
 * @param kernel_start_s Time, in seconds relative to the timeline anchor, at which the GPU kernel started.
 * @param kernel_end_s Time, in seconds relative to the timeline anchor, at which the GPU kernel completed.
 * @param d2h_start_s Time, in seconds relative to the timeline anchor, at which the device-to-host transfer started.
 * @param d2h_end_s Time, in seconds relative to the timeline anchor, at which the device-to-host transfer completed.
 */
static void runner_gpu_write_timeline_row(
    const char *kind,
    const long long step,
    const int runner_id,
    const int substream_id,
    const ticks anchor_tic,
    const ticks sync_toc,
    const double h2d_end_s,
    const double kernel_start_s,
    const double kernel_end_s,
    const double d2h_start_s,
    const double d2h_end_s) {

  const char *path = getenv("SWIFT_GPU_TIMELINE_FILE");

  if (path == NULL || path[0] == '\0')
    path = "gpu_timeline.csv";

  FILE *fp = fopen(path, "a");

  if (fp == NULL)
    return;

  fseek(fp, 0, SEEK_END);

  if (ftell(fp) == 0) {
    fprintf(
        fp,
        "kind,step,runner,substream,anchor_tic,sync_toc,"
        "h2d_end_s,kernel_start_s,kernel_end_s,"
        "d2h_start_s,d2h_end_s\n");
  }

  fprintf(
      fp,
      "%s,%lld,%d,%d,%llu,%llu,"
      "%.9e,%.9e,%.9e,%.9e,%.9e\n",
      kind,
      step,
      runner_id,
      substream_id,
      (unsigned long long)anchor_tic,
      (unsigned long long)sync_toc,
      h2d_end_s,
      kernel_start_s,
      kernel_end_s,
      d2h_start_s,
      d2h_end_s);

  fclose(fp);
}

#endif /* SWIFT_DEBUG_TASKS */

#ifdef SWIFT_GPU_TIMING 
static __thread double runner_gpu_self_pack_time_s = 0.0;
static __thread double runner_gpu_pair_pack_time_s = 0.0;

/**
 * @brief Convert GPU walltime to seconds
 */
static inline double runner_gpu_walltime_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

/**
 * @brief Return the elapsed time between two GPU events in seconds.
 *
 * Synchronizes on the stop event to ensure that the timed GPU work has
 * completed, then computes the elapsed time between the start and stop
 * events. The GPU runtime reports the elapsed time in milliseconds, which
 * is converted to seconds before being returned.
 *
 * @param start The GPU event marking the start of the timed interval.
 * @param stop The GPU event marking the end of the timed interval.
 *
 * @return The elapsed time between @p start and @p stop in seconds.
 */
static double runner_gpu_event_elapsed_s(GPUEvent start, GPUEvent stop) {
  float elapsed_ms = 0.0f;

  GPUEventSynchronize(stop);
  GPUEventElapsedTime(&elapsed_ms, start, stop);

  return 1.0e-3 * (double)elapsed_ms;
}

/**
 * @brief Return the total host-to-device data size for a self-gravity batch.
 *
 * Computes the number of bytes transferred from host to device for the
 * packed self-gravity particle data, active-particle indices, and per-slot
 * metadata associated with a GPU substream.
 *
 * @param substream The GPU runner substream containing the packed self-gravity
 *                  particle and active-particle counts.
 * @param nslots The number of self-gravity slots in the batch.
 *
 * @return The total host-to-device transfer size, in bytes.
 */
static inline size_t runner_gpu_self_h2d_bytes(
    const struct gpu_runner_substream *substream,
    const int nslots) {

  const size_t particle_bytes =
      (size_t)substream->self_total_count *
      (sizeof(float4) + sizeof(float));

  const size_t active_index_bytes =
      (size_t)substream->self_total_active_count * sizeof(int);

  /*
   * Per-slot integer arrays:
   *   self_counts
   *   self_offsets
   *   self_active_counts
   *   self_active_offsets
   *   self_cell_flags
   *   self_use_full
   */
  const size_t slot_metadata_bytes =
      (size_t)nslots * 6u * sizeof(int);

  return particle_bytes +
         active_index_bytes +
         slot_metadata_bytes;
}

/**
 * @brief Return the total host-to-device data size for a pair-gravity batch.
 *
 * Computes the number of bytes transferred from host to device for the
 * packed pair-gravity particle data, active-particle indices, per-slot
 * metadata, and per-pair metadata associated with a GPU substream.
 *
 * @param substream The GPU runner substream containing the packed pair-gravity
 *                  particle and active-particle counts.
 * @param nslots The number of cell slots in the pair-gravity batch.
 * @param npairs The number of cell pairs in the batch.
 *
 * @return The total host-to-device transfer size, in bytes.
 */
static inline size_t runner_gpu_pair_h2d_bytes(
    const struct gpu_runner_substream *substream,
    const int nslots,
    const int npairs) {

  const size_t particle_bytes =
      (size_t)substream->pair_total_count *
      (sizeof(float4) + sizeof(float));

  const size_t active_index_bytes =
      (size_t)substream->pair_total_active_count * sizeof(int);

  /*
   * Per-slot integer arrays:
   *   pair_counts
   *   pair_offsets
   *   pair_active_counts
   *   pair_active_offsets
   *   pair_cell_flags
   */
  const size_t slot_metadata_bytes =
      (size_t)nslots * 5u * sizeof(int);

  /*
   * Per-pair integer arrays:
   *   pair_pair_i
   *   pair_pair_j
   *   pair_use_full
   *   pair_side_active_offsets, with two entries per pair
   *
   * Total: 1 + 1 + 1 + 2 = 5 integers per pair.
   */
  const size_t pair_metadata_bytes =
      (size_t)npairs * 5u * sizeof(int);

  return particle_bytes +
         active_index_bytes +
         slot_metadata_bytes +
         pair_metadata_bytes;
}

/**
 * @brief Append a GPU batch timing entry to the timing CSV file.
 *
 * Writes timing and workload information for a GPU batch to the file
 * specified by the SWIFT_GPU_TIMING_FILE environment variable. If this
 * variable is not set, the default file "gpu_timing_ga008.csv" is used.
 *
 * If the file is empty, a CSV header is written before the first data row.
 * The function silently returns if the output file cannot be opened.
 *
 * @param kind String identifying the type of GPU task, for example "self"
 *             or "pair".
 * @param step Simulation timestep associated with the batch.
 * @param runner_id ID of the runner executing the batch.
 * @param stream_id ID of the GPU stream used for the batch.
 * @param nslots Number of cell slots contained in the batch.
 * @param nparts Number of particles contained in the batch.
 * @param h2d_bytes Number of bytes transferred from host to device.
 * @param pack_s Time spent packing the batch on the host, in seconds.
 * @param h2d_s Time spent transferring data from host to device, in seconds.
 * @param kernel_s Time spent executing the GPU kernel, in seconds.
 * @param d2h_s Time spent transferring data from device to host, in seconds.
 * @param unpack_s Time spent unpacking the batch on the host, in seconds.
 */
static void runner_gpu_write_timing_row(
    const char *kind,
    long long step,
    int runner_id,
    int stream_id,
    int nslots,
    int nparts,
    size_t h2d_bytes,
    double pack_s,
    double h2d_s,
    double kernel_s,
    double d2h_s,
    double unpack_s) {

	#ifdef SWIFT_USE_TIMERS
	const ticks io_tic = getticks();
	#endif

  const char *path = getenv("SWIFT_GPU_TIMING_FILE");

  if (path == NULL || path[0] == '\0')
    path = "gpu_timing_ga008.csv";

  FILE *fp = fopen(path, "a");

  if (fp == NULL)
    return;

  fseek(fp, 0, SEEK_END);

  if (ftell(fp) == 0) {
  fprintf(fp,
          "kind,step,runner,stream,nslots,nparticles,h2d_bytes,"
          "pack_s,h2d_s,kernel_s,d2h_s,unpack_s\n");
	}

  fprintf(fp,
        "%s,%lld,%d,%d,%d,%d,%zu,"
        "%.9e,%.9e,%.9e,%.9e,%.9e\n",
        kind,
        step,
        runner_id,
        stream_id,
        nslots,
        nparts,
        h2d_bytes,
        pack_s,
        h2d_s,
        kernel_s,
        d2h_s,
        unpack_s);

  fclose(fp);

  #ifdef SWIFT_USE_TIMERS
  timers_toc(timer_gpu_grav_timing_io, io_tic);
  #endif
}
#endif

/**
 * @brief Launch the GPU P-P gravity kernel for a packed pair batch.
 *
 * @param periodic Whether periodic boundary conditions are enabled.
 * @param min_trunc Minimum truncation radius for periodic forces.
 * @param r_s_inv Pointer to the inverse splitting scale for periodic mesh
 *                forces.
 * @param pair_use_full_d Device array indicating whether each pair uses the
 *                        full or truncated gravity interaction.
 * @param pair_side_active_offsets_d Device array containing the receive-buffer
 *                                   offsets for each side of every pair.
 * @param pair_counts_d Device array containing the particle count for each
 *                      unique cell.
 * @param pair_offsets_d Device array containing the packed particle offset for
 *                       each unique cell.
 * @param pair_active_counts_d Device array containing the number of active
 *                             particles in each unique cell.
 * @param pair_active_offsets_d Device array containing the offset into the
 *                              active-particle index array for each unique
 *                              cell.
 * @param pair_active_index_d Device array containing the local indices of
 *                            active particles.
 * @param pair_pair_i_d Device array mapping each pair to its first unique-cell
 *                      slot.
 * @param pair_pair_j_d Device array mapping each pair to its second unique-cell
 *                      slot.
 * @param npairs Number of cell pairs in the batch.
 * @param nslots Number of unique cell slots in the batch.
 * @param dim_0 Domain size in the x direction.
 * @param dim_1 Domain size in the y direction.
 * @param dim_2 Domain size in the z direction.
 * @param pair_cell_flags_d Device array indicating which unique cells are
 *                          active and local to this rank.
 * @param send_pair_pos_mass_d Device buffer containing packed particle
 *                             positions and masses.
 * @param send_pair_h_d Device buffer containing packed particle softenings.
 * @param gravity_gpu_values_recv_d Device buffer receiving the calculated
 *                                  gravity results.
 * @param ncells Number of packed cell slots supplied to the kernel.
 * @param max_cell_size Maximum number of particles permitted in a packed cell.
 * @param max_active_count Maximum number of active particles in any packed
 *                         cell.
 * @param stream GPU stream used for the kernel launch.
 */
extern void pair_pp_offload_gpu(
    int periodic, double min_trunc, const float *r_s_inv,
    const int *pair_use_full_d,
    const int *pair_side_active_offsets_d,
    const int *pair_counts_d,
    const int *pair_offsets_d,
    const int *pair_active_counts_d,
    const int *pair_active_offsets_d,
    const int *pair_active_index_d,
    const int *pair_pair_i_d,
    const int *pair_pair_j_d,
    int npairs,
    int nslots,
    float dim_0, float dim_1, float dim_2,
    const int *pair_cell_flags_d,
    const float4 *send_pair_pos_mass_d,
    const float *send_pair_h_d,
    struct gravity_gpu_values_recv *gravity_gpu_values_recv_d,
    int ncells,
    int max_cell_size,
    int max_active_count,
    GPUStream stream);
    
/**
 * @brief Launch the GPU P-P self-gravity kernel for a packed batch.
 *
 * @param periodic Whether periodic boundary conditions are enabled.
 * @param r_s_inv Inverse splitting scale for periodic mesh forces.
 * @param self_cell_flags_d Device array indicating which cells are active.
 * @param self_use_full_d Device array indicating whether each cell uses the
 *                        full or truncated gravity interaction.
 * @param counts_d Device array containing the particle count for each cell.
 * @param offsets_d Device array containing the particle offset for each cell.
 * @param active_counts_d Device array containing the number of active
 *                        particles in each cell.
 * @param active_offsets_d Device array containing the offset into the active
 *                         particle index array for each cell.
 * @param active_index_d Device array containing the local indices of active
 *                       particles.
 * @param send_self_pos_mass_d Device buffer containing packed particle
 *                             positions and masses.
 * @param send_self_h_d Device buffer containing packed particle softenings.
 * @param recv_d Device buffer receiving the calculated gravity results.
 * @param ncells Number of packed cells in the batch.
 * @param max_cell_size Maximum number of particles in a packed cell.
 * @param max_active_count Maximum number of active particles in any cell.
 * @param stream GPU stream used for the kernel launch.
 */    
extern void self_pp_offload_gpu(
    int periodic,
    const float *r_s_inv,
    const int *self_cell_flags_d,
    const int *self_use_full_d,
    const int *counts_d,
    const int *offsets_d,
    const int *active_counts_d,
    const int *active_offsets_d,
    const int *active_index_d,
    const float4 *send_self_pos_mass_d,
    const float *send_self_h_d,
    struct gravity_gpu_values_recv *recv_d,
    int ncells,
    int max_cell_size,
    int max_active_count,
    GPUStream stream);
/**
 * @brief Unpack the GPU parameters from the parameter file.
 *
 * @param e The #engine to unpack the parameters for.
 */
void runner_gpu_params_init(struct engine *e) {

  /* Unpack the number of cells we will pack onto the GPU at a time. */
  e->ncells_per_gpu_grav_pack = parser_get_opt_param_int(
      e->parameter_file, "GPU:ncells_per_gpu_grav_pack", -1);

  if (e->ncells_per_gpu_grav_pack == 0 ||
      e->ncells_per_gpu_grav_pack < -1) {
    error("GPU:ncells_per_gpu_grav_pack must be >= 1 for user selected or -1 for auto");
  }
}

/**
 * @brief Store the values for packing the self cells without gaps.
 *
 * @param substream The stream the task is on.
 * @param slot The number cell in the pack it is
 * @param count The size of the cell.
 */
static inline void append_packed_self_cell(
    struct gpu_runner_substream *substream, int slot, int count) {

  substream->self_offsets_h[slot] = substream->self_total_count;
  substream->self_counts_h[slot] = count;
  substream->self_total_count += count;
}

/**
 * @brief Store the values for packing the pair cells without gaps.
 *
 * @param substream The stream the task is on.
 * @param slot The number cell in the pack it is
 * @param count The size of the cell.
 */
static inline void append_packed_pair_cell(
    struct gpu_runner_substream *substream, int slot, int count) {

  substream->pair_offsets_h[slot] = substream->pair_total_count;
  substream->pair_counts_h[slot] = count;
  substream->pair_total_count += count;
}

/**
 * @brief Find an existing pair-gravity cell slot or pack a new cell.
 *
 * Searches the current GPU pair-gravity batch for the supplied cell. If the
 * cell has already been packed, its existing slot index is returned. Otherwise,
 * a new slot is allocated and the cell's gravity particle data are copied from
 * the gravity cache into the packed host buffers.
 *
 * The function also constructs the list of active particle indices for the
 * cell, updates the per-slot and total active-particle counts, and records the
 * cell in the unique-cell list. An error is raised if the number of unique
 * cells exceeds the available batch slots or if the cell contains more
 * particles than the configured maximum cell size.
 *
 * @param r The runner processing the GPU gravity batch.
 * @param substream The GPU runner substream containing the pair-gravity
 *                  packing buffers and metadata.
 * @param c The cell to find or pack.
 * @param cache The gravity cache containing the particle data for @p c.
 * @param max_cell_size Maximum number of gravity particles allowed in a
 *                      packed cell.
 *
 * @return The slot index corresponding to @p c in the pair-gravity batch.
 */
static int runner_gpu_find_or_pack_pair_cell(
    struct runner *r,
    struct gpu_runner_substream *substream,
    struct cell *c,
    struct gravity_cache *cache,
    int max_cell_size) {

  for (int s = 0; s < substream->pair_unique_cell_count; s++) {
    if (substream->pair_unique_cells[s] == c)
      return s;
  }

  const int slot = substream->pair_unique_cell_count++;
  const int gcount = c->grav.count;

  if (slot >= r->gpu.grav_batch_ncells)
    error("Too many unique pair cells in GPU batch");

  if (gcount > max_cell_size)
    error("Pair unique-cell pack overflow: gcount=%d > max_cell_size=%d",
          gcount, max_cell_size);

  append_packed_pair_cell(substream, slot, gcount);

  const int off = substream->pair_offsets_h[slot];

  for (int i = 0; i < gcount; i++) {
    const int k = off + i;

    substream->send_pair_pos_mass[k].x = cache->x[i];
    substream->send_pair_pos_mass[k].y = cache->y[i];
    substream->send_pair_pos_mass[k].z = cache->z[i];
    substream->send_pair_pos_mass[k].w = cache->m[i];

    substream->send_pair_h[k] = cache->epsilon[i];
  }

  int active_count = 0;
  const int active_base = substream->pair_total_active_count;
  substream->pair_active_offsets_h[slot] = active_base;

  const int local_active_cell =
    (c->nodeID == r->e->nodeID) && cell_is_active_gravity(c, r->e);

  for (int i = 0; i < gcount; i++) {
    if (local_active_cell && cache->active[i] > 0) {
      substream->pair_active_index_h[active_base + active_count] = i;
      active_count++;
    }
  }

  substream->pair_total_active_count += active_count;
  substream->pair_active_counts_h[slot] = active_count;

  if (active_count > substream->pair_max_active_count)
    substream->pair_max_active_count = active_count;

  substream->pair_unique_cells[slot] = c;

  return slot;
}


/**
 * @brief Mark a packed self-gravity task as complete on the scheduler.
 *
 * @param r The #runner owning the task.
 * @param sched The scheduler tracking the task.
 * @param t The task to complete.
 */
static void runner_gpu_complete_self_task(struct runner *r,
                                          struct scheduler *sched,
                                          struct task *t) {

  if (t == NULL)
    error("runner_gpu_complete_self_task got NULL task.");

  if (t->gpu_completed) {
	  error("runner_gpu_complete_self_task called for already-completed task: "
		"task=%p type=%s subtype=%s gpu_counted=%d"
		"self_left=%d qid=%d",
		(void *)t,
		taskID_names[t->type],
		subtaskID_names[t->subtype],
		t->gpu_counted,
		sched->queues[r->qid].gpu_self_tasks_left,
		r->qid);
	}

  t->gpu_completed = 1;

  if (!t->gpu_counted) {
	  error("Completing uncounted GPU self task: task=%p type=%s subtype=%s "
		"gpu_completed=%d qid=%d",
		(void *)t,
		taskID_names[t->type],
		subtaskID_names[t->subtype],
		t->gpu_completed,
		r->qid);
	}

	lock_lock(&sched->queues[r->qid].lock);

	if (sched->queues[r->qid].gpu_self_tasks_left <= 0)
	  error("gpu_self_tasks_left underflow: task=%p type=%s subtype=%s qid=%d",
		(void *)t,
		taskID_names[t->type],
		subtaskID_names[t->subtype],
		r->qid);

	sched->queues[r->qid].gpu_self_tasks_left--;
	runner_gpu_check_queue_counters(r, sched, "runner_gpu_complete_self_task");

	(void)lock_unlock(&sched->queues[r->qid].lock);

	t->gpu_counted = 0;

  scheduler_done(sched, t);
}

/**
 * @brief Wrapper to call runner_gpu_complete_self_task
 *
 * @param r The #runner owning the task.
 * @param sched The scheduler tracking the task.
 * @param t The task to complete.
 */
void runner_gpu_complete_current_self_task(struct runner *r,
                                           struct scheduler *sched,
                                           struct task *t) {
  runner_gpu_complete_self_task(r, sched, t);
}

/**
 * @brief Mark a packed pair-gravity task as complete on the scheduler.
 *
 * @param r The #runner owning the task.
 * @param sched The scheduler tracking the task.
 * @param t The task to complete.
 */
void runner_gpu_complete_pair_task(struct runner *r, struct scheduler *sched,
                                   struct task *t) {

  if (t == NULL)
    error("runner_gpu_complete_pair_task got NULL task.");

  if (t->gpu_completed) {
    error("runner_gpu_complete_pair_task called for already-completed task: "
          "task=%p type=%s subtype=%s gpu_counted=%d "
          "pair_left=%d qid=%d",
          (void *)t,
          taskID_names[t->type],
          subtaskID_names[t->subtype],
          t->gpu_counted,
          sched->queues[r->qid].gpu_pair_tasks_left,
          r->qid);
  }

  t->gpu_completed = 1;

  if (!t->gpu_counted) {
    error("Completing uncounted GPU pair task: task=%p type=%s subtype=%s "
          "gpu_completed=%d qid=%d pair_left=%d",
          (void *)t,
          taskID_names[t->type],
          subtaskID_names[t->subtype],
          t->gpu_completed,
          r->qid,
          sched->queues[r->qid].gpu_pair_tasks_left);
  }

  lock_lock(&sched->queues[r->qid].lock);

  const int before = sched->queues[r->qid].gpu_pair_tasks_left;

  if (before <= 0)
    error("gpu_pair_tasks_left underflow before completing task=%p "
          "type=%s subtype=%s qid=%d pair_left=%d",
          (void *)t,
          taskID_names[t->type],
          subtaskID_names[t->subtype],
          r->qid,
          before);

  sched->queues[r->qid].gpu_pair_tasks_left--;
  runner_gpu_check_queue_counters(r, sched, "runner_gpu_complete_pair_task");

  (void)lock_unlock(&sched->queues[r->qid].lock);

  t->gpu_counted = 0;

  scheduler_done(sched, t);
}

/**
 * @brief Complete all eligible self-gravity tasks in the current GPU batch.
 *
 * Marks each unique counted scheduler task in the batch as complete, excluding
 * the currently executing task when one is supplied. The function then clears
 * the batch task and cell pointers and resets the self-gravity packing
 * metadata for the substream.
 *
 * @param r The #runner owning the batch.
 * @param sched The scheduler tracking the tasks.
 * @param substream The GPU substream containing the completed self-gravity
 *                  batch.
 * @param current_task The task currently being processed by the runner, or
 *                     NULL for a leftover batch flush. This task is not
 *                     completed here.
 */
void runner_gpu_complete_self_batch(struct runner *r, struct scheduler *sched,
                                    struct gpu_runner_substream *substream,
                                    struct task *current_task) {

  const int count = substream->grav_batch_self_count;

  /* Complete each top-level self task at most once. Recursive self walks can
     pack many leaf cells for the same scheduler task, so duplicate task
     pointers in grav_tasks_self[] are normal. */
  for (int i = 0; i < count; i++) {

    struct task *task = substream->grav_tasks_self[i];

    if (task == NULL)
      continue;

    /* The currently-walking task is completed by runner_main when
       flushed_self_task is returned. */
    if (task == current_task)
      continue;

    /* Duplicate slot for a task already completed earlier in this batch. */
    if (task->gpu_completed)
      continue;

    /* Only counted top-level GPU tasks own a queue counter. */
    if (!task->gpu_counted)
      continue;

    runner_gpu_complete_self_task(r, sched, task);
  }

  /* Clear self batch cell/task entries. */
  for (int i = 0; i < count; i++) {
    substream->grav_cells_self[i] = NULL;
    substream->grav_tasks_self[i] = NULL;
  }

  substream->grav_batch_self_count = 0;
  substream->busy = 0;

  for (int i = 0; i < count; i++) {
	  substream->self_counts_h[i] = 0;
	  substream->self_offsets_h[i] = 0;
	  substream->self_active_counts_h[i] = 0;
	  substream->self_active_offsets_h[i] = 0;
	  substream->self_rmax_h[i] = 0.f;

	  substream->self_cell_flags_h[i] = 0;
	  substream->self_use_full_h[i] = 0;
	}

  for (int a = 0; a < substream->self_total_active_count; a++)
    substream->self_active_index_h[a] = 0;

  substream->self_total_count = 0;
  substream->self_max_active_count = 0;
  substream->self_total_active_count = 0;

  runner_gpu_check_queue_counters(r, sched,
                                  "runner_gpu_complete_self_batch:end");
}

/**
 * @brief Pack one leaf pair-gravity interaction into a GPU batch.
 *
 * Populates the gravity caches for both cells, packs any unique cells not
 * already present in the batch, records the pair-to-cell mapping, active-cell
 * metadata, and periodic-force mode, and associates the interaction with its
 * parent scheduler task. No GPU work is launched by this function.
 *
 * @param r The #runner processing the interaction.
 * @param substream The GPU substream into which the interaction is packed.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param symmetric Whether both cells are updated.
 * @param allow_mpole Whether multipole information may be used when populating
 *                    the gravity caches.
 * @param grav_cells_pair Array storing the cell pointers for each packed pair.
 * @param grav_tasks_pair Array storing the scheduler task associated with each
 *                        packed pair.
 * @param grav_pair_internal_from_self Array marking pair interactions generated
 *                                     internally during a self-gravity walk.
 * @param t The top-level #task currently being processed.
 * @param internal_from_self Whether this pair was generated internally from a
 *                           self-gravity task.
 * @param max_cell_size Maximum number of gravity particles permitted in a
 *                      packed cell.
 * @param stream GPU stream associated with the batch.
 */
static void runner_dopair_grav_pp_pack(
    struct runner *r, struct gpu_runner_substream *substream,
    struct cell *ci, struct cell *cj, const int symmetric,
    const int allow_mpole,
    struct cell **grav_cells_pair, struct task **grav_tasks_pair,
    unsigned char *grav_pair_internal_from_self,
    struct task *t, int internal_from_self,
    int max_cell_size, GPUStream stream) {

	#ifdef SWIFT_USE_TIMERS
	const ticks pair_pack_tic = getticks();
	#endif

  (void)symmetric;
  (void)stream;

  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const float dim[3] = {(float)e->mesh->dim[0], (float)e->mesh->dim[1],
                        (float)e->mesh->dim[2]};
  const double min_trunc = e->mesh->r_cut_min;

  #ifdef SWIFT_GPU_TIMING 
  const double pack_time_start = runner_gpu_walltime_s();
  #endif

  const int ci_active =
      cell_is_active_gravity(ci, e) && (ci->nodeID == e->nodeID);
  const int cj_active =
      cell_is_active_gravity(cj, e) && (cj->nodeID == e->nodeID);

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->split || cj->split) error("Running P-P on splitable cells");
  if (!cell_are_gpart_drifted(ci, e)) error("Un-drifted gparts");
  if (!cell_are_gpart_drifted(cj, e)) error("Un-drifted gparts");

  if (cj_active && ci->grav.ti_old_multipole != e->ti_current)
    error("Un-drifted multipole");
  if (ci_active && cj->grav.ti_old_multipole != e->ti_current)
    error("Un-drifted multipole");
#endif

  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;
  struct gravity_cache *const cj_cache = &r->cj_gravity_cache;

  const double shift_i[3] = {0., 0., 0.};
  const double shift_j[3] = {0., 0., 0.};

  const float rmax_i = ci->grav.multipole->r_max;
  const float rmax_j = cj->grav.multipole->r_max;

  const float CoM_i[3] = {
      (float)(ci->grav.multipole->CoM[0] - shift_i[0]),
      (float)(ci->grav.multipole->CoM[1] - shift_i[1]),
      (float)(ci->grav.multipole->CoM[2] - shift_i[2])};

  const float CoM_j[3] = {
      (float)(cj->grav.multipole->CoM[0] - shift_j[0]),
      (float)(cj->grav.multipole->CoM[1] - shift_j[1]),
      (float)(cj->grav.multipole->CoM[2] - shift_j[2])};

  const int gcount_i = ci->grav.count;
  const int gcount_j = cj->grav.count;

  const int gcount_padded_i = gcount_i - (gcount_i % VEC_SIZE) + VEC_SIZE;
  const int gcount_padded_j = gcount_j - (gcount_j % VEC_SIZE) + VEC_SIZE;

  const int allow_multipole_i = allow_mpole && ci->grav.count > 1;
  const int allow_multipole_j = allow_mpole && cj->grav.count > 1;

  if (gcount_i > max_cell_size)
    error("Pair pack overflow: gcount_i=%d > max_cell_size=%d",
          gcount_i, max_cell_size);

  if (gcount_j > max_cell_size)
    error("Pair pack overflow: gcount_j=%d > max_cell_size=%d",
          gcount_j, max_cell_size);

  if (ci->nodeID == e->nodeID) {
    gravity_cache_populate(e->max_active_bin, allow_multipole_j, periodic, dim,
                           ci_cache, ci->grav.parts, gcount_i, gcount_padded_i,
                           shift_i, CoM_j, cj->grav.multipole, ci,
                           e->gravity_properties);
  } else {
    gravity_cache_populate_foreign(periodic, dim, ci_cache,
                                   ci->grav.parts_foreign, gcount_i,
                                   gcount_padded_i, shift_i, ci,
                                   e->gravity_properties);
  }

  if (cj->nodeID == e->nodeID) {
    gravity_cache_populate(e->max_active_bin, allow_multipole_i, periodic, dim,
                           cj_cache, cj->grav.parts, gcount_j, gcount_padded_j,
                           shift_j, CoM_i, ci->grav.multipole, cj,
                           e->gravity_properties);
  } else {
    gravity_cache_populate_foreign(periodic, dim, cj_cache,
                                   cj->grav.parts_foreign, gcount_j,
                                   gcount_padded_j, shift_j, cj,
                                   e->gravity_properties);
  }

  struct cell *a = ci;
  struct cell *b = cj;

  if (a > b) {
    struct cell *tmp = a;
    a = b;
    b = tmp;
  }

  while (cell_glocktree(a)) {
    ;
  }
  while (cell_glocktree(b)) {
    ;
  }
  
  const int pair_capacity = r->gpu.grav_batch_ncells / 2;
	const int cell_capacity = r->gpu.grav_batch_ncells;

	if (substream->grav_batch_pair_count < 0 ||
	    substream->grav_batch_pair_count >= pair_capacity) {
	  error("Pair batch overflow before pack: pair_count=%d pair_capacity=%d",
		substream->grav_batch_pair_count,
		pair_capacity);
	}

	if (substream->pair_unique_cell_count < 0 ||
	    substream->pair_unique_cell_count + 2 > cell_capacity) {
	  error("Pair unique-cell overflow before pack: unique=%d cell_capacity=%d",
		substream->pair_unique_cell_count,
		cell_capacity);
	}

  const int slot_i =
    runner_gpu_find_or_pack_pair_cell(r, substream, ci, ci_cache, max_cell_size);

  const int slot_j =
    runner_gpu_find_or_pack_pair_cell(r, substream, cj, cj_cache, max_cell_size);

  const int pair_id = substream->grav_batch_pair_count;
  
  if (pair_id < 0 || pair_id >= pair_capacity)
  	error("Bad pair_id=%d pair_capacity=%d", pair_id, pair_capacity);

  substream->pair_pair_i_h[pair_id] = slot_i;
  substream->pair_pair_j_h[pair_id] = slot_j;

  grav_cells_pair[2 * pair_id] = ci;
  grav_cells_pair[2 * pair_id + 1] = cj;
  grav_tasks_pair[pair_id] = t;
  grav_pair_internal_from_self[pair_id] =
      (unsigned char)internal_from_self;
      
  if (!internal_from_self)
  	runner_gpu_count_pair_task(r, &r->e->sched, t);

  int use_full = 1;

  if (periodic) {
    double d0 = CoM_j[0] - CoM_i[0];
    double d1 = CoM_j[1] - CoM_i[1];
    double d2 = CoM_j[2] - CoM_i[2];

    d0 = nearest(d0, e->mesh->dim[0]);
    d1 = nearest(d1, e->mesh->dim[1]);
    d2 = nearest(d2, e->mesh->dim[2]);

    const double r2 = d0 * d0 + d1 * d1 + d2 * d2;
    const double max_r = sqrt(r2) + rmax_i + rmax_j;

    use_full = (max_r <= min_trunc);
  }
  
  substream->pair_use_full_h[pair_id] = use_full;

  substream->pair_cell_flags_h[slot_i] =
    (ci->nodeID == e->nodeID && cell_is_active_gravity(ci, e)) ? 1 : 0;

  substream->pair_cell_flags_h[slot_j] =
    (cj->nodeID == e->nodeID && cell_is_active_gravity(cj, e)) ? 1 : 0;
    
  const int side_i = 2 * pair_id;
  const int side_j = side_i + 1;

  substream->pair_side_active_offsets_h[side_i] =
    substream->pair_total_pair_active_count;

  substream->pair_total_pair_active_count +=
    substream->pair_active_counts_h[slot_i];
    
  if (2 * pair_id + 1 >= 2 * pair_capacity)
	  error("Bad pair side offset index: pair_id=%d pair_capacity=%d",
		pair_id, pair_capacity);

  substream->pair_side_active_offsets_h[side_j] =
    substream->pair_total_pair_active_count;

  substream->pair_total_pair_active_count +=
    substream->pair_active_counts_h[slot_j];

  substream->grav_batch_pair_count++;

  gravity_cache_zero_output(ci_cache, gcount_padded_i);
  gravity_cache_zero_output(cj_cache, gcount_padded_j);

  cell_gunlocktree(b);
  cell_gunlocktree(a);

  #ifdef SWIFT_GPU_TIMING 
  runner_gpu_pair_pack_time_s += runner_gpu_walltime_s() - pack_time_start;
  #endif

  #ifdef SWIFT_USE_TIMERS
  timers_toc(timer_gpu_grav_pair_pack_cpu, pair_pack_tic);
  #endif
}

/**
 * @brief Unpack the GPU gravity results for one side of a pair interaction.
 *
 * Adds the accelerations and potentials returned by the GPU to the active
 * particles belonging to the specified cell. Foreign or inactive cells are
 * not updated.
 *
 * @param r The #runner processing the GPU batch.
 * @param substream The GPU substream containing the returned pair results and
 *                  active-particle metadata.
 * @param c The cell whose results are to be unpacked.
 * @param slot The unique-cell slot associated with @p c.
 * @param recv_base Offset into the packed receive buffer for this side of the
 *                  pair interaction.
 */
static inline void runner_gpu_unpack_pair_side(
    struct runner *r,
    struct gpu_runner_substream *substream,
    struct cell *c,
    int slot,
    int recv_base) {

  const struct engine *e = r->e;

  if (c == NULL)
    error("GPU pair unpack received NULL cell.");

  /* MPI safety: never write to a foreign cell. */
  if (c->nodeID != e->nodeID)
    return;

  if (!cell_is_active_gravity(c, e))
    return;

  const int active_count = substream->pair_active_counts_h[slot];
  const int active_base = substream->pair_active_offsets_h[slot];

  if (active_count == 0)
    return;

  while (cell_glocktree(c)) {
    ;
  }

  for (int a = 0; a < active_count; a++) {
    const int local_pid =
        substream->pair_active_index_h[active_base + a];

    const int k = recv_base + a;

#ifdef SWIFT_DEBUG_CHECKS
    if (local_pid < 0 || local_pid >= c->grav.count)
      error("GPU pair unpack local_pid=%d out of range [0,%d).",
            local_pid, c->grav.count);
#endif

    c->grav.parts[local_pid].a_grav[0] +=
        substream->recv_pair_active[k].values_i.x;
    c->grav.parts[local_pid].a_grav[1] +=
        substream->recv_pair_active[k].values_i.y;
    c->grav.parts[local_pid].a_grav[2] +=
        substream->recv_pair_active[k].values_i.z;
    c->grav.parts[local_pid].potential +=
        substream->recv_pair_active[k].values_i.w;
  }

  cell_gunlocktree(c);
}

/**
 * @brief Flush a full pair-gravity GPU batch: H2D copy, kernel launch, D2H
 *        copy, stream synchronisation, result unpacking, scheduler completion,
 *        and metadata reset.
 *
 * current_task is skipped because it may still be in the recursive walk.
 * If this flush happened while walking current_task, runner_dopair_grav_pp_gpu()
 * must return flushed_pair_task so runner_main() completes current_task.
 * If this is a leftover flush, current_task should be NULL and all non-internal
 * tasks in the batch are completed here.
 */
static void runner_dopair_grav_pp_flush(
    struct runner *r,
    struct gpu_runner_substream *substream,
    struct cell **grav_cells_pair,
    struct task **grav_tasks_pair,
    struct task *current_task,
    int ncells,
    int max_cell_size,
    GPUStream stream) {

  runner_gpu_bind_device(r);

  const int npairs = substream->grav_batch_pair_count;
  const int nslots = substream->pair_unique_cell_count;
  const int ncells_flush = nslots;

  if (npairs == 0 || ncells_flush == 0)
    return;

  #ifdef SWIFT_USE_TIMERS
  const ticks pair_flush_tic = getticks();
  #endif
    
  #ifdef SWIFT_DEBUG_TASKS

  /* Host-clock anchors used to align the GPU timeline with the
   * normal SWIFT task timeline. */
  ticks gpu_anchor_tic = 0;
  ticks gpu_sync_toc = 0;

  /* Each runner owns its own array of GPU substreams. */
  const int gpu_substream_id =
      (int)(substream - r->gpu.substreams);

#endif
    
  #ifdef SWIFT_DEBUG_TASKS
  /*
   * Record that an actual pair GPU batch flush occurred while this
   * scheduler task was being processed.
   *
   * current_task is NULL for leftover/end-of-queue flushes, so those
   * are deliberately not attributed to a scheduler task.
   */
  if (current_task != NULL) {
    current_task->gpu_debug_pair_flushes++;
  }
#endif

  if (npairs < 0 || npairs > ncells / 2)
    error("Bad pair flush npairs=%d capacity=%d", npairs, ncells / 2);

  if (nslots < 0 || nslots > ncells)
    error("Bad pair flush nslots=%d ncells=%d", nslots, ncells);

  if (substream->pair_total_count < 0)
    error("Bad pair_total_count=%d", substream->pair_total_count);

  if (substream->pair_total_active_count < 0)
    error("Bad pair_total_active_count=%d",
          substream->pair_total_active_count);

  if (substream->pair_total_pair_active_count < 0)
    error("Bad pair_total_pair_active_count=%d",
          substream->pair_total_pair_active_count);

  #ifdef SWIFT_GPU_TIMING 
  double h2d_s = 0.0;
  double kernel_s = 0.0;
  double d2h_s = 0.0;
  double unpack_s = 0.0;
  #endif

  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const float r_s_inv = e->mesh->r_s_inv;
  const double min_trunc = e->mesh->r_cut_min;

  const float dim_0 = (float)e->mesh->dim[0];
  const float dim_1 = (float)e->mesh->dim[1];
  const float dim_2 = (float)e->mesh->dim[2];

  #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)

GPUEvent h2d_start, h2d_stop;
GPUEvent kernel_start, kernel_stop;
GPUEvent d2h_start, d2h_stop;

GPUEventCreate(&h2d_start);
GPUEventCreate(&h2d_stop);

GPUEventCreate(&kernel_start);
GPUEventCreate(&kernel_stop);

GPUEventCreate(&d2h_start);
GPUEventCreate(&d2h_stop);

#endif

  /* ---- H2D copies ---- */
  {
  #ifdef SWIFT_DEBUG_TASKS
  gpu_anchor_tic = getticks();
#endif

#if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
GPUEventRecord(h2d_start, stream);
#endif

    GPUMemcpyAsync(
        substream->send_pair_pos_mass_d,
        substream->send_pair_pos_mass,
        (size_t)substream->pair_total_count * sizeof(float4),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->send_pair_h_d,
        substream->send_pair_h,
        (size_t)substream->pair_total_count * sizeof(float),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_counts_d,
        substream->pair_counts_h,
        (size_t)nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_offsets_d,
        substream->pair_offsets_h,
        (size_t)nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_active_counts_d,
        substream->pair_active_counts_h,
        (size_t)nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_active_offsets_d,
        substream->pair_active_offsets_h,
        (size_t)nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_active_index_d,
        substream->pair_active_index_h,
        (size_t)substream->pair_total_active_count * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_pair_i_d,
        substream->pair_pair_i_h,
        (size_t)npairs * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_pair_j_d,
        substream->pair_pair_j_h,
        (size_t)npairs * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_use_full_d,
        substream->pair_use_full_h,
        (size_t)npairs * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_side_active_offsets_d,
        substream->pair_side_active_offsets_h,
        (size_t)(2 * npairs) * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);

    GPUMemcpyAsync(
        substream->pair_cell_flags_d,
        substream->pair_cell_flags_h,
        (size_t)nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        stream);
        
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(h2d_stop, stream);
    #endif
  }

  runner_gpu_check_error("runner_dopair_grav_pp_flush H2D");

  /* ---- Kernel launch ---- */
  {
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(kernel_start, stream);
    #endif

    pair_pp_offload_gpu(
        periodic,
        min_trunc,
        &r_s_inv,
        substream->pair_use_full_d,
        substream->pair_side_active_offsets_d,
        substream->pair_counts_d,
        substream->pair_offsets_d,
        substream->pair_active_counts_d,
        substream->pair_active_offsets_d,
        substream->pair_active_index_d,
        substream->pair_pair_i_d,
        substream->pair_pair_j_d,
        npairs,
        nslots,
        dim_0,
        dim_1,
        dim_2,
        substream->pair_cell_flags_d,
        substream->send_pair_pos_mass_d,
        substream->send_pair_h_d,
        substream->recv_pair_active_d,
        ncells_flush,
        max_cell_size,
        substream->pair_max_active_count,
        stream);
    
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(kernel_stop, stream);
    #endif
  }

  runner_gpu_check_error("runner_dopair_grav_pp_flush kernel");

  /* ---- D2H copy ---- */
  {
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(d2h_start, stream);
    #endif

    const size_t recv_pair_active_capacity =
        (size_t)(2 * (r->gpu.grav_batch_ncells / 2)) *
        (size_t)r->gpu.grav_max_cell_size;

    if ((size_t)substream->pair_total_pair_active_count >
        recv_pair_active_capacity) {
      error("GPU pair recv overflow: pair_total_pair_active_count=%d "
            "capacity=%zu",
            substream->pair_total_pair_active_count,
            recv_pair_active_capacity);
    }

    GPUMemcpyAsync(
        substream->recv_pair_active,
        substream->recv_pair_active_d,
        (size_t)substream->pair_total_pair_active_count *
            sizeof(struct gravity_gpu_values_recv),
        GPU_MEMCPY_DEVICE_TO_HOST,
        stream);
        
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(d2h_stop, stream);
    #endif 
    GPUEventRecord(substream->done, stream);
    #ifdef SWIFT_USE_TIMERS
    const ticks gpu_wait_tic = getticks();
    #endif

    GPUEventSynchronize(substream->done);

    #ifdef SWIFT_USE_TIMERS
    timers_toc(timer_gpu_grav_wait, gpu_wait_tic);
    #endif
    
    #ifdef SWIFT_DEBUG_TASKS
  gpu_sync_toc = getticks();
#endif
    
    #ifdef SWIFT_DEBUG_TASKS

const double gpu_h2d_end_s =
    runner_gpu_event_offset_s(h2d_start, h2d_stop);

const double gpu_kernel_start_s =
    runner_gpu_event_offset_s(h2d_start, kernel_start);

const double gpu_kernel_end_s =
    runner_gpu_event_offset_s(h2d_start, kernel_stop);

const double gpu_d2h_start_s =
    runner_gpu_event_offset_s(h2d_start, d2h_start);

const double gpu_d2h_end_s =
    runner_gpu_event_offset_s(h2d_start, d2h_stop);

runner_gpu_write_timeline_row(
    "pair",
    (long long)e->step,
    r->id,
    gpu_substream_id,
    gpu_anchor_tic,
    gpu_sync_toc,
    gpu_h2d_end_s,
    gpu_kernel_start_s,
    gpu_kernel_end_s,
    gpu_d2h_start_s,
    gpu_d2h_end_s);

#endif
    
    #ifdef SWIFT_GPU_TIMING 
    h2d_s = runner_gpu_event_elapsed_s(h2d_start, h2d_stop);
    kernel_s = runner_gpu_event_elapsed_s(kernel_start, kernel_stop);
    d2h_s = runner_gpu_event_elapsed_s(d2h_start, d2h_stop);
    #endif 
  }

  runner_gpu_check_error("runner_dopair_grav_pp_flush D2H");

  /* ---- Unpack results back to particles ---- */
  {
    #ifdef SWIFT_GPU_TIMING 
    const double unpack_t0 = runner_gpu_walltime_s();
    #endif

    TIMER_TIC;

    for (int pair_id = 0; pair_id < npairs; pair_id++) {

      const int slot_i = substream->pair_pair_i_h[pair_id];
      const int slot_j = substream->pair_pair_j_h[pair_id];

      if (slot_i < 0 || slot_i >= nslots)
        error("Bad pair slot_i=%d nslots=%d pair_id=%d",
              slot_i,
              nslots,
              pair_id);

      if (slot_j < 0 || slot_j >= nslots)
        error("Bad pair slot_j=%d nslots=%d pair_id=%d",
              slot_j,
              nslots,
              pair_id);

      const int recv_i = substream->pair_side_active_offsets_h[2 * pair_id];
      const int recv_j =
          substream->pair_side_active_offsets_h[2 * pair_id + 1];

      if (recv_i < 0 ||
          recv_i + substream->pair_active_counts_h[slot_i] >
              substream->pair_total_pair_active_count) {
        error("Bad pair recv_i=%d active_i=%d total_pair_active=%d "
              "pair_id=%d",
              recv_i,
              substream->pair_active_counts_h[slot_i],
              substream->pair_total_pair_active_count,
              pair_id);
      }

      if (recv_j < 0 ||
          recv_j + substream->pair_active_counts_h[slot_j] >
              substream->pair_total_pair_active_count) {
        error("Bad pair recv_j=%d active_j=%d total_pair_active=%d "
              "pair_id=%d",
              recv_j,
              substream->pair_active_counts_h[slot_j],
              substream->pair_total_pair_active_count,
              pair_id);
      }

      runner_gpu_unpack_pair_side(
          r,
          substream,
          grav_cells_pair[2 * pair_id],
          slot_i,
          recv_i);

      runner_gpu_unpack_pair_side(
          r,
          substream,
          grav_cells_pair[2 * pair_id + 1],
          slot_j,
          recv_j);
          
      #ifdef SWIFT_DEBUG_CHECKS
/* Account for direct particle-particle interactions processed on the GPU. */
struct cell *ci_pair = grav_cells_pair[2 * pair_id];
struct cell *cj_pair = grav_cells_pair[2 * pair_id + 1];

const int count_i = ci_pair->grav.count;
const int count_j = cj_pair->grav.count;

/* Only local, active cells receive force contributions. */
if (ci_pair->nodeID == e->nodeID &&
    cell_is_active_gravity(ci_pair, e)) {

  for (int i = 0; i < count_i; i++) {
    accumulate_add_ll(
        &ci_pair->grav.parts[i].num_interacted, count_j);
  }
}

if (cj_pair->nodeID == e->nodeID &&
    cell_is_active_gravity(cj_pair, e)) {

  for (int j = 0; j < count_j; j++) {
    accumulate_add_ll(
        &cj_pair->grav.parts[j].num_interacted, count_i);
  }
}
#endif
    }
    

    TIMER_TOC(timer_dopair_grav_pp);

    #ifdef SWIFT_GPU_TIMING
    unpack_s = runner_gpu_walltime_s() - unpack_t0;
    #endif
  }

  /* ---- Complete scheduler tasks before clearing task pointers ---- */
  {
    struct scheduler *sched = &r->e->sched;

    int printed_current_skip = 0;

	for (int pair_id = 0; pair_id < npairs; pair_id++) {

	  struct task *batch_task = grav_tasks_pair[pair_id];
	  const int internal = substream->grav_pair_internal_from_self[pair_id];

	  if (batch_task == NULL)
	    error("NULL task in GPU pair batch: pair_id=%d npairs=%d",
		  pair_id, npairs);

	  if (!internal && batch_task == current_task) {
	    if (!printed_current_skip) {
	      printed_current_skip = 1;
	    }
	    continue;
	  }

	  if (!internal && !batch_task->gpu_completed && batch_task->gpu_counted)
	    runner_gpu_complete_pair_task(r, sched, batch_task);
	
    }
  }
  
  #ifdef SWIFT_GPU_TIMING

  const double pair_pack_s = runner_gpu_pair_pack_time_s;
  runner_gpu_pair_pack_time_s = 0.0;

  const size_t h2d_bytes =
    runner_gpu_pair_h2d_bytes(
        substream,
        nslots,
        npairs);

/* ---- Timing output ---- */
runner_gpu_write_timing_row(
    "pair",
    (long long)e->step,
    r->id,
    0,
    nslots,
    substream->pair_total_count,
    h2d_bytes,
    pair_pack_s,
    h2d_s,
    kernel_s,
    d2h_s,
    unpack_s);

#endif
      
  #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
  GPUEventDestroy(h2d_start);
  GPUEventDestroy(h2d_stop);
  GPUEventDestroy(kernel_start);
  GPUEventDestroy(kernel_stop);
  GPUEventDestroy(d2h_start);
  GPUEventDestroy(d2h_stop);
  #endif 

  /* ---- Reset pair entries: arrays indexed by pair_id ---- */
  for (int pair_id = 0; pair_id < npairs; pair_id++) {
    grav_cells_pair[2 * pair_id] = NULL;
    grav_cells_pair[2 * pair_id + 1] = NULL;
    grav_tasks_pair[pair_id] = NULL;

    substream->grav_pair_internal_from_self[pair_id] = 0;
    substream->pair_pair_i_h[pair_id] = 0;
    substream->pair_pair_j_h[pair_id] = 0;
    substream->pair_use_full_h[pair_id] = 0;

    substream->pair_side_active_offsets_h[2 * pair_id] = 0;
    substream->pair_side_active_offsets_h[2 * pair_id + 1] = 0;
  }

  /* ---- Reset unique-cell entries: arrays indexed by unique cell slot ---- */
  for (int slot = 0; slot < nslots; slot++) {
    substream->pair_counts_h[slot] = 0;
    substream->pair_offsets_h[slot] = 0;
    substream->pair_active_counts_h[slot] = 0;
    substream->pair_active_offsets_h[slot] = 0;
    substream->pair_cell_flags_h[slot] = 0;
    substream->pair_unique_cells[slot] = NULL;
  }

  /* ---- Reset active particle index buffer only for the span used ---- */
  for (int a = 0; a < substream->pair_total_active_count; a++) {
    substream->pair_active_index_h[a] = 0;
  }

  substream->grav_batch_pair_count = 0;
  substream->busy = 0;

  substream->pair_total_count = 0;
  substream->pair_unique_cell_count = 0;
  substream->pair_total_active_count = 0;
  substream->pair_max_active_count = 0;
  substream->pair_total_pair_active_count = 0;

  runner_gpu_check_queue_counters(
      r,
      &r->e->sched,
      "runner_dopair_grav_pp_flush:end");

  #ifdef SWIFT_USE_TIMERS
  timers_toc(timer_gpu_grav_pair_flush_cpu, pair_flush_tic);
  #endif
}

/**
 * @brief Pack a leaf pair-gravity interaction and flush if the batch is full.
 *
 * This is the main entry point called from runner_dopair_recursive_grav_gpu()
 * when two unsplit leaf cells are reached. It packs the pair, and if the batch
 * is now full it flushes the GPU work before returning.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 * @param symmetric Are we updating both cells (1) or just ci (0) ?
 * @param allow_mpole Are we allowing the use of M2P interactions ?
 */
enum runner_gpu_task_type runner_dopair_grav_pp_gpu(
    struct runner *r, struct gpu_runner_substream *substream, struct cell *ci,
    struct cell *cj, const int symmetric, const int allow_mpole,
    struct cell **grav_cells_pair, struct task **grav_tasks_pair,
    unsigned char *grav_pair_internal_from_self,
    struct task *t, int internal_from_self,
    int ncells, int max_cell_size, GPUStream stream) {

  const int pair_capacity = ncells / 2;
  enum runner_gpu_task_type result = packed_task;

  /* If the existing batch cannot accept this pair, flush it first.
     This flush skips current_task=t, so runner_main must complete t. */
  if (substream->grav_batch_pair_count + 1 > pair_capacity ||
      substream->pair_unique_cell_count + 2 > ncells) {

    runner_dopair_grav_pp_flush(
        r, substream,
        grav_cells_pair, grav_tasks_pair,
        t, ncells, max_cell_size, stream);

    result = flushed_pair_task;
  }

  runner_dopair_grav_pp_pack(
      r, substream, ci, cj, symmetric, allow_mpole,
      grav_cells_pair, grav_tasks_pair,
      grav_pair_internal_from_self,
      t, internal_from_self,
      max_cell_size, stream);

  /* If packing this pair filled the batch, flush it too.
     Again, current_task=t is skipped in the flush, so return flushed_pair_task. */
  if (substream->grav_batch_pair_count >= pair_capacity ||
      substream->pair_unique_cell_count >= ncells) {

    runner_dopair_grav_pp_flush(
        r, substream,
        grav_cells_pair, grav_tasks_pair,
        t, ncells, max_cell_size, stream);

    return flushed_pair_task;
  }

  return result;
}

/**
 * @brief Launch the GPU P-P self-gravity kernel for a packed self batch.
 *
 * Passes the packed self-gravity data and associated metadata for the current
 * substream to the GPU kernel.
 *
 * @param r The #runner processing the batch.
 * @param substream The GPU substream containing the packed self-gravity data.
 * @param nslots Number of packed cells in the batch.
 * @param max_active_count Maximum number of active particles in any packed
 *                         cell.
 * @param max_cell_size Maximum number of particles permitted in a packed cell.
 * @param stream GPU stream used for the kernel launch.
 */
static void runner_doself_grav_pp_flush(
    struct runner *r,
    struct gpu_runner_substream *substream,
    int nslots,
    int max_active_count,
    int max_cell_size,
    GPUStream stream) {

  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const float r_s_inv = e->mesh->r_s_inv;
  const double min_trunc = e->mesh->r_cut_min;

  self_pp_offload_gpu(
    periodic,
    &r_s_inv,
    substream->self_cell_flags_d,
    substream->self_use_full_d,
    substream->self_counts_d,
    substream->self_offsets_d,
    substream->self_active_counts_d,
    substream->self_active_offsets_d,
    substream->self_active_index_d,
    substream->send_self_pos_mass_d,
    substream->send_self_h_d,
    substream->recv_self_active_d,
    nslots,
    max_cell_size,
    max_active_count,
    stream);
}

/**
 * @brief Pack a self-gravity interaction and flush the GPU batch when full.
 *
 * Packs one self-gravity cell into the current GPU substream. If the batch
 * reaches its configured capacity, the function copies the packed data to the
 * device, launches the self-gravity kernel, copies the results back, unpacks
 * them into the particles, completes the corresponding scheduler tasks, and
 * resets the batch state.
 *
 * @param r The #runner processing the task.
 * @param substream The GPU substream used for packing and executing the batch.
 * @param ci The #cell to pack.
 * @param t The top-level #task being executed.
 * @param ncells Maximum number of cells in the GPU batch.
 * @param max_cell_size Maximum number of gravity particles permitted in a
 *                      packed cell.
 *
 * @return packed_task if the interaction was packed without flushing, or
 *         flushed_self_task if the batch was flushed.
 */
  enum runner_gpu_task_type runner_doself_grav_pp_task_gpu(
    struct runner *r,
    struct gpu_runner_substream *substream,
    struct cell *ci,
    struct task *t,
    int ncells,
    int max_cell_size) {
    
    runner_gpu_bind_device(r);

  if (ci->nodeID != r->e->nodeID)
    error("GPU self task attempted to pack a foreign cell.");

  #ifdef SWIFT_GPU_TIMING 
  const double pack_time_start = runner_gpu_walltime_s();
  #endif

  #ifdef SWIFT_USE_TIMERS
  const ticks self_pack_tic = getticks();
  #endif

  const struct engine *e = r->e;
  const int periodic = e->mesh->periodic;
  const double min_trunc = e->mesh->r_cut_min;
  struct gravity_cache *const ci_cache = &r->ci_gravity_cache;

  const int gcount = ci->grav.count;
  const int gcount_padded = gcount - (gcount % VEC_SIZE) + VEC_SIZE;

  if (gcount > max_cell_size)
    error("More particles than allocated memory!");

  const int slot = substream->grav_batch_self_count;

  const double loc[3] = {
      ci->loc[0] + 0.5 * ci->width[0],
      ci->loc[1] + 0.5 * ci->width[1],
      ci->loc[2] + 0.5 * ci->width[2]};

  gravity_cache_populate_no_mpole(
      e->max_active_bin, ci_cache,
      ci->grav.parts,
      gcount, gcount_padded,
      loc, ci,
      e->gravity_properties);

  while (cell_glocktree(ci)) {
    ;
  }

  /*record packed offset/count */
  append_packed_self_cell(substream, slot, gcount);
  substream->self_rmax_h[slot] = 2.f * ci->grav.multipole->r_max;
  const int offset = substream->self_offsets_h[slot];

  /* Build compact active-target list for this cell. */
  int active_count = 0;
  const int active_base = substream->self_total_active_count;
  substream->self_active_offsets_h[slot] = active_base;

  for (int i = 0; i < gcount; i++) {
  	if (ci_cache->active[i] > 0) {
    		substream->self_active_index_h[active_base + active_count] = i;
    		active_count++;
  	}
  }

  substream->self_total_active_count += active_count;
  
  substream->self_active_counts_h[slot] = active_count;
  if (active_count > substream->self_max_active_count)
    substream->self_max_active_count = active_count;

  /* Pack contiguously */
  substream->self_cell_flags_h[slot] =
    cell_is_active_gravity(ci, e) ? 1 : 0;

  substream->self_rmax_h[slot] =
    2.f * ci->grav.multipole->r_max;

  substream->self_use_full_h[slot] =
    (!periodic || substream->self_rmax_h[slot] <= min_trunc) ? 1 : 0;

  for (int i = 0; i < gcount; i++) {
  	const int k = offset + i;

  	substream->send_self_pos_mass[k].x = ci_cache->x[i];
  	substream->send_self_pos_mass[k].y = ci_cache->y[i];
  	substream->send_self_pos_mass[k].z = ci_cache->z[i];
  	substream->send_self_pos_mass[k].w = ci_cache->m[i];

  	substream->send_self_h[k] = ci_cache->epsilon[i];
  }

  substream->grav_cells_self[slot] = ci;
  substream->grav_tasks_self[slot] = t;
  
  runner_gpu_count_self_task(r, &r->e->sched, t);
  
  substream->grav_batch_self_count++;

  gravity_cache_zero_output(ci_cache, gcount_padded);
  cell_gunlocktree(ci);

  #ifdef SWIFT_GPU_TIMING 
  runner_gpu_self_pack_time_s += runner_gpu_walltime_s() - pack_time_start;
  #endif

  #ifdef SWIFT_USE_TIMERS
  timers_toc(timer_gpu_grav_self_pack_cpu, self_pack_tic);
  #endif

#ifdef SWIFT_DEBUG_CHECKS
  for (int j = 0; j < gcount; j++) {
    for (int i = 0; i < gcount; i++) {
      if (i == j) continue;
      accumulate_inc_ll(&ci->grav.parts[j].num_interacted);
    }
  }
#endif

  /* ===================== FLUSH ===================== */

    if (substream->grav_batch_self_count >= ncells) {

    const int nslots = substream->grav_batch_self_count;
    const int total = substream->self_total_count;
    
    #ifdef SWIFT_DEBUG_TASKS

  ticks gpu_anchor_tic = 0;
  ticks gpu_sync_toc = 0;

  const int gpu_substream_id =
      (int)(substream - r->gpu.substreams);

#endif

    #ifdef SWIFT_GPU_TIMING

  double h2d_s = 0.0;
  double kernel_s = 0.0;
  double d2h_s = 0.0;
  double unpack_s = 0.0;

#endif

#if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)

  GPUEvent h2d_start, h2d_stop;
  GPUEvent kernel_start, kernel_stop;
  GPUEvent d2h_start, d2h_stop;

  GPUEventCreate(&h2d_start);
  GPUEventCreate(&h2d_stop);

  GPUEventCreate(&kernel_start);
  GPUEventCreate(&kernel_stop);

  GPUEventCreate(&d2h_start);
  GPUEventCreate(&d2h_stop);

#endif 

    /* copy packed metadata */
    #ifdef SWIFT_DEBUG_TASKS
  gpu_anchor_tic = getticks();
#endif

#if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
  GPUEventRecord(h2d_start, substream->stream);
#endif
    GPUMemcpyAsync(
        substream->self_counts_d,
        substream->self_counts_h,
        nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        substream->stream);

    GPUMemcpyAsync(
        substream->self_offsets_d,
        substream->self_offsets_h,
        nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        substream->stream);

    GPUMemcpyAsync(
        substream->self_active_counts_d,
        substream->self_active_counts_h,
        nslots * sizeof(int),
        GPU_MEMCPY_HOST_TO_DEVICE,
        substream->stream);

    GPUMemcpyAsync(
    substream->self_active_offsets_d,
    substream->self_active_offsets_h,
    nslots * sizeof(int),
    GPU_MEMCPY_HOST_TO_DEVICE,
    substream->stream);

    GPUMemcpyAsync(
    substream->self_active_index_d,
    substream->self_active_index_h,
    (size_t)substream->self_total_active_count * sizeof(int),
    GPU_MEMCPY_HOST_TO_DEVICE,
    substream->stream);

    /* H2D: only live data */
    GPUMemcpyAsync(
    	substream->self_cell_flags_d,
    	substream->self_cell_flags_h,
    	nslots * sizeof(int),
    	GPU_MEMCPY_HOST_TO_DEVICE,
    	substream->stream);

	GPUMemcpyAsync(
	    substream->self_use_full_d,
	    substream->self_use_full_h,
	    nslots * sizeof(int),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

	GPUMemcpyAsync(
	    substream->send_self_pos_mass_d,
	    substream->send_self_pos_mass,
	    total * sizeof(float4),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

	GPUMemcpyAsync(
	    substream->send_self_h_d,
	    substream->send_self_h,
	    total * sizeof(float),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(h2d_stop, substream->stream);
    #endif 

    /* kernel */
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(kernel_start, substream->stream);
    #endif 

    runner_doself_grav_pp_flush(
    r, substream, nslots, substream->self_max_active_count, max_cell_size, substream->stream);

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(kernel_stop, substream->stream);
    #endif

    /* D2H: only live data */
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS) 
    GPUEventRecord(d2h_start, substream->stream);
    #endif

    GPUMemcpyAsync(
	    substream->recv_self_active,
	    substream->recv_self_active_d,
	    (size_t)substream->self_total_active_count *
		sizeof(struct gravity_gpu_values_recv),
	    GPU_MEMCPY_DEVICE_TO_HOST,
	    substream->stream);

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS) 
    GPUEventRecord(d2h_stop, substream->stream);
    #endif
    GPUEventRecord(substream->done, substream->stream);
    #ifdef SWIFT_USE_TIMERS
    const ticks gpu_wait_tic = getticks();
    #endif

    GPUEventSynchronize(substream->done);

    #ifdef SWIFT_USE_TIMERS
    timers_toc(timer_gpu_grav_wait, gpu_wait_tic);
    #endif
    
    #ifdef SWIFT_DEBUG_TASKS
  gpu_sync_toc = getticks();
#endif

    #ifdef SWIFT_GPU_TIMING 
    h2d_s = runner_gpu_event_elapsed_s(h2d_start, h2d_stop);
    kernel_s = runner_gpu_event_elapsed_s(kernel_start, kernel_stop);
    d2h_s = runner_gpu_event_elapsed_s(d2h_start, d2h_stop);
    #endif
    
    #ifdef SWIFT_DEBUG_TASKS

  const double gpu_h2d_end_s =
      runner_gpu_event_offset_s(h2d_start, h2d_stop);

  const double gpu_kernel_start_s =
      runner_gpu_event_offset_s(h2d_start, kernel_start);

  const double gpu_kernel_end_s =
      runner_gpu_event_offset_s(h2d_start, kernel_stop);

  const double gpu_d2h_start_s =
      runner_gpu_event_offset_s(h2d_start, d2h_start);

  const double gpu_d2h_end_s =
      runner_gpu_event_offset_s(h2d_start, d2h_stop);

  runner_gpu_write_timeline_row(
      "self",
      (long long)r->e->step,
      r->id,
      gpu_substream_id,
      gpu_anchor_tic,
      gpu_sync_toc,
      gpu_h2d_end_s,
      gpu_kernel_start_s,
      gpu_kernel_end_s,
      gpu_d2h_start_s,
      gpu_d2h_end_s);

#endif

    /* ===================== UNPACK ===================== */

    #ifdef SWIFT_GPU_TIMING 
    const double unpack_t0 = runner_gpu_walltime_s();
    #endif

    for (int j = 0; j < nslots; j++) {

      struct cell *c_unpack = substream->grav_cells_self[j];
      const int count = substream->self_counts_h[j];
      const int offset = substream->self_offsets_h[j];

      while (cell_glocktree(c_unpack)) {
        ;
      }

      const int active_count = substream->self_active_counts_h[j];
	const int active_base = substream->self_active_offsets_h[j];

	for (int a = 0; a < active_count; a++) {
	  if (c_unpack == NULL)
  		error("GPU self unpack received NULL cell.");

	  if (c_unpack->nodeID != r->e->nodeID)
  		error("GPU self unpack attempted to write a foreign cell.");

	  if (!cell_is_active_gravity(c_unpack, r->e))
  		continue;
  
	  const int local_pid = substream->self_active_index_h[active_base + a];
	  const int k = active_base + a;

	  c_unpack->grav.parts[local_pid].a_grav[0] +=
	      substream->recv_self_active[k].values_i.x;
	  c_unpack->grav.parts[local_pid].a_grav[1] +=
	      substream->recv_self_active[k].values_i.y;
	  c_unpack->grav.parts[local_pid].a_grav[2] +=
	      substream->recv_self_active[k].values_i.z;
	  c_unpack->grav.parts[local_pid].potential +=
	      substream->recv_self_active[k].values_i.w;
	}

      cell_gunlocktree(c_unpack);
    }

    #ifdef SWIFT_GPU_TIMING 
    unpack_s = runner_gpu_walltime_s() - unpack_t0;

    const double self_pack_s = runner_gpu_self_pack_time_s;
    runner_gpu_self_pack_time_s = 0.0;
    
    const size_t h2d_bytes =
    runner_gpu_self_h2d_bytes(
        substream,
        nslots);
        
    runner_gpu_write_timing_row(
	    "self",
	    (long long)r->e->step,
	    r->id,
	    r->qid,
	    nslots,
	    total,
	    h2d_bytes,
	    self_pack_s,
	    h2d_s,
	    kernel_s,
	    d2h_s,
	    unpack_s);
    #endif
    
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)

  GPUEventDestroy(h2d_start);
  GPUEventDestroy(h2d_stop);

  GPUEventDestroy(kernel_start);
  GPUEventDestroy(kernel_stop);

  GPUEventDestroy(d2h_start);
  GPUEventDestroy(d2h_stop);

#endif

    /* ===================== COMPLETE TASKS ===================== */

    runner_gpu_complete_self_batch(r, &r->e->sched, substream, t);

    return flushed_self_task;
  }

  return packed_task;
}

/**
 * @brief Recursively process a pair-gravity interaction, offloading leaf
 *        P-P interactions to the GPU where appropriate.
 *
 * Recurses through the gravity cell hierarchy and selects between truncated
 * interactions, CPU multipole interactions, further cell splitting, and
 * direct particle-particle interactions. Leaf P-P interactions are packed
 * into the supplied GPU substream and may trigger a batch flush.
 *
 * @param r The #runner processing the interaction.
 * @param substream The GPU substream used for packed pair interactions.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Whether to record the sub-pair gravity timer.
 * @param grav_cells_pair Array storing the cell pointers associated with
 *                        packed pairs.
 * @param grav_tasks_pair Array storing the scheduler tasks associated with
 *                        packed pairs.
 * @param grav_pair_internal_from_self Array marking pair interactions generated
 *                                     internally from self-gravity walks.
 * @param t The top-level scheduler #task associated with this recursive walk.
 * @param internal_from_self Whether this pair interaction originates from a
 *                           self-gravity walk.
 * @param ncells Maximum number of cell slots available in the GPU pair batch.
 * @param max_cell_size Maximum number of gravity particles permitted in a
 *                      packed cell.
 * @param stream GPU stream associated with the batch.
 *
 * @return regular_task if no GPU leaf interaction was generated, packed_task
 *         if GPU work remains packed in the batch, or flushed_pair_task if
 *         GPU work was flushed during the recursive walk.
 */
enum runner_gpu_task_type runner_dopair_recursive_grav_gpu(
    struct runner *r, struct gpu_runner_substream *substream, struct cell *ci,
    struct cell *cj, const int gettimer,
    struct cell **grav_cells_pair, struct task **grav_tasks_pair,
    unsigned char *grav_pair_internal_from_self,
    struct task *t, int internal_from_self,
    int ncells, int max_cell_size, GPUStream stream) {

  if (ci == NULL || cj == NULL)
    error("runner_dopair_recursive_grav_gpu got NULL cell");

  const struct engine *e = r->e;

  if (!cell_are_gpart_drifted(ci, e))
    cell_drift_gpart(ci, e, /*force=*/1, /*init_particles=*/0, NULL);
  if (!cell_are_gpart_drifted(cj, e))
    cell_drift_gpart(cj, e, /*force=*/1, /*init_particles=*/0, NULL);

  /* Clear the flags */
  runner_clear_grav_flags(ci, e);
  runner_clear_grav_flags(cj, e);

  /* Some constants */
  const int nodeID = e->nodeID;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const double max_distance = e->mesh->r_cut_max;

  /* Anything to do here? */
  if (!((cell_is_active_gravity(ci, e) && ci->nodeID == nodeID) ||
        (cell_is_active_gravity(cj, e) && cj->nodeID == nodeID)))
    return regular_task;

#ifdef SWIFT_DEBUG_CHECKS

  const int gcount_i = ci->grav.count;
  const int gcount_j = cj->grav.count;

  /* Early abort? */
  if (gcount_i == 0 || gcount_j == 0)
    error("Doing pair gravity on an empty cell !");

  /* Sanity check */
  if (ci == cj) error("Pair interaction between a cell and itself.");

  if (cell_is_active_gravity(ci, e) &&
      ci->grav.ti_old_multipole != e->ti_current)
    error("ci->grav.multipole not drifted.");
  if (cell_is_active_gravity(cj, e) &&
      cj->grav.ti_old_multipole != e->ti_current)
    error("cj->grav.multipole not drifted.");
#endif

  TIMER_TIC;

  /* Recover the multipole information */
  struct gravity_tensors *const multi_i = ci->grav.multipole;
  struct gravity_tensors *const multi_j = cj->grav.multipole;

  /* Get the distance between the CoMs */
  double dx = multi_i->CoM[0] - multi_j->CoM[0];
  double dy = multi_i->CoM[1] - multi_j->CoM[1];
  double dz = multi_i->CoM[2] - multi_j->CoM[2];

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  /* Minimal distance between any 2 particles in the two cells */
  const double r_lr_check = sqrt(r2) - (multi_i->r_max + multi_j->r_max);

  /* Are we beyond the distance where the truncated forces are 0? */
  if (periodic && r_lr_check > max_distance) {

#ifdef SWIFT_DEBUG_CHECKS
    if (cell_is_active_gravity(ci, e))
      accumulate_add_ll(&multi_i->pot.num_interacted,
                        multi_j->m_pole.num_gpart);
    if (cell_is_active_gravity(cj, e))
      accumulate_add_ll(&multi_j->pot.num_interacted,
                        multi_i->m_pole.num_gpart);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
    /* Need to account for the interactions we missed */
    if (cell_is_active_gravity(ci, e))
      accumulate_add_ll(&multi_i->pot.num_interacted_pm,
                        multi_j->m_pole.num_gpart);
    if (cell_is_active_gravity(cj, e))
      accumulate_add_ll(&multi_j->pot.num_interacted_pm,
                        multi_i->m_pole.num_gpart);
#endif
    if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
    return regular_task;
  }

  /* OK, we actually need to compute this pair. Let's find the cheapest
   * option... */

  if (ci->grav.count <= 1 || cj->grav.count <= 1) {

    /* We have two cheap cells. Go P-P. */
    #ifdef SWIFT_USE_TIMERS
    const ticks cpu_fallback_tic = getticks();
    #endif

    runner_dopair_recursive_grav(r, ci, cj, 0);

    #ifdef SWIFT_USE_TIMERS
    timers_toc(timer_gpu_grav_cpu_fallback, cpu_fallback_tic);
    #endif

    if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
    return regular_task;

    /* Can we use M-M interactions ? */
  } else if (gravity_M2L_accept_symmetric(e->gravity_properties, multi_i,
                                          multi_j, r2,
                                          /*use_rebuild_sizes=*/0, periodic)) {

    /* Go M-M */
    #ifdef SWIFT_USE_TIMERS
    const ticks cpu_fallback_tic = getticks();
    #endif

    runner_dopair_recursive_grav(r, ci, cj, 0);

    #ifdef SWIFT_USE_TIMERS
    timers_toc(timer_gpu_grav_cpu_fallback, cpu_fallback_tic);
    #endif

    if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
    return regular_task;

    /* Did we reach the bottom? */
  } else if (!ci->split && !cj->split) {

    /* We have two leaves. Go P-P. */
    const enum runner_gpu_task_type result = runner_dopair_grav_pp_gpu(
    r, substream, ci, cj, 1, 1,
    substream->grav_cells_pair,
    substream->grav_tasks_pair,
    substream->grav_pair_internal_from_self,
    t, internal_from_self,
    ncells, max_cell_size, substream->stream);

    if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
    return result;

  } else {

    enum runner_gpu_task_type task_type = regular_task;

    /* Alright, we'll have to split and recurse. */
    /* We know at least one of ci and cj is splittable */

    const double ri_max = multi_i->r_max;
    const double rj_max = multi_j->r_max;

    /* Split the larger of the two cells and start over again */
    if (ri_max > rj_max) {

      /* Can we actually split that interaction ? */
      if (ci->split) {

        /* Loop over ci's children */
        for (int k = 0; k < 8; k++) {
          if (ci->progeny[k] != NULL) {
            enum runner_gpu_task_type child_type =
    		runner_dopair_recursive_grav_gpu(
        		r, substream, ci->progeny[k], cj, 0,
        		substream->grav_cells_pair, substream->grav_tasks_pair,
        		substream->grav_pair_internal_from_self,
        		t, internal_from_self, ncells, max_cell_size, substream->stream);
            if (child_type > task_type) task_type = child_type;
          }
        }

      } else {
        /* cj is split */

        /* MATTHIEU: This could maybe be replaced by P-M interactions ?  */

        /* Loop over cj's children */
        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL) {
            enum runner_gpu_task_type child_type =
                runner_dopair_recursive_grav_gpu(
                    r, substream, ci, cj->progeny[k], 0,
		    substream->grav_cells_pair, substream->grav_tasks_pair,
		    substream->grav_pair_internal_from_self,
		    t, internal_from_self, ncells, max_cell_size, substream->stream);
            if (child_type > task_type) task_type = child_type;
          }
        }
      }
    } else {

      /* Can we actually split that interaction ? */
      if (cj->split) {

        /* Loop over cj's children */
        for (int k = 0; k < 8; k++) {
          if (cj->progeny[k] != NULL) {
            enum runner_gpu_task_type child_type =
                runner_dopair_recursive_grav_gpu(
                    r, substream, ci, cj->progeny[k], 0,
		    substream->grav_cells_pair, substream->grav_tasks_pair,
		    substream->grav_pair_internal_from_self,
		    t, internal_from_self, ncells, max_cell_size, substream->stream);
            if (child_type > task_type) task_type = child_type;
          }
        }

      } else {
        /* ci is split */

        /* MATTHIEU: This could maybe be replaced by P-M interactions ?  */

        /* Loop over ci's children */
        for (int k = 0; k < 8; k++) {
          if (ci->progeny[k] != NULL) {
            enum runner_gpu_task_type child_type =
                runner_dopair_recursive_grav_gpu(
                    r, substream, ci->progeny[k], cj, 0,
		    substream->grav_cells_pair, substream->grav_tasks_pair,
		    substream->grav_pair_internal_from_self,
		    t, internal_from_self, ncells, max_cell_size, substream->stream);
            if (child_type > task_type) task_type = child_type;
          }
        }
      }
    }

    /* Determine the return type based on whether *this* walk actually
       produced any leaf pairs.  task_type tracks the highest child return
       value: packed_task or flushed_pair_task means at least one leaf
       pair was generated by this walk. */
    enum runner_gpu_task_type final_type;
    if (task_type >= packed_task) {
      /* This walk produced leaf pairs.  If some are still in the buffer
         they will be flushed later; if all were flushed the task is done. */
      if (substream->grav_batch_pair_count > 0) {
        final_type = packed_task;
      } else {
        final_type = flushed_pair_task;
      }
    } else {
      /* No leaf pairs were produced at all (all M-M or truncated). */
      final_type = regular_task;
    }

    if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
    return final_type;
  }

  if (gettimer) TIMER_TOC(timer_dosub_pair_grav);
  return regular_task;
}


/**
 * @brief Choose the number of cell slots to allocate per GPU gravity batch.
 *
 * Estimates a safe batch size from the currently available GPU memory,
 * maximum cell size, number of GPU streams, number of runner threads, and the
 * number of local MPI ranks sharing the device. A user-specified
 * GPU:ncells_per_gpu_grav_pack value is honoured where it fits within the
 * calculated memory budget.
 *
 * @param e The #engine containing the GPU configuration and runtime state.
 * @param max_cell_size Maximum number of gravity particles permitted in a
 *                      packed cell.
 *
 * @return The selected number of cell slots per GPU gravity batch.
 */
static int runner_gpu_choose_batch_ncells(const struct engine *e,
                                          int max_cell_size) {
                                          
  /* User override (<=0 means "not set") */
  const int user_ncells = e->ncells_per_gpu_grav_pack;

  size_t free_bytes = 0, total_bytes = 0;
  GPUMemGetInfo(&free_bytes, &total_bytes);

  /* Leave headroom for CUDA/HIP context/runtime overhead */
  const double usable_fraction = 0.80;
  const size_t usable_bytes = (size_t)(free_bytes * usable_fraction);
  
  const int nstreams = parser_get_opt_param_int(
      e->parameter_file, "GPU:nstreams", 1);

  /* Estimate memory for self and pair gravity buffers. */
  const size_t bytes_per_particle =
    2 * (sizeof(float4) + sizeof(float)) +
    2 * sizeof(struct gravity_gpu_values_recv);

  const size_t bytes_per_cell_per_substream =
    (size_t)max_cell_size * bytes_per_particle;

  const size_t bytes_per_cell_per_runner =
      (size_t)nstreams * bytes_per_cell_per_substream;

  /* Replace with the actual number of runners that allocate GPU buffers */
  const int nr_threads = e->nr_threads > 0 ? e->nr_threads : 1;
  const int nr_gpu_runners =
    nr_threads * runner_gpu_local_ranks_on_device_for_budget;

  const size_t metadata = 64ULL * 1024ULL * 1024ULL; /* 64 MB */

  size_t budget = 0;
  if (usable_bytes > metadata)
    budget = usable_bytes - metadata;

  const size_t per_runner_budget = budget / (size_t)nr_gpu_runners;

  int ncells = (int)(per_runner_budget / bytes_per_cell_per_runner);

  if (ncells < 2) ncells = 2;
  if (ncells > 10000) ncells = 10000;

  if (user_ncells > 0) {

    if (user_ncells > ncells) {
      if (e->verbose) {
        message("GPU:ncells_per_gpu_grav_pack=%d too large, limiting to %d",
                user_ncells, ncells);
      }
      return ncells;
    }

    return user_ncells;
  }

  return ncells;
}


/**
 * @brief Choose a safe number of GPU substreams per runner.
 *
 * Estimates the maximum number of substreams that can be allocated from the
 * available GPU memory after accounting for the configured batch size,
 * maximum cell size, runner threads, and local MPI ranks sharing the device.
 *
 * @param e The #engine containing the GPU configuration and runtime state.
 * @param max_cell_size Maximum number of gravity particles permitted in a
 *                      packed cell.
 * @param ncells Number of cell slots allocated per GPU gravity batch.
 *
 * @return The maximum safe number of GPU substreams per runner.
 */
static int runner_gpu_choose_nstreams(const struct engine *e,
                                      int max_cell_size,
                                      int ncells) {

  size_t free_bytes = 0, total_bytes = 0;
  GPUMemGetInfo(&free_bytes, &total_bytes);

  /* Leave headroom for other allocations */
  const double usable_fraction = 0.80;
  const size_t usable_bytes = (size_t)(free_bytes * usable_fraction);

  const int nr_threads = e->nr_threads > 0 ? e->nr_threads : 1;
  const int nr_gpu_runners =
    nr_threads * runner_gpu_local_ranks_on_device_for_budget;
  const size_t metadata = 64ULL * 1024ULL * 1024ULL; /* 64 MB */

  size_t budget = 0;
  if (usable_bytes > metadata)
    budget = usable_bytes - metadata;

  const size_t per_runner_budget = budget / (size_t)nr_gpu_runners;

  const size_t bytes_per_substream =
      (size_t)ncells * (size_t)max_cell_size *
      (2 * sizeof(struct gravity_gpu_values_send) +
       2 * sizeof(struct gravity_gpu_values_recv));

  if (bytes_per_substream == 0)
    return 1;

  int max_safe_nstreams = (int)(per_runner_budget / bytes_per_substream);

  if (max_safe_nstreams < 1) max_safe_nstreams = 1;

  return max_safe_nstreams;
}

/**
 * @brief Select the GPU device used by a runner.
 *
 * Determines the local MPI rank and assigns it to a visible GPU. If a device
 * is explicitly specified using GPU:device_id, that device is used instead.
 * The function also determines how many local MPI ranks share the selected
 * device for use in GPU memory budgeting.
 *
 * @param r The #runner whose GPU device is to be selected.
 *
 * @return The selected GPU device ID.
 */
static int runner_gpu_select_device(struct runner *r) {

  struct engine *e = r->e;
  struct gpu_runner *gpu = &r->gpu;

  int ngpu = 0;
  GPUGetDeviceCount(&ngpu);

  if (ngpu <= 0)
    error("No CUDA/HIP GPU visible to this MPI rank.");

  int local_rank = 0;
  int local_size = 1;

#ifdef WITH_MPI
  MPI_Comm local_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &local_comm);

  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &local_size);

  MPI_Comm_free(&local_comm);
#endif

  const int user_device = parser_get_opt_param_int(
      e->parameter_file, "GPU:device_id", -1);

  int device_id = 0;

  if (user_device >= 0) {
    if (user_device >= ngpu)
      error("GPU:device_id=%d requested but only %d GPU(s) visible.",
            user_device, ngpu);

    device_id = user_device;

  } else {
    device_id = local_rank % ngpu;
  }

  int local_ranks_on_device = 1;

#ifdef WITH_MPI
  /*
   * Count ranks on this node that will map to the same GPU.
   * This is used only for conservative memory budgeting.
   */
  MPI_Comm local_comm2;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &local_comm2);

  const int my_device_id = device_id;
  int *all_devices = malloc((size_t)local_size * sizeof(int));

  if (all_devices == NULL)
    error("Failed to allocate local GPU mapping array.");

  MPI_Allgather(&my_device_id, 1, MPI_INT,
                all_devices, 1, MPI_INT, local_comm2);

  local_ranks_on_device = 0;

  for (int i = 0; i < local_size; i++) {
    if (all_devices[i] == my_device_id)
      local_ranks_on_device++;
  }

  free(all_devices);
  MPI_Comm_free(&local_comm2);
#endif

  gpu->device_id = device_id;
  gpu->local_mpi_rank = local_rank;
  gpu->local_mpi_size = local_size;
  gpu->local_ranks_on_device = local_ranks_on_device;

  GPUSetDevice(device_id);
  runner_gpu_check_error("GPUSetDevice");

  return device_id;
}

/**
 * @brief Initialise the GPU-specific state attached to a runner.
 *
 * @param r The runner whose GPU state to initialise.
 */
void runner_gpu_init(struct runner *r) {

  struct gpu_runner *gpu = &r->gpu;
  struct engine *e = r->e;
  
  #ifdef WITH_MPI
  MPI_Comm local_comm;
  int local_rank = 0, local_size = 1;

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &local_comm);
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &local_size);
  MPI_Comm_free(&local_comm);
  #else
  int local_rank = 0;
  #endif

  int ngpu = 0;
  GPUGetDeviceCount(&ngpu);

  if (ngpu <= 0)
    error("No GPU visible to MPI rank");

  const int device_id = runner_gpu_select_device(r);
  
  runner_gpu_local_ranks_on_device_for_budget =
    gpu->local_ranks_on_device > 0 ? gpu->local_ranks_on_device : 1;

  GPUDeviceProp prop;
  GPUGetDeviceProperties(&prop, device_id);
  runner_gpu_check_error("GPUGetDeviceProperties");

  if (e->verbose && r->id == 0) {
    message("MPI rank %d local_rank=%d using GPU device %d "
          "(%d visible GPU(s), %d local rank(s) sharing this device)",
          e->nodeID,
          gpu->local_mpi_rank,
          gpu->device_id,
          gpu->local_mpi_size,
          gpu->local_ranks_on_device);
  }

  const int max_cell_size = space_subsize_self_grav + 100;

  gpu->grav_max_cell_size = max_cell_size;
  gpu->grav_batch_ncells =
      runner_gpu_choose_batch_ncells(e, gpu->grav_max_cell_size);

  /* Pair batches consume 2 slots at a time */
  if (gpu->grav_batch_ncells % 2 != 0)
    gpu->grav_batch_ncells--;

  if (gpu->grav_batch_ncells < 2)
    gpu->grav_batch_ncells = 2;

  /* User may request nstreams, otherwise auto-pick based on ncells/max_cell_size */
  const int user_nstreams = parser_get_opt_param_int(
      e->parameter_file, "GPU:nstreams", 1);

  const int auto_nstreams =
      runner_gpu_choose_nstreams(e,
                                 gpu->grav_max_cell_size,
                                 gpu->grav_batch_ncells);

  if (user_nstreams > 0) {
    gpu->nstreams = user_nstreams;
    if (gpu->nstreams > auto_nstreams) {
      if (r->id == 0) {
        message("GPU:nstreams=%d too large for current ncells/max_cell_size, "
                "limiting to %d",
                user_nstreams, auto_nstreams);
      }
      gpu->nstreams = auto_nstreams;
    }
  } else {
    gpu->nstreams = auto_nstreams;
  }

  if (gpu->nstreams < 1) gpu->nstreams = 1;
  if (gpu->nstreams > 8) gpu->nstreams = 8;

  gpu->substreams = malloc((size_t)gpu->nstreams *
                           sizeof(struct gpu_runner_substream));
  if (gpu->substreams == NULL)
    error("Failed to allocate GPU substreams");

  const size_t nparticles =
    (size_t)gpu->grav_batch_ncells *
    (size_t)gpu->grav_max_cell_size;

  const size_t pos_mass_bytes =
    nparticles * sizeof(float4);

  const size_t h_bytes =
    nparticles * sizeof(float);

  const size_t recv_bytes =
    nparticles * sizeof(struct gravity_gpu_values_recv);

  /* Conservative estimate for self and pair gravity buffers. */
  const size_t bytes_per_substream =
    2 * (pos_mass_bytes + h_bytes + recv_bytes);

  size_t free_bytes = 0, total_bytes = 0;
  GPUMemGetInfo(&free_bytes, &total_bytes);

  if (r->id == 0) {
    message("GPU device: %s", prop.name);
    message("GPU free memory: %.2f GB",
            free_bytes / (1024.0 * 1024.0 * 1024.0));
    message("GPU total memory: %.2f GB",
            total_bytes / (1024.0 * 1024.0 * 1024.0));
    message("Max cell size: %i", gpu->grav_max_cell_size);
    message("ncells per pack: %i", gpu->grav_batch_ncells);
    message("Streams per runner: %i", gpu->nstreams);
    message("Per-substream buffer bytes: %zu", bytes_per_substream);
  }

  gpu->next_substream = 0;

  for (int i = 0; i < gpu->nstreams; i++) {

    struct gpu_runner_substream *substream = &gpu->substreams[i];

    GPUStreamCreateWithFlags(&substream->stream, GPUStreamNonBlocking);
    GPUEventCreate(&substream->done);
    substream->busy = 0;

    /* ---------- Pair state ---------- */

    substream->grav_batch_pair_count = 0;
    substream->pair_unique_cell_count = 0;

    substream->pair_unique_cells =
    malloc((size_t)gpu->grav_batch_ncells * sizeof(struct cell *));

    substream->pair_pair_i_h =
    malloc((size_t)gpu->grav_batch_ncells * sizeof(int));

    substream->pair_pair_j_h =
    malloc((size_t)gpu->grav_batch_ncells * sizeof(int));

    GPUMalloc((void **)&substream->pair_pair_i_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));

    GPUMalloc((void **)&substream->pair_pair_j_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));
          
    const int pair_capacity = gpu->grav_batch_ncells / 2;

    substream->pair_total_pair_active_count = 0;

    substream->pair_use_full_h =
    	malloc((size_t)pair_capacity * sizeof(int));

    GPUMalloc((void **)&substream->pair_use_full_d,
          (size_t)pair_capacity * sizeof(int));

    substream->pair_side_active_offsets_h =
        malloc((size_t)(2 * pair_capacity) * sizeof(int));

    GPUMalloc((void **)&substream->pair_side_active_offsets_d,
          (size_t)(2 * pair_capacity) * sizeof(int));
    
    const size_t pos_mass_bytes =
    (size_t)gpu->grav_batch_ncells *
    (size_t)gpu->grav_max_cell_size *
    sizeof(float4);

	const size_t h_bytes =
	    (size_t)gpu->grav_batch_ncells *
	    (size_t)gpu->grav_max_cell_size *
	    sizeof(float);

	GPUMalloc((void **)&substream->send_pair_pos_mass_d, pos_mass_bytes);
	GPUHostMalloc((void **)&substream->send_pair_pos_mass, pos_mass_bytes);

	GPUMalloc((void **)&substream->send_pair_h_d, h_bytes);
	GPUHostMalloc((void **)&substream->send_pair_h, h_bytes);
    
    const size_t active_recv_bytes =
    (size_t)gpu->grav_batch_ncells *
    (size_t)gpu->grav_max_cell_size *
    sizeof(struct gravity_gpu_values_recv);

    GPUMalloc((void **)&substream->recv_pair_active_d, active_recv_bytes);
    GPUHostMalloc((void **)&substream->recv_pair_active, active_recv_bytes);

    substream->grav_cells_pair =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(struct cell *));
    substream->grav_tasks_pair =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(struct task *));
    substream->grav_pair_internal_from_self =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(unsigned char));

    if (substream->grav_pair_internal_from_self != NULL) {
      memset(substream->grav_pair_internal_from_self, 0,
             (size_t)gpu->grav_batch_ncells * sizeof(unsigned char));
    }

    substream->pair_total_count = 0;
    substream->pair_max_active_count = 0;
    substream->pair_total_active_count = 0;
    substream->pair_total_pair_active_count = 0;

    substream->pair_counts_h =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->pair_offsets_h =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->pair_active_counts_h =
    	malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->pair_active_offsets_h =
    	malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->pair_active_index_h =
    	malloc((size_t)gpu->grav_batch_ncells *
           (size_t)gpu->grav_max_cell_size * sizeof(int));

    GPUMalloc((void **)&substream->pair_counts_d,
              (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->pair_offsets_d,
              (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->pair_active_counts_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->pair_active_offsets_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->pair_active_index_d,
          (size_t)gpu->grav_batch_ncells *
          (size_t)gpu->grav_max_cell_size * sizeof(int));
          
    substream->pair_cell_flags_h =
    malloc((size_t)gpu->grav_batch_ncells * sizeof(int));

    GPUMalloc((void **)&substream->pair_cell_flags_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));

    for (int j = 0; j < gpu->grav_batch_ncells; j++) {
      substream->pair_counts_h[j] = 0;
      substream->pair_offsets_h[j] = 0;
      substream->pair_active_counts_h[j] = 0;
      substream->pair_active_offsets_h[j] = 0;
      substream->pair_unique_cells[j] = NULL;
      substream->pair_pair_i_h[j] = 0;
      substream->pair_pair_j_h[j] = 0;
    }

    /* ---------- Self state ---------- */

    substream->grav_batch_self_count = 0;

    GPUMalloc((void **)&substream->send_self_pos_mass_d, pos_mass_bytes);
    GPUHostMalloc((void **)&substream->send_self_pos_mass, pos_mass_bytes);

    GPUMalloc((void **)&substream->send_self_h_d, h_bytes);
    GPUHostMalloc((void **)&substream->send_self_h, h_bytes);
    
    GPUMalloc((void **)&substream->recv_self_active_d, active_recv_bytes);
    GPUHostMalloc((void **)&substream->recv_self_active, active_recv_bytes);

    substream->grav_cells_self =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(struct cell *));
    substream->grav_tasks_self =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(struct task *));

    substream->self_total_count = 0;
    substream->self_max_active_count = 0;
    substream->self_total_active_count = 0;

    substream->self_counts_h =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->self_offsets_h =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->self_active_counts_h =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->self_active_offsets_h =
    	malloc((size_t)gpu->grav_batch_ncells * sizeof(int));
    substream->self_active_index_h =
    	malloc((size_t)gpu->grav_batch_ncells *
           (size_t)gpu->grav_max_cell_size * sizeof(int));
    substream->self_rmax_h =
        malloc((size_t)gpu->grav_batch_ncells * sizeof(float));

    GPUMalloc((void **)&substream->self_counts_d,
              (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->self_offsets_d,
              (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->self_active_counts_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->self_active_offsets_d,
          (size_t)gpu->grav_batch_ncells * sizeof(int));
    GPUMalloc((void **)&substream->self_active_index_d,
          (size_t)gpu->grav_batch_ncells *
          (size_t)gpu->grav_max_cell_size * sizeof(int));
          
    substream->self_cell_flags_h =
    malloc((size_t)gpu->grav_batch_ncells * sizeof(int));

	substream->self_use_full_h =
	    malloc((size_t)gpu->grav_batch_ncells * sizeof(int));

	GPUMalloc((void **)&substream->self_cell_flags_d,
		  (size_t)gpu->grav_batch_ncells * sizeof(int));

	GPUMalloc((void **)&substream->self_use_full_d,
		  (size_t)gpu->grav_batch_ncells * sizeof(int));
		  
    for (int j = 0; j < gpu->grav_batch_ncells; j++) {
	  substream->self_counts_h[j] = 0;
	  substream->self_offsets_h[j] = 0;
	  substream->self_active_counts_h[j] = 0;
	  substream->self_active_offsets_h[j] = 0;
	  substream->self_rmax_h[j] = 0.f;

	  substream->self_cell_flags_h[j] = 0;
	  substream->self_use_full_h[j] = 0;
	}

    /* ---------- checks ---------- */

    if (substream->grav_cells_pair == NULL ||
        substream->grav_tasks_pair == NULL ||
        substream->grav_pair_internal_from_self == NULL ||
        substream->pair_counts_h == NULL ||
        substream->pair_offsets_h == NULL ||
        substream->pair_cell_flags_h == NULL ||
        substream->pair_active_counts_h == NULL ||
        substream->pair_unique_cells == NULL ||
        substream->pair_pair_i_h == NULL ||
        substream->pair_pair_j_h == NULL ||
        substream->pair_use_full_h == NULL ||
	substream->pair_side_active_offsets_h == NULL ||
        substream->pair_active_index_h == NULL)
      error("Failed to allocate runner GPU pair substream metadata arrays.");

    if (substream->grav_cells_self == NULL ||
    substream->grav_tasks_self == NULL ||
    substream->self_counts_h == NULL ||
    substream->self_offsets_h == NULL ||
    substream->self_active_counts_h == NULL ||
    substream->self_active_offsets_h == NULL ||
    substream->self_active_index_h == NULL ||
    substream->self_cell_flags_h == NULL ||
    substream->self_use_full_h == NULL ||
    substream->self_rmax_h == NULL)
  error("Failed to allocate runner GPU self substream metadata arrays.");
  }

  const GPUError err = GPUGetLastError();
  if (err != GPU_SUCCESS)
    error("runner_gpu_init failed: %s", GPUGetErrorString(err));
}

/**
 * @brief Acquire the substream for the GPU work to be launched to
 *
 * @param r The runner
 */
struct gpu_runner_substream *
runner_gpu_acquire_substream(struct runner *r) {

  runner_gpu_bind_device(r);

  struct gpu_runner *gpu = &r->gpu;

  for (int i = 0; i < gpu->nstreams; i++) {
    const int idx = (gpu->next_substream + i) % gpu->nstreams;
    struct gpu_runner_substream *substream = &gpu->substreams[idx];

    if (!substream->busy) {
      substream->busy = 1;
      gpu->next_substream = (idx + 1) % gpu->nstreams;
      return substream;
    }

    const GPUError qerr = GPUEventQuery(substream->done);

    if (qerr == GPU_SUCCESS) {
      substream->busy = 1;
      gpu->next_substream = (idx + 1) % gpu->nstreams;
      return substream;
    }

#if defined(WITH_CUDA)
    if (qerr != cudaErrorNotReady)
      error("runner_gpu_acquire_substream: GPUEventQuery failed: %s",
            GPUGetErrorString(qerr));
#elif defined(WITH_HIP)
    if (qerr != hipErrorNotReady)
      error("runner_gpu_acquire_substream: GPUEventQuery failed: %s",
            GPUGetErrorString(qerr));
#endif
  }

  struct gpu_runner_substream *substream =
      &gpu->substreams[gpu->next_substream];

  #ifdef SWIFT_USE_TIMERS
  const ticks gpu_wait_tic = getticks();
  #endif

  const GPUError serr = GPUEventSynchronize(substream->done);

  #ifdef SWIFT_USE_TIMERS
  timers_toc(timer_gpu_grav_wait, gpu_wait_tic);
  #endif
  if (serr != GPU_SUCCESS)
    error("runner_gpu_acquire_substream: GPUEventSynchronize failed: %s",
          GPUGetErrorString(serr));

  substream->busy = 1;
  gpu->next_substream = (gpu->next_substream + 1) % gpu->nstreams;

  return substream;
}

/**
 * @brief Clean the GPU-specific state attached to a runner.
 *
 * @param r The runner whose GPU state to clean.
 */
void runner_gpu_clean(struct runner *r) {

  runner_gpu_bind_device(r);

  struct gpu_runner *gpu = &r->gpu;

  for (int i = 0; i < gpu->nstreams; i++) {

    struct gpu_runner_substream *substream = &gpu->substreams[i];

    /* Pair buffers */
    GPUFreeHost(substream->send_pair_pos_mass);
    GPUFree(substream->send_pair_pos_mass_d);

    GPUFreeHost(substream->send_pair_h);
    GPUFree(substream->send_pair_h_d);;

    free(substream->grav_cells_pair);
    free(substream->grav_tasks_pair);
    free(substream->grav_pair_internal_from_self);
    
    free(substream->pair_counts_h);
    free(substream->pair_offsets_h);
    free(substream->pair_active_counts_h);
    free(substream->pair_active_index_h);
    GPUFree(substream->pair_counts_d);
    GPUFree(substream->pair_offsets_d);
    GPUFree(substream->pair_active_counts_d);
    GPUFree(substream->pair_active_index_d);
    
    GPUFreeHost(substream->recv_pair_active);
    GPUFree(substream->recv_pair_active_d);
    
    free(substream->pair_cell_flags_h);
    GPUFree(substream->pair_cell_flags_d);
    
    free(substream->pair_unique_cells);
    free(substream->pair_pair_i_h);
    free(substream->pair_pair_j_h);

    GPUFree(substream->pair_pair_i_d);
    GPUFree(substream->pair_pair_j_d);
    
    free(substream->pair_use_full_h);
    GPUFree(substream->pair_use_full_d);

    free(substream->pair_side_active_offsets_h);
    GPUFree(substream->pair_side_active_offsets_d);

    /* Self buffers */
    GPUFreeHost(substream->send_self_pos_mass);
    GPUFree(substream->send_self_pos_mass_d);

    GPUFreeHost(substream->send_self_h);
    GPUFree(substream->send_self_h_d);
    
    free(substream->self_cell_flags_h);
    GPUFree(substream->self_cell_flags_d);

    free(substream->self_use_full_h);
    GPUFree(substream->self_use_full_d);

    free(substream->grav_cells_self);
    free(substream->grav_tasks_self);
    
    free(substream->self_counts_h);
    free(substream->self_offsets_h);
    
    GPUFree(substream->self_counts_d);
    GPUFree(substream->self_offsets_d);
    
    free(substream->self_active_counts_h);
    free(substream->self_active_offsets_h);
    free(substream->self_active_index_h);

    GPUFree(substream->self_active_counts_d);
    GPUFree(substream->self_active_offsets_d);
    GPUFree(substream->self_active_index_d);
    
    GPUFreeHost(substream->recv_self_active);
    GPUFree(substream->recv_self_active_d);
    
    free(substream->self_rmax_h);

    /* Stream/event */
    GPUEventDestroy(substream->done);
    GPUStreamDestroy(substream->stream);

    /* Reset */
    substream->grav_cells_pair = NULL;
    substream->grav_tasks_pair = NULL;
    substream->grav_pair_internal_from_self = NULL;
    substream->grav_batch_pair_count = 0;
    substream->pair_counts_h = NULL;
    substream->pair_offsets_h = NULL;
    substream->pair_active_counts_h = NULL;
    substream->pair_active_index_h = NULL;
    substream->pair_counts_d = NULL;
    substream->pair_offsets_d = NULL;
    substream->pair_total_count = 0;

    substream->send_self_pos_mass = NULL;
    substream->send_self_pos_mass_d = NULL;
    substream->send_self_h = NULL;
    substream->send_self_h_d = NULL;
    substream->grav_cells_self = NULL;
    substream->grav_tasks_self = NULL;
    substream->grav_batch_self_count = 0;
    substream->self_counts_h = NULL;
    substream->self_offsets_h = NULL;
    substream->self_counts_d = NULL;
    substream->self_offsets_d = NULL;
    substream->self_total_count = 0;
    substream->self_rmax_h = NULL;

    substream->busy = 0;
  }

  gpu->next_substream = 0;
  gpu->grav_batch_ncells = 0;
  gpu->grav_max_cell_size = 0;
  
  free(gpu->substreams);
  gpu->substreams = NULL;
  gpu->nstreams = 0;
}

/**
 * @brief Flush any leftover packed self-gravity work owned by a runner.
 *
 * @param r The runner whose GPU batch should be flushed.
 * @return The outcome of the leftover flush attempt.
 */
enum runner_gpu_task_type runner_gpu_flush_leftover_self(struct runner *r) {

  runner_gpu_bind_device(r);

  enum runner_gpu_task_type result = regular_task;

  for (int l = 0; l < r->gpu.nstreams; l++) {
    struct gpu_runner_substream *substream = &r->gpu.substreams[l];
    const int nslots = substream->grav_batch_self_count;
    const int total = substream->self_total_count;
    const int max_cell_size = r->gpu.grav_max_cell_size;

    if (nslots == 0) continue;
    
    #ifdef SWIFT_DEBUG_TASKS

  ticks gpu_anchor_tic = 0;
  ticks gpu_sync_toc = 0;

  /* Here l already is the substream index. */
  const int gpu_substream_id = l;

#endif

    #ifdef SWIFT_GPU_TIMING

  double h2d_s = 0.0;
  double kernel_s = 0.0;
  double d2h_s = 0.0;
  double unpack_s = 0.0;

#endif

#if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)

  GPUEvent h2d_start, h2d_stop;
  GPUEvent kernel_start, kernel_stop;
  GPUEvent d2h_start, d2h_stop;

  GPUEventCreate(&h2d_start);
  GPUEventCreate(&h2d_stop);

  GPUEventCreate(&kernel_start);
  GPUEventCreate(&kernel_stop);

  GPUEventCreate(&d2h_start);
  GPUEventCreate(&d2h_stop);

#endif

    /* ===================== H2D ===================== */

    #ifdef SWIFT_DEBUG_TASKS
  gpu_anchor_tic = getticks();
#endif

#if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
  GPUEventRecord(h2d_start, substream->stream);
#endif

    GPUMemcpyAsync(substream->self_counts_d, substream->self_counts_h,
                   nslots * sizeof(int),
                   GPU_MEMCPY_HOST_TO_DEVICE, substream->stream);

    GPUMemcpyAsync(substream->self_offsets_d, substream->self_offsets_h,
                   nslots * sizeof(int),
                   GPU_MEMCPY_HOST_TO_DEVICE, substream->stream);

    GPUMemcpyAsync(substream->self_active_counts_d, substream->self_active_counts_h,
                   nslots * sizeof(int),
                   GPU_MEMCPY_HOST_TO_DEVICE, substream->stream);

    GPUMemcpyAsync(substream->self_active_offsets_d,
               substream->self_active_offsets_h,
               nslots * sizeof(int),
               GPU_MEMCPY_HOST_TO_DEVICE, substream->stream);

    GPUMemcpyAsync(substream->self_active_index_d,
               substream->self_active_index_h,
               (size_t)substream->self_total_active_count * sizeof(int),
               GPU_MEMCPY_HOST_TO_DEVICE, substream->stream);

    GPUMemcpyAsync(
	    substream->self_cell_flags_d,
	    substream->self_cell_flags_h,
	    nslots * sizeof(int),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

	GPUMemcpyAsync(
	    substream->self_use_full_d,
	    substream->self_use_full_h,
	    nslots * sizeof(int),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

	GPUMemcpyAsync(
	    substream->send_self_pos_mass_d,
	    substream->send_self_pos_mass,
	    total * sizeof(float4),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

	GPUMemcpyAsync(
	    substream->send_self_h_d,
	    substream->send_self_h,
	    total * sizeof(float),
	    GPU_MEMCPY_HOST_TO_DEVICE,
	    substream->stream);

    GPUMemsetAsync(
        substream->recv_self_active_d,
        0,
        (size_t)substream->self_total_active_count *
            sizeof(struct gravity_gpu_values_recv),
        substream->stream);

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS) 
    GPUEventRecord(h2d_stop, substream->stream);
    #endif

    /* ===================== KERNEL ===================== */

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS) 
    GPUEventRecord(kernel_start, substream->stream);
    #endif 

    runner_doself_grav_pp_flush(
        r, substream, nslots, substream->self_max_active_count, max_cell_size,
        substream->stream);

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)
    GPUEventRecord(kernel_stop, substream->stream);
    #endif 

    /* ===================== D2H ===================== */

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS) 
    GPUEventRecord(d2h_start, substream->stream);
    #endif 

    GPUMemcpyAsync(
        substream->recv_self_active,
        substream->recv_self_active_d,
        (size_t)substream->self_total_active_count *
            sizeof(struct gravity_gpu_values_recv),
        GPU_MEMCPY_DEVICE_TO_HOST,
        substream->stream);

    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS) 
    GPUEventRecord(d2h_stop, substream->stream);
    #endif

    GPUEventRecord(substream->done, substream->stream);
    #ifdef SWIFT_USE_TIMERS
    const ticks gpu_wait_tic = getticks();
    #endif

    GPUEventSynchronize(substream->done);

    #ifdef SWIFT_USE_TIMERS
    timers_toc(timer_gpu_grav_wait, gpu_wait_tic);
    #endif
    
    #ifdef SWIFT_DEBUG_TASKS
  gpu_sync_toc = getticks();
#endif

    #ifdef SWIFT_GPU_TIMING 
    h2d_s = runner_gpu_event_elapsed_s(h2d_start, h2d_stop);
    kernel_s = runner_gpu_event_elapsed_s(kernel_start, kernel_stop);
    d2h_s = runner_gpu_event_elapsed_s(d2h_start, d2h_stop);
    #endif
    
    #ifdef SWIFT_DEBUG_TASKS

  const double gpu_h2d_end_s =
      runner_gpu_event_offset_s(h2d_start, h2d_stop);

  const double gpu_kernel_start_s =
      runner_gpu_event_offset_s(h2d_start, kernel_start);

  const double gpu_kernel_end_s =
      runner_gpu_event_offset_s(h2d_start, kernel_stop);

  const double gpu_d2h_start_s =
      runner_gpu_event_offset_s(h2d_start, d2h_start);

  const double gpu_d2h_end_s =
      runner_gpu_event_offset_s(h2d_start, d2h_stop);

  runner_gpu_write_timeline_row(
      "self",
      (long long)r->e->step,
      r->id,
      gpu_substream_id,
      gpu_anchor_tic,
      gpu_sync_toc,
      gpu_h2d_end_s,
      gpu_kernel_start_s,
      gpu_kernel_end_s,
      gpu_d2h_start_s,
      gpu_d2h_end_s);

#endif

    /* ===================== UNPACK ===================== */

    #ifdef SWIFT_GPU_TIMING 
    const double unpack_t0 = runner_gpu_walltime_s();
    #endif 

    for (int j = 0; j < nslots; j++) {
      struct cell *c_unpack = substream->grav_cells_self[j];
      
      if (c_unpack == NULL)
        error("GPU self unpack received NULL cell.");

      if (c_unpack->nodeID != r->e->nodeID)
        error("GPU self unpack attempted to write a foreign cell.");

      if (!cell_is_active_gravity(c_unpack, r->e))
        continue;

      while (cell_glocktree(c_unpack)) {
        ;
      }

      const int active_count = substream->self_active_counts_h[j];
      const int active_base = substream->self_active_offsets_h[j];

      for (int a = 0; a < active_count; a++) {
        const int local_pid = substream->self_active_index_h[active_base + a];
        const int k = active_base + a;

        c_unpack->grav.parts[local_pid].a_grav[0] +=
            substream->recv_self_active[k].values_i.x;
        c_unpack->grav.parts[local_pid].a_grav[1] +=
            substream->recv_self_active[k].values_i.y;
        c_unpack->grav.parts[local_pid].a_grav[2] +=
            substream->recv_self_active[k].values_i.z;
        c_unpack->grav.parts[local_pid].potential +=
            substream->recv_self_active[k].values_i.w;
      }

      cell_gunlocktree(c_unpack);
    }

    #ifdef SWIFT_GPU_TIMING 
    unpack_s = runner_gpu_walltime_s() - unpack_t0;

    const double self_pack_s = runner_gpu_self_pack_time_s;
    runner_gpu_self_pack_time_s = 0.0;

    const size_t h2d_bytes =
    runner_gpu_self_h2d_bytes(
        substream,
        nslots);
        
    runner_gpu_write_timing_row(
	    "self",
	    (long long)r->e->step,
	    r->id,
	    r->qid,
	    nslots,
	    total,
	    h2d_bytes,
	    self_pack_s,
	    h2d_s,
	    kernel_s,
	    d2h_s,
	    unpack_s);
    #endif 
    
    #if defined(SWIFT_GPU_TIMING) || defined(SWIFT_DEBUG_TASKS)

  GPUEventDestroy(h2d_start);
  GPUEventDestroy(h2d_stop);

  GPUEventDestroy(kernel_start);
  GPUEventDestroy(kernel_stop);

  GPUEventDestroy(d2h_start);
  GPUEventDestroy(d2h_stop);

#endif

    runner_gpu_complete_self_batch(r, &r->e->sched, substream, NULL);

    result = flushed_self_task;
  }

  return result;
}

/**
 * @brief Flush any leftover packed pair-gravity work owned by a runner.
 *
 * @param r The runner whose GPU batch should be flushed.
 * @return The outcome of the leftover flush attempt.
 */
enum runner_gpu_task_type runner_gpu_flush_leftover_pair(struct runner *r) {

  runner_gpu_bind_device(r);

  enum runner_gpu_task_type result = regular_task;

for (int l = 0; l < r->gpu.nstreams; l++) {
  struct gpu_runner_substream *ss = &r->gpu.substreams[l];

    if (ss->grav_batch_pair_count == 0) continue;

    runner_dopair_grav_pp_flush(
        r, ss,
        ss->grav_cells_pair,
        ss->grav_tasks_pair,
        NULL,
        r->gpu.grav_batch_ncells,
        r->gpu.grav_max_cell_size,
        ss->stream);

    result = flushed_pair_task;
  }

  return result;
}

/**
 * @brief Computes the interaction of all the particles in a cell.
 *
 * This function will try to recurse as far down the tree as possible and only
 * default to direct summation if there is no better option.
 *
 * @param r The #runner.
 * @param c The first #cell.
 * @param gettimer Are we timing this ?
 */
enum runner_gpu_task_type runner_doself_recursive_grav_gpu(
    struct runner *r,
    struct gpu_runner_substream *substream,
    struct cell *c,
    const int gettimer,
    struct cell **grav_cells_self,
    struct task **grav_tasks_self,
    struct task *t,
    int ncells,
    int max_cell_size,
    GPUStream stream) {

  const struct engine *e = r->e;

  runner_clear_grav_flags(c, e);

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.count == 0) error("Doing self gravity on an empty cell !");
#endif

  TIMER_TIC;

  if (!cell_is_active_gravity(c, e)) {
    if (gettimer) TIMER_TOC(timer_dosub_self_grav);
    return regular_task;
  }

  enum runner_gpu_task_type task_type = regular_task;

  if (c->split) {

    for (int j = 0; j < 8; j++) {
      if (c->progeny[j] == NULL) continue;

      enum runner_gpu_task_type child_self_type =
          runner_doself_recursive_grav_gpu(
		    r,
		    substream,
		    c->progeny[j],
		    0,
		    grav_cells_self,
		    grav_tasks_self,
		    t,
		    ncells,
		    max_cell_size,
		    stream);

      if (child_self_type > task_type) task_type = child_self_type;

      for (int k = j + 1; k < 8; k++) {
        if (c->progeny[k] == NULL) continue;

        enum runner_gpu_task_type child_pair_type =
    	runner_dopair_recursive_grav_gpu(
        	r, substream, c->progeny[j], c->progeny[k], 0,
        	substream->grav_cells_pair, substream->grav_tasks_pair,
        	substream->grav_pair_internal_from_self,
        	t, 1, ncells, max_cell_size, substream->stream);

        if (child_pair_type > task_type) task_type = child_pair_type;
      }
    }

  } else {

    /* Leaf self cell: pack it */
    task_type = runner_doself_grav_pp_task_gpu(
        r, substream, c, t, ncells, max_cell_size);
  }

  /*
	 * A self-gravity scheduler task can generate pair interactions while
	 * recursively walking its children. These internal pair interactions belong
	 * to the parent self task and must be completed before that task can be
	 * released to the scheduler.
	 */
	if (gettimer && substream->grav_batch_pair_count > 0) {

	  runner_dopair_grav_pp_flush(
	      r,
	      substream,
	      substream->grav_cells_pair,
	      substream->grav_tasks_pair,
	      t,  /* current parent self task */
	      ncells,
	      max_cell_size,
	      substream->stream);
	}
      

    enum runner_gpu_task_type final_type = regular_task;

  if (substream->grav_batch_self_count > 0) {
  final_type = packed_task;
  } else if (task_type != regular_task) {
  final_type = flushed_self_task;
  } else {
  final_type = regular_task;
  }
  
  #ifdef SWIFT_DEBUG_CHECKS
  if (gettimer && final_type == packed_task) {
    if (!t->gpu_counted || t->gpu_completed) {
      error("Self task deferred without valid GPU bookkeeping.");
    }
  }
  #endif

  if (gettimer) TIMER_TOC(timer_dosub_self_grav);
  return final_type;
}
