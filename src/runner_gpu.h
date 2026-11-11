/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_RUNNER_GPU_H
#define SWIFT_RUNNER_GPU_H

#include "gpu_mapping.h"
#include "gpu_pack_params.h"
#include "cuda/gpu_offload_data.h"

struct cell;
struct engine;
struct runner;
struct scheduler;
struct task;
struct gravity_gpu_values_recv;
struct gravity_gpu_values_send;

/**
 * @brief Enumeration of the types of operation a GPU task can have performed.
 *
 * `regular_task` means the task completed normally and should be finished with
 * the standard scheduler path.
 * `packed_task` means the task completed packing into a GPU batch but the batch
 * was not flushed yet.
 * `flushed_self_task` means a self-gravity GPU batch was flushed.
 * `flushed_pair_task` means a pair-gravity GPU batch was flushed.
 */
enum runner_gpu_task_type {
  regular_task = 0,
  packed_task = 1,
  flushed_self_task = 2,
  flushed_pair_task = 3
};

/**
 * @brief GPU pair-work state owned by a single runner substream.
 *
 * Each substream corresponds to one explicit GPU stream and owns one in-flight
 * pair batch.
 */
struct gpu_runner_substream {

  /*! Stream used for this substream's GPU pair work. */
  GPUStream stream;

  /*! Completion event for this substream's current async work. */
  GPUEvent done;

  /*! Whether this substream currently owns in-flight work. */
  int busy;
  
  /*SELF OPERATIONS*/
  /*! Number of self cells currently packed in this substream's GPU batch. */
  int grav_batch_self_count;
  
  /*! Host and device buffers for sent self interactions. */
  float4 *send_self_pos_mass, *send_self_pos_mass_d;
  float *send_self_h, *send_self_h_d;

  /*! Host and device buffers for recieved self results. */
  struct gravity_gpu_values_recv *recv_self, *recv_self_d;
  struct gravity_gpu_values_recv *recv_self_active, *recv_self_active_d;

  /*! Packed self-batch cell and task handles. */
  struct cell **grav_cells_self;
  struct task **grav_tasks_self;
  
  /*! Counts for self operations for cell packing*/
  int *self_counts_h, *self_counts_d;
  int *self_offsets_h, *self_offsets_d;
  int self_total_count;
  
  /*! Active-target compaction metadata for self work. */
  int *self_active_counts_h, *self_active_counts_d;
  int *self_active_offsets_h, *self_active_offsets_d;
  int *self_active_index_h, *self_active_index_d;
  int self_total_active_count;
  int self_max_active_count;
  
  int *self_cell_flags_h, *self_cell_flags_d;
  int *self_use_full_h, *self_use_full_d;
  
  /*! rmax values */
  float *self_rmax_h, *self_rmax_d;
  
  /*PAIR OPERATIONS*/
  /*! Number of pair cells currently packed in this substream's GPU batch. */
  int grav_batch_pair_count;
  
  int pair_unique_cell_count;
  struct cell **pair_unique_cells;

  int *pair_pair_i_h, *pair_pair_i_d;
  int *pair_pair_j_h, *pair_pair_j_d;
  
  int *pair_use_full_h, *pair_use_full_d;

  /* size = 2 * max_pairs */
  int *pair_side_active_offsets_h, *pair_side_active_offsets_d;
  int pair_total_pair_active_count;
  /*-----*/
  
  /*! Compact host and device buffers for coalesced pair source-tile loads. */
  float4 *send_pair_pos_mass;
  float4 *send_pair_pos_mass_d;

  float *send_pair_h;
  float *send_pair_h_d;

  /*! Host and device buffers for recieved pair results. */
  struct gravity_gpu_values_recv *recv_pair_active, *recv_pair_active_d;

  /*! Packed pair-batch cell and task handles. */
  struct cell **grav_cells_pair;
  struct task **grav_tasks_pair;
  
  /* Flag for pairs spawned from self recursion */
  unsigned char *grav_pair_internal_from_self;
  
  /*! Counts for pair operations for cell packing*/
  int *pair_counts_h, *pair_counts_d;
  int *pair_offsets_h, *pair_offsets_d;
  int pair_total_count;
  
  int *pair_cell_flags_h, *pair_cell_flags_d;
  
  /*! Active-target compaction metadata for pair work. */
  int *pair_active_counts_h, *pair_active_counts_d;
  int *pair_active_offsets_h, *pair_active_offsets_d;
  int *pair_active_index_h, *pair_active_index_d;
  int pair_total_active_count;
  int pair_max_active_count;
};

struct gpu_hydro_state {
    /*hydro packing params*/
    struct gpu_global_pack_params params;

    /* hydro gpu buffers */
    struct gpu_offload_data density;
    struct gpu_offload_data gradient;
    struct gpu_offload_data force;
    
    /* hydro streams */
    GPUStream *streams;
    int nstreams;

    int initialised;
};


/**
 * @brief GPU-specific state owned by a single runner.
 */
struct gpu_runner {
  /*! Substreams for explicit multi-stream execution. */
  /* Gravity-specific */
  struct gpu_runner_substream *substreams;
  
  /*! Number of GPU streams launched per runner. */
  int nstreams;

  /*! Number of cells that fit in one GPU batch for this runner. */
  int grav_batch_ncells;

  /*! Maximum number of particles packed per cell. */
  int grav_max_cell_size;

  /*! Next substream index to try when acquiring a pair-work substream. */
  int next_substream;
  
  /*! CUDA/HIP device selected by this MPI rank. */
  int device_id;

  /*! Local MPI rank within the node. */
  int local_mpi_rank;

  /*! Number of MPI ranks on this node. */
  int local_mpi_size;

  /*! Number of local MPI ranks sharing this selected GPU. */
  int local_ranks_on_device;
  
  /*Hydro*/
  struct gpu_hydro_state hydro;
};

void runner_gpu_hydro_init(struct runner *r);
void runner_gpu_hydro_clean(struct runner *r);

struct gpu_runner_substream *runner_gpu_acquire_substream(struct runner *r);

/**
 * @brief Initialise GPU parameters from the parameter file.
 *
 * @param e The #engine to unpack parameters for.
 */
void runner_gpu_params_init(struct engine *e);

/**
 * @brief Initialise the GPU-specific state attached to a runner.
 *
 * @param r The runner whose GPU state to initialise.
 */
void runner_gpu_init(struct runner *r);

/**
 * @brief Clean the GPU-specific state attached to a runner.
 *
 * @param r The runner whose GPU state to clean.
 */
void runner_gpu_clean(struct runner *r);

/**
 * @brief Pack, launch, and unpack a batched self-gravity GPU task.
 *
 * @param r The #runner.
 * @param c The #cell to pack.
 * @param t The #task being executed.
 * @param ncells The batch capacity in cells.
 * @param max_cell_size The maximum number of particles per packed cell.
 * @return The outcome of the GPU wrapper for this task.
 */
enum runner_gpu_task_type runner_doself_grav_pp_task_gpu(
    struct runner *r, struct gpu_runner_substream *substream,
    struct cell *c, struct task *t, int ncells, int max_cell_size);

/**
 * @brief Pack one leaf pair-gravity interaction into the substream GPU batch.
 *
 * @return The outcome of the GPU wrapper for this task.
 */
enum runner_gpu_task_type runner_dopair_grav_pp_gpu(
    struct runner *r, struct gpu_runner_substream *substream, struct cell *ci,
    struct cell *cj, const int symmetric, const int allow_mpole,
    struct cell **grav_cells_pair, struct task **grav_tasks_pair,
    unsigned char *grav_pair_internal_from_self,
    struct task *t, int internal_from_self,
    int ncells, int max_cell_size, GPUStream stream);

/**
 * @brief Recursively process a pair-gravity task with GPU batching.
 *
 * @return The outcome of the GPU wrapper for this task.
 */
enum runner_gpu_task_type runner_dopair_recursive_grav_gpu(
    struct runner *r, struct gpu_runner_substream *substream, struct cell *ci,
    struct cell *cj, const int gettimer,
    struct cell **grav_cells_pair, struct task **grav_tasks_pair,
    unsigned char *grav_pair_internal_from_self,
    struct task *t, int internal_from_self,
    int ncells, int max_cell_size, GPUStream stream);
    
    
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
    GPUStream stream);

/**
 * @brief Flush any leftover packed self-gravity work owned by a runner.
 *
 * @param r The runner whose GPU batch should be flushed.
 * @return The outcome of the leftover flush attempt.
 */
enum runner_gpu_task_type runner_gpu_flush_leftover_self(struct runner *r);

/**
 * @brief Flush any leftover packed pair-gravity work owned by all substreams of a
 * runner.
 *
 * @param r The runner whose GPU pair batches should be flushed.
 * @return The outcome of the leftover flush attempt.
 */
enum runner_gpu_task_type runner_gpu_flush_leftover_pair(struct runner *r);

/**
 * @brief Complete all self tasks currently stored in the runner GPU batch.
 */
void runner_gpu_complete_self_batch(struct runner *r, struct scheduler *sched,
                                    struct gpu_runner_substream *substream, struct task *current_task);

/**
 * @brief Complete a single pair task (decrement queue counter and mark done).
 */
void runner_gpu_complete_pair_task(struct runner *r, struct scheduler *sched,
                                   struct task *t);
                                    
void runner_gpu_bind_device(struct runner *r);

void runner_gpu_complete_current_self_task(struct runner *r,
                                           struct scheduler *sched,
                                           struct task *t);
                                           
/**
 * @brief Count a hydro scheduler task offloaded to the GPU.
 */
void runner_gpu_count_hydro_task(struct runner *r,
                                 struct scheduler *sched,
                                 struct task *t);
                                 
void runner_gpu_complete_hydro_task(struct runner *r,
                                    struct scheduler *sched,
                                    struct task *t);                          

#endif /* SWIFT_RUNNER_GPU_H */
