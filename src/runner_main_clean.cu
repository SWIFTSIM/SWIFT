/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#define GPUOFFLOAD_DENSITY 1   // off-load hydro density to GPU
#define GPUOFFLOAD_GRADIENT 1  // off-load hydro gradient to GPU
#define GPUOFFLOAD_FORCE 1     // off-load hydro force to GPU

// #define DO_CORNERS 1 //do corner pair tasks on CPU
// #define DUMP_TIMINGS 1
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "engine.h"
#include "feedback.h"
#include "runner_doiact_sinks.h"
#include "scheduler.h"
#include "space_getsid.h"
#include "timers.h"

/* Import the gravity loop functions. */
#include "runner_doiact_grav.h"

/* Import the density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_hydro.h"
#include "runner_doiact_undef.h"

/* Import the gradient loop functions (if required). */
#ifdef EXTRA_HYDRO_LOOP
#define FUNCTION gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_GRADIENT
#include "runner_doiact_hydro.h"
#include "runner_doiact_undef.h"
#endif

/* Import the force loop functions. */
#define FUNCTION force
#define FUNCTION_TASK_LOOP TASK_LOOP_FORCE
#include "runner_doiact_hydro.h"
#include "runner_doiact_undef.h"

/* Import the limiter loop functions. */
#define FUNCTION limiter
#define FUNCTION_TASK_LOOP TASK_LOOP_LIMITER
#include "runner_doiact_limiter.h"
#include "runner_doiact_undef.h"

/* Import the stars density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_stars.h"
#include "runner_doiact_undef.h"

#ifdef EXTRA_STAR_LOOPS

/* Import the stars prepare1 loop functions. */
#define FUNCTION prep1
#define FUNCTION_TASK_LOOP TASK_LOOP_STARS_PREP1
#include "runner_doiact_stars.h"
#include "runner_doiact_undef.h"

/* Import the stars prepare2 loop functions. */
#define FUNCTION prep2
#define FUNCTION_TASK_LOOP TASK_LOOP_STARS_PREP2
#include "runner_doiact_stars.h"
#include "runner_doiact_undef.h"

#endif /* EXTRA_STAR_LOOPS */

/* Import the stars feedback loop functions. */
#define FUNCTION feedback
#define FUNCTION_TASK_LOOP TASK_LOOP_FEEDBACK
#include "runner_doiact_stars.h"
#include "runner_doiact_undef.h"

/* Import the black hole density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_black_holes.h"
#include "runner_doiact_undef.h"

/* Import the black hole feedback loop functions. */
#define FUNCTION swallow
#define FUNCTION_TASK_LOOP TASK_LOOP_SWALLOW
#include "runner_doiact_black_holes.h"
#include "runner_doiact_undef.h"

/* Import the black hole feedback loop functions. */
#define FUNCTION feedback
#define FUNCTION_TASK_LOOP TASK_LOOP_FEEDBACK
#include "runner_doiact_black_holes.h"
#include "runner_doiact_undef.h"

/* Import the RT gradient loop functions */
#define FUNCTION rt_gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_GRADIENT
#include "runner_doiact_hydro.h"
#include "runner_doiact_undef.h"

/* Import the RT transport (force) loop functions. */
#define FUNCTION rt_transport
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_TRANSPORT
#include "runner_doiact_hydro.h"
#include "runner_doiact_undef.h"

#ifdef __cplusplus
}
#endif
/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 **/

/* CUDA Header */
#ifdef WITH_CUDA
#ifdef __cplusplus
extern "C" {
#endif

#include "cuda/part_gpu.h"
#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include "runner_doiact_functions_hydro_gpu.h"
#include "runner_gpu_pack_functions.h"
#include "cuda/GPU_runner_functions.h"

#ifdef __cplusplus
}
#endif
// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
#define CUDA_DEBUG

inline cudaError_t checkCuda(cudaError_t result) {
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}

void *runner_main2(void *data) {
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  struct space *space = e->s;

  /*pack_vars contain data required for packing tasks destined for the GPU*/
  struct pack_vars_self *pack_vars_self_dens;
  struct pack_vars_self *pack_vars_self_forc;
  struct pack_vars_self *pack_vars_self_grad;

  /*pack_vars contain data required for packing tasks destined for the GPU*/
  struct pack_vars_pair *pack_vars_pair_dens;
  struct pack_vars_pair *pack_vars_pair_forc;
  struct pack_vars_pair *pack_vars_pair_grad;

  cudaMallocHost((void **)&pack_vars_self_dens,
                 sizeof(struct pack_vars_self *));
  cudaMallocHost((void **)&pack_vars_self_forc,
                 sizeof(struct pack_vars_self *));
  cudaMallocHost((void **)&pack_vars_self_grad,
                 sizeof(struct pack_vars_self *));

  cudaMallocHost((void **)&pack_vars_pair_dens,
                 sizeof(struct pack_vars_pair *));
  cudaMallocHost((void **)&pack_vars_pair_forc,
                 sizeof(struct pack_vars_pair *));
  cudaMallocHost((void **)&pack_vars_pair_grad,
                 sizeof(struct pack_vars_pair *));

  int devId = 0;  // find and print gpu device name
  struct cudaDeviceProp prop;
  int nDevices;
  int maxBlocksSM;
  int nSMs;
  cudaGetDeviceCount(&nDevices);
  cudaGetDeviceProperties(&prop, devId);
  cudaDeviceGetAttribute(&maxBlocksSM, cudaDevAttrMaxBlocksPerMultiprocessor,
                         devId);
  cudaDeviceGetAttribute(&nSMs, cudaDevAttrMultiProcessorCount, devId);
  int nPartsPerCell = space->nr_parts / space->tot_cells;
  int mpi_rank = 0;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
  if (r->cpuid == 0 && mpi_rank == 0) {
    fprintf(stderr, "%i devices available device id is %i\n", nDevices, devId);
    fprintf(stderr, "Device : %s\n", prop.name);
    fprintf(stderr, "nSMs %i max blocks per SM %i maxnBlocks per stream %i\n",
            nSMs, maxBlocksSM, nSMs * maxBlocksSM);
    fprintf(stderr, "Target nBlocks per kernel is %i\n",
            N_TASKS_BUNDLE_SELF * nPartsPerCell / BLOCK_SIZE);
    fprintf(stderr, "Target nBlocks per stream is %i\n",
            N_TASKS_PER_PACK_SELF * nPartsPerCell / BLOCK_SIZE);
  }
  if (nDevices == 1) cudaSetDevice(devId);
#ifndef WITH_MPI
  else {
    cudaSetDevice(devId);
  }
#endif
#ifdef WITH_MPI
  else {
    cudaSetDevice(mpi_rank);
    fprintf(stderr, "%i devices available device id is %i\n", nDevices,
            mpi_rank);
  }
#endif
  fprintf(stderr, "after dev select engine_rank %i rank %i\n", engine_rank,
          mpi_rank);

  cudaError_t cu_error;
  size_t free_mem, total_mem;
  cudaMemGetInfo(&free_mem, &total_mem);

  fprintf(stderr, "free mem %lu, total mem %lu\n", free_mem, total_mem);
  // how many tasks do we want for each launch of GPU kernel
  const int target_n_tasks = sched->pack_size;
  const int target_n_tasks_pair = sched->pack_size_pair;
  pack_vars_self_dens->target_n_tasks = target_n_tasks;
  pack_vars_pair_dens->target_n_tasks = target_n_tasks_pair;
  pack_vars_self_forc->target_n_tasks = target_n_tasks;
  pack_vars_pair_forc->target_n_tasks = target_n_tasks_pair;
  pack_vars_self_grad->target_n_tasks = target_n_tasks;
  pack_vars_pair_grad->target_n_tasks = target_n_tasks_pair;
  // how many tasks we want in each bundle (used for launching kernels in
  // different streams)
  const int bundle_size = N_TASKS_BUNDLE_SELF;
  const int bundle_size_pair = N_TASKS_BUNDLE_PAIR;

  pack_vars_self_dens->bundle_size = bundle_size;
  pack_vars_pair_dens->bundle_size = bundle_size_pair;
  pack_vars_self_forc->bundle_size = bundle_size;
  pack_vars_pair_forc->bundle_size = bundle_size_pair;
  pack_vars_self_grad->bundle_size = bundle_size;
  pack_vars_pair_grad->bundle_size = bundle_size_pair;
  // Keep track of first and last particles for each task (particle data is
  // arranged in long arrays containing particles from all the tasks we will
  // work with)

  // Copy of the above residing on the GPU
  int *d_task_first_part_self_dens, *d_task_last_part_self_dens;
  int2 *task_first_part_self_dens_f4;
  int2 *task_first_part_f4;
  int2 *task_first_part_f4_f;
  int2 *task_first_part_f4_g;
  int2 *d_task_first_part_f4;
  int2 *d_task_first_part_f4_f;
  int2 *d_task_first_part_f4_g;
  int *d_task_first_part_self_forc, *d_task_last_part_self_forc;
  int *d_task_first_part_self_grad, *d_task_last_part_self_grad;
  int *d_task_first_parts_pair_dens, *d_task_last_parts_pair_dens;

  int4 *fparti_fpartj_lparti_lpartj_dens;
  int4 *fparti_fpartj_lparti_lpartj_forc, *d_fparti_fpartj_lparti_lpartj_forc;
  int4 *fparti_fpartj_lparti_lpartj_grad, *d_fparti_fpartj_lparti_lpartj_grad;

  int *d_task_first_parts_pair_forc, *d_task_last_parts_pair_forc;
  int *d_task_first_parts_pair_grad, *d_task_last_parts_pair_grad;

  cudaMallocManaged((void **)&task_first_part_self_dens_f4,
                    target_n_tasks * sizeof(int2), cudaMemAttachGlobal);
  cudaMallocHost((void **)&task_first_part_f4, target_n_tasks * sizeof(int2));
  cudaMalloc((void **)&d_task_first_part_f4, target_n_tasks * sizeof(int2));
  cudaMallocHost((void **)&task_first_part_f4_f, target_n_tasks * sizeof(int2));
  cudaMalloc((void **)&d_task_first_part_f4_f, target_n_tasks * sizeof(int2));
  cudaMallocHost((void **)&task_first_part_f4_g, target_n_tasks * sizeof(int2));
  cudaMalloc((void **)&d_task_first_part_f4_g, target_n_tasks * sizeof(int2));

  cudaMallocHost((void **)&fparti_fpartj_lparti_lpartj_dens,
                 target_n_tasks * sizeof(int4));

  cudaMallocHost((void **)&fparti_fpartj_lparti_lpartj_forc,
                 target_n_tasks * sizeof(int4));
  cudaMalloc((void **)&d_fparti_fpartj_lparti_lpartj_forc,
             target_n_tasks * sizeof(int4));

  cudaMallocHost((void **)&fparti_fpartj_lparti_lpartj_grad,
                 target_n_tasks * sizeof(int4));
  cudaMalloc((void **)&d_fparti_fpartj_lparti_lpartj_grad,
             target_n_tasks * sizeof(int4));

  // Arrays keeping track of the row numbers of the first and last particles
  // within each bundle. Required by the GPU code

  cudaMallocHost((void **)&pack_vars_self_dens->task_first_part,
                 target_n_tasks * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_dens->task_last_part,
                 target_n_tasks * sizeof(int));

  cudaMallocHost((void **)&pack_vars_pair_dens->task_first_part,
                 2 * target_n_tasks * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_dens->task_last_part,
                 2 * target_n_tasks * sizeof(int));

  cudaMallocHost((void **)&pack_vars_self_forc->task_first_part,
                 target_n_tasks * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_forc->task_last_part,
                 target_n_tasks * sizeof(int));

  cudaMallocHost((void **)&pack_vars_pair_forc->task_first_part,
                 2 * target_n_tasks * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_forc->task_last_part,
                 2 * target_n_tasks * sizeof(int));

  cudaMallocHost((void **)&pack_vars_self_grad->task_first_part,
                 target_n_tasks * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_grad->task_last_part,
                 target_n_tasks * sizeof(int));

  cudaMallocHost((void **)&pack_vars_pair_grad->task_first_part,
                 2 * target_n_tasks * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_grad->task_last_part,
                 2 * target_n_tasks * sizeof(int));

  /* nBundles is the number of task bundles each
  thread has ==> Used to loop through bundles */
  int nBundles = (target_n_tasks + bundle_size - 1) / bundle_size;
  int nBundles_pair =
      (target_n_tasks_pair + bundle_size_pair - 1) / bundle_size_pair;

  if (r->cpuid == 0) {
    fprintf(stderr, "engine_rank %i cpuid %i nBundles/nStreams %i\n",
            engine_rank, r->cpuid, nBundles);
    fprintf(stderr, "nBundles/nStreams Pair %i\n", nBundles_pair);
  }

  pack_vars_self_dens->nBundles = nBundles;
  pack_vars_pair_dens->nBundles = nBundles_pair;
  pack_vars_self_forc->nBundles = nBundles;
  pack_vars_pair_forc->nBundles = nBundles_pair;
  pack_vars_self_grad->nBundles = nBundles;
  pack_vars_pair_grad->nBundles = nBundles_pair;

  // first part and last part are the first and last particle ids (locally
  // within this thread)

  cudaMallocHost((void **)&pack_vars_self_dens->bundle_first_part,
                 nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_dens->bundle_last_part,
                 nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_dens->bundle_first_task_list,
                 nBundles * sizeof(int));

  cudaMallocHost((void **)&pack_vars_pair_dens->bundle_first_part,
                 2 * nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_dens->bundle_last_part,
                 2 * nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_dens->bundle_first_task_list,
                 2 * nBundles * sizeof(int));

  cudaMallocHost((void **)&pack_vars_self_forc->bundle_first_part,
                 nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_forc->bundle_last_part,
                 nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_forc->bundle_first_task_list,
                 nBundles * sizeof(int));

  cudaMallocHost((void **)&pack_vars_pair_forc->bundle_first_part,
                 2 * nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_forc->bundle_last_part,
                 2 * nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_forc->bundle_first_task_list,
                 2 * nBundles * sizeof(int));

  cudaMallocHost((void **)&pack_vars_self_grad->bundle_first_part,
                 nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_grad->bundle_last_part,
                 nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_self_grad->bundle_first_task_list,
                 nBundles * sizeof(int));

  cudaMallocHost((void **)&pack_vars_pair_grad->bundle_first_part,
                 2 * nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_grad->bundle_last_part,
                 2 * nBundles * sizeof(int));
  cudaMallocHost((void **)&pack_vars_pair_grad->bundle_first_task_list,
                 2 * nBundles * sizeof(int));

  // These I need to keep/////////////////
  cudaMalloc((void **)&d_task_first_part_self_dens,
             target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_first_part_self_forc,
             target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_first_part_self_grad,
             target_n_tasks * sizeof(int));

  cudaMalloc((void **)&d_task_last_part_self_dens,
             target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_last_part_self_forc,
             target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_last_part_self_grad,
             target_n_tasks * sizeof(int));
  // These I need to keep/////////////////
  pack_vars_self_dens->d_task_first_part = d_task_first_part_self_dens;
  pack_vars_self_dens->d_task_last_part = d_task_last_part_self_dens;
  // These I need to keep/////////////////
  pack_vars_self_forc->d_task_first_part = d_task_first_part_self_forc;
  pack_vars_self_forc->d_task_last_part = d_task_last_part_self_forc;
  // These I need to keep/////////////////
  pack_vars_self_grad->d_task_first_part = d_task_first_part_self_grad;
  pack_vars_self_grad->d_task_last_part = d_task_last_part_self_grad;

  // These I need to keep/////////////////
  cudaMalloc((void **)&d_task_first_parts_pair_dens,
             2 * target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_first_parts_pair_forc,
             2 * target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_first_parts_pair_grad,
             2 * target_n_tasks * sizeof(int));

  cudaMalloc((void **)&d_task_last_parts_pair_dens,
             2 * target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_last_parts_pair_forc,
             2 * target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_last_parts_pair_grad,
             2 * target_n_tasks * sizeof(int));
  // These I need to keep/////////////////
  pack_vars_pair_dens->d_task_first_part = d_task_first_parts_pair_dens;
  pack_vars_pair_dens->d_task_last_part = d_task_last_parts_pair_dens;
  pack_vars_pair_forc->d_task_first_part = d_task_first_parts_pair_forc;
  pack_vars_pair_forc->d_task_last_part = d_task_last_parts_pair_forc;
  pack_vars_pair_grad->d_task_first_part = d_task_first_parts_pair_grad;
  pack_vars_pair_grad->d_task_last_part = d_task_last_parts_pair_grad;
  // cell positions for self tasks REMEMBER to remove CPU copies as these are no
  // longer necessary
  double *d_dens_cell_x, *d_dens_cell_y, *d_dens_cell_z;
  float3 *d_dens_f3_cell_x;
  double *d_grad_cell_x, *d_grad_cell_y, *d_grad_cell_z;
  double *d_forc_cell_x, *d_forc_cell_y, *d_forc_cell_z;
  // Shifts for pair tasks REMEMBER to remove CPU copies as these are no longer
  // necessary
  double *d_dens_shift_x, *d_dens_shift_y, *d_dens_shift_z;
  double *d_grad_shift_x, *d_grad_shift_y, *d_grad_shift_z;
  double *d_forc_shift_x, *d_forc_shift_y, *d_forc_shift_z;

  // These I need to keep/////////////////
  cudaMalloc((void **)&d_dens_cell_x, target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_dens_cell_y, target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_dens_cell_z, target_n_tasks * sizeof(double));

  cudaMalloc((void **)&d_dens_f3_cell_x, target_n_tasks * sizeof(float3));

  cudaMalloc((void **)&d_forc_cell_x, target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_forc_cell_y, target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_forc_cell_z, target_n_tasks * sizeof(double));

  cudaMalloc((void **)&d_grad_cell_x, target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_grad_cell_y, target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_grad_cell_z, target_n_tasks * sizeof(double));

  cudaMalloc((void **)&d_dens_shift_x, 2 * target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_dens_shift_y, 2 * target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_dens_shift_z, 2 * target_n_tasks * sizeof(double));

  cudaMalloc((void **)&d_forc_shift_x, 2 * target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_forc_shift_y, 2 * target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_forc_shift_z, 2 * target_n_tasks * sizeof(double));

  cudaMalloc((void **)&d_grad_shift_x, 2 * target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_grad_shift_y, 2 * target_n_tasks * sizeof(double));
  cudaMalloc((void **)&d_grad_shift_z, 2 * target_n_tasks * sizeof(double));
  // These I need to keep/////////////////

  cudaMallocHost((void **)&pack_vars_self_dens->cellx,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost((void **)&pack_vars_self_dens->celly,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost((void **)&pack_vars_self_dens->cellz,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host

  pack_vars_self_dens->d_cellx = d_dens_cell_x;
  pack_vars_self_dens->d_celly = d_dens_cell_y;
  pack_vars_self_dens->d_cellz = d_dens_cell_z;

  cudaMallocHost(
      (void **)&pack_vars_pair_dens->shiftx,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost(
      (void **)&pack_vars_pair_dens->shifty,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost(
      (void **)&pack_vars_pair_dens->shiftz,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host

  pack_vars_pair_dens->d_shiftx = d_dens_shift_x;
  pack_vars_pair_dens->d_shifty = d_dens_shift_y;
  pack_vars_pair_dens->d_shiftz = d_dens_shift_z;

  cudaMallocHost((void **)&pack_vars_self_forc->cellx,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost((void **)&pack_vars_self_forc->celly,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost((void **)&pack_vars_self_forc->cellz,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host

  pack_vars_self_forc->d_cellx = d_forc_cell_x;
  pack_vars_self_forc->d_celly = d_forc_cell_y;
  pack_vars_self_forc->d_cellz = d_forc_cell_z;

  cudaMallocHost(
      (void **)&pack_vars_pair_forc->shiftx,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost(
      (void **)&pack_vars_pair_forc->shifty,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost(
      (void **)&pack_vars_pair_forc->shiftz,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host

  pack_vars_pair_forc->d_shiftx = d_forc_shift_x;
  pack_vars_pair_forc->d_shifty = d_forc_shift_y;
  pack_vars_pair_forc->d_shiftz = d_forc_shift_z;

  cudaMallocHost((void **)&pack_vars_self_grad->cellx,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost((void **)&pack_vars_self_grad->celly,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost((void **)&pack_vars_self_grad->cellz,
                 target_n_tasks * sizeof(double));  // Pinned allocation on host

  pack_vars_self_grad->d_cellx = d_grad_cell_x;
  pack_vars_self_grad->d_celly = d_grad_cell_y;
  pack_vars_self_grad->d_cellz = d_grad_cell_z;

  cudaMallocHost(
      (void **)&pack_vars_pair_grad->shiftx,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost(
      (void **)&pack_vars_pair_grad->shifty,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host
  cudaMallocHost(
      (void **)&pack_vars_pair_grad->shiftz,
      2 * target_n_tasks * sizeof(double));  // Pinned allocation on host

  pack_vars_pair_grad->d_shiftx = d_grad_shift_x;
  pack_vars_pair_grad->d_shifty = d_grad_shift_y;
  pack_vars_pair_grad->d_shiftz = d_grad_shift_z;

  cudaStream_t stream[nBundles];
  cudaStream_t stream_pairs[nBundles_pair];

  cudaEvent_t self_end[nBundles];
  for (int i = 0; i < nBundles; i++) cudaEventCreate(&self_end[i]);

  cudaEvent_t self_end_g[nBundles];
  for (int i = 0; i < nBundles; i++) cudaEventCreate(&self_end_g[i]);

  cudaEvent_t self_end_f[nBundles];
  for (int i = 0; i < nBundles; i++) cudaEventCreate(&self_end_f[i]);

  cudaEvent_t pair_end[nBundles];
  for (int i = 0; i < nBundles; i++) cudaEventCreate(&pair_end[i]);

  cudaEvent_t pair_end_g[nBundles];
  for (int i = 0; i < nBundles; i++) cudaEventCreate(&pair_end_g[i]);

  cudaEvent_t pair_end_f[nBundles];
  for (int i = 0; i < nBundles; i++) cudaEventCreate(&pair_end_f[i]);

  int tasksperbundle = (target_n_tasks + nBundles - 1) / nBundles;
  int tasksperbundle_pair =
      (target_n_tasks_pair + nBundles_pair - 1) / nBundles_pair;

  pack_vars_self_dens->tasksperbundle = tasksperbundle;
  pack_vars_pair_dens->tasksperbundle = tasksperbundle_pair;
  pack_vars_self_forc->tasksperbundle = tasksperbundle;
  pack_vars_pair_forc->tasksperbundle = tasksperbundle_pair;
  pack_vars_self_grad->tasksperbundle = tasksperbundle;
  pack_vars_pair_grad->tasksperbundle = tasksperbundle_pair;

  for (int i = 0; i < nBundles; ++i)
    cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking);

  for (int i = 0; i < nBundles_pair; ++i)
    cudaStreamCreateWithFlags(&stream_pairs[i], cudaStreamNonBlocking);

  pack_vars_self_dens->count_parts = 0;
  pack_vars_pair_dens->count_parts = 0;
  pack_vars_self_forc->count_parts = 0;
  pack_vars_pair_forc->count_parts = 0;
  pack_vars_self_grad->count_parts = 0;
  pack_vars_pair_grad->count_parts = 0;

  /*Estimate how many particles to pack for GPU for each GPU launch
   * instruction*/
  int nr_nodes = 1, res = 0;
#ifdef WITH_MPI
  if ((res = MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes)) != MPI_SUCCESS)
    error("MPI_Comm_size failed with error %i.", res);
#endif
  int parts_per_top_level_cell =
      space->nr_local_cells_with_particles /
      space->nr_parts; /*A. Nasar: What I think is a good approximation for
                                   average N particles in each top level cell*/
  float eta_neighbours = e->s->eta_neighbours;
  int np_per_cell = ceil(2.0 * eta_neighbours);
  np_per_cell *= np_per_cell * np_per_cell;
  /*A. Nasar: Increase parts per recursed task-level cell by buffer to
    ensure we allocate enough memory*/
  int buff = ceil(0.5 * np_per_cell);

  int tot_self_tasks = space->nr_parts / np_per_cell;

  /*A. Nasar: Multiplication by 2 is also to ensure we do not over-run
   *  the allocated memory on buffers and GPU. This can happen if calculated h
   * is larger than cell width and splitting makes bigger than target cells*/
  int count_max_parts_tmp = 64 * 8 * target_n_tasks * (np_per_cell + buff);

  //  message("np per cell %i, max_parts %i, n_tasks_GPU %i\n", np_per_cell,
  //  count_max_parts_tmp, target_n_tasks);
  pack_vars_self_dens->count_max_parts = count_max_parts_tmp;
  pack_vars_pair_dens->count_max_parts = count_max_parts_tmp;
  pack_vars_self_forc->count_max_parts = count_max_parts_tmp;
  pack_vars_pair_forc->count_max_parts = count_max_parts_tmp;
  pack_vars_self_grad->count_max_parts = count_max_parts_tmp;
  pack_vars_pair_grad->count_max_parts = count_max_parts_tmp;

  struct part_aos *parts_aos_dens;
  struct part_aos_f4 *parts_aos_dens_f4;
  struct part_aos_f4_send *parts_aos_f4_send;
  struct part_aos_f4_recv *parts_aos_f4_recv;

  struct part_aos_f *parts_aos_forc;
  struct part_aos_f4_f *parts_aos_forc_f4;
  struct part_aos_f4_f_send *parts_aos_forc_f4_send;
  struct part_aos_f4_f_recv *parts_aos_forc_f4_recv;

  struct part_aos_g *parts_aos_grad;
  struct part_aos_f4_g *parts_aos_grad_f4;
  struct part_aos_f4_g_send *parts_aos_grad_f4_send;
  struct part_aos_f4_g_recv *parts_aos_grad_f4_recv;

  struct part_aos *d_parts_aos_dens;
  struct part_aos_f4 *d_parts_aos_dens_f4;
  struct part_aos_f4_send *d_parts_aos_f4_send;
  struct part_aos_f4_recv *d_parts_aos_f4_recv;

  struct part_aos_f *d_parts_aos_forc;
  struct part_aos_f4_f *d_parts_aos_forc_f4;
  struct part_aos_f4_f_send *d_parts_aos_forc_f4_send;
  struct part_aos_f4_f_recv *d_parts_aos_forc_f4_recv;

  struct part_aos_g *d_parts_aos_grad;
  struct part_aos_f4_g *d_parts_aos_grad_f4;
  struct part_aos_f4_g_send *d_parts_aos_grad_f4_send;
  struct part_aos_f4_g_recv *d_parts_aos_grad_f4_recv;

  struct part_aos *parts_aos_pair_dens;
  struct part_aos_f4_send *parts_aos_pair_f4_send;
  struct part_aos_f4_recv *parts_aos_pair_f4_recv;

  struct part_aos *d_parts_aos_pair_dens;
  struct part_aos_f4_send *d_parts_aos_pair_f4_send;
  struct part_aos_f4_recv *d_parts_aos_pair_f4_recv;

  struct part_aos_f *parts_aos_pair_forc;
  struct part_aos_f4_f_send *parts_aos_pair_f4_f_send;
  struct part_aos_f4_f_recv *parts_aos_pair_f4_f_recv;

  struct part_aos_f *d_parts_aos_pair_forc;
  struct part_aos_f4_f_send *d_parts_aos_pair_f4_f_send;
  struct part_aos_f4_f_recv *d_parts_aos_pair_f4_f_recv;

  struct part_aos_g *parts_aos_pair_grad;
  struct part_aos_f4_g_send *parts_aos_pair_f4_g_send;
  struct part_aos_f4_g_recv *parts_aos_pair_f4_g_recv;

  struct part_aos_g *d_parts_aos_pair_grad;
  struct part_aos_f4_g_send *d_parts_aos_pair_f4_g_send;
  struct part_aos_f4_g_recv *d_parts_aos_pair_f4_g_recv;

  cudaMalloc((void **)&d_parts_aos_f4_send,
             count_max_parts_tmp * sizeof(struct part_aos_f4_send));
  cudaMalloc((void **)&d_parts_aos_f4_recv,
             count_max_parts_tmp * sizeof(struct part_aos_f4_recv));

  cudaMalloc((void **)&d_parts_aos_forc_f4_send,
             count_max_parts_tmp * sizeof(struct part_aos_f4_f_send));
  cudaMalloc((void **)&d_parts_aos_forc_f4_recv,
             count_max_parts_tmp * sizeof(struct part_aos_f4_f_recv));

  cudaMalloc((void **)&d_parts_aos_grad_f4_send,
             count_max_parts_tmp * sizeof(struct part_aos_f4_g_send));
  cudaMalloc((void **)&d_parts_aos_grad_f4_recv,
             count_max_parts_tmp * sizeof(struct part_aos_f4_g_recv));

  cudaMallocHost((void **)&parts_aos_f4_send,
                 count_max_parts_tmp * sizeof(struct part_aos_f4_send));
  cudaMallocHost((void **)&parts_aos_f4_recv,
                 count_max_parts_tmp * sizeof(struct part_aos_f4_recv));

  cudaMallocHost((void **)&parts_aos_forc_f4_send,
                 count_max_parts_tmp * sizeof(struct part_aos_f4_f_send));
  cudaMallocHost((void **)&parts_aos_forc_f4_recv,
                 count_max_parts_tmp * sizeof(struct part_aos_f4_f_recv));

  cudaMallocHost((void **)&parts_aos_grad_f4_send,
                 count_max_parts_tmp * sizeof(struct part_aos_f4_g_send));
  cudaMallocHost((void **)&parts_aos_grad_f4_recv,
                 count_max_parts_tmp * sizeof(struct part_aos_f4_g_recv));

  cudaMalloc((void **)&d_parts_aos_pair_f4_send,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f4_send));
  cudaMalloc((void **)&d_parts_aos_pair_f4_recv,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f4_recv));

  cudaMalloc((void **)&d_parts_aos_pair_f4_f_send,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f4_f_send));
  cudaMalloc((void **)&d_parts_aos_pair_f4_f_recv,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f4_f_recv));

  cudaMalloc((void **)&d_parts_aos_pair_f4_g_send,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f4_g_send));
  cudaMalloc((void **)&d_parts_aos_pair_f4_g_recv,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f4_g_recv));

  ///////////Probably not needed
  /// anymore////////////////////////////////////////////////////////////////
  cudaMalloc((void **)&d_parts_aos_pair_forc,
             2 * count_max_parts_tmp * sizeof(struct part_aos_f));
  cudaMalloc((void **)&d_parts_aos_pair_grad,
             2 * count_max_parts_tmp * sizeof(struct part_aos_g));
  ///////////Probably not needed
  /// anymore////////////////////////////////////////////////////////////////

  cudaMallocHost((void **)&parts_aos_pair_f4_send,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f4_send));
  cudaMallocHost((void **)&parts_aos_pair_f4_recv,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f4_recv));

  cudaMallocHost((void **)&parts_aos_pair_f4_g_send,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f4_g_send));
  cudaMallocHost((void **)&parts_aos_pair_f4_g_recv,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f4_g_recv));

  cudaMallocHost((void **)&parts_aos_pair_f4_f_send,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f4_f_send));
  cudaMallocHost((void **)&parts_aos_pair_f4_f_recv,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f4_f_recv));

  cudaMallocHost((void **)&parts_aos_pair_forc,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_f));
  cudaMallocHost((void **)&parts_aos_pair_grad,
                 2 * count_max_parts_tmp * sizeof(struct part_aos_g));

  /*Declare some global variables*/
  float d_a = e->cosmology->a;
  float d_H = e->cosmology->H;
  int step = 0;

  // a list of the cells and tasks the GPU will work on
  pack_vars_self_dens->task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_self_dens->cell_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));

  pack_vars_pair_dens->task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_pair_dens->top_task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_pair_dens->ci_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));
  pack_vars_pair_dens->cj_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));

  pack_vars_self_forc->task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_self_forc->cell_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));

  pack_vars_pair_forc->task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_pair_forc->ci_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));
  pack_vars_pair_forc->cj_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));

  pack_vars_self_grad->task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_self_grad->cell_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));

  pack_vars_pair_grad->task_list =
      (struct task **)calloc(target_n_tasks, sizeof(struct task *));
  pack_vars_pair_grad->ci_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));
  pack_vars_pair_grad->cj_list =
      (struct cell **)calloc(target_n_tasks, sizeof(struct cell *));

  // number of density self tasks executed
  int tasks_done_cpu = 0;
  int tasks_done_gpu = 0;
  int tasks_done_gpu_inc = 0;

  /* Main loop. */
  while (1) {
    /*Stuff for debugging*/
    int n_full_d_bundles = 0, n_full_g_bundles = 0, n_full_f_bundles = 0;
    int n_full_p_d_bundles = 0, n_full_p_g_bundles = 0, n_full_p_f_bundles = 0;
    int n_partial_d_bundles = 0, n_partial_g_bundles = 0,
        n_partial_f_bundles = 0;
    int n_partial_p_d_bundles = 0, n_partial_p_g_bundles = 0,
        n_partial_p_f_bundles = 0;
    int output = 0;
    int packed_self = 0;
    int packed_pair = 0;
    int packed_self_f = 0;
    int packed_pair_f = 0;
    int packed_self_g = 0;
    int packed_pair_g = 0;
    int density = 0;
    int density_sub = 0;
    int unpacked = 0;
    int unpacked_f = 0;
    int unpacked_g = 0;
    int unpacked_pair = 0;
    int unpacked_pair_f = 0;
    int unpacked_pair_g = 0;
    int ghost_in = 0;
    int cpu_self = 0;
    int cpu_self_f = 0;
    int cpu_self_g = 0;
    int cpu_pair = 0;
    int cpu_pair_f = 0;
    int cpu_pair_g = 0;
    int n_leafs_total = 0;
    //	Initialise timers to zero
    double time_for_density_cpu = 0.0;
    double time_for_density_cpu_pair = 0.0;
    double time_for_cpu_g = 0.0;
    double time_for_cpu_pair_g = 0.0;
    double time_for_cpu_f = 0.0;
    double time_for_cpu_pair_f = 0.0;
    double time_for_density_cpu_sub = 0.0;
    double time_for_density_gpu = 0.0;
    double time_for_density_gpu_pair = 0.0;
    double time_for_gpu_f = 0.0;
    double time_for_gpu_pair_f = 0.0;
    double time_for_gpu_g = 0.0;
    double time_for_gpu_pair_g = 0.0;
    double unpack_time_self_g = 0.0;
    double unpack_time_self_f = 0.0;
    double unpack_time_self = 0.0;
    double time_for_gpu_pair = 0.0;
    int nr_cells = space->nr_cells;
    /* Wait at the barrier. */
    engine_barrier(e);
    // Initialise packing counters
    pack_vars_self_dens->tasks_packed = 0;
    pack_vars_pair_dens->tasks_packed = 0;
    pack_vars_self_dens->count_parts = 0;
    pack_vars_pair_dens->count_parts = 0;
    pack_vars_pair_dens->task_locked = 0;
    pack_vars_pair_dens->top_tasks_packed = 0;
    // Initialise packing counters
    pack_vars_self_forc->tasks_packed = 0;
    pack_vars_pair_forc->tasks_packed = 0;
    pack_vars_self_forc->count_parts = 0;
    pack_vars_pair_forc->count_parts = 0;
    // Initialise packing counters
    pack_vars_self_grad->tasks_packed = 0;
    pack_vars_pair_grad->tasks_packed = 0;
    pack_vars_self_grad->count_parts = 0;
    pack_vars_pair_grad->count_parts = 0;

    int total_tasks_packed_this_time_pair = 0;
    double packing_time = 0.0;
    double packing_time_f = 0.0;
    double packing_time_g = 0.0;
    double unpacking_time = 0.0;
    double unpacking_time_f = 0.0;
    double unpacking_time_g = 0.0;
    double packing_time_pair = 0.0;
    double packing_time_pair_f = 0.0;
    double packing_time_pair_g = 0.0;
    double unpacking_time_pair = 0.0;
    double unpacking_time_pair_f = 0.0;
    double unpacking_time_pair_g = 0.0;
    double time_for_copy_to_struct = 0.0;
    double tot_time_for_hard_memcpys = 0.0;
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done) break;
    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    /*Some bits for output in case of debug*/
    char buf5[20];
    snprintf(buf5, sizeof(buf5), "t%dr%dstep%d", r->cpuid, engine_rank, step);
#ifdef DUMP_TIMINGS
    FILE *fgpu_steps;
    fgpu_steps = fopen(buf5, "w");
#endif
    //    if (step == 0) cudaProfilerStart();
    step++;

    sched->nr_packs_self_dens_done = 0;
    sched->nr_packs_pair_dens_done = 0;
    sched->nr_packs_self_forc_done = 0;
    sched->nr_packs_pair_forc_done = 0;
    sched->nr_packs_self_grad_done = 0;
    sched->nr_packs_pair_grad_done = 0;
    int n_cells_d = 0;
    int n_cells_g = 0;
    int n_cells_f = 0;
    int n_cells_p_d = 0;
    int n_cells_p_g = 0;
    int n_cells_p_f = 0;
    int n_w_prts_gtr_target_d = 0;
    int n_w_prts_gtr_target_g = 0;
    int n_w_prts_gtr_target_f = 0;
    int n_w_prts_gtr_target_p_d = 0;
    int n_w_prts_gtr_target_p_g = 0;
    int n_w_prts_gtr_target_p_f = 0;
    int g100 = 0;
    int l100 = 0;
    int maxcount = 0;
    /* Loop while there are tasks... */
    tasks_done_gpu_inc = 0;
    ticks hang_time = getticks();
    while (1) {
      // A. Nasar: Get qid for re-use later
      int qid = r->qid;
      /* If there's no old task, try to get a new one. */
      if (t == NULL) {
        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, qid, prev);
        TIMER_TOC(timer_gettask);
        /* Did I get anything? */
        if (t == NULL) break;
      }
      /* Get the cells. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      if (ci == NULL && (t->subtype != task_subtype_gpu_unpack_d
    		  && t->subtype != task_subtype_gpu_unpack_g
			  && t->subtype != task_subtype_gpu_unpack_f)) error("This cannot be");

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        struct cell *ci_temp = ci;
        struct cell *cj_temp = cj;
        double shift[3];
        if (t->subtype != task_subtype_gpu_unpack_d &&
            t->subtype != task_subtype_gpu_unpack_g &&
            t->subtype != task_subtype_gpu_unpack_f)
          t->sid = space_getsid_and_swap_cells(e->s, &ci_temp, &cj_temp, shift);
      } else {
        t->sid = -1;
      }
#endif

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      r->t = t;
#endif

      const ticks task_beg = getticks();
      /* Different types of tasks... */
      switch (t->type) {
        case task_type_self:
          if (t->subtype == task_subtype_gpu_unpack_d) {
            unpacked++;
          } else if (t->subtype == task_subtype_gpu_unpack_g) {
            unpacked_g++;
          } else if (t->subtype == task_subtype_gpu_unpack_f) {
            unpacked_f++;
          } else if (t->subtype == task_subtype_density) {
            cpu_self++;
#ifndef GPUOFFLOAD_DENSITY
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_doself1_branch_density(r, ci);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_density_cpu += (t1.tv_sec - t0.tv_sec) +
                                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
            density++;
#endif
            /* GPU WORK */
          } else if (t->subtype == task_subtype_gpu_pack_d) {
            packed_self++;
#ifdef GPUOFFLOAD_DENSITY
            ticks tic_cpu_pack = getticks();

            packing_time +=
                runner_doself1_pack_f4(r, sched, pack_vars_self_dens, ci, t,
                                       parts_aos_f4_send, task_first_part_f4);

            t->total_cpu_pack_ticks += getticks() - tic_cpu_pack;

            /* No pack tasks left in queue, flag that we want to run */
            int launch_leftovers = pack_vars_self_dens->launch_leftovers;
            n_cells_d++;
            maxcount = max(maxcount, ci->hydro.count);
            if (ci->hydro.count > 1.5 * np_per_cell) {
              n_w_prts_gtr_target_d++;
            }
            /*Packed enough tasks. Let's go*/
            int launch = pack_vars_self_dens->launch;
            /* Do we have enough stuff to run the GPU ? */
            if (launch) n_full_d_bundles++;
            if (launch_leftovers) n_partial_d_bundles++;
            if (launch || launch_leftovers) {
              /*Launch GPU tasks*/
              int t_packed = pack_vars_self_dens->tasks_packed;
              runner_doself1_launch_f4(
                  r, sched, pack_vars_self_dens, ci, t, parts_aos_f4_send,
                  parts_aos_f4_recv, d_parts_aos_f4_send, d_parts_aos_f4_recv,
                  stream, d_a, d_H, e, &packing_time, &time_for_density_gpu,
                  &unpack_time_self, task_first_part_self_dens_f4, devId,
                  task_first_part_f4, d_task_first_part_f4, self_end);
            } /*End of GPU work Self*/
#endif
          } /* self / pack */
          else if (t->subtype == task_subtype_gpu_pack_g) {
            packed_self_g++;
#ifdef GPUOFFLOAD_GRADIENT

            ticks tic_cpu_pack = getticks();

            n_cells_g++;
            maxcount = max(maxcount, ci->hydro.count);
            if (ci->hydro.count > 1.5 * np_per_cell) {
              n_w_prts_gtr_target_g++;
            }
            packing_time_g += runner_doself1_pack_f4_g(
                r, sched, pack_vars_self_grad, ci, t, parts_aos_grad_f4_send,
                task_first_part_f4_g);

            t->total_cpu_pack_ticks += getticks() - tic_cpu_pack;

            /* No pack tasks left in queue, flag that we want to run */
            int launch_leftovers = pack_vars_self_grad->launch_leftovers;
            /*Packed enough tasks let's go*/
            int launch = pack_vars_self_grad->launch;

            /* Do we have enough stuff to run the GPU ? */
            if (launch || launch_leftovers) {
              /*Launch GPU tasks*/
              int t_packed = pack_vars_self_grad->tasks_packed;
              //              signal_sleeping_runners(sched, t, t_packed);
              runner_doself1_launch_f4_g(
                  r, sched, pack_vars_self_grad, ci, t, parts_aos_grad_f4_send,
                  parts_aos_grad_f4_recv, d_parts_aos_grad_f4_send,
                  d_parts_aos_grad_f4_recv, stream, d_a, d_H, e,
                  &packing_time_g, &time_for_gpu_g, task_first_part_f4_g,
                  d_task_first_part_f4_g, self_end_g, &unpack_time_self_g);
            } /*End of GPU work Self*/
#endif  // GPUGRADSELF
          } else if (t->subtype == task_subtype_gpu_pack_f) {
            packed_self_f++;
#ifdef GPUOFFLOAD_FORCE
            ticks tic_cpu_pack = getticks();

            n_cells_f++;
            maxcount = max(maxcount, ci->hydro.count);
            if (ci->hydro.count > 1.5 * np_per_cell) {
              n_w_prts_gtr_target_f++;
            }
            packing_time_f += runner_doself1_pack_f4_f(
                r, sched, pack_vars_self_forc, ci, t, parts_aos_forc_f4_send,
                task_first_part_f4_f);

            t->total_cpu_pack_ticks += getticks() - tic_cpu_pack;

            /* No pack tasks left in queue, flag that we want to run */
            int launch_leftovers = pack_vars_self_forc->launch_leftovers;
            /*Packed enough tasks let's go*/
            int launch = pack_vars_self_forc->launch;

            /* Do we have enough stuff to run the GPU ? */
            if (launch || launch_leftovers) {
              /*Launch GPU tasks*/
              int t_packed = pack_vars_self_forc->tasks_packed;
              //              signal_sleeping_runners(sched, t, t_packed);
              runner_doself1_launch_f4_f(
                  r, sched, pack_vars_self_forc, ci, t, parts_aos_forc_f4_send,
                  parts_aos_forc_f4_recv, d_parts_aos_forc_f4_send,
                  d_parts_aos_forc_f4_recv, stream, d_a, d_H, e,
                  &packing_time_f, &time_for_gpu_f, task_first_part_f4_f,
                  d_task_first_part_f4_f, self_end_f, &unpack_time_self_f);
            } /*End of GPU work Self*/
#endif
          }
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient) {
            cpu_self_g++;
#ifndef GPUOFFLOAD_GRADIENT
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_doself1_branch_gradient(r, ci);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_cpu_g += (t1.tv_sec - t0.tv_sec) +
                              (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
#endif
          }
#endif
          else if (t->subtype == task_subtype_force) {
            cpu_self_f++;
#ifndef GPUOFFLOAD_FORCE
            struct timespec t0, t1;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_doself2_branch_force(r, ci);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_cpu_f += (t1.tv_sec - t0.tv_sec) +
                              (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
#endif
          } else if (t->subtype == task_subtype_limiter)
            runner_doself1_branch_limiter(r, ci);
          else if (t->subtype == task_subtype_grav)
            runner_doself_recursive_grav(r, ci, 1);
          else if (t->subtype == task_subtype_external_grav)
            runner_do_grav_external(r, ci, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_doself_branch_stars_density(r, ci);
#ifdef EXTRA_STAR_LOOPS
          else if (t->subtype == task_subtype_stars_prep1)
            runner_doself_branch_stars_prep1(r, ci);
          else if (t->subtype == task_subtype_stars_prep2)
            runner_doself_branch_stars_prep2(r, ci);
#endif
          else if (t->subtype == task_subtype_stars_feedback)
            runner_doself_branch_stars_feedback(r, ci);
          else if (t->subtype == task_subtype_bh_density)
            runner_doself_branch_bh_density(r, ci);
          else if (t->subtype == task_subtype_bh_swallow)
            runner_doself_branch_bh_swallow(r, ci);
          else if (t->subtype == task_subtype_do_gas_swallow)
            runner_do_gas_swallow_self(r, ci, 1);
          else if (t->subtype == task_subtype_do_bh_swallow)
            runner_do_bh_swallow_self(r, ci, 1);
          else if (t->subtype == task_subtype_bh_feedback)
            runner_doself_branch_bh_feedback(r, ci);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_doself1_branch_rt_gradient(r, ci);
          else if (t->subtype == task_subtype_rt_transport)
            runner_doself2_branch_rt_transport(r, ci);
          else if (t->subtype == task_subtype_sink_swallow)
            runner_doself_branch_sinks_swallow(r, ci);
          else if (t->subtype == task_subtype_sink_do_gas_swallow)
            runner_do_sinks_gas_swallow_self(r, ci, 1);
          else if (t->subtype == task_subtype_sink_do_sink_swallow)
            runner_do_sinks_sink_swallow_self(r, ci, 1);
          else
            error("Unknown/invalid task subtype (%s).",
                  subtaskID_names[t->subtype]);
          break;

        case task_type_pair:
          if (t->subtype == task_subtype_density) {
            /* Abouzied: To be commented out when the GPU pairs have been coded
             * up */
            cpu_pair++;
#ifndef GPUOFFLOAD_DENSITY
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_dopair1_branch_density(r, ci, cj);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_density_cpu_pair +=
                (t1.tv_sec - t0.tv_sec) +
                (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
#endif
          }
          /* GPU WORK */
          else if (t->subtype == task_subtype_gpu_pack_d) {
            packed_pair++;
#ifdef GPUOFFLOAD_DENSITY
#ifdef DO_CORNERS
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            double shift[3] = {0.0};
            t->corner_pair = 0;
            int sid = space_getsid_filter(e->s, &ci, &cj, shift);
            clock_gettime(CLOCK_REALTIME, &t1);
            packing_time_pair += (t1.tv_sec - t0.tv_sec) +
                                 (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
            if ((sid == 0 || sid == 2 || sid == 6 || sid == 8) && step > 1) {
              clock_gettime(CLOCK_REALTIME, &t0);
              runner_dopair1_branch_density(r, ci, cj);
              t->corner_pair = 1;
              int qid = r->qid;
              atomic_dec(&(sched->queues[qid].n_packs_pair_left));
              /* Tell the cells they have been packed */
              ci->pack_done++;
              cj->pack_done++;
              t->done = 1;
              int launch = 0, launch_leftovers = 0;
              if ((sched->queues[qid].n_packs_pair_left == 0))
                launch_leftovers = 1;
              /* Tasks done. Release the lock ! */
              task_unlock(t);
              /*schedule my dependencies (Only unpacks really)*/
              enqueue_dependencies(sched, t);
              /*Signal sleeping runners*/
              signal_sleeping_runners(sched, t);
              clock_gettime(CLOCK_REALTIME, &t1);
              packing_time_pair += (t1.tv_sec - t0.tv_sec) +
                                   (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
              if (launch_leftovers) {
                pack_vars_pair_dens->launch_leftovers = 1;
                runner_dopair1_launch_f4_one_memcpy(
                    r, sched, pack_vars_pair_dens, t, parts_aos_pair_f4_send,
                    parts_aos_pair_f4_recv, d_parts_aos_pair_f4_send,
                    d_parts_aos_pair_f4_recv, stream_pairs, d_a, d_H, e,
                    &packing_time_pair, &time_for_density_gpu_pair,
                    &unpacking_time_pair, fparti_fpartj_lparti_lpartj_dens,
                    pair_end);
              }
            } else {
#endif  // DO_CORNERS

              ticks tic_cpu_pack = getticks();
              n_cells_p_d++;
              maxcount = max(maxcount, ci->hydro.count);
              if (ci->hydro.count > 1.5 * np_per_cell) {
                n_w_prts_gtr_target_p_d++;
              }

              /*Call recursion here. This will be a function in runner_doiact_functions_hydro_gpu.h.
               * We are recursing separately to find out how much work we have before offloading*/
              //We need to allocate a list to put cell pointers into for each new task
              int n_expected_tasks = 1024;
              int n_leafs_found = 0;
              int depth = 0;
              struct cell * cells_left[n_expected_tasks];
              struct cell * cells_right[n_expected_tasks];
              runner_recurse_gpu(r, sched, pack_vars_pair_dens, ci, cj, t,
                      parts_aos_pair_f4_send, e, fparti_fpartj_lparti_lpartj_dens, &n_leafs_found, cells_left, cells_right, depth);
              n_leafs_total += n_leafs_found;

              int cstart = 0, cend = n_leafs_found;

              int cid = 0;
              pack_vars_pair_dens->task_locked = 1;
              int top_tasks_packed = pack_vars_pair_dens->top_tasks_packed;
              pack_vars_pair_dens->top_tasks_packed++;
              pack_vars_pair_dens->top_task_list[top_tasks_packed] = t;
              int t_s, t_e;
              t_s = 0;
              while(cid < n_leafs_found){
                //////////////////////////////////////////////////////////////////////////////////
                /*Loop through n_daughters such that the pack_vars_pair_dens counters are updated*/
                for (cid = cstart; pack_vars_pair_dens->tasks_packed < pack_vars_pair_dens->target_n_tasks
                     && cid < n_leafs_found; cid++){

                  packing_time_pair += runner_dopair1_pack_f4(
                  r, sched, pack_vars_pair_dens, cells_left[cid], cells_right[cid], t,
                  parts_aos_pair_f4_send, e, fparti_fpartj_lparti_lpartj_dens);
//                  if (pack_vars_pair_dens->unfinished)
//                	break;
//                message("Packing task %i in recursed tasks\n", cid);
                }
                /* Copies done. Release the lock ! */
                pack_vars_pair_dens->task_locked = 0;
//                if(cid == n_leafs_found){
//                  cell_unlocktree(ci);
//                  cell_unlocktree(cj);
//                  pack_vars_pair_dens->task_locked = 0;
//                }
                cstart = cid + 1;
                t->total_cpu_pack_ticks += getticks() - tic_cpu_pack;
                /* Packed enough tasks or no pack tasks left in queue, flag that
                 * we want to run */
                int launch = pack_vars_pair_dens->launch;
                int launch_leftovers = pack_vars_pair_dens->launch_leftovers;

                /* Do we have enough stuff to run the GPU ? */
                if (launch) n_full_p_d_bundles++;
                if (launch_leftovers) n_partial_p_d_bundles++;

                if (launch || launch_leftovers) {
                  /*Launch GPU tasks*/
                  int t_packed = pack_vars_pair_dens->tasks_packed;
                  //                signal_sleeping_runners(sched, t, t_packed);
                  runner_dopair1_launch_f4_one_memcpy(
                    r, sched, pack_vars_pair_dens, t, parts_aos_pair_f4_send,
                    parts_aos_pair_f4_recv, d_parts_aos_pair_f4_send,
                    d_parts_aos_pair_f4_recv, stream_pairs, d_a, d_H, e,
                    &packing_time_pair, &time_for_density_gpu_pair,
                    &unpacking_time_pair, fparti_fpartj_lparti_lpartj_dens,
                    pair_end);
                  for (int tid = 0; tid < pack_vars_pair_dens->top_tasks_packed -1; tid++){
                    /*schedule my dependencies (Only unpacks really)*/
                	struct task *tii = pack_vars_pair_dens->top_task_list[tid];
                    enqueue_dependencies(sched, tii);
                  }
                  pack_vars_pair_dens->top_tasks_packed = 1;
                  pack_vars_pair_dens->top_task_list[0] = t;
                }
                ///////////////////////////////////////////////////////////////////////
              }
              cell_unlocktree(ci);
              cell_unlocktree(cj);
              pack_vars_pair_dens->task_locked = 0;
              pack_vars_pair_dens->launch_leftovers = 0;

#ifdef DO_CORNERS
            } /* End of GPU work Pairs */
#endif  // DO_CORNERS
#endif  // GPUOFFLOAD_DENSITY
          } /* pair / pack */
          else if (t->subtype == task_subtype_gpu_pack_g) {
            packed_pair_g++;
#ifdef GPUOFFLOAD_GRADIENT
#ifdef DO_CORNERS
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            double shift[3] = {0.0};
            t->corner_pair = 0;
            int sid = space_getsid_filter(e->s, &ci, &cj, shift);
            clock_gettime(CLOCK_REALTIME, &t1);
            packing_time_pair += (t1.tv_sec - t0.tv_sec) +
                                 (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
            if ((sid == 0 || sid == 2 || sid == 6 || sid == 8) && step > 1) {
              clock_gettime(CLOCK_REALTIME, &t0);
              runner_dopair1_branch_gradient(r, ci, cj);
              t->corner_pair = 1;
              int qid = r->qid;
              atomic_dec(&(sched->queues[qid].n_packs_pair_left_g));
              /* Tell the cells they have been packed */
              ci->pack_done++;
              cj->pack_done++;
              t->done = 1;
              int launch = 0, launch_leftovers = 0;
              if ((sched->queues[qid].n_packs_pair_left_g == 0))
                launch_leftovers = 1;
              /* Tasks done. Release the lock ! */
              task_unlock(t);
              /*schedule my dependencies (Only unpacks really)*/
              enqueue_dependencies(sched, t);
              /*Signal sleeping runners*/
              signal_sleeping_runners(sched, t);
              clock_gettime(CLOCK_REALTIME, &t1);
              packing_time_pair_g += (t1.tv_sec - t0.tv_sec) +
                                     (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
              if (launch_leftovers) {
                pack_vars_pair_grad->launch_leftovers = 1;
                runner_dopair1_launch_f4_g_one_memcpy(
                    r, sched, pack_vars_pair_grad, t, parts_aos_pair_f4_g_send,
                    parts_aos_pair_f4_g_recv, d_parts_aos_pair_f4_g_send,
                    d_parts_aos_pair_f4_g_recv, stream_pairs, d_a, d_H, e,
                    &packing_time_pair_g, &time_for_gpu_pair_g,
                    &unpacking_time_pair_g, fparti_fpartj_lparti_lpartj_grad,
                    pair_end_g);
              }
            } else {
#endif  // DO_CORNERS
              ticks tic_cpu_pack = getticks();
              n_cells_p_g++;
              maxcount = max(maxcount, ci->hydro.count);
              if (ci->hydro.count > 1.5 * np_per_cell) {
                n_w_prts_gtr_target_p_g++;
  //              message("count %i target %i", ci->hydro.count, np_per_cell);
              }
              packing_time_pair_g +=
                  runner_dopair1_pack_f4_g(r, sched, pack_vars_pair_grad, ci,
                                           cj, t, parts_aos_pair_f4_g_send, e,
                                           fparti_fpartj_lparti_lpartj_grad);

              t->total_cpu_pack_ticks += getticks() - tic_cpu_pack;

              /* No pack tasks left in queue, flag that we want to run */
              int launch_leftovers = pack_vars_pair_grad->launch_leftovers;
              /*Packed enough tasks, let's go*/
              int launch = pack_vars_pair_grad->launch;

              /* Do we have enough stuff to run the GPU ? */
              if (launch || launch_leftovers) {
                /*Launch GPU tasks*/
                int t_packed = pack_vars_pair_grad->tasks_packed;
                //                signal_sleeping_runners(sched, t, t_packed);
                runner_dopair1_launch_f4_g_one_memcpy(
                    r, sched, pack_vars_pair_grad, t, parts_aos_pair_f4_g_send,
                    parts_aos_pair_f4_g_recv, d_parts_aos_pair_f4_g_send,
                    d_parts_aos_pair_f4_g_recv, stream_pairs, d_a, d_H, e,
                    &packing_time_pair_g, &time_for_gpu_pair_g,
                    &unpacking_time_pair_g, fparti_fpartj_lparti_lpartj_grad,
                    pair_end_g);
              }
              pack_vars_pair_grad->launch_leftovers = 0;
#ifdef DO_CORNERS
            } /* End of GPU work Pairs */
#endif  // DO_CORNERS
#endif  // GPUOFFLOAD_GRADIENT
          } else if (t->subtype == task_subtype_gpu_pack_f) {
            packed_pair_f++;
#ifdef GPUOFFLOAD_FORCE
#ifdef DO_CORNERS
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            double shift[3] = {0.0};
            t->corner_pair = 0;
            int sid = space_getsid_filter(e->s, &ci, &cj, shift);
            clock_gettime(CLOCK_REALTIME, &t1);
            packing_time_pair += (t1.tv_sec - t0.tv_sec) +
                                 (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
            if ((sid == 0 || sid == 2 || sid == 6 || sid == 8) && step > 1) {
              //          if((sid != 4 && sid != 10 && sid == 12) && step > 1){
              runner_dopair1_branch_force(r, ci, cj);
              t->corner_pair = 1;
              int qid = r->qid;
              atomic_dec(&(sched->queues[qid].n_packs_pair_left_f));
              /* Tell the cells they have been packed */
              ci->pack_done++;
              cj->pack_done++;
              t->done = 1;
              int launch = 0, launch_leftovers = 0;
              if ((sched->queues[qid].n_packs_pair_left_f == 0))
                launch_leftovers = 1;
              /* Tasks done. Release the lock ! */
              task_unlock(t);
              /*schedule my dependencies (Only unpacks really)*/
              enqueue_dependencies(sched, t);
              /*Signal sleeping runners*/
              signal_sleeping_runners(sched, t);
              clock_gettime(CLOCK_REALTIME, &t1);
              packing_time_pair_f += (t1.tv_sec - t0.tv_sec) +
                                     (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
              if (launch_leftovers) {
                pack_vars_pair_forc->launch_leftovers = 1;
                runner_dopair1_launch_f4_f_one_memcpy(
                    r, sched, pack_vars_pair_forc, t, parts_aos_pair_f4_f_send,
                    parts_aos_pair_f4_f_recv, d_parts_aos_pair_f4_f_send,
                    d_parts_aos_pair_f4_f_recv, stream_pairs, d_a, d_H, e,
                    &packing_time_pair_f, &time_for_gpu_pair_f,
                    &unpacking_time_pair_f, fparti_fpartj_lparti_lpartj_forc,
                    pair_end_f);
              }
            } else {
#endif  // DO_CORNERS

              ticks tic_cpu_pack = getticks();

              packing_time_pair_f +=
                  runner_dopair1_pack_f4_f(r, sched, pack_vars_pair_forc, ci,
                                           cj, t, parts_aos_pair_f4_f_send, e,
                                           fparti_fpartj_lparti_lpartj_forc);
              n_cells_p_f++;
              maxcount = max(maxcount, ci->hydro.count);
              if (ci->hydro.count > 1.5 * np_per_cell) {
                n_w_prts_gtr_target_p_f++;
  //              message("count %i target %i", ci->hydro.count, np_per_cell);
              }
              t->total_cpu_pack_ticks += getticks() - tic_cpu_pack;

              /* No pack tasks left in queue, flag that we want to run */
              int launch_leftovers = pack_vars_pair_forc->launch_leftovers;
              /*Packed enough tasks let's go*/
              int launch = pack_vars_pair_forc->launch;
              //              if ((sched->p_f_left[qid] < 1)){
              //            	  launch_leftovers = 1;
              //            	  pack_vars_pair_forc->launch_leftovers = 1;
              //              }
              /* Do we have enough stuff to run the GPU ? */
              if (launch || launch_leftovers) {
                /*Launch GPU tasks*/
                int t_packed = pack_vars_pair_forc->tasks_packed;
                //                signal_sleeping_runners(sched, t, t_packed);
                runner_dopair1_launch_f4_f_one_memcpy(
                    r, sched, pack_vars_pair_forc, t, parts_aos_pair_f4_f_send,
                    parts_aos_pair_f4_f_recv, d_parts_aos_pair_f4_f_send,
                    d_parts_aos_pair_f4_f_recv, stream_pairs, d_a, d_H, e,
                    &packing_time_pair_f, &time_for_gpu_pair_f,
                    &unpacking_time_pair_f, fparti_fpartj_lparti_lpartj_forc,
                    pair_end_f);

                pack_vars_pair_forc->launch_leftovers = 0;
              } /* End of GPU work Pairs */
#ifdef DO_CORNERS
            }
#endif  // DO_CORNERS
#endif  // GPUOFFLOAD_FORCE
          } else if (t->subtype == task_subtype_gpu_unpack_d) {
            unpacked_pair++;
          } else if (t->subtype == task_subtype_gpu_unpack_g) {
            unpacked_pair_g++;
          } else if (t->subtype == task_subtype_gpu_unpack_f) {
            unpacked_pair_f++;
          }
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient) {
            int Do_nothing = 0;
#ifndef GPUOFFLOAD_GRADIENT
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_dopair1_branch_gradient(r, ci, cj);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_cpu_pair_g += (t1.tv_sec - t0.tv_sec) +
                                   (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
#endif
          }
#endif  // EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_force) {
            int Do_nothing = 0;
#ifndef GPUOFFLOAD_FORCE
            struct timespec t0, t1, dt;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_dopair2_branch_force(r, ci, cj);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_cpu_pair_f += (t1.tv_sec - t0.tv_sec) +
                                   (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
#endif  // GPUOFFLOAD_FORCE
          } else if (t->subtype == task_subtype_limiter)
            runner_dopair1_branch_limiter(r, ci, cj);
          else if (t->subtype == task_subtype_grav)
            runner_dopair_recursive_grav(r, ci, cj, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_dopair_branch_stars_density(r, ci, cj);
#ifdef EXTRA_STAR_LOOPS
          else if (t->subtype == task_subtype_stars_prep1)
            runner_dopair_branch_stars_prep1(r, ci, cj);
          else if (t->subtype == task_subtype_stars_prep2)
            runner_dopair_branch_stars_prep2(r, ci, cj);
#endif
          else if (t->subtype == task_subtype_stars_feedback)
            runner_dopair_branch_stars_feedback(r, ci, cj);
          else if (t->subtype == task_subtype_bh_density)
            runner_dopair_branch_bh_density(r, ci, cj);
          else if (t->subtype == task_subtype_bh_swallow)
            runner_dopair_branch_bh_swallow(r, ci, cj);
          else if (t->subtype == task_subtype_do_gas_swallow)
            runner_do_gas_swallow_pair(r, ci, cj, 1);
          else if (t->subtype == task_subtype_do_bh_swallow)
            runner_do_bh_swallow_pair(r, ci, cj, 1);
          else if (t->subtype == task_subtype_bh_feedback)
            runner_dopair_branch_bh_feedback(r, ci, cj);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_dopair1_branch_rt_gradient(r, ci, cj);
          else if (t->subtype == task_subtype_rt_transport)
            runner_dopair2_branch_rt_transport(r, ci, cj);
          else if (t->subtype == task_subtype_sink_swallow)
            runner_dopair_branch_sinks_swallow(r, ci, cj);
          else if (t->subtype == task_subtype_sink_do_gas_swallow)
            runner_do_sinks_gas_swallow_pair(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_do_sink_swallow)
            runner_do_sinks_sink_swallow_pair(r, ci, cj, 1);
          else
            error("Unknown/invalid task subtype (%s/%s).",
                  taskID_names[t->type], subtaskID_names[t->subtype]);
          break;

        case task_type_sub_self:
          if (t->subtype == task_subtype_density) {
            struct timespec t0, t1, dt;
            const int count = ci->hydro.count;
            density_sub++;
            clock_gettime(CLOCK_REALTIME, &t0);
            runner_dosub_self1_density(r, ci, 1);
            clock_gettime(CLOCK_REALTIME, &t1);
            tasks_done_cpu++;
            time_for_density_cpu_sub +=
                (t1.tv_sec - t0.tv_sec) +
                (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
          }
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient) {
            runner_dosub_self1_gradient(r, ci, 1);
          }
#endif
          else if (t->subtype == task_subtype_force) {
            runner_dosub_self2_force(r, ci, 1);
          } else if (t->subtype == task_subtype_limiter)
            runner_dosub_self1_limiter(r, ci, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_dosub_self_stars_density(r, ci, 1);
#ifdef EXTRA_STAR_LOOPS
          else if (t->subtype == task_subtype_stars_prep1)
            runner_dosub_self_stars_prep1(r, ci, 1);
          else if (t->subtype == task_subtype_stars_prep2)
            runner_dosub_self_stars_prep2(r, ci, 1);
#endif
          else if (t->subtype == task_subtype_stars_feedback)
            runner_dosub_self_stars_feedback(r, ci, 1);
          else if (t->subtype == task_subtype_bh_density)
            runner_dosub_self_bh_density(r, ci, 1);
          else if (t->subtype == task_subtype_bh_swallow)
            runner_dosub_self_bh_swallow(r, ci, 1);
          else if (t->subtype == task_subtype_do_gas_swallow)
            runner_do_gas_swallow_self(r, ci, 1);
          else if (t->subtype == task_subtype_do_bh_swallow)
            runner_do_bh_swallow_self(r, ci, 1);
          else if (t->subtype == task_subtype_bh_feedback)
            runner_dosub_self_bh_feedback(r, ci, 1);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_dosub_self1_rt_gradient(r, ci, 1);
          else if (t->subtype == task_subtype_rt_transport)
            runner_dosub_self2_rt_transport(r, ci, 1);
          else if (t->subtype == task_subtype_sink_swallow)
            runner_dosub_self_sinks_swallow(r, ci, 1);
          else if (t->subtype == task_subtype_sink_do_gas_swallow)
            runner_do_sinks_gas_swallow_self(r, ci, 1);
          else if (t->subtype == task_subtype_sink_do_sink_swallow)
            runner_do_sinks_sink_swallow_self(r, ci, 1);
          else
            error("Unknown/invalid task subtype (%s/%s).",
                  taskID_names[t->type], subtaskID_names[t->subtype]);
          break;

        case task_type_sub_pair:
          if (t->subtype == task_subtype_density) {
            int nothing = 0;
            runner_dosub_pair1_density(r, ci, cj, 1);
          }
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient) {
            runner_dosub_pair1_gradient(r, ci, cj, 1);
          }
#endif
          else if (t->subtype == task_subtype_force) {
            runner_dosub_pair2_force(r, ci, cj, 1);
          } else if (t->subtype == task_subtype_limiter)
            runner_dosub_pair1_limiter(r, ci, cj, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_dosub_pair_stars_density(r, ci, cj, 1);
#ifdef EXTRA_STAR_LOOPS
          else if (t->subtype == task_subtype_stars_prep1)
            runner_dosub_pair_stars_prep1(r, ci, cj, 1);
          else if (t->subtype == task_subtype_stars_prep2)
            runner_dosub_pair_stars_prep2(r, ci, cj, 1);
#endif
          else if (t->subtype == task_subtype_stars_feedback)
            runner_dosub_pair_stars_feedback(r, ci, cj, 1);
          else if (t->subtype == task_subtype_bh_density)
            runner_dosub_pair_bh_density(r, ci, cj, 1);
          else if (t->subtype == task_subtype_bh_swallow)
            runner_dosub_pair_bh_swallow(r, ci, cj, 1);
          else if (t->subtype == task_subtype_do_gas_swallow)
            runner_do_gas_swallow_pair(r, ci, cj, 1);
          else if (t->subtype == task_subtype_do_bh_swallow)
            runner_do_bh_swallow_pair(r, ci, cj, 1);
          else if (t->subtype == task_subtype_bh_feedback)
            runner_dosub_pair_bh_feedback(r, ci, cj, 1);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_dosub_pair1_rt_gradient(r, ci, cj, 1);
          else if (t->subtype == task_subtype_rt_transport)
            runner_dosub_pair2_rt_transport(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_swallow)
            runner_dosub_pair_sinks_swallow(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_do_gas_swallow)
            runner_do_sinks_gas_swallow_pair(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_do_sink_swallow)
            runner_do_sinks_sink_swallow_pair(r, ci, cj, 1);
          else
            error("Unknown/invalid task subtype (%s/%s).",
                  taskID_names[t->type], subtaskID_names[t->subtype]);
          break;

        case task_type_sort:
          /* Cleanup only if any of the indices went stale. */
          runner_do_hydro_sort(
              r, ci, t->flags,
              ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin,
              cell_get_flag(ci, cell_flag_rt_requests_sort), 1);
          /* Reset the sort flags as our work here is done. */
          t->flags = 0;
          break;
        case task_type_rt_sort:
          /* Cleanup only if any of the indices went stale.
           * NOTE: we check whether we reset the sort flags when the
           * recv tasks are running. Cells without an RT recv task
           * don't have rt_sort tasks. */
          runner_do_hydro_sort(
              r, ci, t->flags,
              ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin, 1, 1);
          /* Reset the sort flags as our work here is done. */
          t->flags = 0;
          break;
        case task_type_stars_sort:
          /* Cleanup only if any of the indices went stale. */
          runner_do_stars_sort(
              r, ci, t->flags,
              ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin, 1);
          /* Reset the sort flags as our work here is done. */
          t->flags = 0;
          break;
        case task_type_init_grav:
          runner_do_init_grav(r, ci, 1);
          break;
        case task_type_ghost:
          runner_do_ghost(r, ci, 1);
          break;
#ifdef EXTRA_HYDRO_LOOP
        case task_type_extra_ghost:
          runner_do_extra_ghost(r, ci, 1);
          break;
#endif
        case task_type_stars_ghost:
          runner_do_stars_ghost(r, ci, 1);
          break;
        case task_type_bh_density_ghost:
          runner_do_black_holes_density_ghost(r, ci, 1);
          break;
        case task_type_bh_swallow_ghost3:
          runner_do_black_holes_swallow_ghost(r, ci, 1);
          break;
        case task_type_drift_part:
          runner_do_drift_part(r, ci, 1);
          break;
        case task_type_drift_spart:
          runner_do_drift_spart(r, ci, 1);
          break;
        case task_type_drift_sink:
          runner_do_drift_sink(r, ci, 1);
          break;
        case task_type_drift_bpart:
          runner_do_drift_bpart(r, ci, 1);
          break;
        case task_type_drift_gpart:
          runner_do_drift_gpart(r, ci, 1);
          break;
        case task_type_kick1:
          runner_do_kick1(r, ci, 1);
          break;
        case task_type_kick2:
          runner_do_kick2(r, ci, 1);
          break;
        case task_type_end_hydro_force:
          runner_do_end_hydro_force(r, ci, 1);
          break;
        case task_type_end_grav_force:
          runner_do_end_grav_force(r, ci, 1);
          break;
        case task_type_csds:
          runner_do_csds(r, ci, 1);
          break;
        case task_type_timestep:
          runner_do_timestep(r, ci, 1);
          break;
        case task_type_timestep_limiter:
          runner_do_limiter(r, ci, 0, 1);
          break;
        case task_type_timestep_sync:
          runner_do_sync(r, ci, 0, 1);
          break;
        case task_type_collect:
          runner_do_timestep_collect(r, ci, 1);
          break;
        case task_type_rt_collect_times:
          runner_do_collect_rt_times(r, ci, 1);
          break;
#ifdef WITH_MPI
        case task_type_send:
          if (t->subtype == task_subtype_tend) {
            free(t->buff);
          } else if (t->subtype == task_subtype_sf_counts) {
            free(t->buff);
          } else if (t->subtype == task_subtype_part_swallow) {
            free(t->buff);
          } else if (t->subtype == task_subtype_bpart_merger) {
            free(t->buff);
          } else if (t->subtype == task_subtype_limiter) {
            free(t->buff);
          }
          break;
        case task_type_recv:
          if (t->subtype == task_subtype_tend) {
            cell_unpack_end_step(ci, (struct pcell_step *)t->buff);
            free(t->buff);
          } else if (t->subtype == task_subtype_sf_counts) {
            cell_unpack_sf_counts(ci, (struct pcell_sf *)t->buff);
            cell_clear_stars_sort_flags(ci, /*clear_unused_flags=*/0);
            free(t->buff);
          } else if (t->subtype == task_subtype_xv) {
            runner_do_recv_part(r, ci, 1, 1);
          } else if (t->subtype == task_subtype_rho) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_gradient) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_rt_gradient) {
            runner_do_recv_part(r, ci, 2, 1);
          } else if (t->subtype == task_subtype_rt_transport) {
            runner_do_recv_part(r, ci, -1, 1);
          } else if (t->subtype == task_subtype_part_swallow) {
            cell_unpack_part_swallow(ci,
                                     (struct black_holes_part_data *)t->buff);
            free(t->buff);
          } else if (t->subtype == task_subtype_bpart_merger) {
            cell_unpack_bpart_swallow(ci,
                                      (struct black_holes_bpart_data *)t->buff);
            free(t->buff);
          } else if (t->subtype == task_subtype_limiter) {
            /* Nothing to do here. Unpacking done in a separate task */
          } else if (t->subtype == task_subtype_gpart) {
            runner_do_recv_gpart(r, ci, 1);
          } else if (t->subtype == task_subtype_spart_density) {
            runner_do_recv_spart(r, ci, 1, 1);
          } else if (t->subtype == task_subtype_part_prep1) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_spart_prep2) {
            runner_do_recv_spart(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_bpart_rho) {
            runner_do_recv_bpart(r, ci, 1, 1);
          } else if (t->subtype == task_subtype_bpart_feedback) {
            runner_do_recv_bpart(r, ci, 0, 1);
          } else {
            error("Unknown/invalid task subtype (%d).", t->subtype);
          }
          break;

        case task_type_pack:
          runner_do_pack_limiter(r, ci, &t->buff, 1);
          task_get_unique_dependent(t)->buff = t->buff;
          break;
        case task_type_unpack:
          runner_do_unpack_limiter(r, ci, t->buff, 1);
          break;
#endif
        case task_type_grav_down:
          runner_do_grav_down(r, t->ci, 1);
          break;
        case task_type_grav_long_range:
          runner_do_grav_long_range(r, t->ci, 1);
          break;
        case task_type_grav_mm:
          runner_dopair_grav_mm_progenies(r, t->flags, t->ci, t->cj);
          break;
        case task_type_cooling:
          runner_do_cooling(r, t->ci, 1);
          break;
        case task_type_star_formation:
          runner_do_star_formation(r, t->ci, 1);
          break;
        case task_type_star_formation_sink:
          runner_do_star_formation_sink(r, t->ci, 1);
          break;
        case task_type_stars_resort:
          runner_do_stars_resort(r, t->ci, 1);
          break;
        case task_type_sink_formation:
          runner_do_sink_formation(r, t->ci);
          break;
        case task_type_fof_self:
          runner_do_fof_search_self(r, t->ci, 1);
          break;
        case task_type_fof_pair:
          runner_do_fof_search_pair(r, t->ci, t->cj, 1);
          break;
        case task_type_fof_attach_self:
          runner_do_fof_attach_self(r, t->ci, 1);
          break;
        case task_type_fof_attach_pair:
          runner_do_fof_attach_pair(r, t->ci, t->cj, 1);
          break;
        case task_type_neutrino_weight:
          runner_do_neutrino_weighting(r, ci, 1);
          break;
        case task_type_rt_ghost1:
          runner_do_rt_ghost1(r, t->ci, 1);
          break;
        case task_type_rt_ghost2:
          runner_do_rt_ghost2(r, t->ci, 1);
          break;
        case task_type_rt_tchem:
          runner_do_rt_tchem(r, t->ci, 1);
          break;
        case task_type_rt_advance_cell_time:
          runner_do_rt_advance_cell_time(r, t->ci, 1);
          break;
        default:
          error("Unknown/invalid task type (%d).", t->type);
      }
      r->active_time += (getticks() - task_beg);

/* Mark that we have run this task on these cells */
#ifdef SWIFT_DEBUG_CHECKS
      if (ci != NULL) {
        ci->tasks_executed[t->type]++;
        ci->subtasks_executed[t->subtype]++;
      }
      if (cj != NULL) {
        cj->tasks_executed[t->type]++;
        cj->subtasks_executed[t->subtype]++;
      }
      /* This runner is not doing a task anymore */
      r->t = NULL;
#endif

      /* We're done with this task, see if we get a next one. */
      prev = t;

      if (t->subtype == task_subtype_gpu_pack_d) {
#ifdef GPUOFFLOAD_DENSITY
        /* Don't enqueue unpacks yet. Just signal the runners */
        t->skip = 1;
        t->toc = getticks();
        t->total_ticks += t->toc - t->tic;
        t = NULL;
#else
        t = scheduler_done(sched, t);
#endif
      }

      else if (t->subtype == task_subtype_gpu_pack_g) {
#ifdef GPUOFFLOAD_GRADIENT
        /* Don't enqueue unpacks yet. Just signal the runners */
        t->skip = 1;
        t->toc = getticks();
        t->total_ticks += t->toc - t->tic;
        t = NULL;
#else
        t = scheduler_done(sched, t);
#endif
      }

      else if (t->subtype == task_subtype_gpu_pack_f) {
#ifdef GPUOFFLOAD_FORCE
        /* Don't enqueue unpacks yet. Just signal the runners */
        t->skip = 1;
        t->toc = getticks();
        t->total_ticks += t->toc - t->tic;
        t = NULL;
#else
        t = scheduler_done(sched, t);
#endif
      }

      else if (t->subtype != task_subtype_gpu_pack_d &&
               t->subtype != task_subtype_gpu_pack_g &&
               t->subtype != task_subtype_gpu_pack_f) {
        t = scheduler_done(sched, t);
      }
    } /* main loop. */

    message("n_leafs found %i", n_leafs_total);
//    message("cpu %i packed %i cells with %i containing more parts than target of %i max_count %i",
//            r->cpuid, n_cells_d, n_w_prts_gtr_target_d, np_per_cell, maxcount);
//    message("cpu %i packed %i cells_G with %i containing more parts than target of %i max_count %i",
//            r->cpuid, n_cells_g, n_w_prts_gtr_target_g, np_per_cell, maxcount);
//    message("cpu %i packed %i cells_F with %i containing more parts than target of %i max_count %i",
//            r->cpuid, n_cells_f, n_w_prts_gtr_target_f, np_per_cell, maxcount);
//    message("cpu %i packed %i pairs_D with %i containing more parts than target of %i max_count %i",
//            r->cpuid, n_cells_p_d, n_w_prts_gtr_target_p_d, np_per_cell, maxcount);
//    message("cpu %i packed %i pairs_G with %i containing more parts than target of %i max_count %i",
//            r->cpuid, n_cells_p_g, n_w_prts_gtr_target_p_g, np_per_cell, maxcount);
//    message("cpu %i packed %i pairs_F with %i containing more parts than target of %i max_count %i",
//            r->cpuid, n_cells_p_f, n_w_prts_gtr_target_p_f, np_per_cell, maxcount);

    //    message("Worked on %i supers w more than 100 parts", g100);
    // Stuff for writing debug data to file for validation
    ////        if (step % 10 == 0 || step == 1) {
    //      if(r->cpuid == 0 && engine_rank == 0)fprintf(fgpu_steps, "x, y, z,
    //      rho, rhodh, v_sig, lap_u, a_visc_max, ax, ay, az\n"); for (int tid
    //      = 0; tid < space->nr_local_cells;
    //           tid++) { /* This should indeed be tasks_done_gpu as they are
    //           the only
    ////                     tasks which have been done*/
    //        struct cell *ctemp = &(space->cells_top[tid]);
    //        for (int i = 0; i < ctemp->hydro.count; i++) {
    //          fprintf(fgpu_steps, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,
    //          %f, %f\n",
    //                  ctemp->hydro.parts[i].x[0],
    //                  ctemp->hydro.parts[i].x[1],
    //                  ctemp->hydro.parts[i].x[2], ctemp->hydro.parts[i].rho,
    //                  ctemp->hydro.parts[i].density.rho_dh,
    //                  ctemp->hydro.parts[i].viscosity.v_sig,
    //                  ctemp->hydro.parts[i].diffusion.laplace_u,
    //                  ctemp->hydro.parts[i].force.alpha_visc_max_ngb,
    //                  ctemp->hydro.parts[i].a_hydro[0],
    //				  ctemp->hydro.parts[i].a_hydro[1],
    //				  ctemp->hydro.parts[i].a_hydro[2]);
    ////          message("wcount %f density %f",
    /// ctemp->hydro.parts[i].density.wcount, ctemp->hydro.parts[i].rho); /
    /// message("wcount is %f\n", ctemp->hydro.parts[i].density.wcount);
    //        }
    //      }
    ////  }
    /*Output compute times to separate files. cat later into one file*/
//    if (step % 11 == 0 || step == 1) {
#ifdef DUMP_TIMINGS
#if defined(GPUOFFLOAD_DENSITY) || defined(GPUOFFLOAD_GRADIENT) || \
    defined(GPUOFFLOAD_FORCE)
    //        char buffer[30];
    //        snprintf(buffer, sizeof(buffer), "t%d_stepnfullbundles%d",
    //        r->cpuid, step); FILE *fullbundles = fopen(buffer, "w");
    //        if(r->cpuid == 0)fprintf(fullbundles, "nfull, npartial,
    //        nfullpair, npartialpair\n"); else fprintf(fullbundles, "%i, %i,
    //        %i, %i\n", 		n_full_d_bundles, n_partial_d_bundles,
    //        n_full_p_d_bundles, n_partial_p_d_bundles); fflush(fullbundles);

    ///////////////////////////////////////////////////////////////
    /// to ooutput timings uncomment this
    ///////////////////////////////////////////////////////////////
    if (r->cpuid == 0 && engine_rank == 0)
      fprintf(fgpu_steps,
              "GPU_SD, P_SD, U_SD, GPU_PD,  P_PD, U_PD, "
              "GPU_SF, P_SF, U_SF, GPU_PF, P_PF, U_PF, GPU_SG, P_SG, U_SG, "
              "GPU_PG, P_PG, U_PG\n "
              "%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, "
              "%e, %e\n",
              time_for_density_gpu, packing_time, unpack_time_self,
              time_for_density_gpu_pair, packing_time_pair, unpacking_time_pair,
              time_for_gpu_f, packing_time_f, unpack_time_self_f,
              time_for_gpu_pair_f, packing_time_pair_f, unpacking_time_pair_f,
              time_for_gpu_g, packing_time_g, unpack_time_self_g,
              time_for_gpu_pair_g, packing_time_pair_g, unpacking_time_pair_f);

    else
      fprintf(fgpu_steps,
              "%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, "
              "%e, %e\n",
              time_for_density_gpu, packing_time, unpack_time_self,
              time_for_density_gpu_pair, packing_time_pair, unpacking_time_pair,
              time_for_gpu_f, packing_time_f, unpack_time_self_f,
              time_for_gpu_pair_f, packing_time_pair_f, unpacking_time_pair_f,
              time_for_gpu_g, packing_time_g, unpack_time_self_g,
              time_for_gpu_pair_g, packing_time_pair_g, unpacking_time_pair_f);
      //////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////

#else  // No GPU offload
    if (r->cpuid == 0 && engine_rank == 0)
      fprintf(fgpu_steps,
              "CPU TIME SELF, CPU TIME PAIR, "
              "CPU TIME SELF F, CPU TIME PAIR F, CPU TIME SELF G, CPU TIME "
              "PAIR G\n "
              "%e, %e, %e, %e, %e, %e\n",
              time_for_density_cpu, time_for_density_cpu_pair, time_for_cpu_f,
              time_for_cpu_pair_f, time_for_cpu_g, time_for_cpu_pair_g);

    else
      fprintf(fgpu_steps, "%e, %e, %e, %e, %e, %e,\n", time_for_density_cpu,
              time_for_density_cpu_pair, time_for_cpu_f, time_for_cpu_pair_f,
              time_for_cpu_g, time_for_cpu_pair_g);
#endif
    //    }
    fflush(fgpu_steps);
    fclose(fgpu_steps);
#endif  // DUMPTIMINGS
    time_for_density_cpu = 0.0;
    time_for_density_gpu = 0.0;
    time_for_density_cpu_pair = 0.0;
    time_for_density_gpu_pair = 0.0;
    time_for_density_cpu_sub = 0.0;
    tot_time_for_hard_memcpys = 0.0;
    tasks_done_gpu = 0;
    tasks_done_cpu = 0;
    tasks_done_gpu_inc = 0;
    if (ghost_in > 0)
      fprintf(stderr, "total tasks not done on GPU %i is %i\n", r->cpuid,
              ghost_in);
    packed_self = 0;
    packed_pair = 0;
    packed_self_f = 0;
    packed_pair_f = 0;
    packed_self_g = 0;
    packed_pair_g = 0;
    density = 0;
    density_sub = 0;
    unpacked = 0;
    //	if(step == 2)cudaProfilerStop();
    //	if(step == 2)exit(0);
    //	  size_t free_byte ;
    //	  size_t total_byte ;
    //	  cudaError_t cuda_status = cudaMemGetInfo( &free_byte,
    //&total_byte ) ; 	  double free = (double)free_byte; 	  double
    // available = (double)total_byte; 	  double used = (available - free);
    // fprintf(stderr, "Used %f GB GPU memory\n", used/1e9);
    /* Wait at the wait barrier. */
    //    swift_barrier_wait(&e->wait_barrier);
  }
  // Free all data
  //  cudaFree(d_tid_p);
  //  cudaFree(d_id);
  //  cudaFree(d_x_p);
  //  cudaFree(d_y_p);
  //  cudaFree(d_z_p);
  //  cudaFree(d_ux);
  //  cudaFree(d_uy);
  //  cudaFree(d_uz);
  //  cudaFree(d_a_hydrox);
  //  cudaFree(d_a_hydroy);
  //  cudaFree(d_a_hydroz);
  //  cudaFree(d_mass);
  //  cudaFree(d_h);
  //  cudaFree(d_u);
  //  cudaFree(d_u_dt);
  //  cudaFree(d_rho);
  //  cudaFree(d_SPH_sum);
  //  cudaFree(d_locx);
  //  cudaFree(d_locy);
  //  cudaFree(d_locz);
  //  cudaFree(d_widthx);
  //  cudaFree(d_widthy);
  //  cudaFree(d_widthz);
  //  cudaFree(d_h_max);
  //  cudaFree(d_count_p);
  //  cudaFree(d_wcount);
  //  cudaFree(d_wcount_dh);
  //  cudaFree(d_rho_dh);
  //  cudaFree(d_rot_ux);
  //  cudaFree(d_rot_uy);
  //  cudaFree(d_rot_uz);
  //  cudaFree(d_div_v);
  //  cudaFree(d_div_v_previous_step);
  //  cudaFree(d_alpha_visc);
  //  cudaFree(d_v_sig);
  //  cudaFree(d_laplace_u);
  //  cudaFree(d_alpha_diff);
  //  cudaFree(d_f);
  //  cudaFree(d_soundspeed);
  //  cudaFree(d_h_dt);
  //  cudaFree(d_balsara);
  //  cudaFree(d_pressure);
  //  cudaFree(d_alpha_visc_max_ngb);
  //  cudaFree(d_time_bin);
  //  cudaFree(d_wakeup);
  //  cudaFree(d_min_ngb_time_bin);
  //  cudaFree(d_to_be_synchronized);
  //  cudaFree(tid_p);
  //  cudaFree(id);
  //  cudaFree(mass);
  //  cudaFree(h);
  //  cudaFree(u);
  //  cudaFree(u_dt);
  //  cudaFree(rho);
  //  cudaFree(SPH_sum);
  //  cudaFree(x_p);
  //  cudaFree(y_p);
  //  cudaFree(z_p);
  //  cudaFree(ux);
  //  cudaFree(uy);
  //  cudaFree(uz);
  //  cudaFree(a_hydrox);
  //  cudaFree(a_hydroy);
  //  cudaFree(a_hydroz);
  //  cudaFree(locx);
  //  cudaFree(locy);
  //  cudaFree(locz);
  //  cudaFree(widthx);
  //  cudaFree(widthy);
  //  cudaFree(widthz);
  //  cudaFree(h_max);
  //  cudaFree(count_p);
  //  cudaFree(wcount);
  //  cudaFree(wcount_dh);
  //  cudaFree(rho_dh);
  //  cudaFree(rot_ux);
  //  cudaFree(rot_uy);
  //  cudaFree(rot_uz);
  //  cudaFree(div_v);
  //  cudaFree(div_v_previous_step);
  //  cudaFree(alpha_visc);
  //  cudaFree(v_sig);
  //  cudaFree(laplace_u);
  //  cudaFree(alpha_diff);
  //  cudaFree(f);
  //  cudaFree(soundspeed);
  //  cudaFree(h_dt);
  //  cudaFree(balsara);
  //  cudaFree(pressure);
  //  cudaFree(alpha_visc_max_ngb);
  //  cudaFree(time_bin);
  //  cudaFree(wakeup);
  //  cudaFree(min_ngb_time_bin);
  //  cudaFree(to_be_synchronized);
  //  cudaFree(partid_p);
  //  cudaFree(d_task_first_part);
  //  cudaFree(d_task_last_part);
  //  cudaFree(task_first_part_self_dens);
  //  cudaFree(task_last_part_self_dens);
  //  cudaFree(task_first_part_pair_ci);
  //  cudaFree(task_last_part_pair_ci);
  //  cudaFree(task_first_part_pair_cj);
  //  cudaFree(task_last_part_pair_cj);
  //  cudaFree(d_bundle_first_part_self_dens);
  //  cudaFree(d_bundle_last_part_self_dens);
  //  cudaFree(bundle_first_part_self_dens);
  //  cudaFree(bundle_last_part_self_dens);
  //  cudaFree(bundle_first_part_pair_ci);
  //  cudaFree(bundle_last_part_pair_ci);
  //  cudaFree(bundle_first_part_pair_cj);
  //  cudaFree(bundle_last_part_pair_cj);
  //  free(ci_list_self_dens);
  //  free(ci_list_pair);
  //  free(cj_list_pair);

  /* Be kind, rewind. */
  return NULL;
}

#endif  // WITH_CUDA

