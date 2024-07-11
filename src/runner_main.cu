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

#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* This object's header. */
#include "runner.h"
/* Local headers. */
#include "active.h"
#include "engine.h"
#include "scheduler.h"
#include "space_getsid.h"
#include "timers.h"

/* Import the gravity loop functions. */
#include "runner_doiact_grav.h"

/* Import the density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the gradient loop functions (if required). */
#ifdef EXTRA_HYDRO_LOOP
#define FUNCTION gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_GRADIENT
#include "runner_doiact_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP
#endif

/* Import the force loop functions. */
#define FUNCTION force
#define FUNCTION_TASK_LOOP TASK_LOOP_FORCE
#include "runner_doiact_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the limiter loop functions. */
#define FUNCTION limiter
#define FUNCTION_TASK_LOOP TASK_LOOP_LIMITER
#include "runner_doiact_limiter.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the stars density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the stars feedback loop functions. */
#define FUNCTION feedback
#define FUNCTION_TASK_LOOP TASK_LOOP_FEEDBACK
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the black hole density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_black_holes.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the black hole feedback loop functions. */
#define FUNCTION swallow
#define FUNCTION_TASK_LOOP TASK_LOOP_SWALLOW
#include "runner_doiact_black_holes.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the black hole feedback loop functions. */
#define FUNCTION feedback
#define FUNCTION_TASK_LOOP TASK_LOOP_FEEDBACK
#include "runner_doiact_black_holes.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import radiative transfer loop functions. */
#define FUNCTION inject
#include "runner_doiact_rt.h"
#undef FUNCTION

/* Import the sink compute formation loop functions. */
#define FUNCTION compute_formation
#define FUNCTION_TASK_LOOP TASK_LOOP_SINK_FORMATION
#include "runner_doiact_sinks.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

#ifdef __cplusplus
}
#endif
/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.

/* CUDA Header */
#ifdef WITH_CUDA
#include "./cuda/BLOCK_SIZE.h"
#include "./cuda/Data_and_GPU_prep_functions.cu"
#include "./cuda/cell_gpu.h"
#include "./cuda/cuda_headers.h"
#include "./cuda/tasks_gpu.h"
#include "runner_gpu_pack_functions.h"
#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline cudaError_t checkCuda(cudaError_t result) {
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}

inline void gpuErrchk(cudaError_t code) {
#define __FILE__ __LINE__
  inline void gpuAssert(cudaError_t code, const char *file, int line) {
    int abort = 0;
    if (code != cudaSuccess) {
      //			fprintf( stderr, "cudaCheckError() failed at
      //%s:%i : %s\n",
      //                 file, line, cudaGetErrorString( code ) );
      abort = 1;
      fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
              line);
      if (abort)
        exit(code);
    }
  }
}

void *runner_main_cuda(void *data) {

  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);

  int counter = 0;
  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);

    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;

    /* Loop while there are tasks... */
    while (1) {

      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, prev);
        TIMER_TOC(timer_gettask);
        /* Did I get anything? */
        /* If yes, can I change it for another task if its not what I want? */
        if (t == NULL)
          break;
      }
      /*if task type is force self*/
      //		if(r->cpuid==0 && t->subtype != task_subtype_force &&
      //t->type != task_type_self){ 			printf("Task ignored \n"); 			continue; //This
      //doesn't work, obviously! Task gets left undone as it is not re-scheduled
      //for another runner!
      //		}
      /* Get the cells. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;
#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        struct cell *ci_temp = ci;
        struct cell *cj_temp = cj;
        double shift[3];
        t->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
      } else {
        t->sid = -1;
      }
#endif
      if (t->type != task_type_self && t->type != task_type_sub_self)
        printf("type is not self or sub-self\n");
      if (t->subtype != task_subtype_density &&
          t->subtype != task_subtype_gradient &&
          t->subtype != task_subtype_force)
        printf("subtype is not force, gradient or density\n");
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
       * i.e., rather than just one r->t, r would have a nested bunch of structs
       * r->t[size=nStreams].*/
      r->t = t;
#endif
      /*****Initialise GPU copies of variables*****/
      const int count = ci->hydro.count;
      const int countj = cj->hydro.count;
      struct cell_gpu *ci_gpu = malloc(sizeof(*ci_gpu));
      struct cell_gpu *cj_gpu = malloc(sizeof(*cj_gpu));
      struct part_gpu *parts_gpu =
          (struct part_gpu *)malloc(count * sizeof(struct part_gpu));
      struct part_gpu *parts_gpuj =
          (struct part_gpu *)malloc(count * sizeof(struct part_gpu));
      /*****Copy variables for cell i (self interaction)*****/
      ci_gpu->hydro.count = count;
      ci_gpu->hydro.h_max = ci->hydro.h_max;
      for (int d = 0; d < 3; d++) {
        ci_gpu->width[d] = ci->width[d];
        ci_gpu->loc[d] = ci->loc[d];
      }
      for (int p = 0; p < count; p++) {
        parts_gpu[p].id = ci->hydro.parts[p].id;
        for (int d = 0; d < 3; d++) {
          parts_gpu[p].x[d] = ci->hydro.parts[p].x[d];
          parts_gpu[p].v[d] = ci->hydro.parts[p].v[d];
          parts_gpu[p].a_hydro[d] = ci->hydro.parts[p].a_hydro[d];
        }
        parts_gpu[p].mass = ci->hydro.parts[p].mass;
        parts_gpu[p].h = ci->hydro.parts[p].h;
        parts_gpu[p].u = ci->hydro.parts[p].u;
        parts_gpu[p].u_dt = ci->hydro.parts[p].u_dt;
        parts_gpu[p].rho = ci->hydro.parts[p].rho;
        parts_gpu[p].div_v = ci->hydro.parts[p].viscosity.div_v;
        parts_gpu[p].div_v_previous_step =
            ci->hydro.parts[p].viscosity.div_v_previous_step;
        parts_gpu[p].alpha_visc = ci->hydro.parts[p].viscosity.alpha;
        parts_gpu[p].v_sig = ci->hydro.parts[p].viscosity.v_sig;
        parts_gpu[p].laplace_u = ci->hydro.parts[p].diffusion.laplace_u;
        parts_gpu[p].alpha_diff = ci->hydro.parts[p].diffusion.alpha;
        parts_gpu[p].f = ci->hydro.parts[p].force.f;
        parts_gpu[p].soundspeed = ci->hydro.parts[p].force.soundspeed;
        parts_gpu[p].h_dt = ci->hydro.parts[p].force.h_dt;
        parts_gpu[p].balsara = ci->hydro.parts[p].force.balsara;
        parts_gpu[p].pressure = ci->hydro.parts[p].force.pressure;
        parts_gpu[p].time_bin = ci->hydro.parts[p].time_bin;
        parts_gpu[p].wakeup = ci->hydro.parts[p].limiter_data.wakeup;
        parts_gpu[p].min_ngb_time_bin =
            ci->hydro.parts[p].limiter_data.min_ngb_time_bin;
        parts_gpu[p].to_be_synchronized =
            ci->hydro.parts[p].limiter_data.to_be_synchronized;
        parts_gpu[p].wcount = ci->hydro.parts[p].density.wcount;
        parts_gpu[p].wcount_dh = ci->hydro.parts[p].density.wcount_dh;
        parts_gpu[p].rho_dh = ci->hydro.parts[p].density.rho_dh;
        parts_gpu[p].div_v = ci->hydro.parts[p].viscosity.div_v;
        parts_gpu[p].rot_v[0] = ci->hydro.parts[p].density.rot_v[0];
        parts_gpu[p].rot_v[1] = ci->hydro.parts[p].density.rot_v[1];
        parts_gpu[p].rot_v[2] = ci->hydro.parts[p].density.rot_v[2];
        parts_gpu[p].SPH_sum = 0.f;
      }
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      switch (t->type) {
      case task_type_self:
        if (t->subtype == task_subtype_density) {
          counter++;
          printf("Counter is %i\n", counter);
          const char *loop_type = "density";
          double t1 = clock();
          //		     	 runner_doself1_branch_density(r, ci);
          double t2 = clock();

          printf("CPU time is %f\n", (t2 - t1) / CLOCKS_PER_SEC);
          double tstart = clock();
          // FILE *fp;
          // fp = fopen("./res_CPU_density.txt", "w");
          // for(int p=0; p<count; p++){
          //	float xx=parts_gpu[p].x[0], yy=parts_gpu[p].x[1],
          //zz=parts_gpu[p].x[2]; 	fprintf(fp, "%f %f %f %f %f %f %f %f %f\n",
          //xx, yy, zz, ci->hydro.parts[p].rho,
          //	 ci->hydro.parts[p].density.rho_dh,
          //ci->hydro.parts[p].density.wcount,
          //	 ci->hydro.parts[p].density.wcount_dh,
          //ci->hydro.parts[p].viscosity.div_v 	 ,
          //ci->hydro.parts[p].density.rot_v[0]);
          //}
          // fclose(fp);
          double tend = clock();
          // printf("time setting up arrays is %f\n",(tend-tstart)/
          // CLOCKS_PER_SEC);

          t1 = clock();
          int numBlocks = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
          /*define device structs*/
          struct cell_gpu *d_ci_gpu;
          struct part_gpu *d_parts;
          int size_parts = count * sizeof(struct part_gpu);
          int size_cell = sizeof(struct cell_gpu);
          cudaMalloc((void **)&d_parts, size_parts);
          cudaMalloc((void **)&d_ci_gpu, size_cell);
          t2 = clock();
          cudaMemcpy(d_ci_gpu, ci_gpu, size_cell, cudaMemcpyHostToDevice);
          cudaMemcpy(d_parts, parts_gpu, size_parts, cudaMemcpyHostToDevice);
          double t3 = clock();
          launch_cuda_kernel(d_ci_gpu, d_parts, numBlocks, d_a, d_H, loop_type);
          double t4 = clock();
          cudaMemcpy(parts_gpu, d_parts, size_parts, cudaMemcpyDeviceToHost);
          cudaMemcpy(ci_gpu, d_ci_gpu, size_cell, cudaMemcpyDeviceToHost);
          double t5 = clock();
          cudaDeviceSynchronize();
          cudaFree(d_parts), cudaFree(d_ci_gpu);
          for (int p = 0; p < count; p++) {
            ci->hydro.parts[p].density.wcount = parts_gpu[p].wcount;
            ci->hydro.parts[p].density.wcount_dh = parts_gpu[p].wcount_dh;
            ci->hydro.parts[p].density.rho_dh = parts_gpu[p].rho_dh;
            ci->hydro.parts[p].rho = parts_gpu[p].rho;
            ci->hydro.parts[p].viscosity.div_v = parts_gpu[p].div_v;
            ci->hydro.parts[p].density.rot_v[0] = parts_gpu[p].rot_v[0];
            ci->hydro.parts[p].density.rot_v[1] = parts_gpu[p].rot_v[1];
            ci->hydro.parts[p].density.rot_v[2] = parts_gpu[p].rot_v[2];
          }
          double t6 = clock();
          printf("time density loop for malloc %f\n",
                 (t2 - t1) / CLOCKS_PER_SEC);
          printf("time for memcpy H2D %f\n", (t3 - t2) / CLOCKS_PER_SEC);
          printf("time for kernel %f\n", (t4 - t3) / CLOCKS_PER_SEC);
          printf("time for memcpy D2H %f\n", (t5 - t4) / CLOCKS_PER_SEC);
          printf("Total GPU time is %f\n\n", (t6 - t1) / CLOCKS_PER_SEC);
          //             GPU_runner_doself1_branch_gradient(ci_gpu, parts_gpu);
          // fp = fopen("./res_GPU_density.txt", "w");
          // for(int p=0; p<count; p++){
          //	float xx=parts_gpu[p].x[0], yy=parts_gpu[p].x[1],
          //zz=parts_gpu[p].x[2]; 	fprintf(fp, "%f %f %f %f %f %f %f %f %f\n",
          //xx, yy, zz, parts_gpu[p].rho, 	parts_gpu[p].rho_dh,
          //parts_gpu[p].wcount, parts_gpu[p].wcount_dh, 	parts_gpu[p].div_v,
          //parts_gpu[p].rot_v[0]);
          //}
          // fclose(fp);
          if (counter == 1000)
            exit(0);
          free(parts_gpu), free(ci_gpu), free(cj_gpu);
        }
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient) {
          //  				 counter++;
          //  				 printf("Counter is %i",counter);
          const char *loop_type = "gradient";
          double t1 = clock();
          //             runner_doself1_branch_gradient(r, ci);
          double t2 = clock();
          int in = 1;
          if (in == 1) {
            //             printf("CPU time is %f\n",(t2-t1)/ CLOCKS_PER_SEC);
            double tstart = clock();
            FILE *fp;
            //             fp = fopen("./res_CPU_gradient.txt", "w");
            //             for(int p=0; p<count; p++){
            //             	float xx=parts_gpu[p].x[0],
            //             yy=parts_gpu[p].x[1], zz=parts_gpu[p].x[2];
            //					fprintf(fp, "%f %f %f %f %f %f\n", xx, yy, zz,
            //ci->hydro.parts[p].viscosity.v_sig,
            //ci->hydro.parts[p].diffusion.laplace_u,
            //ci->hydro.parts[p].force.alpha_visc_max_ngb);
            //             }
            //             fclose(fp);
            double tend = clock();
            //             printf("time setting up arrays is
            //             %f\n",(tend-tstart)/ CLOCKS_PER_SEC);

            t1 = clock();
            int numBlocks = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
            /*define device structs*/
            struct cell_gpu *d_ci_gpu;
            struct part_gpu *d_parts;
            int size_parts = count * sizeof(struct part_gpu);
            int size_cell = sizeof(struct cell_gpu);
            cudaMalloc((void **)&d_parts, size_parts);
            cudaMalloc((void **)&d_ci_gpu, size_cell);
            t2 = clock();
            cudaMemcpy(d_ci_gpu, ci_gpu, size_cell, cudaMemcpyHostToDevice);
            cudaMemcpy(d_parts, parts_gpu, size_parts, cudaMemcpyHostToDevice);
            double t3 = clock();
            launch_cuda_kernel(d_ci_gpu, d_parts, numBlocks, d_a, d_H,
                               loop_type);
            double t4 = clock();
            cudaMemcpy(parts_gpu, d_parts, size_parts, cudaMemcpyDeviceToHost);
            cudaMemcpy(ci_gpu, d_ci_gpu, size_cell, cudaMemcpyDeviceToHost);
            double t5 = clock();
            cudaDeviceSynchronize();
            cudaFree(d_parts), cudaFree(d_ci_gpu);
            for (int p = 0; p < count; p++) {
              ci->hydro.parts[p].viscosity.v_sig = parts_gpu[p].v_sig;
              ci->hydro.parts[p].diffusion.laplace_u = parts_gpu[p].laplace_u;
              ci->hydro.parts[p].force.alpha_visc_max_ngb =
                  parts_gpu[p].alpha_visc_max_ngb;
            }
            double t6 = clock();
            //				 printf("time for malloc %f\n",(t2-t1)/
            //CLOCKS_PER_SEC); 				 printf("time for memcpy H2D %f\n",(t3-t2)/
            //CLOCKS_PER_SEC); 				 printf("time for kernel %f\n",(t4-t3)/
            //CLOCKS_PER_SEC); 				 printf("time for memcpy D2H %f\n",(t5-t4)/
            //CLOCKS_PER_SEC); 				 printf("Total GPU time is %f\n\n",(t6-t1)/
            //CLOCKS_PER_SEC);
            //             GPU_runner_doself1_branch_gradient(ci_gpu,
            //             parts_gpu); fp = fopen("./res_GPU_gradient.txt",
            //             "w"); for(int p=0; p<count; p++){ 	float
            //             xx=parts_gpu[p].x[0], yy=parts_gpu[p].x[1],
            //             zz=parts_gpu[p].x[2];
            //					fprintf(fp, "%f %f %f %f %f %f\n", xx, yy, zz,
            //parts_gpu[p].v_sig, parts_gpu[p].laplace_u,
            //parts_gpu[p].alpha_visc_max_ngb);
            //             }
            //             fclose(fp);
            //             for(int p=0; p<count; p++){
            //             	if(parts_gpu[p].v_sig>10.f){
            //             		printf("vsig in gradient is
            //             %f\n",parts_gpu[p].v_sig); 		exit(0);
            //             	}
            //             }
            free(parts_gpu), free(ci_gpu), free(cj_gpu);
          }
        }
#endif
        else if (t->subtype == task_subtype_force) {
          /*cudafree is used here to initialise cuda and take the hit of mapping
             GPU memory, etc., before we start doing SPH calculations. The
             initialisation takes a long time so we do it here before we start
             timing*/
          cudaFree(0);
          counter++;
          const char *loop_type = "force";
          double t1 = clock();
          //             if (!cell_is_active_hydro(ci, e)) break;
          //             runner_doself2_branch_force(r, ci);
          double t2 = clock();
          //             if(count==0)continue;
          //             printf("CPU time is %f\n",(t2-t1)/ CLOCKS_PER_SEC);
          double tstart = clock();
          FILE *fp;
          //             fp = fopen("./res_CPU_force.txt", "w");
          //             for(int p=0; p<count; p++){
          //             	float xx=parts_gpu[p].x[0],
          //             yy=parts_gpu[p].x[1], zz=parts_gpu[p].x[2];
          //						fprintf(fp, "%f %f %f %f %f %f %f %f %f\n",
          //xx, yy, zz, ci->hydro.parts[p].a_hydro[0], 						ci->hydro.parts[p].u_dt,
          //ci->hydro.parts[p].force.h_dt, 						ci->hydro.parts[p].viscosity.v_sig),
          //ci->hydro.parts[p].limiter_data.min_ngb_time_bin,
          //						ci->hydro.parts[p].time_bin;
          //            }
          //             fclose(fp);
          double tend = clock();
          //             printf("time setting up arrays is %f\n",(tend-tstart)/
          //             CLOCKS_PER_SEC);

          t1 = clock();
          int numBlocks = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
          /*define device structs*/
          struct cell_gpu *d_ci_gpu;
          struct part_gpu *d_parts;
          int size_parts = count * sizeof(struct part_gpu);
          int size_cell = sizeof(struct cell_gpu);
          cudaMalloc((void **)&d_parts, size_parts);
          cudaMalloc((void **)&d_ci_gpu, size_cell);
          t2 = clock();
          cudaMemcpy(d_ci_gpu, ci_gpu, size_cell, cudaMemcpyHostToDevice);
          cudaMemcpy(d_parts, parts_gpu, size_parts, cudaMemcpyHostToDevice);
          double t3 = clock();
          launch_cuda_kernel(d_ci_gpu, d_parts, numBlocks, d_a, d_H, loop_type);
          double t4 = clock();
          cudaMemcpy(parts_gpu, d_parts, size_parts, cudaMemcpyDeviceToHost);
          cudaMemcpy(ci_gpu, d_ci_gpu, size_cell, cudaMemcpyDeviceToHost);
          double t5 = clock();
          cudaDeviceSynchronize();
          cudaFree(d_parts), cudaFree(d_ci_gpu);
          for (int p = 0; p < count; p++) {
            //				 	 float a_hydro;
            for (int d = 0; d < 3; d++) {
              //					 	 a_hydro=parts_gpu[p].a_hydro[d];
              ci->hydro.parts[p].a_hydro[d] = parts_gpu[p].a_hydro[d];
            }
            ci->hydro.parts[p].u_dt = parts_gpu[p].u_dt;
            ci->hydro.parts[p].viscosity.v_sig = parts_gpu[p].v_sig;
            ci->hydro.parts[p].force.h_dt = parts_gpu[p].h_dt;
            ci->hydro.parts[p].time_bin = parts_gpu[p].time_bin;
            ci->hydro.parts[p].limiter_data.min_ngb_time_bin =
                parts_gpu[p].min_ngb_time_bin;
          }
          double t6 = clock();
          //				 printf("time for malloc %f\n",(t2-t1)/
          //CLOCKS_PER_SEC); 				 printf("time for memcpy H2D %f\n",(t3-t2)/
          //CLOCKS_PER_SEC); 				 printf("time for kernel %f\n",(t4-t3)/
          //CLOCKS_PER_SEC); 				 printf("time for memcpy D2H %f\n",(t5-t4)/
          //CLOCKS_PER_SEC); 				 printf("Total GPU time is %f\n\n",(t6-t1)/
          //CLOCKS_PER_SEC);
          //             GPU_runner_doself1_branch_gradient(ci_gpu, parts_gpu);
          //				 printf("Counter is %i", counter);
          //             fp = fopen("./res_GPU_force.txt", "w");
          //             for(int p=0; p<count; p++){
          //            	float xx=ci->hydro.parts[p].x[0],
          //            yy=ci->hydro.parts[p].x[1], zz=ci->hydro.parts[p].x[2];
          ////             	fprintf(fp, "%f %f %f %f %f %f %f %f\n", xx, yy,
          ///zz, ci->hydro.parts[p].a_hydro[0], /
          ///ci->hydro.parts[p].u_dt, ci->hydro.parts[p].force.h_dt, /
          ///ci->hydro.parts[p].viscosity.v_sig),
          ///ci->hydro.parts[p].limiter_data.min_ngb_time_bin;
          //					fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", xx,
          //yy, zz, parts_gpu[p].a_hydro[0] 					, parts_gpu[p].u_dt,
          //parts_gpu[p].h_dt, parts_gpu[p].v_sig,
          //parts_gpu[p].min_ngb_time_bin, 					parts_gpu[p].time_bin);
          //           }
          //           fclose(fp);
          //             for(int p=0; p<count; p++){
          //             	if(parts_gpu[p].v_sig>10.f){
          //             		printf("vsig in force is
          //             %f\n",parts_gpu[p].v_sig); 		exit(0);
          //             	}
          //             }
          free(parts_gpu), free(ci_gpu), free(cj_gpu);

          //             exit(0);
          //             if(counter==1)exit(0);
        } else if (t->subtype == task_subtype_limiter)
          runner_doself1_branch_limiter(r, ci);
        else if (t->subtype == task_subtype_grav)
          runner_doself_recursive_grav(r, ci, 1);
        else if (t->subtype == task_subtype_external_grav)
          runner_do_grav_external(r, ci, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_doself_branch_stars_density(r, ci);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_doself_branch_rt_inject(r, ci, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_doself_branch_sinks_compute_formation(r, ci);
        else
          error("Unknown/invalid task subtype (%s).",
                subtaskID_names[t->subtype]);
        break;

      case task_type_pair:
        /*****Copy variables for cell j (if it is a pair interaction)*****/

        //				cj_gpu->hydro.count=countj;
        //				printf("count is %i %i\n",count,
        //countj); 				cj_gpu->hydro.h_max=cj->hydro.h_max; 				for(int d=0; d<3; d++){
        //					cj_gpu->width[d]=cj->width[d];
        //				   cj_gpu->loc[d]=cj->loc[d];
        //				}
        //				for(int p=0; p<count; p++){
        //					parts_gpuj[p].id=cj->hydro.parts[p].id;
        //					for(int d=0; d<3; d++){
        //				   	parts_gpuj[p].x[d]=cj->hydro.parts[p].x[d];
        //				      parts_gpuj[p].v[d]=cj->hydro.parts[p].v[d];
        //				      parts_gpuj[p].a_hydro[d]=cj->hydro.parts[p].a_hydro[d];
        //				   }
        //				   parts_gpuj[p].mass=cj->hydro.parts[p].mass;
        //				   parts_gpuj[p].h=cj->hydro.parts[p].h;
        //				   parts_gpuj[p].u=cj->hydro.parts[p].u;
        //				   parts_gpuj[p].u_dt=cj->hydro.parts[p].u_dt;
        //				   parts_gpuj[p].rho=1.f;//ci->hydro.parts[p].rho;
        //				   parts_gpuj[p].SPH_sum=0.f;
        //				}
        if (t->subtype == task_subtype_density)
          runner_dopair1_branch_density(r, ci, cj);
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_dopair1_branch_gradient(r, ci, cj);
#endif
        else if (t->subtype == task_subtype_force)
          runner_dopair2_branch_force(r, ci, cj);
        else if (t->subtype == task_subtype_limiter)
          runner_dopair1_branch_limiter(r, ci, cj);
        else if (t->subtype == task_subtype_grav)
          runner_dopair_recursive_grav(r, ci, cj, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_dopair_branch_stars_density(r, ci, cj);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_dopair_branch_rt_inject(r, ci, cj, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_dopair_branch_sinks_compute_formation(r, ci, cj);
        else
          error("Unknown/invalid task subtype (%s/%s).", taskID_names[t->type],
                subtaskID_names[t->subtype]);
        break;

      case task_type_sub_self:
        if (t->subtype == task_subtype_density)
          runner_dosub_self1_density(r, ci, 1);
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_dosub_self1_gradient(r, ci, 1);
#endif
        else if (t->subtype == task_subtype_force)
          runner_dosub_self2_force(r, ci, 1);
        else if (t->subtype == task_subtype_limiter)
          runner_dosub_self1_limiter(r, ci, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_dosub_self_stars_density(r, ci, 1);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_dosub_self_rt_inject(r, ci, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_dosub_self_sinks_compute_formation(r, ci, 1);
        else
          error("Unknown/invalid task subtype (%s/%s).", taskID_names[t->type],
                subtaskID_names[t->subtype]);
        break;

      case task_type_sub_pair:
        /*****Copy variables for cell j (if it is a pair interaction)*****/

        //				cj_gpu->hydro.count=countj;
        //				printf("count is %i %i\n",count,
        //countj); 				cj_gpu->hydro.h_max=cj->hydro.h_max; 				for(int d=0; d<3; d++){
        //					cj_gpu->width[d]=cj->width[d];
        //				   cj_gpu->loc[d]=cj->loc[d];
        //				}
        //				for(int p=0; p<count; p++){
        //					parts_gpuj[p].id=cj->hydro.parts[p].id;
        //					for(int d=0; d<3; d++){
        //				   	parts_gpuj[p].x[d]=cj->hydro.parts[p].x[d];
        //				      parts_gpuj[p].v[d]=cj->hydro.parts[p].v[d];
        //				      parts_gpuj[p].a_hydro[d]=cj->hydro.parts[p].a_hydro[d];
        //				   }
        //				   parts_gpuj[p].mass=cj->hydro.parts[p].mass;
        //				   parts_gpuj[p].h=cj->hydro.parts[p].h;
        //				   parts_gpuj[p].u=cj->hydro.parts[p].u;
        //				   parts_gpuj[p].u_dt=cj->hydro.parts[p].u_dt;
        //				   parts_gpuj[p].rho=1.f;//ci->hydro.parts[p].rho;
        //				   parts_gpuj[p].SPH_sum=0.f;
        //				}
        if (t->subtype == task_subtype_density)
          runner_dosub_pair1_density(r, ci, cj, 1);
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_dosub_pair1_gradient(r, ci, cj, 1);
#endif
        else if (t->subtype == task_subtype_force)
          runner_dosub_pair2_force(r, ci, cj, 1);
        else if (t->subtype == task_subtype_limiter)
          runner_dosub_pair1_limiter(r, ci, cj, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_dosub_pair_stars_density(r, ci, cj, 1);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_dosub_pair_rt_inject(r, ci, cj, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_dosub_pair_sinks_compute_formation(r, ci, cj, 1);
        else
          error("Unknown/invalid task subtype (%s/%s).", taskID_names[t->type],
                subtaskID_names[t->subtype]);
        break;

      case task_type_sort:
        /* Cleanup only if any of the indices went stale. */
        runner_do_hydro_sort(
            r, ci, t->flags,
            ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin, 1);
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
      case task_type_logger:
        runner_do_logger(r, ci, 1);
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
#ifdef WITH_MPI
      case task_type_send:
        if (t->subtype == task_subtype_tend_part) {
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_gpart) {
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_spart) {
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_bpart) {
          free(t->buff);
        } else if (t->subtype == task_subtype_sf_counts) {
          free(t->buff);
        } else if (t->subtype == task_subtype_part_swallow) {
          free(t->buff);
        } else if (t->subtype == task_subtype_bpart_merger) {
          free(t->buff);
        }
        break;
      case task_type_recv:
        if (t->subtype == task_subtype_tend_part) {
          cell_unpack_end_step_hydro(ci, (struct pcell_step_hydro *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_gpart) {
          cell_unpack_end_step_grav(ci, (struct pcell_step_grav *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_spart) {
          cell_unpack_end_step_stars(ci, (struct pcell_step_stars *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_bpart) {
          cell_unpack_end_step_black_holes(
              ci, (struct pcell_step_black_holes *)t->buff);
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
        } else if (t->subtype == task_subtype_part_swallow) {
          cell_unpack_part_swallow(ci, (struct black_holes_part_data *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_bpart_merger) {
          cell_unpack_bpart_swallow(ci,
                                    (struct black_holes_bpart_data *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_limiter) {
          runner_do_recv_part(r, ci, 0, 1);
        } else if (t->subtype == task_subtype_gpart) {
          runner_do_recv_gpart(r, ci, 1);
        } else if (t->subtype == task_subtype_spart) {
          runner_do_recv_spart(r, ci, 1, 1);
        } else if (t->subtype == task_subtype_bpart_rho) {
          runner_do_recv_bpart(r, ci, 1, 1);
        } else if (t->subtype == task_subtype_bpart_swallow) {
          runner_do_recv_bpart(r, ci, 0, 1);
        } else if (t->subtype == task_subtype_bpart_feedback) {
          runner_do_recv_bpart(r, ci, 0, 1);
        } else if (t->subtype == task_subtype_multipole) {
          cell_unpack_multipoles(ci, (struct gravity_tensors *)t->buff);
          free(t->buff);
        } else {
          error("Unknown/invalid task subtype (%d).", t->subtype);
        }
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
      case task_type_stars_resort:
        runner_do_stars_resort(r, t->ci, 1);
        break;
      case task_type_sink_formation:
        runner_do_sink_formation(r, t->ci);
        break;
      case task_type_fof_self:
        runner_do_fof_self(r, t->ci, 1);
        break;
      case task_type_fof_pair:
        runner_do_fof_pair(r, t->ci, t->cj, 1);
        break;
      case task_type_rt_ghost1:
        runner_do_rt_ghost1(r, t->ci, 1);
        break;
      default:
        error("Unknown/invalid task type (%d).", t->type);
      }

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
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////****same as runner_main_cuda but using streams to bunch tasks together b4
///sending to GPU***********/////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void *runner_main_cuda_streams(void *data) {

  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  cudaError_t cu_error;
  //  cudaMemPool_t memPool;
  //  gpuErrchk( cudaDeviceGetDefaultMemPool(&memPool, devId) );
  //  int maxmem = 0.9f*prop.totalGlobalMem;
  //  gpuErrchk(cudaMemPoolSetAttribute(memPool,
  //  cudaMemPoolAttrReleaseThreshold, (void *)&maxmem));
  int counter = 0;

  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    struct task **tasks;
    gpuErrchk(cudaMallocHost(
        (void **)&tasks,
        n_streams * sizeof(struct task))); // Pinned allocation on host
    // tasks=malloc(n_streams*sizeof(struct task));
    for (int i = 0; i < n_streams; i++) {
      tasks[i] = NULL;
    }
    //	 struct task* tasks = malloc(n_streams * sizeof(struct task));

    /* Loop while there are tasks... */
    //    struct tasks_self_gpu tasks;
    // for (int  tid=0; tid<n_streams; tid++){
    //	tasks.t[tid] = NULL;
    //}
    while (1) {
      for (int i = 0; i < n_streams; i++) {
        tasks[i] = NULL;
      }
      // for(int tid=0; tid<n_streams; tid++){
      int count_tasks = 0;
      int count_stale = 0;
      ////////////////////////////////////
      /*Grab a bunch of tasks*/
      while (1) { //(count_tasks<n_streams){
        // there's a potential bug here. Code hangs if n_streams set to > 256 in
        // cuda_headers.h
        /* If there's no old task, try to get a new one. */
        if (tasks[count_tasks] == NULL) {
          /* Get the task. */
          TIMER_TIC
          //        fprintf(stderr, "searching for task %i\n", count_stale);
          tasks[count_tasks] = scheduler_gettask(sched, r->qid, prev);
          //        fprintf(stderr, "got one\n,");
          TIMER_TOC(timer_gettask);
          if (tasks[count_tasks] != NULL) {
            count_tasks++;
            //        	count_stale++;
          }
          //        fprintf(stderr, "got one, count tasks is %i\n",
          //        count_tasks); else{ 	count_stale++;
          //        }
        }
        /* If there's an old task, move onto the next one */
        // else count_tasks++;
        count_stale++;
        //      printf("count stale is: %i, count tasks is: %i, n_streams is:
        //      %i\n", count_stale, count_tasks, n_streams);
        if (count_stale >= n_streams) {
          //      	printf("count stale is: %i, count tasks is: %i,
          //      n_streams is: %i\n", count_stale, count_tasks, n_streams);
          break;
        }
      }
      fprintf(stderr, "Count tasks is %i\n", count_tasks);
      //     printf("Got out of loop\n");
      ////////////////////////////////////

      ////////////////////////////////////
      /*Get cell data for each task*/
      for (int tid = 1; tid < count_tasks; tid++) {
        struct cell *c_gpu = tasks[tid]->ci;
        //     	printf("Cell pos for task %i is %f %f %f\n",tid, c_gpu->loc[0],
        //     c_gpu->loc[1], c_gpu->loc[2]);
      }

      /* Get the cells. */
      struct cell **ci_list;
      int max_count = 0;
      gpuErrchk(cudaMallocHost(
          (void **)&ci_list,
          count_tasks * sizeof(struct cell))); // Pinned allocation on host
      //		ci_list=malloc(count_tasks*(sizeof(struct
      //cell_gpu)+sizeof(struct cell_hydro_gpu)));
      for (int tid = 0; tid < count_tasks; tid++) {
        ci_list[tid] = tasks[tid]->ci;
        const int count = ci_list[tid]->hydro.count;
        max_count = max(count, max_count);
      }
      fprintf(stderr, "max count is %i\n", max_count);
      ////////////////////////////////////

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->type == task_type_pair ||
            tasks[tid]->type == task_type_sub_pair) {
          struct cell *ci_temp = ci_list[tid];
          //		     struct cell *cj_temp = cj;
          double shift[3];
          tasks[tid]->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        } else {
          tasks[tid]->sid = -1;
        }
      }
#endif
//		if(t->type != task_type_self && t->type !=
//task_type_sub_self)printf("type is not self or sub-self\n"); 		if(t->subtype !=
//task_subtype_density && t->subtype != task_subtype_gradient && t->subtype !=
//task_subtype_force)printf("subtype is not force, gradient or density\n");
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
       * i.e., rather than just one r->t, r would have a nested bunch of structs
       * r->t[size=nStreams].*/
      r->t = t;
#endif
      //////////////////////////////////////////////
      double t_reorder_start = clock();
      /*****Initialise GPU copies of variables and malloc pinned memory*****/
      struct cell_gpu **ci_gpu_list;
      struct part_gpu **parts_gpu_list;
      struct part_gpu *CSR_parts;
      int *task_first_part, *task_last_part;
      gpuErrchk(cudaMallocHost((void **)&ci_gpu_list,
                               count_tasks * sizeof(struct cell_gpu *)));
      gpuErrchk(cudaMallocHost((void **)&parts_gpu_list,
                               count_tasks * sizeof(struct part_gpu *)));
      gpuErrchk(
          cudaMallocHost((void **)&CSR_parts,
                         count_tasks * max_count * sizeof(struct part_gpu)));
      gpuErrchk(
          cudaMallocHost((void **)&task_first_part, count_tasks * sizeof(int)));
      gpuErrchk(
          cudaMallocHost((void **)&task_last_part, count_tasks * sizeof(int)));
      int total_num_parts = 0;
      int total_parts = 0;
      int first_part_tmp = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_list[tid]->hydro.count;
        task_first_part[tid] = first_part_tmp;
        task_last_part[tid] = first_part_tmp + count;
        first_part_tmp += count;
        gpuErrchk(cudaMallocHost(
            (void **)&ci_gpu_list[tid],
            sizeof(struct cell_gpu))); // Pinned allocation on host
        //			ci_gpu_list[tid]=malloc(sizeof(struct
        //cell_gpu));
        gpuErrchk(cudaMallocHost(
            (void **)&parts_gpu_list[tid],
            max_count * sizeof(struct part_gpu))); // Pinned allocation on host
        //			parts_gpu_list[tid]= malloc(count*sizeof(struct
        //part_gpu));

        //		struct part_gpu* parts_gpuj = (struct
        //part_gpu*)malloc(count*sizeof(struct part_gpu));
        ///*****Copy variables for cell i (self interaction)*****/
        ci_gpu_list[tid]->hydro.count = count;
        ci_gpu_list[tid]->hydro.h_max = ci_list[tid]->hydro.h_max;
        //			printf("count is %i for cell %i\n",
        //ci_gpu_list[tid]->hydro.count, tid);
        for (int d = 0; d < 3; d++) {
          ci_gpu_list[tid]->width[d] = ci_list[tid]->width[d];
          ci_gpu_list[tid]->loc[d] = ci_list[tid]->loc[d];
        }
        for (int p = 0; p < count; p++) {
          total_num_parts++;
          parts_gpu_list[tid][p].id = ci_list[tid]->hydro.parts[p].id;
          //		   	printf("id is %i for part %i\n",
          //parts_gpu_list[tid][p].id, p);
          for (int d = 0; d < 3; d++) {
            parts_gpu_list[tid][p].x[d] = ci_list[tid]->hydro.parts[p].x[d];
            parts_gpu_list[tid][p].v[d] = ci_list[tid]->hydro.parts[p].v[d];
            parts_gpu_list[tid][p].a_hydro[d] =
                ci_list[tid]->hydro.parts[p].a_hydro[d];
          }
          //		      printf("x is %f, y is %f, z is %f for part %i\n",
          //parts_gpu_list[tid][p].x[0], parts_gpu_list[tid][p].x[1],
          //parts_gpu_list[tid][p].x[2], p);
          parts_gpu_list[tid][p].mass = ci_list[tid]->hydro.parts[p].mass;
          parts_gpu_list[tid][p].h = ci_list[tid]->hydro.parts[p].h;
          parts_gpu_list[tid][p].u = ci_list[tid]->hydro.parts[p].u;
          parts_gpu_list[tid][p].u_dt = ci_list[tid]->hydro.parts[p].u_dt;
          parts_gpu_list[tid][p].rho = ci_list[tid]->hydro.parts[p].rho;
          parts_gpu_list[tid][p].div_v =
              ci_list[tid]->hydro.parts[p].viscosity.div_v;
          parts_gpu_list[tid][p].div_v_previous_step =
              ci_list[tid]->hydro.parts[p].viscosity.div_v_previous_step;
          parts_gpu_list[tid][p].alpha_visc =
              ci_list[tid]->hydro.parts[p].viscosity.alpha;
          parts_gpu_list[tid][p].v_sig =
              ci_list[tid]->hydro.parts[p].viscosity.v_sig;
          parts_gpu_list[tid][p].laplace_u =
              ci_list[tid]->hydro.parts[p].diffusion.laplace_u;
          parts_gpu_list[tid][p].alpha_diff =
              ci_list[tid]->hydro.parts[p].diffusion.alpha;
          parts_gpu_list[tid][p].f = ci_list[tid]->hydro.parts[p].force.f;
          parts_gpu_list[tid][p].soundspeed =
              ci_list[tid]->hydro.parts[p].force.soundspeed;
          parts_gpu_list[tid][p].h_dt = ci_list[tid]->hydro.parts[p].force.h_dt;
          parts_gpu_list[tid][p].balsara =
              ci_list[tid]->hydro.parts[p].force.balsara;
          parts_gpu_list[tid][p].pressure =
              ci_list[tid]->hydro.parts[p].force.pressure;
          parts_gpu_list[tid][p].time_bin =
              ci_list[tid]->hydro.parts[p].time_bin;
          parts_gpu_list[tid][p].wakeup =
              ci_list[tid]->hydro.parts[p].limiter_data.wakeup;
          parts_gpu_list[tid][p].min_ngb_time_bin =
              ci_list[tid]->hydro.parts[p].limiter_data.min_ngb_time_bin;
          parts_gpu_list[tid][p].to_be_synchronized =
              ci_list[tid]->hydro.parts[p].limiter_data.to_be_synchronized;
          parts_gpu_list[tid][p].wcount =
              ci_list[tid]->hydro.parts[p].density.wcount;
          parts_gpu_list[tid][p].wcount_dh =
              ci_list[tid]->hydro.parts[p].density.wcount_dh;
          parts_gpu_list[tid][p].rho_dh =
              ci_list[tid]->hydro.parts[p].density.rho_dh;
          parts_gpu_list[tid][p].div_v =
              ci_list[tid]->hydro.parts[p].viscosity.div_v;
          parts_gpu_list[tid][p].rot_v[0] =
              ci_list[tid]->hydro.parts[p].density.rot_v[0];
          parts_gpu_list[tid][p].rot_v[1] =
              ci_list[tid]->hydro.parts[p].density.rot_v[1];
          parts_gpu_list[tid][p].rot_v[2] =
              ci_list[tid]->hydro.parts[p].density.rot_v[2];
          parts_gpu_list[tid][p].SPH_sum = 0.f;
          CSR_parts[total_parts] = parts_gpu_list[tid][p];
          total_parts++;
        }
      }

      double t_reorder_end = clock();
      fprintf(stderr, "time to re-order is %f\n",
              (t_reorder_end - t_reorder_start) / CLOCKS_PER_SEC);
      //////////////////////////////////////////////

      //////Prepare streams and malloc memory on the GPU
      struct part_gpu *parts_gpu_all;
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      // create events and streams
      cudaProfilerStart();
      double t_stream_create_start = clock();
      cudaEvent_t startEvent, stopEvent, dummyEvent;
      cudaStream_t stream[count_tasks];
      gpuErrchk(cudaEventCreate(&startEvent));
      gpuErrchk(cudaEventCreate(&stopEvent));
      gpuErrchk(cudaEventCreate(&dummyEvent));
      for (int i = 0; i < count_tasks; ++i)
        gpuErrchk(cudaStreamCreate(&stream[i]));
      //			gpuErrchk( cudaStreamCreateWithFlags(&stream[i],
      //cudaStreamNonBlocking) ); 		launch_cuda_print_streams(10, stream[0], 0);
      double t_stream_create_end = clock();
      fprintf(stderr, "time to create streams is %f\n",
              (t_stream_create_end - t_stream_create_start) / CLOCKS_PER_SEC);

      double t_gpu_malloc_begin = clock();
      //		struct cell_gpu *d_ci_gpu[count_tasks];
      //		struct part_gpu *d_parts[count_tasks];
      struct part_gpu *d_parts;
      //		struct part_gpu **d_parts;
      //		gpuErrchk(cudaMalloc((void**)&d_parts,
      //count_tasks*sizeof(struct part_gpu*)));
      //		gpuErrchk(cudaMalloc((void**)&d_ci_gpu, count_tasks));
      FILE *fpcpu, *fpgpu;
      fpcpu = fopen("./res_CPU_density.txt", "w");
      double t_cpu_start = clock();
      double t_gpu_elapsed = 0.f;

      cudaMalloc((void **)&d_parts,
                 count_tasks * max_count * sizeof(struct part_gpu));
      //	   for(int tid=0; tid<count_tasks; tid++){
      //		    const int count=ci_gpu_list[tid]->hydro.count;
      //			const int size_parts=max_count*sizeof(struct
      //part_gpu);
      ////			const int size_cell=sizeof(struct cell_gpu);
      //			cudaMalloc((void**)&d_parts[tid], size_parts);
      ////			cudaMalloc((void**)&d_ci_gpu[tid], size_cell);
      //		}
      double t_gpu_malloc_end = clock();
      fprintf(stderr, "malloc time is %f\n",
              (t_gpu_malloc_end - t_gpu_malloc_begin) / CLOCKS_PER_SEC);
      //////////////////////////////////////////////////////
      /*Copy data to GPU and launch kernels*/
      double t_gpu_start = clock();
      cudaError_t err;
      //		cudaDeviceSynchronize();
      double t_memcpy_kernel = clock();
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_gpu_list[tid]->hydro.count;
        const int size_parts = count * sizeof(struct part_gpu);
        const int size_cell = sizeof(struct cell_gpu);
        int first_part = task_first_part[tid];
        int last_part = task_last_part[tid];
        //			cudaMemcpyAsync(d_parts[tid], parts_gpu_list[tid],
        //size_parts, cudaMemcpyHostToDevice, stream[tid]);
        cudaMemcpyAsync(&d_parts[first_part], &CSR_parts[first_part],
                        size_parts, cudaMemcpyHostToDevice, stream[tid]);
        //			cudaMemcpy(&d_parts[first_part],
        //&CSR_parts[first_part], size_parts, cudaMemcpyHostToDevice);

        //			cudaDeviceSynchronize();
        //			err = cudaPeekAtLastError();//cudaGetLastError();
        //// Get error code 			if ( err != cudaSuccess )
        //			 {
        //			   fprintf(stderr, "CUDA Error first memcpy: %s tid
        //is %i \n ", cudaGetErrorString(err), tid);
        //			 }
        float cellx = ci_gpu_list[tid]->loc[0],
              celly = ci_gpu_list[tid]->loc[1],
              cellz = ci_gpu_list[tid]->loc[2];
        const int numBlocks = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
        const char *loop_type = "density";
        launch_cuda_kernel_streams(d_parts, numBlocks, d_a, d_H, loop_type,
                                   stream[tid], tid, count, max_count, cellx,
                                   celly, cellz, first_part, last_part);
        //			cudaDeviceSynchronize();
        //			err = cudaPeekAtLastError();//cudaGetLastError();
        //// Get error code 			if ( err != cudaSuccess )
        //			 {
        //			   fprintf(stderr, "CUDA Error kernel: %s tid is %i
        //\n ", cudaGetErrorString(err), tid);
        //			 }
        //			cudaMemcpyAsync(parts_gpu_list[tid], d_parts[tid],
        //size_parts, cudaMemcpyDeviceToHost, stream[tid]);
        cudaMemcpyAsync(&CSR_parts[first_part], &d_parts[first_part],
                        size_parts, cudaMemcpyDeviceToHost, stream[tid]);
        //			cudaMemcpy(&CSR_parts[first_part],
        //&d_parts[first_part], size_parts, cudaMemcpyDeviceToHost);
        //			cudaDeviceSynchronize();
        //			err = cudaPeekAtLastError();//cudaGetLastError();
        //// Get error code 			if ( err != cudaSuccess )
        //			 {
        //			   fprintf(stderr, "CUDA Error second memcpy: %s tid
        //is %i \n ", cudaGetErrorString(err), tid);
        //			 }
      }
      cudaDeviceSynchronize();
      t_memcpy_kernel = (clock() - t_memcpy_kernel) / CLOCKS_PER_SEC;
      fprintf(stderr, "time for kernels and memcpys is %f\n", t_memcpy_kernel);
      double t_cpu_elapsed = 0.f;
      for (int tid = 0; tid < count_tasks; tid++) {
        //      	switch (tasks[tid]->type) {
        //		     case task_type_self:
        if (tasks[tid]->subtype == task_subtype_density) {
          double t1 = clock();
          runner_doself1_branch_density(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        //				 if (tasks[tid]->subtype ==
        //task_subtype_gradient){ 				 	 double t1=clock();
        //				 	 runner_doself1_branch_gradient(r,
        //ci_list[tid]); 				 	 double t2=clock(); 				 	 t_cpu_elapsed+=t2-t1;
        //				 }
        //				 if (tasks[tid]->subtype ==
        //task_subtype_force){ 				 	 double t1=clock(); 				 	 runner_doself2_branch_force(r,
        //ci_list[tid]); 				 	 double t2=clock(); 				 	 t_cpu_elapsed+=t2-t1;
        //				 }
        //      	}
      }

      double t_print_start = clock();
      fclose(fpcpu);
      fpcpu = fopen("./res_CPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_gpu_list[tid]->hydro.count;
        for (int p = 0; p < count; p++) {
          float xx = ci_list[tid]->hydro.parts[p].x[0],
                yy = ci_list[tid]->hydro.parts[p].x[1],
                zz = ci_list[tid]->hydro.parts[p].x[2];
          fprintf(fpcpu, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
                  ci_list[tid]->hydro.parts[p].rho,
                  ci_list[tid]->hydro.parts[p].density.rho_dh,
                  ci_list[tid]->hydro.parts[p].density.wcount,
                  ci_list[tid]->hydro.parts[p].density.wcount_dh,
                  ci_list[tid]->hydro.parts[p].viscosity.div_v,
                  ci_list[tid]->hydro.parts[p].density.rot_v[0]);
        }
      }

      fclose(fpcpu);
      double t_print_end = clock();
      double t_stream_destroy = clock();
      for (int tid = 0; tid < count_tasks; ++tid)
        gpuErrchk(cudaStreamDestroy(stream[tid]));
      t_stream_destroy = (clock() - t_stream_destroy) / CLOCKS_PER_SEC;
      fprintf(stderr, "stream destruction takes %f\n", t_stream_destroy);
      double t_gpu_end = clock();
      printf("time for cpu is %f\n", t_cpu_elapsed / CLOCKS_PER_SEC);
      printf("time for gpu is %f\n", (t_gpu_end - t_gpu_start - t_cpu_elapsed -
                                      t_print_end + t_print_start) /
                                         CLOCKS_PER_SEC);
      printf("finished looping through tasks\n");
      FILE *fp;
      fp = fopen("./res_GPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_gpu_list[tid]->hydro.count;
        //			for(int p=0; p<count; p++){
        //				ci_list[tid]->hydro.parts[p].density.wcount=parts_gpu_list[tid][p].wcount;
        //				ci_list[tid]->hydro.parts[p].density.wcount_dh=parts_gpu_list[tid][p].wcount_dh;
        //				ci_list[tid]->hydro.parts[p].density.rho_dh=parts_gpu_list[tid][p].rho_dh;
        //				ci_list[tid]->hydro.parts[p].rho=parts_gpu_list[tid][p].rho;
        //				ci_list[tid]->hydro.parts[p].viscosity.div_v=parts_gpu_list[tid][p].div_v;
        //				ci_list[tid]->hydro.parts[p].density.rot_v[0]=parts_gpu_list[tid][p].rot_v[0];
        //				ci_list[tid]->hydro.parts[p].density.rot_v[1]=parts_gpu_list[tid][p].rot_v[1];
        //				ci_list[tid]->hydro.parts[p].density.rot_v[2]=parts_gpu_list[tid][p].rot_v[2];
        //		  	}
        double t6 = clock();
        // printf("time for malloc %f\n",(t2-t1)/ CLOCKS_PER_SEC);
        // printf("time for memcpy H2D %f\n",(t3-t2)/ CLOCKS_PER_SEC);
        // printf("time for kernel %f\n",(t4-t3)/ CLOCKS_PER_SEC);
        // printf("time for memcpy D2H %f\n",(t5-t4)/ CLOCKS_PER_SEC);
        // printf("Total GPU time is %f\n\n",(t6-t1)/ CLOCKS_PER_SEC);
        //             GPU_runner_doself1_branch_gradient(ci_gpu, parts_gpu);
        for (int p = 0; p < count; p++) {
          parts_gpu_list[tid][p] = CSR_parts[p + task_first_part[tid]];
          float xx = parts_gpu_list[tid][p].x[0],
                yy = parts_gpu_list[tid][p].x[1],
                zz = parts_gpu_list[tid][p].x[2];
          fprintf(
              fp, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
              parts_gpu_list[tid][p].rho, parts_gpu_list[tid][p].rho_dh,
              parts_gpu_list[tid][p].wcount, parts_gpu_list[tid][p].wcount_dh,
              parts_gpu_list[tid][p].div_v, parts_gpu_list[tid][p].rot_v[0]);
        }
      }
      fclose(fp);
      for (int tid = 0; tid < count_tasks; tid++) {
        gpuErrchk(cudaFree(d_parts));
      }
      cudaFreeHost(ci_gpu_list);
      cudaFreeHost(parts_gpu_list);
      cudaProfilerStop();
      cudaDeviceReset();
      exit(0);
      ///* Mark that we have run this task on these cells */
      //#ifdef SWIFT_DEBUG_CHECKS
      //				if (ci_list[tid] != NULL) {
      //				  ci_list[tid]->tasks_executed[tasks[tid]->type]++;
      //				  ci_list[tid]->subtasks_executed[tasks-[tid]>subtype]++;
      //				}
      ////				if (cj != NULL) {
      ////				  cj->tasks_executed[t->type]++;
      ////				  cj->subtasks_executed[t->subtype]++;
      ////				}

      //				/* This runner is not doing a task
      //anymore */ 				r->t = NULL; #endif

      /* We're done with this task, see if we get a next one. */
      //      prev = t;
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////****same as runner_main_cuda but bunching tasks into bundles. Each stream
///works on one bundle***********/////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void *runner_main_cuda_bundles(void *data) {
  cudaDeviceReset();
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  cudaError_t cu_error;
  cudaMemPool_t memPool;
  gpuErrchk(cudaDeviceGetDefaultMemPool(&memPool, devId));
  int maxmem = 0.9f * prop.totalGlobalMem;
  //  gpuErrchk(cudaMemPoolSetAttribute(memPool,
  //  cudaMemPoolAttrReleaseThreshold, (void *)&maxmem));
  int counter = 0;

  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    struct task **tasks;
    gpuErrchk(cudaMallocHost(
        (void **)&tasks,
        n_streams * sizeof(struct task))); // Pinned allocation on host
    for (int i = 0; i < n_streams; i++) {
      tasks[i] = NULL;
    }

    while (1) {
      for (int i = 0; i < n_streams; i++) {
        tasks[i] = NULL;
      }
      // for(int tid=0; tid<n_streams; tid++){
      int count_tasks = 0;
      int count_stale = 0;
      ////////////////////////////////////
      /*Grab a bunch of tasks*/
      while (1) { //(count_tasks<n_streams){
        // there's a potential bug here. Code hangs if n_streams set to > 256 in
        // cuda_headers.h
        /* If there's no old task, try to get a new one. */
        if (tasks[count_tasks] == NULL) {
          /* Get the task. */
          TIMER_TIC
          tasks[count_tasks] = scheduler_gettask(sched, r->qid, prev);
          TIMER_TOC(timer_gettask);
          if (tasks[count_tasks] != NULL) {
            count_tasks++;
          }
        }
        /* If there's an old task, move onto the next one */
        // else count_tasks++;
        count_stale++;
        if (count_stale >= n_streams) {
          break;
        }
      }
      fprintf(stderr, "Count tasks is %i\n", count_tasks);
      ////////////////////////////////////

      ////////////////////////////////////
      /*Get cell data for each task*/
      for (int tid = 1; tid < count_tasks; tid++) {
        struct cell *c_gpu = tasks[tid]->ci;
      }

      /* Get the cells. */
      struct cell **ci_list;
      gpuErrchk(cudaMallocHost(
          (void **)&ci_list,
          count_tasks * sizeof(struct cell))); // Pinned allocation on host
      for (int tid = 0; tid < count_tasks; tid++) {
        ci_list[tid] = tasks[tid]->ci;
        fprintf(stderr, "positon is %d\n", ci_list[tid]->hydro.parts[0].x[0]);
      }
      ////////////////////////////////////
      exit(0);

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->type == task_type_pair ||
            tasks[tid]->type == task_type_sub_pair) {
          struct cell *ci_temp = ci_list[tid];
          //		     struct cell *cj_temp = cj;
          double shift[3];
          tasks[tid]->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        } else {
          tasks[tid]->sid = -1;
        }
      }
#endif
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
      i.e., rather than just one r->t, r would have a nested bunch of structs
       r->t[size=nStreams].*/
      r->t = t;
#endif

      //////////////////////////////////////////////
      /*****Initialise GPU suitable copies of variables and malloc pinned
       * memory*****/
      struct cell_gpu_flat ci_gpu_list[count_tasks];
      struct part_gpu *parts_gpu_list[count_tasks];

      int total_num_parts = 0;
      ///////Start GPU timer
      //		double t_gpu_start=clock();
      ////////////////////////////////////////////

      ////////Malloc on HOST
      int count_parts_total = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_list[tid]->hydro.count;
        count_parts_total += count;
        gpuErrchk(cudaMallocHost(
            (void **)&parts_gpu_list[tid],
            count * sizeof(struct part_gpu))); // Pinned allocation on host
        gpuErrchk(cudaMallocHost(
            (void **)&ci_gpu_list[tid],
            sizeof(struct cell_gpu_flat))); // Pinned allocation on host

        ////////////////////////////////////////////
        ///*****Copy variables for cell i (self interaction)*****/
        ci_gpu_list[tid].count = count;
        ci_gpu_list[tid].h_max = ci_list[tid]->hydro.h_max;

        ci_gpu_list[tid].width[0] = ci_list[tid]->width[0];
        ci_gpu_list[tid].width[1] = ci_list[tid]->width[1];
        ci_gpu_list[tid].width[2] = ci_list[tid]->width[2];
        ci_gpu_list[tid].loc[0] = ci_list[tid]->loc[0];
        ci_gpu_list[tid].loc[1] = ci_list[tid]->loc[1];
        ci_gpu_list[tid].loc[2] = ci_list[tid]->loc[2];

        for (int p = 0; p < count; p++) {
          total_num_parts++;
          parts_gpu_list[tid][p].id = ci_list[tid]->hydro.parts[p].id;

          parts_gpu_list[tid][p].count = count;
          parts_gpu_list[tid][p].h_max = ci_list[tid]->hydro.h_max;

          for (int d = 0; d < 3; d++) {
            parts_gpu_list[tid][p].x[d] = ci_list[tid]->hydro.parts[p].x[d];
            parts_gpu_list[tid][p].v[d] = ci_list[tid]->hydro.parts[p].v[d];
            parts_gpu_list[tid][p].a_hydro[d] =
                ci_list[tid]->hydro.parts[p].a_hydro[d];
            parts_gpu_list[tid][p].loc[d] = ci_list[tid]->loc[d];
          }

          parts_gpu_list[tid][p].mass = ci_list[tid]->hydro.parts[p].mass;
          parts_gpu_list[tid][p].h = ci_list[tid]->hydro.parts[p].h;
          parts_gpu_list[tid][p].u = ci_list[tid]->hydro.parts[p].u;
          parts_gpu_list[tid][p].u_dt = ci_list[tid]->hydro.parts[p].u_dt;
          parts_gpu_list[tid][p].rho = ci_list[tid]->hydro.parts[p].rho;
          parts_gpu_list[tid][p].div_v =
              ci_list[tid]->hydro.parts[p].viscosity.div_v;
          parts_gpu_list[tid][p].div_v_previous_step =
              ci_list[tid]->hydro.parts[p].viscosity.div_v_previous_step;
          parts_gpu_list[tid][p].alpha_visc =
              ci_list[tid]->hydro.parts[p].viscosity.alpha;
          parts_gpu_list[tid][p].v_sig =
              ci_list[tid]->hydro.parts[p].viscosity.v_sig;
          parts_gpu_list[tid][p].laplace_u =
              ci_list[tid]->hydro.parts[p].diffusion.laplace_u;
          parts_gpu_list[tid][p].alpha_diff =
              ci_list[tid]->hydro.parts[p].diffusion.alpha;
          parts_gpu_list[tid][p].f = ci_list[tid]->hydro.parts[p].force.f;
          parts_gpu_list[tid][p].soundspeed =
              ci_list[tid]->hydro.parts[p].force.soundspeed;
          parts_gpu_list[tid][p].h_dt = ci_list[tid]->hydro.parts[p].force.h_dt;
          parts_gpu_list[tid][p].balsara =
              ci_list[tid]->hydro.parts[p].force.balsara;
          parts_gpu_list[tid][p].pressure =
              ci_list[tid]->hydro.parts[p].force.pressure;
          parts_gpu_list[tid][p].time_bin =
              ci_list[tid]->hydro.parts[p].time_bin;
          parts_gpu_list[tid][p].wakeup =
              ci_list[tid]->hydro.parts[p].limiter_data.wakeup;
          parts_gpu_list[tid][p].min_ngb_time_bin =
              ci_list[tid]->hydro.parts[p].limiter_data.min_ngb_time_bin;
          parts_gpu_list[tid][p].to_be_synchronized =
              ci_list[tid]->hydro.parts[p].limiter_data.to_be_synchronized;
          parts_gpu_list[tid][p].wcount =
              ci_list[tid]->hydro.parts[p].density.wcount;
          parts_gpu_list[tid][p].wcount_dh =
              ci_list[tid]->hydro.parts[p].density.wcount_dh;
          parts_gpu_list[tid][p].rho_dh =
              ci_list[tid]->hydro.parts[p].density.rho_dh;
          parts_gpu_list[tid][p].div_v =
              ci_list[tid]->hydro.parts[p].viscosity.div_v;
          parts_gpu_list[tid][p].rot_v[0] =
              ci_list[tid]->hydro.parts[p].density.rot_v[0];
          parts_gpu_list[tid][p].rot_v[1] =
              ci_list[tid]->hydro.parts[p].density.rot_v[1];
          parts_gpu_list[tid][p].rot_v[2] =
              ci_list[tid]->hydro.parts[p].density.rot_v[2];
          parts_gpu_list[tid][p].SPH_sum = 0.f;
        }
      }
      fprintf(stderr, "got here 1\n");
      //////////////////////////////////////////////

      //////Prepare streams and malloc memory on the GPU
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      // create events and streams
      ///////Start GPU timer
      double t_gpu_start = clock();
      ////////////////////////////////////////////

      ////////////////////////////////////////////////
      ///////////Malloc GPU data

      cudaProfilerStart();
      double t_gpu_malloc_begin = clock();

      FILE *fpcpu, *fpgpu;
      fpcpu = fopen("./res_CPU_density.txt", "w");
      double t_cpu_start = clock();
      double t_gpu_elapsed = 0.f;

      // see
      // https://stackoverflow.com/questions/23609770/cuda-double-pointer-memory-copy
      // for details of how this is allocated and copied
      struct part_gpu **d_all_parts;
      struct cell_gpu_flat *d_all_cells;
      gpuErrchk(cudaMalloc(
          (void **)&d_all_parts,
          count_tasks *
              sizeof(struct part_gpu *))); // Pinned allocation on host
      gpuErrchk(cudaMalloc(
          (void **)&d_all_cells,
          count_tasks *
              sizeof(struct cell_gpu_flat))); // Pinned allocation on host

      // For another possible solution, see
      // https://forums.developer.nvidia.com/t/how-do-i-pass-a-double-pointers-array-to-the-device-im-getting-cudaerrorillegaladdress/72518/2

      int size_parts_mmcpy = 0;

      struct part_gpu *h_all_parts[count_tasks];
      struct cell_gpu_flat h_all_cells[count_tasks];

      cudaEvent_t startEvent, stopEvent, dummyEvent;

      gpuErrchk(cudaEventCreate(&startEvent));
      gpuErrchk(cudaEventCreate(&stopEvent));
      gpuErrchk(cudaEventCreate(&dummyEvent));

      int max_parts = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        int count = ci_gpu_list[tid].count;
        max_parts = max(max_parts, count);
        gpuErrchk(cudaMalloc((void **)&h_all_parts[tid],
                             count * sizeof(struct part_gpu)));
        gpuErrchk(cudaMalloc((void **)&h_all_cells[tid],
                             sizeof(struct cell_gpu_flat)));
        // make host pointers point to parts_gpu_list and ci_gpu_list
        h_all_parts[tid] = parts_gpu_list[tid];
        h_all_cells[tid] = ci_gpu_list[tid];
      }
      double t_gpu_malloc_finish = clock();
      fprintf(stderr, "malloc time is %f\n",
              (t_gpu_malloc_finish - t_gpu_malloc_begin) / CLOCKS_PER_SEC);

      int bundle_size = 32;
      int nBundles = (count_tasks + bundle_size - 1) / bundle_size;
      cudaStream_t stream[nBundles];
      int tasksperbundle = (count_tasks + nBundles - 1) / nBundles;
      double t_gpu_copyandkernel_begin = clock();
      for (int i = 0; i < nBundles; ++i)
        //			gpuErrchk( cudaStreamCreate(&stream[i]) );
        gpuErrchk(cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking));
      t_gpu_elapsed = 0.f;
      for (int bid = 0; bid < nBundles; bid++) {
        for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
             tid++) {
          max_parts = 0;
          if (tid < count_tasks) {
            int count = ci_gpu_list[tid].count;
            max_parts = max(max_parts, count);
          }
        }
        int tid = 0;
        int offset = bid * tasksperbundle;
        int sizecells = 0, sizeparts = 0;

        int numBlocks_y = tasksperbundle;
        sizecells = tasksperbundle * sizeof(struct cell_gpu_flat);
        sizeparts = tasksperbundle * sizeof(struct part_gpu *);
        int tasks_left = tasksperbundle;
        if (bid == nBundles - 1) {
          tasks_left = count_tasks - (nBundles - 1) * tasksperbundle;
          sizecells = tasks_left * sizeof(struct cell_gpu_flat);
          sizeparts = tasks_left * sizeof(struct part_gpu *);
          numBlocks_y = tasks_left;
        }

        //			gpuErrchk(cudaMemcpyAsync(&d_all_cells[offset],
        //&h_all_cells[offset], sizecells, cudaMemcpyHostToDevice,
        //stream[bid]));
        gpuErrchk(cudaMemcpyAsync(&d_all_parts[offset], &h_all_parts[offset],
                                  sizeparts, cudaMemcpyHostToDevice,
                                  stream[bid]));

        int block_size = 128;
        int numBlocks; //=(tasks_left + block_size - 1) / block_size;
                       ////////////////////////////////////////////////
        //		launch_cuda_kernel_bundles_test(d_all_cells, d_all_parts,
        //numBlocks, d_a, d_H, count_tasks);
        int numBlocks_x = (max_parts + block_size - 1) / block_size;
        numBlocks = tasks_left * numBlocks_x;
        fprintf(stderr,
                "bid is %i, max parts/cell in bundle %i Ntasks is %i "
                "tasksperbundle is %i\n NBlocks x is %i NBlocksy is %i "
                "numBundles is %i\n",
                bid, max_parts, count_tasks, tasksperbundle, numBlocks_x,
                numBlocks_y, nBundles);
        const char *loop_type = "density";
        launch_cuda_kernel_bundles(d_all_cells, d_all_parts, numBlocks, d_a,
                                   d_H, loop_type, stream[bid], bid, block_size,
                                   count_tasks, tasksperbundle, numBlocks_x,
                                   numBlocks_y, tid, offset);
        //////////////////////////////////////////////////////

        //			gpuErrchk(cudaMemcpyAsync(&h_all_cells[offset],
        //&d_all_cells[offset], sizecells, cudaMemcpyDeviceToHost,
        //stream[bid]));
        gpuErrchk(cudaMemcpyAsync(&h_all_parts[offset], &d_all_parts[offset],
                                  sizeparts, cudaMemcpyDeviceToHost,
                                  stream[bid]));
      }
      cudaDeviceSynchronize();
      double t_gpu_copyandkernel_finish = clock();
      fprintf(stderr, "kernel and cpy time is %f\n",
              (t_gpu_copyandkernel_finish - t_gpu_copyandkernel_begin) /
                  CLOCKS_PER_SEC);
      //		cu_error = cudaPeekAtLastError();//cudaGetLastError();
      //// Get error code 		if ( cu_error != cudaSuccess )
      //		{
      //			fprintf(stderr, "CUDA Error: %s\n",
      //cudaGetErrorString(cu_error));
      //		}

      for (int tid = 0; tid < count_tasks; tid++) {
        ci_gpu_list[tid] = h_all_cells[tid];
        parts_gpu_list[tid] = h_all_parts[tid];
      }

      double t_cpu_elapsed = 0.f;
      for (int tid = 0; tid < count_tasks; tid++) {
        //      	switch (tasks[tid]->type) {
        //		     case task_type_self:
        if (tasks[tid]->subtype == task_subtype_density) {
          double t1 = clock();
          runner_doself1_branch_density(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        if (tasks[tid]->subtype == task_subtype_gradient) {
          double t1 = clock();
          runner_doself1_branch_gradient(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        if (tasks[tid]->subtype == task_subtype_force) {
          double t1 = clock();
          runner_doself2_branch_force(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        //      	}
      }
      double t_cpu_print_start = clock();
      fclose(fpcpu);
      fpcpu = fopen("./res_CPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_gpu_list[tid].count;
        for (int p = 0; p < count; p++) {
          float xx = ci_list[tid]->hydro.parts[p].x[0],
                yy = ci_list[tid]->hydro.parts[p].x[1],
                zz = ci_list[tid]->hydro.parts[p].x[2];
          fprintf(fpcpu, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
                  ci_list[tid]->hydro.parts[p].rho,
                  ci_list[tid]->hydro.parts[p].density.rho_dh,
                  ci_list[tid]->hydro.parts[p].density.wcount,
                  ci_list[tid]->hydro.parts[p].density.wcount_dh,
                  ci_list[tid]->hydro.parts[p].viscosity.div_v,
                  ci_list[tid]->hydro.parts[p].density.rot_v[0]);
        }
      }

      fclose(fpcpu);
      double t_cpu_print_end = clock();
      cudaProfilerStop();
      //		for (int tid = 0; tid < count_tasks; ++tid)
      //			gpuErrchk( cudaStreamDestroy(stream[tid]));
      double t_gpu_end = clock();
      printf("time for cpu is %f\n", t_cpu_elapsed / CLOCKS_PER_SEC);
      printf("time for gpu is %f\n", (t_gpu_end - t_gpu_start - t_cpu_elapsed -
                                      t_cpu_print_end + t_cpu_print_start) /
                                         CLOCKS_PER_SEC);
      printf("finished looping through tasks\n");
      FILE *fp;
      fp = fopen("./res_GPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_gpu_list[tid].count;
        //			for(int p=0; p<count; p++){
        //				ci_list[tid]->hydro.parts[p].density.wcount=parts_gpu_list[tid][p].wcount;
        //				ci_list[tid]->hydro.parts[p].density.wcount_dh=parts_gpu_list[tid][p].wcount_dh;
        //				ci_list[tid]->hydro.parts[p].density.rho_dh=parts_gpu_list[tid][p].rho_dh;
        //				ci_list[tid]->hydro.parts[p].rho=parts_gpu_list[tid][p].rho;
        //				ci_list[tid]->hydro.parts[p].viscosity.div_v=parts_gpu_list[tid][p].div_v;
        //				ci_list[tid]->hydro.parts[p].density.rot_v[0]=parts_gpu_list[tid][p].rot_v[0];
        //				ci_list[tid]->hydro.parts[p].density.rot_v[1]=parts_gpu_list[tid][p].rot_v[1];
        //				ci_list[tid]->hydro.parts[p].density.rot_v[2]=parts_gpu_list[tid][p].rot_v[2];
        //		  	}
        double t6 = clock();
        // printf("time for malloc %f\n",(t2-t1)/ CLOCKS_PER_SEC);
        // printf("time for memcpy H2D %f\n",(t3-t2)/ CLOCKS_PER_SEC);
        // printf("time for kernel %f\n",(t4-t3)/ CLOCKS_PER_SEC);
        // printf("time for memcpy D2H %f\n",(t5-t4)/ CLOCKS_PER_SEC);
        // printf("Total GPU time is %f\n\n",(t6-t1)/ CLOCKS_PER_SEC);
        //             GPU_runner_doself1_branch_gradient(ci_gpu, parts_gpu);
        for (int p = 0; p < count; p++) {
          float xx = parts_gpu_list[tid][p].x[0],
                yy = parts_gpu_list[tid][p].x[1],
                zz = parts_gpu_list[tid][p].x[2];
          fprintf(
              fp, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
              parts_gpu_list[tid][p].rho, parts_gpu_list[tid][p].rho_dh,
              parts_gpu_list[tid][p].wcount, parts_gpu_list[tid][p].wcount_dh,
              parts_gpu_list[tid][p].div_v, parts_gpu_list[tid][p].rot_v[0]);
        }
      }
      fclose(fp);
      //      cudaFreeHost(ci_gpu_list);
      cudaFreeHost(parts_gpu_list);
      exit(0);
      ///* Mark that we have run this task on these cells */
      //#ifdef SWIFT_DEBUG_CHECKS
      //				if (ci_list[tid] != NULL) {
      //				  ci_list[tid]->tasks_executed[tasks[tid]->type]++;
      //				  ci_list[tid]->subtasks_executed[tasks-[tid]>subtype]++;
      //				}
      ////				if (cj != NULL) {
      ////				  cj->tasks_executed[t->type]++;
      ////				  cj->subtasks_executed[t->subtype]++;
      ////				}

      //				/* This runner is not doing a task
      //anymore */ 				r->t = NULL; #endif

      /* We're done with this task, see if we get a next one. */
      //      prev = t;
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////****same as runner_main_cuda but bunching tasks into bundles. Each stream
///works on one bundle***********/////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void *runner_main_cuda_bundles_revised(void *data) {
  cudaDeviceReset();
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  cudaError_t cu_error;
  cudaMemPool_t memPool;
  gpuErrchk(cudaDeviceGetDefaultMemPool(&memPool, devId));
  int maxmem = 0.9f * prop.totalGlobalMem;
  gpuErrchk(cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold,
                                    (void *)&maxmem));
  int counter = 0;

  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    //	 gpuErrchk( cudaMallocHost((void**)&tasks, n_streams*sizeof(struct
    //task)) ); // Pinned allocation on host Typecasting from void** to struct
    // task **//

    ////////////////////////////////////////////////////////////////////////////////
    ///////This is likely where the problem is!!!!
    ////////////////////////////////////////////////////////////////////////////////
    struct task *tasks[n_streams];

    while (1) {
      for (int i = 0; i < n_streams; i++) {
        tasks[i] = NULL;
      }
      // for(int tid=0; tid<n_streams; tid++){
      int count_tasks = 0;
      int count_stale = 0;
      ////////////////////////////////////
      /*Grab a bunch of tasks*/
      while (1) { //(count_tasks<n_streams){
        // there's a potential bug here. Code hangs if n_streams set to > 256 in
        // cuda_headers.h
        /* If there's no old task, try to get a new one. */
        if (tasks[count_tasks] == NULL) {
          /* Get the task. */
          TIMER_TIC
          tasks[count_tasks] = scheduler_gettask(sched, r->qid, prev);
          TIMER_TOC(timer_gettask);
          if (tasks[count_tasks] != NULL) {
            count_tasks++;
          }
        }
        /* If there's an old task, move onto the next one */
        // else count_tasks++;
        count_stale++;
        if (count_stale >= n_streams) {
          break;
        }
      }
      fprintf(stderr, "Count tasks is %i\n", count_tasks);
      ////////////////////////////////////

      ////////////////////////////////////
      /* Get the cells. */
      struct cell **ci_list =
          (struct cell **)malloc(count_tasks * sizeof(struct cell *));
      int count_all_parts = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        //			ci_list[tid]=(struct cell) malloc(sizeof(struct
        //cell*));
        ci_list[tid] = tasks[tid]->ci;
        struct cell *c_tmp = ci_list[tid];
        count_all_parts += ci_list[tid]->hydro.count;
      }
      ////////////////////////////////////

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->type == task_type_pair ||
            tasks[tid]->type == task_type_sub_pair) {
          struct cell *ci_temp = ci_list[tid];
          //		     struct cell *cj_temp = cj;
          double shift[3];
          tasks[tid]->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        } else {
          tasks[tid]->sid = -1;
        }
      }
#endif
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
      i.e., rather than just one r->t, r would have a nested bunch of structs
       r->t[size=nStreams].*/
      r->t = t;
#endif

      //////////////////////////////////////////////
      /*****Initialise GPU suitable copies of variables and malloc pinned
       * memory*****/
      struct part_gpu **parts_gpu_list =
          malloc(count_tasks * sizeof(struct part_gpu *));
      //		struct part_gpu **
      //parts_gpu_list=malloc(count_tasks*sizeof(struct part_gpu *)); 		struct
      //part_gpu **parts_gpu_list; 		gpuErrchk(
      //cudaMallocHost((void**)parts_gpu_list, count_tasks*sizeof(struct
      //part_gpu *)) ); // Pinned allocation on host
      struct part_gpu *parts_gpu_CSR;
      gpuErrchk(cudaMallocHost(
          (void **)&parts_gpu_CSR,
          count_all_parts *
              sizeof(struct part_gpu))); // Pinned allocation on host
      fprintf(stderr, "Size I'm allocating is %i\n",
              count_all_parts * sizeof(struct part_gpu));
      //		exit(0);
      int *task_first_part, *task_last_part;
      int *d_task_first_part, *d_task_last_part;
      gpuErrchk(cudaMallocHost((void **)&task_first_part,
                               count_tasks *
                                   sizeof(int))); // Pinned allocation on host
      gpuErrchk(cudaMallocHost((void **)&task_last_part,
                               count_tasks *
                                   sizeof(int))); // Pinned allocation on host

      int bundle_size = 1;
      int nBundles = (count_tasks + bundle_size - 1) / bundle_size;

      int *bundle_first_part, *bundle_last_part;
      int *d_bundle_first_part, *d_bundle_last_part;
      gpuErrchk(
          cudaMallocHost((void **)&bundle_first_part,
                         nBundles * sizeof(int))); // Pinned allocation on host
      gpuErrchk(
          cudaMallocHost((void **)&bundle_last_part,
                         nBundles * sizeof(int))); // Pinned allocation on host
      int total_num_parts = 0;
      ///////Start GPU timer
      //		double t_gpu_start=clock();
      ////////////////////////////////////////////

      ////////Malloc on HOST
      int first_part_tmp = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        int count = ci_list[tid]->hydro.count;
        parts_gpu_list[tid] = malloc(count * sizeof(struct part_gpu));
        parts_gpu_list[tid][0].id = ci_list[tid]->hydro.parts[0].id;
        populate_parts_list(ci_list[tid], parts_gpu_list[tid]);
        for (int p = 0; p < count; p++) {
          parts_gpu_CSR[p + first_part_tmp] = parts_gpu_list[tid][p];
          parts_gpu_CSR[p + first_part_tmp].tid = tid;
        }
        task_first_part[tid] = first_part_tmp;
        task_last_part[tid] = first_part_tmp + count;
        first_part_tmp += count;
        total_num_parts = count_all_parts;
        if (tid % bundle_size == 0) {
          int bid = tid / bundle_size;
          bundle_first_part[bid] = task_first_part[tid];
        }
      }
      for (int bid = 0; bid < nBundles - 1; bid++) {
        bundle_last_part[bid] = bundle_first_part[bid + 1];
      }
      bundle_last_part[nBundles - 1] =
          count_all_parts - bundle_last_part[nBundles - 1];

      fprintf(stderr, "GOT HERE\n");

      struct part_gpu *d_parts_gpu_CSR;
      cudaMalloc((void **)&d_parts_gpu_CSR,
                 count_all_parts * sizeof(struct part_gpu));
      cudaMalloc((void **)&d_task_first_part, count_tasks * sizeof(int));
      cudaMalloc((void **)&d_task_last_part, count_tasks * sizeof(int));
      cudaMalloc((void **)&d_bundle_first_part, nBundles * sizeof(int));
      cudaMalloc((void **)&d_bundle_last_part, nBundles * sizeof(int));
      //////////////////////////////////////////////

      //////Prepare streams and malloc memory on the GPU
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      // create events and streams
      ///////Start GPU timer
      double t_gpu_start = clock();
      ////////////////////////////////////////////

      ////////////////////////////////////////////////
      ///////////Malloc GPU data

      cudaProfilerStart();
      double t_gpu_malloc_begin = clock();

      FILE *fpcpu, *fpgpu;
      fpcpu = fopen("./res_CPU_density.txt", "w");
      double t_cpu_start = clock();
      double t_gpu_elapsed = 0.f;

      // see
      // https://stackoverflow.com/questions/23609770/cuda-double-pointer-memory-copy
      // for details of how this is allocated and copied

      fprintf(stderr, "Malloc d_all_parts\n");
      cudaError_t err;
      for (int tid = 0; tid < count_tasks; tid++) {
        int first_part_tmp = task_first_part[tid];
        int count = parts_gpu_list[tid][0].count;
        gpuErrchk(cudaMemcpy(
            &d_parts_gpu_CSR[first_part_tmp], &parts_gpu_CSR[first_part_tmp],
            count * sizeof(struct part_gpu), cudaMemcpyHostToDevice));

        err = cudaPeekAtLastError(); // cudaGetLastError();        // Get error
                                     // code
        if (err != cudaSuccess) {
          fprintf(stderr, "CUDA kernel Error first memcpy: %s \n ",
                  cudaGetErrorString(err));
        }
      }
      fprintf(stderr,
              "number of parts %i, first part in last task %i, last part in "
              "last task %i\n",
              count_all_parts, task_first_part[count_tasks - 1],
              task_last_part[count_tasks - 1]);
      fprintf(stderr,
              "number of parts %i, first part in fist bundle %i, last part in "
              "last bundle %i\n",
              count_all_parts, bundle_first_part[nBundles - 1],
              bundle_last_part[nBundles - 1]);
      //	   exit(0);
      cudaDeviceSynchronize();

      gpuErrchk(cudaMemcpy(d_task_first_part, task_first_part,
                           count_tasks * sizeof(int), cudaMemcpyHostToDevice));
      gpuErrchk(cudaMemcpy(d_task_last_part, task_last_part,
                           count_tasks * sizeof(int), cudaMemcpyHostToDevice));
      gpuErrchk(cudaMemcpy(d_bundle_first_part, bundle_first_part,
                           nBundles * sizeof(int), cudaMemcpyHostToDevice));
      gpuErrchk(cudaMemcpy(d_bundle_last_part, bundle_last_part,
                           nBundles * sizeof(int), cudaMemcpyHostToDevice));
      err =
          cudaPeekAtLastError(); // cudaGetLastError();        // Get error code
      if (err != cudaSuccess) {
        fprintf(stderr, "CUDA kernel Error first memcpy: %s \n ",
                cudaGetErrorString(err));
      }
      fprintf(stderr, "Malloc d_all_parts\n");

      // For another possible solution, see
      // https://forums.developer.nvidia.com/t/how-do-i-pass-a-double-pointers-array-to-the-device-im-getting-cudaerrorillegaladdress/72518/2

      int size_parts_mmcpy = 0;

      cudaEvent_t startEvent, stopEvent, dummyEvent;

      gpuErrchk(cudaEventCreate(&startEvent));
      gpuErrchk(cudaEventCreate(&stopEvent));
      gpuErrchk(cudaEventCreate(&dummyEvent));
      fprintf(stderr, "got here in the end\n");
      double t_gpu_malloc_finish = clock();
      fprintf(stderr, "malloc time is %f\n",
              (t_gpu_malloc_finish - t_gpu_malloc_begin) / CLOCKS_PER_SEC);

      cudaStream_t stream[nBundles];
      int tasksperbundle = (count_tasks + nBundles - 1) / nBundles;
      fprintf(stderr, "count_tasks/nBundles %i, tasksperbundle %i\n",
              count_tasks / nBundles, tasksperbundle);
      //		exit(0);
      double t_gpu_copyandkernel_begin = clock();
      for (int i = 0; i < nBundles; ++i)
        //			gpuErrchk( cudaStreamCreate(&stream[i]) );
        gpuErrchk(cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking));
      err =
          cudaPeekAtLastError(); // cudaGetLastError();        // Get error code
      if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Stream Creation Error: %s \n ",
                cudaGetErrorString(err));
      }
      t_gpu_elapsed = 0.f;
      for (int bid = 0; bid < nBundles; bid++) {
        int max_parts = 0;
        int parts_in_bundle = 0;
        for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
             tid++) {
          max_parts = 0;
          if (tid < count_tasks) {
            int count = task_last_part[tid] - task_first_part[tid];
            //					if(tid<count_tasks-1)
            parts_in_bundle += count;
            max_parts = max(max_parts, count);
          }
        }
        int tid = 0;
        int offset = bid * tasksperbundle;
        int sizecells = 0, sizeparts = 0;

        int numBlocks_y = tasksperbundle;
        sizeparts = tasksperbundle * sizeof(struct part_gpu *);
        int tasks_left = tasksperbundle;
        if (bid == nBundles - 1) {
          tasks_left = count_tasks - (nBundles - 1) * tasksperbundle;
          sizeparts = tasks_left * sizeof(struct part_gpu *);
          numBlocks_y = tasks_left;
        }
        int block_size = BLOCK_SIZE;
        int numBlocks_x = (max_parts + block_size - 1) / block_size;
        int numBlocks = (parts_in_bundle + block_size - 1) / block_size;
        //			fprintf(stderr, "numblocks %i, parts in bundle
        //%i\n", numBlocks, parts_in_bundle); 			fprintf(stderr, "bid is %i, max
        //parts/cell in bundle %i Ntasks is %i tasksperbundle is %i\n NBlocks x
        //is %i NBlocksy is %i numBundles is %i\n", bid, max_parts, count_tasks,
        //tasksperbundle, numBlocks_x, numBlocks_y, nBundles);
        const char *loop_type = "density";
        launch_cuda_kernel_bundles_revised(
            d_parts_gpu_CSR, d_task_first_part, d_task_last_part,
            d_bundle_first_part, d_bundle_last_part, numBlocks, d_a, d_H,
            loop_type, stream[bid], bid, block_size, count_tasks,
            tasksperbundle, numBlocks_x, numBlocks_y, tid, offset);
        //////////////////////////////////////////////////////
        err = cudaPeekAtLastError(); // cudaGetLastError();        // Get error
                                     // code
        if (err != cudaSuccess) {
          fprintf(stderr, "CUDA kernel Error: %s \n ", cudaGetErrorString(err));
        }
      }
      cudaDeviceSynchronize();
      //		for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=task_last_part[tid]-task_first_part[tid];
      //		   gpuErrchk(cudaMemcpy(&parts_gpu_CSR[first_part_tmp],
      //&d_parts_gpu_CSR[first_part_tmp], count*sizeof(struct part_gpu),
      //cudaMemcpyDeviceToHost)); 		   err =
      //cudaPeekAtLastError();//cudaGetLastError();        // Get error code 		   if
      //( err != cudaSuccess )
      //		     {
      //		   	   fprintf(stderr, "CUDA Error second memcpy: %s
      //tid is %i\n ", cudaGetErrorString(err), tid);
      //		   	 }
      //		}
      gpuErrchk(cudaMemcpy(parts_gpu_CSR, d_parts_gpu_CSR,
                           count_all_parts * sizeof(struct part_gpu),
                           cudaMemcpyDeviceToHost));
      err =
          cudaPeekAtLastError(); // cudaGetLastError();        // Get error code
      if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error second memcpy: %s \n ",
                cudaGetErrorString(err));
      }

      double t_gpu_copyandkernel_finish = clock();
      fprintf(stderr, "kernel and cpy time is %f\n",
              (t_gpu_copyandkernel_finish - t_gpu_copyandkernel_begin) /
                  CLOCKS_PER_SEC);

      double t_cpu_elapsed = 0.f;
      for (int tid = 0; tid < count_tasks; tid++) {
        //      	switch (tasks[tid]->type) {
        //		     case task_type_self:
        if (tasks[tid]->subtype == task_subtype_density) {
          double t1 = clock();
          runner_doself1_branch_density(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        if (tasks[tid]->subtype == task_subtype_gradient) {
          double t1 = clock();
          runner_doself1_branch_gradient(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        if (tasks[tid]->subtype == task_subtype_force) {
          double t1 = clock();
          runner_doself2_branch_force(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        //      	}
      }
      double t_cpu_print_start = clock();
      fclose(fpcpu);
      fprintf(stderr, "Printing CPU results\n");
      fpcpu = fopen("./res_CPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_list[tid]->hydro.count;
        for (int p = 0; p < count; p++) {
          float xx = ci_list[tid]->hydro.parts[p].x[0],
                yy = ci_list[tid]->hydro.parts[p].x[1],
                zz = ci_list[tid]->hydro.parts[p].x[2];
          fprintf(fpcpu, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
                  ci_list[tid]->hydro.parts[p].rho,
                  ci_list[tid]->hydro.parts[p].density.rho_dh,
                  ci_list[tid]->hydro.parts[p].density.wcount,
                  ci_list[tid]->hydro.parts[p].density.wcount_dh,
                  ci_list[tid]->hydro.parts[p].viscosity.div_v,
                  ci_list[tid]->hydro.parts[p].density.rot_v[0]);
        }
      }

      fclose(fpcpu);
      double t_cpu_print_end = clock();

      fprintf(stderr, "Finished printing GPU results\n");
      //		for (int tid = 0; tid < count_tasks; ++tid)
      //			gpuErrchk( cudaStreamDestroy(stream[tid]));
      double t_gpu_end = clock();
      printf("time for cpu is %f\n", t_cpu_elapsed / CLOCKS_PER_SEC);
      printf("time for gpu is %f\n", (t_gpu_end - t_gpu_start - t_cpu_elapsed -
                                      t_cpu_print_end + t_cpu_print_start) /
                                         CLOCKS_PER_SEC);
      printf("finished looping through tasks\n");
      FILE *fp;
      fp = fopen("./res_GPU_density.txt", "w");
      //		for(int tid=0; tid<count_tasks; tid++){
      //			const int
      //count=task_last_part[tid]-task_first_part[tid];
      double t6 = clock();
      for (int p = 0; p < count_all_parts; p++) {
        float xx = parts_gpu_CSR[p].x[0], yy = parts_gpu_CSR[p].x[1],
              zz = parts_gpu_CSR[p].x[2];
        fprintf(fp, "%f %f %f %f %f %f %f %f %f %i\n", xx, yy, zz,
                parts_gpu_CSR[p].rho, parts_gpu_CSR[p].rho_dh,
                parts_gpu_CSR[p].wcount, parts_gpu_CSR[p].wcount_dh,
                parts_gpu_CSR[p].div_v, parts_gpu_CSR[p].rot_v[0],
                parts_gpu_CSR[p].tid);
      }
      //		}
      fclose(fp);
      //      cudaFreeHost(ci_gpu_list);
      cudaFreeHost(parts_gpu_list);
      cudaProfilerStop();
      cudaDeviceReset();
      exit(0);
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////****same as runner_main_cuda but bunching tasks into bundles. Each stream
///works on one bundle***********/////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void *runner_main_cuda_bundles_revised_soa(void *data) {
  cudaDeviceReset();
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  cudaError_t cu_error;
  cudaMemPool_t memPool;
  gpuErrchk(cudaDeviceGetDefaultMemPool(&memPool, devId));
  int maxmem = 0.9f * prop.totalGlobalMem;
  gpuErrchk(cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold,
                                    (void *)&maxmem));
  int counter = 0;

  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    //	 gpuErrchk( cudaMallocHost((void**)&tasks, n_streams*sizeof(struct
    //task)) ); // Pinned allocation on host Typecasting from void** to struct
    // task **//

    ////////////////////////////////////////////////////////////////////////////////
    ///////This is likely where the problem is!!!!
    ////////////////////////////////////////////////////////////////////////////////
    float time_to_get_tasks = clock();
    struct task *tasks[n_streams];

    while (1) {
      for (int i = 0; i < n_streams; i++) {
        tasks[i] = NULL;
      }
      // for(int tid=0; tid<n_streams; tid++){
      int count_tasks = 0;
      int count_stale = 0;
      ////////////////////////////////////
      /*Grab a bunch of tasks*/
      while (1) { //(count_tasks<n_streams){
        // there's a potential bug here. Code hangs if n_streams set to > 256 in
        // cuda_headers.h
        /* If there's no old task, try to get a new one. */
        if (tasks[count_tasks] == NULL) {
          /* Get the task. */
          TIMER_TIC
          tasks[count_tasks] = scheduler_gettask(sched, r->qid, prev);
          TIMER_TOC(timer_gettask);
          if (tasks[count_tasks] != NULL) {
            count_tasks++;
          }
        }
        /* If there's an old task, move onto the next one */
        // else count_tasks++;
        count_stale++;
        if (count_stale >= n_streams) {
          break;
        }
      }
      fprintf(stderr, "Count tasks is %i\n", count_tasks);
      ////////////////////////////////////

      ////////////////////////////////////
      /* Get the cells. */
      struct cell **ci_list =
          (struct cell **)malloc(count_tasks * sizeof(struct cell *));
      int count_all_parts = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        //			ci_list[tid]=(struct cell) malloc(sizeof(struct
        //cell*));
        ci_list[tid] = tasks[tid]->ci;
        struct cell *c_tmp = ci_list[tid];
        count_all_parts += ci_list[tid]->hydro.count;
      }

      ////////////////////////////////////

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->type == task_type_pair ||
            tasks[tid]->type == task_type_sub_pair) {
          struct cell *ci_temp = ci_list[tid];
          //		     struct cell *cj_temp = cj;
          double shift[3];
          tasks[tid]->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        } else {
          tasks[tid]->sid = -1;
        }
      }
#endif
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
      i.e., rather than just one r->t, r would have a nested bunch of structs
       r->t[size=nStreams].*/
      r->t = t;
#endif
      time_to_get_tasks = (clock() - time_to_get_tasks) / CLOCKS_PER_SEC;
      fprintf(stderr, "time to get tasks is %f\n", time_to_get_tasks);
      float time_to_malloc_variables = clock();
      //////////////////////////////////////////////
      /*****Initialise GPU suitable copies of variables and malloc pinned
       * memory*****/
      //		struct part_gpu **
      //parts_gpu_list=malloc(count_tasks*sizeof(struct part_gpu *)); 		struct
      //part_gpu * parts_gpu_CSR;
      part_soa parts_soa;
      part_soa d_parts_soa;
      ////////All data contained within parts_gpu_CSR_soa struct. A lot, innit?
      int *d_tid_p;
      long long *d_id;
      double *d_x_p;
      double *d_y_p;
      double *d_z_p;
      float *d_ux;
      float *d_uy;
      float *d_uz;
      float *d_a_hydrox;
      float *d_a_hydroy;
      float *d_a_hydroz;
      float *d_mass;
      float *d_h;
      float *d_u;
      float *d_u_dt;
      float *d_rho;
      float *d_SPH_sum;
      float *d_locx;
      float *d_locy;
      float *d_locz;
      float *d_widthx;
      float *d_widthy;
      float *d_widthz;
      float *d_h_max;
      int *d_count_p;
      float *d_wcount;
      float *d_wcount_dh;
      float *d_rho_dh;
      float *d_rot_ux;
      float *d_rot_uy;
      float *d_rot_uz;
      float *d_div_v;
      float *d_div_v_previous_step;
      float *d_alpha_visc;
      float *d_v_sig;
      float *d_laplace_u;
      float *d_alpha_diff;
      float *d_f;
      float *d_soundspeed;
      float *d_h_dt;
      float *d_balsara;
      float *d_pressure;
      float *d_alpha_visc_max_ngb;
      /* timestep stuff */
      timebin_t *d_time_bin;
      timebin_t *d_wakeup;
      timebin_t *d_min_ngb_time_bin;
      char *d_to_be_synchronized;
      ////////now malloc all this data. Sheesh
      fprintf(stderr, "before mallocing soas\n");
      cudaMalloc((void **)&(d_tid_p), sizeof(int) * count_all_parts);
      cudaMalloc((void **)&(d_id), sizeof(long long) * count_all_parts);
      cudaMalloc((void **)&(d_x_p), sizeof(double) * count_all_parts);
      cudaMalloc((void **)&(d_y_p), sizeof(double) * count_all_parts);
      cudaMalloc((void **)&(d_z_p), sizeof(double) * count_all_parts);
      cudaMalloc((void **)&(d_ux), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_uy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_uz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_a_hydrox), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_a_hydroy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_a_hydroz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_mass), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_h), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_u), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_u_dt), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rho), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_SPH_sum), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_locx), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_locy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_locz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_widthx), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_widthy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_widthz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_h_max), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_count_p), sizeof(int) * count_all_parts);
      cudaMalloc((void **)&(d_wcount), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_wcount_dh), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rho_dh), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rot_ux), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rot_uy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rot_uz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_div_v), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_div_v_previous_step),
                 sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_alpha_visc), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_v_sig), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_laplace_u), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_alpha_diff), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_f), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_soundspeed), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_h_dt), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_balsara), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_pressure), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_alpha_visc_max_ngb),
                 sizeof(float) * count_all_parts);
      /* timestep stuff */
      cudaMalloc((void **)&(d_time_bin), sizeof(timebin_t) * count_all_parts);
      cudaMalloc((void **)&(d_wakeup), sizeof(timebin_t) * count_all_parts);
      cudaMalloc((void **)&(d_min_ngb_time_bin),
                 sizeof(timebin_t) * count_all_parts);
      cudaMalloc((void **)&(d_to_be_synchronized),
                 sizeof(char) * count_all_parts);
      fprintf(stderr, "after mallocing soas\n");
      ///////////Host arrays
      //		int *tid_p=malloc(count_all_parts*sizeof(int));
      int *tid_p;
      cudaMallocHost((void **)&tid_p, count_all_parts * sizeof(int));
      long long *id;
      cudaMallocHost((void **)&id, count_all_parts * sizeof(long long));
      float *mass;
      cudaMallocHost((void **)&mass, count_all_parts * sizeof(float));
      float *h;
      cudaMallocHost((void **)&h, count_all_parts * sizeof(float));
      float *u;
      cudaMallocHost((void **)&u, count_all_parts * sizeof(float));
      float *u_dt;
      cudaMallocHost((void **)&u_dt, count_all_parts * sizeof(float));
      float *rho;
      cudaMallocHost((void **)&rho, count_all_parts * sizeof(float));
      float *SPH_sum;
      cudaMallocHost((void **)&SPH_sum, count_all_parts * sizeof(float));
      //		double x_p[count_all_parts];
      double *x_p;
      cudaMallocHost((void **)&x_p, count_all_parts * sizeof(double));
      double *y_p;
      cudaMallocHost((void **)&y_p, count_all_parts * sizeof(double));
      double *z_p;
      cudaMallocHost((void **)&z_p, count_all_parts * sizeof(double));
      float *ux;
      cudaMallocHost((void **)&ux, count_all_parts * sizeof(float));
      float *uy;
      cudaMallocHost((void **)&uy, count_all_parts * sizeof(float));
      float *uz;
      cudaMallocHost((void **)&uz, count_all_parts * sizeof(float));
      float *a_hydrox;
      cudaMallocHost((void **)&a_hydrox, count_all_parts * sizeof(float));
      float *a_hydroy;
      cudaMallocHost((void **)&a_hydroy, count_all_parts * sizeof(float));
      float *a_hydroz;
      cudaMallocHost((void **)&a_hydroz, count_all_parts * sizeof(float));
      float *locx;
      cudaMallocHost((void **)&locx, count_all_parts * sizeof(float));
      float *locy;
      cudaMallocHost((void **)&locy, count_all_parts * sizeof(float));
      float *locz;
      cudaMallocHost((void **)&locz, count_all_parts * sizeof(float));
      float *widthx;
      cudaMallocHost((void **)&widthx, count_all_parts * sizeof(float));
      float *widthy;
      cudaMallocHost((void **)&widthy, count_all_parts * sizeof(float));
      float *widthz;
      cudaMallocHost((void **)&widthz, count_all_parts * sizeof(float));
      float *h_max;
      cudaMallocHost((void **)&h_max, count_all_parts * sizeof(float));
      int *count_p;
      cudaMallocHost((void **)&count_p, count_all_parts * sizeof(int));
      float *wcount;
      cudaMallocHost((void **)&wcount, count_all_parts * sizeof(float));
      float *wcount_dh;
      cudaMallocHost((void **)&wcount_dh, count_all_parts * sizeof(float));
      float *rho_dh;
      cudaMallocHost((void **)&rho_dh, count_all_parts * sizeof(float));
      float *rot_ux;
      cudaMallocHost((void **)&rot_ux, count_all_parts * sizeof(float));
      float *rot_uy;
      cudaMallocHost((void **)&rot_uy, count_all_parts * sizeof(float));
      float *rot_uz;
      cudaMallocHost((void **)&rot_uz, count_all_parts * sizeof(float));
      float *div_v;
      cudaMallocHost((void **)&div_v, count_all_parts * sizeof(float));
      float *div_v_previous_step;
      cudaMallocHost((void **)&div_v_previous_step,
                     count_all_parts * sizeof(float));
      float *alpha_visc;
      cudaMallocHost((void **)&alpha_visc, count_all_parts * sizeof(float));
      float *v_sig;
      cudaMallocHost((void **)&v_sig, count_all_parts * sizeof(float));
      float *laplace_u;
      cudaMallocHost((void **)&laplace_u, count_all_parts * sizeof(float));
      float *alpha_diff;
      cudaMallocHost((void **)&alpha_diff, count_all_parts * sizeof(float));
      float *f;
      cudaMallocHost((void **)&f, count_all_parts * sizeof(float));
      float *soundspeed;
      cudaMallocHost((void **)&soundspeed, count_all_parts * sizeof(float));
      float *h_dt;
      cudaMallocHost((void **)&h_dt, count_all_parts * sizeof(float));
      float *balsara;
      cudaMallocHost((void **)&balsara, count_all_parts * sizeof(float));
      float *pressure;
      cudaMallocHost((void **)&pressure, count_all_parts * sizeof(float));
      float *alpha_visc_max_ngb;
      cudaMallocHost((void **)&alpha_visc_max_ngb,
                     count_all_parts * sizeof(float));
      /* timestep stuff */
      timebin_t *time_bin;
      cudaMallocHost((void **)&time_bin, count_all_parts * sizeof(timebin_t));
      timebin_t *wakeup;
      cudaMallocHost((void **)&wakeup, count_all_parts * sizeof(timebin_t));
      timebin_t *min_ngb_time_bin;
      cudaMallocHost((void **)&min_ngb_time_bin,
                     count_all_parts * sizeof(timebin_t));
      char *to_be_synchronized;
      cudaMallocHost((void **)&to_be_synchronized,
                     count_all_parts * sizeof(char));
      ///////////////////////////////////////////////////

      /////////////Host arrays
      ////		int *tid_p=malloc(count_all_parts*sizeof(int));
      //		int *tid_p;
      //		gpuErrchk( cudaMallocHost((void**)&tid_p,
      //count_all_parts*sizeof(int)) ); 		long long
      //*id=malloc(count_all_parts*sizeof(long long)); 		float
      //*mass=malloc(count_all_parts*sizeof(float)); 		float
      //*h=malloc(count_all_parts*sizeof(float)); 		float
      //*u=malloc(count_all_parts*sizeof(float)); 		float
      //*u_dt=malloc(count_all_parts*sizeof(float)); 		float
      //*rho=malloc(count_all_parts*sizeof(float)); 		float
      //*SPH_sum=malloc(count_all_parts*sizeof(float));
      ////		double x_p[count_all_parts];
      //		double *x_p=malloc(count_all_parts*sizeof(double));
      //		double *y_p=malloc(count_all_parts*sizeof(double));
      //		double *z_p=malloc(count_all_parts*sizeof(double));
      //		float *ux=malloc(count_all_parts*sizeof(float));
      //		float *uy=malloc(count_all_parts*sizeof(float));
      //		float *uz=malloc(count_all_parts*sizeof(float));
      //		float *a_hydrox=malloc(count_all_parts*sizeof(float));
      //		float *a_hydroy=malloc(count_all_parts*sizeof(float));
      //		float *a_hydroz=malloc(count_all_parts*sizeof(float));
      //		float *locx=malloc(count_all_parts*sizeof(float));
      //		float *locy=malloc(count_all_parts*sizeof(float));
      //		float *locz=malloc(count_all_parts*sizeof(float));
      //		float *widthx=malloc(count_all_parts*sizeof(float));
      //		float *widthy=malloc(count_all_parts*sizeof(float));
      //		float *widthz=malloc(count_all_parts*sizeof(float));
      //		float *h_max=malloc(count_all_parts*sizeof(float));
      //		int *count_p=malloc(count_all_parts*sizeof(int));
      //		float *wcount=malloc(count_all_parts*sizeof(float));
      //		float *wcount_dh=malloc(count_all_parts*sizeof(float));
      //		float *rho_dh=malloc(count_all_parts*sizeof(float));
      //		float *rot_ux=malloc(count_all_parts*sizeof(float));
      //		float *rot_uy=malloc(count_all_parts*sizeof(float));
      //		float *rot_uz=malloc(count_all_parts*sizeof(float));
      //		float *div_v=malloc(count_all_parts*sizeof(float));
      //		float
      //*div_v_previous_step=malloc(count_all_parts*sizeof(float)); 		float
      //*alpha_visc=malloc(count_all_parts*sizeof(float)); 		float
      //*v_sig=malloc(count_all_parts*sizeof(float)); 		float
      //*laplace_u=malloc(count_all_parts*sizeof(float)); 		float
      //*alpha_diff=malloc(count_all_parts*sizeof(float)); 		float
      //*f=malloc(count_all_parts*sizeof(float)); 		float
      //*soundspeed=malloc(count_all_parts*sizeof(float)); 		float
      //*h_dt=malloc(count_all_parts*sizeof(float)); 		float
      //*balsara=malloc(count_all_parts*sizeof(float)); 		float
      //*pressure=malloc(count_all_parts*sizeof(float)); 		float
      //*alpha_visc_max_ngb=malloc(count_all_parts*sizeof(float));
      //		/* timestep stuff */
      //		timebin_t
      //*time_bin=malloc(count_all_parts*sizeof(timebin_t)); 		timebin_t
      //*wakeup=malloc(count_all_parts*sizeof(timebin_t)); 		timebin_t
      //*min_ngb_time_bin=malloc(count_all_parts*sizeof(timebin_t)); 		char
      //*to_be_synchronized=malloc(count_all_parts*sizeof(char));
      /////////////////////////////////////////////////////

      fprintf(stderr, "Size I'm allocating is %i\n",
              count_all_parts * sizeof(struct part_gpu));
      int *task_first_part, *task_last_part;
      int *d_task_first_part, *d_task_last_part;
      cudaMallocHost((void **)&task_first_part,
                     count_tasks * sizeof(int)); // Pinned allocation on host
      cudaMallocHost((void **)&task_last_part,
                     count_tasks * sizeof(int)); // Pinned allocation on host

      int bundle_size = 32;
      int nBundles = (count_tasks + bundle_size - 1) / bundle_size;

      int *bundle_first_part, *bundle_last_part;
      int *d_bundle_first_part, *d_bundle_last_part;
      cudaMallocHost((void **)&bundle_first_part,
                     nBundles * sizeof(int)); // Pinned allocation on host
      cudaMallocHost((void **)&bundle_last_part,
                     nBundles * sizeof(int)); // Pinned allocation on host

      cudaMalloc((void **)&d_task_first_part, count_tasks * sizeof(int));
      cudaMalloc((void **)&d_task_last_part, count_tasks * sizeof(int));
      cudaMalloc((void **)&d_bundle_first_part, nBundles * sizeof(int));
      cudaMalloc((void **)&d_bundle_last_part, nBundles * sizeof(int));
      int total_num_parts = 0;
      ///////Start GPU timer
      //		double t_gpu_start=clock();
      ////////////////////////////////////////////
      time_to_malloc_variables =
          (clock() - time_to_malloc_variables) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_to_malloc_variables %f\n",
              time_to_malloc_variables);
      ////////Malloc on HOST
      float time_for_data_arrangement = clock();
      int first_part_tmp = 0;
      //		fprintf(stderr,"Before populate parts\n");
      int count_prev = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        int count = ci_list[tid]->hydro.count;
        //			fprintf(stderr, "Just before popluator\n");
        ///////////////////////////////////////////////////////////////////////
        //////////////This should probably be done on the GPU??////////////////
        ///////////////////////////////////////////////////////////////////////
        //			struct cell *ci_temp = ci_list[tid];
        populate_parts_list_soa(
            count_all_parts, ci_list[tid], first_part_tmp, count, tid, tid_p,
            id, x_p, y_p, z_p, ux, uy, uz, a_hydrox, a_hydroy, a_hydroz, mass,
            h, u, u_dt, rho, SPH_sum, locx, locy, locz, widthx, widthy, widthz,
            h_max, count_p, wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
            div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
            alpha_diff, f, soundspeed, h_dt, balsara, pressure,
            alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
            to_be_synchronized);
        //			fprintf(stderr, "Just after popluator\n");
        //			populate_parts_list_soa(count_all_parts,
        //ci_list[tid], first_part_tmp, count, &tid_p, &id, &x_p, &y_p, &z_p,
        //&ux, &uy, &uz, &a_hydrox, 									&a_hydroy, &a_hydroz, &mass, &h, &u, &u_dt,
        //&rho, &SPH_sum, &locx, &locy, &locz, 									&widthx, &widthy, &widthz,
        //&h_max, &count_p, &wcount, &wcount_dh, &rho_dh, &rot_u, &rot_v,
        //									&rot_w, &div_v,
        //&div_v_previous_step, &alpha_visc, &v_sig, &laplace_u, &alpha_diff,
        //&f, &soundspeed, 									&h_dt, &balsara, &pressure, &alpha_visc_max_ngb,
        //&time_bin, &wakeup, &min_ngb_time_bin, 									&to_be_synchronized);
        //			fprintf(stderr, "Just after popluator\n");
        //			for(int p=0; p<count; p++){
        //				int p_gid=p+first_part_tmp;
        //				x_p[p_gid]=ci_list[tid]->hydro.parts[p].x[0];
        //				fprintf(stderr,"part is %i x is %f x from cell
        //struct is %f\n",p, x_p[p_gid], ci_list[tid]->hydro.parts[p].x[0]);
        //		    }
        task_first_part[tid] = first_part_tmp;
        task_last_part[tid] = first_part_tmp + count;
        first_part_tmp += count;
        if (tid % bundle_size == 0) {
          int bid = tid / bundle_size;
          bundle_first_part[bid] = task_first_part[tid];
          //				fprintf(stderr, "bundle id %i nBundles %i\n", bid,
          //nBundles);
        }
      }
      for (int bid = 0; bid < nBundles - 1; bid++) {
        bundle_last_part[bid] = bundle_first_part[bid + 1];
      }
      bundle_last_part[nBundles - 1] =
          count_all_parts; //-bundle_last_part[nBundles-1];
                           //		for(int bid=0; bid<nBundles; bid++){
                           //			fprintf(stderr, "count_all_parts %i, bid %i,
      //nBundles %i, first part %i, last part %i\n", count_all_parts, bid,
      //nBundles, bundle_first_part[bid], bundle_last_part[bid]);
      //		}
      //		exit(0);
      //		fprintf(stderr, "GOT HERE\n");
      time_for_data_arrangement =
          (clock() - time_for_data_arrangement) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_data_arrangement %f\n",
              time_for_data_arrangement);
      //		parts_soa.tid_p=tid_p;
      //	    parts_soa.locx=locx;
      //	    parts_soa.locy=locy;
      //	    parts_soa.locz=locz;
      //	    parts_soa.h=h;
      //	    parts_soa.mass=mass;
      //	    parts_soa.x_p=x_p;
      //	    parts_soa.y_p=y_p;
      //	    parts_soa.z_p=z_p;
      //	    parts_soa.rho=rho;
      //	    parts_soa.rho_dh=rho_dh;
      //	    parts_soa.wcount=wcount;
      //	    parts_soa.wcount_dh=wcount_dh;
      //	    parts_soa.ux=ux;
      //	    parts_soa.uy=uy;
      //	    parts_soa.uz=uz;
      //	    parts_soa.div_v=div_v;
      //	    parts_soa.rot_ux=rot_ux;
      //	    parts_soa.rot_uy=rot_uy;
      //	    parts_soa.rot_uz=rot_uz;
      //	    parts_soa.count_p=count_p;
      //	    for(int p=0; p<count_all_parts; p++){
      //	    	parts_soa.count_p[p]=0;
      //	    }
      //	    fprintf(stderr,"got here 2\n");
      //		struct part_gpu * d_parts_gpu_CSR;
      //		cudaMalloc((void**)&d_parts_gpu_CSR,
      //count_all_parts*sizeof(struct part_gpu));

      //////////////////////////////////////////////

      //////Prepare streams and malloc memory on the GPU
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      // create events and streams
      ///////Start GPU timer
      double t_gpu_start = clock();
      ////////////////////////////////////////////

      ////////////////////////////////////////////////
      ///////////Malloc GPU data

      cudaProfilerStart();
      double t_gpu_malloc_begin = clock();

      FILE *fpcpu, *fpgpu;
      fpcpu = fopen("./res_CPU_density.txt", "w");
      double t_cpu_start = clock();
      double t_gpu_elapsed = 0.f;

      // see
      // https://stackoverflow.com/questions/23609770/cuda-double-pointer-memory-copy
      // for details of how this is allocated and copied

      //	   fprintf(stderr,"Malloc d_all_parts\n");
      cudaError_t err;
      //	   for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=ci_list[tid]->hydro.count;
      //
      //		   gpuErrchk(cudaMemcpy(&d_tid_p[first_part_tmp],
      //&tid_p[first_part_tmp], count*sizeof(int), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_locx[first_part_tmp],
      //&locx[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_locy[first_part_tmp],
      //&locy[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_locz[first_part_tmp],
      //&locz[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_h[first_part_tmp],
      //&h[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_mass[first_part_tmp],
      //&mass[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_x_p[first_part_tmp],
      //&x_p[first_part_tmp], count*sizeof(double), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_y_p[first_part_tmp],
      //&y_p[first_part_tmp], count*sizeof(double), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_z_p[first_part_tmp],
      //&z_p[first_part_tmp], count*sizeof(double), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rho[first_part_tmp],
      //&rho[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rho_dh[first_part_tmp],
      //&rho_dh[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_wcount[first_part_tmp],
      //&wcount[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_wcount_dh[first_part_tmp],
      //&wcount_dh[first_part_tmp], count*sizeof(float),
      //cudaMemcpyHostToDevice)); 		   gpuErrchk(cudaMemcpy(&d_ux[first_part_tmp],
      //&ux[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_uy[first_part_tmp],
      //&uy[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_uz[first_part_tmp],
      //&uz[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_div_v[first_part_tmp],
      //&div_v[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rot_ux[first_part_tmp],
      //&rot_ux[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rot_uy[first_part_tmp],
      //&rot_uy[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rot_uz[first_part_tmp],
      //&rot_uz[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_count_p[first_part_tmp],
      //&count_p[first_part_tmp], count*sizeof(int), cudaMemcpyHostToDevice));
      //
      //		   err = cudaPeekAtLastError();//cudaGetLastError(); //
      //Get error code 		   if ( err != cudaSuccess )
      //			{
      //			   fprintf(stderr, "CUDA kernel Error first
      //memcpy11: %s \n ", cudaGetErrorString(err));
      //			}
      //	   }

      //	   for(int p=0; p<count_all_parts; p++){
      //		   fprintf(stderr,"x position is
      //%f\n",parts_soa.x_p[p]); 		   fprintf(stderr,"tid is
      //%i\n",parts_soa.tid_p[p]);
      //	   }

      //	   FILE *fptest = fopen("./res_GPU_density.txt", "w");
      //	   	   for(int p=0; p<count_all_parts; p++){
      //	   			float xx=parts_soa.x_p[p],
      //yy=parts_soa.y_p[p], zz=parts_soa.z_p[p];
      //	   //			fprintf(stderr,"x is %f y is %f z is
      //%f\n", xx, yy, zz); 	   			fprintf(fptest, "%f %f %f %f %f %f %f %f %f %i\n",
      //xx, yy, zz, parts_soa.rho[p], 	   					parts_soa.rho_dh[p], parts_soa.wcount[p],
      //parts_soa.wcount_dh[p], 	   					parts_soa.div_v[p], parts_soa.rot_ux[p],
      //parts_soa.tid_p[p]);
      //	   	   }
      //	   	   fprintf(stderr,"Finished printing GPU results\n");
      //	   	   //		}
      //	   	   fclose(fptest);
      //	   	FILE *fptestcpu = fopen("./res_CPU_density.txt", "w");
      //	   		   for(int tid=0; tid<count_tasks; tid++){
      //	   			int count = ci_list[tid]->hydro.count;
      //	   			for (int p=0; p<count; p++){
      //	   				float
      //xx=ci_list[tid]->hydro.parts[p].x[0],
      //yy=ci_list[tid]->hydro.parts[p].x[1],
      //zz=ci_list[tid]->hydro.parts[p].x[2]; 					fprintf(fptestcpu, "%f %f %f %f %f
      //%f %f %f %f\n", xx, yy, zz, ci_list[tid]->hydro.parts[p].rho,
      //					ci_list[tid]->hydro.parts[p].density.rho_dh,
      //ci_list[tid]->hydro.parts[p].density.wcount,
      //					ci_list[tid]->hydro.parts[p].density.wcount_dh,
      //ci_list[tid]->hydro.parts[p].viscosity.div_v 					,
      //ci_list[tid]->hydro.parts[p].density.rot_v[0]);
      //	   			}
      //	   		   }
      //	   		   fprintf(stderr,"Finished printing GPU
      //results\n");
      //	   		   //		}
      //	   		   fclose(fptestcpu);
      //	   exit(0);
      float time_for_extraneous_memcpys = clock();
      parts_soa.tid_p = d_tid_p;
      parts_soa.locx = d_locx;
      parts_soa.locy = d_locy;
      parts_soa.locz = d_locz;
      parts_soa.h = d_h;
      parts_soa.mass = d_mass;
      parts_soa.x_p = d_x_p;
      parts_soa.y_p = d_y_p;
      parts_soa.z_p = d_z_p;
      parts_soa.rho = d_rho;
      parts_soa.rho_dh = d_rho_dh;
      parts_soa.wcount = d_wcount;
      parts_soa.wcount_dh = d_wcount_dh;
      parts_soa.ux = d_ux;
      parts_soa.uy = d_uy;
      parts_soa.uz = d_uz;
      parts_soa.div_v = d_div_v;
      parts_soa.rot_ux = d_rot_ux;
      parts_soa.rot_uy = d_rot_uy;
      parts_soa.rot_uz = d_rot_uz;
      parts_soa.count_p = d_count_p;

      //	   fprintf(stderr,"number of parts %i, first part in last task
      //%i, last part in last task %i\n", count_all_parts,
      //task_first_part[count_tasks-1], task_last_part[count_tasks-1]);
      //	   fprintf(stderr,"number of parts %i, first part in fist bundle
      //%i, last part in last bundle %i\n", count_all_parts,
      //bundle_first_part[nBundles-1], bundle_last_part[nBundles-1]);

      cudaMemcpy(d_task_first_part, task_first_part, count_tasks * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy1: %s
      //\n ", cudaGetErrorString(err));
      //		}
      cudaMemcpy(d_task_last_part, task_last_part, count_tasks * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy2: %s
      //\n ", cudaGetErrorString(err));
      //		}
      cudaMemcpy(d_bundle_first_part, bundle_first_part, nBundles * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy3: %s
      //\n ", cudaGetErrorString(err));
      //		}
      cudaMemcpy(d_bundle_last_part, bundle_last_part, nBundles * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy4: %s
      //\n ", cudaGetErrorString(err));
      //		}
      time_for_extraneous_memcpys =
          (clock() - time_for_extraneous_memcpys) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_extraneous_memcpys is %f\n",
              time_for_extraneous_memcpys);
      // For another possible solution, see
      // https://forums.developer.nvidia.com/t/how-do-i-pass-a-double-pointers-array-to-the-device-im-getting-cudaerrorillegaladdress/72518/2

      int size_parts_mmcpy = 0;
      float time_for_stream_creation = clock();
      cudaEvent_t startEvent, stopEvent, dummyEvent;

      cudaEventCreate(&startEvent);
      cudaEventCreate(&stopEvent);
      cudaEventCreate(&dummyEvent);
      //		fprintf(stderr,"got here in the end\n");
      double t_gpu_malloc_finish = clock();
      //		fprintf(stderr,"malloc time is
      //%f\n",(t_gpu_malloc_finish-t_gpu_malloc_begin)/CLOCKS_PER_SEC);

      cudaStream_t stream[nBundles];
      int tasksperbundle = (count_tasks + nBundles - 1) / nBundles;
      //		fprintf(stderr,"count_tasks/nBundles %i, tasksperbundle
      //%i\n", count_tasks/nBundles, tasksperbundle); 		exit(0);
      for (int i = 0; i < nBundles; ++i)
        //			gpuErrchk( cudaStreamCreate(&stream[i]) );
        cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking);
      //		err = cudaPeekAtLastError();//cudaGetLastError(); // Get
      //error code 		if ( err != cudaSuccess )
      //		{
      //		  fprintf(stderr, "CUDA Stream Creation Error: %s \n ",
      //cudaGetErrorString(err));
      //		}
      time_for_stream_creation =
          (clock() - time_for_stream_creation) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_stream_creation is %f\n",
              time_for_stream_creation);
      t_gpu_elapsed = 0.f;
      float time_for_memcpys_and_kernel = clock();
      //		for(int bid=0; bid<nBundles; bid++){
      //		   int first_part_tmp=bundle_first_part[bid];
      //		   int bundle_size=bundle_last_part[bid]-first_part_tmp;
      //
      //		   gpuErrchk(cudaMemcpyAsync(&d_tid_p[first_part_tmp],
      //&tid_p[first_part_tmp], bundle_size*sizeof(int), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_locx[first_part_tmp],
      //&locx[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_locy[first_part_tmp],
      //&locy[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_locz[first_part_tmp],
      //&locz[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_h[first_part_tmp],
      //&h[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_mass[first_part_tmp],
      //&mass[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_x_p[first_part_tmp],
      //&x_p[first_part_tmp], bundle_size*sizeof(double),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_y_p[first_part_tmp],
      //&y_p[first_part_tmp], bundle_size*sizeof(double),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_z_p[first_part_tmp],
      //&z_p[first_part_tmp], bundle_size*sizeof(double),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rho[first_part_tmp],
      //&rho[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_rho_dh[first_part_tmp],
      //&rho_dh[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_wcount[first_part_tmp],
      //&wcount[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_wcount_dh[first_part_tmp],
      //&wcount_dh[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_ux[first_part_tmp],
      //&ux[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_uy[first_part_tmp],
      //&uy[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_uz[first_part_tmp],
      //&uz[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_div_v[first_part_tmp],
      //&div_v[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rot_ux[first_part_tmp],
      //&rot_ux[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rot_uy[first_part_tmp],
      //&rot_uy[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rot_uz[first_part_tmp],
      //&rot_uz[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_count_p[first_part_tmp],
      //&count_p[first_part_tmp], bundle_size*sizeof(int),
      //cudaMemcpyHostToDevice, stream[bid]));
      //
      //		   err = cudaPeekAtLastError();//cudaGetLastError(); //
      //Get error code 		   if ( err != cudaSuccess )
      //			{
      //			   fprintf(stderr, "CUDA kernel Error first
      //memcpy11: %s \n ", cudaGetErrorString(err));
      //			}
      //		}
      int shared_memory_used = 0;
      double t_gpu_copyandkernel_begin = clock();
      for (int bid = 0; bid < nBundles; bid++) {
        int max_parts = 0;
        int parts_in_bundle = 0;
        for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
             tid++) {
          if (tid < count_tasks) {
            int count = task_last_part[tid] - task_first_part[tid];
            //					if(tid<count_tasks-1)
            parts_in_bundle += count;
            max_parts = max(max_parts, count);
          }
        }

        int first_part_tmp = bundle_first_part[bid];
        int bundle_size = bundle_last_part[bid] - first_part_tmp;
        cudaMemcpyAsync(&d_tid_p[first_part_tmp], &tid_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_locx[first_part_tmp], &locx[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_locy[first_part_tmp], &locy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_locz[first_part_tmp], &locz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_h[first_part_tmp], &h[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_mass[first_part_tmp], &mass[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_x_p[first_part_tmp], &x_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_y_p[first_part_tmp], &y_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_z_p[first_part_tmp], &z_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rho[first_part_tmp], &rho[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rho_dh[first_part_tmp], &rho_dh[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_wcount[first_part_tmp], &wcount[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_wcount_dh[first_part_tmp],
                        &wcount_dh[first_part_tmp], bundle_size * sizeof(float),
                        cudaMemcpyHostToDevice, stream[bid]);
        cudaMemcpyAsync(&d_ux[first_part_tmp], &ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_uy[first_part_tmp], &uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_uz[first_part_tmp], &uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_div_v[first_part_tmp], &div_v[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rot_ux[first_part_tmp], &rot_ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rot_uy[first_part_tmp], &rot_uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rot_uz[first_part_tmp], &rot_uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_count_p[first_part_tmp], &count_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyHostToDevice,
                        stream[bid]);
        int tid = 0;
        int offset = bid * tasksperbundle;
        int sizecells = 0, sizeparts = 0;

        sizeparts = tasksperbundle * sizeof(struct part_gpu *);
        int tasks_left = tasksperbundle;
        if (bid == nBundles - 1) {
          tasks_left = count_tasks - (nBundles - 1) * tasksperbundle;
          sizeparts = tasks_left * sizeof(struct part_gpu *);
        }
        int numBlocks_y = tasks_left;
        int block_size = BLOCK_SIZE;
        //			fprintf(stderr, "block size in runner_main is %i\n",
        //block_size);
        int numBlocks_x = (max_parts + block_size - 1) / block_size;
        int numBlocks = (parts_in_bundle + block_size - 1) / block_size;
        //			fprintf(stderr, "numblocks %i, parts in bundle
        //%i\n", numBlocks, parts_in_bundle); 			fprintf(stderr, "bid is %i, max
        //parts/cell in bundle %i Ntasks is %i tasksperbundle is %i\n NBlocks x
        //is %i NBlocksy is %i numBundles is %i\n", bid, max_parts, count_tasks,
        //tasksperbundle, numBlocks_x, numBlocks_y, nBundles);
        const char *loop_type = "density";
        int shared_memory_block =
            BLOCK_SIZE * (8 * sizeof(float) + sizeof(float *));
        int bundle_part_0 = bundle_first_part[bid];
        int bundle_first_task = tid_p[bundle_part_0];
        //			fprintf(stderr, "size required per block is %i,
        //total requested for kernel is %i\n", shared_memory_block,
        //shared_memory_block*numBlocks);
        launch_cuda_kernel_bundles_revised_soa(
            parts_soa, d_task_first_part, d_task_last_part, d_bundle_first_part,
            d_bundle_last_part, numBlocks, d_a, d_H, loop_type, stream[bid],
            bid, block_size, count_tasks, tasksperbundle, numBlocks_x,
            numBlocks_y, tid, offset, bundle_first_task, max_parts);
        //			cudaDeviceSynchronize();
        //			err = cudaPeekAtLastError();//cudaGetLastError();
        //// Get error code 			if ( err != cudaSuccess )
        //			{
        //			  fprintf(stderr, "CUDA kernel launch Error: %s \n
        //", cudaGetErrorString(err));
        //			}
        cudaMemcpyAsync(&tid_p[first_part_tmp], &d_tid_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&locx[first_part_tmp], &d_locx[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&locy[first_part_tmp], &d_locy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&locz[first_part_tmp], &d_locz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&h[first_part_tmp], &d_h[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&mass[first_part_tmp], &d_mass[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&x_p[first_part_tmp], &d_x_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&y_p[first_part_tmp], &d_y_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&z_p[first_part_tmp], &d_z_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rho[first_part_tmp], &d_rho[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rho_dh[first_part_tmp], &d_rho_dh[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&wcount[first_part_tmp], &d_wcount[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(
            &wcount_dh[first_part_tmp], &d_wcount_dh[first_part_tmp],
            bundle_size * sizeof(float), cudaMemcpyDeviceToHost, stream[bid]);
        cudaMemcpyAsync(&ux[first_part_tmp], &d_ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&uy[first_part_tmp], &d_uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&uz[first_part_tmp], &d_uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&div_v[first_part_tmp], &d_div_v[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rot_ux[first_part_tmp], &d_rot_ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rot_uy[first_part_tmp], &d_rot_uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rot_uz[first_part_tmp], &d_rot_uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&count_p[first_part_tmp], &d_count_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyDeviceToHost,
                        stream[bid]);
        //		    if(bid%10==0)cudaDeviceSynchronize();
      }

      cudaDeviceSynchronize();
      time_for_memcpys_and_kernel =
          (clock() - time_for_memcpys_and_kernel) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_memcpys_and_kernel is %f\n",
              time_for_memcpys_and_kernel);
      //		err = cudaPeekAtLastError();//cudaGetLastError(); // Get
      //error code 		if ( err != cudaSuccess )
      //		{
      //		  fprintf(stderr, "CUDA kernel launch Error: %s \n ",
      //cudaGetErrorString(err));
      //		}
      //		for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=task_last_part[tid]-task_first_part[tid];
      //		   gpuErrchk(cudaMemcpy(&parts_gpu_CSR[first_part_tmp],
      //&d_parts_gpu_CSR[first_part_tmp], count*sizeof(struct part_gpu),
      //cudaMemcpyDeviceToHost)); 		   err =
      //cudaPeekAtLastError();//cudaGetLastError();        // Get error code 		   if
      //( err != cudaSuccess )
      //		     {
      //		   	   fprintf(stderr, "CUDA Error second memcpy: %s
      //tid is %i\n ", cudaGetErrorString(err), tid);
      //		   	 }
      //		}
      //		for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=ci_list[tid]->hydro.count;
      //		   gpuErrchk(cudaMemcpy(&tid_p[first_part_tmp],
      //&d_tid_p[first_part_tmp], count*sizeof(int), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&locx[first_part_tmp],
      //&d_locx[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&locy[first_part_tmp],
      //&d_locy[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&locz[first_part_tmp],
      //&d_locz[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&h[first_part_tmp],
      //&d_h[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&mass[first_part_tmp],
      //&d_mass[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&x_p[first_part_tmp],
      //&d_x_p[first_part_tmp], count*sizeof(double), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&y_p[first_part_tmp],
      //&d_y_p[first_part_tmp], count*sizeof(double), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&z_p[first_part_tmp],
      //&d_z_p[first_part_tmp], count*sizeof(double), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&rho[first_part_tmp],
      //&d_rho[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&rho_dh[first_part_tmp],
      //&d_rho_dh[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&wcount[first_part_tmp],
      //&d_wcount[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&wcount_dh[first_part_tmp],
      //&d_wcount_dh[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&ux[first_part_tmp],
      //&d_ux[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&uy[first_part_tmp],
      //&d_uy[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&uz[first_part_tmp],
      //&d_uz[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&div_v[first_part_tmp],
      //&d_div_v[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&rot_ux[first_part_tmp],
      //&d_rot_ux[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&rot_uy[first_part_tmp],
      //&d_rot_uy[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&rot_uz[first_part_tmp],
      //&d_rot_uz[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&count_p[first_part_tmp],
      //&d_count_p[first_part_tmp], count*sizeof(int), cudaMemcpyDeviceToHost));
      //	   }
      //		cudaDeviceSynchronize();
      //	    err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	    if ( err != cudaSuccess )
      //		 {
      //		   fprintf(stderr, "CUDA Error second memcpy: %s \n ",
      //cudaGetErrorString(err));
      //		 }
      parts_soa.tid_p = tid_p;
      parts_soa.locx = locx;
      parts_soa.locy = locy;
      parts_soa.locz = locz;
      parts_soa.h = h;
      parts_soa.mass = mass;
      parts_soa.x_p = x_p;
      parts_soa.y_p = y_p;
      parts_soa.z_p = z_p;
      parts_soa.rho = rho;
      parts_soa.rho_dh = rho_dh;
      parts_soa.wcount = wcount;
      parts_soa.wcount_dh = wcount_dh;
      parts_soa.ux = ux;
      parts_soa.uy = uy;
      parts_soa.uz = uz;
      parts_soa.div_v = div_v;
      parts_soa.rot_ux = rot_ux;
      parts_soa.rot_uy = rot_uy;
      parts_soa.rot_uz = rot_uz;
      parts_soa.count_p = count_p;
      //	    for(int p=0; p<count_all_parts; p++){
      //			parts_soa.x_p[p]=x_p[p];
      //parts_soa.y_p[p]=y_p[p]; parts_soa.z_p[p]=z_p[p];
      //			parts_soa.rho[p]=rho[p];
      //parts_soa.rho_dh[p]=rho_dh[p]; parts_soa.wcount[p]=wcount[p];
      //parts_soa.wcount_dh[p]=wcount_dh[p]; 			parts_soa.div_v[p]=div_v[p];
      //parts_soa.rot_ux[p]=ux[p]; parts_soa.tid_p[p]=tid_p[p];
      //	    }
      double t_gpu_copyandkernel_finish = clock();
      fprintf(stderr, "kernel, cpy and re-assignment to struct time is %f\n",
              (t_gpu_copyandkernel_finish - t_gpu_copyandkernel_begin) /
                  CLOCKS_PER_SEC);

      double t_cpu_elapsed = 0.f;
      for (int tid = 0; tid < count_tasks; tid++) {
        //      	switch (tasks[tid]->type) {
        //		     case task_type_self:
        if (tasks[tid]->subtype == task_subtype_density) {
          double t1 = clock();
          runner_doself1_branch_density(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        //				 if (tasks[tid]->subtype ==
        //task_subtype_gradient){ 				 	 double t1=clock();
        //				 	 runner_doself1_branch_gradient(r,
        //ci_list[tid]); 				 	 double t2=clock(); 				 	 t_cpu_elapsed+=t2-t1;
        //				 }
        //				 if (tasks[tid]->subtype ==
        //task_subtype_force){ 				 	 double t1=clock(); 				 	 runner_doself2_branch_force(r,
        //ci_list[tid]); 				 	 double t2=clock(); 				 	 t_cpu_elapsed+=t2-t1;
        //				 }
        //      	}
      }
      double t_cpu_print_start = clock();
      fclose(fpcpu);
      fprintf(stderr, "Printing CPU results\n");
      fpcpu = fopen("./res_CPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_list[tid]->hydro.count;
        for (int p = 0; p < count; p++) {
          float xx = ci_list[tid]->hydro.parts[p].x[0],
                yy = ci_list[tid]->hydro.parts[p].x[1],
                zz = ci_list[tid]->hydro.parts[p].x[2];
          fprintf(fpcpu, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
                  ci_list[tid]->hydro.parts[p].rho,
                  ci_list[tid]->hydro.parts[p].density.rho_dh,
                  ci_list[tid]->hydro.parts[p].density.wcount,
                  ci_list[tid]->hydro.parts[p].density.wcount_dh,
                  ci_list[tid]->hydro.parts[p].viscosity.div_v,
                  ci_list[tid]->hydro.parts[p].density.rot_v[0]);
        }
      }

      fclose(fpcpu);
      double t_cpu_print_end = clock();

      //		for (int tid = 0; tid < count_tasks; ++tid)
      //			gpuErrchk( cudaStreamDestroy(stream[tid]));
      double t_gpu_end = clock();
      fprintf(stderr, "time for cpu is %f\n", t_cpu_elapsed / CLOCKS_PER_SEC);
      fprintf(stderr, "time for gpu is %f\n",
              (t_gpu_end - t_gpu_copyandkernel_begin - t_cpu_elapsed -
               t_cpu_print_end + t_cpu_print_start) /
                  CLOCKS_PER_SEC);
      fprintf(stderr, "finished looping through tasks\n");
      FILE *fp;
      fp = fopen("./res_GPU_density.txt", "w");
      //		for(int tid=0; tid<count_tasks; tid++){
      //			const int
      //count=task_last_part[tid]-task_first_part[tid];
      double t6 = clock();
      for (int p = 0; p < count_all_parts; p++) {
        float xx = parts_soa.x_p[p], yy = parts_soa.y_p[p],
              zz = parts_soa.z_p[p];
        fprintf(fp, "%f %f %f %f %f %f %f %f %f %i %i %f\n", xx, yy, zz,
                parts_soa.rho[p], parts_soa.rho_dh[p], parts_soa.wcount[p],
                parts_soa.wcount_dh[p], parts_soa.div_v[p], parts_soa.rot_ux[p],
                parts_soa.tid_p[p], parts_soa.count_p[p], parts_soa.h[p]);
      }
      fprintf(stderr, "Finished printing GPU results\n");
      //		}
      fclose(fp);
      //      cudaFreeHost(ci_gpu_list);
      cudaProfilerStop();
      cudaDeviceReset();
      exit(0);
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
//////////////////////////////////////Adding pair interactions
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////****same as runner_main_cuda but bunching tasks into bundles. Each stream
///works on one bundle***********/////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void *runner_main_self_pair_sum_loops(void *data) {
  cudaDeviceReset();
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  cudaError_t cu_error;
  cudaMemPool_t memPool;
  gpuErrchk(cudaDeviceGetDefaultMemPool(&memPool, devId));
  int maxmem = 0.9f * prop.totalGlobalMem;
  gpuErrchk(cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold,
                                    (void *)&maxmem));
  int counter = 0;

  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    //	 gpuErrchk( cudaMallocHost((void**)&tasks, n_streams*sizeof(struct
    //task)) ); // Pinned allocation on host Typecasting from void** to struct
    // task **//

    ////////////////////////////////////////////////////////////////////////////////
    ///////This is likely where the problem is!!!!
    ////////////////////////////////////////////////////////////////////////////////
    float time_to_get_tasks = clock();
    struct task *tasks[n_streams];

    while (1) {
      for (int i = 0; i < n_streams; i++) {
        tasks[i] = NULL;
      }
      // for(int tid=0; tid<n_streams; tid++){
      int count_tasks = 0;
      int count_stale = 0;
      ////////////////////////////////////
      /*Grab a bunch of tasks*/
      int count_self = 0, count_sub_self = 0, count_pair = 0, count_self_d = 0,
          count_self_f = 0, count_self_g = 0, count_pair_d = 0,
          count_pair_f = 0, count_pair_g = 0;
      while (1) { //(count_tasks<n_streams){
        // there's a potential bug here. Code hangs if n_streams set to > 256 in
        // cuda_headers.h
        /* If there's no old task, try to get a new one. */
        if (tasks[count_tasks] == NULL) {
          /* Get the task. */
          TIMER_TIC
          tasks[count_tasks] = scheduler_gettask(sched, r->qid, prev);
          TIMER_TOC(timer_gettask);
          if (tasks[count_tasks] != NULL) {
            t = tasks[count_tasks];
            count_tasks++;
            if ((t->type == task_type_self) &&
                (t->subtype == task_subtype_density)) {
              count_self++;
              count_self_d++;
            }
            if ((t->type == task_type_self) &&
                (t->subtype == task_subtype_gradient)) {
              count_self++;
              count_self_g++;
            }
            if ((t->type == task_type_self) &&
                (t->subtype == task_subtype_force)) {
              count_self++;
              count_self_f++;
            }
            if ((t->type == task_type_sub_self) &&
                (t->subtype == task_subtype_density)) {
              count_pair++;
              count_pair_d++;
            }
            if ((t->type == task_type_sub_self) &&
                (t->subtype == task_subtype_gradient)) {
              count_pair++;
              count_pair_g++;
            }
            if ((t->type == task_type_sub_self) &&
                (t->subtype == task_subtype_force)) {
              count_pair++;
              count_pair_f++;
            }
          }
        }

        /* If there's an old task, move onto the next one */
        // else count_tasks++;
        count_stale++;
        if (count_stale >= n_streams) {
          break;
        }
      }

      fprintf(stderr, "Count tasks is %i\n", count_tasks);
      fprintf(stderr, "Count tasks self is %i\n", count_self);
      fprintf(stderr, "Count tasks pair is %i\n", count_pair);
      fprintf(stderr, "Count self d is %i\n", count_self_d);
      fprintf(stderr, "Count self g is %i\n", count_self_g);
      fprintf(stderr, "Count self f is %i\n", count_self_f);
      fprintf(stderr, "Count pair d is %i\n", count_pair_d);
      fprintf(stderr, "Count pair g is %i\n", count_pair_g);
      fprintf(stderr, "Count pair f is %i\n", count_pair_f);
      exit(0);
      ////////////////////////////////////

      ////////////////////////////////////
      /* Get the cells. */
      struct cell **ci_list =
          (struct cell **)malloc(count_tasks * sizeof(struct cell *));
      int count_all_parts = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        //			ci_list[tid]=(struct cell) malloc(sizeof(struct
        //cell*));
        ci_list[tid] = tasks[tid]->ci;
        struct cell *c_tmp = ci_list[tid];
        count_all_parts += ci_list[tid]->hydro.count;
      }

      ////////////////////////////////////

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->type == task_type_pair ||
            tasks[tid]->type == task_type_sub_pair) {
          struct cell *ci_temp = ci_list[tid];
          //		     struct cell *cj_temp = cj;
          double shift[3];
          tasks[tid]->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        } else {
          tasks[tid]->sid = -1;
        }
      }
#endif
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
      i.e., rather than just one r->t, r would have a nested bunch of structs
       r->t[size=nStreams].*/
      r->t = t;
#endif
      time_to_get_tasks = (clock() - time_to_get_tasks) / CLOCKS_PER_SEC;
      fprintf(stderr, "time to get tasks is %f\n", time_to_get_tasks);
      float time_to_malloc_variables = clock();
      //////////////////////////////////////////////
      /*****Initialise GPU suitable copies of variables and malloc pinned
       * memory*****/
      //		struct part_gpu **
      //parts_gpu_list=malloc(count_tasks*sizeof(struct part_gpu *)); 		struct
      //part_gpu * parts_gpu_CSR;
      part_soa parts_soa;
      part_soa d_parts_soa;
      ////////All data contained within parts_gpu_CSR_soa struct. A lot, innit?
      int *d_tid_p;
      long long *d_id;
      double *d_x_p;
      double *d_y_p;
      double *d_z_p;
      float *d_ux;
      float *d_uy;
      float *d_uz;
      float *d_a_hydrox;
      float *d_a_hydroy;
      float *d_a_hydroz;
      float *d_mass;
      float *d_h;
      float *d_u;
      float *d_u_dt;
      float *d_rho;
      float *d_SPH_sum;
      float *d_locx;
      float *d_locy;
      float *d_locz;
      float *d_widthx;
      float *d_widthy;
      float *d_widthz;
      float *d_h_max;
      int *d_count_p;
      float *d_wcount;
      float *d_wcount_dh;
      float *d_rho_dh;
      float *d_rot_ux;
      float *d_rot_uy;
      float *d_rot_uz;
      float *d_div_v;
      float *d_div_v_previous_step;
      float *d_alpha_visc;
      float *d_v_sig;
      float *d_laplace_u;
      float *d_alpha_diff;
      float *d_f;
      float *d_soundspeed;
      float *d_h_dt;
      float *d_balsara;
      float *d_pressure;
      float *d_alpha_visc_max_ngb;
      /* timestep stuff */
      timebin_t *d_time_bin;
      timebin_t *d_wakeup;
      timebin_t *d_min_ngb_time_bin;
      char *d_to_be_synchronized;
      ////////now malloc all this data. Sheesh
      fprintf(stderr, "before mallocing soas\n");
      cudaMalloc((void **)&(d_tid_p), sizeof(int) * count_all_parts);
      cudaMalloc((void **)&(d_id), sizeof(long long) * count_all_parts);
      cudaMalloc((void **)&(d_x_p), sizeof(double) * count_all_parts);
      cudaMalloc((void **)&(d_y_p), sizeof(double) * count_all_parts);
      cudaMalloc((void **)&(d_z_p), sizeof(double) * count_all_parts);
      cudaMalloc((void **)&(d_ux), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_uy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_uz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_a_hydrox), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_a_hydroy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_a_hydroz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_mass), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_h), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_u), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_u_dt), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rho), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_SPH_sum), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_locx), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_locy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_locz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_widthx), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_widthy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_widthz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_h_max), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_count_p), sizeof(int) * count_all_parts);
      cudaMalloc((void **)&(d_wcount), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_wcount_dh), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rho_dh), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rot_ux), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rot_uy), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_rot_uz), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_div_v), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_div_v_previous_step),
                 sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_alpha_visc), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_v_sig), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_laplace_u), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_alpha_diff), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_f), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_soundspeed), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_h_dt), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_balsara), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_pressure), sizeof(float) * count_all_parts);
      cudaMalloc((void **)&(d_alpha_visc_max_ngb),
                 sizeof(float) * count_all_parts);
      /* timestep stuff */
      cudaMalloc((void **)&(d_time_bin), sizeof(timebin_t) * count_all_parts);
      cudaMalloc((void **)&(d_wakeup), sizeof(timebin_t) * count_all_parts);
      cudaMalloc((void **)&(d_min_ngb_time_bin),
                 sizeof(timebin_t) * count_all_parts);
      cudaMalloc((void **)&(d_to_be_synchronized),
                 sizeof(char) * count_all_parts);
      fprintf(stderr, "after mallocing soas\n");
      ///////////Host arrays
      //		int *tid_p=malloc(count_all_parts*sizeof(int));
      int *tid_p;
      cudaMallocHost((void **)&tid_p, count_all_parts * sizeof(int));
      long long *id;
      cudaMallocHost((void **)&id, count_all_parts * sizeof(long long));
      float *mass;
      cudaMallocHost((void **)&mass, count_all_parts * sizeof(float));
      float *h;
      cudaMallocHost((void **)&h, count_all_parts * sizeof(float));
      float *u;
      cudaMallocHost((void **)&u, count_all_parts * sizeof(float));
      float *u_dt;
      cudaMallocHost((void **)&u_dt, count_all_parts * sizeof(float));
      float *rho;
      cudaMallocHost((void **)&rho, count_all_parts * sizeof(float));
      float *SPH_sum;
      cudaMallocHost((void **)&SPH_sum, count_all_parts * sizeof(float));
      //		double x_p[count_all_parts];
      double *x_p;
      cudaMallocHost((void **)&x_p, count_all_parts * sizeof(double));
      double *y_p;
      cudaMallocHost((void **)&y_p, count_all_parts * sizeof(double));
      double *z_p;
      cudaMallocHost((void **)&z_p, count_all_parts * sizeof(double));
      float *ux;
      cudaMallocHost((void **)&ux, count_all_parts * sizeof(float));
      float *uy;
      cudaMallocHost((void **)&uy, count_all_parts * sizeof(float));
      float *uz;
      cudaMallocHost((void **)&uz, count_all_parts * sizeof(float));
      float *a_hydrox;
      cudaMallocHost((void **)&a_hydrox, count_all_parts * sizeof(float));
      float *a_hydroy;
      cudaMallocHost((void **)&a_hydroy, count_all_parts * sizeof(float));
      float *a_hydroz;
      cudaMallocHost((void **)&a_hydroz, count_all_parts * sizeof(float));
      float *locx;
      cudaMallocHost((void **)&locx, count_all_parts * sizeof(float));
      float *locy;
      cudaMallocHost((void **)&locy, count_all_parts * sizeof(float));
      float *locz;
      cudaMallocHost((void **)&locz, count_all_parts * sizeof(float));
      float *widthx;
      cudaMallocHost((void **)&widthx, count_all_parts * sizeof(float));
      float *widthy;
      cudaMallocHost((void **)&widthy, count_all_parts * sizeof(float));
      float *widthz;
      cudaMallocHost((void **)&widthz, count_all_parts * sizeof(float));
      float *h_max;
      cudaMallocHost((void **)&h_max, count_all_parts * sizeof(float));
      int *count_p;
      cudaMallocHost((void **)&count_p, count_all_parts * sizeof(int));
      float *wcount;
      cudaMallocHost((void **)&wcount, count_all_parts * sizeof(float));
      float *wcount_dh;
      cudaMallocHost((void **)&wcount_dh, count_all_parts * sizeof(float));
      float *rho_dh;
      cudaMallocHost((void **)&rho_dh, count_all_parts * sizeof(float));
      float *rot_ux;
      cudaMallocHost((void **)&rot_ux, count_all_parts * sizeof(float));
      float *rot_uy;
      cudaMallocHost((void **)&rot_uy, count_all_parts * sizeof(float));
      float *rot_uz;
      cudaMallocHost((void **)&rot_uz, count_all_parts * sizeof(float));
      float *div_v;
      cudaMallocHost((void **)&div_v, count_all_parts * sizeof(float));
      float *div_v_previous_step;
      cudaMallocHost((void **)&div_v_previous_step,
                     count_all_parts * sizeof(float));
      float *alpha_visc;
      cudaMallocHost((void **)&alpha_visc, count_all_parts * sizeof(float));
      float *v_sig;
      cudaMallocHost((void **)&v_sig, count_all_parts * sizeof(float));
      float *laplace_u;
      cudaMallocHost((void **)&laplace_u, count_all_parts * sizeof(float));
      float *alpha_diff;
      cudaMallocHost((void **)&alpha_diff, count_all_parts * sizeof(float));
      float *f;
      cudaMallocHost((void **)&f, count_all_parts * sizeof(float));
      float *soundspeed;
      cudaMallocHost((void **)&soundspeed, count_all_parts * sizeof(float));
      float *h_dt;
      cudaMallocHost((void **)&h_dt, count_all_parts * sizeof(float));
      float *balsara;
      cudaMallocHost((void **)&balsara, count_all_parts * sizeof(float));
      float *pressure;
      cudaMallocHost((void **)&pressure, count_all_parts * sizeof(float));
      float *alpha_visc_max_ngb;
      cudaMallocHost((void **)&alpha_visc_max_ngb,
                     count_all_parts * sizeof(float));
      /* timestep stuff */
      timebin_t *time_bin;
      cudaMallocHost((void **)&time_bin, count_all_parts * sizeof(timebin_t));
      timebin_t *wakeup;
      cudaMallocHost((void **)&wakeup, count_all_parts * sizeof(timebin_t));
      timebin_t *min_ngb_time_bin;
      cudaMallocHost((void **)&min_ngb_time_bin,
                     count_all_parts * sizeof(timebin_t));
      char *to_be_synchronized;
      cudaMallocHost((void **)&to_be_synchronized,
                     count_all_parts * sizeof(char));
      ///////////////////////////////////////////////////

      /////////////Host arrays
      ////		int *tid_p=malloc(count_all_parts*sizeof(int));
      //		int *tid_p;
      //		gpuErrchk( cudaMallocHost((void**)&tid_p,
      //count_all_parts*sizeof(int)) ); 		long long
      //*id=malloc(count_all_parts*sizeof(long long)); 		float
      //*mass=malloc(count_all_parts*sizeof(float)); 		float
      //*h=malloc(count_all_parts*sizeof(float)); 		float
      //*u=malloc(count_all_parts*sizeof(float)); 		float
      //*u_dt=malloc(count_all_parts*sizeof(float)); 		float
      //*rho=malloc(count_all_parts*sizeof(float)); 		float
      //*SPH_sum=malloc(count_all_parts*sizeof(float));
      ////		double x_p[count_all_parts];
      //		double *x_p=malloc(count_all_parts*sizeof(double));
      //		double *y_p=malloc(count_all_parts*sizeof(double));
      //		double *z_p=malloc(count_all_parts*sizeof(double));
      //		float *ux=malloc(count_all_parts*sizeof(float));
      //		float *uy=malloc(count_all_parts*sizeof(float));
      //		float *uz=malloc(count_all_parts*sizeof(float));
      //		float *a_hydrox=malloc(count_all_parts*sizeof(float));
      //		float *a_hydroy=malloc(count_all_parts*sizeof(float));
      //		float *a_hydroz=malloc(count_all_parts*sizeof(float));
      //		float *locx=malloc(count_all_parts*sizeof(float));
      //		float *locy=malloc(count_all_parts*sizeof(float));
      //		float *locz=malloc(count_all_parts*sizeof(float));
      //		float *widthx=malloc(count_all_parts*sizeof(float));
      //		float *widthy=malloc(count_all_parts*sizeof(float));
      //		float *widthz=malloc(count_all_parts*sizeof(float));
      //		float *h_max=malloc(count_all_parts*sizeof(float));
      //		int *count_p=malloc(count_all_parts*sizeof(int));
      //		float *wcount=malloc(count_all_parts*sizeof(float));
      //		float *wcount_dh=malloc(count_all_parts*sizeof(float));
      //		float *rho_dh=malloc(count_all_parts*sizeof(float));
      //		float *rot_ux=malloc(count_all_parts*sizeof(float));
      //		float *rot_uy=malloc(count_all_parts*sizeof(float));
      //		float *rot_uz=malloc(count_all_parts*sizeof(float));
      //		float *div_v=malloc(count_all_parts*sizeof(float));
      //		float
      //*div_v_previous_step=malloc(count_all_parts*sizeof(float)); 		float
      //*alpha_visc=malloc(count_all_parts*sizeof(float)); 		float
      //*v_sig=malloc(count_all_parts*sizeof(float)); 		float
      //*laplace_u=malloc(count_all_parts*sizeof(float)); 		float
      //*alpha_diff=malloc(count_all_parts*sizeof(float)); 		float
      //*f=malloc(count_all_parts*sizeof(float)); 		float
      //*soundspeed=malloc(count_all_parts*sizeof(float)); 		float
      //*h_dt=malloc(count_all_parts*sizeof(float)); 		float
      //*balsara=malloc(count_all_parts*sizeof(float)); 		float
      //*pressure=malloc(count_all_parts*sizeof(float)); 		float
      //*alpha_visc_max_ngb=malloc(count_all_parts*sizeof(float));
      //		/* timestep stuff */
      //		timebin_t
      //*time_bin=malloc(count_all_parts*sizeof(timebin_t)); 		timebin_t
      //*wakeup=malloc(count_all_parts*sizeof(timebin_t)); 		timebin_t
      //*min_ngb_time_bin=malloc(count_all_parts*sizeof(timebin_t)); 		char
      //*to_be_synchronized=malloc(count_all_parts*sizeof(char));
      /////////////////////////////////////////////////////

      fprintf(stderr, "Size I'm allocating is %i\n",
              count_all_parts * sizeof(struct part_gpu));
      int *task_first_part, *task_last_part;
      int *d_task_first_part, *d_task_last_part;
      cudaMallocHost((void **)&task_first_part,
                     count_tasks * sizeof(int)); // Pinned allocation on host
      cudaMallocHost((void **)&task_last_part,
                     count_tasks * sizeof(int)); // Pinned allocation on host

      int bundle_size = 32;
      int nBundles = (count_tasks + bundle_size - 1) / bundle_size;

      int *bundle_first_part, *bundle_last_part;
      int *d_bundle_first_part, *d_bundle_last_part;
      cudaMallocHost((void **)&bundle_first_part,
                     nBundles * sizeof(int)); // Pinned allocation on host
      cudaMallocHost((void **)&bundle_last_part,
                     nBundles * sizeof(int)); // Pinned allocation on host

      cudaMalloc((void **)&d_task_first_part, count_tasks * sizeof(int));
      cudaMalloc((void **)&d_task_last_part, count_tasks * sizeof(int));
      cudaMalloc((void **)&d_bundle_first_part, nBundles * sizeof(int));
      cudaMalloc((void **)&d_bundle_last_part, nBundles * sizeof(int));
      int total_num_parts = 0;
      ///////Start GPU timer
      //		double t_gpu_start=clock();
      ////////////////////////////////////////////
      time_to_malloc_variables =
          (clock() - time_to_malloc_variables) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_to_malloc_variables %f\n",
              time_to_malloc_variables);
      ////////Malloc on HOST
      float time_for_data_arrangement = clock();
      int first_part_tmp = 0;
      //		fprintf(stderr,"Before populate parts\n");
      int count_prev = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        int count = ci_list[tid]->hydro.count;
        //			fprintf(stderr, "Just before popluator\n");
        ///////////////////////////////////////////////////////////////////////
        //////////////This should probably be done on the GPU??////////////////
        ///////////////////////////////////////////////////////////////////////
        //			struct cell *ci_temp = ci_list[tid];
        populate_parts_list_soa(
            count_all_parts, ci_list[tid], first_part_tmp, count, tid, tid_p,
            id, x_p, y_p, z_p, ux, uy, uz, a_hydrox, a_hydroy, a_hydroz, mass,
            h, u, u_dt, rho, SPH_sum, locx, locy, locz, widthx, widthy, widthz,
            h_max, count_p, wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
            div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
            alpha_diff, f, soundspeed, h_dt, balsara, pressure,
            alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
            to_be_synchronized);
        //			fprintf(stderr, "Just after popluator\n");
        //			populate_parts_list_soa(count_all_parts,
        //ci_list[tid], first_part_tmp, count, &tid_p, &id, &x_p, &y_p, &z_p,
        //&ux, &uy, &uz, &a_hydrox, 									&a_hydroy, &a_hydroz, &mass, &h, &u, &u_dt,
        //&rho, &SPH_sum, &locx, &locy, &locz, 									&widthx, &widthy, &widthz,
        //&h_max, &count_p, &wcount, &wcount_dh, &rho_dh, &rot_u, &rot_v,
        //									&rot_w, &div_v,
        //&div_v_previous_step, &alpha_visc, &v_sig, &laplace_u, &alpha_diff,
        //&f, &soundspeed, 									&h_dt, &balsara, &pressure, &alpha_visc_max_ngb,
        //&time_bin, &wakeup, &min_ngb_time_bin, 									&to_be_synchronized);
        //			fprintf(stderr, "Just after popluator\n");
        //			for(int p=0; p<count; p++){
        //				int p_gid=p+first_part_tmp;
        //				x_p[p_gid]=ci_list[tid]->hydro.parts[p].x[0];
        //				fprintf(stderr,"part is %i x is %f x from cell
        //struct is %f\n",p, x_p[p_gid], ci_list[tid]->hydro.parts[p].x[0]);
        //		    }
        task_first_part[tid] = first_part_tmp;
        task_last_part[tid] = first_part_tmp + count;
        first_part_tmp += count;
        if (tid % bundle_size == 0) {
          int bid = tid / bundle_size;
          bundle_first_part[bid] = task_first_part[tid];
          //				fprintf(stderr, "bundle id %i nBundles %i\n", bid,
          //nBundles);
        }
      }
      for (int bid = 0; bid < nBundles - 1; bid++) {
        bundle_last_part[bid] = bundle_first_part[bid + 1];
      }
      bundle_last_part[nBundles - 1] =
          count_all_parts; //-bundle_last_part[nBundles-1];
                           //		for(int bid=0; bid<nBundles; bid++){
                           //			fprintf(stderr, "count_all_parts %i, bid %i,
      //nBundles %i, first part %i, last part %i\n", count_all_parts, bid,
      //nBundles, bundle_first_part[bid], bundle_last_part[bid]);
      //		}
      //		exit(0);
      //		fprintf(stderr, "GOT HERE\n");
      time_for_data_arrangement =
          (clock() - time_for_data_arrangement) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_data_arrangement %f\n",
              time_for_data_arrangement);
      //		parts_soa.tid_p=tid_p;
      //	    parts_soa.locx=locx;
      //	    parts_soa.locy=locy;
      //	    parts_soa.locz=locz;
      //	    parts_soa.h=h;
      //	    parts_soa.mass=mass;
      //	    parts_soa.x_p=x_p;
      //	    parts_soa.y_p=y_p;
      //	    parts_soa.z_p=z_p;
      //	    parts_soa.rho=rho;
      //	    parts_soa.rho_dh=rho_dh;
      //	    parts_soa.wcount=wcount;
      //	    parts_soa.wcount_dh=wcount_dh;
      //	    parts_soa.ux=ux;
      //	    parts_soa.uy=uy;
      //	    parts_soa.uz=uz;
      //	    parts_soa.div_v=div_v;
      //	    parts_soa.rot_ux=rot_ux;
      //	    parts_soa.rot_uy=rot_uy;
      //	    parts_soa.rot_uz=rot_uz;
      //	    parts_soa.count_p=count_p;
      //	    for(int p=0; p<count_all_parts; p++){
      //	    	parts_soa.count_p[p]=0;
      //	    }
      //	    fprintf(stderr,"got here 2\n");
      //		struct part_gpu * d_parts_gpu_CSR;
      //		cudaMalloc((void**)&d_parts_gpu_CSR,
      //count_all_parts*sizeof(struct part_gpu));

      //////////////////////////////////////////////

      //////Prepare streams and malloc memory on the GPU
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      // create events and streams
      ///////Start GPU timer
      double t_gpu_start = clock();
      ////////////////////////////////////////////

      ////////////////////////////////////////////////
      ///////////Malloc GPU data

      cudaProfilerStart();
      double t_gpu_malloc_begin = clock();

      FILE *fpcpu, *fpgpu;
      fpcpu = fopen("./res_CPU_density.txt", "w");
      double t_cpu_start = clock();
      double t_gpu_elapsed = 0.f;

      // see
      // https://stackoverflow.com/questions/23609770/cuda-double-pointer-memory-copy
      // for details of how this is allocated and copied

      //	   fprintf(stderr,"Malloc d_all_parts\n");
      cudaError_t err;
      //	   for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=ci_list[tid]->hydro.count;
      //
      //		   gpuErrchk(cudaMemcpy(&d_tid_p[first_part_tmp],
      //&tid_p[first_part_tmp], count*sizeof(int), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_locx[first_part_tmp],
      //&locx[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_locy[first_part_tmp],
      //&locy[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_locz[first_part_tmp],
      //&locz[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_h[first_part_tmp],
      //&h[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_mass[first_part_tmp],
      //&mass[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_x_p[first_part_tmp],
      //&x_p[first_part_tmp], count*sizeof(double), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_y_p[first_part_tmp],
      //&y_p[first_part_tmp], count*sizeof(double), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_z_p[first_part_tmp],
      //&z_p[first_part_tmp], count*sizeof(double), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rho[first_part_tmp],
      //&rho[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rho_dh[first_part_tmp],
      //&rho_dh[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_wcount[first_part_tmp],
      //&wcount[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_wcount_dh[first_part_tmp],
      //&wcount_dh[first_part_tmp], count*sizeof(float),
      //cudaMemcpyHostToDevice)); 		   gpuErrchk(cudaMemcpy(&d_ux[first_part_tmp],
      //&ux[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_uy[first_part_tmp],
      //&uy[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_uz[first_part_tmp],
      //&uz[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_div_v[first_part_tmp],
      //&div_v[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rot_ux[first_part_tmp],
      //&rot_ux[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rot_uy[first_part_tmp],
      //&rot_uy[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_rot_uz[first_part_tmp],
      //&rot_uz[first_part_tmp], count*sizeof(float), cudaMemcpyHostToDevice));
      //		   gpuErrchk(cudaMemcpy(&d_count_p[first_part_tmp],
      //&count_p[first_part_tmp], count*sizeof(int), cudaMemcpyHostToDevice));
      //
      //		   err = cudaPeekAtLastError();//cudaGetLastError(); //
      //Get error code 		   if ( err != cudaSuccess )
      //			{
      //			   fprintf(stderr, "CUDA kernel Error first
      //memcpy11: %s \n ", cudaGetErrorString(err));
      //			}
      //	   }

      //	   for(int p=0; p<count_all_parts; p++){
      //		   fprintf(stderr,"x position is
      //%f\n",parts_soa.x_p[p]); 		   fprintf(stderr,"tid is
      //%i\n",parts_soa.tid_p[p]);
      //	   }

      //	   FILE *fptest = fopen("./res_GPU_density.txt", "w");
      //	   	   for(int p=0; p<count_all_parts; p++){
      //	   			float xx=parts_soa.x_p[p],
      //yy=parts_soa.y_p[p], zz=parts_soa.z_p[p];
      //	   //			fprintf(stderr,"x is %f y is %f z is
      //%f\n", xx, yy, zz); 	   			fprintf(fptest, "%f %f %f %f %f %f %f %f %f %i\n",
      //xx, yy, zz, parts_soa.rho[p], 	   					parts_soa.rho_dh[p], parts_soa.wcount[p],
      //parts_soa.wcount_dh[p], 	   					parts_soa.div_v[p], parts_soa.rot_ux[p],
      //parts_soa.tid_p[p]);
      //	   	   }
      //	   	   fprintf(stderr,"Finished printing GPU results\n");
      //	   	   //		}
      //	   	   fclose(fptest);
      //	   	FILE *fptestcpu = fopen("./res_CPU_density.txt", "w");
      //	   		   for(int tid=0; tid<count_tasks; tid++){
      //	   			int count = ci_list[tid]->hydro.count;
      //	   			for (int p=0; p<count; p++){
      //	   				float
      //xx=ci_list[tid]->hydro.parts[p].x[0],
      //yy=ci_list[tid]->hydro.parts[p].x[1],
      //zz=ci_list[tid]->hydro.parts[p].x[2]; 					fprintf(fptestcpu, "%f %f %f %f %f
      //%f %f %f %f\n", xx, yy, zz, ci_list[tid]->hydro.parts[p].rho,
      //					ci_list[tid]->hydro.parts[p].density.rho_dh,
      //ci_list[tid]->hydro.parts[p].density.wcount,
      //					ci_list[tid]->hydro.parts[p].density.wcount_dh,
      //ci_list[tid]->hydro.parts[p].viscosity.div_v 					,
      //ci_list[tid]->hydro.parts[p].density.rot_v[0]);
      //	   			}
      //	   		   }
      //	   		   fprintf(stderr,"Finished printing GPU
      //results\n");
      //	   		   //		}
      //	   		   fclose(fptestcpu);
      //	   exit(0);
      float time_for_extraneous_memcpys = clock();
      parts_soa.tid_p = d_tid_p;
      parts_soa.locx = d_locx;
      parts_soa.locy = d_locy;
      parts_soa.locz = d_locz;
      parts_soa.h = d_h;
      parts_soa.mass = d_mass;
      parts_soa.x_p = d_x_p;
      parts_soa.y_p = d_y_p;
      parts_soa.z_p = d_z_p;
      parts_soa.rho = d_rho;
      parts_soa.rho_dh = d_rho_dh;
      parts_soa.wcount = d_wcount;
      parts_soa.wcount_dh = d_wcount_dh;
      parts_soa.ux = d_ux;
      parts_soa.uy = d_uy;
      parts_soa.uz = d_uz;
      parts_soa.div_v = d_div_v;
      parts_soa.rot_ux = d_rot_ux;
      parts_soa.rot_uy = d_rot_uy;
      parts_soa.rot_uz = d_rot_uz;
      parts_soa.count_p = d_count_p;

      //	   fprintf(stderr,"number of parts %i, first part in last task
      //%i, last part in last task %i\n", count_all_parts,
      //task_first_part[count_tasks-1], task_last_part[count_tasks-1]);
      //	   fprintf(stderr,"number of parts %i, first part in fist bundle
      //%i, last part in last bundle %i\n", count_all_parts,
      //bundle_first_part[nBundles-1], bundle_last_part[nBundles-1]);

      cudaMemcpy(d_task_first_part, task_first_part, count_tasks * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy1: %s
      //\n ", cudaGetErrorString(err));
      //		}
      cudaMemcpy(d_task_last_part, task_last_part, count_tasks * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy2: %s
      //\n ", cudaGetErrorString(err));
      //		}
      cudaMemcpy(d_bundle_first_part, bundle_first_part, nBundles * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy3: %s
      //\n ", cudaGetErrorString(err));
      //		}
      cudaMemcpy(d_bundle_last_part, bundle_last_part, nBundles * sizeof(int),
                 cudaMemcpyHostToDevice);
      //	   err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	   if ( err != cudaSuccess )
      //		{
      //		   fprintf(stderr, "CUDA kernel Error first memcpy4: %s
      //\n ", cudaGetErrorString(err));
      //		}
      time_for_extraneous_memcpys =
          (clock() - time_for_extraneous_memcpys) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_extraneous_memcpys is %f\n",
              time_for_extraneous_memcpys);
      // For another possible solution, see
      // https://forums.developer.nvidia.com/t/how-do-i-pass-a-double-pointers-array-to-the-device-im-getting-cudaerrorillegaladdress/72518/2

      int size_parts_mmcpy = 0;
      float time_for_stream_creation = clock();
      cudaEvent_t startEvent, stopEvent, dummyEvent;

      cudaEventCreate(&startEvent);
      cudaEventCreate(&stopEvent);
      cudaEventCreate(&dummyEvent);
      //		fprintf(stderr,"got here in the end\n");
      double t_gpu_malloc_finish = clock();
      //		fprintf(stderr,"malloc time is
      //%f\n",(t_gpu_malloc_finish-t_gpu_malloc_begin)/CLOCKS_PER_SEC);

      cudaStream_t stream[nBundles];
      int tasksperbundle = (count_tasks + nBundles - 1) / nBundles;
      //		fprintf(stderr,"count_tasks/nBundles %i, tasksperbundle
      //%i\n", count_tasks/nBundles, tasksperbundle); 		exit(0);
      for (int i = 0; i < nBundles; ++i)
        //			gpuErrchk( cudaStreamCreate(&stream[i]) );
        cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking);
      //		err = cudaPeekAtLastError();//cudaGetLastError(); // Get
      //error code 		if ( err != cudaSuccess )
      //		{
      //		  fprintf(stderr, "CUDA Stream Creation Error: %s \n ",
      //cudaGetErrorString(err));
      //		}
      time_for_stream_creation =
          (clock() - time_for_stream_creation) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_stream_creation is %f\n",
              time_for_stream_creation);
      t_gpu_elapsed = 0.f;
      float time_for_memcpys_and_kernel = clock();
      //		for(int bid=0; bid<nBundles; bid++){
      //		   int first_part_tmp=bundle_first_part[bid];
      //		   int bundle_size=bundle_last_part[bid]-first_part_tmp;
      //
      //		   gpuErrchk(cudaMemcpyAsync(&d_tid_p[first_part_tmp],
      //&tid_p[first_part_tmp], bundle_size*sizeof(int), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_locx[first_part_tmp],
      //&locx[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_locy[first_part_tmp],
      //&locy[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_locz[first_part_tmp],
      //&locz[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_h[first_part_tmp],
      //&h[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_mass[first_part_tmp],
      //&mass[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_x_p[first_part_tmp],
      //&x_p[first_part_tmp], bundle_size*sizeof(double),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_y_p[first_part_tmp],
      //&y_p[first_part_tmp], bundle_size*sizeof(double),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_z_p[first_part_tmp],
      //&z_p[first_part_tmp], bundle_size*sizeof(double),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rho[first_part_tmp],
      //&rho[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_rho_dh[first_part_tmp],
      //&rho_dh[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_wcount[first_part_tmp],
      //&wcount[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_wcount_dh[first_part_tmp],
      //&wcount_dh[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_ux[first_part_tmp],
      //&ux[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_uy[first_part_tmp],
      //&uy[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_uz[first_part_tmp],
      //&uz[first_part_tmp], bundle_size*sizeof(float), cudaMemcpyHostToDevice,
      //stream[bid])); 		   gpuErrchk(cudaMemcpyAsync(&d_div_v[first_part_tmp],
      //&div_v[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rot_ux[first_part_tmp],
      //&rot_ux[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rot_uy[first_part_tmp],
      //&rot_uy[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_rot_uz[first_part_tmp],
      //&rot_uz[first_part_tmp], bundle_size*sizeof(float),
      //cudaMemcpyHostToDevice, stream[bid]));
      //		   gpuErrchk(cudaMemcpyAsync(&d_count_p[first_part_tmp],
      //&count_p[first_part_tmp], bundle_size*sizeof(int),
      //cudaMemcpyHostToDevice, stream[bid]));
      //
      //		   err = cudaPeekAtLastError();//cudaGetLastError(); //
      //Get error code 		   if ( err != cudaSuccess )
      //			{
      //			   fprintf(stderr, "CUDA kernel Error first
      //memcpy11: %s \n ", cudaGetErrorString(err));
      //			}
      //		}
      int shared_memory_used = 0;
      double t_gpu_copyandkernel_begin = clock();
      for (int bid = 0; bid < nBundles; bid++) {
        int max_parts = 0;
        int parts_in_bundle = 0;
        for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
             tid++) {
          if (tid < count_tasks) {
            int count = task_last_part[tid] - task_first_part[tid];
            //					if(tid<count_tasks-1)
            parts_in_bundle += count;
            max_parts = max(max_parts, count);
          }
        }

        int first_part_tmp = bundle_first_part[bid];
        int bundle_size = bundle_last_part[bid] - first_part_tmp;
        cudaMemcpyAsync(&d_tid_p[first_part_tmp], &tid_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_locx[first_part_tmp], &locx[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_locy[first_part_tmp], &locy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_locz[first_part_tmp], &locz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_h[first_part_tmp], &h[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_mass[first_part_tmp], &mass[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_x_p[first_part_tmp], &x_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_y_p[first_part_tmp], &y_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_z_p[first_part_tmp], &z_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rho[first_part_tmp], &rho[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rho_dh[first_part_tmp], &rho_dh[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_wcount[first_part_tmp], &wcount[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_wcount_dh[first_part_tmp],
                        &wcount_dh[first_part_tmp], bundle_size * sizeof(float),
                        cudaMemcpyHostToDevice, stream[bid]);
        cudaMemcpyAsync(&d_ux[first_part_tmp], &ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_uy[first_part_tmp], &uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_uz[first_part_tmp], &uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_div_v[first_part_tmp], &div_v[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rot_ux[first_part_tmp], &rot_ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rot_uy[first_part_tmp], &rot_uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_rot_uz[first_part_tmp], &rot_uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyHostToDevice,
                        stream[bid]);
        cudaMemcpyAsync(&d_count_p[first_part_tmp], &count_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyHostToDevice,
                        stream[bid]);
        int tid = 0;
        int offset = bid * tasksperbundle;
        int sizecells = 0, sizeparts = 0;

        sizeparts = tasksperbundle * sizeof(struct part_gpu *);
        int tasks_left = tasksperbundle;
        if (bid == nBundles - 1) {
          tasks_left = count_tasks - (nBundles - 1) * tasksperbundle;
          sizeparts = tasks_left * sizeof(struct part_gpu *);
        }
        int numBlocks_y = tasks_left;
        int block_size = BLOCK_SIZE;
        //			fprintf(stderr, "block size in runner_main is %i\n",
        //block_size);
        int numBlocks_x = (max_parts + block_size - 1) / block_size;
        int numBlocks = (parts_in_bundle + block_size - 1) / block_size;
        //			fprintf(stderr, "numblocks %i, parts in bundle
        //%i\n", numBlocks, parts_in_bundle); 			fprintf(stderr, "bid is %i, max
        //parts/cell in bundle %i Ntasks is %i tasksperbundle is %i\n NBlocks x
        //is %i NBlocksy is %i numBundles is %i\n", bid, max_parts, count_tasks,
        //tasksperbundle, numBlocks_x, numBlocks_y, nBundles);
        const char *loop_type = "density";
        int shared_memory_block =
            BLOCK_SIZE * (8 * sizeof(float) + sizeof(float *));
        int bundle_part_0 = bundle_first_part[bid];
        int bundle_first_task = tid_p[bundle_part_0];
        //			fprintf(stderr, "size required per block is %i,
        //total requested for kernel is %i\n", shared_memory_block,
        //shared_memory_block*numBlocks);
        launch_cuda_kernel_bundles_revised_soa(
            parts_soa, d_task_first_part, d_task_last_part, d_bundle_first_part,
            d_bundle_last_part, numBlocks, d_a, d_H, loop_type, stream[bid],
            bid, block_size, count_tasks, tasksperbundle, numBlocks_x,
            numBlocks_y, tid, offset, bundle_first_task, max_parts);
        //			cudaDeviceSynchronize();
        //			err = cudaPeekAtLastError();//cudaGetLastError();
        //// Get error code 			if ( err != cudaSuccess )
        //			{
        //			  fprintf(stderr, "CUDA kernel launch Error: %s \n
        //", cudaGetErrorString(err));
        //			}
        cudaMemcpyAsync(&tid_p[first_part_tmp], &d_tid_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&locx[first_part_tmp], &d_locx[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&locy[first_part_tmp], &d_locy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&locz[first_part_tmp], &d_locz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&h[first_part_tmp], &d_h[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&mass[first_part_tmp], &d_mass[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&x_p[first_part_tmp], &d_x_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&y_p[first_part_tmp], &d_y_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&z_p[first_part_tmp], &d_z_p[first_part_tmp],
                        bundle_size * sizeof(double), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rho[first_part_tmp], &d_rho[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rho_dh[first_part_tmp], &d_rho_dh[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&wcount[first_part_tmp], &d_wcount[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(
            &wcount_dh[first_part_tmp], &d_wcount_dh[first_part_tmp],
            bundle_size * sizeof(float), cudaMemcpyDeviceToHost, stream[bid]);
        cudaMemcpyAsync(&ux[first_part_tmp], &d_ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&uy[first_part_tmp], &d_uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&uz[first_part_tmp], &d_uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&div_v[first_part_tmp], &d_div_v[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rot_ux[first_part_tmp], &d_rot_ux[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rot_uy[first_part_tmp], &d_rot_uy[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&rot_uz[first_part_tmp], &d_rot_uz[first_part_tmp],
                        bundle_size * sizeof(float), cudaMemcpyDeviceToHost,
                        stream[bid]);
        cudaMemcpyAsync(&count_p[first_part_tmp], &d_count_p[first_part_tmp],
                        bundle_size * sizeof(int), cudaMemcpyDeviceToHost,
                        stream[bid]);
        //		    if(bid%10==0)cudaDeviceSynchronize();
      }

      cudaDeviceSynchronize();
      time_for_memcpys_and_kernel =
          (clock() - time_for_memcpys_and_kernel) / CLOCKS_PER_SEC;
      fprintf(stderr, "time_for_memcpys_and_kernel is %f\n",
              time_for_memcpys_and_kernel);
      //		err = cudaPeekAtLastError();//cudaGetLastError(); // Get
      //error code 		if ( err != cudaSuccess )
      //		{
      //		  fprintf(stderr, "CUDA kernel launch Error: %s \n ",
      //cudaGetErrorString(err));
      //		}
      //		for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=task_last_part[tid]-task_first_part[tid];
      //		   gpuErrchk(cudaMemcpy(&parts_gpu_CSR[first_part_tmp],
      //&d_parts_gpu_CSR[first_part_tmp], count*sizeof(struct part_gpu),
      //cudaMemcpyDeviceToHost)); 		   err =
      //cudaPeekAtLastError();//cudaGetLastError();        // Get error code 		   if
      //( err != cudaSuccess )
      //		     {
      //		   	   fprintf(stderr, "CUDA Error second memcpy: %s
      //tid is %i\n ", cudaGetErrorString(err), tid);
      //		   	 }
      //		}
      //		for(int tid=0; tid<count_tasks; tid++){
      //		   int first_part_tmp=task_first_part[tid];
      //		   int count=ci_list[tid]->hydro.count;
      //		   gpuErrchk(cudaMemcpy(&tid_p[first_part_tmp],
      //&d_tid_p[first_part_tmp], count*sizeof(int), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&locx[first_part_tmp],
      //&d_locx[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&locy[first_part_tmp],
      //&d_locy[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&locz[first_part_tmp],
      //&d_locz[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&h[first_part_tmp],
      //&d_h[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&mass[first_part_tmp],
      //&d_mass[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&x_p[first_part_tmp],
      //&d_x_p[first_part_tmp], count*sizeof(double), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&y_p[first_part_tmp],
      //&d_y_p[first_part_tmp], count*sizeof(double), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&z_p[first_part_tmp],
      //&d_z_p[first_part_tmp], count*sizeof(double), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&rho[first_part_tmp],
      //&d_rho[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&rho_dh[first_part_tmp],
      //&d_rho_dh[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&wcount[first_part_tmp],
      //&d_wcount[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&wcount_dh[first_part_tmp],
      //&d_wcount_dh[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&ux[first_part_tmp],
      //&d_ux[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&uy[first_part_tmp],
      //&d_uy[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&uz[first_part_tmp],
      //&d_uz[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&div_v[first_part_tmp],
      //&d_div_v[first_part_tmp], count*sizeof(float), cudaMemcpyDeviceToHost));
      //		   gpuErrchk(cudaMemcpy(&rot_ux[first_part_tmp],
      //&d_rot_ux[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&rot_uy[first_part_tmp],
      //&d_rot_uy[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&rot_uz[first_part_tmp],
      //&d_rot_uz[first_part_tmp], count*sizeof(float),
      //cudaMemcpyDeviceToHost)); 		   gpuErrchk(cudaMemcpy(&count_p[first_part_tmp],
      //&d_count_p[first_part_tmp], count*sizeof(int), cudaMemcpyDeviceToHost));
      //	   }
      //		cudaDeviceSynchronize();
      //	    err = cudaPeekAtLastError();//cudaGetLastError();        //
      //Get error code 	    if ( err != cudaSuccess )
      //		 {
      //		   fprintf(stderr, "CUDA Error second memcpy: %s \n ",
      //cudaGetErrorString(err));
      //		 }
      parts_soa.tid_p = tid_p;
      parts_soa.locx = locx;
      parts_soa.locy = locy;
      parts_soa.locz = locz;
      parts_soa.h = h;
      parts_soa.mass = mass;
      parts_soa.x_p = x_p;
      parts_soa.y_p = y_p;
      parts_soa.z_p = z_p;
      parts_soa.rho = rho;
      parts_soa.rho_dh = rho_dh;
      parts_soa.wcount = wcount;
      parts_soa.wcount_dh = wcount_dh;
      parts_soa.ux = ux;
      parts_soa.uy = uy;
      parts_soa.uz = uz;
      parts_soa.div_v = div_v;
      parts_soa.rot_ux = rot_ux;
      parts_soa.rot_uy = rot_uy;
      parts_soa.rot_uz = rot_uz;
      parts_soa.count_p = count_p;
      //	    for(int p=0; p<count_all_parts; p++){
      //			parts_soa.x_p[p]=x_p[p];
      //parts_soa.y_p[p]=y_p[p]; parts_soa.z_p[p]=z_p[p];
      //			parts_soa.rho[p]=rho[p];
      //parts_soa.rho_dh[p]=rho_dh[p]; parts_soa.wcount[p]=wcount[p];
      //parts_soa.wcount_dh[p]=wcount_dh[p]; 			parts_soa.div_v[p]=div_v[p];
      //parts_soa.rot_ux[p]=ux[p]; parts_soa.tid_p[p]=tid_p[p];
      //	    }
      double t_gpu_copyandkernel_finish = clock();
      fprintf(stderr, "kernel, cpy and re-assignment to struct time is %f\n",
              (t_gpu_copyandkernel_finish - t_gpu_copyandkernel_begin) /
                  CLOCKS_PER_SEC);

      double t_cpu_elapsed = 0.f;
      for (int tid = 0; tid < count_tasks; tid++) {
        //      	switch (tasks[tid]->type) {
        //		     case task_type_self:
        if (tasks[tid]->subtype == task_subtype_density) {
          double t1 = clock();
          runner_doself1_branch_density(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        //				 if (tasks[tid]->subtype ==
        //task_subtype_gradient){ 				 	 double t1=clock();
        //				 	 runner_doself1_branch_gradient(r,
        //ci_list[tid]); 				 	 double t2=clock(); 				 	 t_cpu_elapsed+=t2-t1;
        //				 }
        //				 if (tasks[tid]->subtype ==
        //task_subtype_force){ 				 	 double t1=clock(); 				 	 runner_doself2_branch_force(r,
        //ci_list[tid]); 				 	 double t2=clock(); 				 	 t_cpu_elapsed+=t2-t1;
        //				 }
        //      	}
      }
      double t_cpu_print_start = clock();
      fclose(fpcpu);
      fprintf(stderr, "Printing CPU results\n");
      fpcpu = fopen("./res_CPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = ci_list[tid]->hydro.count;
        for (int p = 0; p < count; p++) {
          float xx = ci_list[tid]->hydro.parts[p].x[0],
                yy = ci_list[tid]->hydro.parts[p].x[1],
                zz = ci_list[tid]->hydro.parts[p].x[2];
          fprintf(fpcpu, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
                  ci_list[tid]->hydro.parts[p].rho,
                  ci_list[tid]->hydro.parts[p].density.rho_dh,
                  ci_list[tid]->hydro.parts[p].density.wcount,
                  ci_list[tid]->hydro.parts[p].density.wcount_dh,
                  ci_list[tid]->hydro.parts[p].viscosity.div_v,
                  ci_list[tid]->hydro.parts[p].density.rot_v[0]);
        }
      }

      fclose(fpcpu);
      double t_cpu_print_end = clock();

      //		for (int tid = 0; tid < count_tasks; ++tid)
      //			gpuErrchk( cudaStreamDestroy(stream[tid]));
      double t_gpu_end = clock();
      fprintf(stderr, "time for cpu is %f\n", t_cpu_elapsed / CLOCKS_PER_SEC);
      fprintf(stderr, "time for gpu is %f\n",
              (t_gpu_end - t_gpu_copyandkernel_begin - t_cpu_elapsed -
               t_cpu_print_end + t_cpu_print_start) /
                  CLOCKS_PER_SEC);
      fprintf(stderr, "finished looping through tasks\n");
      FILE *fp;
      fp = fopen("./res_GPU_density.txt", "w");
      //		for(int tid=0; tid<count_tasks; tid++){
      //			const int
      //count=task_last_part[tid]-task_first_part[tid];
      double t6 = clock();
      for (int p = 0; p < count_all_parts; p++) {
        float xx = parts_soa.x_p[p], yy = parts_soa.y_p[p],
              zz = parts_soa.z_p[p];
        fprintf(fp, "%f %f %f %f %f %f %f %f %f %i %i %f\n", xx, yy, zz,
                parts_soa.rho[p], parts_soa.rho_dh[p], parts_soa.wcount[p],
                parts_soa.wcount_dh[p], parts_soa.div_v[p], parts_soa.rot_ux[p],
                parts_soa.tid_p[p], parts_soa.count_p[p], parts_soa.h[p]);
      }
      fprintf(stderr, "Finished printing GPU results\n");
      //		}
      fclose(fp);
      //      cudaFreeHost(ci_gpu_list);
      cudaProfilerStop();
      cudaDeviceReset();
      exit(0);
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////****same as runner_main_cuda but bunching tasks into bundles. Each stream
///works on one bundle***********/////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void *runner_main_cuda_bundles_unified_mem(void *data) {
  cudaDeviceReset();
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  //#ifdef WITH_CUDA
  //  Initialise_GPU();
  //#endif
  /* Main loop. */
  int devId = 0; // find and print device name
  struct cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  cudaError_t cu_error;
  cudaMemPool_t memPool;
  gpuErrchk(cudaDeviceGetDefaultMemPool(&memPool, devId));
  int maxmem = 0.9f * prop.totalGlobalMem;
  gpuErrchk(cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold,
                                    (void *)&maxmem));
  int counter = 0;
  int n_tasks_per_thread = 32;

  while (1) {
    /* Wait at the barrier. */
    engine_barrier(e);
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    struct task **tasks;
    gpuErrchk(cudaMallocHost(
        (void **)&tasks,
        n_streams * sizeof(struct task))); // Pinned allocation on host
    for (int i = 0; i < n_streams; i++) {
      tasks[i] = NULL;
    }

    while (1) {
      for (int i = 0; i < n_streams; i++) {
        tasks[i] = NULL;
      }
      // for(int tid=0; tid<n_streams; tid++){
      int count_tasks = 0;
      int count_stale = 0;
      ////////////////////////////////////
      /*Grab a bunch of tasks*/
      while (1) { //(count_tasks<n_streams){
        // there's a potential bug here. Code hangs if n_streams set to > 256 in
        // cuda_headers.h
        /* If there's no old task, try to get a new one. */
        if (tasks[count_tasks] == NULL) {
          /* Get the task. */
          TIMER_TIC
          tasks[count_tasks] = scheduler_gettask(sched, r->qid, prev);
          TIMER_TOC(timer_gettask);
          if (tasks[count_tasks] != NULL) {
            count_tasks++;
          }
        }
        /* If there's an old task, move onto the next one */
        // else count_tasks++;
        count_stale++;
        if (count_stale >= n_streams) {
          break;
        }
      }
      fprintf(stderr, "Count tasks is %i\n", count_tasks);
      ////////////////////////////////////

      ////////////////////////////////////
      /*Get cell data for each task*/
      for (int tid = 1; tid < count_tasks; tid++) {
        struct cell *c_gpu = tasks[tid]->ci;
      }

      /* Get the cells. */
      struct cell **ci_list;
      struct cell **ci_list_mgd;

      gpuErrchk(cudaMallocHost(
          (void **)&ci_list,
          count_tasks * sizeof(struct cell *))); // Pinned allocation on host

      for (int tid = 0; tid < count_tasks; tid++) {
        ci_list[tid] = tasks[tid]->ci;
      }

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->type == task_type_pair ||
            tasks[tid]->type == task_type_sub_pair) {
          struct cell *ci_temp = ci_list[tid];
          //		     struct cell *cj_temp = cj;
          double shift[3];
          tasks[tid]->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        } else {
          tasks[tid]->sid = -1;
        }
      }
#endif
#ifdef SWIFT_DEBUG_CHECKS
      /* Check that we haven't scheduled an inactive task */
      t->ti_run = e->ti_current;
      /* Store the task that will be running (for debugging only) */
      /*COULD THE t SUBSTRUCT BE USED TO GROUP A BUNCH OF TASKS FOR A RUNNER r?
      i.e., rather than just one r->t, r would have a nested bunch of structs
       r->t[size=nStreams].*/
      r->t = t;
#endif

      //////////////////////////////////////////////
      /*****Initialise GPU suitable copies of variables and malloc pinned
       * memory*****/
      double t_gpu_start = clock();
      double data_rearrangement_begin = clock();
      struct part_gpu **parts_gpu_list;
      cudaMallocManaged((void **)&parts_gpu_list,
                        count_tasks * sizeof(struct part_gpu *),
                        cudaMemAttachGlobal);

      int total_num_parts = 0;
      ///////Start GPU timer
      ////////////////////////////////////////////

      ////////Malloc on HOST
      fprintf(stderr, "Got to before loop\n");
      int count_parts_total = 0;
      for (int tid = 0; tid < count_tasks; tid++) {
        int count = ci_list[tid]->hydro.count;
        //			fprintf(stderr,"Postion is
        //%f\n",ci_list[tid]->hydro.parts[0].x[0]);
        count_parts_total += ci_list[tid]->hydro.count;
        cudaMallocManaged((void **)&parts_gpu_list[tid],
                          count * sizeof(struct part_gpu), cudaMemAttachGlobal);
        //			populate_parts_list(ci_list[tid],
        //parts_gpu_list[tid]);
        total_num_parts = count_parts_total;
      }
      double rearrangement_malloc_end = clock();
      for (int tid = 0; tid < count_tasks; tid++) {
        populate_parts_list(ci_list[tid], parts_gpu_list[tid]);
      }
      //////////////////////////////////////////////
      //////////////////////////////////////////////
      double data_rearrangement_end = clock();
      printf("rearrangement malloc takes %f s\n",
             (rearrangement_malloc_end - data_rearrangement_begin) /
                 CLOCKS_PER_SEC);
      printf("data_rearrangement takes %f s\n",
             (data_rearrangement_end - data_rearrangement_begin) /
                 CLOCKS_PER_SEC);
      //////Prepare streams and malloc memory on the GPU
      float d_a = e->cosmology->a;
      float d_H = e->cosmology->H;
      /* Different types of tasks... */
      // create events and streams
      ///////Start GPU timer
      t_gpu_start = clock();
      ////////////////////////////////////////////

      ////////////////////////////////////////////////
      ///////////Malloc GPU data

      cudaProfilerStart();
      double t_gpu_malloc_begin = clock();

      FILE *fpcpu, *fpgpu;
      fpcpu = fopen("./res_CPU_density.txt", "w");
      double t_cpu_start = clock();
      double t_gpu_elapsed = 0.f;

      int size_parts_mmcpy = 0;

      double t_gpu_malloc_finish = clock();
      fprintf(stderr, "malloc time is %f\n",
              (t_gpu_malloc_finish - t_gpu_malloc_begin) / CLOCKS_PER_SEC);

      int bundle_size = 8;
      int nBundles = (count_tasks + bundle_size - 1) / bundle_size;
      // create events and streams
      double stream_setup_begin = clock();
      cudaStream_t stream[nBundles];
      int tasksperbundle = (count_tasks + nBundles - 1) / nBundles;
      fprintf(stderr, "n bundles is %i, tasks per bundle is %i\n", nBundles,
              tasksperbundle);

      cudaEvent_t startEvent, stopEvent, dummyEvent;
      gpuErrchk(cudaEventCreate(&startEvent));
      gpuErrchk(cudaEventCreate(&stopEvent));
      gpuErrchk(cudaEventCreate(&dummyEvent));
      for (int i = 0; i < nBundles; ++i)
        gpuErrchk(cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking));
      //			gpuErrchk( cudaStreamCreate(&stream[i]) );
      double stream_setup_end = clock();
      printf("stream setup takes %f seconds\n",
             (stream_setup_end - stream_setup_begin) / CLOCKS_PER_SEC);
      t_gpu_elapsed = 0.f;
      int max_parts = 0;
      double gpu_time_start_no_setup = clock();
      double t_gpu_copyandkernel_begin = clock();

      for (int bid = 0; bid < nBundles; bid++) {
        max_parts = 0;
        int countperbundle = 0;
        for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
             tid++) {
          if (tid < count_tasks) {
            int count = parts_gpu_list[tid][0].count;
            countperbundle += count;
            max_parts = max(max_parts, count);
          }
        }
        //			fprintf(stderr,"bundle id %i, max_parts is %i\n",
        //bid, max_parts);
        int tid = 0;
        int offset = bid * tasksperbundle;
        int sizecells = 0, sizeparts = 0;

        int numBlocks_y = tasksperbundle;
        sizecells = tasksperbundle * sizeof(struct cell_gpu_flat);
        sizeparts = tasksperbundle * sizeof(struct part_gpu);
        int tasks_left = tasksperbundle;
        if (bid == nBundles - 1) {
          tasks_left = count_tasks - (nBundles - 1) * tasksperbundle;
          sizecells = tasks_left * sizeof(struct cell_gpu_flat);
          sizeparts = tasks_left * sizeof(struct part_gpu);
          numBlocks_y = tasks_left;
        }
        sizeparts = countperbundle * sizeof(struct part_gpu) +
                    tasks_left * sizeof(struct part_gpu *);

        int block_size = 32;
        int numBlocks; //=(tasks_left + block_size - 1) / block_size;
                       ////////////////////////////////////////////////
        int numBlocks_x = (max_parts + block_size - 1) / block_size;
        const char *loop_type = "density";
        int device = -1, host = -1;
        cudaGetDevice(&device);
        //			cudaMemPrefetchAsync(parts_gpu_list[bid*bundle_size],
        //sizeparts, device, stream[bid]);
        mgd_mem_cuda_kernel_bundles(parts_gpu_list, numBlocks, d_a, d_H,
                                    loop_type, stream, bid, block_size,
                                    count_tasks, tasksperbundle, numBlocks_x,
                                    numBlocks_y, tid, offset);
        //			cudaMemPrefetchAsync(parts_gpu_list[bid*bundle_size],
        //sizeparts, -1, stream[bid]);
        //////////////////////////////////////////////////////
      }
      cudaDeviceSynchronize();
      double t_gpu_copyandkernel_finish = clock();
      double gpu_time_end_no_setup = clock();
      printf("GPU time no setup is %f\n",
             (gpu_time_end_no_setup - gpu_time_start_no_setup) /
                 CLOCKS_PER_SEC);
      fprintf(stderr, "kernel and cpy time is %f\n",
              (t_gpu_copyandkernel_finish - t_gpu_copyandkernel_begin) /
                  CLOCKS_PER_SEC);
      cu_error =
          cudaPeekAtLastError(); // cudaGetLastError();        // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(cu_error));
      }

      double t_cpu_elapsed = 0.f;
      for (int tid = 0; tid < count_tasks; tid++) {
        if (tasks[tid]->subtype == task_subtype_density) {
          double t1 = clock();
          runner_doself1_branch_density(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        if (tasks[tid]->subtype == task_subtype_gradient) {
          double t1 = clock();
          runner_doself1_branch_gradient(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        if (tasks[tid]->subtype == task_subtype_force) {
          double t1 = clock();
          runner_doself2_branch_force(r, ci_list[tid]);
          double t2 = clock();
          t_cpu_elapsed += t2 - t1;
        }
        //      	}
      }
      double t_cpu_print_start = clock();
      fclose(fpcpu);
      fpcpu = fopen("./res_CPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = parts_gpu_list[tid][0].count;
        for (int p = 0; p < count; p++) {
          float xx = ci_list[tid]->hydro.parts[p].x[0],
                yy = ci_list[tid]->hydro.parts[p].x[1],
                zz = ci_list[tid]->hydro.parts[p].x[2];
          fprintf(fpcpu, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
                  ci_list[tid]->hydro.parts[p].rho,
                  ci_list[tid]->hydro.parts[p].density.rho_dh,
                  ci_list[tid]->hydro.parts[p].density.wcount,
                  ci_list[tid]->hydro.parts[p].density.wcount_dh,
                  ci_list[tid]->hydro.parts[p].viscosity.div_v,
                  ci_list[tid]->hydro.parts[p].density.rot_v[0]);
        }
      }

      fclose(fpcpu);
      double t_cpu_print_end = clock();
      double t_gpu_end = clock();
      double stream_destruction_start = clock();
      for (int bid = 0; bid < nBundles; ++bid)
        gpuErrchk(cudaStreamDestroy(stream[bid]));

      double stream_destruction_end = clock();
      printf("stream destruction time is %f seconds \n",
             (stream_destruction_end - stream_destruction_start) /
                 CLOCKS_PER_SEC);
      printf("time for cpu is %f\n", t_cpu_elapsed / CLOCKS_PER_SEC);
      printf("time for gpu is %f\n", (t_gpu_end - t_gpu_start - t_cpu_elapsed -
                                      t_cpu_print_end + t_cpu_print_start) /
                                         CLOCKS_PER_SEC);
      printf("finished looping through tasks\n");
      FILE *fp;
      fp = fopen("./res_GPU_density.txt", "w");
      for (int tid = 0; tid < count_tasks; tid++) {
        const int count = parts_gpu_list[tid][0].count;
        double t6 = clock();
        for (int p = 0; p < count; p++) {
          float xx = parts_gpu_list[tid][p].x[0],
                yy = parts_gpu_list[tid][p].x[1],
                zz = parts_gpu_list[tid][p].x[2];
          fprintf(
              fp, "%f %f %f %f %f %f %f %f %f\n", xx, yy, zz,
              parts_gpu_list[tid][p].rho, parts_gpu_list[tid][p].rho_dh,
              parts_gpu_list[tid][p].wcount, parts_gpu_list[tid][p].wcount_dh,
              parts_gpu_list[tid][p].div_v, parts_gpu_list[tid][p].rot_v[0]);
        }
      }
      fclose(fp);

      cudaProfilerStop();
      cudaFreeHost(parts_gpu_list);
      cudaDeviceReset();
      exit(0);
      ///* Mark that we have run this task on these cells */
      //#ifdef SWIFT_DEBUG_CHECKS
      //				if (ci_list[tid] != NULL) {
      //				  ci_list[tid]->tasks_executed[tasks[tid]->type]++;
      //				  ci_list[tid]->subtasks_executed[tasks-[tid]>subtype]++;
      //				}
      ////				if (cj != NULL) {
      ////				  cj->tasks_executed[t->type]++;
      ////				  cj->subtasks_executed[t->subtype]++;
      ////				}

      //				/* This runner is not doing a task
      //anymore */ 				r->t = NULL; #endif

      /* We're done with this task, see if we get a next one. */
      //      prev = t;
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
//#ifdef __cplusplus
// extern "C" {
//#endif

//#include <vector>

//#ifdef __cplusplus
//}
//#endif
void *runner_main2(void *data) {
  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  int nr_packs = sched->nr_pack_tasks;
  //  fprintf(stderr,"nr_packs is %i\n", sched->nr_pack_tasks);
  sched->nr_packs_done = 0;
  struct space *space = e->s;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);

  int launch_leftovers = 0;

  int packed = 0;
  int density = 0;
  int density_sub = 0;
  int unpacked = 0;
  int ghost_in = 0;
  cudaDeviceReset();
  int devId = 0; // find and print gpu device name
  struct cudaDeviceProp prop;

  cudaProfilerStart();
  cudaGetDeviceProperties(&prop, devId);
  if (r->cpuid == 0)
    fprintf(stderr, "Device : %s\n", prop.name);
  cudaSetDevice(devId);
  cudaFree(0);
  //  cudaProfilerStart();
  cudaError_t cu_error;
  // how many tasks do we want for each launch of GPU kernel
  const int target_n_tasks = 256;
  // how many tasks we want in each bundle (used for launching kernels in
  // different streams)
  const int bundle_size = 64;
  // array to keep track of GPU tasks
  struct task *tasks[target_n_tasks];
  for (int i = 0; i < target_n_tasks; i++) {
    tasks[i] = NULL;
  }
  // Keep track of first and last particles for each task (particle data is
  // arranged in long arrays containing particles from all the tasks we will
  // work with)
  int *task_first_part, *task_last_part;
  // Copy of the above residing on the GPU
  int *d_task_first_part, *d_task_last_part;

  // Arrays keeping track of the row numbers of the first and last particles
  // within each bundle. Required by the GPU code
  int *bundle_first_part, *bundle_last_part;
  int *d_bundle_first_part, *d_bundle_last_part;

  cudaMallocHost((void **)&task_first_part,
                 target_n_tasks * sizeof(int)); // Pinned allocation on host
  cudaMallocHost((void **)&task_last_part,
                 target_n_tasks * sizeof(int)); // Pinned allocation on host

  int nBundles = (target_n_tasks + bundle_size - 1) /
                 bundle_size; // nBundeles is the number of task bundles each
                              // thread has ==> Used to loop through bundles
  //  fprintf(stderr, "thread %i will try to pack %i task bundles\n", r->cpuid,
  //  nBundles);
  // keep track of how many bundles we pack
  int bundles_packed = 0;
  // first part and last part are the first and last particle ids (locally
  // within this thread)
  cudaMallocHost((void **)&bundle_first_part,
                 nBundles * sizeof(int)); // Pinned allocation on host
  cudaMallocHost((void **)&bundle_last_part,
                 nBundles * sizeof(int)); // Pinned allocation on host

  cudaMalloc((void **)&d_task_first_part, target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_task_last_part, target_n_tasks * sizeof(int));
  cudaMalloc((void **)&d_bundle_first_part, nBundles * sizeof(int));
  cudaMalloc((void **)&d_bundle_last_part, nBundles * sizeof(int));

  cudaEvent_t startEvent, stopEvent, dummyEvent;
  cudaEventCreate(&startEvent);
  cudaEventCreate(&stopEvent);
  cudaEventCreate(&dummyEvent);
  cudaStream_t stream[nBundles];
  int tasksperbundle = (target_n_tasks + nBundles - 1) / nBundles;

  for (int i = 0; i < nBundles; ++i)
    cudaStreamCreateWithFlags(&stream[i], cudaStreamNonBlocking);
  // pack length is the length of vector holding task data for GPU
  int pack_length = 0;
  int total_parts_packed = 0;
  int total_parts_done_cpu = 0;

  // Pointer referencing the length (pack_length) which will be incremented
  // depending on how many tasks we want before sending to GPU
  int *pack_lengthp = &pack_length;

  /*Allocate particle positions using arbitrary large array*/

  //  fprintf(stderr, "nr_cells is %i\n", space->tot_cells);
  //  fprintf(stderr, "nr_cells is %i\n", space->nr_cells);
  //  fprintf(stderr, "nr_parts is %i\n", space->nr_parts);
  int count_max_parts_tmp =
      2 * target_n_tasks * space->nr_parts /
      space->nr_cells; /*(1 << 18); which is 2^18=262144 particles*/
  //  fprintf(stderr, "count_max_parts_tmp is %i\n", count_max_parts_tmp);
  int tid = 0;
  //////////////////////////////////////////////
  /*****Initialise GPU suitable copies of variables and malloc pinned
   * memory*****/
  part_soa parts_soa;
  part_soa d_parts_soa;
  ////////All data contained within parts_soa struct. A lot, innit?
  int *d_tid_p;
  long long *d_id;
  double *d_x_p;
  double *d_y_p;
  double *d_z_p;
  float *d_ux;
  float *d_uy;
  float *d_uz;
  float *d_a_hydrox;
  float *d_a_hydroy;
  float *d_a_hydroz;
  float *d_mass;
  float *d_h;
  float *d_u;
  float *d_u_dt;
  float *d_rho;
  float *d_SPH_sum;
  float *d_locx;
  float *d_locy;
  float *d_locz;
  float *d_widthx;
  float *d_widthy;
  float *d_widthz;
  float *d_h_max;
  int *d_count_p;
  float *d_wcount;
  float *d_wcount_dh;
  float *d_rho_dh;
  float *d_rot_ux;
  float *d_rot_uy;
  float *d_rot_uz;
  float *d_div_v;
  float *d_div_v_previous_step;
  float *d_alpha_visc;
  float *d_v_sig;
  float *d_laplace_u;
  float *d_alpha_diff;
  float *d_f;
  float *d_soundspeed;
  float *d_h_dt;
  float *d_balsara;
  float *d_pressure;
  float *d_alpha_visc_max_ngb;
  /* timestep stuff */
  timebin_t *d_time_bin;
  timebin_t *d_wakeup;
  timebin_t *d_min_ngb_time_bin;
  char *d_to_be_synchronized;
  struct cell *mycell[count_max_parts_tmp];
  for (int i = 0; i < count_max_parts_tmp; i++)
    mycell[i] = (struct cell *)malloc(sizeof(struct cell));
  ////////now malloc space for all this data on the GPU. Sheesh

  cudaMalloc((void **)&(d_tid_p), sizeof(int) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_id), sizeof(long long) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_x_p), sizeof(double) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_y_p), sizeof(double) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_z_p), sizeof(double) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_ux), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_uy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_uz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_a_hydrox), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_a_hydroy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_a_hydroz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_mass), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_h), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_u), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_u_dt), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_rho), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_SPH_sum), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_locx), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_locy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_locz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_widthx), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_widthy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_widthz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_h_max), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_count_p), sizeof(int) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_wcount), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_wcount_dh), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_rho_dh), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_rot_ux), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_rot_uy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_rot_uz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_div_v), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_div_v_previous_step),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_alpha_visc), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_v_sig), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_laplace_u), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_alpha_diff), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_f), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_soundspeed), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_h_dt), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_balsara), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_pressure), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_alpha_visc_max_ngb),
             sizeof(float) * count_max_parts_tmp);
  /* timestep stuff */
  cudaMalloc((void **)&(d_time_bin), sizeof(timebin_t) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_wakeup), sizeof(timebin_t) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_min_ngb_time_bin),
             sizeof(timebin_t) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_to_be_synchronized),
             sizeof(char) * count_max_parts_tmp);
  ///////////Malloc Host arrays
  int *tid_p;
  cudaMallocHost((void **)&tid_p, count_max_parts_tmp * sizeof(int));
  long long *id;
  cudaMallocHost((void **)&id, count_max_parts_tmp * sizeof(long long));
  float *mass;
  cudaMallocHost((void **)&mass, count_max_parts_tmp * sizeof(float));
  float *h;
  cudaMallocHost((void **)&h, count_max_parts_tmp * sizeof(float));
  float *u;
  cudaMallocHost((void **)&u, count_max_parts_tmp * sizeof(float));
  float *u_dt;
  cudaMallocHost((void **)&u_dt, count_max_parts_tmp * sizeof(float));
  float *rho;
  cudaMallocHost((void **)&rho, count_max_parts_tmp * sizeof(float));
  float *SPH_sum;
  cudaMallocHost((void **)&SPH_sum, count_max_parts_tmp * sizeof(float));
  double *x_p;
  cudaMallocHost((void **)&x_p, count_max_parts_tmp * sizeof(double));
  double *y_p;
  cudaMallocHost((void **)&y_p, count_max_parts_tmp * sizeof(double));
  double *z_p;
  cudaMallocHost((void **)&z_p, count_max_parts_tmp * sizeof(double));
  float *ux;
  cudaMallocHost((void **)&ux, count_max_parts_tmp * sizeof(float));
  float *uy;
  cudaMallocHost((void **)&uy, count_max_parts_tmp * sizeof(float));
  float *uz;
  cudaMallocHost((void **)&uz, count_max_parts_tmp * sizeof(float));
  float *a_hydrox;
  cudaMallocHost((void **)&a_hydrox, count_max_parts_tmp * sizeof(float));
  float *a_hydroy;
  cudaMallocHost((void **)&a_hydroy, count_max_parts_tmp * sizeof(float));
  float *a_hydroz;
  cudaMallocHost((void **)&a_hydroz, count_max_parts_tmp * sizeof(float));
  float *locx;
  cudaMallocHost((void **)&locx, count_max_parts_tmp * sizeof(float));
  float *locy;
  cudaMallocHost((void **)&locy, count_max_parts_tmp * sizeof(float));
  float *locz;
  cudaMallocHost((void **)&locz, count_max_parts_tmp * sizeof(float));
  float *widthx;
  cudaMallocHost((void **)&widthx, count_max_parts_tmp * sizeof(float));
  float *widthy;
  cudaMallocHost((void **)&widthy, count_max_parts_tmp * sizeof(float));
  float *widthz;
  cudaMallocHost((void **)&widthz, count_max_parts_tmp * sizeof(float));
  float *h_max;
  cudaMallocHost((void **)&h_max, count_max_parts_tmp * sizeof(float));
  int *count_p;
  cudaMallocHost((void **)&count_p, count_max_parts_tmp * sizeof(int));
  float *wcount;
  cudaMallocHost((void **)&wcount, count_max_parts_tmp * sizeof(float));
  float *wcount_dh;
  cudaMallocHost((void **)&wcount_dh, count_max_parts_tmp * sizeof(float));
  float *rho_dh;
  cudaMallocHost((void **)&rho_dh, count_max_parts_tmp * sizeof(float));
  float *rot_ux;
  cudaMallocHost((void **)&rot_ux, count_max_parts_tmp * sizeof(float));
  float *rot_uy;
  cudaMallocHost((void **)&rot_uy, count_max_parts_tmp * sizeof(float));
  float *rot_uz;
  cudaMallocHost((void **)&rot_uz, count_max_parts_tmp * sizeof(float));
  float *div_v;
  cudaMallocHost((void **)&div_v, count_max_parts_tmp * sizeof(float));
  float *div_v_previous_step;
  cudaMallocHost((void **)&div_v_previous_step,
                 count_max_parts_tmp * sizeof(float));
  float *alpha_visc;
  cudaMallocHost((void **)&alpha_visc, count_max_parts_tmp * sizeof(float));
  float *v_sig;
  cudaMallocHost((void **)&v_sig, count_max_parts_tmp * sizeof(float));
  float *laplace_u;
  cudaMallocHost((void **)&laplace_u, count_max_parts_tmp * sizeof(float));
  float *alpha_diff;
  cudaMallocHost((void **)&alpha_diff, count_max_parts_tmp * sizeof(float));
  float *f;
  cudaMallocHost((void **)&f, count_max_parts_tmp * sizeof(float));
  float *soundspeed;
  cudaMallocHost((void **)&soundspeed, count_max_parts_tmp * sizeof(float));
  float *h_dt;
  cudaMallocHost((void **)&h_dt, count_max_parts_tmp * sizeof(float));
  float *balsara;
  cudaMallocHost((void **)&balsara, count_max_parts_tmp * sizeof(float));
  float *pressure;
  cudaMallocHost((void **)&pressure, count_max_parts_tmp * sizeof(float));
  float *alpha_visc_max_ngb;
  cudaMallocHost((void **)&alpha_visc_max_ngb,
                 count_max_parts_tmp * sizeof(float));
  /* timestep stuff */
  timebin_t *time_bin;
  cudaMallocHost((void **)&time_bin, count_max_parts_tmp * sizeof(timebin_t));
  timebin_t *wakeup;
  cudaMallocHost((void **)&wakeup, count_max_parts_tmp * sizeof(timebin_t));
  timebin_t *min_ngb_time_bin;
  cudaMallocHost((void **)&min_ngb_time_bin,
                 count_max_parts_tmp * sizeof(timebin_t));
  char *to_be_synchronized;
  cudaMallocHost((void **)&to_be_synchronized,
                 count_max_parts_tmp * sizeof(char));
  long long *partid_p;
  cudaMallocHost((void **)&partid_p, count_max_parts_tmp * sizeof(long long));

  float d_a = e->cosmology->a;
  float d_H = e->cosmology->H;

  // Stuff for writing debug data to file for validation
  char buf[20], buf2[20], buf3[20];
  snprintf(buf, sizeof(buf), "gpu_out%d", r->cpuid);
  snprintf(buf2, sizeof(buf2), "cpu_out%d", r->cpuid);
  snprintf(buf3, sizeof(buf3), "gpu_out_in%d", r->cpuid);
  FILE *fgpu = fopen(buf, "w");
  FILE *fcpu = fopen(buf2, "w");
  FILE *fgpuin = fopen(buf3, "w");

  FILE *fbefore = fopen("./original.txt", "w");
  FILE *fafter = fopen("./modified.txt", "w");
  FILE *fbefore_cpu = fopen("./original_cpu.txt", "w");
  FILE *fafter_cpu = fopen("./modified_cpu.txt", "w");
  // a list of the cells the GPU will work on
  struct cell **ci_list =
      (struct cell **)malloc(target_n_tasks * sizeof(struct cell *));
  struct cell *ci_list_2[target_n_tasks];
  struct cell *ci_list_cpu[sched->nr_pack_tasks];
  struct cell *ci_list_gpu[sched->nr_pack_tasks];
  // number of density self tasks executed
  int tasks_done_cpu = 0;
  int tasks_done_gpu = 0;
  int tasks_done_gpu_inc = 0;

  /* Main loop. */
  while (1) {
    //	Initialise timers to zero
    double time_for_density_cpu = 0.0;
    double time_for_density_cpu_sub = 0.0;
    double time_for_density_gpu = 0.0;
    //	tasks_done_cpu=0;
    //	tasks_done_gpu=0;
    /* Wait at the barrier. */
    engine_barrier(e);
    // Initialise packing counters
    pack_length = 0;
    int tasks_packed = 0;
    int total_tasks_packed_this_time = 0;
    double packing_time = 0.0;
    double time_for_copy_to_struct = 0.0;
    double tot_time_for_hard_memcpys = 0.0;
    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done)
      break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;
    /* Loop while there are tasks... */
    while (1) {
      // get the number of pack tasks for this scheduler
      nr_packs = sched->nr_pack_tasks;
      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, prev);
        TIMER_TOC(timer_gettask);

        /* Did I get anything? */
        if (t == NULL)
          break;
      }

      /* Get the cells. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        struct cell *ci_temp = ci;
        struct cell *cj_temp = cj;
        double shift[3];
        t->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
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
      /* Different types of tasks... */
      switch (t->type) {
      case task_type_gpu_unpack:
        if (t->subtype == task_subtype_gpu_unpack) {
          int do_nothing = 0;
          //			if(r->cpuid==0)printf("unpacking\n");
          if (unpacked < 1)
            fprintf(stderr, "unpacked cpu %i\n", r->cpuid);
          unpacked++;
        }
        break;
      case task_type_self:
        if (t->subtype == task_subtype_density) {
          const int count = ci->hydro.count;
          for (int i = 0; i < count; i++) {
            if (part_is_inhibited(&ci->hydro.parts[i], e))
              fprintf(stderr, "part inhibited\n");
          }
          double tstart = clock();
          runner_doself1_branch_density(r, ci);
          double tend = clock();
          tasks_done_cpu++;
          time_for_density_cpu += (tend - tstart) / CLOCKS_PER_SEC;
          if (density < 1)
            fprintf(stderr, "density self cpu %i\n", r->cpuid);
          density++;
          //        	    if(ci->hydro.gpu_pack->skip==0||ci->hydro.gpu_unpack->skip==0){
          //        	    	fprintf(stderr,"Doing density before skipping
          //        pack\n"); 	    	exit(0);
          //        	    }
          //        	    for(int i=0; i<ci->hydro.count; i++){
          //        	    	fprintf(fafter_cpu, "%f %f %f %f %f\n",
          //        ci->hydro.parts[i].x[0], ci->hydro.parts[i].x[1],
          //        ci->hydro.parts[i].x[2], ci->hydro.parts[i].rho,
          //        ci->hydro.parts[i].density.rho_dh);
          //        	    }
        } else if (t->subtype == task_subtype_gpu_pack) {
          if (packed < 1)
            fprintf(stderr, "packed by cpu %i\n", r->cpuid);
          packed++;
          ci_list_2[tasks_packed] = ci;
          ci_list_gpu[tasks_done_gpu] = ci;
          // timers for how long this all takes
          struct timespec t0, t1, tp0, tp1, dt;
          clock_gettime(CLOCK_REALTIME, &t0);
          double tstart = clock();
          clock_gettime(CLOCK_REALTIME, &tp0);
          double pack_start = clock();
          // create a list of pointers to the tasks
          //        	tasks[tasks_packed]=t;
          //        	struct cell *ci2_p=t->ci;
          //        	fprintf(stderr,"Before\n");
          //        	memcpy(&ci_list_2[tasks_packed], &ci, sizeof(ci));

          //        	*ci_list_2[tasks_packed]=*ci;
          //        	fprintf(stderr,"After\n");
          //        	struct cell *ci2=&ci_list_2[tasks_packed];
          //        	struct cell *ci2=&ci_list_2[tasks_packed];
          // identisy row in particle data arrays where this task starts
          task_first_part[tasks_packed] = pack_length;
          // id for the task
          tid = tasks_packed;

          // This re-arranges the particle data from cell->hydro->parts into a
          // series of long arrays
          runner_doself1_gpu_pack(
              r, ci, 0, pack_lengthp, x_p, y_p, z_p, tid, tid_p, id, ux, uy, uz,
              a_hydrox, a_hydroy, a_hydroz, mass, h, u, u_dt, rho, SPH_sum,
              locx, locy, locz, widthx, widthy, widthz, h_max, count_p, wcount,
              wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz, div_v,
              div_v_previous_step, alpha_visc, v_sig, laplace_u, alpha_diff, f,
              soundspeed, h_dt, balsara, pressure, alpha_visc_max_ngb, time_bin,
              wakeup, min_ngb_time_bin, to_be_synchronized, count_max_parts_tmp,
              fgpuin);
          /*Print out particle data before launching kernel*/
          //            if(r->cpuid==0){
          //            for(int i=pack_length-ci->hydro.count; i<pack_length;
          //            i++){
          ////            for(int i=0; i<ci2->hydro.count; i++){
          //            	fprintf(fgpuin, "%f %f %f %f %f %f\n", x_p[i],
          //            y_p[i], z_p[i], rho[i], rho_dh[i], div_v[i]);
          //            }
          //            }
          //
          //            for(int i=0; i<ci->hydro.count; i++){
          //				fprintf(fbefore_cpu, "%f %f %f %f %f %f\n",
          //ci->hydro.parts[i].x[0], ci->hydro.parts[i].x[1],
          //ci->hydro.parts[i].x[2], ci->hydro.parts[i].rho,
          //ci->hydro.parts[i].density.rho_dh,
          //ci->hydro.parts[i].viscosity.div_v);
          //			}

          // identify the row in the array where this task ends (row id of its
          // last particle)
          task_last_part[tasks_packed] = pack_length;
          // Identify first particle for each bundle of tasks
          if (tasks_packed % bundle_size == 0) {
            int bid = tasks_packed / bundle_size;
            bundle_first_part[bid] = task_first_part[tasks_packed];
          }
          tasks_packed++;
          tasks_done_gpu_inc++;
          total_tasks_packed_this_time++;
          //            atomic_inc(&sched->nr_packs_done);

          // Currently designed to work in serial (one thread) with the number
          // of gpu tasks created divisible by target_n_tasks
          //            fprintf(stderr, "n_pack_tasks_left is %i and
          //            n_threads*tasrget_n_tasks is %i\n",
          //            (nr_packs-sched->nr_packs_done),
          //            (e->nr_threads*target_n_tasks));
          // clumsy workaround for threads being unleashed whilst no tasks have
          // been generated
          //            if(nr_packs==0)nr_packs=1<<18;
          //            if((nr_packs-sched->nr_packs_done)<(e->nr_threads*target_n_tasks)){
          int qid = r->qid;
          atomic_dec(&sched->queues[qid].n_packs_self_left);
          //            fprintf(stderr, "in qid %i npacks_left is %i\n", qid,
          //            sched->queues[qid].n_packs_left);
          if ((sched->queues[qid].n_packs_self_left == 0)) {
            //            if((nr_packs-sched->nr_packs_done)==0){
            launch_leftovers = 1;
            //            	fprintf(stderr,"thread %i signalled tasks left
            //            lower than splittable. n left is %i, ntotal is %i\n",
            //            r->cpuid, nr_packs-sched->nr_packs_done, nr_packs);

            //            	engine_barrier(e);
            //            	int waiter=
            //            swift_barrier_wait(&e->wait_barrier);
          }
          double pack_end = clock();
          clock_gettime(CLOCK_REALTIME, &tp1);
          //            packing_time+=(pack_end-pack_start)/CLOCKS_PER_SEC;
          packing_time += (tp1.tv_sec - tp0.tv_sec) +
                          (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
          if (tasks_packed == target_n_tasks || launch_leftovers) {
            //            	if(r->cpuid==0){
            //            	fprintf(stderr, "packed %i tasks\n",
            //            tasks_packed);
            //				for(int i=0; i<pack_length; i++){
            //	//            for(int i=0; i<ci2->hydro.count; i++){
            //					fprintf(fgpuin, "%f %f %f %f %f %f\n", x_p[i], y_p[i],
            //z_p[i], rho[i], rho_dh[i], div_v[i]);
            //				}
            //				}
            // Printout particles before packing
            //            	for(int i=0; i<pack_length; i++){
            //					fprintf(fbefore, "%f %f %f %f %f %f\n", x_p[i], y_p[i],
            //z_p[i], rho[i], rho_dh[i], div_v[i]);
            //				}
            ///////////Might need tweaking!!!!//////////
            int nBundles_temp = nBundles;
            //            	if((nr_packs-sched->nr_packs_done)<(e->nr_threads*target_n_tasks)){
            if (launch_leftovers) {
              int bid = tasks_packed / bundle_size;
              bundle_first_part[bid] = task_first_part[tasks_packed - 1];
              //            		fprintf(stderr,"CPUid %i Tasks packed %i
              //            bundle id %i task_first_part %i tasks_first_part-1
              //            %i pack length %i\n", r->cpuid, tasks_packed, bid,
              //            task_first_part[tasks_packed],
              //            task_first_part[tasks_packed-1], pack_length);
              nBundles_temp = bid;
              //            		exit(0);
            }
            tasks_done_gpu += tasks_packed;
            double t_gpu_elapsed = clock();
            double t_gpu_copy_to_struct = clock();
            // Identify last particle for each bundle of tasks
            for (int bid = 0; bid < nBundles_temp - 1; bid++) {
              bundle_last_part[bid] = bundle_first_part[bid + 1];
              //					fprintf(stderr,"CPUid %
              //i bundle last part %i\n", r->cpuid, bundle_last_part[bid]);
            }
            bundle_last_part[nBundles_temp - 1] = pack_length;
            // Point parts_soa to the GPU data
            parts_soa.tid_p = d_tid_p;
            parts_soa.locx = d_locx;
            parts_soa.locy = d_locy;
            parts_soa.locz = d_locz;
            parts_soa.h = d_h;
            parts_soa.mass = d_mass;
            parts_soa.x_p = d_x_p;
            parts_soa.y_p = d_y_p;
            parts_soa.z_p = d_z_p;
            parts_soa.rho = d_rho;
            parts_soa.rho_dh = d_rho_dh;
            parts_soa.wcount = d_wcount;
            parts_soa.wcount_dh = d_wcount_dh;
            parts_soa.ux = d_ux;
            parts_soa.uy = d_uy;
            parts_soa.uz = d_uz;
            parts_soa.div_v = d_div_v;
            parts_soa.rot_ux = d_rot_ux;
            parts_soa.rot_uy = d_rot_uy;
            parts_soa.rot_uz = d_rot_uz;
            parts_soa.count_p = d_count_p;

            time_for_copy_to_struct +=
                (clock() - t_gpu_copy_to_struct) / CLOCKS_PER_SEC;
            // Copies data for start and end of each tasks and bundle to GPU
            double time_for_hard_memcpys = clock();
            cudaMemcpy(d_task_first_part, task_first_part,
                       tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_task_last_part, task_last_part,
                       tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_bundle_first_part, bundle_first_part,
                       nBundles_temp * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_bundle_last_part, bundle_last_part,
                       nBundles_temp * sizeof(int), cudaMemcpyHostToDevice);
            tot_time_for_hard_memcpys +=
                (clock() - time_for_hard_memcpys) / CLOCKS_PER_SEC;

            int shared_memory_used = 0;
            // time how long cpu/gpu copies and kernel take
            double t_gpu_copyandkernel_begin = clock();

            for (int bid = 0; bid < nBundles_temp; bid++) {
              int max_parts = 0;
              int parts_in_bundle = 0;

              for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
                   tid++) {
                if (tid < tasks_packed) {
                  int count = task_last_part[tid] - task_first_part[tid];
                  parts_in_bundle += count;
                  max_parts = max(max_parts, count);
                }
              }

              int first_part_tmp = bundle_first_part[bid];
              int bundle_size = bundle_last_part[bid] - first_part_tmp;
              // copy data to GPU
              cudaMemcpyAsync(&d_tid_p[first_part_tmp], &tid_p[first_part_tmp],
                              bundle_size * sizeof(int), cudaMemcpyHostToDevice,
                              stream[bid]);
              cudaMemcpyAsync(&d_locx[first_part_tmp], &locx[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_locy[first_part_tmp], &locy[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_locz[first_part_tmp], &locz[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_h[first_part_tmp], &h[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_mass[first_part_tmp], &mass[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_x_p[first_part_tmp], &x_p[first_part_tmp],
                              bundle_size * sizeof(double),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_y_p[first_part_tmp], &y_p[first_part_tmp],
                              bundle_size * sizeof(double),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_z_p[first_part_tmp], &z_p[first_part_tmp],
                              bundle_size * sizeof(double),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_rho[first_part_tmp], &rho[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_rho_dh[first_part_tmp],
                              &rho_dh[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_wcount[first_part_tmp],
                              &wcount[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_wcount_dh[first_part_tmp],
                              &wcount_dh[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_ux[first_part_tmp], &ux[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_uy[first_part_tmp], &uy[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_uz[first_part_tmp], &uz[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_div_v[first_part_tmp], &div_v[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_rot_ux[first_part_tmp],
                              &rot_ux[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_rot_uy[first_part_tmp],
                              &rot_uy[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_rot_uz[first_part_tmp],
                              &rot_uz[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyHostToDevice, stream[bid]);
              cudaMemcpyAsync(&d_count_p[first_part_tmp],
                              &count_p[first_part_tmp],
                              bundle_size * sizeof(int), cudaMemcpyHostToDevice,
                              stream[bid]);
              //					cu_error =
              //cudaPeekAtLastError();//cudaGetLastError();        // Get error
              //code 					if ( cu_error != cudaSuccess )
              //					{
              //						fprintf(stderr,
              //"CUDA kernel Error host 2 device memcpy: %s cpuid id is: %i\n ",
              //cudaGetErrorString(cu_error), r->cpuid); 						exit(0);
              //					}
              int tid = 0;
              int offset = bid * tasksperbundle;
              int sizecells = 0, sizeparts = 0;
              int tasks_left = tasksperbundle;

              if (bid == nBundles_temp - 1) {
                tasks_left =
                    tasks_packed - (nBundles_temp - 1) * tasksperbundle;
                sizeparts = tasks_left * sizeof(struct part_gpu *);
              }
              // Will launch a 2d grid of GPU thread blocks (number of tasks is
              // the y dimension and max_parts is the x dimension
              int numBlocks_y = tasks_left;
              int block_size = BLOCK_SIZE;
              int numBlocks_x = (max_parts + block_size - 1) / block_size;
              int numBlocks = (parts_in_bundle + block_size - 1) / block_size;
              int shared_memory_block =
                  BLOCK_SIZE * (8 * sizeof(float) + sizeof(float *));
              int bundle_part_0 = bundle_first_part[bid];
              int bundle_first_task = tid_p[bundle_part_0];
              const char *loop_type = "density";
              // Launch the kernel
              launch_cuda_kernel_bundles_revised_soa(
                  parts_soa, d_task_first_part, d_task_last_part,
                  d_bundle_first_part, d_bundle_last_part, numBlocks, d_a, d_H,
                  loop_type, stream[bid], bid, block_size, tasks_packed,
                  tasksperbundle, numBlocks_x, numBlocks_y, tid, offset,
                  bundle_first_task, max_parts);
              //					cu_error =
              //cudaPeekAtLastError();//cudaGetLastError();        // Get error
              //code 					if ( cu_error != cudaSuccess )
              //					{
              //					  fprintf(stderr, "CUDA
              //kernel Error kernel: %s cpuid is %i\n ",
              //cudaGetErrorString(cu_error), r->cpuid);
              //					}
              // Copy data back to CPU
              //					cudaMemcpyAsync(&tid_p[first_part_tmp],
              //&d_tid_p[first_part_tmp], bundle_size*sizeof(int),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&locx[first_part_tmp],
              //&d_locx[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&locy[first_part_tmp],
              //&d_locy[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&locz[first_part_tmp],
              //&d_locz[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&h[first_part_tmp],
              //&d_h[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&mass[first_part_tmp],
              //&d_mass[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&x_p[first_part_tmp],
              //&d_x_p[first_part_tmp], bundle_size*sizeof(double),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&y_p[first_part_tmp],
              //&d_y_p[first_part_tmp], bundle_size*sizeof(double),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&z_p[first_part_tmp],
              //&d_z_p[first_part_tmp], bundle_size*sizeof(double),
              //cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&rho[first_part_tmp], &d_rho[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&rho_dh[first_part_tmp],
                              &d_rho_dh[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&wcount[first_part_tmp],
                              &d_wcount[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&wcount_dh[first_part_tmp],
                              &d_wcount_dh[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&ux[first_part_tmp],
              //&d_ux[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&uy[first_part_tmp],
              //&d_uy[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&uz[first_part_tmp],
              //&d_uz[first_part_tmp], bundle_size*sizeof(float),
              //cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&div_v[first_part_tmp], &d_div_v[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&rot_ux[first_part_tmp],
                              &d_rot_ux[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&rot_uy[first_part_tmp],
                              &d_rot_uy[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              cudaMemcpyAsync(&rot_uz[first_part_tmp],
                              &d_rot_uz[first_part_tmp],
                              bundle_size * sizeof(float),
                              cudaMemcpyDeviceToHost, stream[bid]);
              //					cudaMemcpyAsync(&count_p[first_part_tmp],
              //&d_count_p[first_part_tmp], bundle_size*sizeof(int),
              //cudaMemcpyDeviceToHost, stream[bid]); 					if ( cu_error !=
              //cudaSuccess )
              //					{
              //					  fprintf(stderr, "CUDA
              //kernel Error device2host: %s \n ",
              //cudaGetErrorString(cu_error));
              //					}
            }
            cudaDeviceSynchronize();
            //                parts_soa.tid_p=tid_p;
            //				parts_soa.locx=locx;
            //				parts_soa.locy=locy;
            //				parts_soa.locz=locz;
            //				parts_soa.h=h;
            //				parts_soa.mass=mass;
            //				parts_soa.x_p=x_p;
            //				parts_soa.y_p=y_p;
            //				parts_soa.z_p=z_p;
            //				parts_soa.rho=rho;
            //				parts_soa.rho_dh=rho_dh;
            //				parts_soa.wcount=wcount;
            //				parts_soa.wcount_dh=wcount_dh;
            //				parts_soa.ux=ux;
            //				parts_soa.uy=uy;
            //				parts_soa.uz=uz;
            //				parts_soa.div_v=div_v;
            //				parts_soa.rot_ux=rot_ux;
            //				parts_soa.rot_uy=rot_uy;
            //				parts_soa.rot_uz=rot_uz;
            //				parts_soa.count_p=count_p;
            t_gpu_elapsed = (clock() - t_gpu_elapsed) / CLOCKS_PER_SEC;
            //				fprintf(stderr,"time for gpu to calculate %i tasks is
            //%f\n", tasks_packed, t_gpu_elapsed);
            total_parts_packed += pack_length;
            //				if(r->cpuid==0){
            //                  for(int i=0; i<pack_length; i++){
            //                	  fprintf(fgpuin, "%f %f %f %f %f\n", x_p[i],
            //                y_p[i], z_p[i], rho[i], rho_dh[i]);
            //                  }
            ////                  fprintf(stderr,"pack length %i\n",
            ///pack_length);
            //				}
            pack_length = 0;
            /////////////////Unpack
            ////////////////////////////////////////////////////////////
            for (int tid = 0; tid < tasks_packed; tid++) {
              //					ci2=&ci_list_2[tid];
              //					int
              //count=ci2->hydro.count; 					struct *cell ci2=ci_list_2[tid]; 					struct
              //cell *ci_ttemp=ci_list_2[tid]; 					struct part
              //*parts_tmp=ci_ttemp->hydro.parts; 					int
              //count=ci_ttemp->hydro.count; 					if(r->cpuid==0){ 					for(int i=0;
              //i<count; i++){ 					  fprintf(fgpuin, "%f %f %f %f %f\n",
              //parts_tmp[i].x[0], parts_tmp[i].x[1], parts_tmp[i].x[2],
              //parts_tmp[i].rho, parts_tmp[i].density.rho_dh);
              //					}
              //					}
              //					if(r->cpuid==0){
              //						fprintf(stderr,"count
              //%i packlength %i firstpart %i lastpart %i\n",count, pack_length,
              //task_first_part[tid], task_last_part[tid]);
              ////						exit(0);
              //					}
              clock_gettime(CLOCK_REALTIME, &tp0);

              //					runner_doself1_gpu_unpack(r,
              //ci_list_2[tid], 0, pack_lengthp, 					            x_p, y_p, z_p, tid, tid_p, id,
              //ux, uy, uz, a_hydrox, 					            a_hydroy, a_hydroz, mass, h, u, u_dt, rho,
              //SPH_sum, locx, locy, locz, 					            widthx, widthy, widthz, h_max,
              //count_p, wcount, wcount_dh, rho_dh, rot_ux, rot_uy, 					            rot_uz,
              //div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
              //alpha_diff, f, soundspeed, 					            h_dt, balsara, pressure,
              //alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
              //					            to_be_synchronized,
              //count_max_parts_tmp, fgpuin);

              clock_gettime(CLOCK_REALTIME, &tp1);

              packing_time += (tp1.tv_sec - tp0.tv_sec) +
                              (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
              //					if(r->cpuid==0){
              //						for(int i=0;
              //i<count; i++){ 							fprintf(fgpu, "%f %f %f %f %f\n",
              //parts_tmp[i].x[0], parts_tmp[i].x[1], parts_tmp[i].x[2],
              //parts_tmp[i].rho, parts_tmp[i].density.rho_dh);
              //						}
              //					}
            }
            //				if(r->cpuid==0){
            //					fprintf(stderr,"pack length %i\n",
            //*pack_lengthp);
            //				}
            ///////////////////////////////////////////////////////////////////////////////////
            //				    if(r->cpuid==0){
            //						for(int i=0; i<total_parts_packed;
            //i++){ 							fprintf(fgpu, "%f %f %f %f %f\n", x_p[i], y_p[i], z_p[i],
            //rho[i], rho_dh[i]);
            //						}
            //						for(int tid=0; tid<tasks_packed;
            //tid++){ 							for(int i=0; i<ci_list_2[tid]->hydro.count; i++){
            //								fprintf(fafter, "%f %f %f %f %f\n",
            //ci_list_2[tid]->hydro.parts[i].x[0],
            //ci_list_2[tid]->hydro.parts[i].x[1],
            //ci_list_2[tid]->hydro.parts[i].x[2],
            //ci_list_2[tid]->hydro.parts[i].rho,
            //ci_list_2[tid]->hydro.parts[i].density.rho_dh);
            //							}
            //						}
            //				    }
            //					count=tasks[tid]->ci->hydro.count;
            //					for(int i=0; i<tasks[tid]->ci->hydro.count;
            //i++){ 						double diff_rho=tasks[tid]->ci->hydro.parts[i].rho -
            //ci2->hydro.parts[i].rho; 						double
            //diff_rho_dh=tasks[tid]->ci->hydro.parts[i].density.rho_dh -
            //ci2->hydro.parts[i].density.rho_dh;
            ////						fprintf(fafter_cpu, "%f %f %f %f %f\n",
            ///tasks[tid]->ci->hydro.parts[i].x[0],
            ///tasks[tid]->ci->hydro.parts[i].x[1],
            ///tasks[tid]->ci->hydro.parts[i].x[2], diff_rho, diff_rho_dh);
            //						fprintf(fafter_cpu, "%f %f %f %f %f\n",
            //tasks[tid]->ci->hydro.parts[i].x[0],
            //tasks[tid]->ci->hydro.parts[i].x[1],
            //tasks[tid]->ci->hydro.parts[i].x[2],
            //tasks[tid]->ci->hydro.parts[i].rho,
            //tasks[tid]->ci->hydro.parts[i].density.rho_dh);
            //					}
            //					if(r->cpuid==0){
            //						fprintf(stderr,"cpu 0 output its
            //data\n");
            //					}
            pack_length = 0;
            //				}

            //////////////////////////////////////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////////////////////////////////
            /*Ignore this. Purely for outputting data to check the density sums
             * were done correctly*/

            //				for(int tid=0; tid<tasks_packed; tid++){
            //					const int
            //count=ci_list[tid]->hydro.count; 					for(int p=0; p<count; p++){ 						float
            //xx=ci_list[tid]->hydro.parts[p].x[0],
            //yy=ci_list[tid]->hydro.parts[p].x[1],
            //zz=ci_list[tid]->hydro.parts[p].x[2]; 						fprintf(fcpu, "%f %f %f %f
            //%f %f %f %f %f\n", xx, yy, zz, ci_list[tid]->hydro.parts[p].rho,
            //						ci_list[tid]->hydro.parts[p].density.rho_dh,
            //ci_list[tid]->hydro.parts[p].density.wcount,
            //						ci_list[tid]->hydro.parts[p].density.wcount_dh,
            //ci_list[tid]->hydro.parts[p].viscosity.div_v 						,
            //ci_list[tid]->hydro.parts[p].density.rot_v[0]);
            //					}
            //				}
            //				double t_cpu_elapsed=clock();
            // This was to re-run the tasks on the CPU to compare timings
            //				for(int tid=0; tid<tasks_packed; tid++){
            //					 runner_doself1_branch_density(r,
            //ci_list[tid]);
            //				}
            //				t_cpu_elapsed=(clock()-t_cpu_elapsed)/CLOCKS_PER_SEC;
            //				fprintf(stderr,"packing time is %f gpu time is %f cpu time
            //is %f cpuid is %i\n", packing_time, t_gpu_elapsed, t_cpu_elapsed,
            //r->cpuid); 				for(int tid=0; tid<tasks_packed; tid++){ 					const int
            //count=ci_list[tid]->hydro.count; 					for(int p=0; p<count; p++){ 						float
            //xx=ci_list[tid]->hydro.parts[p].x[0],
            //yy=ci_list[tid]->hydro.parts[p].x[1],
            //zz=ci_list[tid]->hydro.parts[p].x[2]; 						fprintf(fcpu, "%f %f %f %f
            //%f %f %f %f %f\n", xx, yy, zz, ci_list[tid]->hydro.parts[p].rho,
            //						ci_list[tid]->hydro.parts[p].density.rho_dh,
            //ci_list[tid]->hydro.parts[p].density.wcount,
            //						ci_list[tid]->hydro.parts[p].density.wcount_dh,
            //ci_list[tid]->hydro.parts[p].viscosity.div_v 						,
            //ci_list[tid]->hydro.parts[p].density.rot_v[0]);
            //					}
            //				}
            //				for(int p=0; p<pack_length; p++){
            //					float xx=parts_soa.x_p[p], yy=parts_soa.y_p[p],
            //zz=parts_soa.z_p[p]; 					fprintf(fgpu, "%f %f %f %f %f %f %f %f %f %i
            //%i %f\n", xx, yy, zz, parts_soa.rho[p], 					parts_soa.rho_dh[p],
            //parts_soa.wcount[p], parts_soa.wcount_dh[p], 					parts_soa.div_v[p],
            //parts_soa.rot_ux[p], parts_soa.tid_p[p], parts_soa.count_p[p],
            //parts_soa.h[p]);
            //				}
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //				cudaProfilerStop();
            tasks_packed = 0;
            launch_leftovers = 0;
          }
          double tend = clock();
          // increment time taken for density calcs on gpu
          clock_gettime(CLOCK_REALTIME, &t1);
          //            packing_time+=(pack_end-pack_start)/CLOCKS_PER_SEC;
          time_for_density_gpu += (t1.tv_sec - t0.tv_sec) +
                                  (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
          break;
        }
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_doself1_branch_gradient(r, ci);
#endif
        else if (t->subtype == task_subtype_force)
          runner_doself2_branch_force(r, ci);
        else if (t->subtype == task_subtype_limiter)
          runner_doself1_branch_limiter(r, ci);
        else if (t->subtype == task_subtype_grav)
          runner_doself_recursive_grav(r, ci, 1);
        else if (t->subtype == task_subtype_external_grav)
          runner_do_grav_external(r, ci, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_doself_branch_stars_density(r, ci);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_doself_branch_rt_inject(r, ci, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_doself_branch_sinks_compute_formation(r, ci);
        else
          error("Self Unknown/invalid task subtype (%s).",
                subtaskID_names[t->subtype]);
        break;

      case task_type_pair:
        if (t->subtype == task_subtype_density)
          runner_dopair1_branch_density(r, ci, cj);
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_dopair1_branch_gradient(r, ci, cj);
#endif
        else if (t->subtype == task_subtype_force)
          runner_dopair2_branch_force(r, ci, cj);
        else if (t->subtype == task_subtype_limiter)
          runner_dopair1_branch_limiter(r, ci, cj);
        else if (t->subtype == task_subtype_grav)
          runner_dopair_recursive_grav(r, ci, cj, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_dopair_branch_stars_density(r, ci, cj);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_dopair_branch_rt_inject(r, ci, cj, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_dopair_branch_sinks_compute_formation(r, ci, cj);
        else
          error("Pair Unknown/invalid task subtype (%s/%s).",
                taskID_names[t->type], subtaskID_names[t->subtype]);
        break;

      case task_type_sub_self:
        if (t->subtype == task_subtype_density) {
          struct timespec t0, t1, dt;

          const int count = ci->hydro.count;
          for (int i = 0; i < count; i++) {
            if (part_is_inhibited(&ci->hydro.parts[i], e))
              fprintf(stderr, "part inhibited\n");
          }
          if (density_sub < 1)
            fprintf(stderr, "density sub cpu %i\n", r->cpuid);
          density_sub++;
          clock_gettime(CLOCK_REALTIME, &t0);
          runner_dosub_self1_density(r, ci, 1);
          clock_gettime(CLOCK_REALTIME, &t1);
          ci_list_cpu[tasks_done_cpu] = ci;
          tasks_done_cpu++;
          time_for_density_cpu_sub +=
              (t1.tv_sec - t0.tv_sec) +
              (t1.tv_nsec - t0.tv_nsec) /
                  1000000000.0; // difftime(&t1,
                                // &t0);//(tend-tstart)/CLOCKS_PER_SEC;
          total_parts_done_cpu += ci->hydro.count;
          //            if(ci->hydro.gpu_pack->skip==0||ci->hydro.gpu_unpack->skip==0){
          //              fprintf(stderr,"Doing density before skipping
          //              pack\n"); exit(0);
          //            }

          //            fprintf(stderr,"sub_density interaction\n");
        }
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_dosub_self1_gradient(r, ci, 1);
#endif
        else if (t->subtype == task_subtype_force)
          runner_dosub_self2_force(r, ci, 1);
        else if (t->subtype == task_subtype_limiter)
          runner_dosub_self1_limiter(r, ci, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_dosub_self_stars_density(r, ci, 1);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_dosub_self_rt_inject(r, ci, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_dosub_self_sinks_compute_formation(r, ci, 1);
        else
          error("Sub Self Unknown/invalid task subtype (%s/%s).",
                taskID_names[t->type], subtaskID_names[t->subtype]);
        break;

      case task_type_sub_pair:
        if (t->subtype == task_subtype_density)
          runner_dosub_pair1_density(r, ci, cj, 1);
#ifdef EXTRA_HYDRO_LOOP
        else if (t->subtype == task_subtype_gradient)
          runner_dosub_pair1_gradient(r, ci, cj, 1);
#endif
        else if (t->subtype == task_subtype_force)
          runner_dosub_pair2_force(r, ci, cj, 1);
        else if (t->subtype == task_subtype_limiter)
          runner_dosub_pair1_limiter(r, ci, cj, 1);
        else if (t->subtype == task_subtype_stars_density)
          runner_dosub_pair_stars_density(r, ci, cj, 1);
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
        else if (t->subtype == task_subtype_rt_inject)
          runner_dosub_pair_rt_inject(r, ci, cj, 1);
        else if (t->subtype == task_subtype_sink_compute_formation)
          runner_dosub_pair_sinks_compute_formation(r, ci, cj, 1);
        else
          error("Sub Pair Unknown/invalid task subtype (%s/%s).",
                taskID_names[t->type], subtaskID_names[t->subtype]);
        break;

      case task_type_sort:
        /* Cleanup only if any of the indices went stale. */
        runner_do_hydro_sort(
            r, ci, t->flags,
            ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin, 1);
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
        if (ghost_in < 1)
          fprintf(stderr, "ghosted cpu %i\n", r->cpuid);
        ghost_in++;
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
      case task_type_logger:
        runner_do_logger(r, ci, 1);
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
#ifdef WITH_MPI
      case task_type_send:
        if (t->subtype == task_subtype_tend_part) {
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_gpart) {
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_spart) {
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_bpart) {
          free(t->buff);
        } else if (t->subtype == task_subtype_sf_counts) {
          free(t->buff);
        } else if (t->subtype == task_subtype_part_swallow) {
          free(t->buff);
        } else if (t->subtype == task_subtype_bpart_merger) {
          free(t->buff);
        }
        break;
      case task_type_recv:
        if (t->subtype == task_subtype_tend_part) {
          cell_unpack_end_step_hydro(ci, (struct pcell_step_hydro *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_gpart) {
          cell_unpack_end_step_grav(ci, (struct pcell_step_grav *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_spart) {
          cell_unpack_end_step_stars(ci, (struct pcell_step_stars *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_tend_bpart) {
          cell_unpack_end_step_black_holes(
              ci, (struct pcell_step_black_holes *)t->buff);
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
        } else if (t->subtype == task_subtype_part_swallow) {
          cell_unpack_part_swallow(ci, (struct black_holes_part_data *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_bpart_merger) {
          cell_unpack_bpart_swallow(ci,
                                    (struct black_holes_bpart_data *)t->buff);
          free(t->buff);
        } else if (t->subtype == task_subtype_limiter) {
          runner_do_recv_part(r, ci, 0, 1);
        } else if (t->subtype == task_subtype_gpart) {
          runner_do_recv_gpart(r, ci, 1);
        } else if (t->subtype == task_subtype_spart) {
          runner_do_recv_spart(r, ci, 1, 1);
        } else if (t->subtype == task_subtype_bpart_rho) {
          runner_do_recv_bpart(r, ci, 1, 1);
        } else if (t->subtype == task_subtype_bpart_swallow) {
          runner_do_recv_bpart(r, ci, 0, 1);
        } else if (t->subtype == task_subtype_bpart_feedback) {
          runner_do_recv_bpart(r, ci, 0, 1);
        } else if (t->subtype == task_subtype_multipole) {
          cell_unpack_multipoles(ci, (struct gravity_tensors *)t->buff);
          free(t->buff);
        } else {
          error("Unknown/invalid task subtype (%d).", t->subtype);
        }
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
      case task_type_stars_resort:
        runner_do_stars_resort(r, t->ci, 1);
        break;
      case task_type_sink_formation:
        runner_do_sink_formation(r, t->ci);
        break;
      case task_type_fof_self:
        runner_do_fof_self(r, t->ci, 1);
        break;
      case task_type_fof_pair:
        runner_do_fof_pair(r, t->ci, t->cj, 1);
        break;
      case task_type_rt_ghost1:
        runner_do_rt_ghost1(r, t->ci, 1);
        break;
      default:
        error("Unknown/invalid task type (%d).", t->type);
      }

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
      t = scheduler_done(sched, t);

      //      if(t->subtype==task_subtype_gpu_pack){
      //    	  //This is a pointer to the cell's unpack self task.
      //    Incremented to wait if task type is packing. Wait increment will be
      //    removed after the task is executed on the GPU
      ////    	  ci->hydro.gpu_unpack->wait++;
      //    	  //This is a pointer to the cell's density self task.
      //    Incremented to wait if task type is packing. Wait increment will be
      //    removed after the task is executed on the GPU
      ////    	  t->ci->hydro.density_self->wait++;
      //    	  t = scheduler_done(sched, t);
      //      }
      //      else{
      //		  t = scheduler_done(sched, t);
      //      }
      //	  if((nr_packs-sched->nr_packs_done)<(e->nr_threads*target_n_tasks)
      //&& sched->nr_pack_tasks > 0){ 		sched->nr_packs_done=0; 		fprintf(stderr,
      //"sched->nr_packs set to zero by thread %i\n", r->cpuid);
      //	  }
    } /* main loop. */
      //    if(r->cpuid==0){
      //		sched->nr_packs_done=0;
      ////		fprintf(stderr,"sched->nr_packs set to zero\n");
      //	}
      //    if(r->cpuid==0)
    fprintf(stderr,
            "cpu time for %i tasks and %i parts %f for sub %f Gpu did %i tasks"
            " and %i parts in %f packing time %f total_tasks_packed %i\n",
            tasks_done_cpu, total_parts_done_cpu, time_for_density_cpu,
            time_for_density_cpu_sub, tasks_done_gpu, total_parts_packed,
            time_for_density_gpu, packing_time, total_tasks_packed_this_time);
    for (int tid = 0; tid < total_parts_packed; tid++) {
      //    	if(ci_list_gpu[tid]->hydro.count==0||ci_list_gpu[tid]==NULL){
      //    		fprintf(stderr,"a cell has no particles or has a NULL
      //    value!\n");
      //    	}
      for (int i = 0; i < ci_list_gpu[tid]->hydro.count; i++) {
        fprintf(fgpu, "%f %f %f %f %f\n", ci_list_gpu[tid]->hydro.parts[i].x[0],
                ci_list_gpu[tid]->hydro.parts[i].x[1],
                ci_list_gpu[tid]->hydro.parts[i].x[2],
                ci_list_gpu[tid]->hydro.parts[i].rho,
                ci_list_gpu[tid]->hydro.parts[i].density.rho_dh);
      }
    }

    //    for(int tid=0;tid<tasks_done_cpu; tid++){
    //      for(int i=0; i<ci_list_cpu[tid]->hydro.count; i++){
    //        fprintf(fgpu, "%f %f %f %f %f\n",
    //        ci_list_cpu[tid]->hydro.parts[i].x[0],
    //        ci_list_cpu[tid]->hydro.parts[i].x[1],
    //        ci_list_cpu[tid]->hydro.parts[i].x[2],
    //        ci_list_cpu[tid]->hydro.parts[i].rho,
    //        ci_list_cpu[tid]->hydro.parts[i].density.rho_dh);
    //      }
    //    }
    time_for_density_cpu = 0.0;
    total_parts_done_cpu = 0;
    total_parts_packed = 0;
    tasks_done_gpu = 0;
    tasks_done_cpu = 0;

    packed = 0;
    density = 0;
    density_sub = 0;
    unpacked = 0;
    fprintf(stderr, "cpuid %i has %i packs self left\n", r->cpuid,
            sched->queues[r->qid].n_packs_self_left);
    /* Wait at the wait barrier. */
    //      swift_barrier_wait(&e->wait_barrier);
    //    if(r->cpuid==0)exit(0);
    cudaProfilerStop();
    //    exit(0);
  }
  // Free all data
  fprintf(stderr, "got to end of runner_main\n");
  cudaFree(d_tid_p);
  cudaFree(d_id);
  cudaFree(d_x_p);
  cudaFree(d_y_p);
  cudaFree(d_z_p);
  cudaFree(d_ux);
  cudaFree(d_uy);
  cudaFree(d_uz);
  cudaFree(d_a_hydrox);
  cudaFree(d_a_hydroy);
  cudaFree(d_a_hydroz);
  cudaFree(d_mass);
  cudaFree(d_h);
  cudaFree(d_u);
  cudaFree(d_u_dt);
  cudaFree(d_rho);
  cudaFree(d_SPH_sum);
  cudaFree(d_locx);
  cudaFree(d_locy);
  cudaFree(d_locz);
  cudaFree(d_widthx);
  cudaFree(d_widthy);
  cudaFree(d_widthz);
  cudaFree(d_h_max);
  cudaFree(d_count_p);
  cudaFree(d_wcount);
  cudaFree(d_wcount_dh);
  cudaFree(d_rho_dh);
  cudaFree(d_rot_ux);
  cudaFree(d_rot_uy);
  cudaFree(d_rot_uz);
  cudaFree(d_div_v);
  cudaFree(d_div_v_previous_step);
  cudaFree(d_alpha_visc);
  cudaFree(d_v_sig);
  cudaFree(d_laplace_u);
  cudaFree(d_alpha_diff);
  cudaFree(d_f);
  cudaFree(d_soundspeed);
  cudaFree(d_h_dt);
  cudaFree(d_balsara);
  cudaFree(d_pressure);
  cudaFree(d_alpha_visc_max_ngb);
  cudaFree(d_time_bin);
  cudaFree(d_wakeup);
  cudaFree(d_min_ngb_time_bin);
  cudaFree(d_to_be_synchronized);
  cudaFree(tid_p);
  cudaFree(id);
  cudaFree(mass);
  cudaFree(h);
  cudaFree(u);
  cudaFree(u_dt);
  cudaFree(rho);
  cudaFree(SPH_sum);
  cudaFree(x_p);
  cudaFree(y_p);
  cudaFree(z_p);
  cudaFree(ux);
  cudaFree(uy);
  cudaFree(uz);
  cudaFree(a_hydrox);
  cudaFree(a_hydroy);
  cudaFree(a_hydroz);
  cudaFree(locx);
  cudaFree(locy);
  cudaFree(locz);
  cudaFree(widthx);
  cudaFree(widthy);
  cudaFree(widthz);
  cudaFree(h_max);
  cudaFree(count_p);
  cudaFree(wcount);
  cudaFree(wcount_dh);
  cudaFree(rho_dh);
  cudaFree(rot_ux);
  cudaFree(rot_uy);
  cudaFree(rot_uz);
  cudaFree(div_v);
  cudaFree(div_v_previous_step);
  cudaFree(alpha_visc);
  cudaFree(v_sig);
  cudaFree(laplace_u);
  cudaFree(alpha_diff);
  cudaFree(f);
  cudaFree(soundspeed);
  cudaFree(h_dt);
  cudaFree(balsara);
  cudaFree(pressure);
  cudaFree(alpha_visc_max_ngb);
  cudaFree(time_bin);
  cudaFree(wakeup);
  cudaFree(min_ngb_time_bin);
  cudaFree(to_be_synchronized);
  cudaFree(partid_p);
  cudaFree(d_task_first_part);
  cudaFree(d_task_last_part);
  cudaFree(task_first_part);
  cudaFree(task_last_part);
  cudaFree(d_bundle_first_part);
  cudaFree(d_bundle_last_part);
  cudaFree(bundle_first_part);
  cudaFree(bundle_last_part);
  free(ci_list);

  /* Be kind, rewind. */
  return NULL;
}

#endif
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>

// uint64_t time_used ( ) {
//    struct rusage ru;
//    struct timeval t;
//    getrusage(RUSAGE_THREAD,&ru);
//    t = ru.ru_utime;
//    return (uint64_t) t.tv_sec*1000 + t.tv_usec/1000;
// }
