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

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "engine.h"
#include "gravity_io.h"
#include "infinity_wrapper.h"
#include "mpiuse.h"
#include "feedback.h"
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

#ifdef EXTRA_STAR_LOOPS

/* Import the stars prepare1 loop functions. */
#define FUNCTION prep1
#define FUNCTION_TASK_LOOP TASK_LOOP_STARS_PREP1
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the stars prepare2 loop functions. */
#define FUNCTION prep2
#define FUNCTION_TASK_LOOP TASK_LOOP_STARS_PREP2
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

#endif /* EXTRA_STAR_LOOPS */

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
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_INJECT
#include "runner_doiact_rt.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the RT gradient loop functions */
#define FUNCTION rt_gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_GRADIENT
#include "runner_doiact_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the RT transport (force) loop functions. */
#define FUNCTION rt_transport
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_TRANSPORT
#include "runner_doiact_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the sink compute formation loop functions. */
#define FUNCTION compute_formation
#define FUNCTION_TASK_LOOP TASK_LOOP_SINK_FORMATION
#include "runner_doiact_sinks.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the sink compute formation loop functions. */
#define FUNCTION accretion
#define FUNCTION_TASK_LOOP TASK_LOOP_SINK_ACCRETION
#include "runner_doiact_sinks.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the sink merger loop functions. */
#define FUNCTION merger
#define FUNCTION_TASK_LOOP TASK_LOOP_SINK_MERGER
#include "runner_doiact_sinks_merger.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 */
void *runner_main(void *data) {

  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  /* Main loop. */
  while (1) {

    /* Wait at the barrier. */
    engine_barrier(e);

    /* Can we go home yet? */
    if (e->step_props & engine_step_prop_done) break;

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;

    /* Loop while there are tasks... */
    while (1) {

      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, r->cpuid, prev);
        TIMER_TOC(timer_gettask);

        /* Did I get anything? */
        if (t == NULL) break;
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

      const ticks task_beg = getticks();
      /* Different types of tasks... */
      switch (t->type) {
        case task_type_self:
          if (t->subtype == task_subtype_density)
            runner_doself1_branch_density(r, ci);
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
          else if (t->subtype == task_subtype_rt_inject)
            runner_doself_branch_rt_inject(r, ci, 1);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_doself1_branch_rt_gradient(r, ci);
          else if (t->subtype == task_subtype_rt_transport)
            runner_doself2_branch_rt_transport(r, ci);
          else if (t->subtype == task_subtype_sink_compute_formation)
            runner_doself_branch_sinks_compute_formation(r, ci);
          else if (t->subtype == task_subtype_sink_accretion)
            runner_doself_branch_sinks_accretion(r, ci);
          else if (t->subtype == task_subtype_sink_merger)
            runner_doself_sinks_merger(r, ci);
          else
            error("Unknown/invalid task subtype (%s).",
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
          else if (t->subtype == task_subtype_rt_inject)
            runner_dopair_branch_rt_inject(r, ci, cj, 1);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_dopair1_branch_rt_gradient(r, ci, cj);
          else if (t->subtype == task_subtype_rt_transport)
            runner_dopair2_branch_rt_transport(r, ci, cj);
          else if (t->subtype == task_subtype_sink_compute_formation)
            runner_dopair_branch_sinks_compute_formation(r, ci, cj);
          else if (t->subtype == task_subtype_sink_accretion)
            runner_dopair_branch_sinks_accretion(r, ci, cj);
          else if (t->subtype == task_subtype_sink_merger)
            runner_do_sym_pair_sinks_merger(r, ci, cj);
          else
            error("Unknown/invalid task subtype (%s/%s).",
                  taskID_names[t->type], subtaskID_names[t->subtype]);
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
          else if (t->subtype == task_subtype_rt_inject)
            runner_dosub_self_rt_inject(r, ci, 1);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_dosub_self1_rt_gradient(r, ci, 1);
          else if (t->subtype == task_subtype_rt_transport)
            runner_dosub_self2_rt_transport(r, ci, 1);
          else if (t->subtype == task_subtype_sink_compute_formation)
            runner_dosub_self_sinks_compute_formation(r, ci, 1);
          else if (t->subtype == task_subtype_sink_accretion)
            runner_dosub_self_sinks_accretion(r, ci, 1);
          else if (t->subtype == task_subtype_sink_merger)
            runner_dosub_self_sinks_merger(r, ci);
          else
            error("Unknown/invalid task subtype (%s/%s).",
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
          else if (t->subtype == task_subtype_rt_inject)
            runner_dosub_pair_rt_inject(r, ci, cj, 1);
          else if (t->subtype == task_subtype_rt_gradient)
            runner_dosub_pair1_rt_gradient(r, ci, cj, 1);
          else if (t->subtype == task_subtype_rt_transport)
            runner_dosub_pair2_rt_transport(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_compute_formation)
            runner_dosub_pair_sinks_compute_formation(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_accretion)
            runner_dosub_pair_sinks_accretion(r, ci, cj, 1);
          else if (t->subtype == task_subtype_sink_merger)
            runner_dosub_pair_sinks_merger(r, ci, cj);
          else
            error("Unknown/invalid task subtype (%s/%s).",
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
#ifdef WITH_MPI
        case task_type_send:
          {
            if (t->flags != -1) {
              // XXX log actual time used in call, not how long we are queued.
              mpiuse_log_allocation(t->type, t->subtype, &t->buff, 1, t->win_size,
                                    cj->nodeID, t->flags, t->sendfull);
              
              /* Need space for the data and the headers. */
              scheduler_rdma_blocktype datasize =
                scheduler_rdma_toblocks(t->win_size) + scheduler_rdma_header_size;

              /* Access the registered memory for transferring this data. */
              scheduler_rdma_blocktype *dataptr =
                infinity_get_send_buffer(sched->send_handle, cj->nodeID);

              /* First element is marked as LOCKED, so only we can update. */
              dataptr[t->win_offset] = scheduler_rdma_locked;
              dataptr[t->win_offset + 1] = t->win_size;
              dataptr[t->win_offset + 2] = t->flags;
              dataptr[t->win_offset + 3] = engine_rank;

              // For this type the data should already be copied and ready to go.
              if (t->subtype != task_subtype_gpart && t->subtype != task_subtype_xv) {
                memcpy(&dataptr[t->win_offset + scheduler_rdma_header_size],
                       t->buff, t->win_size);
              }

#ifdef SWIFT_DEBUG_CHECKS
              if (e->verbose)
                message(
                        "Sending message to %d from %d subtype %d tag %zu size %zu"
                        " (cf %lld %zu) offsets %zu",
                        cj->nodeID, ci->nodeID, t->subtype, dataptr[t->win_offset+2],
                        dataptr[t->win_offset+1], t->flags, t->win_size, t->win_offset);
#endif
              infinity_send_data(sched->send_handle, cj->nodeID,
                                 scheduler_rdma_tobytes(datasize),
                                 scheduler_rdma_tobytes(t->win_offset),
                                 scheduler_rdma_tobytes(t->win_offset));
#ifdef SWIFT_DEBUG_CHECKS
              if (e->verbose) {
                message(
                        "Sent message to %d subtype %d tag %zu size %zu offset %zu"
                        " (cf %lld %zu)",
                        cj->nodeID, t->subtype, dataptr[t->win_offset+2],
                        dataptr[t->win_offset+1], t->win_offset, t->flags,
                        t->win_size);
              }
#endif
            }
          }

          if (t->subtype == task_subtype_subgpart) {

            /* Get the send buffer for this remote and copy our data chunk. */
            char *rdmaptr = infinity_get_send_buffer(sched->send_handle, cj->nodeID);
            size_t offset = scheduler_rdma_tobytes(scheduler_rdma_header_size) +
              scheduler_rdma_tobytes(t->main_task->win_offset);

            if (t->sendfull) {
              //message("sending sendfull");
              offset += (t->sub_offset * sizeof(struct gpart));
              memcpy(rdmaptr + offset, t->buff, t->sub_size * sizeof(struct gpart));
            } else {
              //message("sending not sendfull");

              /* Only sending position and mass, so need to pack them into
               * suitable strides. */
              offset += (t->sub_offset * sizeof(struct reduced_gpart));

              struct gpart *gp = (struct gpart*) t->buff;
              struct reduced_gpart *rgp = (struct reduced_gpart *) (rdmaptr + offset);
              for (int i = 0; i < t->sub_size; i++) {
                memcpy(rgp[i].x, gp[i].x, sizeof(struct reduced_gpart));
                //rgp[i].x[0] = gp[i].x[0];
                //rgp[i].x[1] = gp[i].x[1];
                //rgp[i].x[2] = gp[i].x[2];
                //rgp[i].mass = gp[i].mass;
#ifdef SWIFT_DEBUG_CHECKS
                rgp[i].ti_drift = gp[i].ti_drift;
                rgp[i].id = i;
#endif
              }
            }

          }
          else if (t->subtype == task_subtype_subxv) {

            /* Get the send buffer for this remote and copy our data chunk. */
            char *dataptr = infinity_get_send_buffer(sched->send_handle, cj->nodeID);
            size_t offset = scheduler_rdma_tobytes(scheduler_rdma_header_size) +
              scheduler_rdma_tobytes(t->main_task->win_offset) +
              (t->sub_offset * sizeof(struct part));
            memcpy(dataptr + offset, t->buff, t->sub_size * sizeof(struct part));
            //message("running subxv copy");
          }
          else if (t->subtype == task_subtype_tend_part) {
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
          } else if (t->subtype == task_subtype_limiter) {
            free(t->buff);
          }

          /* And log deactivation, if logging enabled. */
          if (t->flags != -1) {
            mpiuse_log_allocation(t->type, t->subtype, &t->buff, 0, 0, 0, 0, t->sendfull);
          }

          /* Only sendfull once, reduced_gparts parts next time. */
          t->sendfull = 0;

          break;

        case task_type_recv:
          /* Ready to process. So copy to local buffer, unless we are doing
           * this in parts. */
          if (t->flags != -1) {
            if (t->subtype != task_subtype_gpart && t->subtype != task_subtype_xv) {
              if (t->rdmabuff == NULL) error("No RDMA data yet!");
              memcpy(t->buff, t->rdmabuff, t->win_size);
            }
          }

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

          } else if (t->subtype == task_subtype_doxv) {
            runner_do_recv_part(r, ci, 1, 1);
          } else if (t->subtype == task_subtype_xv) {
            /* Does nothing, just unlocks subxv and doxv recvs, already
             * accepted remove data send. */

          } else if (t->subtype == task_subtype_subxv) {

            /* Only move our part of the data. */
            char *ptr1 = ((char *)t->buff);
            char *ptr2 = ((char *)t->main_task->rdmabuff) + (t->sub_offset * sizeof(struct part));
            memcpy(ptr1, ptr2, t->sub_size * sizeof(struct part));

          } else if (t->subtype == task_subtype_rho) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_gradient) {
            runner_do_recv_part(r, ci, 0, 1);
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

          } else if (t->subtype == task_subtype_dogpart) {
            /* Handle complete data movement. */
            //message("completing recv with the dogpart");
            runner_do_recv_gpart(r, ci, 1);

          } else if (t->subtype == task_subtype_gpart) {
            /* Does nothing, just unlocks subgpart and dogpart recvs, already
             * accepted remove data send. */

          } else if (t->subtype == task_subtype_subgpart) {

            char *rdmaptr = t->main_task->rdmabuff;

            if (t->sendfull) {
              //message("expecting sendfull");
              /* Being sent full gparts, but only our chunk... */
              char *dataptr = rdmaptr + (t->sub_offset * sizeof(struct gpart));
              memcpy(t->buff, dataptr, t->sub_size * sizeof(struct gpart));

            } else {
              //message("expecting not sendfull");
              /* Expecting partial updates. */
              struct gpart *gp = (struct gpart*) t->buff;
              struct reduced_gpart *rgp = (struct reduced_gpart *)
                (rdmaptr + (t->sub_offset * sizeof(struct reduced_gpart)));
              for (int i = 0; i < t->sub_size; i++) {
                memcpy(gp[i].x, rgp[i].x, sizeof(struct reduced_gpart));
                //gp[i].x[0] = rgp[i].x[0];
                //gp[i].x[1] = rgp[i].x[1];
                //gp[i].x[2] = rgp[i].x[2];
                //gp[i].mass = rgp[i].mass;
#ifdef SWIFT_DEBUG_CHECKS
                gp[i].ti_drift = rgp[i].ti_drift;
                message("ti_drift = %lld", gp[i].ti_drift);
                if (rgp[i].id != i)
                  message("index misaligned: %d != %d", rgp[i].id, i);
#endif
              }
            }

          } else if (t->subtype == task_subtype_spart_density) {
            runner_do_recv_spart(r, ci, 1, 1);
          } else if (t->subtype == task_subtype_part_prep1) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_spart_prep2) {
            runner_do_recv_spart(r, ci, 0, 1);
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

          /* And log deactivation, if logging enabled. */
          if (t->flags != -1) {
            mpiuse_log_allocation(t->type, t->subtype, &t->buff, 0, 0, 0, 0, t->sendfull);
          }

          /* Not expecting full sends again until next rebuild. */
          t->sendfull = 0;

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
          runner_do_fof_self(r, t->ci, 1);
          break;
        case task_type_fof_pair:
          runner_do_fof_pair(r, t->ci, t->cj, 1);
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
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}

ticks runner_get_active_time(const struct runner *restrict r) {
  return r->active_time;
}

void runner_reset_active_time(struct runner *restrict r) { r->active_time = 0; }
