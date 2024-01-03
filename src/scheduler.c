/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2016 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Some standard headers. */
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "scheduler.h"

/* Local headers. */
#include "atomic.h"
#include "cycle.h"
#include "engine.h"
#include "error.h"
#include "intrinsics.h"
#include "kernel_hydro.h"
#include "memuse.h"
#include "mpiuse.h"
#include "queue.h"
#include "sort_part.h"
#include "space.h"
#include "space_getsid.h"
#include "task.h"
#include "threadpool.h"
#include "timers.h"
#include "version.h"

/**
 * @brief Re-set the list of active tasks.
 */
void scheduler_clear_active(struct scheduler *s) { s->active_count = 0; }

/**
 * @brief Increase the space available for unlocks. Only call when
 *        current index == s->size_unlock;
 */
static void scheduler_extend_unlocks(struct scheduler *s) {
  /* Allocate the new buffer. */
  const int size_unlocks_new = s->size_unlocks * 2;
  struct task **unlocks_new = (struct task **)swift_malloc(
      "unlocks", sizeof(struct task *) * size_unlocks_new);
  int *unlock_ind_new =
      (int *)swift_malloc("unlock_ind", sizeof(int) * size_unlocks_new);
  if (unlocks_new == NULL || unlock_ind_new == NULL)
    error("Failed to re-allocate unlocks.");

  /* Wait for all writes to the old buffer to complete. */
  while (s->completed_unlock_writes < s->size_unlocks)
    ;

  /* Copy the buffers. */
  memcpy(unlocks_new, s->unlocks, sizeof(struct task *) * s->size_unlocks);
  memcpy(unlock_ind_new, s->unlock_ind, sizeof(int) * s->size_unlocks);
  swift_free("unlocks", s->unlocks);
  swift_free("unlock_ind", s->unlock_ind);
  s->unlocks = unlocks_new;
  s->unlock_ind = unlock_ind_new;

  /* Publish the new buffer size. */
  s->size_unlocks = size_unlocks_new;
}

/**
 * @brief Add an unlock_task to the given task.
 *
 * @param s The #scheduler.
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.

 */
void scheduler_addunlock(struct scheduler *s, struct task *ta,
                         struct task *tb) {
#ifdef SWIFT_DEBUG_CHECKS
  if (ta == NULL) error("Unlocking task is NULL.");
  if (tb == NULL) error("Unlocked task is NULL.");
#endif

  /* Get an index at which to store this unlock. */
  const int ind = atomic_inc(&s->nr_unlocks);

  /* Does the buffer need to be grown? */
  if (ind == s->size_unlocks) scheduler_extend_unlocks(s);

#ifdef SWIFT_DEBUG_CHECKS
  if (ind > s->size_unlocks * 2)
    message("unlocks guard enabled: %d / %d", ind, s->size_unlocks);
#endif

  /* Wait for there to actually be space at my index. */
  while (ind > s->size_unlocks)
    ;

  /* Guard against case when more than (old) s->size_unlocks unlocks
   * are now pending. */
  if (ind == s->size_unlocks) scheduler_extend_unlocks(s);

  /* Write the unlock to the scheduler. */
  s->unlocks[ind] = tb;
  s->unlock_ind[ind] = ta - s->tasks;
  atomic_inc(&s->completed_unlock_writes);
}

/* Conservative number of dependencies per task type */
#define MAX_NUMBER_DEP 128

/**
 * @brief Describe the level at which the task are done.
 * WARNING: the order is supposed to be sorted from the root
 * to the leaf.
 */
enum task_dependency_level {
  task_dependency_level_top = 0,
  task_dependency_level_super,
  task_dependency_level_super_hydro,
  task_dependency_level_super_grav,
  task_dependency_level_none,
};

/**
 * @brief Informations about all the task dependencies of
 *   a single task.
 */
struct task_dependency {
  /* Main task */
  /* ID of the task */
  int type_in;

  /* ID of the subtask */
  int subtype_in;

  /* Is the task implicit */
  int implicit_in;

  /* Is the taks_in at the top level? */
  int task_in_is_top;

  /* Is the taks_in at the grav.super level? */
  int task_in_is_grav_super;

  /* Is the taks_in at the hydro.super level? */
  int task_in_is_hydro_super;

  /* Dependent task */
  /* ID of the dependent task */
  int type_out[MAX_NUMBER_DEP];

  /* ID of the dependent subtask */
  int subtype_out[MAX_NUMBER_DEP];

  /* Is the dependent task implicit */
  int implicit_out[MAX_NUMBER_DEP];

  /* Is the taks_out at the top level? */
  int task_out_is_top[MAX_NUMBER_DEP];

  /* Is the taks_out at the grav.super level? */
  int task_out_is_grav_super[MAX_NUMBER_DEP];

  /* Is the taks_out at the hydro.super level? */
  int task_out_is_hydro_super[MAX_NUMBER_DEP];

  /* Statistics */
  /* number of link between the two task type */
  int number_link[MAX_NUMBER_DEP];

  /* number of ranks having this relation */
  int number_rank[MAX_NUMBER_DEP];
};

#ifdef WITH_MPI

/**
 * @brief Define the #task_dependency for MPI
 *
 * @param tstype The MPI_Datatype to initialize
 */
void task_dependency_define(MPI_Datatype *tstype) {
  /* Define the variables */
  const int count = 14;
  int blocklens[count];
  MPI_Datatype types[count];
  MPI_Aint disps[count];

  /* all the type are int */
  for (int i = 0; i < count; i++) {
    types[i] = MPI_INT;
  }

  /* Task in */
  disps[0] = offsetof(struct task_dependency, type_in);
  blocklens[0] = 1;
  disps[1] = offsetof(struct task_dependency, subtype_in);
  blocklens[1] = 1;
  disps[2] = offsetof(struct task_dependency, implicit_in);
  blocklens[2] = 1;
  disps[3] = offsetof(struct task_dependency, task_in_is_top);
  blocklens[3] = 1;
  disps[4] = offsetof(struct task_dependency, task_in_is_hydro_super);
  blocklens[4] = 1;
  disps[5] = offsetof(struct task_dependency, task_in_is_grav_super);
  blocklens[5] = 1;

  /* Task out */
  disps[6] = offsetof(struct task_dependency, type_out);
  blocklens[6] = MAX_NUMBER_DEP;
  disps[7] = offsetof(struct task_dependency, subtype_out);
  blocklens[7] = MAX_NUMBER_DEP;
  disps[8] = offsetof(struct task_dependency, implicit_out);
  blocklens[8] = MAX_NUMBER_DEP;
  disps[9] = offsetof(struct task_dependency, task_out_is_top);
  blocklens[9] = MAX_NUMBER_DEP;
  disps[10] = offsetof(struct task_dependency, task_out_is_hydro_super);
  blocklens[10] = MAX_NUMBER_DEP;
  disps[11] = offsetof(struct task_dependency, task_out_is_grav_super);
  blocklens[11] = MAX_NUMBER_DEP;

  /* statistics */
  disps[12] = offsetof(struct task_dependency, number_link);
  blocklens[12] = MAX_NUMBER_DEP;
  disps[13] = offsetof(struct task_dependency, number_rank);
  blocklens[13] = MAX_NUMBER_DEP;

  /* define it for MPI */
  MPI_Type_create_struct(count, blocklens, disps, types, tstype);
  MPI_Type_commit(tstype);
}

/**
 * @brief Sum operator of #task_dependency for MPI
 *
 * @param in_p The #task_dependency to add
 * @param out_p The #task_dependency where in_p is added
 * @param len The length of the arrays
 * @param type The MPI datatype
 */
void task_dependency_sum(void *in_p, void *out_p, int *len,
                         MPI_Datatype *type) {
  /* change pointer type */
  struct task_dependency *in = (struct task_dependency *)in_p;
  struct task_dependency *out = (struct task_dependency *)out_p;

  /* Loop over all the current objects */
  for (int i = 0; i < *len; i++) {
    /* loop over all the object set in invals */
    for (int j = 0; j < MAX_NUMBER_DEP; j++) {

      /* Have we reached the end of the links? */
      if (in[i].number_link[j] == -1) {
        break;
      }

      /* get a few variables */
      int tb_type = in[i].type_out[j];
      int tb_subtype = in[i].subtype_out[j];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check tasks */
      if (tb_type >= task_type_count) {
        error("Unknown task type %i", tb_type);
      }

      if (tb_subtype >= task_subtype_count) {
        error("Unknown subtask type %i", tb_subtype);
      }
#endif

      /* find the corresponding id */
      int k = 0;
      while (k < MAX_NUMBER_DEP) {
        /* have we reached the end of the links? */
        if (out[i].number_link[k] == -1) {
          /* reset the counter in order to be safe */
          out[i].number_link[k] = 0;
          out[i].number_rank[k] = 0;

          /* set the relation */
          out[i].type_in = in[i].type_in;
          out[i].subtype_in = in[i].subtype_in;
          out[i].implicit_in = in[i].implicit_in;

          out[i].type_out[k] = in[i].type_out[j];
          out[i].subtype_out[k] = in[i].subtype_out[j];
          out[i].implicit_out[k] = in[i].implicit_out[j];
          break;
        }

        /* do we have the same relation? */
        if (out[i].type_out[k] == tb_type &&
            out[i].subtype_out[k] == tb_subtype) {
          break;
        }

        k++;
      }

      /* Check if we are still in the memory */
      if (k == MAX_NUMBER_DEP) {
        error("Not enough memory, please increase MAX_NUMBER_DEP");
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check if correct relation */
      if (out[i].type_in != in[i].type_in ||
          out[i].subtype_in != in[i].subtype_in ||
          out[i].implicit_in != in[i].implicit_in ||
          out[i].type_out[k] != in[i].type_out[j] ||
          out[i].subtype_out[k] != in[i].subtype_out[j] ||
          out[i].implicit_out[k] != in[i].implicit_out[j]) {
        error("Tasks do not correspond");
      }
#endif
      /* sum the contributions */
      out[i].number_link[k] += in[i].number_link[j];
      out[i].number_rank[k] += in[i].number_rank[j];

      /* Get the task in level */
      out[i].task_in_is_top = min(out[i].task_in_is_top, in[i].task_in_is_top);
      out[i].task_in_is_hydro_super =
          min(out[i].task_in_is_hydro_super, in[i].task_in_is_hydro_super);
      out[i].task_in_is_grav_super =
          min(out[i].task_in_is_grav_super, in[i].task_in_is_grav_super);

      /* Get the task out level */
      out[i].task_out_is_top[j] =
          min(out[i].task_out_is_top[j], in[i].task_out_is_top[j]);
      out[i].task_out_is_hydro_super[j] = min(out[i].task_out_is_hydro_super[j],
                                              in[i].task_out_is_hydro_super[j]);
      out[i].task_out_is_grav_super[j] = min(out[i].task_out_is_grav_super[j],
                                             in[i].task_out_is_grav_super[j]);
    }
  }

  return;
}

#endif  // WITH_MPI

/**
 * @brief Write a csv file with the task dependencies.
 *
 * Run plot_task_dependencies.py for an example of how to use it
 * to generate the figure.
 *
 * @param s The #scheduler we are working in.
 * @param verbose Are we verbose about this?
 * @param step The current step number.
 */
void scheduler_write_dependencies(struct scheduler *s, int verbose, int step) {
  const ticks tic = getticks();

  /* Number of possible relations between tasks */
  const int nber_tasks = task_type_count * task_subtype_count;

  /* To get the table for a task:
   * ind = (ta * task_subtype_count + sa)
   * where ta is the value of task_type and sa is the value of
   * task_subtype  */
  struct task_dependency *task_dep = (struct task_dependency *)malloc(
      nber_tasks * sizeof(struct task_dependency));
  /* keep track of whether a task exists in this run */
  int *task_exists = (int *)malloc(nber_tasks * sizeof(int));
  /* keep track of whether a task has a dependency or an unlock,
   * and hence will be drawn in the task graph */
  int *task_has_deps = (int *)malloc(nber_tasks * sizeof(int));

  /* Special marker for tasks with no dependencies. */
  const int no_dependency = -3;

  if (task_dep == NULL)
    error("Error allocating memory for task-dependency graph (table).");

  /* Reset counter */
  for (int i = 0; i < nber_tasks; i++) {
    /* Assume that the tasks are at all levels and correct later. */
    task_dep[i].task_in_is_top = 1;
    task_dep[i].task_in_is_grav_super = 1;
    task_dep[i].task_in_is_hydro_super = 1;
    const int tt = i / task_subtype_count;
    const int tst = i % task_subtype_count;
    task_dep[i].type_in = tt;
    task_dep[i].subtype_in = tst;

    for (int j = 0; j < MAX_NUMBER_DEP; j++) {
      /* Use number_link as indicator of the existance of a relation */
      task_dep[i].number_link[j] = -1;

      /* Assume that the tasks are at all levels and correct later. */
      task_dep[i].task_out_is_top[j] = 1;
      task_dep[i].task_out_is_grav_super[j] = 1;
      task_dep[i].task_out_is_hydro_super[j] = 1;
    }
    task_exists[i] = 0;
    task_has_deps[i] = 0;
  }

  /* loop over all tasks */
  for (int i = 0; i < s->nr_tasks; i++) {
    const struct task *ta = &s->tasks[i];

    /* Are we using this task?
     * For the 0-step, we wish to show all the tasks (even the inactives). */
    if (step != 0 && ta->skip) continue;

    /* Current index */
    const int ind = ta->type * task_subtype_count + ta->subtype;

    struct task_dependency *cur = &task_dep[ind];
    task_exists[ind]++;

#ifdef SWIFT_DEBUG_CHECKS
    if (cur->type_in != ta->type)
      error("wrong indexing for task %d: Expect type %d got %d", i,
            cur->type_in, ta->type);
    if (cur->subtype_in != ta->subtype)
      error("wrong indexing for task %d: Expect subtype %d got %d", i,
            cur->subtype_in, ta->subtype);
#endif
    /* Is ta implicit? */
    cur->implicit_in = ta->implicit;

    /* Set the task level. */
    const struct cell *ci = ta->ci;
    const struct cell *cj = ta->cj;
    const int is_ci_top = ci == NULL || (ci != NULL && ci == ci->top);
    const int is_cj_top = cj == NULL || (cj != NULL && cj == cj->top);
    const int is_hydro_super =
        (cj == NULL || (cj != NULL && cj == cj->hydro.super)) &&
        (ci == NULL || (ci != NULL && ci == ci->hydro.super));
    const int is_grav_super =
        (cj == NULL || (cj != NULL && cj == cj->grav.super)) &&
        (ci == NULL || (ci != NULL && ci == ci->grav.super));

    /* Are we dealing with a task at the top level? */
    if (!(is_ci_top && is_cj_top)) {
      cur->task_in_is_top = 0;
    }
    /* At the hydro level? */
    if (!is_hydro_super) {
      cur->task_in_is_hydro_super = 0;
    }
    /* At the gravity level? */
    if (!is_grav_super) {
      cur->task_in_is_grav_super = 0;
    }

    /* If this task unlocks nothing, make a note of it. */
    if (ta->nr_unlock_tasks == 0) {
      int k = 0;
      while (k < MAX_NUMBER_DEP) {
        /* not written yet */
        if (cur->number_link[k] == -1) {
          cur->type_out[k] = no_dependency;
          cur->subtype_out[k] = no_dependency;
          cur->implicit_out[k] = no_dependency;

          /* statistics */
          cur->number_link[k] = 0;
          cur->number_rank[k] = 1;

          /* Are we dealing with a task at the top level? */
          cur->task_out_is_top[k] = no_dependency;
          cur->task_out_is_hydro_super[k] = no_dependency;
          cur->task_out_is_grav_super[k] = no_dependency;

          break;
        }

        /* already written */
        if (cur->type_out[k] == no_dependency &&
            cur->subtype_out[k] == no_dependency) {
          break;
        }
        k += 1;
      }

      /* if this task unlocks nothing, we have nothing left to do.
       * go to next. */
      continue;
    }

    /* This task unlocks stuff, so check the dependencies */
    for (int j = 0; j < ta->nr_unlock_tasks; j++) {
      const struct task *tb = ta->unlock_tasks[j];

      /* Are we using this task?
       * For the 0-step, we wish to show all the tasks (even the inactive). */
      if (step != 0 && tb->skip) continue;

      int indj = tb->type * task_subtype_count + tb->subtype;

#ifdef SWIFT_DEBUG_CHECKS
      const struct task_dependency *target = &task_dep[indj];
      if (target->type_in != tb->type)
        error("wrong indexing for task %d: Expect type %d got %d", i,
              target->type_in, tb->type);
      if (target->subtype_in != tb->subtype)
        error("wrong indexing for task %d: Expect subtype %d got %d", i,
              target->subtype_in, tb->subtype);
#endif
      task_exists[indj]++;

      const struct cell *ci_b = tb->ci;
      const struct cell *cj_b = tb->cj;
      const int is_ci_b_top =
          ci_b == NULL || (ci_b != NULL && ci_b == ci_b->top);
      const int is_cj_b_top =
          cj_b == NULL || (cj_b != NULL && cj_b == cj_b->top);
      const int is_b_hydro_super =
          (cj_b == NULL || (cj_b != NULL && cj_b == cj_b->hydro.super)) &&
          (ci_b == NULL || (ci_b != NULL && ci_b == ci_b->hydro.super));
      const int is_b_grav_super =
          (cj_b == NULL || (cj_b != NULL && cj_b == cj_b->grav.super)) &&
          (ci_b == NULL || (ci_b != NULL && ci_b == ci_b->grav.super));

      int k = 0;
      while (k < MAX_NUMBER_DEP) {
        /* not written yet */
        if (cur->number_link[k] == -1) {
          /* set tb */
          cur->type_out[k] = tb->type;
          cur->subtype_out[k] = tb->subtype;
          cur->implicit_out[k] = tb->implicit;

          /* statistics */
          cur->number_link[k] = 1;
          cur->number_rank[k] = 1;

          /* Are we dealing with a task at the top level? */
          if (!(is_ci_b_top && is_cj_b_top)) {
            cur->task_out_is_top[k] = 0;
          }
          /* At the hydro level? */
          if (!is_b_hydro_super) {
            cur->task_out_is_hydro_super[k] = 0;
          }
          /* At the gravity level? */
          if (!is_b_grav_super) {
            cur->task_out_is_grav_super[k] = 0;
          }

          break;
        }

        /* already written */
        if (cur->type_out[k] == tb->type &&
            cur->subtype_out[k] == tb->subtype) {

          /* Increase the number of link. */
          cur->number_link[k] += 1;

          /* Are we dealing with a task at the top level? */
          if (!(is_ci_b_top && is_cj_b_top)) {
            cur->task_out_is_top[k] = 0;
          }
          /* At the hydro level? */
          if (!is_b_hydro_super) {
            cur->task_out_is_hydro_super[k] = 0;
          }
          /* At the gravity level? */
          if (!is_b_grav_super) {
            cur->task_out_is_grav_super[k] = 0;
          }
          break;
        }

        k += 1;
      }

      /* MAX_NUMBER_DEP is too small */
      if (k == MAX_NUMBER_DEP)
        error("Not enough memory, please increase MAX_NUMBER_DEP");
    }
  }

#ifdef WITH_MPI
  /* create MPI operator */
  MPI_Datatype dependency_data_type;
  task_dependency_define(&dependency_data_type);

  MPI_Op dependency_sum;
  MPI_Op_create(task_dependency_sum, /* commute */ 1, &dependency_sum);

  /* create recv buffer */
  struct task_dependency *recv = NULL;
  int *recv_exists = NULL;

  if (s->nodeID == 0) {
    recv = (struct task_dependency *)malloc(nber_tasks *
                                            sizeof(struct task_dependency));
    recv_exists = (int *)malloc(nber_tasks * sizeof(int));

    /* reset counter */
    for (int i = 0; i < nber_tasks; i++) {
      for (int j = 0; j < MAX_NUMBER_DEP; j++) {
        /* Use number_link as indicator of the existance of a relation */
        recv[i].number_link[j] = -1;
      }
      recv_exists[i] = 0;
    }
  }

  /* Do the reduction */
  int test = MPI_Reduce(task_dep, recv, nber_tasks, dependency_data_type,
                        dependency_sum, 0, MPI_COMM_WORLD);
  if (test != MPI_SUCCESS) error("MPI reduce failed");

  test = MPI_Reduce(task_exists, recv_exists, nber_tasks, MPI_INT, MPI_SUM, 0,
                    MPI_COMM_WORLD);
  if (test != MPI_SUCCESS) error("MPI reduce failed");

  /* free some memory */
  if (s->nodeID == 0) {
    free(task_dep);
    task_dep = recv;
    free(task_exists);
    task_exists = recv_exists;
  }
#endif

  if (s->nodeID == 0) {
    /* Create file */
    char filename[50];
    sprintf(filename, "dependency_graph_%i.csv", step);
    FILE *f = fopen(filename, "w");
    if (f == NULL) error("Error opening dependency graph file.");

    /* Write header */
    fprintf(f, "# %s\n", git_revision());
    fprintf(
        f,
        "task_in,task_out,implicit_in,implicit_out,mpi_in,mpi_out,cluster_in,"
        "cluster_out,number_link,number_rank,task_in_is_top,task_in_is_hydro_"
        "super,task_in_is_grav_super,task_out_is_top,task_out_is_hydro_super,"
        "task_out_is_grav_super,cell_has_active_task\n");

    for (int i = 0; i < nber_tasks; i++) {
      for (int j = 0; j < MAX_NUMBER_DEP; j++) {
        /* Does this link exists */
        if (task_dep[i].number_link[j] == -1) continue;
        /* Don't write tasks without dependencies (yet) */
        if (task_dep[i].type_out[j] == no_dependency) continue;

        /* Define a few variables */
        const int ta_type = task_dep[i].type_in;
        const int ta_subtype = task_dep[i].subtype_in;
        const int ta_implicit = task_dep[i].implicit_in;

        const int tb_type = task_dep[i].type_out[j];
        const int tb_subtype = task_dep[i].subtype_out[j];
        const int tb_implicit = task_dep[i].implicit_out[j];

        const int count = task_dep[i].number_link[j];
        const int number_rank = task_dep[i].number_rank[j];

        const int task_in_is_top = task_dep[i].task_in_is_top;
        const int task_in_is_grav_super = task_dep[i].task_in_is_grav_super;
        const int task_in_is_hydro_super = task_dep[i].task_in_is_hydro_super;

        const int task_out_is_top = task_dep[i].task_out_is_top[j];
        const int task_out_is_grav_super =
            task_dep[i].task_out_is_grav_super[j];
        const int task_out_is_hydro_super =
            task_dep[i].task_out_is_hydro_super[j];

        /* text to write */
        char ta_name[200];
        char tb_name[200];

        /* take note that these tasks have dependencies and unlocks */
        task_has_deps[i]++;
        int indj = tb_type * task_subtype_count + tb_subtype;
        task_has_deps[indj]++;

        /* construct line */
        task_get_full_name(ta_type, ta_subtype, ta_name);
        task_get_full_name(tb_type, tb_subtype, tb_name);

        /* Check if MPI */
        int ta_mpi = 0;
        if (ta_type == task_type_send || ta_type == task_type_recv) ta_mpi = 1;

        int tb_mpi = 0;
        if (tb_type == task_type_send || tb_type == task_type_recv) tb_mpi = 1;

        /* Get group name */
        char ta_cluster[20];
        char tb_cluster[20];
        task_get_group_name(ta_type, ta_subtype, ta_cluster);
        task_get_group_name(tb_type, tb_subtype, tb_cluster);

        fprintf(f, "%s,%s,%d,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
                ta_name, tb_name, ta_implicit, tb_implicit, ta_mpi, tb_mpi,
                ta_cluster, tb_cluster, count, number_rank, task_in_is_top,
                task_in_is_hydro_super, task_in_is_grav_super, task_out_is_top,
                task_out_is_hydro_super, task_out_is_grav_super,
                /*cell_has_active_task=*/1);
      }
    }

    /* Now write the tasks without dependencies */
    for (int i = 0; i < nber_tasks; i++) {

      /* There may be several tasks that don't unlock anything,
       * e.g. timestep_collect or kick1 tasks. Those are covered
       * in the graph by them being unlocked by some other task.
       * If a task however doesn't unlock anything, nor is unlocked
       * by any other task, it needs special treatment, which it
       * receives now. The condition to be written down is a) the
       * task must exist, i.e. we must have encountered it in the
       * list of tasks, and b) it must not have been unlocked by
       * anyting. */

      if (task_exists[i] && !task_has_deps[i]) {

        /* Define a few variables */
        const int ta_type = task_dep[i].type_in;
        const int ta_subtype = task_dep[i].subtype_in;
        const int ta_implicit = task_dep[i].implicit_in;

        const int tb_implicit = 0;

        const int count = 0;
        const int number_rank = 1;

        const int task_in_is_top = task_dep[i].task_in_is_top;
        const int task_in_is_grav_super = task_dep[i].task_in_is_grav_super;
        const int task_in_is_hydro_super = task_dep[i].task_in_is_hydro_super;

        const int task_out_is_top = -1;
        const int task_out_is_grav_super = -1;
        const int task_out_is_hydro_super = -1;

        /* text to write */
        char ta_name[200];
        task_get_full_name(ta_type, ta_subtype, ta_name);
        char *tb_name = "task_unlocks_nothing\0";

        /* Check if MPI */
        int ta_mpi = 0;
        if (ta_type == task_type_send || ta_type == task_type_recv) ta_mpi = 1;

        int tb_mpi = 0;

        /* Get group name */
        char ta_cluster[20];
        task_get_group_name(ta_type, ta_subtype, ta_cluster);
        char *tb_cluster = "None\0";

        fprintf(f, "%s,%s,%d,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
                ta_name, tb_name, ta_implicit, tb_implicit, ta_mpi, tb_mpi,
                ta_cluster, tb_cluster, count, number_rank, task_in_is_top,
                task_in_is_hydro_super, task_in_is_grav_super, task_out_is_top,
                task_out_is_hydro_super, task_out_is_grav_super,
                /*cell_has_active_task=*/1);
      }
    }

    /* Close the file */
    fclose(f);
  }

#if defined(SWIFT_DEBUG_CHECKS)
  /* Check if we have the correct number of dependencies. */
  if (step == 0) {
    int count_total = 0;
    for (int i = 0; i < nber_tasks; i++) {
      for (int j = 0; j < MAX_NUMBER_DEP; j++) {
        if (task_dep[i].number_link[j] != -1)
          count_total += task_dep[i].number_link[j];
      }
    }

    /* Get the number of unlocks from all the ranks */
    int nr_unlocks = s->nr_unlocks;
#ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &nr_unlocks, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    if (s->nodeID == 0 && count_total != nr_unlocks) {
      error("Not all the dependencies were found: %i != %i", count_total,
            nr_unlocks);
    }
  }
#endif

  /* Be clean */
  free(task_dep);
  free(task_exists);
  free(task_has_deps);
#ifdef WITH_MPI
  MPI_Type_free(&dependency_data_type);
  MPI_Op_free(&dependency_sum);
#endif

  if (verbose)
    message("Printing task graph took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

/**
 * @brief Write a csv file with the task dependencies for a single cell.
 *
 * Run plot_task_dependencies.py for an example of how to use it
 * to generate the figure.
 *
 * @param s The #scheduler we are working in.
 * @param verbose Are we verbose about this?
 * @param step The current step number.
 */
void scheduler_write_cell_dependencies(struct scheduler *s, int verbose,
                                       int step) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)

  const ticks tic = getticks();

  const long long cellID = s->dependency_graph_cellID;
  if (cellID == 0LL) return;

  /* Number of possible relations between tasks */
  const int nber_tasks = task_type_count * task_subtype_count;

  /* To get the table for a task:
   * ind = (ta * task_subtype_count + sa)
   * where ta is the value of task_type and sa is the value of
   * task_subtype  */
  struct task_dependency *task_dep = (struct task_dependency *)malloc(
      nber_tasks * sizeof(struct task_dependency));

  /* Keep track whether the requested cell is also involved in the
   * dependency */
  int cell_involved[nber_tasks][MAX_NUMBER_DEP];

  if (task_dep == NULL)
    error("Error allocating memory for task-dependency graph (table).");

  /* Reset counter */
  for (int i = 0; i < nber_tasks; i++) {
    /* Assume that the tasks are at all levels and correct later. */
    task_dep[i].task_in_is_top = 1;
    task_dep[i].task_in_is_grav_super = 1;
    task_dep[i].task_in_is_hydro_super = 1;
    const int tt = i / task_subtype_count;
    const int tst = i % task_subtype_count;
    task_dep[i].type_in = tt;
    task_dep[i].subtype_in = tst;

    for (int j = 0; j < MAX_NUMBER_DEP; j++) {
      /* Use number_link as indicator of the existance of a relation */
      task_dep[i].number_link[j] = -1;

      /* Assume that the tasks are at all levels and correct later. */
      task_dep[i].task_out_is_top[j] = 1;
      task_dep[i].task_out_is_grav_super[j] = 1;
      task_dep[i].task_out_is_hydro_super[j] = 1;
      cell_involved[i][j] = 0;
    }
  }

  /* loop over all tasks */
  int local_count = 0;
  for (int i = 0; i < s->nr_tasks; i++) {
    const struct task *ta = &s->tasks[i];

    /* Are we using this task?
     * For the 0-step, we wish to show all the tasks (even the inactives). */
    if (step != 0 && ta->skip) continue;
    /* Note: task_type_none may have t->ci==NULL too */
    if (!(((ta->ci != NULL) && ta->ci->cellID == cellID) ||
          ((ta->cj != NULL) && ta->cj->cellID == cellID)))
      continue;

    /* Current index */
    const int ind = ta->type * task_subtype_count + ta->subtype;
    struct task_dependency *cur = &task_dep[ind];

#ifdef SWIFT_DEBUG_CHECKS
    if (cur->type_in != ta->type)
      error("wrong indexing for task %d: Expect type %d got %d", i,
            cur->type_in, ta->type);
    if (cur->subtype_in != ta->subtype)
      error("wrong indexing for task %d: Expect subtype %d got %d", i,
            cur->subtype_in, ta->subtype);
#endif
    /* Is ta implicit? */
    cur->implicit_in = ta->implicit;

    /* Set the task level. */
    const struct cell *ci = ta->ci;
    const struct cell *cj = ta->cj;
    const int is_ci_top = ci == NULL || (ci != NULL && ci == ci->top);
    const int is_cj_top = cj == NULL || (cj != NULL && cj == cj->top);
    const int is_hydro_super =
        (cj == NULL || (cj != NULL && cj == cj->hydro.super)) &&
        (ci == NULL || (ci != NULL && ci == ci->hydro.super));
    const int is_grav_super =
        (cj == NULL || (cj != NULL && cj == cj->grav.super)) &&
        (ci == NULL || (ci != NULL && ci == ci->grav.super));

    /* Are we dealing with a task at the top level? */
    if (!(is_ci_top && is_cj_top)) {
      cur->task_in_is_top = 0;
    }
    /* At the hydro level? */
    if (!is_hydro_super) {
      cur->task_in_is_hydro_super = 0;
    }
    /* At the gravity level? */
    if (!is_grav_super) {
      cur->task_in_is_grav_super = 0;
    }

    /* and their dependencies */
    for (int j = 0; j < ta->nr_unlock_tasks; j++) {
      const struct task *tb = ta->unlock_tasks[j];

      /* Are we using this task?
       * For the 0-step, we wish to show all the tasks (even the inactive). */
      if (step != 0 && tb->skip) continue;

      /* Found a task with a dependency. */
      local_count++;

      const struct cell *ci_b = tb->ci;
      const struct cell *cj_b = tb->cj;
      const int is_ci_b_top =
          ci_b == NULL || (ci_b != NULL && ci_b == ci_b->top);
      const int is_cj_b_top =
          cj_b == NULL || (cj_b != NULL && cj_b == cj_b->top);
      const int is_b_hydro_super =
          (cj_b == NULL || (cj_b != NULL && cj_b == cj_b->hydro.super)) &&
          (ci_b == NULL || (ci_b != NULL && ci_b == ci_b->hydro.super));
      const int is_b_grav_super =
          (cj_b == NULL || (cj_b != NULL && cj_b == cj_b->grav.super)) &&
          (ci_b == NULL || (ci_b != NULL && ci_b == ci_b->grav.super));

      int k = 0;
      while (k < MAX_NUMBER_DEP) {
        /* not written yet */
        if (cur->number_link[k] == -1) {
          /* set tb */
          cur->type_out[k] = tb->type;
          cur->subtype_out[k] = tb->subtype;
          cur->implicit_out[k] = tb->implicit;

          /* statistics */
          cur->number_link[k] = 1;
          cur->number_rank[k] = 1;

          /* Are we dealing with a task at the top level? */
          if (!(is_ci_b_top && is_cj_b_top)) {
            cur->task_out_is_top[k] = 0;
          }
          /* At the hydro level? */
          if (!is_b_hydro_super) {
            cur->task_out_is_hydro_super[k] = 0;
          }
          /* At the gravity level? */
          if (!is_b_grav_super) {
            cur->task_out_is_grav_super[k] = 0;
          }

          if (ci_b->cellID == cellID) cell_involved[ind][k]++;
          if (cj_b != NULL && cj_b->cellID == cellID) cell_involved[ind][k]++;

          break;
        }

        /* already written */
        if (cur->type_out[k] == tb->type &&
            cur->subtype_out[k] == tb->subtype) {

          /* Increase the number of link. */
          cur->number_link[k] += 1;

          /* Are we dealing with a task at the top level? */
          if (!(is_ci_b_top && is_cj_b_top)) {
            cur->task_out_is_top[k] = 0;
          }
          /* At the hydro level? */
          if (!is_b_hydro_super) {
            cur->task_out_is_hydro_super[k] = 0;
          }
          /* At the gravity level? */
          if (!is_b_grav_super) {
            cur->task_out_is_grav_super[k] = 0;
          }

          if (ci_b->cellID == cellID) cell_involved[ind][k]++;
          if (cj_b != NULL && cj_b->cellID == cellID) cell_involved[ind][k]++;

          break;
        }

        k += 1;
      }

      /* MAX_NUMBER_DEP is too small */
      if (k == MAX_NUMBER_DEP)
        error("Not enough memory, please increase MAX_NUMBER_DEP");
    }

    /* Some tasks might not unlock anything, like the kick1 tasks. This is
     * expected, and they should turn up in the task graph because they are
     * being unlocked by some other task. However, if a dev missed a
     * dependency and has tasks with no unlocks nor dependencies, they
     * wouldn't show up in the graph. So we write these tasks with no
     * unlocks down too, but as a special case. */
    if (ta->nr_unlock_tasks == 0) {
      cur->number_link[0] = 0;
      cur->type_out[0] = -1;
      cur->subtype_out[0] = -1;
      cur->implicit_out[0] = -1;
      cur->number_rank[0] = 1;
      cell_involved[ind][0] = 1;
    }
  }

  if (local_count > 0) {
    /* We have tasks involving the requested cell on this node */

    /* Create file */
    char filename[50];
    sprintf(filename, "dependency_graph_cell_%lld_step_%i_rank_%i.csv", cellID,
            step, engine_rank);
    FILE *f = fopen(filename, "w");
    if (f == NULL) error("Error opening dependency graph file.");

    /* Write header */
    fprintf(f, "# %s\n", git_revision());
    fprintf(
        f,
        "task_in,task_out,implicit_in,implicit_out,mpi_in,mpi_out,cluster_in,"
        "cluster_out,number_link,number_rank,task_in_is_top,task_in_is_hydro_"
        "super,task_in_is_grav_super,task_out_is_top,task_out_is_hydro_super,"
        "task_out_is_grav_super,cell_has_active_task\n");

    for (int i = 0; i < nber_tasks; i++) {
      for (int j = 0; j < MAX_NUMBER_DEP; j++) {
        /* Does this link exists */
        if (task_dep[i].number_link[j] == -1) {
          continue;
        }

        /* Define a few variables */
        const int ta_type = task_dep[i].type_in;
        const int ta_subtype = task_dep[i].subtype_in;
        const int ta_implicit = task_dep[i].implicit_in;

        const int tb_type = task_dep[i].type_out[j];
        const int tb_subtype = task_dep[i].subtype_out[j];
        const int tb_implicit = task_dep[i].implicit_out[j];

        const int count = task_dep[i].number_link[j];
        const int number_rank = task_dep[i].number_rank[j];

        const int task_in_is_top = task_dep[i].task_in_is_top;
        const int task_in_is_grav_super = task_dep[i].task_in_is_grav_super;
        const int task_in_is_hydro_super = task_dep[i].task_in_is_hydro_super;

        const int task_out_is_top = task_dep[i].task_out_is_top[j];
        const int task_out_is_grav_super =
            task_dep[i].task_out_is_grav_super[j];
        const int task_out_is_hydro_super =
            task_dep[i].task_out_is_hydro_super[j];

        /* text to write */
        char ta_name[200];
        char tb_name[200];

        /* construct line */
        task_get_full_name(ta_type, ta_subtype, ta_name);
        if (tb_type == -1) {
          /* special handling of tasks which have no unlocks */
          strcpy(tb_name, "task_unlocks_nothing");
        } else {
          task_get_full_name(tb_type, tb_subtype, tb_name);
        }

        /* Check if MPI */
        int ta_mpi = 0;
        if (ta_type == task_type_send || ta_type == task_type_recv) ta_mpi = 1;

        int tb_mpi = 0;
        if (tb_type == task_type_send || tb_type == task_type_recv) tb_mpi = 1;

        /* Get group name */
        char ta_cluster[20];
        char tb_cluster[20];
        task_get_group_name(ta_type, ta_subtype, ta_cluster);
        task_get_group_name(tb_type, tb_subtype, tb_cluster);

        fprintf(f, "%s,%s,%d,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
                ta_name, tb_name, ta_implicit, tb_implicit, ta_mpi, tb_mpi,
                ta_cluster, tb_cluster, count, number_rank, task_in_is_top,
                task_in_is_hydro_super, task_in_is_grav_super, task_out_is_top,
                task_out_is_hydro_super, task_out_is_grav_super,
                cell_involved[i][j]);
      }
    }
    /* Close the file */
    fclose(f);
  }

  /* Clean up after yourself */
  free(task_dep);

  if (verbose)
    message("Printing task graph took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#endif /* defined SWIFT_DEBUG_CHECKS || defined CELL_GRAPH */
}

/**
 * @brief Split a hydrodynamic task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_hydro(struct task *t, struct scheduler *s) {
  /* Are we considering both stars and hydro when splitting? */
  /* Note this is not very clean as the scheduler should not really
     access the engine... */
  const int with_feedback = (s->space->e->policy & engine_policy_feedback);
  const int with_stars = (s->space->e->policy & engine_policy_stars);
  const int with_sinks = (s->space->e->policy & engine_policy_sinks);
  const int with_black_holes =
      (s->space->e->policy & engine_policy_black_holes);

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Is this a non-empty self-task? */
    const int is_self =
        (t->type == task_type_self) && (t->ci != NULL) &&
        ((t->ci->hydro.count > 0) || (with_stars && t->ci->stars.count > 0) ||
         (with_sinks && t->ci->sinks.count > 0) ||
         (with_black_holes && t->ci->black_holes.count > 0));

    /* Is this a non-empty pair-task? */
    const int is_pair = (t->type == task_type_pair) && (t->ci != NULL) &&
                        (t->cj != NULL) &&
                        ((t->ci->hydro.count > 0) ||
                         (with_feedback && t->ci->stars.count > 0) ||
                         (with_sinks && t->ci->sinks.count > 0) ||
                         (with_black_holes && t->ci->black_holes.count > 0)) &&
                        ((t->cj->hydro.count > 0) ||
                         (with_feedback && t->cj->stars.count > 0) ||
                         (with_sinks && t->cj->sinks.count > 0) ||
                         (with_black_holes && t->cj->black_holes.count > 0));

    /* Empty task? */
    if (!is_self && !is_pair) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split and the task does not violate h ? */
      if (cell_can_split_self_hydro_task(ci)) {
        /* Make a sub? */
        if (scheduler_dosub && (ci->hydro.count < space_subsize_self_hydro) &&
            (ci->stars.count < space_subsize_self_stars)) {
          /* convert to a self-subtask. */
          t->type = task_type_sub_self;

          /* Otherwise, make tasks explicitly. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[k] != NULL &&
                (ci->progeny[k]->hydro.count ||
                 (with_stars && ci->progeny[k]->stars.count))) {
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);
            }
          }

          /* Make a task for each pair of progeny */
          for (int j = 0; j < 8; j++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[j] != NULL &&
                (ci->progeny[j]->hydro.count ||
                 (with_feedback && ci->progeny[j]->stars.count))) {
              for (int k = j + 1; k < 8; k++) {
                /* Do we have a second non-empty progenitor? */
                if (ci->progeny[k] != NULL &&
                    (ci->progeny[k]->hydro.count ||
                     (with_feedback && ci->progeny[k]->stars.count))) {
                  scheduler_splittask_hydro(
                      scheduler_addtask(s, task_type_pair, t->subtype,
                                        sub_sid_flag[j][k], 0, ci->progeny[j],
                                        ci->progeny[k]),
                      s);
                }
              }
            }
          }
        }

      } /* Cell is split */

    } /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Get the sort ID, use space_getsid and not t->flags
         to make sure we get ci and cj swapped if needed. */
      double shift[3];
      const int sid = space_getsid(s->space, &ci, &cj, shift);

#ifdef SWIFT_DEBUG_CHECKS
      if (sid != t->flags)
        error("Got pair task with incorrect flags: sid=%d flags=%lld", sid,
              t->flags);
#endif

      /* Should this task be split-up? */
      if (cell_can_split_pair_hydro_task(ci) &&
          cell_can_split_pair_hydro_task(cj)) {

        const int h_count_i = ci->hydro.count;
        const int h_count_j = cj->hydro.count;

        const int s_count_i = ci->stars.count;
        const int s_count_j = cj->stars.count;

        int do_sub_hydro = 1;
        int do_sub_stars_i = 1;
        int do_sub_stars_j = 1;
        if (h_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_hydro =
              h_count_i * sid_scale[sid] < space_subsize_pair_hydro / h_count_j;
        }
        if (s_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_stars_i =
              s_count_i * sid_scale[sid] < space_subsize_pair_stars / h_count_j;
        }
        if (s_count_j > 0 && h_count_i > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_stars_j =
              s_count_j * sid_scale[sid] < space_subsize_pair_stars / h_count_i;
        }

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            (do_sub_hydro && do_sub_stars_i && do_sub_stars_j) &&
            !sort_is_corner(sid)) {

          /* Make this task a sub task. */
          t->type = task_type_sub_pair;

          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Loop over the sub-cell pairs for the current sid and add new tasks
           * for them. */
          struct cell_split_pair *csp = &cell_split_pairs[sid];

          t->ci = ci->progeny[csp->pairs[0].pid];
          t->cj = cj->progeny[csp->pairs[0].pjd];
          if (t->ci != NULL) cell_set_flag(t->ci, cell_flag_has_tasks);
          if (t->cj != NULL) cell_set_flag(t->cj, cell_flag_has_tasks);

          t->flags = csp->pairs[0].sid;
          for (int k = 1; k < csp->count; k++) {
            scheduler_splittask_hydro(
                scheduler_addtask(s, task_type_pair, t->subtype,
                                  csp->pairs[k].sid, 0,
                                  ci->progeny[csp->pairs[k].pid],
                                  cj->progeny[csp->pairs[k].pjd]),
                s);
          }
        }

        /* Otherwise, break it up if it is too large? */
      } else if (scheduler_doforcesplit && ci->split && cj->split &&
                 (ci->hydro.count > space_maxsize / cj->hydro.count)) {
        // message( "force splitting pair with %i and %i parts." ,
        // ci->hydro.count , cj->hydro.count );

        /* Replace the current task. */
        t->type = task_type_none;

        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->hydro.count)
            for (int k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL && cj->progeny[k]->hydro.count) {
                struct task *tl =
                    scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                      ci->progeny[j], cj->progeny[k]);
                scheduler_splittask_hydro(tl, s);
                tl->flags = space_getsid(s->space, &t->ci, &t->cj, shift);
              }
      }
    } /* pair interaction? */
  }   /* iterate over the current task. */
}

/**
 * @brief Split a gravity task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_gravity(struct task *t, struct scheduler *s) {
  const struct space *sp = s->space;
  struct engine *e = sp->e;

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_pair && t->cj == NULL)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      const struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Should we split this task? */
      if (cell_can_split_self_gravity_task(ci)) {
        if (scheduler_dosub && ci->grav.count < space_subsize_self_grav) {
          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              scheduler_splittask_gravity(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);

          /* Make a task for each pair of progeny */
          if (t->subtype != task_subtype_external_grav) {
            for (int j = 0; j < 8; j++)
              if (ci->progeny[j] != NULL)
                for (int k = j + 1; k < 8; k++)
                  if (ci->progeny[k] != NULL)
                    scheduler_splittask_gravity(
                        scheduler_addtask(s, task_type_pair, t->subtype,
                                          sub_sid_flag[j][k], 0, ci->progeny[j],
                                          ci->progeny[k]),
                        s);

          } /* Self-gravity only */
        }   /* Make tasks explicitly */
      }     /* Cell is split */
    }       /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Should this task be split-up? */
      if (cell_can_split_pair_gravity_task(ci) &&
          cell_can_split_pair_gravity_task(cj)) {
        const long long gcount_i = ci->grav.count;
        const long long gcount_j = cj->grav.count;

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            gcount_i * gcount_j < ((long long)space_subsize_pair_grav)) {
          /* Otherwise, split it. */
        } else {
          /* Turn the task into a M-M task that will take care of all the
           * progeny pairs */
          t->type = task_type_grav_mm;
          t->subtype = task_subtype_none;
          t->flags = 0;

          /* Make a task for every other pair of progeny */
          for (int i = 0; i < 8; i++) {
            if (ci->progeny[i] != NULL) {
              for (int j = 0; j < 8; j++) {
                if (cj->progeny[j] != NULL) {
                  /* Can we use a M-M interaction here? */
                  if (cell_can_use_pair_mm(ci->progeny[i], cj->progeny[j], e,
                                           sp, /*use_rebuild_data=*/1,
                                           /*is_tree_walk=*/1)) {

                    /* Flag this pair as being treated by the M-M task.
                     * We use the 64 bits in the task->flags field to store
                     * this information. The corresponding taks will unpack
                     * the information and operate according to the choices
                     * made here. */
                    const int flag = i * 8 + j;
                    t->flags |= (1ULL << flag);

                  } else {
                    /* Ok, we actually have to create a task */
                    scheduler_splittask_gravity(
                        scheduler_addtask(s, task_type_pair, task_subtype_grav,
                                          0, 0, ci->progeny[i], cj->progeny[j]),
                        s);
                  }
                }
              }
            }
          }

          /* Can none of the progenies use M-M calculations? */
          if (t->flags == 0) {
            t->type = task_type_none;
            t->subtype = task_subtype_none;
            t->ci = NULL;
            t->cj = NULL;
            t->skip = 1;
          }

        } /* Split the pair */
      }
    } /* pair interaction? */
  }   /* iterate over the current task. */
}

/**
 * @brief Split a FOF task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_fof(struct task *t, struct scheduler *s) {

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {

    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_fof_pair && t->cj == NULL) ||
        t->ci->grav.count == 0 || (t->cj != NULL && t->cj->grav.count == 0)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_fof_self) {

      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split? */
      if (cell_can_split_self_fof_task(ci)) {

        /* Take a step back (we're going to recycle the current task)... */
        redo = 1;

        /* Add the self tasks. */
        int first_child = 0;
        while (ci->progeny[first_child] == NULL) first_child++;
        t->ci = ci->progeny[first_child];
        for (int k = first_child + 1; k < 8; k++)
          if (ci->progeny[k] != NULL && ci->progeny[k]->grav.count)
            scheduler_splittask_fof(
                scheduler_addtask(s, task_type_fof_self, t->subtype, 0, 0,
                                  ci->progeny[k], NULL),
                s);

        /* Make a task for each pair of progeny */
        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->grav.count)
            for (int k = j + 1; k < 8; k++)
              if (ci->progeny[k] != NULL && ci->progeny[k]->grav.count)
                scheduler_splittask_fof(
                    scheduler_addtask(s, task_type_fof_pair, t->subtype, 0, 0,
                                      ci->progeny[j], ci->progeny[k]),
                    s);
      } /* Cell is split */

    } /* Self interaction */

  } /* iterate over the current task. */
}

/**
 * @brief Mapper function to split FOF tasks that may be too large.
 *
 * @param map_data the tasks to process
 * @param num_elements the number of tasks.
 * @param extra_data The #scheduler we are working in.
 */
void scheduler_splittasks_fof_mapper(void *map_data, int num_elements,
                                     void *extra_data) {
  /* Extract the parameters. */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *tasks = (struct task *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Invoke the correct splitting strategy */
    if (t->type == task_type_fof_self || t->type == task_type_fof_pair) {
      scheduler_splittask_fof(t, s);
    }
  }
}

/**
 * @brief Mapper function to split non-FOF tasks that may be too large.
 *
 * @param map_data the tasks to process
 * @param num_elements the number of tasks.
 * @param extra_data The #scheduler we are working in.
 */
void scheduler_splittasks_mapper(void *map_data, int num_elements,
                                 void *extra_data) {
  /* Extract the parameters. */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *tasks = (struct task *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Invoke the correct splitting strategy */
    if (t->subtype == task_subtype_density) {
      scheduler_splittask_hydro(t, s);
    } else if (t->subtype == task_subtype_external_grav) {
      scheduler_splittask_gravity(t, s);
    } else if (t->subtype == task_subtype_grav) {
      scheduler_splittask_gravity(t, s);
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Unexpected task sub-type %s/%s", taskID_names[t->type],
            subtaskID_names[t->subtype]);
#endif
    }
  }
}

/**
 * @brief Splits all the tasks in the scheduler that are too large.
 *
 * @param s The #scheduler.
 * @param fof_tasks Are we splitting the FOF tasks (1)? Or the regular tasks
 * (0)?
 * @param verbose Are we talkative?
 */
void scheduler_splittasks(struct scheduler *s, const int fof_tasks,
                          const int verbose) {

  if (verbose) {
    message("space_subsize_self_hydro= %d", space_subsize_self_hydro);
    message("space_subsize_pair_hydro= %d", space_subsize_pair_hydro);
    message("space_subsize_self_stars= %d", space_subsize_self_stars);
    message("space_subsize_pair_stars= %d", space_subsize_pair_stars);
    message("space_subsize_self_grav= %d", space_subsize_self_grav);
    message("space_subsize_pair_grav= %d", space_subsize_pair_grav);
  }

  if (fof_tasks) {
    /* Call the mapper on each current task. */
    threadpool_map(s->threadpool, scheduler_splittasks_fof_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   s);

  } else {
    /* Call the mapper on each current task. */
    threadpool_map(s->threadpool, scheduler_splittasks_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   s);
  }
}

/**
 * @brief Add a #task to the #scheduler.
 *
 * @param s The #scheduler we are working in.
 * @param type The type of the task.
 * @param subtype The sub-type of the task.
 * @param flags The flags of the task.
 * @param implicit If true, only use this task to unlock dependencies, i.e.
 *        this task is never enqueued.
 * @param ci The first cell to interact.
 * @param cj The second cell to interact.
 */
struct task *scheduler_addtask(struct scheduler *s, enum task_types type,
                               enum task_subtypes subtype, long long flags,
                               int implicit, struct cell *ci, struct cell *cj) {
  /* Get the next free task. */
  const int ind = atomic_inc(&s->tasks_next);

  /* Overflow? */
  if (ind >= s->size)
    error(
        "Task list overflow (%d). Need to increase "
        "Scheduler:tasks_per_cell.",
        ind);

  /* Get a pointer to the new task. */
  struct task *t = &s->tasks[ind];

  /* Copy the data. */
  t->type = type;
  t->subtype = subtype;
  t->flags = flags;
  t->wait = 0;
  t->ci = ci;
  t->cj = cj;
  t->skip = 1; /* Mark tasks as skip by default. */
  t->implicit = implicit;
  t->weight = 0;
  t->rank = 0;
  t->nr_unlock_tasks = 0;
#ifdef SWIFT_DEBUG_TASKS
  t->rid = -1;
#endif
  t->tic = 0;
  t->toc = 0;
  t->total_ticks = 0;

  if (ci != NULL) cell_set_flag(ci, cell_flag_has_tasks);
  if (cj != NULL) cell_set_flag(cj, cell_flag_has_tasks);

  /* Add an index for it. */
  // lock_lock( &s->lock );
  s->tasks_ind[atomic_inc(&s->nr_tasks)] = ind;
  // lock_unlock_blind( &s->lock );

  /* Return a pointer to the new task. */
  return t;
}

/**
 * @brief Set the unlock pointers in each task.
 *
 * @param s The #scheduler.
 */
void scheduler_set_unlocks(struct scheduler *s) {
  /* Store the counts for each task. */
  int *counts;
  if ((counts = (int *)swift_malloc("counts", sizeof(int) * s->nr_tasks)) ==
      NULL)
    error("Failed to allocate temporary counts array.");
  bzero(counts, sizeof(int) * s->nr_tasks);
  for (int k = 0; k < s->nr_unlocks; k++) {
    counts[s->unlock_ind[k]] += 1;

    /* Check that we are not overflowing */
    if (counts[s->unlock_ind[k]] < 0)
      error(
          "Task (type=%s/%s) unlocking more than %lld other tasks!\n"
          "This likely a result of having tasks at vastly different levels"
          "in the tree.\nYou may want to play with the 'Scheduler' "
          "parameters to modify the task splitting strategy and reduce"
          "the difference in task depths.",
          taskID_names[s->tasks[s->unlock_ind[k]].type],
          subtaskID_names[s->tasks[s->unlock_ind[k]].subtype],
          (1LL << (8 * sizeof(int) - 1)) - 1);
  }

  /* Compute the offset for each unlock block. */
  int *offsets;
  if ((offsets = (int *)swift_malloc("offsets",
                                     sizeof(int) * (s->nr_tasks + 1))) == NULL)
    error("Failed to allocate temporary offsets array.");
  offsets[0] = 0;
  for (int k = 0; k < s->nr_tasks; k++) {
    offsets[k + 1] = offsets[k] + counts[k];

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that we are not overflowing */
    if (offsets[k + 1] < 0) error("Task unlock offset array overflowing");
#endif
  }

  /* Create and fill a temporary array with the sorted unlocks. */
  struct task **unlocks;
  if ((unlocks = (struct task **)swift_malloc(
           "unlocks", sizeof(struct task *) * s->size_unlocks)) == NULL)
    error("Failed to allocate temporary unlocks array.");
  for (int k = 0; k < s->nr_unlocks; k++) {
    const int ind = s->unlock_ind[k];
    unlocks[offsets[ind]] = s->unlocks[k];
    offsets[ind] += 1;
  }

  /* Swap the unlocks. */
  swift_free("unlocks", s->unlocks);
  s->unlocks = unlocks;

  /* Re-set the offsets. */
  offsets[0] = 0;
  for (int k = 1; k < s->nr_tasks; k++)
    offsets[k] = offsets[k - 1] + counts[k - 1];

  /* Set the unlocks in the tasks. */
  for (int k = 0; k < s->nr_tasks; k++) {
    struct task *t = &s->tasks[k];
    t->nr_unlock_tasks = counts[k];
    t->unlock_tasks = &s->unlocks[offsets[k]];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that there are no duplicate unlocks. */
  for (int k = 0; k < s->nr_tasks; k++) {
    struct task *t = &s->tasks[k];
    for (int i = 0; i < t->nr_unlock_tasks; i++) {
      for (int j = i + 1; j < t->nr_unlock_tasks; j++) {
        if (t->unlock_tasks[i] == t->unlock_tasks[j])
          error("duplicate unlock! t->type=%s/%s unlocking type=%s/%s",
                taskID_names[t->type], subtaskID_names[t->subtype],
                taskID_names[t->unlock_tasks[i]->type],
                subtaskID_names[t->unlock_tasks[i]->subtype]);
      }
    }
  }
#endif

  /* Clean up. */
  swift_free("counts", counts);
  swift_free("offsets", offsets);
}

/**
 * @brief Sort the tasks in topological order over all queues.
 *
 * @param s The #scheduler.
 */
void scheduler_ranktasks(struct scheduler *s) {
  struct task *tasks = s->tasks;
  int *tid = s->tasks_ind;
  const int nr_tasks = s->nr_tasks;

  /* Run through the tasks and get all the waits right. */
  for (int i = 0; i < nr_tasks; i++) {
    struct task *t = &tasks[i];

    // Increment the waits of the dependances
    for (int k = 0; k < t->nr_unlock_tasks; k++) {
      t->unlock_tasks[k]->wait++;
    }
  }

  /* Load the tids of tasks with no waits. */
  int left = 0;
  for (int k = 0; k < nr_tasks; k++)
    if (tasks[k].wait == 0) {
      tid[left] = k;
      left += 1;
    }

  /* Main loop. */
  for (int j = 0, rank = 0; j < nr_tasks; rank++) {
    /* Did we get anything? */
    if (j == left) error("Unsatisfiable task dependencies detected.");

    /* Unlock the next layer of tasks. */
    const int left_old = left;
    for (; j < left_old; j++) {
      struct task *t = &tasks[tid[j]];
      t->rank = rank;
      /* message( "task %i of type %s has rank %i." , i ,
          (t->type == task_type_self) ? "self" : (t->type == task_type_pair) ?
         "pair" : "sort" , rank ); */
      for (int k = 0; k < t->nr_unlock_tasks; k++) {
        struct task *u = t->unlock_tasks[k];
        if (--u->wait == 0) {
          tid[left] = u - tasks;
          left += 1;
        }
      }
    }

    /* Move back to the old left (like Sanders!). */
    j = left_old;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the tasks were ranked correctly. */
  for (int k = 1; k < s->nr_tasks; k++)
    if (tasks[tid[k - 1]].rank > tasks[tid[k]].rank)
      error("Task ranking failed.");
#endif
}

/**
 * @brief (Re)allocate the task arrays.
 *
 * @param s The #scheduler.
 * @param size The maximum number of tasks in the #scheduler.
 */
void scheduler_reset(struct scheduler *s, int size) {

  /* Do we need to re-allocate? */
  if (size > s->size) {
    /* Free existing task lists if necessary. */
    scheduler_free_tasks(s);

    /* Allocate the new lists. */
    if (swift_memalign("tasks", (void **)&s->tasks, task_align,
                       size * sizeof(struct task)) != 0)
      error("Failed to allocate task array.");

    if ((s->tasks_ind = (int *)swift_malloc("tasks_ind", sizeof(int) * size)) ==
        NULL)
      error("Failed to allocate task lists.");

    if ((s->tid_active =
             (int *)swift_malloc("tid_active", sizeof(int) * size)) == NULL)
      error("Failed to allocate aactive task lists.");
  }

  /* Reset the counters. */
  s->size = size;
  s->nr_tasks = 0;
  s->tasks_next = 0;
  s->waiting = 0;
  s->nr_unlocks = 0;
  s->completed_unlock_writes = 0;
  s->active_count = 0;
  s->total_ticks = 0;

  /* Set the task pointers in the queues. */
  for (int k = 0; k < s->nr_queues; k++) s->queues[k].tasks = s->tasks;
}

/**
 * @brief Compute the task weights
 *
 * @param s The #scheduler.
 * @param verbose Are we talkative?
 */
void scheduler_reweight(struct scheduler *s, int verbose) {
  const int nr_tasks = s->nr_tasks;
  int *tid = s->tasks_ind;
  struct task *tasks = s->tasks;
  const int nodeID = s->nodeID;
  const float wscale = 0.001f;
  const ticks tic = getticks();

  /* Run through the tasks backwards and set their weights. */
  for (int k = nr_tasks - 1; k >= 0; k--) {
    struct task *t = &tasks[tid[k]];
    float cost = 0.f;
    t->weight = 0.f;

    for (int j = 0; j < t->nr_unlock_tasks; j++)
      t->weight += t->unlock_tasks[j]->weight;

    const float count_i = (t->ci != NULL) ? t->ci->hydro.count : 0.f;
    const float count_j = (t->cj != NULL) ? t->cj->hydro.count : 0.f;
    const float gcount_i = (t->ci != NULL) ? t->ci->grav.count : 0.f;
    const float gcount_j = (t->cj != NULL) ? t->cj->grav.count : 0.f;
    const float scount_i = (t->ci != NULL) ? t->ci->stars.count : 0.f;
    const float scount_j = (t->cj != NULL) ? t->cj->stars.count : 0.f;
    const float sink_count_i = (t->ci != NULL) ? t->ci->sinks.count : 0.f;
    const float sink_count_j = (t->cj != NULL) ? t->cj->sinks.count : 0.f;
    const float bcount_i = (t->ci != NULL) ? t->ci->black_holes.count : 0.f;
    const float bcount_j = (t->cj != NULL) ? t->cj->black_holes.count : 0.f;

    switch (t->type) {
      case task_type_sort:
      case task_type_rt_sort:
        cost = wscale * intrinsics_popcount(t->flags) * count_i *
               (sizeof(int) * 8 - (count_i ? intrinsics_clz(count_i) : 0));
        break;

      case task_type_stars_sort:
        cost = wscale * intrinsics_popcount(t->flags) * scount_i *
               (sizeof(int) * 8 - (scount_i ? intrinsics_clz(scount_i) : 0));
        break;

      case task_type_stars_resort:
        cost = wscale * intrinsics_popcount(t->flags) * scount_i *
               (sizeof(int) * 8 - (scount_i ? intrinsics_clz(scount_i) : 0));
        break;

      case task_type_self:
        if (t->subtype == task_subtype_grav) {
          cost = 1.f * (wscale * gcount_i) * gcount_i;
        } else if (t->subtype == task_subtype_external_grav)
          cost = 1.f * wscale * gcount_i;
        else if (t->subtype == task_subtype_stars_density ||
                 t->subtype == task_subtype_stars_prep1 ||
                 t->subtype == task_subtype_stars_prep2 ||
                 t->subtype == task_subtype_stars_feedback)
          cost = 1.f * wscale * scount_i * count_i;
        else if (t->subtype == task_subtype_sink_swallow ||
                 t->subtype == task_subtype_sink_do_gas_swallow)
          cost = 1.f * wscale * count_i * sink_count_i;
        else if (t->subtype == task_subtype_sink_do_sink_swallow)
          cost = 1.f * wscale * sink_count_i * sink_count_i;
        else if (t->subtype == task_subtype_bh_density ||
                 t->subtype == task_subtype_bh_swallow ||
                 t->subtype == task_subtype_bh_feedback)
          cost = 1.f * wscale * bcount_i * count_i;
        else if (t->subtype == task_subtype_do_gas_swallow)
          cost = 1.f * wscale * count_i;
        else if (t->subtype == task_subtype_do_bh_swallow)
          cost = 1.f * wscale * bcount_i;
        else if (t->subtype == task_subtype_density ||
                 t->subtype == task_subtype_gradient ||
                 t->subtype == task_subtype_force ||
                 t->subtype == task_subtype_limiter)
          cost = 1.f * (wscale * count_i) * count_i;
        else if (t->subtype == task_subtype_rt_gradient)
          cost = 1.f * wscale * count_i * count_i;
        else if (t->subtype == task_subtype_rt_transport)
          cost = 1.f * wscale * count_i * count_i;
        else
          error("Untreated sub-type for selfs: %s",
                subtaskID_names[t->subtype]);
        break;

      case task_type_pair:
        if (t->subtype == task_subtype_grav) {
          if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID)
            cost = 3.f * (wscale * gcount_i) * gcount_j;
          else
            cost = 2.f * (wscale * gcount_i) * gcount_j;

        } else if (t->subtype == task_subtype_stars_density ||
                   t->subtype == task_subtype_stars_prep1 ||
                   t->subtype == task_subtype_stars_prep2 ||
                   t->subtype == task_subtype_stars_feedback) {
          if (t->ci->nodeID != nodeID)
            cost = 3.f * wscale * count_i * scount_j * sid_scale[t->flags];
          else if (t->cj->nodeID != nodeID)
            cost = 3.f * wscale * scount_i * count_j * sid_scale[t->flags];
          else
            cost = 2.f * wscale * (scount_i * count_j + scount_j * count_i) *
                   sid_scale[t->flags];

        } else if (t->subtype == task_subtype_sink_swallow ||
                   t->subtype == task_subtype_sink_do_gas_swallow) {
          if (t->ci->nodeID != nodeID)
            cost = 3.f * wscale * count_i * sink_count_j * sid_scale[t->flags];
          else if (t->cj->nodeID != nodeID)
            cost = 3.f * wscale * sink_count_i * count_j * sid_scale[t->flags];
          else
            cost = 2.f * wscale *
                   (sink_count_i * count_j + sink_count_j * count_i) *
                   sid_scale[t->flags];

        } else if (t->subtype == task_subtype_sink_do_sink_swallow) {
          if (t->ci->nodeID != nodeID)
            cost = 3.f * wscale * sink_count_i * sink_count_j *
                   sid_scale[t->flags];
          else if (t->cj->nodeID != nodeID)
            cost = 3.f * wscale * sink_count_i * sink_count_j *
                   sid_scale[t->flags];
          else
            cost = 2.f * wscale *
                   (sink_count_i * sink_count_j + sink_count_j * sink_count_i) *
                   sid_scale[t->flags];

        } else if (t->subtype == task_subtype_bh_density ||
                   t->subtype == task_subtype_bh_swallow ||
                   t->subtype == task_subtype_bh_feedback) {
          if (t->ci->nodeID != nodeID)
            cost = 3.f * wscale * count_i * bcount_j * sid_scale[t->flags];
          else if (t->cj->nodeID != nodeID)
            cost = 3.f * wscale * bcount_i * count_j * sid_scale[t->flags];
          else
            cost = 2.f * wscale * (bcount_i * count_j + bcount_j * count_i) *
                   sid_scale[t->flags];

        } else if (t->subtype == task_subtype_do_gas_swallow) {
          cost = 1.f * wscale * (count_i + count_j);

        } else if (t->subtype == task_subtype_do_bh_swallow) {
          cost = 1.f * wscale * (bcount_i + bcount_j);

        } else if (t->subtype == task_subtype_density ||
                   t->subtype == task_subtype_gradient ||
                   t->subtype == task_subtype_force ||
                   t->subtype == task_subtype_limiter) {
          if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID)
            cost = 3.f * (wscale * count_i) * count_j * sid_scale[t->flags];
          else
            cost = 2.f * (wscale * count_i) * count_j * sid_scale[t->flags];

        } else if (t->subtype == task_subtype_rt_gradient) {
          cost = 1.f * wscale * count_i * count_j;
        } else if (t->subtype == task_subtype_rt_transport) {
          cost = 1.f * wscale * count_i * count_j;
        } else {
          error("Untreated sub-type for pairs: %s",
                subtaskID_names[t->subtype]);
        }
        break;

      case task_type_sub_pair:
#ifdef SWIFT_DEBUG_CHECKS
        if (t->flags < 0) error("Negative flag value!");
#endif
        if (t->subtype == task_subtype_stars_density ||
            t->subtype == task_subtype_stars_prep1 ||
            t->subtype == task_subtype_stars_prep2 ||
            t->subtype == task_subtype_stars_feedback) {
          if (t->ci->nodeID != nodeID) {
            cost = 3.f * (wscale * count_i) * scount_j * sid_scale[t->flags];
          } else if (t->cj->nodeID != nodeID) {
            cost = 3.f * (wscale * scount_i) * count_j * sid_scale[t->flags];
          } else {
            cost = 2.f * wscale * (scount_i * count_j + scount_j * count_i) *
                   sid_scale[t->flags];
          }

        } else if (t->subtype == task_subtype_sink_swallow ||
                   t->subtype == task_subtype_sink_do_gas_swallow) {
          if (t->ci->nodeID != nodeID) {
            cost =
                3.f * (wscale * count_i) * sink_count_j * sid_scale[t->flags];
          } else if (t->cj->nodeID != nodeID) {
            cost =
                3.f * (wscale * sink_count_i) * count_j * sid_scale[t->flags];
          } else {
            cost = 2.f * wscale *
                   (sink_count_i * count_j + sink_count_j * count_i) *
                   sid_scale[t->flags];
          }

        } else if (t->subtype == task_subtype_sink_do_sink_swallow) {
          if (t->ci->nodeID != nodeID) {
            cost = 3.f * (wscale * sink_count_i) * sink_count_j *
                   sid_scale[t->flags];
          } else if (t->cj->nodeID != nodeID) {
            cost = 3.f * (wscale * sink_count_i) * sink_count_j *
                   sid_scale[t->flags];
          } else {
            cost = 2.f * wscale *
                   (sink_count_i * sink_count_j + sink_count_j * sink_count_i) *
                   sid_scale[t->flags];
          }
        } else if (t->subtype == task_subtype_bh_density ||
                   t->subtype == task_subtype_bh_swallow ||
                   t->subtype == task_subtype_bh_feedback) {
          if (t->ci->nodeID != nodeID) {
            cost = 3.f * (wscale * count_i) * bcount_j * sid_scale[t->flags];
          } else if (t->cj->nodeID != nodeID) {
            cost = 3.f * (wscale * bcount_i) * count_j * sid_scale[t->flags];
          } else {
            cost = 2.f * wscale * (bcount_i * count_j + bcount_j * count_i) *
                   sid_scale[t->flags];
          }

        } else if (t->subtype == task_subtype_do_gas_swallow) {
          cost = 1.f * wscale * (count_i + count_j);

        } else if (t->subtype == task_subtype_do_bh_swallow) {
          cost = 1.f * wscale * (bcount_i + bcount_j);

        } else if (t->subtype == task_subtype_density ||
                   t->subtype == task_subtype_gradient ||
                   t->subtype == task_subtype_force ||
                   t->subtype == task_subtype_limiter) {
          if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID) {
            cost = 3.f * (wscale * count_i) * count_j * sid_scale[t->flags];
          } else {
            cost = 2.f * (wscale * count_i) * count_j * sid_scale[t->flags];
          }
        } else if (t->subtype == task_subtype_rt_gradient) {
          cost = 1.f * wscale * count_i * count_j;
        } else if (t->subtype == task_subtype_rt_transport) {
          cost = 1.f * wscale * count_i * count_j;
        } else {
          error("Untreated sub-type for sub-pairs: %s",
                subtaskID_names[t->subtype]);
        }
        break;

      case task_type_sub_self:
        if (t->subtype == task_subtype_stars_density ||
            t->subtype == task_subtype_stars_prep1 ||
            t->subtype == task_subtype_stars_prep2 ||
            t->subtype == task_subtype_stars_feedback) {
          cost = 1.f * (wscale * scount_i) * count_i;
        } else if (t->subtype == task_subtype_sink_swallow ||
                   t->subtype == task_subtype_sink_do_gas_swallow) {
          cost = 1.f * (wscale * sink_count_i) * count_i;
        } else if (t->subtype == task_subtype_sink_do_sink_swallow) {
          cost = 1.f * (wscale * sink_count_i) * sink_count_i;
        } else if (t->subtype == task_subtype_bh_density ||
                   t->subtype == task_subtype_bh_swallow ||
                   t->subtype == task_subtype_bh_feedback) {
          cost = 1.f * (wscale * bcount_i) * count_i;
        } else if (t->subtype == task_subtype_do_gas_swallow) {
          cost = 1.f * wscale * count_i;
        } else if (t->subtype == task_subtype_do_bh_swallow) {
          cost = 1.f * wscale * bcount_i;
        } else if (t->subtype == task_subtype_density ||
                   t->subtype == task_subtype_gradient ||
                   t->subtype == task_subtype_force ||
                   t->subtype == task_subtype_limiter) {
          cost = 1.f * (wscale * count_i) * count_i;
        } else if (t->subtype == task_subtype_rt_gradient) {
          cost = 1.f * wscale * scount_i * count_i;
        } else if (t->subtype == task_subtype_rt_transport) {
          cost = 1.f * wscale * scount_i * count_i;
        } else {
          error("Untreated sub-type for sub-selfs: %s",
                subtaskID_names[t->subtype]);
        }
        break;
      case task_type_ghost:
        if (t->ci == t->ci->hydro.super) cost = wscale * count_i;
        break;
      case task_type_extra_ghost:
        if (t->ci == t->ci->hydro.super) cost = wscale * count_i;
        break;
      case task_type_stars_ghost:
        if (t->ci == t->ci->hydro.super) cost = wscale * scount_i;
        break;
      case task_type_bh_density_ghost:
        if (t->ci == t->ci->hydro.super) cost = wscale * bcount_i;
        break;
      case task_type_bh_swallow_ghost2:
        if (t->ci == t->ci->hydro.super) cost = wscale * bcount_i;
        break;
      case task_type_drift_part:
        cost = wscale * count_i;
        break;
      case task_type_drift_gpart:
        cost = wscale * gcount_i;
        break;
      case task_type_drift_spart:
        cost = wscale * scount_i;
        break;
      case task_type_drift_sink:
        cost = wscale * sink_count_i;
        break;
      case task_type_drift_bpart:
        cost = wscale * bcount_i;
        break;
      case task_type_init_grav:
        cost = wscale * gcount_i;
        break;
      case task_type_grav_down:
        cost = wscale * gcount_i;
        break;
      case task_type_grav_long_range:
        cost = wscale * gcount_i;
        break;
      case task_type_grav_mm:
        cost = wscale * (gcount_i + gcount_j);
        break;
      case task_type_end_hydro_force:
        cost = wscale * count_i;
        break;
      case task_type_end_grav_force:
        cost = wscale * gcount_i;
        break;
      case task_type_cooling:
        cost = wscale * count_i;
        break;
      case task_type_star_formation:
        cost = wscale * (count_i + scount_i);
        break;
      case task_type_star_formation_sink:
        cost = wscale * (sink_count_i + scount_i);
        break;
      case task_type_sink_formation:
        cost = wscale * (count_i + sink_count_i);
        break;
      case task_type_rt_ghost1:
        cost = wscale * count_i;
        break;
      case task_type_rt_ghost2:
        cost = wscale * count_i;
        break;
      case task_type_rt_tchem:
        cost = wscale * count_i;
        break;
      case task_type_rt_advance_cell_time:
      case task_type_rt_collect_times:
        cost = wscale;
        break;
      case task_type_csds:
        cost =
            wscale * (count_i + gcount_i + scount_i + sink_count_i + bcount_i);
        break;
      case task_type_kick1:
        cost =
            wscale * (count_i + gcount_i + scount_i + sink_count_i + bcount_i);
        break;
      case task_type_kick2:
        cost =
            wscale * (count_i + gcount_i + scount_i + sink_count_i + bcount_i);
        break;
      case task_type_timestep:
        cost =
            wscale * (count_i + gcount_i + scount_i + sink_count_i + bcount_i);
        break;
      case task_type_timestep_limiter:
        cost = wscale * count_i;
        break;
      case task_type_timestep_sync:
        cost = wscale * count_i;
        break;
      case task_type_send:
        if (count_i < 1e5)
          cost = 10.f * (wscale * count_i) * count_i;
        else
          cost = 2e9;
        break;
      case task_type_recv:
        if (count_i < 1e5)
          cost = 5.f * (wscale * count_i) * count_i;
        else
          cost = 1e9;
        break;
      default:
        cost = 0;
        break;
    }
    t->weight += cost;
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* int min = tasks[0].weight, max = tasks[0].weight;
  for ( int k = 1 ; k < nr_tasks ; k++ )
      if ( tasks[k].weight < min )
          min = tasks[k].weight;
      else if ( tasks[k].weight > max )
          max = tasks[k].weight;
  message( "task weights are in [ %i , %i ]." , min , max ); */
}

/**
 * @brief #threadpool_map function which runs through the task
 *        graph and re-computes the task wait counters.
 */
void scheduler_rewait_mapper(void *map_data, int num_elements,
                             void *extra_data) {
  struct scheduler *s = (struct scheduler *)extra_data;
  const int *tid = (int *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &s->tasks[tid[ind]];

    /* Ignore skipped tasks. */
    if (t->skip) continue;

    /* Increment the task's own wait counter for the enqueueing. */
    atomic_inc(&t->wait);

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that we don't have more waits that what can be stored. */
    if (t->wait < 0)
      error("Task (type=%s/%s) unlocked by more than %lld tasks!",
            taskID_names[t->type], subtaskID_names[t->subtype],
            (1LL << (8 * sizeof(t->wait) - 1)) - 1);
#endif

    /* Sets the waits of the dependances */
    for (int k = 0; k < t->nr_unlock_tasks; k++) {
      struct task *u = t->unlock_tasks[k];
      atomic_inc(&u->wait);
    }
  }
}

void scheduler_enqueue_mapper(void *map_data, int num_elements,
                              void *extra_data) {
  struct scheduler *s = (struct scheduler *)extra_data;
  const int *tid = (int *)map_data;
  struct task *tasks = s->tasks;
  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[tid[ind]];
    if (atomic_dec(&t->wait) == 1 && !t->skip) {
      scheduler_enqueue(s, t);
    }
  }
  pthread_cond_broadcast(&s->sleep_cond);
}

/**
 * @brief Start the scheduler, i.e. fill the queues with ready tasks.
 *
 * @param s The #scheduler.
 */
void scheduler_start(struct scheduler *s) {

  /* Re-wait the tasks. */
  if (s->active_count > 1000) {
    threadpool_map(s->threadpool, scheduler_rewait_mapper, s->tid_active,
                   s->active_count, sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    scheduler_rewait_mapper(s->tid_active, s->active_count, s);
  }

  /* Loop over the tasks and enqueue whoever is ready. */
  if (s->active_count > 1000) {
    threadpool_map(s->threadpool, scheduler_enqueue_mapper, s->tid_active,
                   s->active_count, sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    scheduler_enqueue_mapper(s->tid_active, s->active_count, s);
  }

  /* Clear the list of active tasks. */
  s->active_count = 0;

  /* To be safe, fire of one last sleep_cond in a safe way. */
  pthread_mutex_lock(&s->sleep_mutex);
  pthread_cond_broadcast(&s->sleep_cond);
  pthread_mutex_unlock(&s->sleep_mutex);
}

/**
 * @brief Put a task on one of the queues.
 *
 * @param s The #scheduler.
 * @param t The #task.
 */
void scheduler_enqueue(struct scheduler *s, struct task *t) {

  /* Ignore skipped tasks */
  if (t->skip) return;

  /* If this is an implicit task, just pretend it's done. */
  if (t->implicit) {
#ifdef SWIFT_DEBUG_CHECKS
    t->ti_run = s->space->e->ti_current;

    /* Mark that we have run this task on these cells */
    if (t->ci != NULL) {
      t->ci->tasks_executed[t->type]++;
      t->ci->subtasks_executed[t->subtype]++;
    }
    if (t->cj != NULL) {
      t->cj->tasks_executed[t->type]++;
      t->cj->subtasks_executed[t->subtype]++;
    }
#endif
    t->skip = 1;
    for (int j = 0; j < t->nr_unlock_tasks; j++) {
      struct task *t2 = t->unlock_tasks[j];
      if (atomic_dec(&t2->wait) == 1) scheduler_enqueue(s, t2);
    }
  }

  /* Otherwise, look for a suitable queue. */
  else {
#ifdef WITH_MPI
    int err = MPI_SUCCESS;
#endif

    /* Find the previous owner for each task type, and do
     * any pre-processing needed. */
    short int qid = -1;
    short int *owner = NULL;
    switch (t->type) {
      case task_type_self:
      case task_type_sub_self:
        if (t->subtype == task_subtype_grav ||
            t->subtype == task_subtype_external_grav) {
          qid = t->ci->grav.super->owner;
          owner = &t->ci->grav.super->owner;
        } else {
          qid = t->ci->hydro.super->owner;
          owner = &t->ci->hydro.super->owner;
        }
        break;
      case task_type_sort:
      case task_type_ghost:
      case task_type_drift_part:
        qid = t->ci->hydro.super->owner;
        owner = &t->ci->hydro.super->owner;
        break;
      case task_type_drift_gpart:
        qid = t->ci->grav.super->owner;
        owner = &t->ci->grav.super->owner;
        break;
      case task_type_kick1:
      case task_type_kick2:
      case task_type_stars_ghost:
      case task_type_csds:
      case task_type_stars_sort:
      case task_type_timestep:
        qid = t->ci->super->owner;
        owner = &t->ci->super->owner;
        break;
      case task_type_pair:
      case task_type_sub_pair:
        qid = t->ci->super->owner;
        owner = &t->ci->super->owner;
        if ((qid < 0) ||
            ((t->cj->super->owner > -1) &&
             (s->queues[qid].count > s->queues[t->cj->super->owner].count))) {
          qid = t->cj->super->owner;
          owner = &t->cj->super->owner;
        }
        break;
      case task_type_recv:
#ifdef WITH_MPI
      {
        size_t size = 0;              /* Size in bytes. */
        size_t count = 0;             /* Number of elements to receive */
        MPI_Datatype type = MPI_BYTE; /* Type of the elements */
        void *buff = NULL;            /* Buffer to accept elements */

        if (t->subtype == task_subtype_tend) {

          count = size = t->ci->mpi.pcell_size * sizeof(struct pcell_step);
          buff = t->buff = malloc(count);

        } else if (t->subtype == task_subtype_part_swallow) {

          count = size =
              t->ci->hydro.count * sizeof(struct black_holes_part_data);
          buff = t->buff = malloc(count);

        } else if (t->subtype == task_subtype_bpart_merger) {
          count = size =
              sizeof(struct black_holes_bpart_data) * t->ci->black_holes.count;
          buff = t->buff = malloc(count);

        } else if (t->subtype == task_subtype_xv ||
                   t->subtype == task_subtype_rho ||
                   t->subtype == task_subtype_gradient ||
                   t->subtype == task_subtype_rt_gradient ||
                   t->subtype == task_subtype_rt_transport ||
                   t->subtype == task_subtype_part_prep1) {

          count = t->ci->hydro.count;
          size = count * sizeof(struct part);
          type = part_mpi_type;
          buff = t->ci->hydro.parts;

        } else if (t->subtype == task_subtype_limiter) {

          size = count = t->ci->hydro.count * sizeof(timebin_t);
          if (posix_memalign((void **)&buff, SWIFT_CACHE_ALIGNMENT, count) != 0)
            error("Error allocating timebin recv buffer");
          type = MPI_BYTE;
          t->buff = buff;
          task_get_unique_dependent(t)->buff = buff;

        } else if (t->subtype == task_subtype_gpart) {

          count = t->ci->grav.count;
          size = count * sizeof(struct gpart);
          type = gpart_mpi_type;
          buff = t->ci->grav.parts;

        } else if (t->subtype == task_subtype_spart_density ||
                   t->subtype == task_subtype_spart_prep2) {

          count = t->ci->stars.count;
          size = count * sizeof(struct spart);
          type = spart_mpi_type;
          buff = t->ci->stars.parts;

        } else if (t->subtype == task_subtype_bpart_rho ||
                   t->subtype == task_subtype_bpart_feedback) {

          count = t->ci->black_holes.count;
          size = count * sizeof(struct bpart);
          type = bpart_mpi_type;
          buff = t->ci->black_holes.parts;

        } else if (t->subtype == task_subtype_sf_counts) {

          count = size = t->ci->mpi.pcell_size * sizeof(struct pcell_sf);
          buff = t->buff = malloc(count);

        } else {
          error("Unknown communication sub-type");
        }

        err = MPI_Irecv(buff, count, type, t->ci->nodeID, t->flags,
                        subtaskMPI_comms[t->subtype], &t->req);

        if (err != MPI_SUCCESS) {
          mpi_error(err, "Failed to emit irecv for particle data.");
        }

        /* And log, if logging enabled. */
        mpiuse_log_allocation(t->type, t->subtype, &t->req, 1, size,
                              t->ci->nodeID, t->flags);

        qid = 1 % s->nr_queues;
      }
#else
        error("SWIFT was not compiled with MPI support.");
#endif
      break;
      case task_type_send:
#ifdef WITH_MPI
      {
        size_t size = 0;              /* Size in bytes. */
        size_t count = 0;             /* Number of elements to send */
        MPI_Datatype type = MPI_BYTE; /* Type of the elements */
        void *buff = NULL;            /* Buffer to send */

        if (t->subtype == task_subtype_tend) {

          size = count = t->ci->mpi.pcell_size * sizeof(struct pcell_step);
          buff = t->buff = malloc(size);
          cell_pack_end_step(t->ci, (struct pcell_step *)buff);

        } else if (t->subtype == task_subtype_part_swallow) {

          size = count =
              t->ci->hydro.count * sizeof(struct black_holes_part_data);
          buff = t->buff = malloc(size);
          cell_pack_part_swallow(t->ci, (struct black_holes_part_data *)buff);

        } else if (t->subtype == task_subtype_bpart_merger) {

          size = count =
              sizeof(struct black_holes_bpart_data) * t->ci->black_holes.count;
          buff = t->buff = malloc(size);
          cell_pack_bpart_swallow(t->ci,
                                  (struct black_holes_bpart_data *)t->buff);

        } else if (t->subtype == task_subtype_xv ||
                   t->subtype == task_subtype_rho ||
                   t->subtype == task_subtype_gradient ||
                   t->subtype == task_subtype_rt_gradient ||
                   t->subtype == task_subtype_rt_transport ||
                   t->subtype == task_subtype_part_prep1) {

          count = t->ci->hydro.count;
          size = count * sizeof(struct part);
          type = part_mpi_type;
          buff = t->ci->hydro.parts;

        } else if (t->subtype == task_subtype_limiter) {

          size = count = t->ci->hydro.count * sizeof(timebin_t);
          type = MPI_BYTE;
          buff = t->buff;

        } else if (t->subtype == task_subtype_gpart) {

          count = t->ci->grav.count;
          size = count * sizeof(struct gpart);
          type = gpart_mpi_type;
          buff = t->ci->grav.parts;

        } else if (t->subtype == task_subtype_spart_density ||
                   t->subtype == task_subtype_spart_prep2) {

          count = t->ci->stars.count;
          size = count * sizeof(struct spart);
          type = spart_mpi_type;
          buff = t->ci->stars.parts;

        } else if (t->subtype == task_subtype_bpart_rho ||
                   t->subtype == task_subtype_bpart_feedback) {

          count = t->ci->black_holes.count;
          size = count * sizeof(struct bpart);
          type = bpart_mpi_type;
          buff = t->ci->black_holes.parts;

        } else if (t->subtype == task_subtype_sf_counts) {

          size = count = t->ci->mpi.pcell_size * sizeof(struct pcell_sf);
          buff = t->buff = malloc(size);
          cell_pack_sf_counts(t->ci, (struct pcell_sf *)t->buff);

        } else {
          error("Unknown communication sub-type");
        }

        if (size > s->mpi_message_limit) {
          err = MPI_Isend(buff, count, type, t->cj->nodeID, t->flags,
                          subtaskMPI_comms[t->subtype], &t->req);
        } else {
          err = MPI_Issend(buff, count, type, t->cj->nodeID, t->flags,
                           subtaskMPI_comms[t->subtype], &t->req);
        }

        if (err != MPI_SUCCESS) {
          mpi_error(err, "Failed to emit isend for particle data.");
        }

        /* And log, if logging enabled. */
        mpiuse_log_allocation(t->type, t->subtype, &t->req, 1, size,
                              t->cj->nodeID, t->flags);

        qid = 0;
      }
#else
        error("SWIFT was not compiled with MPI support.");
#endif
      break;
      default:
        qid = -1;
    }

    if (qid >= s->nr_queues) error("Bad computed qid.");

    /* If no qid, pick a random queue. */
    if (qid < 0) qid = rand() % s->nr_queues;

    /* Save qid as owner for next time a task accesses this cell. */
    if (owner != NULL) *owner = qid;

    /* Increase the waiting counter. */
    atomic_inc(&s->waiting);

    /* Insert the task into that queue. */
    queue_insert(&s->queues[qid], t);
  }
}

/**
 * @brief Take care of a tasks dependencies.
 *
 * @param s The #scheduler.
 * @param t The finished #task.
 *
 * @return A pointer to the next task, if a suitable one has
 *         been identified.
 */
struct task *scheduler_done(struct scheduler *s, struct task *t) {
  /* Release whatever locks this task held. */
  if (!t->implicit) task_unlock(t);

  /* Loop through the dependencies and add them to a queue if
     they are ready. */
  for (int k = 0; k < t->nr_unlock_tasks; k++) {
    struct task *t2 = t->unlock_tasks[k];
    if (t2->skip) continue;

    const int res = atomic_dec(&t2->wait);
    if (res < 1) {
      error("Negative wait!");
    } else if (res == 1) {
      scheduler_enqueue(s, t2);
    }
  }

  /* Task definitely done, signal any sleeping runners. */
  if (!t->implicit) {
    t->toc = getticks();
    t->total_ticks += t->toc - t->tic;
    pthread_mutex_lock(&s->sleep_mutex);
    atomic_dec(&s->waiting);
    pthread_cond_broadcast(&s->sleep_cond);
    pthread_mutex_unlock(&s->sleep_mutex);
  }

  /* Mark the task as skip. */
  t->skip = 1;

  /* Return the next best task. Note that we currently do not
     implement anything that does this, as getting it to respect
     priorities is too tricky and currently unnecessary. */
  return NULL;
}

/**
 * @brief Resolve a single dependency by hand.
 *
 * @param s The #scheduler.
 * @param t The dependent #task.
 *
 * @return A pointer to the next task, if a suitable one has
 *         been identified.
 */
struct task *scheduler_unlock(struct scheduler *s, struct task *t) {
  /* Loop through the dependencies and add them to a queue if
     they are ready. */
  for (int k = 0; k < t->nr_unlock_tasks; k++) {
    struct task *t2 = t->unlock_tasks[k];
    const int res = atomic_dec(&t2->wait);
    if (res < 1) {
      error("Negative wait!");
    } else if (res == 1) {
      scheduler_enqueue(s, t2);
    }
  }

  /* Task definitely done. */
  if (!t->implicit) {
    t->toc = getticks();
    t->total_ticks += t->toc - t->tic;
    pthread_mutex_lock(&s->sleep_mutex);
    atomic_dec(&s->waiting);
    pthread_cond_broadcast(&s->sleep_cond);
    pthread_mutex_unlock(&s->sleep_mutex);
  }

  /* Return the next best task. Note that we currently do not
     implement anything that does this, as getting it to respect
     priorities is too tricky and currently unnecessary. */
  return NULL;
}

/**
 * Take note of the time at which a task was successfully fetched from the
 * queue.
 *
 * @param s The #scheduler.
 */
void scheduler_mark_last_fetch(struct scheduler *s) {

#if defined(SWIFT_DEBUG_CHECKS)
  if (s->deadlock_waiting_time_ms <= 0.f) return;

  ticks now = getticks();
  ticks last = s->last_successful_task_fetch;
  while (atomic_cas(&s->last_successful_task_fetch, last, now) != last) {
    now = getticks();
    last = s->last_successful_task_fetch;
  }
#endif
}

/**
 * Abort the run if you're stuck doing nothing for too long.
 * This function is intended to abort the mission if you're
 * deadlocked somewhere and somehow. You might get core dumps
 * this way. Alternatively, you might manually set a breakpoint
 * with gdb when this function is called.
 *
 * @param s The #scheduler.
 */
void scheduler_check_deadlock(struct scheduler *s) {

#if defined(SWIFT_DEBUG_CHECKS)
  if (s->deadlock_waiting_time_ms <= 0.f) return;

  /* lock_lock(&s->last_task_fetch_lock); */
  ticks now = getticks();
  ticks last = s->last_successful_task_fetch;

  if (last == 0LL) {
    /* Ensure that the first check each engine_launch doesn't fail. There is no
     * guarantee how long it will take from the point where
     * last_successful_task_fetch was reset to get to this point. A poorly
     * chosen scheduler->deadlock_waiting_time_ms may abort a big run in places
     * where there is no deadlock. Better safe than sorry, so at start-up, the
     * last successful task fetch time is marked as 0. So we just exit without
     * checking the time. */
    while (atomic_cas(&s->last_successful_task_fetch, last, now) != last) {
      now = getticks();
      last = s->last_successful_task_fetch;
    }
    return;
  }

  /* ticks on different CPUs may disagree a bit. So we may end up
   * with last > now, and consequently negative idle time, which
   * then overflows unsigned long longs and gives false positives. */
  const ticks big = max(now, last);
  const ticks small = min(now, last);
  const double idle_time = clocks_diff_ticks(big, small);

  if (idle_time > s->deadlock_waiting_time_ms) {
    message(
        "Detected what looks like a deadlock after %g ms of no new task being "
        "fetched from queues. Dumping diagnostic data.",
        idle_time);
    engine_dump_diagnostic_data(s->e);
    error("Aborting now.");
  }
#endif
}

/**
 * @brief Get a task, preferably from the given queue.
 *
 * @param s The #scheduler.
 * @param qid The ID of the preferred #queue.
 * @param prev the previous task that was run.
 *
 * @return A pointer to a #task or @c NULL if there are no available tasks.
 */
struct task *scheduler_gettask(struct scheduler *s, int qid,
                               const struct task *prev) {
  struct task *res = NULL;
  const int nr_queues = s->nr_queues;
  unsigned int seed = qid;

  /* Check qid. */
  if (qid >= nr_queues || qid < 0) error("Bad queue ID.");

  /* Loop as long as there are tasks... */
  while (s->waiting > 0 && res == NULL) {
    /* Try more than once before sleeping. */
    for (int tries = 0; res == NULL && s->waiting && tries < scheduler_maxtries;
         tries++) {
      /* Try to get a task from the suggested queue. */
      if (s->queues[qid].count > 0 || s->queues[qid].count_incoming > 0) {
        TIMER_TIC
        res = queue_gettask(&s->queues[qid], prev, 0);
        TIMER_TOC(timer_qget);
        if (res != NULL) break;
      }

      /* If unsuccessful, try stealing from the other queues. */
      if (s->flags & scheduler_flag_steal) {
        int count = 0, qids[nr_queues];
        for (int k = 0; k < nr_queues; k++)
          if (s->queues[k].count > 0 || s->queues[k].count_incoming > 0) {
            qids[count++] = k;
          }
        for (int k = 0; k < scheduler_maxsteal && count > 0; k++) {
          const int ind = rand_r(&seed) % count;
          TIMER_TIC
          res = queue_gettask(&s->queues[qids[ind]], prev, 0);
          TIMER_TOC(timer_qsteal);
          if (res != NULL) {
            break;
          } else {
            qids[ind] = qids[--count];
          }
        }
        if (res != NULL) break;
      }
    }

/* If we failed, take a short nap. */
#ifdef WITH_MPI
    if (res == NULL && qid > 1)
#else
    if (res == NULL)
#endif
    {
      pthread_mutex_lock(&s->sleep_mutex);
      res = queue_gettask(&s->queues[qid], prev, 1);
      if (res == NULL && s->waiting > 0) {
        pthread_cond_wait(&s->sleep_cond, &s->sleep_mutex);
      }
      pthread_mutex_unlock(&s->sleep_mutex);
    }

    scheduler_check_deadlock(s);
  }

  if (res != NULL) {
    scheduler_mark_last_fetch(s);
    /* Start the timer on this task, if we got one. */
    res->tic = getticks();
#ifdef SWIFT_DEBUG_TASKS
    res->rid = qid;
#endif
  }

  /* No milk today. */
  return res;
}

/**
 * @brief Initialize the #scheduler.
 *
 * @param s The #scheduler.
 * @param space The #space we are working with
 * @param nr_tasks The number of tasks to allocate initially.
 * @param nr_queues The number of queues in this scheduler.
 * @param flags The #scheduler flags.
 * @param nodeID The MPI rank
 * @param tp Parallel processing threadpool.
 */
void scheduler_init(struct scheduler *s, struct space *space, int nr_tasks,
                    int nr_queues, unsigned int flags, int nodeID,
                    struct threadpool *tp) {
  /* Init the lock. */
  lock_init(&s->lock);

  /* Allocate the queues. */
  if (swift_memalign("queues", (void **)&s->queues, queue_struct_align,
                     sizeof(struct queue) * nr_queues) != 0)
    error("Failed to allocate queues.");

  /* Initialize each queue. */
  for (int k = 0; k < nr_queues; k++) queue_init(&s->queues[k], NULL);

  /* Init the sleep mutex and cond. */
  if (pthread_cond_init(&s->sleep_cond, NULL) != 0 ||
      pthread_mutex_init(&s->sleep_mutex, NULL) != 0)
    error("Failed to initialize sleep barrier.");

  /* Init the unlocks. */
  if ((s->unlocks = (struct task **)swift_malloc(
           "unlocks", sizeof(struct task *) * scheduler_init_nr_unlocks)) ==
          NULL ||
      (s->unlock_ind = (int *)swift_malloc(
           "unlock_ind", sizeof(int) * scheduler_init_nr_unlocks)) == NULL)
    error("Failed to allocate unlocks.");
  s->nr_unlocks = 0;
  s->size_unlocks = scheduler_init_nr_unlocks;

  /* Set the scheduler variables. */
  s->nr_queues = nr_queues;
  s->flags = flags;
  s->space = space;
  s->nodeID = nodeID;
  s->threadpool = tp;

  /* Init the tasks array. */
  s->size = 0;
  s->tasks = NULL;
  s->tasks_ind = NULL;
  scheduler_reset(s, nr_tasks);

#if defined(SWIFT_DEBUG_CHECKS)
  s->e = space->e;
  s->last_successful_task_fetch = 0LL;
#endif
}

/**
 * @brief Prints the list of tasks to a file
 *
 * @param s The #scheduler
 * @param fileName Name of the file to write to
 */
void scheduler_print_tasks(const struct scheduler *s, const char *fileName) {
  const int nr_tasks = s->nr_tasks, *tid = s->tasks_ind;
  struct task *t, *tasks = s->tasks;

  FILE *file = fopen(fileName, "w");
  if (file == NULL) error("Could not create file '%s'.", fileName);

  fprintf(file, "# Rank  Name  Subname  unlocks  waits\n");

  for (int k = nr_tasks - 1; k >= 0; k--) {
    t = &tasks[tid[k]];
    if (t->skip) continue;
    fprintf(file, "%d %s %s %d %d\n", k, taskID_names[t->type],
            subtaskID_names[t->subtype], t->nr_unlock_tasks, t->wait);
  }

  fclose(file);
}

/**
 * @brief Frees up the memory allocated for this #scheduler
 */
void scheduler_clean(struct scheduler *s) {
  scheduler_free_tasks(s);
  swift_free("unlocks", s->unlocks);
  swift_free("unlock_ind", s->unlock_ind);
  for (int i = 0; i < s->nr_queues; ++i) queue_clean(&s->queues[i]);
  swift_free("queues", s->queues);
}

/**
 * @brief Free the task arrays allocated by this #scheduler.
 */
void scheduler_free_tasks(struct scheduler *s) {
  if (s->tasks != NULL) {
    swift_free("tasks", s->tasks);
    s->tasks = NULL;
  }
  if (s->tasks_ind != NULL) {
    swift_free("tasks_ind", s->tasks_ind);
    s->tasks_ind = NULL;
  }
  if (s->tid_active != NULL) {
    swift_free("tid_active", s->tid_active);
    s->tid_active = NULL;
  }
  s->size = 0;
  s->nr_tasks = 0;
}

/**
 * @brief write down the levels and the number of tasks at that level.
 *
 * Run plot_task_level.py for an example of how to use it
 * to generate the figure.
 *
 * @param s The #scheduler we are working in.
 * @param step The current step number.
 */
void scheduler_write_task_level(const struct scheduler *s, int step) {

  /* init */
  const int max_depth = 30;
  const struct task *tasks = s->tasks;
  int nr_tasks = s->nr_tasks;

  /* Init counter */
  int size = task_type_count * task_subtype_count * max_depth;
  int *count = (int *)malloc(size * sizeof(int));
  if (count == NULL) error("Failed to allocate memory");

  for (int i = 0; i < size; i++) count[i] = 0;

  /* Count tasks */
  for (int i = 0; i < nr_tasks; i++) {
    const struct task *t = &tasks[i];
    if (t->ci) {
      if ((int)t->ci->depth >= max_depth)
        error("Cell is too deep, you need to increase max_depth");

      int ind = t->type * task_subtype_count * max_depth;
      ind += t->subtype * max_depth;
      ind += (int)t->ci->depth;

      count[ind] += 1;
    }
  }

  /* Generate filename */
  char filename[200] = "task_level_\0";
#ifdef WITH_MPI
  char rankstr[6];
  sprintf(rankstr, "%04d_", s->nodeID);
  strcat(filename, rankstr);
#endif
  char stepstr[100];
  sprintf(stepstr, "%d.txt", step);
  strcat(filename, stepstr);

  /* Open file */
  FILE *f = fopen(filename, "w");
  if (f == NULL) error("Error opening task level file.");

  /* Print header */
  fprintf(f, "# task_type, task_subtype, depth, count\n");

  /* Print tasks level */
  for (int i = 0; i < size; i++) {
    if (count[i] == 0) continue;

    int type = i / (task_subtype_count * max_depth);
    int subtype = i - task_subtype_count * max_depth * type;
    subtype /= max_depth;
    int depth = i - task_subtype_count * max_depth * type;
    depth -= subtype * max_depth;
    fprintf(f, "%s %s %i %i\n", taskID_names[type], subtaskID_names[subtype],
            depth, count[i]);
  }

  /* clean up */
  fclose(f);
  free(count);
}
/**
 * @brief dump all the active queues of all the known schedulers into files.
 *
 * @param e the #scheduler
 */
void scheduler_dump_queues(struct engine *e) {

  struct scheduler *s = &e->sched;
  char dumpfile[35];

#ifdef WITH_MPI
  /* Open a file per rank and write the header. Use per rank to avoid MPI
   * calls that can interact with other blocking ones.  */
  snprintf(dumpfile, sizeof(dumpfile), "queue_dump_MPI-step%d.dat_%d", e->step,
           e->nodeID);
#else
  snprintf(dumpfile, sizeof(dumpfile), "queue_dump-step%d.dat", e->step);
#endif

  FILE *file_thread = fopen(dumpfile, "w");
  if (file_thread == NULL) error("Could not create file '%s'.", dumpfile);
  fprintf(file_thread, "# rank queue index type subtype weight\n");
  for (int l = 0; l < s->nr_queues; l++) {
    queue_dump(engine_rank, l, file_thread, &s->queues[l]);
  }
  fclose(file_thread);
}

void scheduler_report_task_times_mapper(void *map_data, int num_elements,
                                        void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  float time_local[task_category_count] = {0};
  float *time_global = (float *)extra_data;

  /* Gather the times spent in the different task categories */
  for (int i = 0; i < num_elements; ++i) {

    const struct task *t = &tasks[i];
    const float total_time = clocks_from_ticks(t->total_ticks);
    const enum task_categories cat = task_get_category(t);
    time_local[cat] += total_time;
  }

  /* Update the global counters */
  for (int i = 0; i < task_category_count; ++i) {
    atomic_add_f(&time_global[i], time_local[i]);
  }
}

/**
 * @brief Display the time spent in the different task categories.
 *
 * @param s The #scheduler.
 * @param nr_threads The number of threads used in the engine.
 */
void scheduler_report_task_times(const struct scheduler *s,
                                 const int nr_threads) {

  const ticks tic = getticks();

  /* Total CPU time spent in engine_launch() */
  const float total_tasks_time = clocks_from_ticks(s->total_ticks) * nr_threads;

  if (total_tasks_time > 0.) {

    /* Initialise counters */
    float time[task_category_count] = {0};
    threadpool_map(s->threadpool, scheduler_report_task_times_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   time);

    /* Compute the dead time */
    float total_time = 0.;
    for (int i = 0; i < task_category_count; ++i) {
      total_time += time[i];
    }
    const float dead_time = total_tasks_time - total_time;

    message("*** CPU time spent in different task categories:");
    for (int i = 0; i < task_category_count; ++i) {
      message("*** %20s: %8.2f %s (%.2f %%)", task_category_names[i], time[i],
              clocks_getunit(), time[i] / total_tasks_time * 100.);
    }
    message("*** %20s: %8.2f %s (%.2f %%)", "dead time", dead_time,
            clocks_getunit(), dead_time / total_tasks_time * 100.);
    message("*** %20s: %8.2f %s (%.2f %%)", "total", total_tasks_time,
            clocks_getunit(), total_tasks_time / total_tasks_time * 100.);
  }

  /* Done. Report the time spent doing this analysis */
  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());
}

void scheduler_collect_task_times_this_step_mapper(void *map_data,
                                                   int num_elements,
                                                   void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  double time_local[task_category_count] = {0};
  double *time_global = (double *)extra_data;

  /* Gather the times spent in the different task categories */
  for (int i = 0; i < num_elements; ++i) {

    struct task *t = &tasks[i];
    const double dt = clocks_diff_ticks(t->toc, t->tic);
    const enum task_categories cat = task_get_category(t);
    time_local[cat] += dt;
    /* Here we want task times of each step, not throughout the
     * global runtime. So we need to reset the counters. */
    t->tic = 0;
    t->toc = 0;
  }

  /* Update the global counters */
  for (int i = 0; i < task_category_count; ++i) {
    atomic_add_d(&time_global[i], time_local[i]);
  }
}

/**
 * @brief Display the time spent in the different task categories.
 *
 * @param s The #scheduler.
 * @param e The #engine
 * @param nr_threads The number of threads used in the engine.
 * @param sub_cycle Whether this is called for an RT sub-cycle
 */
void scheduler_collect_task_times_this_step(const struct scheduler *s,
                                            struct engine *e,
                                            const int nr_threads,
                                            const int sub_cycle) {

  const ticks tic = getticks();

  /* Total CPU time spent in engine_launch() */
  const float total_tasks_time = clocks_from_ticks(s->total_ticks) * nr_threads;

  if (total_tasks_time > 0.) {

    /* Write data into the engine arrays. */
    /* Initialise counters */
    double time[task_category_count] = {0};
    threadpool_map(s->threadpool, scheduler_collect_task_times_this_step_mapper,
                   s->tasks, s->nr_tasks, sizeof(struct task),
                   threadpool_auto_chunk_size, time);

    /* First we write everything into the sub-cycle array.
     * We transfer it later to the total array. */
    for (int i = 0; i < task_category_count; i++) {
      e->local_task_timings_sub_cycle[i] = time[i];
    }

#ifdef WITH_MPI
    float task_timings_buf[task_category_count] = {0.f};
    int test =
        MPI_Reduce(e->local_task_timings_sub_cycle, &task_timings_buf,
                   task_category_count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (test != MPI_SUCCESS) error("MPI reduce failed");

    /* Write result back into correct place. */
    for (int i = 0; i < task_category_count; i++) {
      e->local_task_timings_sub_cycle[i] = task_timings_buf[i];
    }
#endif

    /* Add sub-cycle times to total step times. */
    for (int i = 0; i < task_category_count; i++) {
      e->local_task_timings[i] += e->local_task_timings_sub_cycle[i];
    }
  }

  /* Done. Report the time spent doing this analysis */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
