/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2026 Will J. Roper (w.roper@sussex.ac.uk)
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "scheduler.h"

/* Local headers. */
#include "cycle.h"
#include "engine.h"
#include "error.h"
#include "queue.h"
#include "version.h"

/* Conservative number of dependencies per task type. */
#define MAX_NUMBER_DEP 128

/**
 * @brief Information about all the task dependencies of a single task.
 */
struct task_dependency {
  /* Main task */
  int type_in;
  int subtype_in;
  int implicit_in;
  int task_in_is_top;
  int task_in_is_grav_super;
  int task_in_is_hydro_super;

  /* Dependent task */
  int type_out[MAX_NUMBER_DEP];
  int subtype_out[MAX_NUMBER_DEP];
  int implicit_out[MAX_NUMBER_DEP];
  int task_out_is_top[MAX_NUMBER_DEP];
  int task_out_is_grav_super[MAX_NUMBER_DEP];
  int task_out_is_hydro_super[MAX_NUMBER_DEP];

  /* Statistics */
  int number_link[MAX_NUMBER_DEP];
  int number_rank[MAX_NUMBER_DEP];
};

#ifdef WITH_MPI

/**
 * @brief Define the #task_dependency for MPI.
 *
 * @param tstype The MPI_Datatype to initialize.
 */
static void task_dependency_define(MPI_Datatype *tstype) {
  const int count = 14;
  int blocklens[count];
  MPI_Datatype types[count];
  MPI_Aint disps[count];

  for (int i = 0; i < count; i++) {
    types[i] = MPI_INT;
  }

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

  disps[12] = offsetof(struct task_dependency, number_link);
  blocklens[12] = MAX_NUMBER_DEP;
  disps[13] = offsetof(struct task_dependency, number_rank);
  blocklens[13] = MAX_NUMBER_DEP;

  MPI_Type_create_struct(count, blocklens, disps, types, tstype);
  MPI_Type_commit(tstype);
}

/**
 * @brief Sum operator of #task_dependency for MPI.
 *
 * @param in_p The #task_dependency to add.
 * @param out_p The #task_dependency where in_p is added.
 * @param len The length of the arrays.
 * @param type The MPI datatype.
 */
static void task_dependency_sum(void *in_p, void *out_p, int *len,
                                MPI_Datatype *type) {
  struct task_dependency *in = (struct task_dependency *)in_p;
  struct task_dependency *out = (struct task_dependency *)out_p;

  for (int i = 0; i < *len; i++) {
    for (int j = 0; j < MAX_NUMBER_DEP; j++) {
      if (in[i].number_link[j] == -1) {
        break;
      }

      const int tb_type = in[i].type_out[j];
      const int tb_subtype = in[i].subtype_out[j];

#ifdef SWIFT_DEBUG_CHECKS
      if (tb_type >= task_type_count) {
        error("Unknown task type %i", tb_type);
      }

      if (tb_subtype >= task_subtype_count) {
        error("Unknown subtask type %i", tb_subtype);
      }
#endif

      int k = 0;
      while (k < MAX_NUMBER_DEP) {
        if (out[i].number_link[k] == -1) {
          out[i].number_link[k] = 0;
          out[i].number_rank[k] = 0;

          out[i].type_in = in[i].type_in;
          out[i].subtype_in = in[i].subtype_in;
          out[i].implicit_in = in[i].implicit_in;

          out[i].type_out[k] = in[i].type_out[j];
          out[i].subtype_out[k] = in[i].subtype_out[j];
          out[i].implicit_out[k] = in[i].implicit_out[j];
          break;
        }

        if (out[i].type_out[k] == tb_type &&
            out[i].subtype_out[k] == tb_subtype) {
          break;
        }

        k++;
      }

      if (k == MAX_NUMBER_DEP) {
        error("Not enough memory, please increase MAX_NUMBER_DEP");
      }

#ifdef SWIFT_DEBUG_CHECKS
      if (out[i].type_in != in[i].type_in ||
          out[i].subtype_in != in[i].subtype_in ||
          out[i].implicit_in != in[i].implicit_in ||
          out[i].type_out[k] != in[i].type_out[j] ||
          out[i].subtype_out[k] != in[i].subtype_out[j] ||
          out[i].implicit_out[k] != in[i].implicit_out[j]) {
        error("Tasks do not correspond");
      }
#endif

      out[i].number_link[k] += in[i].number_link[j];
      out[i].number_rank[k] += in[i].number_rank[j];

      out[i].task_in_is_top = min(out[i].task_in_is_top, in[i].task_in_is_top);
      out[i].task_in_is_hydro_super =
          min(out[i].task_in_is_hydro_super, in[i].task_in_is_hydro_super);
      out[i].task_in_is_grav_super =
          min(out[i].task_in_is_grav_super, in[i].task_in_is_grav_super);

      out[i].task_out_is_top[j] =
          min(out[i].task_out_is_top[j], in[i].task_out_is_top[j]);
      out[i].task_out_is_hydro_super[j] = min(out[i].task_out_is_hydro_super[j],
                                              in[i].task_out_is_hydro_super[j]);
      out[i].task_out_is_grav_super[j] = min(out[i].task_out_is_grav_super[j],
                                             in[i].task_out_is_grav_super[j]);
    }
  }

  (void)type;
}

#endif

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
