/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will J. Roper (w.roper@sussex.ac.uk)
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

/* This object's header. */
#include "zoom.h"

/* Local headers. */
#include "cell.h"
#include "scheduler.h"
#include "task.h"

/**
 * @brief Check if a task involves only zoom cells.
 *
 * @param t The task to check.
 *
 * @return 1 if the task involves only zoom cells, 0 otherwise.
 */
static int zoom_scheduler_zoom_only_task(const struct task *t) {
  return (t->ci->type == cell_type_zoom &&
          (t->cj == NULL || t->cj->type == cell_type_zoom));
}

/**
 * @brief Check if a task involves any zoom cells.
 *
 * @param t The task to check.
 *
 * @return 1 if the task involves any zoom cells, 0 otherwise.
 */
static int zoom_scheduler_zoom_task(const struct task *t) {
  return (t->ci->type == cell_type_zoom ||
          (t->cj != NULL && t->cj->type == cell_type_zoom));
}

/**
 * @brief Check if a task involves only background cells.
 *
 * @param t The task to check.
 *
 * @return 1 if the task involves only background cells, 0 otherwise.
 */
static int zoom_scheduler_bkg_only_task(const struct task *t) {
  return (t->ci->type == cell_type_bkg &&
          (t->cj == NULL || t->cj->type == cell_type_bkg));
}

/**
 * @brief Check if a task involves any background cells.
 *
 * @param t The task to check.
 *
 * @return 1 if the task involves any background cells, 0 otherwise.
 */
static int zoom_scheduler_bkg_task(const struct task *t) {
  return (t->ci->type == cell_type_bkg ||
          (t->cj != NULL && t->cj->type == cell_type_bkg));
}

/**
 * @brief #threadpool_map function which runs through the task graph and
 *        gathers the time spent in the different task categories where tasks
 *        involve only zoom cells.
 *
 * @param map_data the index of tasks in this pool thread.
 * @param num_elements the number of indexes in this pool thread
 * @param extra_data The array to store the times in.
 */
void zoom_scheduler_report_zoom_task_times_mapper(void *map_data,
                                                  int num_elements,
                                                  void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  float time_local[task_category_count] = {0};
  float *time_global = (float *)extra_data;

  /* Gather the times spent in the different task categories */
  for (int i = 0; i < num_elements; ++i) {

    /* Get the task */
    const struct task *t = &tasks[i];

    /* Skip, skipped tasks (but we can't use the skip flag because all tasks
     * have been executed at this point and the skip flag is True for tasks
     * that were executed). */
    if (t->type == task_type_none) continue;

    /* Did it only involve zoom cells? */
    if (zoom_scheduler_zoom_only_task(t)) {

      const float total_time = clocks_from_ticks(t->total_ticks);
      const enum task_categories cat = task_get_category(t);
      time_local[cat] += total_time;
    }
  }

  /* Update the global counters */
  for (int i = 0; i < task_category_count; ++i) {
    atomic_add_f(&time_global[i], time_local[i]);
  }
}

/**
 * @brief #threadpool_map function which runs through the task graph and
 *       gathers the time spent in the different task categories where tasks
 *       involve only background cells.
 *
 * @param map_data the index of tasks in this pool thread.
 * @param num_elements the number of indexes in this pool thread
 * @param extra_data The array to store the times in.
 */
void zoom_scheduler_report_bkg_task_times_mapper(void *map_data,
                                                 int num_elements,
                                                 void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  float time_local[task_category_count] = {0};
  float *time_global = (float *)extra_data;

  /* Gather the times spent in the different task categories */
  for (int i = 0; i < num_elements; ++i) {

    /* Get the task */
    const struct task *t = &tasks[i];

    /* Skip, skipped tasks (but we can't use the skip flag because all tasks
     * have been executed at this point and the skip flag is True for tasks
     * that were executed). */
    if (t->type == task_type_none) continue;

    /* Did it only involve background cells? */
    if (zoom_scheduler_bkg_only_task(t)) {

      const float total_time = clocks_from_ticks(t->total_ticks);
      const enum task_categories cat = task_get_category(t);
      time_local[cat] += total_time;
    }
  }

  /* Update the global counters */
  for (int i = 0; i < task_category_count; ++i) {
    atomic_add_f(&time_global[i], time_local[i]);
  }
}

/**
 * @brief #threadpool_map function which runs through the task graph and
 *      gathers the time spent in the different task categories where tasks
 *      involve zoom and background cells.
 *
 * @param map_data the index of tasks in this pool thread.
 * @param num_elements the number of indexes in this pool thread
 * @param extra_data The array to store the times in.
 */
void zoom_scheduler_report_mixed_task_times_mapper(void *map_data,
                                                   int num_elements,
                                                   void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  float time_local[task_category_count] = {0};
  float *time_global = (float *)extra_data;

  /* Gather the times spent in the different task categories */
  for (int i = 0; i < num_elements; ++i) {

    /* Get the task */
    const struct task *t = &tasks[i];

    /* Skip, skipped tasks (but we can't use the skip flag because all tasks
     * have been executed at this point and the skip flag is True for tasks
     * that were executed). */
    if (t->type == task_type_none) continue;

    /* Did it involve zoom and background cells? */
    if (zoom_scheduler_zoom_task(t) && zoom_scheduler_bkg_task(t)) {

      const float total_time = clocks_from_ticks(t->total_ticks);
      const enum task_categories cat = task_get_category(t);
      time_local[cat] += total_time;
    }
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
void zoom_scheduler_report_task_times(const struct scheduler *s,
                                      const int nr_threads) {

  const ticks tic = getticks();

  /* Total CPU time spent in engine_launch() */
  const float total_tasks_time = clocks_from_ticks(s->total_ticks) * nr_threads;

  if (total_tasks_time > 0.) {

    /* Initialise counters */
    float zoom_time[task_category_count] = {0};
    float bkg_time[task_category_count] = {0};
    float mixed_time[task_category_count] = {0};

    /* Gather the times spent in the different task categories */
    threadpool_map(s->threadpool, zoom_scheduler_report_zoom_task_times_mapper,
                   s->tasks, s->nr_tasks, sizeof(struct task),
                   threadpool_auto_chunk_size, zoom_time);
    threadpool_map(s->threadpool, zoom_scheduler_report_bkg_task_times_mapper,
                   s->tasks, s->nr_tasks, sizeof(struct task),
                   threadpool_auto_chunk_size, bkg_time);
    threadpool_map(s->threadpool, zoom_scheduler_report_mixed_task_times_mapper,
                   s->tasks, s->nr_tasks, sizeof(struct task),
                   threadpool_auto_chunk_size, mixed_time);

    /* Compute the total time spent in zoom only task, bkg only tasks, mixed
     * tasks. We ignore deadtime since that is reported in the standard
     * scheduler report. */
    float total_zoom_time = 0.f;
    float total_bkg_time = 0.f;
    float total_mixed_time = 0.f;
    for (int i = 0; i < task_category_count; ++i) {
      total_zoom_time += zoom_time[i];
      total_bkg_time += bkg_time[i];
      total_mixed_time += mixed_time[i];
    }

    message(
        "*** CPU time spent in different zoom task categories "
        "(within %% / total %%):");
    for (int i = 0; i < task_category_count; ++i) {
      if (zoom_time[i] == 0.f) continue;
      message("*** %20s: %12.2f %s (%6.2f %% / %5.2f %%)",
              task_category_names[i], zoom_time[i], clocks_getunit(),
              zoom_time[i] / total_zoom_time * 100.,
              zoom_time[i] / total_tasks_time * 100.);
    }
    message(
        "*** CPU time spent in different background task categories "
        "(within %% / total %%):");
    for (int i = 0; i < task_category_count; ++i) {
      if (bkg_time[i] == 0.f) continue;
      message("*** %20s: %12.2f %s (%6.2f %% / %5.2f %%)",
              task_category_names[i], bkg_time[i], clocks_getunit(),
              bkg_time[i] / total_bkg_time * 100.,
              bkg_time[i] / total_tasks_time * 100.);
    }
    message(
        "*** CPU time spent in different mixed (zoom + background) task "
        "categories (within %% / total %%):");
    for (int i = 0; i < task_category_count; ++i) {
      if (mixed_time[i] == 0.f) continue;
      message("*** %20s: %12.2f %s (%6.2f %% / %5.2f %%)",
              task_category_names[i], mixed_time[i], clocks_getunit(),
              mixed_time[i] / total_mixed_time * 100.,
              mixed_time[i] / total_tasks_time * 100.);
    }
    message("*** Fraction of total CPU time spent in different combinations:");
    message("*** %20s: %12.2f %s (%5.2f %%)", "zoom total", total_zoom_time,
            clocks_getunit(), total_zoom_time / total_tasks_time * 100.);
    message("*** %20s: %12.2f %s (%5.2f %%)", "background total",
            total_bkg_time, clocks_getunit(),
            total_bkg_time / total_tasks_time * 100.);
    message("*** %20s: %12.2f %s (%5.2f %%)", "mixed total", total_mixed_time,
            clocks_getunit(), total_mixed_time / total_tasks_time * 100.);
  }

  /* Done. Report the time spent doing this analysis */
  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());
}
