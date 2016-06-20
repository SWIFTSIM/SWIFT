/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_TASK_H
#define SWIFT_TASK_H

/* Includes. */
#include "cell.h"
#include "cycle.h"

/* Some constants. */
#define task_maxwait 3
#define task_maxunlock 15

/* The different task types. */
enum task_types {
  task_type_none = 0,
  task_type_sort,
  task_type_self,
  task_type_pair,
  task_type_sub_self,
  task_type_sub_pair,
  task_type_init,
  task_type_ghost,
  task_type_drift,
  task_type_kick,
  task_type_kick_fixdt,
  task_type_send,
  task_type_recv,
  task_type_grav_pp,
  task_type_grav_mm,
  task_type_grav_up,
  task_type_grav_down,
  task_type_grav_external,
  task_type_part_sort,
  task_type_gpart_sort,
  task_type_split_cell,
  task_type_rewait,
  task_type_count
};

extern const char *taskID_names[];

/* The different task sub-types. */
enum task_subtypes {
  task_subtype_none = 0,
  task_subtype_density,
  task_subtype_force,
  task_subtype_grav,
  task_subtype_tend,
  task_subtype_count
};

extern const char *subtaskID_names[];

/* Data of a task. */
struct task {

  enum task_types type;
  enum task_subtypes subtype;
  char skip, tight, implicit;
  int flags, wait, rank, weight;

  struct cell *ci, *cj;

  void *buff;

#ifdef WITH_MPI
  MPI_Request req;
#endif

  int rid, last_rid;
  ticks tic, toc;

  int nr_unlock_tasks;
  struct task **unlock_tasks;
};

/* Function prototypes. */
void task_unlock(struct task *t);
float task_overlap(const struct task *ta, const struct task *tb);
int task_lock(struct task *t);
void task_print_mask(unsigned int mask);
void task_print_submask(unsigned int submask);
void task_do_rewait(struct task *t);

#endif /* SWIFT_TASK_H */
