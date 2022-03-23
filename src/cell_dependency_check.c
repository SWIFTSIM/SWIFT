/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#include "cell.h"
#include "error.h"

void cell_recursively_count_tasks(struct cell *c, const struct task *t) {
#ifdef SWIFT_DEBUG_CHECKS
  const int itask = t->type * task_subtype_count + t->subtype;
  c->tasks_executed[itask]++;
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_recursively_count_tasks(c->progeny[k], t);
      }
    }
  }
#endif
}

void cell_check_dependencies(const struct cell *c, const struct task *t,
                             const int *task_mask) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int itype = 0; itype < task_type_count; itype++) {
    for (int istype = 0; istype < task_subtype_count; istype++) {
      const int itask = itype * task_subtype_count + istype;
      if (task_mask[itask]) {
        if (c->tasks_executed[itask] > 0) {
          char nameA[100], nameB[100];
          task_get_full_name(itype, istype, nameA);
          task_get_full_name(t->type, t->subtype, nameB);
          error("Task dependency violated (%s has run before %s)!", nameA,
                nameB);
        }
      }
    }
  }
#endif
}

void cell_recursively_check_task_mask(const struct cell *c,
                                      const struct task *t,
                                      const int *task_mask) {
#ifdef SWIFT_DEBUG_CHECKS
  cell_check_dependencies(c, t, task_mask);
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_recursively_check_task_mask(c->progeny[k], t, task_mask);
      }
    }
  }
#endif
}
