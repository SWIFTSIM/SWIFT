/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Unit test of the sink drift flags in a super-cell that holds no sink.
 *
 * A pair loop activates the sink drift of a gas-only sub-cell. The drift of
 * the super-cell then stops early because it has no sink. The test checks
 * that no sink drift flag is left set in the sub-tree, and that activating
 * the same sub-cell on the next step activates the drift of the super-cell
 * again, so that the sub-cell is drifted to the current time.
 */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <string.h>

/* Local headers. */
#include "swift.h"

static int sink_drift_flags_set(const struct cell *c) {
  return cell_get_flag(c, cell_flag_do_sink_drift) ||
         cell_get_flag(c, cell_flag_do_sink_sub_drift);
}

int main(int argc, char *argv[]) {

#ifdef SINK_NONE
  /* Without a sink model the sink fields of a cell share one union. */
  message("No sink model: nothing to test.");
  return 0;
#endif

  /* A super-cell with one gas-only leaf below an intermediate cell. */
  struct cell super, mid, leaf;
  memset(&super, 0, sizeof(struct cell));
  memset(&mid, 0, sizeof(struct cell));
  memset(&leaf, 0, sizeof(struct cell));

  super.split = 1;
  super.progeny[0] = &mid;
  mid.split = 1;
  mid.parent = &super;
  mid.progeny[0] = &leaf;
  leaf.parent = &mid;
  super.hydro.super = &super;
  mid.hydro.super = &super;
  leaf.hydro.super = &super;

  /* The sink drift task of the super-cell. */
  struct task tasks[1];
  memset(tasks, 0, sizeof(tasks));
  tasks[0].type = task_type_drift_sink;
  tasks[0].ci = &super;
  tasks[0].skip = 1;
  super.sinks.drift = &tasks[0];

  int tid_active[1] = {0};
  struct scheduler sched;
  memset(&sched, 0, sizeof(struct scheduler));
  sched.tasks = tasks;
  sched.tid_active = tid_active;

  struct space s;
  memset(&s, 0, sizeof(struct space));
  s.periodic = 1;
  s.dim[0] = s.dim[1] = s.dim[2] = 1.;

  struct hydro_props hydro_props;
  memset(&hydro_props, 0, sizeof(struct hydro_props));
  hydro_props.h_max = 1.f;

  struct engine e;
  memset(&e, 0, sizeof(struct engine));
  e.s = &s;
  e.hydro_properties = &hydro_props;
  e.time_base = 1.;
  s.e = &e;

  /* Step 1: activate and drift. */
  e.ti_current = 1;
  cell_activate_drift_sink(&leaf, &sched);
  if (tasks[0].skip)
    error("Step 1: the sink drift of the super-cell was not activated.");
  cell_drift_sink(&super, &e, /*force=*/0, /*init_particles=*/0);

  if (sink_drift_flags_set(&super) || sink_drift_flags_set(&mid) ||
      sink_drift_flags_set(&leaf))
    error(
        "Step 1: a sink drift flag is still set after the drift (super %d, "
        "mid %d, leaf %d).",
        sink_drift_flags_set(&super), sink_drift_flags_set(&mid),
        sink_drift_flags_set(&leaf));
  if (leaf.sinks.ti_old_part != e.ti_current)
    error("Step 1: the leaf was not drifted.");

  /* Step 2: the task is skipped again, as at the start of a new step. */
  tasks[0].skip = 1;
  sched.active_count = 0;
  e.ti_current = 2;
  cell_activate_drift_sink(&leaf, &sched);
  if (tasks[0].skip)
    error("Step 2: the sink drift of the super-cell was not activated.");
  cell_drift_sink(&super, &e, /*force=*/0, /*init_particles=*/0);
  if (leaf.sinks.ti_old_part != e.ti_current)
    error("Step 2: the leaf was not drifted to the current time.");

  message("The sink drift flags are cleared below an empty super-cell.");

  return 0;
}
