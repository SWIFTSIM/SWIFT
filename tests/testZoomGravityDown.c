/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Will J. Roper (w.roper@sussex.ac.uk).
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

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

/**
 * @brief Check that a void down-pass only updates local, active zoom children.
 *
 * Mock two owners without MPI communication. The foreign child has a valid
 * multipole but no locally initialised field tensor; its activity metadata
 * later becomes stale. Neither may be accessed by this rank's down-pass.
 * The local child is its own super, so its field is updated here but its
 * separate down task is not executed. Expected fields are exact constants.
 * Runs through Automake with TESTS=testZoomGravityDown; no ICs or launcher
 * needed.
 */
int main(int argc, char *argv[]) {
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  struct engine e = {0};
  struct cosmology cosmo = {0};
  struct runner r = {0};
  struct cell parent = {0}, local = {0}, foreign = {0};
  struct gravity_tensors parent_multi = {0}, local_multi = {0},
                         foreign_multi = {0};
  e.ti_current = 1;
  e.cosmology = &cosmo;
  r.e = &e;
  parent.type = cell_type_bkg;
  parent.subtype = cell_subtype_void;
  parent.grav.multipole = &parent_multi;
  parent.grav.ti_old_multipole = e.ti_current;
  parent.progeny[0] = &local;
  parent.progeny[1] = &foreign;
  local.type = foreign.type = cell_type_zoom;
  foreign.nodeID = e.nodeID + 1;
  local.grav.multipole = &local_multi;
  foreign.grav.multipole = &foreign_multi;
  local.grav.super = &local;
  foreign.grav.super = &foreign;
  local.grav.count = foreign.grav.count = 1;
  local_multi.m_pole.M_000 = foreign_multi.m_pole.M_000 = 1.f;
  local.grav.ti_old_multipole = foreign.grav.ti_old_multipole = e.ti_current;
  local.grav.ti_end_min = foreign.grav.ti_end_min = e.ti_current;
  gravity_field_tensors_init(&parent_multi.pot, e.ti_current);
  gravity_field_tensors_init(&local_multi.pot, e.ti_current);
  parent_multi.pot.F_000 = 1.f;
  parent_multi.pot.interacted = 1;
#ifdef SWIFT_DEBUG_CHECKS
  parent_multi.pot.num_interacted = 1;
#endif

  runner_do_grav_down(&r, &parent, /*timer=*/0);
  if (local_multi.pot.F_000 != 1.f)
    error("Void down-pass did not update the local zoom child's field");
  if (foreign_multi.pot.F_000 != 0.f || foreign_multi.pot.interacted)
    error("Void down-pass updated a foreign zoom child's field");

  /* An inactive local child must stay unchanged; the stale foreign activity
   * must not even be queried, as it would trigger a debug time-zone error. */
  local.grav.ti_end_min = e.ti_current + 1;
  foreign.grav.ti_end_min = 0;
  runner_do_grav_down(&r, &parent, /*timer=*/0);
  if (local_multi.pot.F_000 != 1.f || foreign_multi.pot.F_000 != 0.f ||
      foreign_multi.pot.interacted)
    error("Void down-pass changed an inactive or foreign child's field");

  message("Local, inactive and foreign zoom-child down-pass checks passed.");
  return 0;
}
