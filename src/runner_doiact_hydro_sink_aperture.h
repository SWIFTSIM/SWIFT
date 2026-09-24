/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_RUNNER_DOIACT_HYDRO_SINK_APERTURE_H
#define SWIFT_RUNNER_DOIACT_HYDRO_SINK_APERTURE_H

/* Gas-vs-existing-sink neighbour loop for sink formation's "overlapping
 * sink" criterion (task subtype sink_formation_sink).
 *
 * Finds, for every active gas particle, the existing sink particles that
 * could make it ineligible to form a new sink (its would-be accretion
 * sphere would overlap an existing sink's accretion sphere), and calls the
 * existing sink_prepare_part_sink_formation_sink_criteria() (src/sink/GEAR/
 * sink.h) to record the result on the gas particle. This is the task-graph
 * replacement for the O(N_gas_eligible * N_sink_in_space) brute-force loop
 * that runner_do_prepare_part_sink_formation() (runner_sinks.c) falls back
 * to when this loop is not active.
 *
 * Unlike the gas-gas fixed-aperture loop (runner_doiact_functions_hydro_
 * aperture.h), this is a genuinely mixed-type search (gas queries sinks, not
 * gas queries gas), so it does not reuse that file's macro-expansion
 * machinery. It also does not need that file's active-particle index-list
 * optimisation (indt_stack): the cost here is O(N_gas_active * N_sink_in_
 * cell), already linear in the gas count because the inner loop walks the
 * handful of sink particles (real or reserved-but-unformed) local to a
 * cell, not a second gas array. Sorted pair scan bounds are likewise not
 * worth it for the same reason, so the pair function is naive.
 *
 * Geometric caveat (see runner_sinks.c): the true "no overlap" test is a SUM
 * of two independent radii (the candidate's own would-be accretion radius,
 * r_cut, plus the existing sink's actual accretion radius, sink->h *
 * kernel_gamma, which under sink_props->use_fixed_r_cut -- the only mode
 * this loop is active for -- is always exactly r_cut too, hence a true
 * reach of 2 * r_cut). DOSUB recursion accounts for this (it stops at
 * 2 * r_cut < 0.5 * dmin, not r_cut, so it never prunes a sub-cell pair
 * still within the true reach -- the leaf is a plain double loop with no
 * further pruning, so recursing less loses nothing). What DOSUB cannot fix
 * is which *hydro.super* cells get a pair task connecting them at all: that
 * stencil is fixed by scheduler_splittasks.c's cell_can_split_{self,pair}_
 * hydro_task, which folds in r_cut (not 2 * r_cut) to decide when a cell is
 * small enough to stop splitting, so the smallest hydro.super cells end up
 * with dmin only a little above r_cut. A sink two hydro.super cells away
 * from the candidate gas particle can then be within the true 2 * r_cut
 * reach with no task connecting the two cells, and is not visited. Closing
 * that gap needs scheduler_splittasks.c to fold in 2 * r_cut instead of
 * r_cut, which would coarsen the shared hydro decomposition in exactly the
 * dense regions where sinks form; scheduler_splittasks.c is shared with the
 * gas-gas loop and every other hydro/stars/sink task, so that tradeoff
 * applies globally, not just to this search. There is no fallback for that
 * residual case once sink_props->use_fixed_r_cut retires the brute-force
 * scan (see runner_sinks.c). */

void runner_doself1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *c, const float r_cut);

void runner_dopair1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *ci, struct cell *cj);

void runner_dosub_self1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *c, const float r_cut, const int gettimer);

void runner_dosub_pair1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *ci, struct cell *cj, const float r_cut,
    const int gettimer);

#endif /* SWIFT_RUNNER_DOIACT_HYDRO_SINK_APERTURE_H */
