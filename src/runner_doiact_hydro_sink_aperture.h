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

/* Gas-sink neighbour loop for sink formation (task subtype
 * sink_formation_sink).
 *
 * A gas particle cannot form a sink if the new sink would overlap the
 * accretion sphere of an existing sink. This loop finds the existing sinks
 * near each active gas particle and calls
 * sink_prepare_part_sink_formation_sink_criteria() to record the overlap on
 * the gas particle. It replaces the old scan of all sinks in the space in
 * runner_do_prepare_part_sink_formation() (runner_sinks.c), which is too slow
 * with many reserved sink slots. The old scan is still used when the sink
 * cut-off radius is not fixed.
 *
 * The loop is a plain double loop over the gas and the sinks of a cell. It
 * does not use the sorted indices.
 *
 * Two limits:
 * - The overlap test adds two radii, so the reach is 2 * r_cut. The cell
 *   pairs come from the hydro task splitting, which only uses r_cut. A sink
 *   between r_cut and 2 * r_cut away in a cell that has no pair task with the
 *   gas cell is not seen.
 * - The loop sees the sinks that exist at the start of the step. Two gas
 *   particles can form sinks in the same step even if the new sinks overlap.
 *   The sink merging removes the overlap later.
 */

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
