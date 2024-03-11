/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_SWIFT_H
#define SWIFT_SWIFT_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "active.h"
#include "adaptive_softening.h"
#include "atomic.h"
#include "black_holes_properties.h"
#include "cache.h"
#include "cell.h"
#include "chemistry.h"
#include "clocks.h"
#include "common_io.h"
#include "const.h"
#include "cooling.h"
#include "cooling_properties.h"
#include "cosmology.h"
#include "csds.h"
#include "csds_io.h"
#include "cycle.h"
#include "debug.h"
#include "engine.h"
#include "entropy_floor.h"
#include "error.h"
#include "extra_io.h"
#include "feedback.h"
#include "feedback_properties.h"
#include "fof.h"
#include "forcing.h"
#include "gravity.h"
#include "gravity_derivatives.h"
#include "gravity_properties.h"
#include "hashmap.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "ic_info.h"
#include "lightcone/lightcone_array.h"
#include "line_of_sight.h"
#include "lock.h"
#include "map.h"
#include "memuse.h"
#include "mesh_gravity.h"
#include "minmax.h"
#include "mpiuse.h"
#include "multipole.h"
#include "neutrino.h"
#include "output_list.h"
#include "output_options.h"
#include "parallel_io.h"
#include "parser.h"
#include "part.h"
#include "partition.h"
#include "periodic.h"
#include "physical_constants.h"
#include "potential.h"
#include "power_spectrum.h"
#include "pressure_floor.h"
#include "pressure_floor_iact.h"
#include "profiler.h"
#include "queue.h"
#include "random.h"
#include "restart.h"
#include "rt.h"
#include "rt_properties.h"
#include "runner.h"
#include "scheduler.h"
#include "serial_io.h"
#include "single_io.h"
#include "sink_properties.h"
#include "space.h"
#include "star_formation.h"
#include "star_formation_iact.h"
#include "star_formation_logger.h"
#include "stars.h"
#include "stars_io.h"
#include "task.h"
#include "threadpool.h"
#include "timeline.h"
#include "timers.h"
#include "tools.h"
#include "units.h"
#include "velociraptor_interface.h"
#include "version.h"

#endif /* SWIFT_SWIFT_H */
