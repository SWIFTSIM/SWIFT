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
#include "../config.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "cache.h"
#include "cell.h"
#include "chemistry.h"
#include "clocks.h"
#include "common_io.h"
#include "const.h"
#include "cooling.h"
#include "cooling_struct.h"
#include "cosmology.h"
#include "cycle.h"
#include "debug.h"
#include "dump.h"
#include "engine.h"
#include "entropy_floor.h"
#include "error.h"
#include "gravity.h"
#include "gravity_derivatives.h"
#include "gravity_properties.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "lock.h"
#include "logger.h"
#include "logger_io.h"
#include "map.h"
#include "memuse.h"
#include "mesh_gravity.h"
#include "multipole.h"
#include "outputlist.h"
#include "parallel_io.h"
#include "parser.h"
#include "part.h"
#include "partition.h"
#include "periodic.h"
#include "physical_constants.h"
#include "potential.h"
#include "profiler.h"
#include "queue.h"
#include "random.h"
#include "restart.h"
#include "runner.h"
#include "scheduler.h"
#include "serial_io.h"
#include "single_io.h"
#include "space.h"
#include "star_formation.h"
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
