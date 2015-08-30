/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include "const.h"
#include "error.h"
#include "cycle.h"
#include "timers.h"
#include "const.h"
#include "atomic.h"
#include "lock.h"
#include "task.h"
#include "scheduler.h"
#include "part.h"
#include "multipole.h"
#include "cell.h"
#include "space.h"
#include "queue.h"
#include "runner.h"
#include "engine.h"
#include "units.h"
#include "single_io.h"
#include "serial_io.h"
#include "parallel_io.h"
#include "debug.h"
#include "version.h"

#ifdef LEGACY_GADGET2_SPH
#include "runner_iact_legacy.h"
#else
#include "runner_iact.h"
#endif
#include "runner_iact_grav.h"

#endif /* SWIFT_SWIFT_H */
