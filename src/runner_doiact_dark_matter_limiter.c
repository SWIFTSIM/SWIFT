/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "runner.h"
#include "space_getsid.h"
#include "timers.h"
#include "timestep_limiter_iact.h"

/* Import the limiter loop functions. */
#define FUNCTION dm_limiter
#include "runner_doiact_functions_dark_matter_limiter.h"
#undef FUNCTION

