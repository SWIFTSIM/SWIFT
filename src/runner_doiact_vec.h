/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
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

#ifndef SWIFT_RUNNER_VEC_H
#define SWIFT_RUNNER_VEC_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "hydro.h"
#include "part.h"
#include "runner.h"
#include "timers.h"
#include "vector.h"

/* Function prototypes. */
void runner_doself1_density_vec(struct runner *r, struct cell *restrict c);
void runner_dopair1_density_vec(struct runner *r, struct cell *restrict ci,
                                struct cell *restrict cj, const int sid,
                                const double *shift);

#endif /* SWIFT_RUNNER_VEC_H */
