/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_SINK_H
#define SWIFT_SINK_H

/* Config parameters. */
#include <config.h>

/* Select the correct sink model */
#if defined(SINK_NONE)
#include "./sink/Default/sink.h"
#elif defined(SINK_BASIC)
#include "./sink/Basic/sink.h"
#elif defined(SINK_GEAR)
#include "./sink/GEAR/sink.h"
#else
#error "Invalid choice of sink model"
#endif

struct engine;
struct space;

void sink_exact_density_compute(struct space *s, const struct engine *e);
void sink_exact_density_check(struct space *s, const struct engine *e,
                              const double rel_tol);

#endif /* SWIFT_SINK_H */
