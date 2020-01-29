/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TIMESTEP_SYNC_PART_H
#define SWIFT_TIMESTEP_SYNC_PART_H

/* Config parameters. */
#include "../config.h"

/**
 * @brief Flag a particle for synchronization on the time-line.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void timestep_sync_part(
    struct part *p) {

  p->limiter_data.to_be_synchronized = 1;
}

#endif /* SWIFT_TIMESTEP_SYNC_PART_H */
