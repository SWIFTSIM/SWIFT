/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
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
#ifndef SWIFT_SINK_BASIC_DEBUG_H
#define SWIFT_SINK_BASIC_DEBUG_H

__attribute__((always_inline)) INLINE static void sink_debug_particle(
    const struct part* p, const struct xpart* xp) {

  warning("[PID%lld] sink_part_data:", p->id);
  warning("[PID%lld] swallow_id = %lld", p->id, p->sink_data.swallow_id);
}

#endif /* SWIFT_SINK_BASIC_DEBUG_H */
