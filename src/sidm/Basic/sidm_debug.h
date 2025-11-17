/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_SIDM_BASIC_DEBUG_H
#define SWIFT_SIDM_BASIC_DEBUG_H

__attribute__((always_inline)) INLINE static void sidm_debug_particle(
    const struct sipart *sip) {
  warning(
      "[PID%lld] x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e] sip->mass=%.3e",
      sip->id, sip->x[0], sip->x[1], sip->x[2], sip->v[0], sip->v[1], sip->v[2],
      sip->mass);
}
#endif /* SWIFT_SIDM_BASIC_DEBUG_H */
