/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_COOLING_PS2020_DEBUG_H
#define SWIFT_COOLING_PS2020_DEBUG_H

__attribute__((always_inline)) INLINE static void cooling_debug_particle(
    const struct part* p, const struct xpart* xp) {

  warning("[PID%lld] cooling_part_data:", p->id);
  warning("[PID%lld] subgrid_temp=%.3e, subgrid_dens=%.3e", p->id,
          p->cooling_data.subgrid_temp, p->cooling_data.subgrid_dens);
  if (xp != NULL) {
    warning("[PID%lld] cooling_xpart_data:", p->id);
    warning("[PID%lld] radiated_energy=%.3e", p->id,
            xp->cooling_data.radiated_energy);
  }
}

#endif /* SWIFT_COOLING_PS2020_DEBUG_H */
