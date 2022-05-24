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
#ifndef SWIFT_CHEMISTRY_GEAR_DEBUG_H
#define SWIFT_CHEMISTRY_GEAR_DEBUG_H

__attribute__((always_inline)) INLINE static void chemistry_debug_particle(
    const struct part* p, const struct xpart* xp) {

  warning("[PID%lld] chemistry_part_data:", p->id);
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    warning("[PID%lld] metal_mass[%i]: %.3e", p->id, i,
            p->chemistry_data.metal_mass[i]);
  }
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    warning("[PID%lld] smoothed_metal_mass_fraction[%i]: %.3e", p->id, i,
            p->chemistry_data.smoothed_metal_mass_fraction[i]);
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_DEBUG_H */
