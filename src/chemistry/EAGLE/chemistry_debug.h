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
#ifndef SWIFT_CHEMISTRY_EAGLE_DEBUG_H
#define SWIFT_CHEMISTRY_EAGLE_DEBUG_H

__attribute__((always_inline)) INLINE static void chemistry_debug_particle(
    const struct part* p, const struct xpart* xp) {

  warning("[PID%lld chemistry_part_data:", p->id);
  for (int i = 0; i < chemistry_element_count; i++) {
    warning("[PID%lld metal_mass_fraction[%i]: %.3e", p->id, i,
            p->chemistry_data.metal_mass_fraction[i]);
  }
  warning("[PID%lld metal_mass_fraction_total: %.3e", p->id,
          p->chemistry_data.metal_mass_fraction_total);
  for (int i = 0; i < chemistry_element_count; i++) {
    warning("[PID%lld smoothed_metal_mass_fraction[%i]: %.3e", p->id, i,
            p->chemistry_data.smoothed_metal_mass_fraction[i]);
  }
  warning("[PID%lld metal_mass_fraction_total: %.3e", p->id,
          p->chemistry_data.smoothed_metal_mass_fraction_total);
  warning(
      "[PID%lld mass_from_SNIa: %.3e, metal_mass_fraction_from_SNIa: %.3e, "
      "mass_from_AGB: %.3e, "
      "metal_mass_fraction_from_AGB: %.3e, "
      "mass_from_SNII: %.3e, "
      "metal_mass_fraction_from_SNII: %.3e, "
      "iron_mass_fraction_from_SNIa: %.3e, "
      "smoothed_iron_mass_fraction_from_SNIa: %.3e",
      p->id, p->chemistry_data.mass_from_SNIa,
      p->chemistry_data.metal_mass_fraction_from_SNIa,
      p->chemistry_data.mass_from_AGB,
      p->chemistry_data.metal_mass_fraction_from_AGB,
      p->chemistry_data.mass_from_SNII,
      p->chemistry_data.metal_mass_fraction_from_SNII,
      p->chemistry_data.iron_mass_fraction_from_SNIa,
      p->chemistry_data.smoothed_iron_mass_fraction_from_SNIa);
}

#endif /* SWIFT_CHEMISTRY_EAGLE_DEBUG_H */
