/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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

/* Include header */
#include "radiation.h"

float radiation_get_star_ionisation_rate(const struct spart* sp) {


  return 0.0;
}

float radiation_get_part_rate_to_fully_ionize(const struct part* p, const struct xpart* xp) {

  return 0.0;
}

void radiation_consume_ionizing_photons(struct spart* sp, float Delta_dot_N_ion) {
  sp->feedback_data.radiation.dot_N_ion -= Delta_dot_N_ion;
  return;
}

void radiation_tag_part_as_ionized(struct part* p, struct xpart* xpj) {
  xpj->feedback_data.radiation.is_ionized = 1;
  return;
}

int radiation_is_part_ionized(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp) {

  /* Is T > 10^4 K ? */
  const float T =  cooling_get_temperature(phys_const, hydro_props, us,
					   cosmo, cooling, p, xp);
  const float ten_to_four_kelvin = 1e4 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Is the particle ionized ? */
  return (T > ten_to_four_kelvin || xp->feedback_data.radiation.is_ionized);
}

float radiation_consume_ionizing_photons(struct spart* sp, float Delta_dot_N_ion) {
  return 0.0;
}

int radiation_is_part_ionized(const struct part* p, const struct xpart* xpj) {
  return 0;
}
