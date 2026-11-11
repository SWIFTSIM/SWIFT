/*******************************************************************************
 * This file is part of SWIFT.
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
#include <config.h>

/* Local includes. */
#include "swift.h"

#ifdef COOLING_GEAJAY

/**
 * @brief Check that a split divides the radiated energy and keeps the rest.
 */
int main(int argc, char **argv) {

  struct part p;
  struct xpart xp;
  bzero(&p, sizeof(p));
  bzero(&xp, sizeof(xp));

  xp.cooling_data.radiated_energy = 8.f;
  xp.cooling_data.time_last_event = 3.;
#if COOLING_GEAJAY_MODE >= 1
  xp.cooling_data.HI_frac = 0.25f;
  xp.cooling_data.cosmic_ray_value = 7.f;
#endif

  cooling_split_part(&p, &xp, 4.);

  if (xp.cooling_data.radiated_energy != 2.f)
    error("The radiated energy is %g, expected 2.",
          xp.cooling_data.radiated_energy);
  if (xp.cooling_data.time_last_event != 3.)
    error("The time of the last event changed.");
#if COOLING_GEAJAY_MODE >= 1
  if (xp.cooling_data.HI_frac != 0.25f)
    error("A mass fraction changed in the split.");
  if (xp.cooling_data.cosmic_ray_value != 7.f)
    error("The cosmic-ray value changed in the split.");
#endif

  return 0;
}

#else

int main(int argc, char **argv) { return 0; }

#endif
