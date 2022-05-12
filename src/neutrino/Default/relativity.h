/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_NEUTRINO_RELATIVITY_H
#define SWIFT_DEFAULT_NEUTRINO_RELATIVITY_H

/**
 * @brief Calculate the relativistic correction to the 'drift' timestep
 *
 * @param v Spatial part of the 4-velocity
 * @param a Scale-factor
 * @param c Speed of light
 */
INLINE static float relativistic_drift_factor(const float *v, float a,
                                              float c) {
  const float v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  const float ac = a * c;
  const float ac2 = ac * ac;

  return ac / sqrtf(ac2 + v2);
}

#endif /* SWIFT_DEFAULT_NEUTRINO_RELATIVITY_H */
