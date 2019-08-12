/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
#ifndef SWIFT_NONE_PRESSURE_FLOOR_IACT_H
#define SWIFT_NONE_PRESSURE_FLOOR_IACT_H

/**
 * @file NONE/pressure_floor_iact.h
 * @brief Density computation
 */

/**
 * @brief do pressure_floor computation after the runner_iact_density (symmetric
 * version)
 *
 * Compute the velocity dispersion follow eq. 2 in Revaz & Jablonka 2018.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_pressure_floor(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {}

/**
 * @brief do pressure_floor computation after the runner_iact_density (non
 * symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_pressure_floor(float r2, const float *dx, float hi, float hj,
                                  struct part *restrict pi,
                                  const struct part *restrict pj, float a,
                                  float H) {}

#endif /* SWIFT_NONE_PRESSURE_FLOOR_IACT_H */
