/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_STARS_IACT_H
#define SWIFT_NONE_STARS_IACT_H

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * Empty imlementation
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_density(const float r2, const float *dx,
                                 const float hi, const float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, const float a,
                                 const float H) {}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 *
 * Empty implementation.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_feedback(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  struct spart *restrict si,
                                  struct part *restrict pj, const float a,
                                  const float H) {}

#endif /* SWIFT_NONE_STARS_IACT_H */
