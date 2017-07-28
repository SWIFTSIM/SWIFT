
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_HYDRO_FLUX_LIMITERS_H
#define SWIFT_HYDRO_FLUX_LIMITERS_H

#ifdef GIZMO_FLUX_LIMITER

#define HYDRO_FLUX_LIMITER_IMPLEMENTATION "GIZMO flux limiter"

/**
 * @brief Limit the flux between two particles.
 *
 * @param flux Unlimited flux between the particles.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_flux_limiters_apply(
    float* flux, struct part* pi, struct part* pj) {

  float flux_limit_factor = 1.;
  const float timefac = max(pi->force.dt, pj->force.dt);
  const float areafac = max(pi->geometry.Atot, pj->geometry.Atot);
  const float totfac = timefac * areafac;
  if (flux[0] * totfac > pi->conserved.mass) {
    flux_limit_factor = pi->conserved.mass / (flux[0] * totfac);
  }
  if (flux[0] * totfac > pj->conserved.mass) {
    flux_limit_factor =
        min(pj->conserved.mass / (flux[0] * totfac), flux_limit_factor);
  }
  if (flux[4] * totfac > pi->conserved.energy) {
    flux_limit_factor =
        min(pi->conserved.energy / (flux[4] * totfac), flux_limit_factor);
  }
  if (flux[4] * totfac > pj->conserved.energy) {
    flux_limit_factor =
        min(pj->conserved.energy / (flux[4] * totfac), flux_limit_factor);
  }

  flux[0] *= flux_limit_factor;
  flux[1] *= flux_limit_factor;
  flux[2] *= flux_limit_factor;
  flux[3] *= flux_limit_factor;
  flux[4] *= flux_limit_factor;
}

#else

#define HYDRO_FLUX_LIMITER_IMPLEMENTATION "No flux limiter"

/**
 * @brief Limit the flux between two particles.
 *
 * @param flux Unlimited flux between the particles.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_flux_limiters_apply(
    float* flux, struct part* pi, struct part* pj) {}

#endif

#endif  // SWIFT_HYDRO_FLUX_LIMITERS_H
