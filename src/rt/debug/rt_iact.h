/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IACT_DEBUG_H
#define SWIFT_RT_IACT_DEBUG_H

/**
 * @file src/rt/debug/rt_iact.h
 * @brief Main header file for the debug radiative transfer scheme particle
 * interactions.
 */

/**
 * @brief Injection step interaction between star and hydro particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si Star particle.
 * @param xpj Hydro particle extra data.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_inject(
    const float r2, float* dx, const float hi, const float hj,
    struct spart* restrict si, struct xpart* restrict xpj) {

  struct rt_spart_data* restrict sd = &(si->rt_data);
  struct rt_xpart_data* restrict pd = &(xpj->rt_data);

  sd->iact_hydro += 1;
  sd->calls_tot += 1;
  sd->calls_per_step += 1;

  pd->iact_stars += 1;
  pd->calls_tot += 1;
  sd->calls_per_step += 1;
}

#endif /* SWIFT_RT_IACT_DEBUG_H */
