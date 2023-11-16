/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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

#ifndef SWIFT_GEAR_CHEMISTRY_ADDITIONS_H
#define SWIFT_GEAR_CHEMISTRY_ADDITIONS_H

/**
 * @brief update metal mass fluxes between two interacting particles during
 * hydro_iact_(non)sym(...) calls.
 *
 * @param pi first interacting particle
 * @param pj second interacting particle
 * @param mass_flux the mass flux between these two particles.
 * @param flux_dt the time-step over which the fluxes are exchanged
 * @param mode 0: non-symmetric interaction, update i only. 1: symmetric
 * interaction.
 **/
__attribute__((always_inline)) INLINE static void runner_iact_chemistry_fluxes(
    struct part *restrict pi, struct part *restrict pj, float mass_flux,
    float flux_dt, int mode) {
#ifdef HYDRO_DOES_MASS_FLUX
  error("Unimplemented!");
#endif
}

#endif  // SWIFT_GEAR_CHEMISTRY_ADDITIONS_H
