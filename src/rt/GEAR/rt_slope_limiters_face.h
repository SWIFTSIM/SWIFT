/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H
#define SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H

/**
 * @file src/rt/GEAR/rt_slope_limiters_face.h
 * @brief File containing routines concerning the face slope
 * limiter for the GEAR RT scheme. (= second slope limiting
 * step done during actual particle interactions.) */

/**
 * @brief Slope limit the slopes at the interface between two particles
 *
 * @param Qi RT quantities (photon energies and flux) of particle i.
 * @param Qj RT quantities (photon energies and flux) of particle j.
 * @param dQi Difference between the RT quantities of particle i at the
 * position of particle i and at the interface position.
 * @param dQj Difference between the RT quantities of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_face(
    const float Qi[4], const float Qj[4], float dQi[4], float dQj[4],
    const float xij_i[3], const float xij_j[3], const float r) {}

#endif /* SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H */
