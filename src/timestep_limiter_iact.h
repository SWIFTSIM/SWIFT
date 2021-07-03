/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TIMESTEP_LIMITER_IACT_H
#define SWIFT_TIMESTEP_LIMITER_IACT_H

/**
 * @brief Force interaction between two particles.
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
__attribute__((always_inline)) INLINE static void runner_iact_timebin(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  /* Update the minimal time-bin */
  if (pj->time_bin > 0)
    pi->limiter_data.min_ngb_time_bin =
        min(pi->limiter_data.min_ngb_time_bin, pj->time_bin);

  if (pi->time_bin > 0)
    pj->limiter_data.min_ngb_time_bin =
        min(pj->limiter_data.min_ngb_time_bin, pi->time_bin);
}

/**
 * @brief Timebin interaction between two particles (non-symmetric).
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_timebin(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  /* Update the minimal time-bin */
  if (pj->time_bin > 0)
    pi->limiter_data.min_ngb_time_bin =
        min(pi->limiter_data.min_ngb_time_bin, pj->time_bin);
}

/**
 * @brief Timestep limiter loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_limiter(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  /* Nothing to do here if both particles are active */
}

/**
 * @brief Timestep limiter loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_limiter(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  /* Wake up the neighbour? */
  if (pj->time_bin > pi->time_bin + time_bin_neighbour_max_delta_bin) {

    /* Store the smallest time bin that woke up this particle */
    accumulate_max_c(&pj->limiter_data.wakeup, -pi->time_bin);
  }
}

#endif /* SWIFT_TIMESTEP_LIMITER_IACT_H */
