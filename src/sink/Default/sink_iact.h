/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_SINKS_IACT_H
#define SWIFT_DEFAULT_SINKS_IACT_H

/**
 * @brief Compute formation interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_compute_formation(const float r2, const float *dx,
                                           const float hi, const float hj,
                                           struct sink *restrict si,
                                           const struct part *restrict pj,
                                           const float a, const float H) {

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_formation < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_formation[si->num_ngb_formation] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_formation;
#endif
}

/**
 * @brief Compute the sink merger interaction.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param rj Comoving cut off radius of particle j.
 * @param si First sink particle.
 * @param sj Second sink particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 * @param Which particle should be removed?
 * Possible value: (sink_merger_remove_none/first/second)
 */
__attribute__((always_inline)) INLINE static int runner_iact_sym_sinks_merger(
    const float r2, const float *dx, const float hi, const float hj,
    struct sink *restrict si, struct sink *restrict sj, const float a,
    const float H) {

  return sink_merger_remove_none;

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_merger < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_merger[si->num_ngb_merger] = sj->id;

  /* Update ngb counters */
  ++si->num_ngb_merger;
#endif
}

#endif
