/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_SINKS_IACT_H
#define SWIFT_GEAR_SINKS_IACT_H



/**
 * @brief do sink computation after the runner_iact_density (symmetric
 * version)
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
__attribute__((always_inline)) INLINE static void runner_iact_sink(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const struct sink_props *sink_props) {

}


/**
 * @brief do sink computation after the runner_iact_density (non symmetric
 * version)
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_sink(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H, const struct sink_props *sink_props) {


  /* In order to prevent the formation of two sink particles at a distance
   * smaller than the sink cutoff radius, we keep only gas particles with
   * the smallest potential. */

  const float r = sqrtf(r2);
  
  if (r < sink_props->cut_off_radius) {
    
    float potential_i = pi->gpart->potential;
    float potential_j = pj->gpart->potential;
        
    /* if the potential is larger
     * prevent the particle to form a sink */
    if(potential_i > potential_j)
      pi->sink_data.can_form_sink = 0;
  } 
  
}







/**
 * @brief Compute sink swallow interaction between two particles (non-symmetric).
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
runner_iact_nonsym_sinks_swallow(const float r2, const float *dx,
                                           const float hi, const float hj,
                                           struct sink *restrict si,
                                           struct part *restrict pj,
                                           const float a, const float H) {

//message("%lld",pj->sink_data.swallow_id);


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

  message(">>");
  if (sqrt(r2) < si->r_cut || sqrt(r2) < sj->r_cut){
      if (si->id > sj->id){
          sj->mass = sj->mass + si->mass;       // should also add momentum and metals
          return sink_merger_remove_first;
        }  
      else {
          si->mass = si->mass + sj->mass;       // should also add momentum and metals
          return sink_merger_remove_second;  
        }
    }
  
  return sink_merger_remove_none;

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_merger < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_merger[si->num_ngb_merger] = sj->id;

  /* Update ngb counters */
  ++si->num_ngb_merger;
#endif
}

/**
 * @brief Accretion interaction between two particles (non-symmetric).
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
runner_iact_nonsym_sinks_accretion(const float r2, const float *dx,
                                   const float hi, const float hj,
                                   struct sink *restrict si,
                                   const struct part *restrict pj,
                                   const float a, const float H) {

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_accretion < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_accretion[si->num_ngb_accretion] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_accretion;
#endif
}

#endif
