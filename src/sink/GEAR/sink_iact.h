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

/* Local includes */
#include "gravity.h"
#include "gravity_iact.h"


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
    const float H, const float cut_off_radius) {

  /* In order to prevent the formation of two sink particles at a distance
   * smaller than the sink cutoff radius, we keep only gas particles with
   * the smallest potential. */

  const float r = sqrtf(r2);

  /* if the distance is less than the cut off radius */
  if (r < cut_off_radius) {

    float potential_i = pi->gpart->potential;
    float potential_j = pj->gpart->potential;

    /* prevent the particle with the largest potential to form a sink */
    if (potential_i > potential_j) {
      pi->sink_data.can_form_sink = 0;
      return;
    }

    if (potential_j > potential_i) {
      pj->sink_data.can_form_sink = 0;
      return;
    }
  }
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
    const float H, const float cut_off_radius) {

  /* In order to prevent the formation of two sink particles at a distance
   * smaller than the sink cutoff radius, we keep only gas particles with
   * the smallest potential. */

  const float r = sqrtf(r2);

  if (r < cut_off_radius) {

    float potential_i = pi->gpart->potential;
    float potential_j = pj->gpart->potential;

    /* if the potential is larger
     * prevent the particle to form a sink */
    if (potential_i > potential_j) pi->sink_data.can_form_sink = 0;
  }
}

/**
 * @brief Compute sink-sink swallow interaction (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param rj Comoving cut off radius of particle j.
 * @param si First sink particle.
 * @param sj Second sink particle.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_sink_swallow(const float r2, const float dx[3],
                                      const float ri, const float rj,
                                      struct sink *restrict si,
                                      struct sink *restrict sj) {

  /* See runner_iact_nonsym_bh_bh_swallow.
   * The sink with the smaller mass will be merged onto the one with the
   * larger mass.
   * To avoid rounding issues, we additionally check for IDs if the sink
   * have the exact same mass. */

  /* We should check the relative energy */

  if ((sj->mass < si->mass) || (sj->mass == si->mass && sj->id < si->id)) {

    /* This particle is swallowed by the sink with the largest mass of all the
     * candidates wanting to swallow it (we use IDs to break ties)*/
    if ((sj->merger_data.swallow_mass < si->mass) ||
        (sj->merger_data.swallow_mass == si->mass &&
         sj->merger_data.swallow_id < si->id)) {

      // message("sink %lld wants to swallow sink particle %lld", si->id,
      // sj->id);

      sj->merger_data.swallow_id = si->id;
      sj->merger_data.swallow_mass = si->mass;
    }
  }

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_formation < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_formation[si->num_ngb_formation] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_formation;
#endif
}

/**
 * @brief Compute sink-gas swallow interaction (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(const float r2, const float dx[3],
                                     const float ri, const float hj,
                                     struct sink *restrict si,
                                     struct part *restrict pj,
				     const int with_cosmology,
				     const struct cosmology *cosmo,
				     const struct gravity_props *grav_props,
				     const struct sink_props* sink_properties) {

  /* See see runner_iact_nonsym_bh_gas_swallow.
   * We first check if a gas particle has not been already marked to
   * be swallowed by another sink particle. */

  /* We should check the relative energy */

  // message("sink %lld wants to swallow gas particle %lld", si->id, pj->id);

  const float r = sqrtf(r2);
  const float f_acc_r_acc = sink_properties->f_acc * si->r_cut;

  /* If the gas falls within f_acc*r_acc, it is accreted without further check */
  if (r < f_acc_r_acc) {
    if (pj->sink_data.swallow_id < si->id) {
      pj->sink_data.swallow_id = si->id;
    }
  } else /*Otherwise, we perform other checks */ {
    
    /* Compute the physical relative velocity between the particles */
    const float dv[3] = {(pj->v[0] - si->v[0]) * cosmo->a_inv,
			 (pj->v[1] - si->v[1]) * cosmo->a_inv,
			 (pj->v[2] - si->v[2]) * cosmo->a_inv};

    /* Compute the physical distance between the particles */
    const float dx_physical[3] = {dx[0] * cosmo->a,
				  dx[1] * cosmo->a,
				  dx[2] * cosmo->a};
    const float r_physical = r * cosmo->a;


    /* Momentum check */
    /* Relative momentum of the gas */
    const float specific_angular_momentum_gas[3] = {dx_physical[1] * dv[2] - dx_physical[2] * dv[1],
						    dx_physical[2] * dv[0] - dx_physical[0] * dv[2],
						    dx_physical[0] * dv[1] - dx_physical[1] * dv[0]};
    const float L2_gas = specific_angular_momentum_gas[0]*specific_angular_momentum_gas[0] + specific_angular_momentum_gas[1]*specific_angular_momentum_gas[1] + specific_angular_momentum_gas[2]*specific_angular_momentum_gas[2];

    /* Keplerian angular speed squared */
    const float omega_acc_2 = grav_props->G_Newton*si->mass / (r_physical*r_physical*r_physical);

    /*Keplerian angular momentum squared */
    const float L2_acc = (si->r_cut*si->r_cut*si->r_cut*si->r_cut)*omega_acc_2;

    /* To be accreted, the gas momentum shoulb lower than the keplerian orbit momentum. */
    if (L2_gas > L2_acc) {
      return ;
    }


  }
  



#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_formation < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_formation[si->num_ngb_formation] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_formation;
#endif
}

#endif
