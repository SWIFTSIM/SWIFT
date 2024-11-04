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
#ifndef SWIFT_GEARBONDIHOYLE_SINKS_IACT_H
#define SWIFT_GEARBONDIHOYLE_SINKS_IACT_H

/* Local includes */
#include "gravity.h"
#include "gravity_iact.h"
#include "sink_properties.h"

/**
 * @brief do sink computation after the runner_iact_density (symmetric
 * version)
 *
 * Note: This functions breaks MPI.
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

    /*
     * NOTE: Those lines break MPI
     */
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
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First particle (sink).
 * @param pj Second particle (gas, not updated).
 * @param xpj The extended data of the second particle (not updated).
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_density(
    const float r2, const float dx[3], const float ri, const float hj,
    struct sink *si, const struct part *pj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct sink_props *sink_props,
    const integertime_t ti_current, const double time) {

  const float r = sqrtf(r2);

  if (r < ri) {

    /* Neighbour gas mass */
    const float mj = hydro_get_mass(pj);

    /* Contribution to the number of neighbours & mass in cutoff radius */
    /* We'll get the gas density from this in the ghost task */
    si->num_ngbs++;
    si->ngb_mass += mj;

    /* Neighbour's sound speed */
    /* COLIBRE BHs do something more complex than this, keep it simple for now */
    float cj = hydro_get_comoving_soundspeed(pj);

    /* Contribution to the smoothed sound speed */
    si->sound_speed_gas += mj * cj;

    /* Neighbour's (drifted) velocity in the frame of the black hole
    * (we do include a Hubble term) */
    const float dv[3] = {pj->v[0] - si->v[0], pj->v[1] - si->v[1],
                        pj->v[2] - si->v[2]};

    /* Contribution to the smoothed velocity (gas w.r.t. sink) */
    si->velocity_gas[0] += mj * dv[0];
    si->velocity_gas[1] += mj * dv[1];
    si->velocity_gas[2] += mj * dv[2];
  }
}

/**
 * @brief Compute sink-sink swallow interaction (non-symmetric).
 *
 * Note: Energies are computed with physical quantities, not the comoving ones.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param rj Comoving cut off radius of particle j.
 * @param si First sink particle.
 * @param sj Second sink particle.
 * @param with_cosmology if we run with cosmology.
 * @param cosmo The cosmological parameters and properties.
 * @param grav_props The gravity scheme parameters and properties.
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_sink_swallow(
    const float r2, const float dx[3], const float ri, const float rj,
    struct sink *restrict si, struct sink *restrict sj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct sink_props *sink_properties,
    const integertime_t ti_current, const double time) {

  const float r = sqrtf(r2);
  const float f_acc_r_acc_i = sink_properties->f_acc * ri;

  /* If the sink j falls within f_acc*r_acc of sink i, then the
     lightest is accreted on the most massive without further check.
     Note that this is a non-symmetric interaction. So, we do not need to check
     for the f_acc*r_acc_j case here. */
  if (r < f_acc_r_acc_i) {
    /* The sink with the smaller mass will be merged onto the one with the
     * larger mass.
     * To avoid rounding issues, we additionally check for IDs if the sink
     * have the exact same mass. */
    if ((sj->mass < si->mass) || (sj->mass == si->mass && sj->id < si->id)) {
      /* This particle is swallowed by the sink with the largest mass of all the
       * candidates wanting to swallow it (we use IDs to break ties)*/
      if ((sj->merger_data.swallow_mass < si->mass) ||
          (sj->merger_data.swallow_mass == si->mass &&
           sj->merger_data.swallow_id < si->id)) {
        sj->merger_data.swallow_id = si->id;
        sj->merger_data.swallow_mass = si->mass;
      }
    }
  } else {

    /* Relative velocity between the sinks */
    const float dv[3] = {sj->v[0] - si->v[0], sj->v[1] - si->v[1],
                         sj->v[2] - si->v[2]};

    const float a = cosmo->a;
    const float H = cosmo->H;
    const float a2H = a * a * H;

    /* Calculate the velocity with the Hubble flow */
    const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                    a2H * dx[2] + dv[2]};

    /* Binding energy check */
    /* Compute the physical relative velocity between the particles */
    const float dv_physical[3] = {v_plus_H_flow[0] * cosmo->a_inv,
                                  v_plus_H_flow[1] * cosmo->a_inv,
                                  v_plus_H_flow[2] * cosmo->a_inv};

    const float dv_physical_squared = dv_physical[0] * dv_physical[0] +
                                      dv_physical[1] * dv_physical[1] +
                                      dv_physical[2] * dv_physical[2];

    /* Kinetic energy per unit mass of the gas */
    const float E_kin_rel = 0.5f * dv_physical_squared;

    /* Compute the Newtonian or softened potential the sink exherts onto the
       gas particle */
    const float eps = gravity_get_softening(si->gpart, grav_props);
    const float eps2 = eps * eps;
    const float eps_inv = 1.f / eps;
    const float eps_inv3 = eps_inv * eps_inv * eps_inv;
    const float si_mass = si->mass;
    const float sj_mass = sj->mass;

    float dummy, pot_ij, pot_ji;
    runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, si_mass, &dummy,
                             &pot_ij);
    runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, sj_mass, &dummy,
                             &pot_ji);

    /* Compute the physical potential energies per unit mass :
                           E_pot_phys = G*pot_grav*a^(-1) + c(a).
       The normalization is c(a) = 0. */
    const float E_pot_ij = grav_props->G_Newton * pot_ij * cosmo->a_inv;
    const float E_pot_ji = grav_props->G_Newton * pot_ji * cosmo->a_inv;

    /* Mechanical energy per unit mass of the pair i-j and j-i */
    const float E_mec_si = E_kin_rel + E_pot_ij;
    const float E_mec_sj = E_kin_rel + E_pot_ji;

    /* Now, check if one is bound to the other */
    if ((E_mec_si > 0) || (E_mec_sj > 0)) {
      return;
    }

    /* The sink with the smaller mass will be merged onto the one with the
     * larger mass.
     * To avoid rounding issues, we additionally check for IDs if the sink
     * have the exact same mass. */
    if ((sj->mass < si->mass) || (sj->mass == si->mass && sj->id < si->id)) {
      /* This particle is swallowed by the sink with the largest mass of all the
       * candidates wanting to swallow it (we use IDs to break ties)*/
      if ((sj->merger_data.swallow_mass < si->mass) ||
          (sj->merger_data.swallow_mass == si->mass &&
           sj->merger_data.swallow_id < si->id)) {
        sj->merger_data.swallow_id = si->id;
        sj->merger_data.swallow_mass = si->mass;
      }
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
 * Note: Energies are computed with physical quantities, not the comoving ones.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(const float r2, const float dx[3],
                                     const float ri, const float hj,
                                     struct sink *restrict si,
                                     struct part *restrict pj,
                                     const int with_cosmology,
                                     const struct cosmology *cosmo,
                                     const struct gravity_props *grav_props,
                                     const struct sink_props *sink_properties,
                                     const integertime_t ti_current, const double time) {


  /* Get r. */
  const float r = sqrtf(r2);

  /* Check if within cutoff radius. If not, we're done. */
  if (r > ri) return;

  /* Check if the sink needs to be fed. If not, we're done here */
  if (si->mass_to_accrete <= 0) return;

  if (sink_properties->use_nibbling) {

    /* If we do nibbling, things are quite straightforward. We transfer
     * the mass and all associated quantities right here. */

    const float si_mass_orig = si->mass;
    const float pj_mass_orig = hydro_get_mass(pj);

    /* Don't nibble from particles that are too small already */
    if (pj_mass_orig < sink_properties->min_gas_mass_for_nibbling) return;

    /* Sinks don't have a kernel, nibble equally from all particles in cutoff radius*/
    float nibble_mass = si->mass_to_accrete / si->num_ngbs;

    /* Need to check whether nibbling would push gas mass below minimum
     * allowed mass */
    float new_gas_mass = pj_mass_orig - nibble_mass;
    if (new_gas_mass < sink_properties->min_gas_mass_for_nibbling) {
      new_gas_mass = sink_properties->min_gas_mass_for_nibbling;
      nibble_mass = pj_mass_orig - sink_properties->min_gas_mass_for_nibbling;
    }

    /* Transfer (dynamical) mass from the gas particle to the sink */
    si->mass += nibble_mass;
    hydro_set_mass(pj, new_gas_mass);

    /* Track total mass accreted */
    si->total_accreted_gas_mass += nibble_mass;

    /* Add the angular momentum of the accreted gas to the sink total.
     * Note no change to gas here. The cosmological conversion factors for
     * velocity (a^-1) and distance (a) cancel out, so the angular momentum
     * is already in physical units. */
    const float dv[3] = {si->v[0] - pj->v[0], si->v[1] - pj->v[1],
                         si->v[2] - pj->v[2]};
    si->swallowed_angular_momentum[0] +=
        nibble_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
    si->swallowed_angular_momentum[1] +=
        nibble_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
    si->swallowed_angular_momentum[2] +=
        nibble_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

    /* Update the BH momentum and velocity. Again, no change to gas here. */
    const float si_mom[3] = {si_mass_orig * si->v[0] + nibble_mass * pj->v[0],
                             si_mass_orig * si->v[1] + nibble_mass * pj->v[1],
                             si_mass_orig * si->v[2] + nibble_mass * pj->v[2]};

    si->v[0] = si_mom[0] / si->mass;
    si->v[1] = si_mom[1] / si->mass;
    si->v[2] = si_mom[2] / si->mass;

    const float nibbled_fraction = nibble_mass / pj_mass_orig;

    /* Update the BH and also gas metal & dust masses */
    struct chemistry_sink_data *si_chem = &si->chemistry_data;
    struct chemistry_part_data *pj_chem = &pj->chemistry_data;

    /* This should ultimately be in a function called chemistry_transfer_part_to_sink in the chemistry module */
    for (int k = 0; k < GEAR_CHEMISTRY_ELEMENT_COUNT; k++) {
      double mk = si_chem->metal_mass_fraction[k] * si_mass_orig +
                  pj_chem->smoothed_metal_mass_fraction[k] * pj_mass_orig;

      si_chem->metal_mass_fraction[k] = mk / si->mass;
    }

    /* We will need this for the GCs */
    // struct dust_bpart_data *si_dust = &si->dust_data;
    // struct dust_part_data *pj_dust = &pj->dust_data;
    // dust_transfer_part_to_bpart(si_dust, pj_dust, nibble_mass,
    //                             nibbled_fraction);

  } else { /* ends nibbling section, below comes swallowing */

    /* Probability to swallow this particle
     * Recall that in SWIFT the SPH kernel is recovered by computing
     * kernel_eval() and muliplying by (1/h^d) */
    
    /* Sinks don't currently have a kernel - all particles in cutoff radius have equal probability */
    const float prob = si->mass_to_accrete / si->ngb_mass;
    
    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(si->id + pj->id, ti_current,
                                            random_number_BH_swallow);
    
    /* Are we lucky? */
    if (rand < prob) {

      /* This particle is swallowed by the BH with the largest ID of all the
       * candidates wanting to swallow it */
      if (pj->sink_data.swallow_id < si->id) {

        message("Sink %lld wants to swallow gas particle %lld", si->id, pj->id);

        pj->sink_data.swallow_id = si->id;

      } else {

        message(
            "Sink %lld wants to swallow gas particle %lld BUT CANNOT (old "
            "swallow id=%lld)",
            si->id, pj->id, pj->sink_data.swallow_id);
      }
    }
  } /* ends section for swallowing */
}

#endif
