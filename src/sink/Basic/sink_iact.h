/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
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
#ifndef SWIFT_BASIC_SINKS_IACT_H
#define SWIFT_BASIC_SINKS_IACT_H

/* Local includes */
#include "gravity.h"
#include "gravity_iact.h"
#include "random.h"
#include "sink_properties.h"

/**
 * @brief Gas particle interactions relevant for sinks, to run in hydro density
 * interaction
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
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {}

/**
 * @brief Gas particle interactions relevant for sinks, to run in hydro density
 * interaction (non symmetric version)
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
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {}

/**
 * @brief Density interaction between sinks and gas (non-symmetric).
 *
 * The gas particle cannot be touched.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing length or cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First particle (sink).
 * @param pj Second particle (gas, not updated).
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sink* si, const struct part* pj, const int with_cosmology,
    const struct cosmology* cosmo, const struct gravity_props* grav_props,
    const struct sink_props* sink_props, const integertime_t ti_current,
    const double time) {

  float wi, wi_dx;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Contribution to the number of neighbours */
  si->num_ngbs++;

  /* Neighbour gas mass */
  const float mj = hydro_get_mass(pj);

  /* Contribution to the BH gas density */
  si->rho_gas += mj * wi;

  /* Contribution to the total neighbour mass */
  si->ngb_mass += mj;

  /* Neighbour's sound speed */
  const float cj = hydro_get_comoving_soundspeed(pj);

  /* Contribution to the smoothed sound speed */
  si->sound_speed_gas += mj * cj * wi;

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we do include a Hubble term) */
  const float dv[3] = {pj->v[0] - si->v[0], pj->v[1] - si->v[1],
                       pj->v[2] - si->v[2]};

  /* Contribution to the smoothed velocity (gas w.r.t. black hole) */
  si->velocity_gas[0] += mj * dv[0] * wi;
  si->velocity_gas[1] += mj * dv[1] * wi;
  si->velocity_gas[2] += mj * dv[2] * wi;

#ifdef SWIFT_SINK_DENSITY_CHECKS
  si->rho_check += mj * wi;
  si->n_check += wi;
  si->N_check_density++;
#endif
}

/**
 * @brief Compute sink-sink swallow interaction (non-symmetric).
 *
 * Note: Energies are computed with physical quantities, not the comoving ones.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing length or cut off radius of particle i.
 * @param hj Comoving smoothing length or cut off radius of particle j.
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
    const float r2, const float dx[3], const float hi, const float hj,
    struct sink* restrict si, struct sink* restrict sj,
    const int with_cosmology, const struct cosmology* cosmo,
    const struct gravity_props* grav_props,
    const struct sink_props* sink_properties, const integertime_t ti_current,
    const double time) {

  /* Simpler version of GEAR as a placeholder. Sinks bound to each other are
   * merged. */

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
  /* TODO: needs updating for MPI safety. We don't have access to foreign gparts
   * here. */
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
 * @param hi Comoving smoothing length or cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sink* restrict si, struct part* restrict pj,
    const int with_cosmology, const struct cosmology* cosmo,
    const struct gravity_props* grav_props,
    const struct sink_props* sink_properties, const integertime_t ti_current,
    const double time) {

  float wi;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Check if the sink needs to be fed. If not, we're done here */
  const float sink_mass_deficit = si->subgrid_mass - si->mass_at_start_of_step;
  if (sink_mass_deficit <= 0) return;

  if (sink_properties->use_nibbling) {

    /* If we do nibbling, things are quite straightforward. We transfer
     * the mass and all associated quantities right here. */

    const float si_mass_orig = si->mass;
    const float pj_mass_orig = hydro_get_mass(pj);

    /* Don't nibble from particles that are too small already */
    if (pj_mass_orig < sink_properties->min_gas_mass_for_nibbling) return;

    /* Next line is equivalent to w_ij * m_j / Sum_j (w_ij * m_j) */
    const float particle_weight = hi_inv_dim * wi * pj_mass_orig / si->rho_gas;
    float nibble_mass = sink_mass_deficit * particle_weight;

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

  } else { /* ends nibbling section, below comes swallowing */

    /* Probability to swallow this particle
     * Recall that in SWIFT the SPH kernel is recovered by computing
     * kernel_eval() and muliplying by (1/h^d) */

    const float prob =
        (si->subgrid_mass - si->mass) * hi_inv_dim * wi / si->rho_gas;

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval_two_IDs(si->id, pj->id, ti_current,
                                                    random_number_sink_swallow);

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
