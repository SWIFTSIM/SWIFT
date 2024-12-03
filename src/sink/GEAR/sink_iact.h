/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#include "sink.h"
#include "sink_getters.h"
#include "sink_properties.h"

/**
 * @brief do sink computation after the runner_iact_density (symmetric
 * version)
 *
 * In GEAR: This function deactivates the sink formation ability of #part not
 * at a potential minimum.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param cut_off_radius Sink cut off radius.
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
    float potential_i = pi->sink_data.potential;
    float potential_j = pj->sink_data.potential;

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
 * In GEAR: This function deactivates the sink formation ability of #part not
 * at a potential minimum.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param cut_off_radius Sink cut off radius.
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
    float potential_i = pi->sink_data.potential;
    float potential_j = pj->sink_data.potential;

    /* if the potential is larger
     * prevent the particle to form a sink */
    if (potential_i > potential_j) pi->sink_data.can_form_sink = 0;
  }
}

/* @brief Density interaction between two particles (non-symmetric).
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
    const float r2, const float dx[3], const float ri, const float hj,
    struct sink *si, const struct part *pj, const int with_cosmology,
    const struct cosmology *cosmo, const struct gravity_props *grav_props,
    const struct sink_props *sink_props, const integertime_t ti_current,
    const double time) {

  /* Contribution to the number of neighbours in cutoff radius */
  si->num_ngbs++;

  float wi, wi_dx;

  /* Compute the kernel function */
  const float r = sqrtf(r2);
  const float hi = ri / kernel_gamma;
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Neighbour gas mass */
  const float mj = hydro_get_mass(pj);

  /* Minimum smoothing length accros the neighbours */
  /* AND the sink smoothing length */
  si->to_collect.minimal_h_gas = min(hj, si->to_collect.minimal_h_gas);

  /* Contribution to the BH gas density */
  si->to_collect.rho_gas += mj * wi;

  /* Contribution to the smoothed sound speed */
  si->to_collect.sound_speed_gas += mj * wi * hydro_get_comoving_soundspeed(pj);

  /* Neighbour's (drifted) velocity in the frame of the sink
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the sink) */
  const float dv[3] = {pj->v[0] - si->v[0], pj->v[1] - si->v[1],
                       pj->v[2] - si->v[2]};

  /* Contribution to the smoothed velocity (gas w.r.t. black hole) */
  si->to_collect.velocity_gas[0] += mj * dv[0] * wi;
  si->to_collect.velocity_gas[1] += mj * dv[1] * wi;
  si->to_collect.velocity_gas[2] += mj * dv[2] * wi;
}

/**
 * @brief  Update the properties of a sink particles from its sink neighbours.
 *
 * Warning: No symmetric interaction for timesteps since we are in
 * runner_iact_nonsym_sinks_sink_swallow() --> pay attention to not break MPI.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param rj Comoving cut off radius of particle j.
 * @param si First sink particle.
 * @param sj Second sink particle.
 */
__attribute__((always_inline)) INLINE static void
sink_collect_properties_from_sink(const float r2, const float dx[3],
                                  const float ri, const float rj,
                                  struct sink *restrict si,
                                  struct sink *restrict sj,
                                  const struct gravity_props *grav_props) {

  /* Neighbour's (drifted) velocity in the frame of the sink i
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the sink) */
  const float dv[3] = {sj->v[0] - si->v[0], sj->v[1] - si->v[1],
                       sj->v[2] - si->v[2]};
  const float dv_norm = sqrtf(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);

  /* Get the gravitional softening */
  const float eps = gravity_get_softening(si->gpart, grav_props);
  const float eps2 = eps * eps;
  const float eps_inv = 1.f / eps;
  const float eps_inv3 = eps_inv * eps_inv * eps_inv;

  /* Compute the kernel potential and force with mass = 1.0. We multiply by
     the mass below if needed. */
  float dphi_dr, pot;
  runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, 1.0, &dphi_dr, &pot);

  /* From Grudic et al. (2021) eq 6, we replace the plummer functionnal form
     sqrt(r^2 + eps^2) by the kernel 1.0/|phi(r,H=3*eps)| */
  const float t_c = 1.0 / (fabsf(pot) * dv_norm);
  si->to_collect.minimal_sink_t_c = min(t_c, si->to_collect.minimal_sink_t_c);

  /* From Grudic et al. (2021) eq 7, we replace the plummer functionnal form
     (r^2 + eps^2)^{3/2} by the kernel |(d phi(r,H=3*eps)/ dr)^{-1}| */
  const float denominator = grav_props->G_Newton * (si->mass + sj->mass);
  const float numerator = 1.0 / fabsf(dphi_dr);
  const float t_dyn = sqrt(numerator / denominator);
  si->to_collect.minimal_sink_t_dyn =
      min(t_dyn, si->to_collect.minimal_sink_t_dyn);
}

/**
 * @brief Compute sink-sink swallow interaction (non-symmetric).
 *
 * Note: Energies are computed with physical quantities, not the comoving ones.
 *
 * MPI note: This functions invokes the gpart. Hence, it must be performed only
 * on the local node (similarly to runner_iact_nonsym_bh_bh_repos()).
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
    const struct sink_props *sink_properties, const integertime_t ti_current,
    const double time) {

  const float r = sqrtf(r2);
  const float f_acc_r_acc_i = sink_properties->f_acc * ri;

  /* Determine if the sink is dead, i.e. if its age is bigger than the
     age_threshold_unlimited */
  const int si_age = sink_get_sink_age(si, with_cosmology, cosmo, time);
  char si_is_dead = si_age > sink_properties->age_threshold_unlimited;

  const int sj_age = sink_get_sink_age(sj, with_cosmology, cosmo, time);
  char sj_is_dead = sj_age > sink_properties->age_threshold_unlimited;

  /* Collect the properties for 2-body interactions if one is sink is alive. If
     they are both dead, we do not want to restrict the timesteps for 2-body
     encounters since they won't merge. */
  if (!si_is_dead || !sj_is_dead) {
    sink_collect_properties_from_sink(r2, dx, ri, rj, si, sj, grav_props);
  }

  /* If si is dead, do not swallow sj. However, sj can swallow si if it alive.
   */
  if (si_is_dead) {
    return;
  }

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

    /* Compute the physical relative velocity between the particles */
    const float dv_physical[3] = {v_plus_H_flow[0] * cosmo->a_inv,
                                  v_plus_H_flow[1] * cosmo->a_inv,
                                  v_plus_H_flow[2] * cosmo->a_inv};

    const float dv_physical_squared = dv_physical[0] * dv_physical[0] +
                                      dv_physical[1] * dv_physical[1] +
                                      dv_physical[2] * dv_physical[2];

    /* Momentum check------------------------------------------------------- */
    float L2_j = 0.0;      /* Relative momentum of the sink j */
    float L2_kepler = 0.0; /* Keplerian angular momentum squared */
    sink_compute_angular_momenta_criterion(dx, v_plus_H_flow, r, si->r_cut,
                                           si->mass, cosmo, grav_props,
                                           &L2_kepler, &L2_j);

    /* To be accreted, the sink momentum should lower than the keplerian orbit
     * momentum. */
    if (L2_j > L2_kepler) {
      return;
    }

    /* Binding energy check------------------------------------------------- */
    /* Kinetic energy per unit mass of the sink */
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

    /* Swallowed mass threshold--------------------------------------------- */
    si->to_collect.mass_eligible_swallow += sj->mass;

    /* Maximal mass that can be swallowed within a single timestep */
    const float mass_swallow_limit = sink_properties->n_IMF * si->mass_IMF;

    /* If the mass exceeds the threshold, do not swallow. Make sure you can at
       least swallow a particle to avoid running into the problem of never being
       able to spawn a star.
       If n_IMF <= 0, then disable this criterion */
    if (sink_properties->n_IMF > 0 &&
        si->to_collect.mass_swallowed >= mass_swallow_limit &&
        si->to_collect.mass_eligible_swallow != 0) {
      return;
    }

    /* Increment the swallowd mass */
    si->to_collect.mass_swallowed += sj->mass;

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
 * MPI note: This functions invokes the gpart. Hence, it must be performed only
 * on the local node (similarly to runner_iact_nonsym_bh_gas_repos()).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 * @param with_cosmology if we run with cosmology.
 * @param cosmo The cosmological parameters and properties.
 * @param grav_props The gravity scheme parameters and properties.
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(
    const float r2, const float dx[3], const float ri, const float hj,
    struct sink *restrict si, struct part *restrict pj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct sink_props *sink_properties, const integertime_t ti_current,
    const double time) {

  const float r = sqrtf(r2);
  const float f_acc_r_acc = sink_properties->f_acc * ri;

  /* Determine if the sink is dead, i.e. if its age is bigger than the
     age_threshold_unlimited */
  const int sink_age = sink_get_sink_age(si, with_cosmology, cosmo, time);
  char is_dead = sink_age > sink_properties->age_threshold_unlimited;

  /* If si is dead, do not swallow pj. */
  if (is_dead) {
    return;
  }

  /* If the gas falls within f_acc*r_acc, it is accreted without further check
   */
  if (r < f_acc_r_acc) {
    warning("Gas %lld within sink %lld inner accretion radius", pj->id, si->id);
    /* Check if a gas particle has not been already marked to be swallowed by
       another sink particle. */
    if (pj->sink_data.swallow_id < si->id) {
      pj->sink_data.swallow_id = si->id;
    }

    /* f_acc*r_acc <= r <= r_acc, we perform other checks */
  } else if ((r >= f_acc_r_acc) && (r < ri)) {

    /* Relative velocity between the sinks */
    const float dv[3] = {pj->v[0] - si->v[0], pj->v[1] - si->v[1],
                         pj->v[2] - si->v[2]};

    const float a = cosmo->a;
    const float H = cosmo->H;
    const float a2H = a * a * H;

    /* Calculate the velocity with the Hubble flow */
    const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                    a2H * dx[2] + dv[2]};

    /* Compute the physical relative velocity between the particles */
    const float dv_physical[3] = {v_plus_H_flow[0] * cosmo->a_inv,
                                  v_plus_H_flow[1] * cosmo->a_inv,
                                  v_plus_H_flow[2] * cosmo->a_inv};

    const float dv_physical_squared = dv_physical[0] * dv_physical[0] +
                                      dv_physical[1] * dv_physical[1] +
                                      dv_physical[2] * dv_physical[2];

    /* Momentum check------------------------------------------------------- */
    float L2_gas_j = 0.0;  /* Relative momentum of the gas */
    float L2_kepler = 0.0; /* Keplerian angular momentum squared */
    sink_compute_angular_momenta_criterion(dx, v_plus_H_flow, r, si->r_cut,
                                           si->mass, cosmo, grav_props,
                                           &L2_kepler, &L2_gas_j);

    /* To be accreted, the gas momentum should lower than the keplerian orbit
     * momentum. */
    if (L2_gas_j > L2_kepler) {
      return;
    }

    /* Energy check--------------------------------------------------------- */
    /* Kinetic energy per unit mass of the gas */
    float E_kin_relative_gas = 0.5f * dv_physical_squared;

    /* Compute the Newtonian or softened potential the sink exherts onto the
       gas particle */
    const float eps = gravity_get_softening(si->gpart, grav_props);
    const float eps2 = eps * eps;
    const float eps_inv = 1.f / eps;
    const float eps_inv3 = eps_inv * eps_inv * eps_inv;
    const float sink_mass = si->mass;
    float dummy, pot_ij;
    runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, sink_mass, &dummy,
                             &pot_ij);

    /* Compute the physical potential energy per unit mass  that the sink
       exerts in the gas :
                       E_pot_phys = G*pot_grav*a^(-1) + c(a).
       The normalization is c(a) = 0. */
    const float E_pot_gas = grav_props->G_Newton * pot_ij * cosmo->a_inv;

    /* Update: Add thermal energy per unit mass  to avoid the sink to swallow
       hot gas regions */
    const float E_therm = hydro_get_drifted_physical_internal_energy(pj, cosmo);

    /* Energy per unit mass of the pair sink-gas */
    const float E_mec_sink_part = E_kin_relative_gas + E_pot_gas + E_therm;

    /* To be accreted, the gas must be gravitationally bound to the sink. */
    if (E_mec_sink_part >= 0) return;

    /* To be accreted, the gas smoothing length must be smaller than the sink
       accretion radius. This is similar to AMR codes requesting the maximum
       refinement level close to the sink. */
    if (sink_properties->sink_formation_smoothing_length_criterion &&
        (pj->h * kernel_gamma >= si->r_cut))
      return;

    /* Most bound pair check------------------------------------------------ */
    /* The pair gas-sink must be the most bound among all sinks */
    if (E_mec_sink_part >= pj->sink_data.E_mec_bound) {
      return;
    }

    /* Swallowed mass threshold--------------------------------------------- */
    si->to_collect.mass_eligible_swallow += hydro_get_mass(pj);

    /* Maximal mass that can be swallowed within a single timestep */
    const float mass_swallow_limit = sink_properties->n_IMF * si->mass_IMF;

    /* If the mass exceeds the threshold, do not swallow. Make sure you can at
       least swallow a particle to avoid running into the problem of never being
       able to spawn a star.
       If n_IMF <= 0, then disable this criterion */
    if (sink_properties->n_IMF > 0 &&
        si->to_collect.mass_swallowed >= mass_swallow_limit &&
        si->to_collect.mass_eligible_swallow != 0) {
      return;
    }

    /* Increment the swallowd mass */
    si->to_collect.mass_swallowed += hydro_get_mass(pj);

    /* --------------------------------------------------------------------- */
    /* Since this pair gas-sink is the most bound, keep track of the
       E_mec_bound and set the swallow_id accordingly */
    pj->sink_data.E_mec_bound = E_mec_sink_part;
    pj->sink_data.swallow_id = si->id;
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
