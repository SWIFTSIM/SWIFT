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
#include "hydro.h"

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
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_sink_swallow(const float r2, const float dx[3],
                                      const float ri, const float rj,
                                      struct sink *restrict si,
                                      struct sink *restrict sj,
                                      const int with_cosmology,
                                      const struct cosmology *cosmo,
                                      const struct gravity_props *grav_props,
				      const struct sink_props *sink_properties) {

  const float r = sqrtf(r2);
  const float f_acc_r_acc_i = sink_properties->f_acc * ri;
  const float f_acc_r_acc_j = sink_properties->f_acc * rj;

  /* If the sink i falls within f_acc*r_acc_i, it is accreted without further
     check */
  if (r < f_acc_r_acc_i) {
    /* Check if a sink particle has not been already marked to be swallowed by
       another sink particle. */
    if (sj->merger_data.swallow_id < si->id) {
      sj->merger_data.swallow_id = si->id;
    }

  /* If the sink j falls within f_acc*r_acc_j, it is accreted without further
     check */
  } else if (r < f_acc_r_acc_j) {

    /* Check if a sink particle has not been already marked to be swallowed by
       another sink particle. */
    if (si->merger_data.swallow_id < sj->id) {
      si->merger_data.swallow_id = sj->id;
    }

  } else {

    /* Relative velocity between th sinks */
    const float dv[3] = {sj->v[0] - si->v[0],  sj->v[1] - si->v[1],
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

    const float dv_physical_squared =  dv_physical[0] * dv_physical[0]
      + dv_physical[1] * dv_physical[1]
      + dv_physical[2] * dv_physical[2];

    /* Kinetic energy of the gas */
    const float E_kin_rel = 0.5f * dv_physical_squared;

    /* Compute the Newtonian or truncated potential the sink exherts onto the
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

    const double Omega_r = cosmo->Omega_r + cosmo->Omega_nu;
    const double Omega_m = cosmo->Omega_cdm + cosmo->Omega_b;
    const double Omega_l = cosmo->Omega_lambda;
    const double w0 = cosmo->w_0;
    const double wa = cosmo->w_a;
    const double a_inv = cosmo->a_inv;

    const double w_DE =  w0 + wa * (1. - a) ; //cosmology_dark_energy_EoS(a, w0, wa);
    const double w_tilde = (a - 1.) * wa - (1. + w0 + wa) * log(a);
    const double density_sum = Omega_m * a_inv * a_inv * a_inv
      + 2.0*Omega_r* a_inv * a_inv * a_inv * a_inv
      + Omega_l * exp(3. * w_tilde) * (1 + w_DE);
    //w_tilde(a, w0, wa)
    const double a_dot_dot = - H*H/2.0 * density_sum;

    /* Compute the physical potential energies :
       E_pot_phys = G*pot_grav*a^(-1) + c(a). */
    /* The normalization is c(a) = -a_dot*a*r^2/2.0. */
    const double constant = - a_dot_dot*a*r2/2.0;
    const float E_pot_ij = grav_props->G_Newton * pot_ij * cosmo->a_inv + constant;
    const float E_pot_ji = grav_props->G_Newton * pot_ji * cosmo->a_inv + constant;

    /* Mechanical energy of the pair i-j and j-i */
    const float E_mec_si = E_kin_rel + E_pot_ij;
    const float E_mec_sj = E_kin_rel + E_pot_ji;

    /* Now, check if one is bound to the other */
    if ((E_mec_si > 0) && (E_mec_sj > 0)) {
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
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(const float r2, const float dx[3],
                                     const float ri, const float hj,
                                     struct sink *restrict si,
                                     struct part *restrict pj,
                                     const int with_cosmology,
                                     const struct cosmology *cosmo,
                                     const struct gravity_props *grav_props,
                                     const struct sink_props *sink_properties) {

  const float r = sqrtf(r2);
  const float f_acc_r_acc = sink_properties->f_acc * ri;

  /* If the gas falls within f_acc*r_acc, it is accreted without further check
   */
  if (r < f_acc_r_acc) {
    /* Check if a gas particle has not been already marked to be swallowed by
       another sink particle. */
    if (pj->sink_data.swallow_id < si->id) {
      pj->sink_data.swallow_id = si->id;
    }

  /* f_acc*r_acc <= r <= r_acc, we perform other checks */
  } else if ((r >= f_acc_r_acc) && (r < ri))  {

     /* Relative velocity between th sinks */
    const float dv[3] = {pj->v[0] - si->v[0],  pj->v[1] - si->v[1],
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

    const float dv_physical_squared =  dv_physical[0] * dv_physical[0]
      + dv_physical[1] * dv_physical[1]
      + dv_physical[2] * dv_physical[2];

    /* Compute the physical distance between the particles */
    const float dx_physical[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a,
                                  dx[2] * cosmo->a};
    const float r_physical = r * cosmo->a;

    /* Momentum check */
    /* Relative momentum of the gas */
    const float specific_angular_momentum_gas[3] = {
        dx_physical[1] * dv_physical[2] - dx_physical[2] * dv_physical[1],
        dx_physical[2] * dv_physical[0] - dx_physical[0] * dv_physical[2],
        dx_physical[0] * dv_physical[1] - dx_physical[1] * dv_physical[0]};
    const float L2_gas =
        specific_angular_momentum_gas[0] * specific_angular_momentum_gas[0] +
        specific_angular_momentum_gas[1] * specific_angular_momentum_gas[1] +
        specific_angular_momentum_gas[2] * specific_angular_momentum_gas[2];

    /* Keplerian angular speed squared */
    const float omega_acc_2 = grav_props->G_Newton * si->mass /
                              (r_physical * r_physical * r_physical);

    /*Keplerian angular momentum squared */
    const float L2_acc =
        (si->r_cut * si->r_cut * si->r_cut * si->r_cut) * omega_acc_2;

    /* To be accreted, the gas momentum shoulb lower than the keplerian orbit
     * momentum. */
    if (L2_gas > L2_acc) {
      return;
    }

    /* Energy check */
    /* Kinetic energy of the gas */
    float E_kin_relative_gas =
        0.5f * dv_physical_squared;

    /* Compute the Newtonian or truncated potential the sink exherts onto the
       gas particle */
    const float eps = gravity_get_softening(si->gpart, grav_props);
    const float eps2 = eps * eps;
    const float eps_inv = 1.f / eps;
    const float eps_inv3 = eps_inv * eps_inv * eps_inv;
    const float sink_mass = si->mass;
    float dummy, pot_ij;
    runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, sink_mass, &dummy,
                             &pot_ij);

    const double Omega_r = cosmo->Omega_r + cosmo->Omega_nu;
    const double Omega_m = cosmo->Omega_cdm + cosmo->Omega_b;
    const double Omega_l = cosmo->Omega_lambda;
    const double w0 = cosmo->w_0;
    const double wa = cosmo->w_a;
    const double a_inv = cosmo->a_inv;

    const double w_DE =  w0 + wa * (1. - a) ; //cosmology_dark_energy_EoS(a, w0, wa);
    const double w_tilde = (a - 1.) * wa - (1. + w0 + wa) * log(a);
    const double density_sum = Omega_m * a_inv * a_inv * a_inv
                               + 2.0*Omega_r* a_inv * a_inv * a_inv * a_inv
                               + Omega_l * exp(3. * w_tilde) * (1 + w_DE);
                               //w_tilde(a, w0, wa)
    const double a_dot_dot = - H*H/2.0 * density_sum;

    /* Compute the physical potential energy that the sink exerts in the gas :
                       E_pot_phys = G*pot_grav*a^(-1) + c(a). */
    /* The normalization is c(a) = -G*a_dot*a*r^2/2.0. */
    const float E_pot_gas = grav_props->G_Newton * pot_ij * cosmo->a_inv
                       - grav_props->G_Newton * a_dot_dot*a*r2/2.0;

    /* Update: Add thermal energy to avoid the sink to swallow hot gas regions */
    const float E_therm = hydro_get_drifted_physical_internal_energy(pj, cosmo);

    /* Energy of the pair sink-gas */
    const float E_mec_sink_part = E_kin_relative_gas + E_pot_gas + E_therm;

    /* To be accreted, the gas must be gravitationally bound to the sink. */
    if (E_mec_sink_part >= 0) return;

    /* Most bound pair check */
    /* The pair gas-sink must be the most bound among all sinks */
    if (E_mec_sink_part >= pj->sink_data.E_mec_bound) {
      return;
    }

    /* Since this pair gas-sink is the most bounf, keep track of the
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
