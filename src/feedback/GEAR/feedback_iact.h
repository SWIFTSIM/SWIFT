/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_FEEDBACK_IACT_H
#define SWIFT_GEAR_FEEDBACK_IACT_H

/* Local includes */
#include "feedback.h"
#include "hydro.h"
#include "random.h"
#include "timestep_sync_part.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *si, const struct part *pj,
                                    const struct xpart *xpj,
                                    const struct cosmology *cosmo,
                                    const struct feedback_props *fb_props,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  /* The normalization by 1 / h^d is done in feedback.h */
  si->feedback_data.enrichment_weight += mj * wi;
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const integertime_t ti_current) {

  const double e_sn = si->feedback_data.energy_ejected;
  const double e_preSN = si->feedback_data.preSN.energy_ejected;

  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);

  /* Get the kernel for hi. */
  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;

  /* Compute inverse enrichment weight */
  const double si_inv_weight = si->feedback_data.enrichment_weight == 0
                                   ? 0.
                                   : 1. / si->feedback_data.enrichment_weight;

  const double weight = mj * wi * si_inv_weight;
  double m_ej = 0.0;
  double new_mass = 0.0;

  /* Distribute pre-SN */
  if (e_preSN != 0.0 && weight > 0.0) {

    /* Mass received by Stellar Winds */
    /* For physical consistency, we consider that the pre-SN feedback occurs
       before the SN feedback, not at the same time ! Thus, we have to perform
       the calculation first with the mass ejected by the winds, otherwise the
       effects of the winds would be artificially boosted by the SNe ! (It is
       similar to not do the pre-SN feedback in the same step as SN, in that
       case the mass of the gas is not impacted by the mass ejected by SNe.) */
    m_ej = si->feedback_data.preSN.mass_ejected;
    const double dm_SW = m_ej * weight;
    new_mass = mj + dm_SW;
    /* If the distance is null, no need to use calculation ressources. */
    if (r2 > 0.0 && new_mass > 0.0) {
      /* -------------------- set to physical quantities -------------------- */
      /* Cosmology constant */
      const float a = cosmo->a;
      const float a_inv = cosmo->a_inv;
      const float H = cosmo->H;
      const float a_dot = a * H;

      /* physical velocities of the star particle i */
      const float v_i_p[3] = {a_dot * si->x[0] + si->v[0] * a_inv,
                              a_dot * si->x[1] + si->v[1] * a_inv,
                              a_dot * si->x[2] + si->v[2] * a_inv};

      /* physical velocities of the gas particle j */
      const float v_j_p[3] = {a_dot * pj->x[0] + xpj->v_full[0] * a_inv,
                              a_dot * pj->x[1] + xpj->v_full[1] * a_inv,
                              a_dot * pj->x[2] + xpj->v_full[2] * a_inv};

      const float r_p = sqrtf(r2) * a;
      const float dx_p[3] = {dx[0] * a, dx[1] * a, dx[2] * a};

      /* --------------- Compute physical momentum received
       * ---------------------------- */
      /* Total momentum ejected by the winds during the timestep from the star
       * particle i */
      const double p_ej = sqrt(2.0 * si->feedback_data.preSN.mass_ejected *
                               si->feedback_data.preSN.energy_ejected);

      /* norm of physical velocities of the gas particle j */
      const float norm2_v_p =
          v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2];

      double dp_lab_frame[3];

      for (int i = 0; i < 3; i++) {
        /* the unit direction from the gas particle j to the star particle i */
        const double unit_direction = dx_p[i] / r_p;
        /* the additional momentum due to change of frame of reference (from
         * star particle frame to lab frame) */
        const double change_of_frame_dp =
            si->feedback_data.preSN.mass_ejected * v_i_p[i];
        /* momentum in lab frame due to the ejecta */
        dp_lab_frame[i] = weight * (p_ej + change_of_frame_dp) * unit_direction;

        /* Give the comoving momentum to the gas particle */
        xpj->feedback_data.delta_p[i] -=
            dp_lab_frame[i] *
            a; /* The minus sign comes from the direction of dx (si - pj) */
      }

      const double norm2_dp_lab_frame = dp_lab_frame[0] * dp_lab_frame[0] +
                                        dp_lab_frame[1] * dp_lab_frame[1] +
                                        dp_lab_frame[2] * dp_lab_frame[2];
      const double norm2_dp = weight * weight * p_ej * p_ej;

      /* ------------------ calculate physical Energy and internal Energy
       * received -------------------------- */

      /* The energy ejected from the star particle i by stellar wind that is
       * actually received by the gas particle j */
      const double weighted_energy = weight * e_preSN;
      /* The additional energy received by the gas particle j due to the
       * momentum of the star particle i */
      const double dE_change_of_frame =
          0.5 * (norm2_dp_lab_frame - norm2_dp) /
          (weight * si->feedback_data.preSN.mass_ejected);
      /* The total energy received from the gas particle j in the laboratory
       * frame of reference */
      const double dE_lab_frame = weighted_energy + dE_change_of_frame;

      /* The momentum of the gas particle j after receiving the momentum from
       * stellar wind */
      const double p_new[3] = {pj->mass * v_j_p[0] + dp_lab_frame[0],
                               pj->mass * v_j_p[1] + dp_lab_frame[1],
                               pj->mass * v_j_p[2] + dp_lab_frame[2]};
      const double norm2_p_new = {p_new[0] * p_new[0] + p_new[1] * p_new[1] +
                                  p_new[2] * p_new[2]};

      /* The new and old kinetic energy of the gas particle j */
      const double new_kinetic_energy = 0.5 * norm2_p_new / new_mass;
      const double old_kinetic_energy = 0.5 * pj->mass * norm2_v_p;

      /* The additional specific internal energy of the gas particle j.
        Ekin_new + U_new = Ekin_old + U_old + dEtot
        -> du = (U_new - U_old) / new_mass = (Ekin_old + dEtot - Ekin_new) /
        new_mass */
      const float du =
          (old_kinetic_energy + dE_lab_frame - new_kinetic_energy) / new_mass;
      xpj->feedback_data.delta_u += du;

      /* Only used in non-cosmological simulations. Has to be
         investigated in cosmological simulations*/
      if (a == 1.0 && a_inv == 1.0 && cosmo->z == 0.0) {
        /* Calculate the velocity without Hubble flow for signal velocity */
        const float v_i_without_Hubble_flow[3] = {
            si->v[0] * a_inv, si->v[1] * a_inv, si->v[2] * a_inv};
        double dp_without_Hubble[3];
        for (int i = 0; i < 3; i++) {
          /* the unit direction from the gas particle j to the star particle i
           */
          const double unit_direction = dx_p[i] / r_p;
          /* the additional momentum due to change of frame of reference (from
           * star particle frame to lab frame) */
          const double change_of_frame_without_Hubble =
              si->feedback_data.preSN.mass_ejected * v_i_without_Hubble_flow[i];
          /* momentum in lab frame due to the ejecta */
          dp_without_Hubble[i] =
              weight * (p_ej + change_of_frame_without_Hubble) * unit_direction;
        }
        /* The norm of the momentum without the Hubble flow participation */
        const double norm2_dp_without_Hubble = {
            dp_without_Hubble[0] * dp_without_Hubble[0] +
            dp_without_Hubble[1] * dp_without_Hubble[1] +
            dp_without_Hubble[2] * dp_without_Hubble[2]};

        /* Update the signal velocity of the gas particle receiving a kick.
           We want to subtract the Hubble flow participation in the signal
           velocity.*/
        const float dv_phys = sqrtf(norm2_dp_without_Hubble) / new_mass;
        hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, dv_phys);
      }

      xpj->feedback_data.hit_by_preSN = 1;
    }
  }

  /* Distribute SN */
  if (e_sn != 0.0) {

    /* Mass received by SN */
    /* For the conservation of mass and energy, we perform the calculation only
     * with the mass actually ejected by the SN (not the combination of pre-SN
     * and SN) */
    m_ej = si->feedback_data.mass_ejected;
    const double dm_SN = m_ej * weight;
    /* But we are considering that the stellar wind occurs before the SN, so the
       total new mass to take into account is the combination of both. It is
       similar to not do pre-SN feedback in the same step as SN, in that case
       the mass of the gas is the one after receiving mass from the pre-SN
       feedback at the previous step. */
    new_mass += dm_SN;

    /* Energy received */
    const double du = (e_sn)*weight / new_mass;
    xpj->feedback_data.delta_u += du;

    /* Compute momentum received. */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.delta_p[i] += dm * (si->v[i] - xpj->v_full[i]);
    }

    /* Add the metals */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      pj->chemistry_data.metal_mass[i] +=
          weight * si->feedback_data.metal_mass_ejected[i];
    }

    /* Set the indication of SN event for cooling*/
    xpj->feedback_data.hit_by_SN = 1;
  }

  if (xpj->feedback_data.hit_by_preSN || xpj->feedback_data.hit_by_SN) {
    /* Update the mass of the gas particle */
    xpj->feedback_data.delta_mass += dm_SW + dm_SN;
  }

  /* Impose maximal viscosity (only for SN) */
  if (xpj->feedback_data.hit_by_SN) {
    hydro_diffusive_feedback_reset(pj);
  }

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}

#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
