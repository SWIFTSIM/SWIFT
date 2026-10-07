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
#include "../GEAR/radiation_iact.h"
#include "../GEAR/radiation_propagation_iact.h"
#include "feedback.h"
#include "feedback_tracers.h"
#include "hydro.h"
#include "random.h"
#include "timestep_sync_part.h"
#include "tracers.h"

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
runner_iact_nonsym_feedback_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, const struct part *pj, const struct xpart *xpj,
    const struct cosmology *cosmo, const struct feedback_props *fb_props,
    const struct hydro_props *hydro_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
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

  /* Contribution to the number of neighbours */
  si->feedback_data.num_ngbs += 1;

  /*****************************************/
  /* Radiation */
  radiation_iact_nonsym_feedback_density(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                         fb_props, hydro_props, phys_const, us,
                                         cooling, ti_current);
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
 * @param hydro_props The properties of the hydro scheme.
 * @param fb_props Properties of the feedback scheme.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The properties of the cooling scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 * @param time_base The time base used to compute integer times.
 * @param with_cosmology Are we running with cosmology on?
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current, const double time_base,
    const int with_cosmology) {

  const double e_sn = si->feedback_data.supernovae.energy_ejected;
  const double e_winds = si->feedback_data.winds.energy_ejected;

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
  double new_mass = mj;
  double dm_SW = 0.0;
  double dm_SN = 0.0;

  /*****************************************/
  /* Radiation */
  /* TODO: Add hit by radiation */
  radiation_iact_nonsym_feedback_apply(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                       hydro_props, fb_props, phys_const, us,
                                       cooling, ti_current);

  /* Distribute pre-SN */
  if (e_winds != 0.0 && weight > 0.0) {

    /* Mass received by Stellar Winds */
    /* For physical consistency, we consider that the pre-SN feedback occurs
       before the SN feedback, not at the same time ! Thus, we have to perform
       the calculation first with the mass ejected by the winds, otherwise the
       effects of the winds would be artificially boosted by the SNe ! (It is
       similar to not do the pre-SN feedback in the same step as SN, in that
       case the mass of the gas is not impacted by the mass ejected by SNe.) */
    m_ej = si->feedback_data.winds.mass_ejected;
    dm_SW = m_ej * weight;
    new_mass += dm_SW;

    /* If the distance is null, no need to use calculation ressources. */
    if (r2 > 0.0 && new_mass > 0.0) {
      /* -------------------- set to physical quantities -------------------- */
      /* Cosmology constant */
      const float a = cosmo->a;
      const float a_inv = cosmo->a_inv;
      const float H = cosmo->H;
      const float a_dot = a * H;

      /* Physical velocities of the star particle i. The Hubble flow term is
         relative wrt the star particle. Hence, for the star, dx = 0 ; for the
         gas, we use -dx = pj - si. */
      const float v_i_p[3] = {si->v[0] * a_inv, si->v[1] * a_inv,
                              si->v[2] * a_inv};

      /* Physical peculiar velocity of the gas particle j. */
      const float v_j_pec[3] = {xpj->v_full[0] * a_inv, xpj->v_full[1] * a_inv,
                                xpj->v_full[2] * a_inv};

      /* Physical velocities of the gas particle j, with the Hubble flow. */
      const float v_j_p[3] = {-a_dot * dx[0] + v_j_pec[0],
                              -a_dot * dx[1] + v_j_pec[1],
                              -a_dot * dx[2] + v_j_pec[2]};

      const float r_p = sqrtf(r2) * a;
      const float dx_p[3] = {dx[0] * a, dx[1] * a, dx[2] * a};

      /* --------------- Compute physical momentum received ---------------- */
      /* Total momentum ejected by the winds during the timestep from the star
       * particle i */
      const float p_ej = sqrt(2.0 * si->feedback_data.winds.mass_ejected *
                              si->feedback_data.winds.energy_ejected);

      /* Ejecta momentum in the rest frame of the star (away from it, dx points
         to the star), then in the lab frame (the ejected mass carries the star
         velocity), then in the frame of the gas, which is the one applied
         (v += delta_p / m_f): dm_SW also has to be accelerated from v_j. */
      double delta_p_star_frame[3];
      double delta_p_lab_frame[3];
      double delta_p_gas_frame[3];

      for (int i = 0; i < 3; i++) {
        delta_p_star_frame[i] = -weight * p_ej * dx_p[i] / r_p;
        delta_p_lab_frame[i] = delta_p_star_frame[i] + dm_SW * v_i_p[i];
        delta_p_gas_frame[i] = delta_p_lab_frame[i] - dm_SW * v_j_pec[i];

        /* Give the comoving momentum to the gas particle */
        xpj->feedback_data.delta_p[i] += delta_p_gas_frame[i] * a;
      }

      const double norm2_delta_p_gas_frame =
          delta_p_gas_frame[0] * delta_p_gas_frame[0] +
          delta_p_gas_frame[1] * delta_p_gas_frame[1] +
          delta_p_gas_frame[2] * delta_p_gas_frame[2];

      /* ----- Calculate physical Energy and internal Energy received ------ */

      /* The thermal energy of the gas particle j is the energy of the ejecta
         (wind energy and kinetic energy of the frame change) plus its kinetic
         energy, minus the kinetic energy after the momentum is shared:
           dE_th = Ekin_old + dE_lab - Ekin_new
                 = 0.5 m dm / m_f |u_rel|^2,
         with u_rel the velocity of the ejecta relative to the gas (the gas has
         the Hubble flow around the star). This form has no cancellation of
         large terms. The update divides the sum over the events by the final
         mass. Without ejected mass, only the wind energy is given. */
      double dE_th = weight * e_winds;
      if (dm_SW > 0.0) {
        double norm2_u_rel = 0.0;
        for (int i = 0; i < 3; i++) {
          const double u_rel = delta_p_lab_frame[i] / dm_SW - v_j_p[i];
          norm2_u_rel += u_rel * u_rel;
        }
        dE_th = 0.5 * mj * dm_SW / new_mass * norm2_u_rel;
      }
      xpj->feedback_data.delta_E_th += dE_th;

      /* Only used in non-cosmological simulations. Has to be
         investigated in cosmological simulations*/
      if (a == 1.0 && a_inv == 1.0 && cosmo->z == 0.0) {
        /* Update the signal velocity of the gas particle receiving a kick. The
           momentum applied has no Hubble flow in it. */
        const float dv_phys = sqrt(norm2_delta_p_gas_frame) / new_mass;
        hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, dv_phys);
      }

      /* Lifetime-cumulative tracer, using this branch's own locally-computed
         momentum/energy (not the shared feedback_data.delta_p/delta_E_th,
         which the SN branch below can also add to this same step). The
         momentum is the one applied, m_f |dv|, in the frame of the gas as for
         the SN. The GEAR tracers get the sums at the update, with the final
         mass. */
      const float delta_p_mag_winds = (float)sqrt(norm2_delta_p_gas_frame);
      feedback_tracers_event_SW(xpj, delta_p_mag_winds, dE_th, new_mass);

      xpj->feedback_data.hit_by_winds = 1;
    }
  }

  /* Distribute SN. The mass is a condition in its own right: with zero SN
     energy the ejected mass, already removed from the star, must still reach
     the gas. */
  if (e_sn != 0.0 || si->feedback_data.supernovae.mass_ejected != 0.0) {

    /* Mass received by SN */
    /* For the conservation of mass and energy, we perform the calculation only
     * with the mass actually ejected by the SN (not the combination of pre-SN
     * and SN) */
    m_ej = si->feedback_data.supernovae.mass_ejected;
    dm_SN = m_ej * weight;

    /* But we are considering that the stellar wind occurs before the SN, so the
       total new mass to take into account is the combination of both. It is
       similar to not do pre-SN feedback in the same step as SN, in that case
       the mass of the gas is the one after receiving mass from the pre-SN
       feedback at the previous step. */
    new_mass += dm_SN;

    /* Energy received */
    const double dE_th = e_sn * weight;
    xpj->feedback_data.delta_E_th += dE_th;

    /* Compute momentum received. */
    float delta_p_supernovae[3];
    for (int i = 0; i < 3; i++) {
      delta_p_supernovae[i] = dm_SN * (si->v[i] - xpj->v_full[i]);
      xpj->feedback_data.delta_p[i] += delta_p_supernovae[i];
    }

    /* Add the metals */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      xpj->feedback_data.delta_metal_mass[i] +=
          weight * si->feedback_data.metal_mass_ejected[i];

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
      pj->feedback_data.metal_mass[i] +=
          weight * si->feedback_data.metal_mass_ejected[i];
#endif
    }

    /* delta_p_supernovae is comoving; a_inv gives the physical momentum
       actually applied (matches feedback_update_part()'s v_full += p/m). */
    const float delta_p_mag_supernovae_comoving =
        sqrtf(delta_p_supernovae[0] * delta_p_supernovae[0] +
              delta_p_supernovae[1] * delta_p_supernovae[1] +
              delta_p_supernovae[2] * delta_p_supernovae[2]);
    const float delta_p_mag_supernovae =
        delta_p_mag_supernovae_comoving * cosmo->a_inv;
    feedback_tracers_event_SN(xpj, delta_p_mag_supernovae, dE_th, new_mass);

    /* Flag the thermal event for cooling: it tracks the injected energy,
       not the mass. */
    if (e_sn != 0.0) xpj->feedback_data.hit_by_SN = 1;
  }

  /* Must not depend on hit_by_SN: mass can arrive with no energy. */
  xpj->feedback_data.delta_mass += dm_SW + dm_SN;

  /* Impose maximal viscosity (only for SN) */
  if (xpj->feedback_data.hit_by_SN) {
    hydro_diffusive_feedback_reset(pj);
  }

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}

#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
