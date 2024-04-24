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

  /* The normalization by 1 / h^d is done in feedback.h */
  /* si->feedback_data.enrichment_weight += mj * wi; */

  /* Now we can compute f_plus and f_minus for the star */
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj, f_plus_i,
						f_minus_i, w_j);

  /* Assign the f_plus and f_minus to the star */
  si->feedback_data.f_plus[0] = f_plus_i[0];
  si->feedback_data.f_plus[1] = f_plus_i[1];
  si->feedback_data.f_plus[2] = f_plus_i[2];

  si->feedback_data.f_minus[0] = f_minus_i[0];
  si->feedback_data.f_minus[1] = f_minus_i[1];
  si->feedback_data.f_minus[2] = f_minus_i[2];

  /* Accumulate w_j norm for later */
  const double w_j_norm_2 = w_j[0]*w_j[0] + w_j[1]*w_j[1] + w_j[2]*w_j[2];
  si->feedback_data.sum_vector_weight_norm += sqrt(w_j_norm_2);
}


/* __attribute__((always_inline)) INLINE static void */
/* runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3], */
/*                                   const float hi, const float hj, */
/*                                   const struct spart *si, struct part *pj, */
/*                                   const struct xpart *xpj, */
/*                                   const struct cosmology *cosmo, */
/*                                   const integertime_t ti_current) {} */

/* __attribute__((always_inline)) INLINE static void */
/* runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3], */
/*                                   const float hi, const float hj, */
/*                                   struct spart *si, const struct part *pj, */
/*                                   const struct xpart *xpj, */
/*                                   const struct cosmology *cosmo, */
/*                                   const integertime_t ti_current) {} */

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

  /* Do we have supernovae? */
  if (e_sn == 0) {
    return;
  }

  /* Finally, we can compute the w_j_bar. */
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj, f_plus_i, f_minus_i, w_j);

  double w_j_bar[3] =  {w_j[0]/si->feedback_data.sum_vector_weight_norm,
			w_j[1]/si->feedback_data.sum_vector_weight_norm,
			w_j[2]/si->feedback_data.sum_vector_weight_norm};

  const double w_j_bar_norm_2 = w_j_bar[0]*w_j_bar[0] + w_j_bar[1]*w_j_bar[1] + w_j_bar[2]*w_j_bar[2];
  const double w_j_bar_norm = sqrt(w_j_bar_norm_2);

  /* Update the xpart properties */
  const float mj = hydro_get_mass(pj);
  const double m_ej = si->feedback_data.mass_ejected;
  const double p_ej = sqrt(2*m_ej*e_sn) ;

  /* Energy received */
  double dm = w_j_bar_norm * m_ej;
  xpj->feedback_data.delta_mass += dm;

  /* Add the metals */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    pj->chemistry_data.metal_mass[i] +=
	w_j_bar_norm * si->feedback_data.metal_mass_ejected[i];
  }

  /* Update Momentum */
  const double dp[3] = {w_j_bar[0]*p_ej, w_j_bar[1]*p_ej, w_j_bar[2]*p_ej};
  const double dp_prime[3] = {dp[0] + dm*si->v[0], dp[1] + dm*si->v[1], dp[2] + dm*si->v[2] };
  const double dp_norm_2 = dp[0]*dp[0] +  dp[2]*dp[2] +  dp[2]*dp[2];
  const double dp_prime_norm_2 = dp_prime[0]*dp_prime[0] +  dp_prime[2]*dp_prime[2]
				 +  dp_prime[2]*dp_prime[2];

  for (int i = 0; i < 3; i++) {
    xpj->feedback_data.delta_p[i] += dp_prime[i];
  }

  /* Update energy */
  const double dE = w_j_bar_norm * e_sn;
  const double dE_prime = dE + 1.0/dm * (dp_prime_norm_2 - dp_norm_2);

  const double new_mass = mj + dm;

  xpj->feedback_data.delta_u += dE_prime/new_mass;

  /* Verify conservation things */
  si->feedback_data.delta_m_check += dm;
  si->feedback_data.delta_E_check += dE;
  si->feedback_data.delta_p_norm_check += sqrt(dp_norm_2);

  si->feedback_data.delta_p_check[0] += dp[0];
  si->feedback_data.delta_p_check[1] += dp[1];
  si->feedback_data.delta_p_check[2] += dp[2];

  message("Conservation check (star %lld): Sum dm_i = %e, Sum dE_i = %e, Sum |dp_i| = %e, Sum dp_i = (%e, %e, %e)", si->id, si->feedback_data.delta_m_check, si->feedback_data.delta_E_check, si->feedback_data.delta_p_norm_check, si->feedback_data.delta_p_check[0], si->feedback_data.delta_p_check[1], si->feedback_data.delta_p_check[2]);

  /* si->feedback_data.delta_p_tot[0] += mj*pj->v[0] + dp_prime[0]; */
  /* si->feedback_data.delta_p_tot[1] += mj*pj->v[1] + dp_prime[1] ; */
  /* si->feedback_data.delta_p_tot[2] += mj*pj->v[1] + dp_prime[2] ; */

  /* const double v_norm_2 = pj->v[0]*sp->v[0] + pj->v[1]*sp->v[1] + sp->v[2]*sp->v[2]; */
  /* sp->feedback_data.E_tot = 0.5*sp->mass*v_norm_2; */

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}


#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
