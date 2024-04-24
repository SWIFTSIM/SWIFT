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
  /* const float mj = hydro_get_mass(pj); */

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  /* const float hi_inv = 1.0f / hi; */
  /* const float ui = r * hi_inv; */
  /* float wi; */
  /* kernel_eval(ui, &wi); */

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  /* The normalization by 1 / h^d is done in feedback.h */
  /* si->feedback_data.enrichment_weight += mj * wi; */


  /* Recompute the scalar_weight */
  const float r_j2 = pj->x[0]*si->x[0] + pj->x[1]*si->x[1] + pj->x[1]*si->x[1];

  float dW_ij_dr_j;
  float dW_jj_dr_j;
  const float u_ij = sqrt(r_j2)/hi;
  const float u_jj = sqrt(r_j2)/hj;
  float dummy_W;

  kernel_deval(u_ij, &dummy_W, &dW_ij_dr_j);
  kernel_deval(u_jj, &dummy_W, &dW_jj_dr_j);

  dW_ij_dr_j *= pow_dimension(1.0/hi); /* 1/h_i^d */
  dW_jj_dr_j *= pow_dimension(1.0/hj); /* 1/h_j^d */

  /* Compute the projection vectors */
  const double dx_ij_plus[3] = {max(dx[0], 0.0)/r,
			  max(dx[1], 0.0)/r,
			  max(dx[2], 0.0)/r};

  const double dx_ij_minus[3] = {min(dx[0], 0.0)/r,
				 min(dx[1], 0.0)/r,
				 min(dx[2], 0.0)/r};

  /* Correct in stars_iact the 1/r factor (remove it) */
  const double dx_ij_hat[3] = {dx_ij_plus[0] + dx_ij_minus[0],
			       dx_ij_plus[1] + dx_ij_minus[1],
			       dx_ij_plus[2] + dx_ij_minus[2]};

  /* I also need the x_ij_hat... which must be calculated before... */
  double n_bar_i_2_inv = 1.0/(si->density.wcount*si->density.wcount);
  double n_bar_j_2_inv = 1.0/(pj->density.wcount*pj->density.wcount);

  double A_j[3] = {(n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[0],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[1],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[2]};

  double number_1 = A_j[0]*dx_ij_hat[0] *  A_j[1]*dx_ij_hat[1] *  A_j[2]*dx_ij_hat[2];
  double number_2 = M_PI * r2;
  double denom = sqrt(1 + number_1/number_2);

  /* Store at least this value. I would not be against storing dx_ij and A to
     avoid computing them many times and making mistakes... */
  double scalar_weight_j = 0.5*(1.0 + 1.0/denom) ;


  /* Now, that we have accumulated the sums, we can compute the f_plus and
     f_minus */

  const double value_plus[3] = {1 + (si->feedback_data.f_plus_num[0]*si->feedback_data.f_plus_num[0])/(si->feedback_data.f_plus_denom[0]*si->feedback_data.f_plus_denom[0]),
		     1 + (si->feedback_data.f_plus_num[1]*si->feedback_data.f_plus_num[1])/(si->feedback_data.f_plus_denom[1]*si->feedback_data.f_plus_denom[1]),
		     1 + (si->feedback_data.f_plus_num[2]*si->feedback_data.f_plus_num[2])/(si->feedback_data.f_plus_denom[2]*si->feedback_data.f_plus_denom[2])};

  const double f_plus_i[3] = {sqrt(0.5*value_plus[0]), sqrt(0.5*value_plus[1]), sqrt(0.5*value_plus[2])};

  const double value_minus[3] = {1 + (si->feedback_data.f_minus_num[0]*si->feedback_data.f_minus_num[0])/(si->feedback_data.f_minus_denom[0]*si->feedback_data.f_minus_denom[0]),
		     1 + (si->feedback_data.f_minus_num[1]*si->feedback_data.f_minus_num[1])/(si->feedback_data.f_minus_denom[1]*si->feedback_data.f_minus_denom[1]),
		     1 + (si->feedback_data.f_minus_num[2]*si->feedback_data.f_minus_num[2])/(si->feedback_data.f_minus_denom[2]*si->feedback_data.f_minus_denom[2])};

  const double f_minus_i[3] = {sqrt(0.5*value_minus[0]), sqrt(0.5*value_minus[1]), sqrt(0.5*value_minus[2])};

  /* Assign the f_plus and f_minus to the star */
  si->feedback_data.f_plus[0] = f_plus_i[0];
  si->feedback_data.f_plus[1] = f_plus_i[1];
  si->feedback_data.f_plus[2] = f_plus_i[2];

  si->feedback_data.f_minus[0] = f_minus_i[0];
  si->feedback_data.f_minus[1] = f_minus_i[1];
  si->feedback_data.f_minus[2] = f_minus_i[2];

  /* Now compute the w_j (vector) */
  double w_j[3];

  /* Dans le texte, on somme sur +/- et sur alpha=x,y,z. Mais ca ne donne pas
  un vecteur. Je pense que la somme sur alpha est fausse et c'est juste
  multiplier par les composantes */
  w_j[0] = scalar_weight_j*(dx_ij_plus[0]*f_plus_i[0] + dx_ij_minus[0]*f_minus_i[0]);
  w_j[1] = scalar_weight_j*(dx_ij_plus[1]*f_plus_i[1] + dx_ij_minus[1]*f_minus_i[1]);
  w_j[2] = scalar_weight_j*(dx_ij_plus[2]*f_plus_i[2] + dx_ij_minus[1]*f_minus_i[2]);

  /* Accumulate w_j norm for later */
  const double w_j_norm_2 = w_j[0]*w_j[0] + w_j[1]*w_j[1] + w_j[2]*w_j[2];
  si->feedback_data.sum_vector_weight_norm += sqrt(w_j_norm_2);
}


__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  const struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {
  message("HELLO");
}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {}

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

  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);

  /* Get the kernel for hi. */
  /* float hi_inv = 1.0f / hi; */
  /* float hi_inv_dim = pow_dimension(hi_inv); /\* 1/h^d *\/ */
  /* float xi = r * hi_inv; */
  /* float wi, wi_dx; */
  /* kernel_deval(xi, &wi, &wi_dx); */
  /* wi *= hi_inv_dim; */

  /* Compute inverse enrichment weight */
  /* const double si_inv_weight = si->feedback_data.enrichment_weight == 0 */
  /*                                  ? 0. */
  /*                                  : 1. / si->feedback_data.enrichment_weight; */

  /* /\* Mass received *\/ */
  /* const double m_ej = si->feedback_data.mass_ejected; */
  /* const double weight = mj * wi * si_inv_weight; */
  /* const double dm = m_ej * weight; */
  /* const double new_mass = mj + dm; */

  /* /\* Energy received *\/ */
  /* const double du = e_sn * weight / new_mass; */

  /* xpj->feedback_data.delta_mass += dm; */
  /* xpj->feedback_data.delta_u += du; */

  /* /\* Compute momentum received. *\/ */
  /* for (int i = 0; i < 3; i++) { */
  /*   xpj->feedback_data.delta_p[i] += dm * (si->v[i] - xpj->v_full[i]); */
  /* } */

  /* /\* Add the metals *\/ */
  /* for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) { */
  /*   pj->chemistry_data.metal_mass[i] += */
  /*       weight * si->feedback_data.metal_mass_ejected[i]; */
  /* } */



  const float r_j2 = pj->x[0]*si->x[0] + pj->x[1]*si->x[1] + pj->x[1]*si->x[1];

  float dW_ij_dr_j;
  float dW_jj_dr_j;
  const float u_ij = sqrt(r_j2)/hi;
  const float u_jj = sqrt(r_j2)/hj;
  float dummy_W;

  kernel_deval(u_ij, &dummy_W, &dW_ij_dr_j);
  kernel_deval(u_jj, &dummy_W, &dW_jj_dr_j);

  dW_ij_dr_j *= pow_dimension(1.0/hi); /* 1/h_i^d */
  dW_jj_dr_j *= pow_dimension(1.0/hj); /* 1/h_j^d */

  /* Compute the projection vectors */
  const double dx_ij_plus[3] = {max(dx[0], 0.0)/r,
			  max(dx[1], 0.0)/r,
			  max(dx[2], 0.0)/r};

  const double dx_ij_minus[3] = {min(dx[0], 0.0)/r,
				 min(dx[1], 0.0)/r,
				 min(dx[2], 0.0)/r};

  /* Correct in stars_iact the 1/r factor (remove it) */
  const double dx_ij_hat[3] = {dx_ij_plus[0] + dx_ij_minus[0],
			       dx_ij_plus[1] + dx_ij_minus[1],
			       dx_ij_plus[2] + dx_ij_minus[2]};

  /* I also need the x_ij_hat... which must be calculated before... */
  double n_bar_i_2_inv = 1.0/(si->density.wcount*si->density.wcount);
  double n_bar_j_2_inv = 1.0/(pj->density.wcount*pj->density.wcount);

  double A_j[3] = {(n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[0],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[1],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[2]};

  double number_1 = A_j[0]*dx_ij_hat[0] *  A_j[1]*dx_ij_hat[1] *  A_j[2]*dx_ij_hat[2];
  double number_2 = M_PI * r2;
  double denom = sqrt(1 + number_1/number_2);

  /* Store at least this value. I would not be against storing dx_ij and A to
     avoid computing them many times and making mistakes... */
  double scalar_weight_j = 0.5*(1.0 + 1.0/denom) ;


  /* Now, that we have accumulated the sums, we can compute the f_plus and
     f_minus */
  const double value_plus[3] = {1 + (si->feedback_data.f_plus_num[0]*si->feedback_data.f_plus_num[0])/(si->feedback_data.f_plus_denom[0]*si->feedback_data.f_plus_denom[0]),
		     1 + (si->feedback_data.f_plus_num[1]*si->feedback_data.f_plus_num[1])/(si->feedback_data.f_plus_denom[1]*si->feedback_data.f_plus_denom[1]),
		     1 + (si->feedback_data.f_plus_num[2]*si->feedback_data.f_plus_num[2])/(si->feedback_data.f_plus_denom[2]*si->feedback_data.f_plus_denom[2])};

  const double f_plus_i[3] = {sqrt(0.5*value_plus[0]), sqrt(0.5*value_plus[1]), sqrt(0.5*value_plus[2])};

  const double value_minus[3] = {1 + (si->feedback_data.f_minus_num[0]*si->feedback_data.f_minus_num[0])/(si->feedback_data.f_minus_denom[0]*si->feedback_data.f_minus_denom[0]),
		     1 + (si->feedback_data.f_minus_num[1]*si->feedback_data.f_minus_num[1])/(si->feedback_data.f_minus_denom[1]*si->feedback_data.f_minus_denom[1]),
		     1 + (si->feedback_data.f_minus_num[2]*si->feedback_data.f_minus_num[2])/(si->feedback_data.f_minus_denom[2]*si->feedback_data.f_minus_denom[2])};

  const double f_minus_i[3] = {sqrt(0.5*value_minus[0]), sqrt(0.5*value_minus[1]), sqrt(0.5*value_minus[2])};

  /* Now compute the w_j (vector) */
  double w_j[3];

  /* Dans le texte, on somme sur +/- et sur alpha=x,y,z. Mais ca ne donne pas
  un vecteur. Je pense que la somme sur alpha est fausse et c'est juste
  multiplier par les composantes */
  w_j[0] = scalar_weight_j*(dx_ij_plus[0]*f_plus_i[0] + dx_ij_minus[0]*f_minus_i[0]);
  w_j[1] = scalar_weight_j*(dx_ij_plus[1]*f_plus_i[1] + dx_ij_minus[1]*f_minus_i[1]);
  w_j[2] = scalar_weight_j*(dx_ij_plus[2]*f_plus_i[2] + dx_ij_minus[1]*f_minus_i[2]);

  double w_j_bar[3] =  {w_j[0]/si->feedback_data.sum_vector_weight_norm,
			w_j[1]/si->feedback_data.sum_vector_weight_norm,
			w_j[2]/si->feedback_data.sum_vector_weight_norm};

  const double w_j_bar_norm_2 = w_j_bar[0]*w_j_bar[0] + w_j_bar[1]*w_j_bar[1] + w_j_bar[2]*w_j_bar[2];
  const double w_j_bar_norm = sqrt(w_j_bar_norm_2);


  /* */
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

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}


#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
