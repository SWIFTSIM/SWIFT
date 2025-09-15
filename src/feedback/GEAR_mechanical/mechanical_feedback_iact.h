/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_GEAR_MECHANICAL_MECHANICAL_FEEDBACK_IACT_H
#define SWIFT_GEAR_MECHANICAL_MECHANICAL_FEEDBACK_IACT_H

#include "feedback.h"
#include "hydro.h"

__attribute__((always_inline)) INLINE static void
mechanical_feedback_accumulate_fluxes_for_conservation_check(
    struct spart *si, const double dm, const double dp[3], const double m_ej,
    const double p_ej, const double E_ej) {
#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  const float dp_norm_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  si->feedback_data.fluxes_conservation_check.delta_m += dm;
  si->feedback_data.fluxes_conservation_check.delta_p_norm += sqrt(dp_norm_2);

  si->feedback_data.fluxes_conservation_check.delta_p[0] += dp[0];
  si->feedback_data.fluxes_conservation_check.delta_p[1] += dp[1];
  si->feedback_data.fluxes_conservation_check.delta_p[2] += dp[2];

  message(
      "Conservation check (star %lld): Sum dm_i = %e (m_ej), Sum |dp_i| = %e "
      "(p_ej), Sum dp_i = (%e, %e, %e) (0), m_ej = %e, E_ej = %e, p_ej = %e",
      si->id, si->feedback_data.fluxes_conservation_check.delta_m,
      si->feedback_data.fluxes_conservation_check.delta_p_norm,
      si->feedback_data.fluxes_conservation_check.delta_p[0],
      si->feedback_data.fluxes_conservation_check.delta_p[1],
      si->feedback_data.fluxes_conservation_check.delta_p[2], m_ej, E_ej, p_ej);
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */
}

#if FEEDBACK_GEAR_MECHANICAL_MODE == 1
/**
 * @brief Feedback interaction between two particles (non-symmetric) for the
 * mechanical feedback mode 1.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param w_j_bar Normalized vector weight.
 * @param w_j_bar_norm Norm of the normalized vector weight.
 * @param v_i_p Physical velocity of particle i (star).
 * @param v_j_p Physical velocity of particle j (gas).
 * @param E_ej Physical ejected energy by the supernova.
 * @param m_ej Ejected mass by the supernova.
 * @param mj Mass of particle j.
 * @param dm Distributed mass from the SN to gas particle j.
 * @param new_mass Mass of particle j after the SN feedback.
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param constants The physical constants (in internal units).
 * @param us The internal system of units.
 * @param dU (return) Internal energy to distribute to gas particle j.
 * @param dKE (return) Kinetic energy variation before and after the SN
 * feedback.
 * @param dp_prime (return) Momemtum to distribute to gas particle j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_1_feedback_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const double dm, const double new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp_prime[3]) {

  // TODO: Convert to an ifdef constant
  const float internal_energy_snowplow_exponent = -6.5;

  /* ... physical momentum */
  const double p_ej = sqrt(2 * m_ej * E_ej);
  const double dp[3] = {w_j_bar[0] * p_ej, w_j_bar[1] * p_ej,
                        w_j_bar[2] * p_ej};
  const double dE = w_j_bar_norm * E_ej;

  /* Now boost to the 'laboratory' frame */
  dp_prime[0] = dp[0] + dm * v_i_p[0];
  dp_prime[1] = dp[1] + dm * v_i_p[1];
  dp_prime[2] = dp[2] + dm * v_i_p[2];

  /* ... physical total energy */
  const double dp_norm_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  const double dp_prime_norm_2 = dp_prime[0] * dp_prime[0] +
                                 dp_prime[1] * dp_prime[1] +
                                 dp_prime[2] * dp_prime[2];
  const double dE_prime = dE + 1.0 / (2.0 * dm) * (dp_prime_norm_2 - dp_norm_2);

  /* ... physical internal energy */
  /* Compute kinetic energy difference before and after SN */
  const double p_old_norm_2 =
      mj * mj *
      (v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2]);
  const double p_new[3] = {mj * v_j_p[0] + dp_prime[0],
                           mj * v_j_p[1] + dp_prime[1],
                           mj * v_j_p[2] + dp_prime[2]};
  const double p_new_norm_2 =
      p_new[0] * p_new[0] + p_new[1] * p_new[1] + p_new[2] * p_new[2];

  const double E_kin_old = p_old_norm_2 / (2.0 * mj);
  const double E_kin_new = p_new_norm_2 / (2.0 * new_mass);
  *dKE = E_kin_new - E_kin_old;

  const double U_old = hydro_get_physical_internal_energy(pj, xpj, cosmo);
  const double E_old = U_old + E_kin_old;
  const double E_new = E_old + dE_prime;
  const double U_new = E_new - E_kin_new;

  /* Compute the physical internal energy */
  *dU = U_new - U_old;

  /* --Now, we take into account for potentially unresolved energy-conserving
     phase of the SN explosion-- */

  const double PdV_work_fraction = sqrt(1 + mj / dm);
  const double p_terminal = feedback_get_physical_SN_terminal_momentum(
      si, pj, xpj, phys_const, us, cosmo);

  /* If we can resolve the Taylor Sedov, then we give the right coupled
     momentum (which is by definition <= p_terminal). If we cannot resolve it,
     then we have reached the p_terminal. This is the upper limit of momentum,
     since afterwards the cooling is efficient and the thermal energy is
     radiated away. Thus, the factor to multiply dp_prime is: */
  const double p_factor = min(PdV_work_fraction, p_terminal / p_ej);

  dp_prime[0] *= p_factor;
  dp_prime[1] *= p_factor;
  dp_prime[2] *= p_factor;

  /* Compute the comoving cooling radius */
  const float r_cool = cosmo->a_inv * feedback_get_physical_SN_cooling_radius(
                                          si, p_ej, p_terminal, cosmo);

  /* If we do not resolve the Taylor-Sedov, we rescale the internal energy */
  if (r2 > r_cool * r_cool) {
    const float r = sqrt(r2);
    *dU *= pow(r / r_cool, internal_energy_snowplow_exponent);
#ifdef SWIFT_DEBUG_CHECKS
    message("We do not resolve the Sedov-Taylor (r_cool = %e). Rescaling dU.",
            r_cool);
#endif /* SWIFT_DEBUG_CHECKS */
  } /* else we do not change dU */

  mechanical_feedback_accumulate_fluxes_for_conservation_check(si, dm, dp, m_ej,
                                                               p_ej, E_ej);
}

#elif FEEDBACK_GEAR_MECHANICAL_MODE == 2
/**
 @brief Feedback interaction between two particles (non-symmetric) for the
 * mechanical feedback mode 2.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param w_j_bar Normalized vector weight.
 * @param w_j_bar_norm Norm of the normalized vector weight.
 * @param v_i_p Physical velocity of particle i (star).
 * @param v_j_p Physical velocity of particle j (gas).
 * @param E_ej Physical ejected energy by the supernova.
 * @param m_ej Ejected mass by the supernova.
 * @param mj Mass of particle j.
 * @param dm Distributed mass from the SN to gas particle j.
 * @param new_mass Mass of particle j after the SN feedback.
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param constants The physical constants (in internal units).
 * @param us The internal system of units.
 * @param dU (return) Internal energy to distribute to gas particle j.
 * @param dKE (return) Kinetic energy variation before and after the SN
 * feedback.
 * @param dp_prime (return) Momemtum to distribute to gas particle j.

 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_2_feedback_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const double dm, const double new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp_prime[3]) {

  const double f_kin_0 = fb_props->f_kin_0;

  /* ... momentum */
  const double E_tot =
      E_ej + 0.5 * m_ej * si->feedback_data.accumulator.E_total;
  const double epsilon = f_kin_0 * E_tot; /* coupled kinetic energy */
  const double beta_1 =
      sqrt(m_ej / (2.0 * epsilon)) * si->feedback_data.accumulator.beta_1;
  const double beta_2 = m_ej * si->feedback_data.accumulator.beta_2;

  /* Compute the PdV work, taking into account gas in/outflows */
  const double psi = (sqrt(fabs(beta_2 + beta_1 * beta_1)) - beta_1) / beta_2;

  /* Now, we take into account for potentially unresolved energy-conserving
     phase of the SN explosion (xsi != 1). If we cannot resolve this phase, we
     give mostly momentum. The thermal energy will be radiated away because of
     cooling. */
  const double p_available = sqrt(2.0 * epsilon * m_ej);
  const double p_terminal = feedback_get_physical_SN_terminal_momentum(
      si, pj, xpj, phys_const, us, cosmo);
  const double xsi = min(1, p_terminal / (psi * p_available));

  /* Finally, the ejected velocity is */
  const double p_ej = psi * xsi * p_available;

  /* Now, we can compute dp */
  const double dp[3] = {w_j_bar[0] * p_ej, w_j_bar[1] * p_ej,
                        w_j_bar[2] * p_ej};

  /* Now boost to the 'laboratory' frame */
  dp_prime[0] = dp[0] + dm * v_i_p[0];
  dp_prime[1] = dp[1] + dm * v_i_p[1];
  dp_prime[2] = dp[2] + dm * v_i_p[2];

  /* ... physical internal energy */
  const double factor =
      (psi * psi * xsi * xsi) * beta_2 + 2.0 * (psi * xsi) * beta_1;
  const double f_therm = 1.0 - factor * epsilon / E_tot;
  const double U_tot = f_therm * E_tot;
  *dU = w_j_bar_norm * U_tot;

  /* Compute kinetic energy difference before and after SN */
  const double p_old_norm_2 =
      mj * mj *
      (v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2]);
  const double p_new[3] = {mj * v_j_p[0] + dp_prime[0],
                           mj * v_j_p[1] + dp_prime[1],
                           mj * v_j_p[2] + dp_prime[2]};
  const double p_new_norm_2 =
      p_new[0] * p_new[0] + p_new[1] * p_new[1] + p_new[2] * p_new[2];

  const double E_kin_old = p_old_norm_2 / (2.0 * mj);
  const double E_kin_new = p_new_norm_2 / (2.0 * new_mass);
  *dKE = E_kin_new - E_kin_old;

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  message(
      "beta_1 = %e, beta_2 = %e, psi = %e, psi*p_available = %e, p_available = "
      "%e",
      beta_1, beta_2, psi, psi * p_available, p_available);
  message("xsi = %e, p_t = %e", xsi, p_terminal);
  message(
      "E_ej = %e, E_tot = %e, U_tot = %e, E_kin_tot = %e, p_ej = %e, "
      "p_terminal = %e, dU = %e, f_therm = %e",
      E_ej, E_tot, U_tot, epsilon, p_ej, p_terminal, *dU, f_therm);
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */

  /* Now we accumulate to verify the conservation of the fluxes. */
  mechanical_feedback_accumulate_fluxes_for_conservation_check(si, dm, dp, m_ej,
                                                               p_ej, E_ej);
}

#else
#error "Mechanical feedback only supports mode 1 and 2"
#endif /* FEEDBACK_GEAR_MECHANICAL_MODE */

/**
 * @brief Mechanical feedback interaction between two particles (non-symmetric)
 * for supernovae.
 *
 * This function is used as a switch between the two modes to not pollute the
 * main code with too much ifdef.
 *
 *
 * @param r2 Comoving square distance between the two particles.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param w_j_bar Normalized vector weight.
 * @param w_j_bar_norm Norm of the normalized vector weight.
 * @param v_i_p Physical velocity of particle i (star).
 * @param v_j_p Physical velocity of particle j (gas).
 * @param E_ej Physical ejected energy by the supernova.
 * @param m_ej Ejected mass by the supernova.
 * @param mj Mass of particle j.
 * @param dm Distributed mass from the SN to gas particle j.
 * @param new_mass Mass of particle j after the SN feedback.
 * @param cosmo The cosmological model.
 * @param constants The physical constants (in internal units).
 * @param us The internal system of units.
 * @param dU (return) Internal energy to distribute to gas particle j.
 * @param dKE (return) Kinetic energy variation before and after the SN
 * feedback.
 * @param dp_prime (return) Momemtum to distribute to gas particle j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_feedback_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const double dm, const double new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp_prime[3]) {

#if FEEDBACK_GEAR_MECHANICAL_MODE == 1
  runner_iact_nonsym_mechanical_1_feedback_apply(
      r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj, dm,
      new_mass, cosmo, fb_props, phys_const, us, dU, dKE, dp_prime);
#elif FEEDBACK_GEAR_MECHANICAL_MODE == 2
  runner_iact_nonsym_mechanical_2_feedback_apply(
      r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj, dm,
      new_mass, cosmo, fb_props, phys_const, us, dU, dKE, dp_prime);
#else
#error "Mechanical feedback only supports mode 1 and 2"
#endif
}

#endif /* SWIFT_GEAR_MECHANICAL_FEEDBACK_IACT_H */
