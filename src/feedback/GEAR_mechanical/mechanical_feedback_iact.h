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

#define internal_energy_snowplow_exponent -6.5

__attribute__((always_inline)) INLINE static void
mechanical_feedback_accumulate_fluxes_for_conservation_check(
    struct spart *si, const double dm, const double dp[3], const float m_ej,
    const float p_ej, const float E_ej) {
  /* Reminder: This is not the SWIFT_DEBUG_CHECKS */
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
 @brief Stellar wind feedback interaction between two particles (non-symmetric)
 * for the mechanical feedback mode 2.
 *
 * Note that contrary to SN feedback, we do not assume winds do PdV work or end
 * with a terminal momentum.
 * Reference: https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.1578H
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
 * @param dp (return) Momemtum from mechanical feedback to distribute j.
 * @param dp_ejecta (return) Momemtum from ejecta to distribute to j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_1_stellar_winds_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const float dm, const float new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp[3],
    double dp_ejecta[3]) {

  const double new_mass_inv = 1.0 / new_mass;

  /* --------------- Compute physical momentum received ---------------- */
  /* Total momentum ejected by the winds during the timestep from the star
   * particle i */
  const float p_ej = sqrt(2.0 * m_ej * E_ej);
  double dp_prime[3] = {0.0};

  for (int i = 0; i < 3; i++) {
    /* Now, we can compute dp */
    dp[i] = w_j_bar[i] * p_ej;

    /* And the boost to the 'laboratory' frame */
    dp_ejecta[i] = dm * v_i_p[i];

    /* Gather all in a variable */
    dp_prime[i] = dp_ejecta[i] + dp[i];
  }
  /* Norm of physical velocities of the gas particle j */
  const float norm2_v_p =
      v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2];

  /* ----- Calculate physical Energy and internal Energy received ------ */
  const double norm2_delta_p_lab_frame =
      dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  const double norm2_delta_p = w_j_bar_norm * w_j_bar_norm * p_ej * p_ej;

  /* The energy ejected from the star particle i by stellar wind that is
   * actually received by the gas particle j */
  const double weighted_energy = w_j_bar_norm * E_ej;

  /* The additional energy received by the gas particle j due to the
   * momentum of the star particle i */
  const double dE_change_of_frame =
      0.5 * (norm2_delta_p_lab_frame - norm2_delta_p) / dm;

  /* The total energy received from the gas particle j in the laboratory
   * frame of reference */
  const double dE_lab_frame = weighted_energy + dE_change_of_frame;

  /* The momentum of the gas particle j after receiving the momentum from
   * stellar wind */
  const double p_new[3] = {mj * v_j_p[0] + dp[0], mj * v_j_p[1] + dp[1],
                           mj * v_j_p[2] + dp[2]};
  const double norm2_p_new = {p_new[0] * p_new[0] + p_new[1] * p_new[1] +
                              p_new[2] * p_new[2]};

  /* The new and old kinetic energy of the gas particle j */
  const double E_kin_old = 0.5 * mj * norm2_v_j_p;
  const double E_kin_new = 0.5 * norm2_p_new * new_mass_inv;

  /* The additional internal energy of the gas particle j.
     Ekin_new + U_new = Ekin_old + U_old + dEtot
     -> dU = (U_new - U_old) = (Ekin_old + dEtot - Ekin_new) */
  *dU = E_kin_old + dE_lab_frame - E_kin_new;

  /* In the frame of the gas particle. Note that v_j_p = 0 in this frame,
     which simplifies the formulas. */
  *dKE = 0.5 * dp_norm2 * new_mass_inv;

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  const double dE_kin = E_kin_new - E_kin_old;
  const double U_old = hydro_get_physical_internal_energy(pj, xpj, cosmo);
  message(
      "E_ej = %e, p_ej = %e | E_new = %e, U_new = %e, E_kin_new = %e, "
      "E_old = %e, U_old = %e, E_kin_old = %e | dE_prime = %e, dU = %e, "
      "dE_kin = %e, dKE_gas_frame = %e",
      E_ej, p_ej, E_new, U_new, E_kin_new, E_old, U_old, E_kin_old,
      dE_change_lab_frame, *dU, dE_kin, *dKE);
#endif /* SWIFT_FEEDBACK_DEBUG CHECKS */
}

/**
 @brief Supernovae feedback interaction between two particles (non-symmetric)
 * for the mechanical feedback mode 1.
 *
 * Reference: https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.1578H
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
 * @param dp (return) Momemtum from mechanical feedback to distribute j.
 * @param dp_ejecta (return) Momemtum from ejecta to distribute to j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_1_supernovae_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const float dm, const float new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp[3],
    double dp_ejecta[3]) {

  const float dm_inv = 1.0/dm;

  /* ... physical momentum */
  const double p_ej = sqrt(2 * m_ej * E_ej);
  double dp_prime[3] = {0.0};

  for (int i = 0; i < 3; i++) {
    /* Now, we can compute dp */
    dp[i] = w_j_bar[i] * p_ej;

    /* And the boost to the 'laboratory' frame */
    dp_ejecta[i] = dm * v_i_p[i];

    /* Gather all in a variable */
    dp_prime[i] = dp_ejecta[i] + dp[i];
  }

  /* Norm of physical velocities of the gas particle j */
  const float norm2_v_j_p =
      v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2];

  /* ... physical total energy */
  const double dp_norm_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  const double dp_prime_norm_2 = dp_prime[0] * dp_prime[0] +
                                 dp_prime[1] * dp_prime[1] +
                                 dp_prime[2] * dp_prime[2];
  const double dE = w_j_bar_norm * E_ej;
  const double dE_prime = dE + 0.5 * (dp_prime_norm_2 - dp_norm_2) *dm_inv;

  /* Now, we take into account for potentially unresolved energy-conserving
     phase of the SN explosion. If we cannot resolve this phase, we give mostly
     momentum. The thermal energy will be radiated away because of cooling. */
  const double PdV_work_fraction = sqrt(1 + mj * dm_inv);
  const double p_terminal = feedback_get_physical_SN_terminal_momentum(
      si, pj, xpj, phys_const, us, fb_props, cosmo);

  /* If we can resolve the Taylor Sedov, then we give the right coupled
     momentum (which is by definition <= p_terminal). If we cannot resolve it,
     then we have reached the p_terminal. This is the upper limit of momentum,
     since afterwards the cooling is efficient and the thermal energy is
     radiated away. Thus, the factor to multiply dp_prime is: */
  const double p_factor = min(PdV_work_fraction, p_terminal / p_ej);

  for (int i = 0; i < 3; i++) {
    /* Now, we can compute dp */
    dp[i] *= p_factor;

    /* And the boost to the 'laboratory' frame */
    dp_ejecta[i] *= p_factor;

    /* Gather all in a variable */
    dp_prime[i] *= p_factor;
  }

  /* Note: Changing dp_prime after computing the total energy is correct. The
     total energy does not specify how much energy goes into kinetic and
     internal form. The p_factor changes the ratio: the energy that goes into
     the PdV work (thus increasing the momentum and kinetic energy) is
     automatically transformed from thermal to kinetic energy.
     However, to compute the total kinetic energy, we need to use the updated
     dp_prime. */

  /* ... physical internal energy */
  /* Compute kinetic energy difference before and after SN */
  const double p_new[3] = {mj * v_j_p[0] + dp_prime[0],
                           mj * v_j_p[1] + dp_prime[1],
                           mj * v_j_p[2] + dp_prime[2]};
  const double p_new_norm_2 =
      p_new[0] * p_new[0] + p_new[1] * p_new[1] + p_new[2] * p_new[2];

  const double E_kin_old = 0.5 * norm2_v_j_p * mj;
  const double E_kin_new = 0.5 * p_new_norm_2 / new_mass;

  const double U_old = hydro_get_physical_internal_energy(pj, xpj, cosmo);
  const double E_old = U_old + E_kin_old;
  const double E_new = E_old + dE_prime;
  const double U_new = E_new - E_kin_new;

  /* Compute the physical internal energy */
  *dU = U_new - U_old;

  /* In the frame of the gas particle. Note that v_j_p = 0 in this frame,
     which simplifies the formulas. */
  *dKE = 0.5 * dp_norm_2 * new_mass_inv;

  /* Compute the comoving cooling radius */
  const float r_cool = cosmo->a_inv * feedback_get_physical_SN_cooling_radius(
                                          si, p_ej, p_terminal, cosmo);
  const float r_cool_2 = r_cool * r_cool;

  /* If we do not resolve the Taylor-Sedov, we rescale the internal energy */
  if (r2 > r_cool_2) {
    const float r = sqrt(r2);
    *dU *= pow(r / r_cool, internal_energy_snowplow_exponent);
#ifdef SWIFT_DEBUG_CHECKS
    message("We do not resolve the Sedov-Taylor (r_cool = %e). Rescaling dU.",
            r_cool);
#endif /* SWIFT_DEBUG_CHECKS */
  } /* else we do not change dU */

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  message(
      "E_ej = %e, p_ej = %e | p_terminal = %e, p_factor = %e | E_new = %e, "
      "U_new = %e, E_kin_new = %e, E_old = %e, U_old = %e, E_kin_old = %e "
      "| dE_prime = %e, dU = %e",
      E_ej, p_ej, p_terminal, p_factor, E_new, U_new, E_kin_new, E_old, U_old,
      E_kin_old, dE_prime, *dU);
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */

  mechanical_feedback_accumulate_fluxes_for_conservation_check(si, dm, dp, m_ej,
                                                               p_ej, E_ej);
}

#elif FEEDBACK_GEAR_MECHANICAL_MODE == 2

/**
 @brief Stellar wind feedback interaction between two particles (non-symmetric)
 * for the mechanical feedback mode 2.
 *
 * Note that contrary to SN feedback, we do not assume winds do PdV work or end
 * with a terminal momentum. Hence, we have psi = xsi = 1.0. The coupled
 * kinetic energy is thus assumed to be the wind energy:
 *               epsilon = E_ej = 0.5 * m_winds * v_winds^2.
 * From these assumptions, we can derive the internal energy using the
 * equations in appendix A from
 * https://ui.adsabs.harvard.edu/abs/2024arXiv240416987H
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
 * @param dp (return) Momemtum from mechanical feedback to distribute j.
 * @param dp_ejecta (return) Momemtum from ejecta to distribute to j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_2_stellar_winds_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const float dm, const float new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp[3],
    double dp_ejecta[3]) {

  const double new_mass_inv = 1.0 / new_mass;

  /* Compute the relevant variables from the accumulators. */
  const double E_tot =
      E_ej + 0.5 * m_ej * si->feedback_data.accumulator_winds.E_total;
  const double epsilon = E_ej; /* coupled kinetic energy */
  const double beta_1 =
      sqrt(0.5 * m_ej / epsilon) * si->feedback_data.accumulator_winds.beta_1;
  const double beta_2 = m_ej * si->feedback_data.accumulator_winds.beta_2;

  /* Winds do not assume any PdV work or terminal momentum. Therefore,
     the ejected velocity is */
  const double p_ej = sqrt(2.0 * epsilon * m_ej);

  for (int i = 0; i < 3; i++) {
    /* Now, we can compute dp */
    dp[i] = w_j_bar[i] * p_ej;

    /* And the boost to the 'laboratory' frame */
    dp_ejecta[i] = dm * v_i_p[i];
  }

  /* ... physical internal energy */
  const double U_tot = E_tot - epsilon * (beta_2 + 2.0 * beta_1);
  *dU = w_j_bar_norm * U_tot;

  /* In the frame of the gas particle. Note that v_j_p = 0 in this frame,
     which simplifies the formulas. */
  const double dp_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  *dKE = 0.5 * dp_2 * new_mass_inv;

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  /* Compute kinetic energy difference before and after SN */
  const double mj_2 = mj * mj;
  const double mj_inv = 1.0 / mj;
  const double p_old_norm_2 =
      mj_2 * (v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2]);
  const double p_new[3] = {mj * v_j_p[0] + dp_ejecta[0] + dp[0],
                           mj * v_j_p[1] + dp_ejecta[1] + dp[1],
                           mj * v_j_p[2] + dp_ejecta[2] + dp[2]};
  const double p_new_norm_2 =
      p_new[0] * p_new[0] + p_new[1] * p_new[1] + p_new[2] * p_new[2];

  const double E_kin_old = 0.5 * p_old_norm_2 * mj_inv;
  const double E_kin_new = 0.5 * p_new_norm_2 * new_mass_inv;
  const double dE_kin = E_kin_new - E_kin_old;
  message(
      "beta_1 = %e, beta_2 = %e | E_ej = %e, p_ej = %e | E_tot = %e, U_tot = "
      "%e, E_kin_tot = %e | E_kin_old = %e, E_kin_new = %e | dU = %e, "
      "f_therm = %e, dE_kin = %e, dKE_gas_frame = %e",
      beta_1, beta_2, E_ej, p_ej, E_tot, U_tot, epsilon, E_kin_old, E_kin_new,
      *dU, f_therm, dE_kin, *dKE);
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */
}

/**
 @brief Supernovae feedback interaction between two particles (non-symmetric)
 * for the mechanical feedback mode 2.
 *
 * Reference: Appendix A in
 * https://ui.adsabs.harvard.edu/abs/2024arXiv240416987H
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
 * @param dp (return) Momemtum from mechanical feedback to distribute j.
 * @param dp_ejecta (return) Momemtum from ejecta to distribute to j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_2_supernovae_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const float dm, const float new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp[3],
    double dp_ejecta[3]) {

  const double new_mass_inv = 1.0 / new_mass;
  const double f_kin_0 = fb_props->f_kin_0;

  /* Compute the relevant variables from the accumulators. */
  const double E_tot =
      E_ej + 0.5 * m_ej * si->feedback_data.accumulator_sn.E_total;
  const double epsilon = f_kin_0 * E_tot; /* coupled kinetic energy */
  const double beta_1 =
      sqrt(0.5 * m_ej / epsilon) * si->feedback_data.accumulator_sn.beta_1;
  const double beta_2 = m_ej * si->feedback_data.accumulator_sn.beta_2;

  /* Compute the PdV work, taking into account gas in/outflows */
  const double psi = (sqrt(fabs(beta_2 + beta_1 * beta_1)) - beta_1) / beta_2;

  /* Now, we take into account for potentially unresolved energy-conserving
     phase of the SN explosion (xsi != 1). If we cannot resolve this phase, we
     give mostly momentum. The thermal energy will be radiated away because of
     cooling. */
  const double p_available = sqrt(2.0 * epsilon * m_ej);
  const double p_terminal = feedback_get_physical_SN_terminal_momentum(
      si, pj, xpj, phys_const, us, fb_props, cosmo);
  const double xsi = min(1.0, p_terminal / (psi * p_available));

  /* Finally, the ejected velocity is */
  const double p_ej = psi * xsi * p_available;

  for (int i = 0; i < 3; i++) {
    /* Now, we can compute dp */
    dp[i] = w_j_bar[i] * p_ej;

    /* And the boost to the 'laboratory' frame */
    dp_ejecta[i] = dm * v_i_p[i];
  }

  /* ... physical internal energy */
  const double factor =
      (psi * psi * xsi * xsi) * beta_2 + 2.0 * (psi * xsi) * beta_1;
  const double f_therm = 1.0 - factor * epsilon / E_tot;
  const double U_tot = f_therm * E_tot;
  *dU = w_j_bar_norm * U_tot;

  /* In the frame of the gas particle. Note that v_j_p = 0 in this frame,
     which simplifies the formulas. */
  const double dp_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  *dKE = 0.5 * dp_2 * new_mass_inv;

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  /* Compute kinetic energy difference before and after SN */
  const double mj_2 = mj * mj;
  const double mj_inv = 1.0 / mj;
  const double p_old_norm_2 =
      mj_2 * (v_j_p[0] * v_j_p[0] + v_j_p[1] * v_j_p[1] + v_j_p[2] * v_j_p[2]);
  const double p_new[3] = {mj * v_j_p[0] + dp_ejecta[0] + dp[0],
                           mj * v_j_p[1] + dp_ejecta[1] + dp[1],
                           mj * v_j_p[2] + dp_ejecta[2] + dp[2]};
  const double p_new_norm_2 =
      p_new[0] * p_new[0] + p_new[1] * p_new[1] + p_new[2] * p_new[2];

  const double E_kin_old = 0.5 * p_old_norm_2 * mj_inv;
  const double E_kin_new = 0.5 * p_new_norm_2 * new_mass_inv;
  const double dE_kin = E_kin_new - E_kin_old;

  message(
      "beta_1 = %e, beta_2 = %e, psi = %e, psi*p_available = %e, p_available = "
      "%e | xsi = %e, p_t = %e | E_ej = %e, p_ej = %e | E_tot = %e, U_tot = %e"
      ", E_kin_tot = %e | E_kin_old = %e, E_kin_new = %e | p_terminal = %e, "
      "dU = %e, f_therm = %e, dE_kin = %e, dKE_gas_frame = %e",
      beta_1, beta_2, psi, psi * p_available, p_available, xsi, p_terminal,
      E_ej, p_ej, E_tot, U_tot, epsilon, E_kin_old, E_kin_new, p_terminal, *dU,
      f_therm, dE_kin, *dKE);
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
 * @param dp (return) Momemtum from mechanical feedback to distribute j.
 * @param dp_ejecta (return) Momemtum from ejecta to distribute to j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_stellar_winds_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const float dm, const float new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp[3],
    double dp_ejecta[3]) {

#if FEEDBACK_GEAR_MECHANICAL_MODE == 1
  runner_iact_nonsym_mechanical_1_stellar_winds_apply(
      r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj, dm,
      new_mass, cosmo, fb_props, phys_const, us, dU, dKE, dp, dp_ejecta);
#elif FEEDBACK_GEAR_MECHANICAL_MODE == 2
  runner_iact_nonsym_mechanical_2_stellar_winds_apply(
      r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj, dm,
      new_mass, cosmo, fb_props, phys_const, us, dU, dKE, dp, dp_ejecta);
#else
#error "Mechanical feedback only supports mode 1 and 2"
#endif
}

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
 * @param dp (return) Momemtum from mechanical feedback to distribute j.
 * @param dp_ejecta (return) Momemtum from ejecta to distribute to j.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mechanical_feedback_apply(
    const float r2, struct spart *si, struct part *pj, struct xpart *xpj,
    const double w_j_bar[3], const double w_j_bar_norm, const float v_i_p[3],
    const float v_j_p[3], const float E_ej, const float m_ej, const float mj,
    const float dm, const float new_mass, const struct cosmology *cosmo,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, double *dU, double *dKE, double dp[3],
    double dp_ejecta[3]) {

#if FEEDBACK_GEAR_MECHANICAL_MODE == 1
  runner_iact_nonsym_mechanical_1_supernovae_apply(
      r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj, dm,
      new_mass, cosmo, fb_props, phys_const, us, dU, dKE, dp, dp_ejecta);
#elif FEEDBACK_GEAR_MECHANICAL_MODE == 2
  runner_iact_nonsym_mechanical_2_supernovae_apply(
      r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj, dm,
      new_mass, cosmo, fb_props, phys_const, us, dU, dKE, dp, dp_ejecta);
#else
#error "Mechanical feedback only supports mode 1 and 2"
#endif
}

#endif /* SWIFT_GEAR_MECHANICAL_FEEDBACK_IACT_H */
