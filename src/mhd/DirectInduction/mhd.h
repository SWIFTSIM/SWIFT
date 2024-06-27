/*******************************************************************************
 * This file is part o SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DIRECT_INDUCTION_MHD_H
#define SWIFT_DIRECT_INDUCTION_MHD_H
#include "minmax.h"

#include <float.h>

/**
 * @brief Returns the magnetic energy contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const float mu_0) {

  const float rho = p->rho;
  const float B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return 0.5f * p->mass * B_over_rho2 * rho / mu_0;
}
/**
 * @brief Returns the magnetic field squared contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */

__attribute__((always_inline)) INLINE static float mhd_get_Bms(
    const struct part *p, const struct xpart *xp) {

  const float rho = p->rho;
  const float B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return B_over_rho2 * rho * rho;
}
/**
 * @brief Returns the magnetic field divergence of a particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */

__attribute__((always_inline)) INLINE static float mhd_get_magnetic_divergence(
    const struct part *p, const struct xpart *xp) {

  return p->mhd_data.divB;
}

/**
 * @brief Returns the magnetic helicity contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_helicity(
    const struct part *p, const struct xpart *xp) {

  return 0.f;
}

__attribute__((always_inline)) INLINE static float mhd_get_cross_helicity(
    const struct part *p, const struct xpart *xp) {

  const float rho = p->rho;
  return (p->v[0] * p->mhd_data.B_over_rho[0] +
          p->v[1] * p->mhd_data.B_over_rho[1] +
          p->v[2] * p->mhd_data.B_over_rho[2]) *
         rho;
}

/**
 * @brief Returns the magnetic field divergence error of the particle.
 *
 * This is (div B) / (B / h) and is hence dimensionless.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_divB_error(
    const struct part *p, const struct xpart *xp) {

  const float rho = p->rho;
  const float B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return fabs(p->mhd_data.divB) * p->h /
         (sqrtf(B_over_rho2 * rho * rho) + 1.e-18);
}

/**
 * @brief Computes the MHD time-step of a given particle
 *
 * This function returns the time-step of a particle given its hydro-dynamical
 * state. A typical time-step calculation would be the use of the CFL condition.
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The SPH parameters
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float mhd_compute_timestep(
    const struct part *p, const struct xpart *xp,
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo,
    const float mu_0) {

  /* Dt from 1/DivOperator(Alfven speed) */

  const float divB = p->mhd_data.divB;

  const float dt_B_factor = fabsf(divB);

  const float dt_B_derivatives =
      dt_B_factor != 0.f
          ? hydro_properties->CFL_condition * cosmo->a /
                cosmo->a_factor_sound_speed *
                sqrtf(p->rho * mu_0 / (dt_B_factor * dt_B_factor))
          : FLT_MAX;

  const float dt_eta = p->mhd_data.resistive_eta != 0.f
                           ? hydro_properties->CFL_condition * cosmo->a *
                                 cosmo->a * p->h * p->h /
                                 p->mhd_data.resistive_eta
                           : FLT_MAX;

  return fminf(dt_B_derivatives, dt_eta);
}

/**
 * @brief Compute magnetosonic speed
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetosonic_speed(
    const struct part *restrict p, const float a, const float mu_0) {

  /* Recover some data */
  const float rho = p->rho;
  float B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;

  /* B squared */
  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

  const float permeability_inv = 1 / mu_0;

  /* Compute effective sound speeds */
  const float cs = p->force.soundspeed;
  const float cs2 = cs * cs;
  const float v_A2 = permeability_inv * B2 / rho;
  const float c_ms2 = cs2 + v_A2;

  return sqrtf(c_ms2);
}

/**
 * @brief Compute fast magnetosonic wave phase veolcity
 */
__attribute__((always_inline)) INLINE static float
mhd_get_fast_magnetosonic_wave_phase_velocity(const float dx[3],
                                              const struct part *restrict p,
                                              const float a, const float mu_0) {

  /* Get r and 1/r. */
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float rho = p->rho;
  float B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;

  /* B dot r. */
  const float Br = B[0] * dx[0] + B[1] * dx[1] + B[2] * dx[2];
  const float permeability_inv = 1 / mu_0;

  /* Compute effective sound speeds */
  const float cs = p->force.soundspeed;
  const float cs2 = cs * cs;
  const float c_ms = mhd_get_magnetosonic_speed(p, a, mu_0);
  const float c_ms2 = c_ms * c_ms;
  const float projection_correction = c_ms2 * c_ms2 - 4.0f * permeability_inv *
                                                          cs2 * Br * r_inv *
                                                          Br * r_inv / rho;

  const float v_fmsw2 = 0.5f * (c_ms2 + sqrtf(projection_correction));

  return sqrtf(v_fmsw2);
}

/**
 * @brief Compute the MHD signal velocity between two gas particles
 *
 * This is eq. (131) of Price D., JCoPh, 2012, Vol. 231, Issue 3
 *
 * Warning ONLY to be called just after preparation of the force loop.
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float mhd_signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta,
    const float a, const float mu_0) {

  const float v_sigi = mhd_get_magnetosonic_speed(pi, a, mu_0);
  const float v_sigj = mhd_get_magnetosonic_speed(pj, a, mu_0);

  /*
  const float v_sigi =
      mhd_get_fast_magnetosonic_wave_phase_velocity(dx, pi, a, mu_0);
  const float v_sigj =
      mhd_get_fast_magnetosonic_wave_phase_velocity(dx, pj, a, mu_0);
  */

  const float v_sig = v_sigi + v_sigj - beta * mu_ij;

  return v_sig;
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density loop over neighbours. Typically, all fields of the
 * density sub-structure of a particle get zeroed in here.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_init_part(
    struct part *p) {}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 * Additional quantities such as velocity gradients will also get the final
 * terms added to them here.
 *
 * Also adds/multiplies the cosmological terms if need be.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_end_density(
    struct part *p, const struct cosmology *cosmo) {}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void mhd_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props) {

  p->force.balsara = 1.f;
}

/**
 * @brief Resets the variables that are required for a gradient calculation.
 *
 * This function is called after mhd_prepare_gradient.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_reset_gradient(
    struct part *p) {

  /* Zero the fields updated by the mhd gradient loop */
  p->mhd_data.divB = 0.0f;
  p->mhd_data.curl_B[0] = 0.0f;
  p->mhd_data.curl_B[1] = 0.0f;
  p->mhd_data.curl_B[2] = 0.0f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->mhd_data.grad_B_tensor[i][j] = 0.0f;
    }
  }

  /* SPH error*/
  p->mhd_data.mean_SPH_err = 0.f;
  for (int k = 0; k < 3; k++) {
    p->mhd_data.mean_grad_SPH_err[k] = 0.f;
  }
}

/**
 * @brief Finishes the gradient calculation.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void mhd_end_gradient(
    struct part *p) {
  /* Add self contribution */
  p->mhd_data.mean_SPH_err += p->mass * kernel_root;
  /* Finish SPH_1 calculation*/
  p->mhd_data.mean_SPH_err *= pow_dimension(1.f / (p->h)) / p->rho;
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * In the desperate case where a particle has no neighbours (likely because
 * of the h_max ceiling), set the particle fields to something sensible to avoid
 * NaNs in the next calculations.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_part_has_no_neighbours(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {

  p->mhd_data.divB = 0.0f;
  p->mhd_data.curl_B[0] = 0.0f;
  p->mhd_data.curl_B[1] = 0.0f;
  p->mhd_data.curl_B[2] = 0.0f;
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * This function is called in the ghost task to convert some quantities coming
 * from the density loop over neighbours into quantities ready to be used in the
 * force loop over neighbours. Quantities are typically read from the density
 * sub-structure and written to the force sub-structure.
 * Examples of calculations done here include the calculation of viscosity term
 * constants, thermal conduction terms, hydro conversions, etc.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param hydro_props Hydrodynamic properties.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 */
__attribute__((always_inline)) INLINE static void mhd_prepare_force(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float dt_alpha,
    const float mu_0) {

  const float h = p->h;
  const float rho = p->rho;
  float B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;

  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
  const float normB = sqrtf(B2);

  float grad_B_mean_square = 0.0f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      grad_B_mean_square +=
          p->mhd_data.grad_B_tensor[i][j] * p->mhd_data.grad_B_tensor[i][j];
    }
  }

  p->mhd_data.alpha_AR =
      normB ? fminf(1.0f, h * sqrtf(grad_B_mean_square) / normB) : 0.0f;
}

/**
 * @brief Calculate time derivative of dedner scalar
 *
 * @param p The particle to act upon
 * @param a The current value of the cosmological scale factor
 */
__attribute__((always_inline)) INLINE static float mhd_get_psi_over_ch_dt(
    struct part *p, const float a, const float a_factor_sound_speed,
    const float H, const struct hydro_props *hydro_props, const float mu_0) {

  /* Retrieve inverse of smoothing length. */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  /* Compute Dedner cleaning speed. */
  const float ch = mhd_get_magnetosonic_speed(p, a, mu_0);

  /* Compute Dedner cleaning scalar time derivative. */
  const float hyp = hydro_props->mhd.hyp_dedner;
  const float hyp_divv = hydro_props->mhd.hyp_dedner_divv;
  const float par = hydro_props->mhd.par_dedner;

  const float div_B = p->mhd_data.divB;
  const float div_v = a * a * hydro_get_div_v(p) - 3.0f * a * a * H;
  const float psi_over_ch = p->mhd_data.psi_over_ch;

  const float hyperbolic_term =
      -hyp * a * a * a_factor_sound_speed * a_factor_sound_speed * ch * div_B;
  const float hyperbolic_divv_term = -hyp_divv * psi_over_ch * div_v;
  const float parabolic_term =
      -par * a * a_factor_sound_speed * psi_over_ch * ch * h_inv;
  const float Hubble_term = a * a * H * psi_over_ch;

  return hyperbolic_term + hyperbolic_divv_term + parabolic_term + Hubble_term;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_reset_acceleration(
    struct part *restrict p) {

  /* Zero the fields updated by the mhd force loop */
  p->mhd_data.B_over_rho_dt[0] = 0.0f;
  p->mhd_data.B_over_rho_dt[1] = 0.0f;
  p->mhd_data.B_over_rho_dt[2] = 0.0f;

  p->mhd_data.B_over_rho_dt_AR[0] = 0.0f;
  p->mhd_data.B_over_rho_dt_AR[1] = 0.0f;
  p->mhd_data.B_over_rho_dt_AR[2] = 0.0f;

  p->mhd_data.u_dt_AR = 0.0f;
  
  /* Save forces*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.tot_mag_F[k] = 0.0f;
  }
  /* Save induction sources*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.Adv_B_source[k] = 0.0f;
    p->mhd_data.Diff_B_source[k] = 0.0f;
    p->mhd_data.Delta_B[k] = 0.0f;
  }
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model
 */
__attribute__((always_inline)) INLINE static void mhd_reset_predicted_values(
    struct part *p, const struct xpart *xp, const struct cosmology *cosmo) {

  /* Re-set the predicted magnetic flux densities */
  p->mhd_data.B_over_rho[0] = xp->mhd_data.B_over_rho_full[0];
  p->mhd_data.B_over_rho[1] = xp->mhd_data.B_over_rho_full[1];
  p->mhd_data.B_over_rho[2] = xp->mhd_data.B_over_rho_full[2];
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * Note the different time-step sizes used for the different quantities as they
 * include cosmological factors.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 * @param mu_0 The vacuum magnetic permeability.
 */
__attribute__((always_inline)) INLINE static void mhd_predict_extra(
    struct part *p, const struct xpart *xp, const float dt_drift,
    const float dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props, const float mu_0) {

  /* Predict the magnetic flux density */
  p->mhd_data.B_over_rho[0] += p->mhd_data.B_over_rho_dt[0] * dt_therm;
  p->mhd_data.B_over_rho[1] += p->mhd_data.B_over_rho_dt[1] * dt_therm;
  p->mhd_data.B_over_rho[2] += p->mhd_data.B_over_rho_dt[2] * dt_therm;

  p->mhd_data.psi_over_ch_dt = mhd_get_psi_over_ch_dt(
      p, cosmo->a, cosmo->a_factor_sound_speed, cosmo->H, hydro_props, mu_0);
  p->mhd_data.psi_over_ch +=
      mhd_get_psi_over_ch_dt(p, cosmo->a, cosmo->a_factor_sound_speed, cosmo->H,
                             hydro_props, mu_0) *
      dt_therm;
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the force and accelerations by the appropiate constants
 * and add the self-contribution term. In most cases, there is little
 * to do here.
 *
 * Cosmological terms are also added/multiplied here.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_end_force(
    struct part *p, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float mu_0) {

  /* Hubble expansion contribution to induction equation */

  const float Hubble_induction_pref =
      cosmo->a * cosmo->a * cosmo->H * (1.5f * hydro_gamma - 2.f);
  p->mhd_data.B_over_rho_dt[0] +=
      Hubble_induction_pref * p->mhd_data.B_over_rho[0];
  p->mhd_data.B_over_rho_dt[1] +=
      Hubble_induction_pref * p->mhd_data.B_over_rho[1];
  p->mhd_data.B_over_rho_dt[2] +=
      Hubble_induction_pref * p->mhd_data.B_over_rho[2];
}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantities are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void mhd_kick_extra(
    struct part *p, struct xpart *xp, const float dt_therm, const float dt_grav,
    const float dt_hydro, const float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the magnetic flux density forward in time */
  const float delta_Bx = p->mhd_data.B_over_rho_dt[0] * dt_therm;
  const float delta_By = p->mhd_data.B_over_rho_dt[1] * dt_therm;
  const float delta_Bz = p->mhd_data.B_over_rho_dt[2] * dt_therm;

  /* Do not decrease the magnetic flux density by more than a factor of 2*/
  xp->mhd_data.B_over_rho_full[0] = xp->mhd_data.B_over_rho_full[0] + delta_Bx;
  xp->mhd_data.B_over_rho_full[1] = xp->mhd_data.B_over_rho_full[1] + delta_By;
  xp->mhd_data.B_over_rho_full[2] = xp->mhd_data.B_over_rho_full[2] + delta_Bz;
}

/**
 * @brief Converts MHD quantities of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 * This can be used to convert internal energy into entropy in the case
 * of hydro for instance.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 */
__attribute__((always_inline)) INLINE static void mhd_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props) {
  /* Set Restitivity Eta */
  p->mhd_data.resistive_eta = hydro_props->mhd.mhd_eta;
  /* Set Monopole subtraction factor */
  p->mhd_data.monopole_beta = hydro_props->mhd.monopole_subtraction;
  /* Set Artificial Difussion */
  p->mhd_data.art_diff_beta = hydro_props->mhd.art_diffusion;

  /* Convert B into B/rho */
  p->mhd_data.B_over_rho[0] /= p->rho;
  p->mhd_data.B_over_rho[1] /= p->rho;
  p->mhd_data.B_over_rho[2] /= p->rho;

  /* Convert to co-moving B/rho */
  p->mhd_data.B_over_rho[0] *= pow(cosmo->a, 1.5f * hydro_gamma);
  p->mhd_data.B_over_rho[1] *= pow(cosmo->a, 1.5f * hydro_gamma);
  p->mhd_data.B_over_rho[2] *= pow(cosmo->a, 1.5f * hydro_gamma);

  xp->mhd_data.B_over_rho_full[0] = p->mhd_data.B_over_rho[0];
  xp->mhd_data.B_over_rho_full[1] = p->mhd_data.B_over_rho[1];
  xp->mhd_data.B_over_rho_full[2] = p->mhd_data.B_over_rho[2];
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_first_init_part(
    struct part *restrict p, struct xpart *restrict xp,
    const struct mhd_global_data *mhd_data, const double Lsize) {

  mhd_reset_acceleration(p);
  mhd_init_part(p);
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_H */
