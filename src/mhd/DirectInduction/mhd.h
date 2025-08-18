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
__attribute__((always_inline)) INLINE static double mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const double mu_0,
    const double a) {

  const double rho = p->rho;
  // convert to physical
  const double a_fact = pow(a, -3.0 * hydro_gamma);
  const double B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return 0.5 * a_fact * p->mass * B_over_rho2 * rho / mu_0;
}
/**
 * @brief Returns the magnetic field squared contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */

__attribute__((always_inline)) INLINE static double mhd_get_Bms(
    const struct part *p, const struct xpart *xp) {

  const double rho = p->rho;
  const double B_over_rho2 =
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

__attribute__((always_inline)) INLINE static double mhd_get_magnetic_divergence(
    const struct part *p, const struct xpart *xp) {

  return p->mhd_data.divB;
}

/**
 * @brief Returns the magnetic helicity contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static double mhd_get_magnetic_helicity(
    const struct part *p, const struct xpart *xp) {

  return 0.0;
}

__attribute__((always_inline)) INLINE static double mhd_get_cross_helicity(
    const struct part *p, const struct xpart *xp) {

  const double rho = p->rho;
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
__attribute__((always_inline)) INLINE static double mhd_get_divB_error(
    const struct part *p, const struct xpart *xp) {

  const double rho = p->rho;
  const double B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];

  const double error = B_over_rho2 != 0.0 ? fabs(p->mhd_data.divB) * p->h /
                                                sqrt(B_over_rho2 * rho * rho)
                                          : 0.0;

  return error;
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
__attribute__((always_inline)) INLINE static double mhd_compute_timestep(
    const struct part *p, const struct xpart *xp,
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo,
    const double mu_0) {

  const double dt_eta = p->mhd_data.resistive_eta != 0.0
                           ? hydro_properties->CFL_condition * cosmo->a *
                                 cosmo->a * p->h * p->h /
                                 p->mhd_data.resistive_eta
                           : DBL_MAX;

  return dt_eta;
}

/**
 * @brief Compute Alfven speed
 */
__attribute__((always_inline)) INLINE static double
mhd_get_comoving_Alfven_speed(const struct part *restrict p, const double mu_0) {

  /* Recover some data */
  const double rho = p->rho;
  double B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;

  /* B squared */
  const double B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

  /* Square of Alfven speed */
  const double vA2 = B2 / (mu_0 * rho);

  return sqrt(vA2);
}

/**
 * @brief Compute magnetosonic speed
 */
__attribute__((always_inline)) INLINE static double
mhd_get_comoving_magnetosonic_speed(const struct part *restrict p) {

  /* Compute square of fast magnetosonic speed */
  const double cs = hydro_get_comoving_soundspeed(p);
  const double cs2 = cs * cs;

  const double vA = p->mhd_data.Alfven_speed;
  const double vA2 = vA * vA;

  const double cms2 = cs2 + vA2;

  return sqrt(cms2);
}

/**
 * @brief Compute fast magnetosonic wave phase veolcity
 */
__attribute__((always_inline)) INLINE static double
mhd_get_comoving_fast_magnetosonic_wave_phase_velocity(
    const double dx[3], const struct part *restrict p, const double a,
    const double mu_0) {

  /* Get r and 1/r. */
  const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double rho = p->rho;
  double B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;

  /* B dot r. */
  const double Br = B[0] * dx[0] + B[1] * dx[1] + B[2] * dx[2];
  const double permeability_inv = 1.0 / mu_0;

  /* Compute effective sound speeds */
  const double cs = p->force.soundspeed;
  const double cs2 = cs * cs;
  const double c_ms = mhd_get_comoving_magnetosonic_speed(p);
  const double c_ms2 = c_ms * c_ms;
  const double projection_correction = c_ms2 * c_ms2 - 4.0 * permeability_inv *
                                                          cs2 * Br * r_inv *
                                                          Br * r_inv / rho;

  const double v_fmsw2 = 0.5 * (c_ms2 + sqrt(projection_correction));

  return sqrt(v_fmsw2);
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
__attribute__((always_inline)) INLINE static double mhd_signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double mu_ij, const double beta,
    const double a, const double mu_0) {

  const double v_sigi = mhd_get_comoving_magnetosonic_speed(pi);
  const double v_sigj = mhd_get_comoving_magnetosonic_speed(pj);

  const double v_sig = v_sigi + v_sigj - beta * mu_ij;

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
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const double mu_0) {

  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
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
  p->mhd_data.curl_B[0] = 0.0;
  p->mhd_data.curl_B[1] = 0.0;
  p->mhd_data.curl_B[2] = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->mhd_data.grad_B_tensor[i][j] = 0.0;
    }
  }

  /* SPH error*/
  p->mhd_data.mean_SPH_err = 0.0;
  for (int k = 0; k < 3; k++) {
    p->mhd_data.mean_grad_SPH_err[k] = 0.0;
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
  p->mhd_data.mean_SPH_err *= pow_dimension(1.0 / (p->h)) / p->rho;
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

  p->mhd_data.divB = 0.0;
  p->mhd_data.curl_B[0] = 0.0;
  p->mhd_data.curl_B[1] = 0.0;
  p->mhd_data.curl_B[2] = 0.0;
}

/**
 * @brief Calculate time derivative of dedner scalar
 *
 * @param p The particle to act upon
 * @param a The current value of the cosmological scale factor
 */
__attribute__((always_inline)) INLINE static double mhd_get_psi_over_ch_dt(
    struct part *p, const double a, const double a_factor_sound_speed,
    const double H, const struct hydro_props *hydro_props, const double mu_0) {

  /* Retrieve inverse of smoothing length. */
  const double h = p->h;
  const double h_inv = 1.0 / h;

  /* Compute Dedner cleaning speed. */
  const double ch = 0.5 * p->viscosity.v_sig;

  /* Compute Dedner cleaning scalar time derivative. */
  const double hyp = hydro_props->mhd.hyp_dedner;
  const double hyp_divv = hydro_props->mhd.hyp_dedner_divv;
  const double par = hydro_props->mhd.par_dedner;

  const double div_B = p->mhd_data.divB;
  const double div_v = a * a * hydro_get_div_v(p) - 3.0 * a * a * H;
  const double psi_over_ch = p->mhd_data.psi_over_ch;

  const double cp = ch;
  const double tau_inv = par * cp * h_inv;

  const double hyperbolic_term =
      -hyp * a * a * a_factor_sound_speed * a_factor_sound_speed * ch * div_B;
  const double hyperbolic_divv_term = -hyp_divv * psi_over_ch * div_v;
  const double parabolic_term =
      -a * a_factor_sound_speed * psi_over_ch * tau_inv;
  const double Hubble_term = a * a * H * psi_over_ch;

  return hyperbolic_term + hyperbolic_divv_term + parabolic_term + Hubble_term;
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
    const struct hydro_props *hydro_props, const double dt_alpha,
    const double mu_0) {

  const double h = p->h;
  const double rho = p->rho;
  double B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;

  const double B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
  const double normB = sqrt(B2);

  double grad_B_mean_square = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      grad_B_mean_square +=
          p->mhd_data.grad_B_tensor[i][j] * p->mhd_data.grad_B_tensor[i][j];
    }
  }

  const double alpha_AR_max = p->mhd_data.art_diff_beta;

  p->mhd_data.alpha_AR =
      normB ? fmin(alpha_AR_max, h * sqrt(grad_B_mean_square) / normB) : 0.0;
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
  p->mhd_data.divB = 0.0;

  p->mhd_data.B_over_rho_dt[0] = 0.0;
  p->mhd_data.B_over_rho_dt[1] = 0.0;
  p->mhd_data.B_over_rho_dt[2] = 0.0;

  p->mhd_data.B_over_rho_dt_AR[0] = 0.0;
  p->mhd_data.B_over_rho_dt_AR[1] = 0.0;
  p->mhd_data.B_over_rho_dt_AR[2] = 0.0;

  p->mhd_data.u_dt_AR = 0.0;

  /* Save forces*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.tot_mag_F[k] = 0.0;
  }
  /* Save induction sources*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.Adv_B_source[k] = 0.0;
    p->mhd_data.Diff_B_source[k] = 0.0;
    p->mhd_data.Delta_B[k] = 0.0;
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
    struct part *p, const struct xpart *xp, const struct cosmology *cosmo,
    const double mu_0) {

  /* Re-set the predicted magnetic flux densities */
  p->mhd_data.B_over_rho[0] = xp->mhd_data.B_over_rho_full[0];
  p->mhd_data.B_over_rho[1] = xp->mhd_data.B_over_rho_full[1];
  p->mhd_data.B_over_rho[2] = xp->mhd_data.B_over_rho_full[2];

  p->mhd_data.psi_over_ch = xp->mhd_data.psi_over_ch_full;

  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
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
    struct part *p, const struct xpart *xp, const double dt_drift,
    const double dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props, const double mu_0) {

  /* Predict the magnetic flux density */
  p->mhd_data.B_over_rho[0] += p->mhd_data.B_over_rho_dt[0] * dt_therm;
  p->mhd_data.B_over_rho[1] += p->mhd_data.B_over_rho_dt[1] * dt_therm;
  p->mhd_data.B_over_rho[2] += p->mhd_data.B_over_rho_dt[2] * dt_therm;

  p->mhd_data.psi_over_ch += p->mhd_data.psi_over_ch_dt * dt_therm;

  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
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
    const struct hydro_props *hydro_props, const double mu_0) {

  /* Get time derivative of Dedner scalar */
  p->mhd_data.psi_over_ch_dt = mhd_get_psi_over_ch_dt(
      p, cosmo->a, cosmo->a_factor_sound_speed, cosmo->H, hydro_props, mu_0);

  /* Hubble expansion contribution to induction equation */
  const double Hubble_induction_pref =
      cosmo->a * cosmo->a * cosmo->H * (1.5 * hydro_gamma - 2.0);
  p->mhd_data.B_over_rho_dt[0] +=
      Hubble_induction_pref * p->mhd_data.B_over_rho[0];
  p->mhd_data.B_over_rho_dt[1] +=
      Hubble_induction_pref * p->mhd_data.B_over_rho[1];
  p->mhd_data.B_over_rho_dt[2] +=
      Hubble_induction_pref * p->mhd_data.B_over_rho[2];

  /* Save forces*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.tot_mag_F[k] *= p->mass;
  }
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
    struct part *p, struct xpart *xp, const double dt_therm, const double dt_grav,
    const double dt_hydro, const double dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the magnetic flux density forward in time */
  const double delta_Bx = p->mhd_data.B_over_rho_dt[0] * dt_therm;
  const double delta_By = p->mhd_data.B_over_rho_dt[1] * dt_therm;
  const double delta_Bz = p->mhd_data.B_over_rho_dt[2] * dt_therm;

  /* Do not decrease the magnetic flux density by more than a factor of 2*/
  xp->mhd_data.B_over_rho_full[0] = xp->mhd_data.B_over_rho_full[0] + delta_Bx;
  xp->mhd_data.B_over_rho_full[1] = xp->mhd_data.B_over_rho_full[1] + delta_By;
  xp->mhd_data.B_over_rho_full[2] = xp->mhd_data.B_over_rho_full[2] + delta_Bz;

  xp->mhd_data.psi_over_ch_full += p->mhd_data.psi_over_ch_dt * dt_therm;
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
    const struct hydro_props *hydro_props, const double mu_0) {
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
  p->mhd_data.B_over_rho[0] *= pow(cosmo->a, 1.5 * hydro_gamma);
  p->mhd_data.B_over_rho[1] *= pow(cosmo->a, 1.5 * hydro_gamma);
  p->mhd_data.B_over_rho[2] *= pow(cosmo->a, 1.5 * hydro_gamma);

  /* Instantiate full step magnetic field */
  xp->mhd_data.B_over_rho_full[0] = p->mhd_data.B_over_rho[0];
  xp->mhd_data.B_over_rho_full[1] = p->mhd_data.B_over_rho[1];
  xp->mhd_data.B_over_rho_full[2] = p->mhd_data.B_over_rho[2];

  /* Instantiate full step magnetic Dedner scalar */
  xp->mhd_data.psi_over_ch_full = p->mhd_data.psi_over_ch;

  /* Instantiate Alfven speed */
  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
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
