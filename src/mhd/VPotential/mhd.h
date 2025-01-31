/*******************************************************************************
 * This file is part of SWIFT.
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
#ifndef SWIFT_VECTOR_POTENTIAL_MHD_H
#define SWIFT_VECTOR_POTENTIAL_MHD_H

#include "hydro.h"
#include "mhd_parameters.h"
#include "space.h"

#include <float.h>

/**
 * @brief Returns the magnetic energy contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const float mu_0,
    const float a) {

  const float afact = pow(a, 2.f * mhd_comoving_factor);

  const float b2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  return 0.5f * afact * b2 / mu_0 * p->mass / p->rho;
}
/**
 * @brief Returns the magnetic field squared contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */

__attribute__((always_inline)) INLINE static float mhd_get_Bms(
    const struct part *p, const struct xpart *xp) {

  const float b2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  return b2;
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

  return p->mhd_data.APred[0] * p->mhd_data.BPred[0] +
         p->mhd_data.APred[1] * p->mhd_data.BPred[1] +
         p->mhd_data.APred[2] * p->mhd_data.BPred[2];
}

__attribute__((always_inline)) INLINE static float mhd_get_cross_helicity(
    const struct part *p, const struct xpart *xp) {

  return p->v[0] * p->mhd_data.BPred[0] + p->v[1] * p->mhd_data.BPred[1] +
         p->v[2] * p->mhd_data.BPred[2];
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

  const float b2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  float error = b2 != 0.f ? sqrt(p->mhd_data.divB * p->mhd_data.divB / b2) *
  			p->h : 0.0f;

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
__attribute__((always_inline)) INLINE static float mhd_compute_timestep(
    const struct part *p, const struct xpart *xp,
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo,
    const float mu_0) {

  //const float afac_divB = pow(cosmo->a, -mhd_comoving_factor - 0.5f); //(cosmo->a / cosmo->a_factor_sound_speed);
  const float afac_divB = (cosmo->a / cosmo->a_factor_sound_speed);
  const float afac_resistive = cosmo->a * cosmo->a;
  /* Dt from 1/DivOperator(Alfven speed) */

  float dt_divA =
      p->mhd_data.divA != 0.0f
          ? afac_divB * hydro_properties->CFL_condition * p->h *
                sqrtf(p->rho / (p->mhd_data.divA * p->mhd_data.divA) * mu_0)
          : FLT_MAX;
  float dt_divB =
      p->mhd_data.divB != 0.0f
          ? afac_divB * hydro_properties->CFL_condition *
                sqrtf(p->rho / (p->mhd_data.divB * p->mhd_data.divB) * mu_0)
          : FLT_MAX;
  dt_divB = fminf(dt_divA, dt_divB);
  const float resistive_eta = p->mhd_data.resistive_eta;
  const float dt_eta = resistive_eta != 0.0f
                           ? afac_resistive * hydro_properties->CFL_condition *
                                 p->h * p->h / resistive_eta * 0.5
                           : FLT_MAX;

  return fminf(dt_eta, dt_divB);
}

/**
 * @brief Compute magnetosonic speed
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetosonic_speed(
    const struct part *restrict p, const float a, const float mu_0) {

  /* Recover some data */
  const float rho = p->rho;

  /* B squared */
  const float B2 = (p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                    p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                    p->mhd_data.BPred[2] * p->mhd_data.BPred[2]);

  const float permeability_inv = 1.0f / mu_0;

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
  B[0] = p->mhd_data.BPred[0];
  B[1] = p->mhd_data.BPred[1];
  B[2] = p->mhd_data.BPred[2];

  /* B dot r. */
  const float Br = B[0] * dx[0] + B[1] * dx[1] + B[2] * dx[2];
  const float permeability_inv = 1.0f / mu_0;

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
 * @brief Compute the MHD signal velocity between two gas particles,
 *
 * This is eq. (131) of Price D., JCoPh, 2012, Vol. 231, Issue 3.
 *
 * Warning ONLY to be called just after preparation of the force loop.
 * Warning BPred is used XXX/TODO
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
    struct part *p) 
 
    {

      p->mhd_data.divA = 0.f;
      // XXX todo, really is not the predicted variable will be the full step
      p->mhd_data.BPred[0] = 0.f;
      p->mhd_data.BPred[1] = 0.f;
      p->mhd_data.BPred[2] = 0.f;
  }
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
    struct part *p, const struct cosmology *cosmo) {

  const float h_inv_dim_plus_one =
      pow_dimension_plus_one(1.f / p->h); /*1/h^(d+1) */
  const float rho_inv = 1.f / p->rho;
  //p->mhd_data.divA *= h_inv_dim_plus_one * rho_inv;
  for (int i = 0; i < 3; i++)
    p->mhd_data.BPred[i] *= h_inv_dim_plus_one * rho_inv;
}

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
  p->mhd_data.divA = 0.f;
  p->mhd_data.divB = 0.f;

  p->mhd_data.BSmooth[0] = 0.f;
  p->mhd_data.BSmooth[1] = 0.f;
  p->mhd_data.BSmooth[2] = 0.f;
  p->mhd_data.Q0 = 0.f;  // XXX make union for clarification

  /* Curl B*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.curlB[k] = 0.f;
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

  // Self Contribution
  for (int i = 0; i < 3; i++)
    p->mhd_data.BSmooth[i] += p->mass * kernel_root * p->mhd_data.BPred[i];
  p->mhd_data.Q0 += p->mass * kernel_root;

  for (int i = 0; i < 3; i++)
    p->mhd_data.BPred[i] = p->mhd_data.BSmooth[i] / p->mhd_data.Q0;

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
  p->mhd_data.divA = 0.0f;
  p->mhd_data.curlB[0] = 0.0f;
  p->mhd_data.curlB[1] = 0.0f;
  p->mhd_data.curlB[2] = 0.0f;
}

/**
 * @brief Returns the Gauge Scalar Phi evolution
 * time the particle. Gauge all variables in full step
 *
 * @param p The particle of interest
 * @param Gauge Gauge
 */
__attribute__((always_inline)) INLINE static float hydro_get_dGau_dt(
    const struct part *restrict p, const float Gauge, const struct cosmology *c)
    //,  const float mu_0) 
    {

  //const float v_sig = mhd_get_magnetosonic_speed(p, c->a, mu_0);
  //const float v_sig = p->viscosity.v_sig;
  const float v_sig = hydro_get_signal_velocity(p);
  const float afac1 = pow(c->a_factor_sound_speed, 2.f) * c->a * c->a;
  const float afac2 = c->a_factor_sound_speed * c->a;
  // Everything is with a2 for the time integration
  /* Hyperbolic term */
  const float Source_Term = 1.f * afac1 * p->mhd_data.divA * (v_sig * v_sig);
  /* Parabolic evolution term */
  const float Damping_Term = 4.f * afac2 * v_sig * Gauge / p->h;
  /* Density change term */
  const float DivV_Term =  hydro_get_div_v(p) * Gauge * c->a * c->a ;
  /* Cosmological term */
  const float Hubble_Term = (2.f + mhd_comoving_factor + 3.f) * c->H * c->a * c->a * Gauge;

  return (- Source_Term - Damping_Term - DivV_Term + Hubble_Term);
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

//  const float mu_0_1 = 1.f / mu_0;
//  const float pressure = hydro_get_comoving_pressure(p);
//  const float b2 = (p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
//                    p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
//                    p->mhd_data.BPred[2] * p->mhd_data.BPred[2]);
  /* Estimation of the tensile instability due divB */
//  p->mhd_data.Q0 = pressure / (b2 / 2.0f * mu_0_1);  // Plasma Beta
//  p->mhd_data.Q0 =
//      p->mhd_data.Q0 < 10.0f ? 1.0f : 0.0f;  // No correction if not magnetized
  /* divB contribution */
//  const float ACC_corr = fabs(p->mhd_data.divB * sqrt(b2) * mu_0_1);
  // this should go with a /p->h, but I
  // take simplify becasue of ACC_mhd also.
  /* isotropic magnetic presure */
  // add the correct hydro acceleration?
//  const float ACC_mhd = (b2 / p->h) * mu_0_1;
  /* Re normalize the correction in the momentum from the DivB errors*/
//  p->mhd_data.Q0 =
//      ACC_corr > ACC_mhd ? p->mhd_data.Q0 * ACC_mhd / ACC_corr : p->mhd_data.Q0;
  //float delim = 1.f;
  //float change_Gau =
  //    hydro_get_dGau_dt(p, xp->mhd_data.Gau_full, cosmo) * dt_alpha;
  //delim = fabs(change_Gau/xp->mhd_data.Gau_full) < 0.25f ?  1.f : 0.25f * fabs(xp->mhd_data.Gau_full / change_Gau);
  
  //xp->mhd_data.Gau_full += change_Gau * delim;
  
  //p->mhd_data.Gau = xp->mhd_data.Gau_full;
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
  /* Zeroes Induction equation */
  p->mhd_data.dAdt[0] = 0.0f;
  p->mhd_data.dAdt[1] = 0.0f;
  p->mhd_data.dAdt[2] = 0.0f;

  /* Save forces*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.tot_mag_F[k] = 0.0f;
  }
  /* Save induction sources*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.Adv_A_source[k] = 0.0f;
    p->mhd_data.Diff_A_source[k] = 0.0f;
    p->mhd_data.Delta_A[k] = 0.0f;
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
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {
  /* Re-set the predicted magnetic flux densities */
  p->mhd_data.APred[0] = xp->mhd_data.Afull[0];
  p->mhd_data.APred[1] = xp->mhd_data.Afull[1];
  p->mhd_data.APred[2] = xp->mhd_data.Afull[2];
  p->mhd_data.Gau = xp->mhd_data.Gau_full;
  ///xp->mhd_data.Afull[0] = p->mhd_data.APred[0] ;
  ///xp->mhd_data.Afull[1] = p->mhd_data.APred[1] ;
  ///xp->mhd_data.Afull[2] = p->mhd_data.APred[2] ;
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
    struct part *p, struct xpart *xp, const float dt_drift,
    const float dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props, const float mu_0) {

  /* Predict the VP magnetic field */  ///// May we need to predict B? XXX
  p->mhd_data.APred[0] += p->mhd_data.dAdt[0] * dt_therm;
  p->mhd_data.APred[1] += p->mhd_data.dAdt[1] * dt_therm;
  p->mhd_data.APred[2] += p->mhd_data.dAdt[2] * dt_therm;
  ///xp->mhd_data.Afull[0] += p->mhd_data.dAdt[0] * dt_therm;
  ///xp->mhd_data.Afull[1] += p->mhd_data.dAdt[1] * dt_therm;
  ///xp->mhd_data.Afull[2] += p->mhd_data.dAdt[2] * dt_therm;
  float delim = 1.f;
  float change_Gau =
      hydro_get_dGau_dt(p, xp->mhd_data.Gau_full, cosmo) * dt_therm;
  //float delim = fabs(change_Gau/p->mhd_data.Gau) < 0.25f ?  1.f : 0.25f * fabs(p->mhd_data.Gau / change_Gau);
  p->mhd_data.Gau += change_Gau * delim;
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

  float a_fac = (mhd_comoving_factor + 2.f) * cosmo->a * cosmo->a * cosmo->H;
  p->mhd_data.dAdt[0] -= a_fac * p->mhd_data.APred[0];
  p->mhd_data.dAdt[1] -= a_fac * p->mhd_data.APred[1];
  p->mhd_data.dAdt[2] -= a_fac * p->mhd_data.APred[2];
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
  // Slope delimiter
  /*const float dAdt_2 = 
	  p->mhd_data.dAdt[0] * p->mhd_data.dAdt[0] +
	  p->mhd_data.dAdt[1] * p->mhd_data.dAdt[1] +
	  p->mhd_data.dAdt[2] * p->mhd_data.dAdt[2] ;
  const float Afull_2 = 
	  xp->mhd_data.Afull[0] * xp->mhd_data.Afull[0] +
	  xp->mhd_data.Afull[1] * xp->mhd_data.Afull[1] +
	  xp->mhd_data.Afull[2] * xp->mhd_data.Afull[2];
  float delim = dAdt_2 * dt_therm * dt_therm / Afull_2 < 1.0f  ? 1.f : sqrt(Afull_2 / dAdt_2) / dt_therm ;
  */
  float delim = 1.f;
  xp->mhd_data.Afull[0] += p->mhd_data.dAdt[0] * dt_therm * delim;
  xp->mhd_data.Afull[1] += p->mhd_data.dAdt[1] * dt_therm * delim;
  xp->mhd_data.Afull[2] += p->mhd_data.dAdt[2] * dt_therm * delim;
  ///p->mhd_data.APred[0] += p->mhd_data.dAdt[0] * dt_therm * delim;
  ///p->mhd_data.APred[1] += p->mhd_data.dAdt[1] * dt_therm * delim;
  ///p->mhd_data.APred[2] += p->mhd_data.dAdt[2] * dt_therm * delim;

  float change_Gau =
      hydro_get_dGau_dt(p, xp->mhd_data.Gau_full, cosmo) * dt_therm;
  //delim = fabs(change_Gau/xp->mhd_data.Gau_full) < 0.25f ?  1.f : 0.25f * fabs(xp->mhd_data.Gau_full / change_Gau);
  xp->mhd_data.Gau_full += change_Gau * delim;
   
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
	
  const float a_fact = pow(cosmo->a, - mhd_comoving_factor - 1.f);
  /* Set Restitivity Eta */
  p->mhd_data.resistive_eta = hydro_props->mhd.mhd_eta;

  p->mhd_data.APred[0] *= a_fact;
  p->mhd_data.APred[1] *= a_fact;
  p->mhd_data.APred[2] *= a_fact;

  xp->mhd_data.Afull[0] = p->mhd_data.APred[0];
  xp->mhd_data.Afull[1] = p->mhd_data.APred[1];
  xp->mhd_data.Afull[2] = p->mhd_data.APred[2];
  
  xp->mhd_data.Gau_full = 0.f;
  p->mhd_data.Gau = 0.f;

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

#endif /* SWIFT_VECTOR_POTENTIAL_MHD_H */
