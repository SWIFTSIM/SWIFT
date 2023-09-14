/*******************************************************************************
 * Thisfile is part o SWIFT.
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

extern float diffusion_eta;

#include <float.h>

__attribute__((always_inline)) INLINE static float mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const float mu_0) {

  const float rho = p->rho;
  const float B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return 0.5f * p->mass * B_over_rho2 * rho / mu_0;
}

__attribute__((always_inline)) INLINE static float mhd_get_Bms(
    const struct part *p, const struct xpart *xp) {

  const float rho = p->rho;
  const float B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return B_over_rho2 * rho * rho;
}

__attribute__((always_inline)) INLINE static float mhd_get_magnetic_divergence(
    const struct part *p, const struct xpart *xp) {

  return p->mhd_data.B_mon;
}

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

__attribute__((always_inline)) INLINE static float mhd_get_divB_error(
    const struct part *p, const struct xpart *xp) {

  const float rho = p->rho;
  const float B_over_rho2 =
      p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
      p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
      p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return fabs(p->mhd_data.B_mon) * p->h /
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

  const float divB = p->mhd_data.B_mon;

  float curl_B[3];
  curl_B[0] = p->mhd_data.curl_B[0];
  curl_B[1] = p->mhd_data.curl_B[1];
  curl_B[2] = p->mhd_data.curl_B[2];

  const float curl_B_norm = sqrtf(
      curl_B[0] * curl_B[0] + curl_B[1] * curl_B[1] + curl_B[2] * curl_B[2]);
 
  const float dt_B_factor = fmax(fabs(divB), fabs(curl_B_norm));
  
  const float dt_B_derivatives = dt_B_factor != 0.f
             ? hydro_properties->CFL_condition * 
                   sqrt(p->rho * mu_0/ (dt_B_factor * dt_B_factor))
             : FLT_MAX;

  const float dt_eta = diffusion_eta != 0.f 
             ? hydro_properties->CFL_condition * p->h * p->h / diffusion_eta
             : FLT_MAX;

  return fmin(dt_B_derivatives, dt_eta);
}

/**
 * @brief Compute the signal velocity between two gas particles
 *
 * This is eq. (131) of Price D., JCoPh, 2012, Vol. 231, Issue 3
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

  /* Get r and 1/r. */
  const float r2 = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  float Bi[3];
  float Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  /* B squared */
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float c2i = ci * ci;
  const float c2j = cj * cj;
  const float v_A2i = B2i / (rhoi * mu_0);
  const float v_A2j = B2j / (rhoj * mu_0);
  const float c2effi = c2i + v_A2i;
  const float c2effj = c2j + v_A2j;

  const float termi =
      max(0.0f, c2effi * c2effi -
                    4.0f * c2i * (Bri * r_inv) * (Bri * r_inv) / (mu_0 * rhoi));
  const float termj =
      max(0.0f, c2effj * c2effj -
                    4.0f * c2j * (Brj * r_inv) * (Brj * r_inv) / (mu_0 * rhoj));

  const float v_sig2i = 0.5f * (c2effi + sqrtf(termi));
  const float v_sig2j = 0.5f * (c2effj + sqrtf(termj));

  const float v_sigi = sqrtf(v_sig2i);
  const float v_sigj = sqrtf(v_sig2j);

  const float v_sig = v_sigi + v_sigj - const_viscosity_beta * mu_ij;

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
    const struct cosmology *cosmo, const struct hydro_props *hydro_props) {}

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
  p->mhd_data.B_mon = 0.0f;
  p->mhd_data.curl_B[0] = 0.0f;
  p->mhd_data.curl_B[1] = 0.0f;
  p->mhd_data.curl_B[2] = 0.0f;
}

/**
 * @brief Finishes the gradient calculation.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void mhd_end_gradient(
    struct part *p) {}

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
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {}

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
    const struct hydro_props *hydro_props, const float dt_alpha) {}

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
 */
__attribute__((always_inline)) INLINE static void mhd_predict_extra(
    struct part *p, const struct xpart *xp, const float dt_drift,
    const float dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Predict the magnetic flux density */
  p->mhd_data.B_over_rho[0] += p->mhd_data.B_over_rho_dt[0] * dt_therm;
  p->mhd_data.B_over_rho[1] += p->mhd_data.B_over_rho_dt[1] * dt_therm;
  p->mhd_data.B_over_rho[2] += p->mhd_data.B_over_rho_dt[2] * dt_therm;
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

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  /*
  const float rho = p->rho;
  const float rho_inv = 1.0f / rho;
    
  float curlB[3];
  curlB[0] = p->mhd_data.curl_B[0];
  curlB[1] = p->mhd_data.curl_B[1];
  curlB[2] = p->mhd_data.curl_B[2];

  const float normCurlB2 =
      curlB[0] * curlB[0] + curlB[1] * curlB[1] + curlB[2] * curlB[2];
  */

  const float ch = p->viscosity.v_sig;

  /* Physical resistivity contribution to energy equation */
  // p->u_dt += hydro_props->mhd.mhd_eta * normCurlB2 * rho_inv / mu_0;

  /* Dedner cleaning scalar time derivative */
  const float hyp = hydro_props->mhd.hyp_dedner;
  const float hyp_divv = hydro_props->mhd.hyp_dedner_divv;
  const float par = hydro_props->mhd.par_dedner;

  const float div_B = p->mhd_data.B_mon;
  const float div_v = hydro_get_div_v(p);
  const float psi_over_ch = p->mhd_data.psi_over_ch;

  p->mhd_data.psi_over_ch_dt = -hyp * ch * div_B -
                               hyp_divv * psi_over_ch * div_v -
                               par * psi_over_ch * ch * h_inv;
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

  /* Integrate the Dedner scalar forward in time */
  const float delta_psi_over_ch = p->mhd_data.psi_over_ch_dt * dt_therm;

  /* Do not decrease the magnetic flux density by more than a factor of 2*/
  xp->mhd_data.B_over_rho_full[0] = xp->mhd_data.B_over_rho_full[0] + delta_Bx;
  xp->mhd_data.B_over_rho_full[1] = xp->mhd_data.B_over_rho_full[1] + delta_By;
  xp->mhd_data.B_over_rho_full[2] = xp->mhd_data.B_over_rho_full[2] + delta_Bz;

  /* Integrate Dedner scalar in time */
  p->mhd_data.psi_over_ch = p->mhd_data.psi_over_ch + delta_psi_over_ch;
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

  /* Convert B into B/rho */
  p->mhd_data.B_over_rho[0] /= p->rho;
  p->mhd_data.B_over_rho[1] /= p->rho;
  p->mhd_data.B_over_rho[2] /= p->rho;

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
