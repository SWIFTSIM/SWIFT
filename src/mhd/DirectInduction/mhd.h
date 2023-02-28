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
#ifndef SWIFT_DIRECT_INDUCTION_MHD_H
#define SWIFT_DIRECT_INDUCTION_MHD_H
#include "minmax.h"

#include <float.h>

__attribute__((always_inline)) INLINE static float mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const float mu_0) {

  const float rho = p->rho;
  const float b2 = p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
                   p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
                   p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return 0.5f * p->mass * b2 * rho * rho / mu_0;
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
  return (p->v[0] * p->mhd_data.B_over_rho[0] + p->v[1] * p->mhd_data.B_over_rho[1] +
          p->v[2] * p->mhd_data.B_over_rho[2]) * rho;
}

__attribute__((always_inline)) INLINE static float mhd_get_divB_error(
    const struct part *p, const struct xpart *xp) {

  const float rho = p->rho;
  const float B2  = p->mhd_data.B_over_rho[0] * p->mhd_data.B_over_rho[0] +
                    p->mhd_data.B_over_rho[1] * p->mhd_data.B_over_rho[1] +
                    p->mhd_data.B_over_rho[2] * p->mhd_data.B_over_rho[2];
  return fabs(p->mhd_data.divB * p->h / sqrtf(B2 * rho * rho + 1.e-18));
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
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo) {

  return FLT_MAX;
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
    const float dx[3], const struct part *pi, const struct part *pj,
    const float mu_ij, const float beta, const float a, const float mu_0) {

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
  
  const float termi = max(0.0f, c2effi * c2effi - 4.0f * c2i * (Bri * r_inv) * (Bri * r_inv) / (mu_0 * rhoi));
  const float termj = max(0.0f, c2effj * c2effj - 4.0f * c2j * (Brj * r_inv) * (Brj * r_inv) / (mu_0 * rhoj));
    
  const float v_sig2i =
      0.5f * (c2effi + sqrtf(termi));
  const float v_sig2j =
      0.5f * (c2effj + sqrtf(termj));

  const float v_sigi = sqrtf(v_sig2i);
  //pi->mhd_data.v_fm = v_sigi;
  const float v_sigj = sqrtf(v_sig2j);
  //pj->mhd_data.v_fm = v_sigj;

  const float v_sig =
      v_sigi + v_sigj - const_viscosity_beta * mu_ij;

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
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props) {}

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
    struct part *p) {

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
		struct part *p, const struct cosmology *cosmo, const float mu_0) {

  // const float mu_0 = 1.0f;
    
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;    
    
  /* Recover some data */
  const float rho = p->rho;
  float B[3];
  B[0] = p->mhd_data.B_over_rho[0] * rho;
  B[1] = p->mhd_data.B_over_rho[1] * rho;
  B[2] = p->mhd_data.B_over_rho[2] * rho;
  
  /* B squared */
  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
  
  /* Compute sound speeds and signal velocity */
  const float cs = p->force.soundspeed;
  const float cs2 = cs * cs;
  const float v_A2 = B2 / (rho * mu_0);
  const float ch = sqrtf(cs2 + v_A2);
  
  /* Dedner cleaning scalar time derivative */
  // const float v_sig = hydro_get_signal_velocity(p);
  // const float v_sig2 = v_sig * v_sig;
  const float div_B = p->mhd_data.B_mon;
  const float div_v = hydro_get_div_v(p);
  const float psi_over_ch = p->mhd_data.psi_over_ch;
  p->mhd_data.psi_over_ch_dt =
      - ch * div_B - dedner_gamma * psi_over_ch * div_v - psi_over_ch * ch * h_inv;
}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
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

  /*
  if (fabs(delta_Bx) < 0.5f * fabs(xp->mhd_data.B_over_rho_full[0])) {
    xp->mhd_data.B_over_rho_full[0] = xp->mhd_data.B_over_rho_full[0] + delta_Bx;
  }
  else {
    xp->mhd_data.B_over_rho_full[0] = 0.5f * xp->mhd_data.B_over_rho_full[0];
  }

  if (fabs(delta_By) < 0.5f * fabs(xp->mhd_data.B_over_rho_full[1])) {
    xp->mhd_data.B_over_rho_full[1] = xp->mhd_data.B_over_rho_full[1] + delta_By;
  }
  else {
    xp->mhd_data.B_over_rho_full[1] = 0.5f * xp->mhd_data.B_over_rho_full[1];
  }

  if (fabs(delta_Bz) < 0.5f * fabs(xp->mhd_data.B_over_rho_full[2])) {
    xp->mhd_data.B_over_rho_full[2] = xp->mhd_data.B_over_rho_full[2] + delta_Bz;
  }
  else {
    xp->mhd_data.B_over_rho_full[2] = 0.5f * xp->mhd_data.B_over_rho_full[2];
  }
  */

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
    struct part *p, struct xpart *xp, const struct mhd_global_data *mhd_data,
    const double Lsize) {

  mhd_reset_acceleration(p);
  mhd_init_part(p);
}

/**
 * @brief Print out the mhd fields of a particle.
 *
 * Function used for debugging purposes.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_debug_particle(
    const struct part *p, const struct xpart *xp) {

  warning("B/rho= [%.3e,%.3e,%.3e] d(B/rho)/dt= [%.3e,%.3e,%.3e]\n"
	  "divB= %.3e psi/ch= %.3e psi/ch_dt= %.3e",
	  p->mhd_data.B_over_rho[0],
	  p->mhd_data.B_over_rho[1],
	  p->mhd_data.B_over_rho[2],
	  p->mhd_data.B_over_rho_dt[0],
	  p->mhd_data.B_over_rho_dt[1],
	  p->mhd_data.B_over_rho_dt[2],
	  p->mhd_data.divB,
	  p->mhd_data.psi_over_ch,
	  p->mhd_data.psi_over_ch_dt);
  
  /*
  warning("[PID%lld] part:", p->id);
  warning(
      "[PID%lld] "
      "x=[%.6g, %.6g, %.6g], v=[%.3g, %.3g, %.3g], "
      "a=[%.3g, %.3g, %.3g], "
      "m=%.3g, u=%.3g, du/dt=%.3g, P=%.3g, c_s=%.3g, "
      "h=%.3g, dh/dt=%.3g, wcount=%.3g, rho=%.3g, "
      "dh_drho=%.3g, time_bin=%d wakeup=%d",
      p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
      p->a_hydro[0], p->a_hydro[1], p->a_hydro[2], p->mass, p->u, p->u_dt,
      hydro_get_comoving_pressure(p), p->force.soundspeed, p->h,
      p->force.h_dt, p->density.wcount, p->rho, p->density.rho_dh, p->time_bin,
      p->limiter_data.wakeup);
  warning("B/rho= [%.3e,%.3e,%.3e] d(B/rho)/dt= [%.3e,%.3e,%.3e]\n"                                                  
          "divB= %.3e psi= %.3e psi_dt= %.3e",                                                                                 p->mhd_data.B_over_rho[0],                                                                                           p->mhd_data.B_over_rho[1],                                                                                           p->mhd_data.B_over_rho[2],                                                                                           p->mhd_data.B_over_rho_dt[0],                                                                                        p->mhd_data.B_over_rho_dt[1],                                                                                        p->mhd_data.B_over_rho_dt[2],                                                                                        p->mhd_data.divB,                                                                                                    p->mhd_data.psi,                                                                                                     p->mhd_data.psi_dt); 
  */

  if (xp != NULL) {
    warning("[PID%lld] xpart:", p->id);
    warning("[PID%lld] v_full=[%.3g, %.3g, %.3g]", p->id, xp->v_full[0],
            xp->v_full[1], xp->v_full[2]);
  }

}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_H */
