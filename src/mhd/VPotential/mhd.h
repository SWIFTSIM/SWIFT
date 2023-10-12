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
    const struct part *p, const struct xpart *xp, const float mu_0) {

  const float b2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  return 0.5f * b2 / mu_0 * p->mass / p->rho;
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
  return fabs(p->mhd_data.divB * p->h / sqrt(b2 + 1.e-18));
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

  float dt_divB =
      p->mhd_data.divB != 0.0f
          ? cosmo->a * hydro_properties->CFL_condition *
                sqrtf(p->rho / (p->mhd_data.divB * p->mhd_data.divB) * mu_0)
          : FLT_MAX;
  const float Reta = p->mhd_data.Reta;
  const float dt_eta = Reta != 0.0f
                           ? cosmo->a * hydro_properties->CFL_condition * p->h *
                                 p->h / Reta * 0.5
                           : FLT_MAX;

  return fminf(dt_eta, dt_divB);
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

  const float a_fac =
      2.f * mhd_comoving_factor + 3.f + 3.f * (hydro_gamma - 1.f);
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float r2 = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r2 ? 1.f / sqrt(r2) : 0.0f;

  const float b2_i = (pi->mhd_data.BPred[0] * pi->mhd_data.BPred[0] +
                      pi->mhd_data.BPred[1] * pi->mhd_data.BPred[1] +
                      pi->mhd_data.BPred[2] * pi->mhd_data.BPred[2]);
  const float b2_j = (pj->mhd_data.BPred[0] * pj->mhd_data.BPred[0] +
                      pj->mhd_data.BPred[1] * pj->mhd_data.BPred[1] +
                      pj->mhd_data.BPred[2] * pj->mhd_data.BPred[2]);
  const float vcsa2_i = ci * ci + pow(a, a_fac) * b2_i / pi->rho * 0.5 / mu_0;
  const float vcsa2_j = cj * cj + pow(a, a_fac) * b2_j / pj->rho * 0.5 / mu_0;
  float Bpro2_i =
      (pi->mhd_data.BPred[0] * dx[0] + pi->mhd_data.BPred[1] * dx[1] +
       pi->mhd_data.BPred[2] * dx[2]) *
      r_inv;
  Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrtf(
      0.5 * (vcsa2_i +
             sqrtf(max((vcsa2_i * vcsa2_i - 4.f * ci * ci * pow(a, a_fac) *
                                                Bpro2_i / pi->rho * 0.5 / mu_0),
                       0.f))));
  float Bpro2_j =
      (pj->mhd_data.BPred[0] * dx[0] + pj->mhd_data.BPred[1] * dx[1] +
       pj->mhd_data.BPred[2] * dx[2]) *
      r_inv;
  Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrtf(
      0.5 * (vcsa2_j +
             sqrtf(max((vcsa2_j * vcsa2_j - 4.f * cj * cj * pow(a, a_fac) *
                                                Bpro2_j / pj->rho * 0.5 / mu_0),
                       0.f))));

  return (mag_speed_i + mag_speed_j - beta / 2. * mu_ij);
}

/**
 * @brief Returns the Gauge Scalar Phi evolution
 * time the particle. Gauge all variables in full step
 *
 * @param p The particle of interest
 * @param Gauge Gauge
 */
__attribute__((always_inline)) INLINE static float hydro_get_dGau_dt(
    const struct part *restrict p, const float Gauge,
    const struct cosmology *c) {

  const float v_sig = hydro_get_signal_velocity(p);
  const float afac1 = pow(c->a, 2.f * mhd_comoving_factor);
  const float afac2 = pow(c->a, 1.f + mhd_comoving_factor);

  return (-p->mhd_data.divA * v_sig * v_sig * 0.1 * afac1 -
          2.0f * v_sig * Gauge / p->h * afac2 -
          (2.f + mhd_comoving_factor) * c->a * c->a * c->H * Gauge) /
         2.f;
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
    struct part *p) {

  p->mhd_data.divA = 0.f;
  for (int i = 0; i < 3; i++) p->mhd_data.BPred[i] = 0.f;
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

  const float h_inv_dim_plus_one = pow_dimension(1.f / p->h) / p->h;
  const float rho_inv = 1.f / p->rho;
  p->mhd_data.divA *= h_inv_dim_plus_one * rho_inv;
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

  p->mhd_data.divB = 0.f;

  p->mhd_data.BSmooth[0] = 0.f;
  p->mhd_data.BSmooth[1] = 0.f;
  p->mhd_data.BSmooth[2] = 0.f;
  //  p->mhd_data.GauSmooth = 0.f;
  p->mhd_data.Q0 = 0.f;  // XXX make union for clarification
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
    const struct hydro_props *hydro_props, const float dt_alpha,
    const float mu_0) {

  const float mu_0_1 = 1.f / mu_0;
  const float pressure = hydro_get_comoving_pressure(p);
  const float b2 = (p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                    p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                    p->mhd_data.BPred[2] * p->mhd_data.BPred[2]);
  /* Estimation of the tensile instability due divB */
  p->mhd_data.Q0 = pressure / (b2 / 2.0f * mu_0_1);  // Plasma Beta
  p->mhd_data.Q0 =
      p->mhd_data.Q0 < 10.0f ? 1.0f : 0.0f;  // No correction if not magnetized
  /* divB contribution */
  const float ACC_corr = fabs(p->mhd_data.divB * sqrt(b2) * mu_0_1);
  // this should go with a /p->h, but I
  // take simplify becasue of ACC_mhd also.
  /* isotropic magnetic presure */
  // add the correct hydro acceleration?
  const float ACC_mhd = (b2 / p->h) * mu_0_1;
  /* Re normalize the correction in the momentum from the DivB errors*/
  p->mhd_data.Q0 =
      ACC_corr > ACC_mhd ? p->mhd_data.Q0 * ACC_mhd / ACC_corr : p->mhd_data.Q0;
  //     p->mhd_data.Q0 = 0.0f;
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

  p->mhd_data.Gau = xp->mhd_data.Gau_full;
  // p->mhd_data.Gau = p->mhd_data.GauSmooth;

  p->mhd_data.APred[0] = xp->mhd_data.APot_full[0];
  p->mhd_data.APred[1] = xp->mhd_data.APot_full[1];
  p->mhd_data.APred[2] = xp->mhd_data.APot_full[2];
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

  /* Predict the VP magnetic field */  ///// May we need to predict B? XXX
  p->mhd_data.APred[0] += p->mhd_data.dAdt[0] * dt_therm;
  p->mhd_data.APred[1] += p->mhd_data.dAdt[1] * dt_therm;
  p->mhd_data.APred[2] += p->mhd_data.dAdt[2] * dt_therm;

  float change_Gau = hydro_get_dGau_dt(p, p->mhd_data.Gau, cosmo) * dt_therm;
  p->mhd_data.Gau += change_Gau;
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
  // p->mhd_data.dAdt[0] = 0.0f;
  // p->mhd_data.dAdt[1] = 0.0f;
  // p->mhd_data.dAdt[2] = 0.0f;
  float a_fac = (2.f + mhd_comoving_factor) * cosmo->a * cosmo->a * cosmo->H;
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

  /* Integrate the magnetic field */
  xp->mhd_data.APot_full[0] += p->mhd_data.dAdt[0] * dt_therm;
  xp->mhd_data.APot_full[1] += p->mhd_data.dAdt[1] * dt_therm;
  xp->mhd_data.APot_full[2] += p->mhd_data.dAdt[2] * dt_therm;
  float change_Gau = hydro_get_dGau_dt(p, p->mhd_data.Gau, cosmo) * dt_therm;
  xp->mhd_data.Gau_full += change_Gau;
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
  p->mhd_data.Reta = hydro_props->mhd.mhd_eta;
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

  // const float Lsize = s->dims[0];
  const float define_Bfield_in_ics = mhd_data->define_Bfield_in_ics;
  const float define_Afield_in_ics = mhd_data->define_Afield_in_ics;
  const float Nvort = 10;
  const float Bini = define_Bfield_in_ics;
  const float Aini = define_Afield_in_ics / (2 * M_PI * Nvort) * Lsize;
  if (define_Bfield_in_ics) {
    p->mhd_data.BPred[0] = Bini * (sin(2 * M_PI * p->x[2] / Lsize * Nvort) +
                                   cos(2 * M_PI * p->x[1] / Lsize * Nvort));
    p->mhd_data.BPred[1] = Bini * (sin(2 * M_PI * p->x[0] / Lsize * Nvort) +
                                   cos(2 * M_PI * p->x[2] / Lsize * Nvort));
    p->mhd_data.BPred[2] = Bini * (sin(2 * M_PI * p->x[1] / Lsize * Nvort) +
                                   cos(2 * M_PI * p->x[0] / Lsize * Nvort));
    p->mhd_data.APred[0] = Aini * (sin(2 * M_PI * p->x[2] / Lsize * Nvort) +
                                   cos(2 * M_PI * p->x[1] / Lsize * Nvort));
    p->mhd_data.APred[1] = Aini * (sin(2 * M_PI * p->x[0] / Lsize * Nvort) +
                                   cos(2 * M_PI * p->x[2] / Lsize * Nvort));
    p->mhd_data.APred[2] = Aini * (sin(2 * M_PI * p->x[1] / Lsize * Nvort) +
                                   cos(2 * M_PI * p->x[0] / Lsize * Nvort));
  }

  // p->mhd_data.BPred[0] = 0.0;
  // p->mhd_data.BPred[1] = 0.0;
  // p->mhd_data.BPred[2] = 0.0;

  // p->mhd_data.APred[0] = 0.0;
  // p->mhd_data.APred[1] = 0.0;
  // p->mhd_data.APred[2] = 0.0;
  xp->mhd_data.APot_full[0] = p->mhd_data.APred[0];
  xp->mhd_data.APot_full[1] = p->mhd_data.APred[1];
  xp->mhd_data.APot_full[2] = p->mhd_data.APred[2];
  xp->mhd_data.Gau_full = 0.0f * p->mhd_data.Gau;
  p->mhd_data.divB = 0.0f;

  mhd_reset_acceleration(p);
  mhd_init_part(p);
}

#endif /* SWIFT_VECTOR_POTENTIAL_MHD_H */
