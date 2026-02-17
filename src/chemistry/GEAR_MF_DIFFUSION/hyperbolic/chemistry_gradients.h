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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GRADIENTS_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GRADIENTS_H

#include "../chemistry_getters.h"
#include "../chemistry_setters.h"
#include "../chemistry_slope_limiters_cell.h"
#include "../chemistry_slope_limiters_face.h"
#include "../chemistry_unphysical.h"

/**
 * @file src/chemistry/GEAR_MFM_diffusion/hyperbolic/chemistry_gradients.h
 * @brief Header file for hyperbolic diffusion gradient functions.
 */

/**
 * @brief Gradients reconstruction. Predict the value at point x_ij given
 * current values at particle positions and gradients at particle positions.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param half_dt Time-step to integrate.
 *
 * @param cosmo The #cosmology.
 * @param chem_data The global properties of the chemistry scheme.
 * @param dUi (return) Resulting time-predicted diffusion state of particle i
(in physical units).
 * @param Uj (return) Resulting time-predicted and limited diffusion state of
 * particle j (in physical units).
 */
__attribute__((always_inline)) INLINE static void
chemistry_gradients_time_extrapolate(
    const struct part *restrict pi, const struct part *restrict pj,
    float half_dt, double grad_qi[3], double grad_qj[3], double dFx_i[3],
    double dFy_i[3], double dFz_i[3], double dFx_j[3], double dFy_j[3],
    double dFz_j[3], const struct cosmology *cosmo,
    const struct chemistry_global_data *chem_data, double dUi[4],
    double dUj[4]) {

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get drhoZ/dt = - div Flux */
  const double drhoZ_dt_i = -(dFx_i[0] + dFy_i[1] + dFz_i[2]);
  const double drhoZ_dt_j = -(dFx_j[0] + dFy_j[1] + dFz_j[2]);

  /* Compute dF/dt = - F/tau - K/tau * grad q */
  const double tau_inv_i = 1.0 / chi->tau;
  const double tau_inv_j = 1.0 / chj->tau;

  double Ki[3][3];
  chemistry_get_physical_matrix_K(pi, chem_data, cosmo, Ki);
  double Kj[3][3];
  chemistry_get_physical_matrix_K(pj, chem_data, cosmo, Kj);

  double dF_dt_homogeneous_i[3] = {0.0};
  double dF_dt_homogeneous_j[3] = {0.0};

  /* If tau == 0, then tau*dF/dt = 0 and the diffusion is parabolic. Therefore,
     we cannot extrapolate the flux equation in time. */
  if (chi->tau != 0) {
    if (chem_data->diffusion_mode == anisotropic_gradient) {
      for (int i = 0; i < 3; ++i) {
	for (int j = 0; j < 3; ++j) {
	  dF_dt_homogeneous_i[i] -= tau_inv_i*Ki[i][j] * grad_qi[j];
	}
      }
    } else {
      dF_dt_homogeneous_i[0] = -chi->kappa * grad_qi[0];
      dF_dt_homogeneous_i[1] = -chi->kappa * grad_qi[1];
      dF_dt_homogeneous_i[2] = -chi->kappa * grad_qi[2];
    }
  }

  /* TODO: Time integrate with the source term. Use the ODE solution */

  if (chj->tau != 0) {
    if (chem_data->diffusion_mode == anisotropic_gradient) {
      for (int i = 0; i < 3; ++i) {
	for (int j = 0; j < 3; ++j) {
	  dF_dt_homogeneous_j[i] -= tau_inv_j*Kj[i][j] * grad_qj[j];
	}
      }
    } else {
      dF_dt_homogeneous_j[0] = -chj->kappa * grad_qj[0];
      dF_dt_homogeneous_j[1] = -chj->kappa * grad_qj[1];
      dF_dt_homogeneous_j[2] = -chj->kappa * grad_qj[2];
    }
  }
  /* TODO: Time integrate with the source term. Use the ODE solution */

  dUi[0] = drhoZ_dt_i * half_dt;
  dUi[1] = dF_dt_homogeneous_i[0] * half_dt;
  dUi[2] = dF_dt_homogeneous_i[1] * half_dt;
  dUi[3] = dF_dt_homogeneous_i[2] * half_dt;

  dUj[0] = drhoZ_dt_j * half_dt;
  dUj[1] = dF_dt_homogeneous_j[0] * half_dt;
  dUj[2] = dF_dt_homogeneous_j[1] * half_dt;
  dUj[3] = dF_dt_homogeneous_j[2] * half_dt;
}

/**
 * @brief Check and correct unphysical states due to extrapolation.
 *
 * This function also handles particle with negative metal masses.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param metal Metal specie to update
 * @param cosmo The #cosmology.
 * @param Ui (return) Resulting corrected diffusion state of particle i (in physical units).
 * @param Uj (return) Resulting corrected diffusion state of particle j (in physical units).
 */
__attribute__((always_inline)) INLINE static void
chemistry_gradients_correct_unphysical_states(
    const struct part *restrict pi, const struct part *restrict pj, int metal,
    const struct cosmology *cosmo,
    double Ui[4], double Uj[4]) {

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;
  const double mi = hydro_get_mass(pi);
  const double mj = hydro_get_mass(pj);
  const double m_Zi_not_extrapolated =
      chemistry_get_metal_mass_fraction(pi, metal) * mi;
  const double m_Zj_not_extrapolated =
      chemistry_get_metal_mass_fraction(pj, metal) * mj;
  double m_Zi = Ui[0] * mi / hydro_get_physical_density(pi, cosmo);
  double m_Zj = Uj[0] * mj / hydro_get_physical_density(pj, cosmo);

  unsigned int dumb;
  chemistry_check_unphysical_state(&m_Zi, m_Zi_not_extrapolated, mi,
				   /*callloc=*/1, /*element*/ metal, pi->id,
				   /*neg_counter*/ &dumb);
  chemistry_check_unphysical_state(&m_Zj, m_Zj_not_extrapolated, mj,
				   /*callloc=*/1, /*element*/ metal, pj->id,
				   &dumb);

  /* Check that we have meaningful fluxes */
  double flux_i[3] = {Ui[1], Ui[2], Ui[3]};
  double flux_j[3] = {Uj[1], Uj[2], Uj[3]};
  chemistry_check_unphysical_diffusion_flux(flux_i);
  chemistry_check_unphysical_diffusion_flux(flux_j);

  /* If the new masses have been changed, do not extrapolate, use 0th order
     reconstruction and update the state vectors */
  if (m_Zi == m_Zi_not_extrapolated) {
    Ui[0] = m_Zi_not_extrapolated * hydro_get_physical_density(pi, cosmo) / mi;
    flux_i[0] = chi->diffusion_flux[metal][0];
    flux_i[1] = chi->diffusion_flux[metal][1];
    flux_i[2] = chi->diffusion_flux[metal][2];
  }
  if (m_Zj == m_Zj_not_extrapolated) {
    Uj[0] = m_Zj_not_extrapolated * hydro_get_physical_density(pj, cosmo) / mj;
    flux_j[0] = chj->diffusion_flux[metal][0];
    flux_j[1] = chj->diffusion_flux[metal][1];
    flux_j[2] = chj->diffusion_flux[metal][2];
  }

  if (m_Zi_not_extrapolated < 0.0) {
    Ui[0] = 0.0;
    flux_i[0] = 0.0;
    flux_i[1] = 0.0;
    flux_i[2] = 0.0;
  }

  if (m_Zj_not_extrapolated < 0.0) {
    Uj[0] = 0.0;
    flux_j[0] = 0.0;
    flux_j[1] = 0.0;
    flux_j[2] = 0.0;
  }

  /* Something went wrong if we get this one! */
  if (m_Zi_not_extrapolated > mi) {
    Ui[0] = hydro_get_physical_density(pi, cosmo);
    flux_i[0] = 0.0;
    flux_i[1] = 0.0;
    flux_i[2] = 0.0;
  }
  if (m_Zj_not_extrapolated > mj) {
    Uj[0] = hydro_get_physical_density(pj, cosmo);
    flux_j[0] = 0.0;
    flux_j[1] = 0.0;
    flux_j[2] = 0.0;
  }

  /* Now assign the fluxes */
  Ui[1] = flux_i[0];
  Ui[2] = flux_i[1];
  Ui[3] = flux_i[2];
  Uj[1] = flux_j[0];
  Uj[2] = flux_j[1];
  Uj[3] = flux_j[2];
}

/**
 * @brief Gradients reconstruction. Predict the value at point x_ij given
 * current values at particle positions and gradients at particle positions.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param metal Metal specie to update
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i.
 * @param cosmo The #cosmology.
 * @param chem_data The global properties of the chemistry scheme.
 * @param Ui (return) Resulting predicted and limited diffusion state of
 * particle i (in physical units).
 * @param Uj (return) Resulting predicted and limited diffusion state of
 * particle j (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict(
    const struct part *restrict pi, const struct part *restrict pj, int metal,
    const float dx[3], const float r, const float xij_i[3],
    const struct cosmology *cosmo,
    const struct chemistry_global_data *chem_data,
    double Ui[4], double Uj[4]) {

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get the metal density (comoving) and the hyperbolic flux (physical) */
  Ui[0] = chemistry_get_physical_metal_density(pi, metal, cosmo);
  Ui[1] = pi->chemistry_data.diffusion_flux[metal][0];
  Ui[2] = pi->chemistry_data.diffusion_flux[metal][1];
  Ui[3] = pi->chemistry_data.diffusion_flux[metal][2];

  Uj[0] = chemistry_get_physical_metal_density(pj, metal, cosmo);
  Uj[1] = pj->chemistry_data.diffusion_flux[metal][0];
  Uj[2] = pj->chemistry_data.diffusion_flux[metal][1];
  Uj[3] = pj->chemistry_data.diffusion_flux[metal][2];
  /* No need to check unphysical state here: they haven't been touched since
     the call to chemistry_end_density() */

  double grad_rhoZ_i[3];
  double grad_rhoZ_j[3];
  chemistry_get_metal_density_gradients(pi, metal, grad_rhoZ_i);
  chemistry_get_metal_density_gradients(pj, metal, grad_rhoZ_j);

  /* Cell limit the metal density gradients now. The fluxes were already cell
     limited. */
  const double alpha_i = chemistry_slope_limit_quantity(
      /*gradient=*/grad_rhoZ_i,
      /*maxr=    */ chi->limiter.maxr,
      /*value=   */ chemistry_get_comoving_metal_density(pi, metal),
      /*valmin=  */ chi->limiter.rhoZ[metal][0],
      /*valmax=  */ chi->limiter.rhoZ[metal][1],
      /*condition_number*/ pi->geometry.condition_number,
      /*pos_preserve*/ 1,
      /*shoot_tol*/ GEAR_FVPM_DIFFUSION_CELL_LIMITER_SHOOT_TOLERANGE);
  const double alpha_j = chemistry_slope_limit_quantity(
      /*gradient=*/grad_rhoZ_j,
      /*maxr=    */ chj->limiter.maxr,
      /*value=   */ chemistry_get_comoving_metal_density(pj, metal),
      /*valmin=  */ chj->limiter.rhoZ[metal][0],
      /*valmax=  */ chj->limiter.rhoZ[metal][1],
      /*condition_number*/ pj->geometry.condition_number,
      /*pos_preserve*/ 1,
      /*shoot_tol*/ GEAR_FVPM_DIFFUSION_CELL_LIMITER_SHOOT_TOLERANGE);
  chemistry_slope_limit_quantity_apply(grad_rhoZ_i, alpha_i);
  chemistry_slope_limit_quantity_apply(grad_rhoZ_j, alpha_j);

  /* Since the fluxes are stored in physical units, we need to convert the
     gradients to physical units */
  double dFx_i[3], dFy_i[3], dFz_i[3];
  double dFx_j[3], dFy_j[3], dFz_j[3];
  chemistry_get_hyperbolic_flux_gradients(pi, metal, dFx_i, dFy_i, dFz_i);
  chemistry_get_hyperbolic_flux_gradients(pj, metal, dFx_j, dFy_j, dFz_j);

  /* Compute interface position (relative to pj, since we don't need the
     actual position) eqn. (8)
     Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L (rho*Z)
     Note for the flux extrapolation: grad_p = a^{-1} grad_c and xij_i/j_p = a
     * xij_i/j_p. Hence, the scale factors compensate and dUi/j[1-3] end up in
     physical units.
     Physical unit conversion for grad_rhoZ: multiply by a^-4 (a^-3 for density
     and a^-1 for gradient). But we multiply by a distance with a^1. So the
     result is a factor a^3. */
  const double a3_inv = cosmo->a_inv * cosmo->a_inv * cosmo->a_inv;
  double dUi[4], dUj[4];
  dUi[0] = chemistry_gradients_extrapolate_double(grad_rhoZ_i, xij_i)*a3_inv;
  dUi[1] = chemistry_gradients_extrapolate_double(dFx_i, xij_i);
  dUi[2] = chemistry_gradients_extrapolate_double(dFy_i, xij_i);
  dUi[3] = chemistry_gradients_extrapolate_double(dFz_i, xij_i);

  dUj[0] = chemistry_gradients_extrapolate_double(grad_rhoZ_j, xij_j)*a3_inv;
  dUj[1] = chemistry_gradients_extrapolate_double(dFx_j, xij_j);
  dUj[2] = chemistry_gradients_extrapolate_double(dFy_j, xij_j);
  dUj[3] = chemistry_gradients_extrapolate_double(dFz_j, xij_j);

  chemistry_slope_limit_face(Ui, Uj, dUi, dUj, xij_i, xij_j, r);

  /* Now, let's time-extrapolate! */
  double dUi_time_extrapolated[4] = {0.0};
  double dUj_time_extrapolated[4] = {0.0};

  const float mindt =
      (chj->flux.dt > 0.f) ? fminf(chi->flux.dt, chj->flux.dt) : chi->flux.dt;
  const float half_mindt = 0.5 * mindt;

  /* Get the physical diffusion driver gradients */
  double grad_qi[3];
  double grad_qj[3];
  chemistry_get_physical_diffusion_driver_gradients(pi, metal, cosmo, chem_data,
                                                    grad_qi);
  chemistry_get_physical_diffusion_driver_gradients(pj, metal, cosmo, chem_data,
                                                    grad_qj);

  chemistry_gradients_time_extrapolate(
      pi, pj, half_mindt, grad_qi, grad_qj, dFx_i, dFy_i, dFz_i,
      dFx_j, dFy_j, dFz_j, cosmo, chem_data, dUi_time_extrapolated,
      dUj_time_extrapolated);

  /* Apply the extrapolation */
  Ui[0] += dUi[0] + dUi_time_extrapolated[0];
  Ui[1] += dUi[1] + dUi_time_extrapolated[1];
  Ui[2] += dUi[2] + dUi_time_extrapolated[2];
  Ui[3] += dUi[3] + dUi_time_extrapolated[3];

  Uj[0] += dUj[0] + dUj_time_extrapolated[0];
  Uj[1] += dUj[1] + dUj_time_extrapolated[1];
  Uj[2] += dUj[2] + dUj_time_extrapolated[2];
  Uj[3] += dUj[3] + dUj_time_extrapolated[3];

  /* Check that we have physical masses and that we are not overshooting the
     particle's mass */
  chemistry_gradients_correct_unphysical_states(pi, pj, metal, cosmo, Ui, Uj);
}

#endif /* SWIFT_CHEMISTRY_GEAR_CHEMISTRY_HYPERBOLIC_GRADIENTS_H */
