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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H

#include "../chemistry_riemann_HLL.h"
#include "../chemistry_struct.h"
#include "../chemistry_unphysical.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/hyperbolic/chemistry_flux.h
 * @brief Header dealing with hyperbolic diffusion fluxes.
 */

/**
 * @brief Reset the riemann fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_part_reset_fluxes(
    struct part *restrict p) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    p->chemistry_data.metal_mass_riemann[i] = 0.0;
    p->chemistry_data.flux_riemann[i][0] = 0.0;
    p->chemistry_data.flux_riemann[i][1] = 0.0;
    p->chemistry_data.flux_riemann[i][2] = 0.0;
  }
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * The sign convention is that a positive total flux is subtracted from the
 * left state, and added to the right state, based on how we chose the unit
 * vector.
 *
 * @param p Particle.
 * @param metal Metal specie.
 * @param fluxes Fluxes accross the interface.
 * @param dt Time step for the flux exchange.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_fluxes_left(struct part *restrict p, const int metal,
                                  const double fluxes[4], const float dt) {
  p->chemistry_data.metal_mass_riemann[metal] -= fluxes[0] * dt;
  p->chemistry_data.flux_riemann[metal][0] -= fluxes[1] * dt;
  p->chemistry_data.flux_riemann[metal][1] -= fluxes[2] * dt;
  p->chemistry_data.flux_riemann[metal][2] -= fluxes[3] * dt;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the right of the interparticle interface.
 *
 * The sign convention is that a positive total flux is subtracted from the
 * left state, and added to the right state, based on how we chose the unit
 * vector.
 *
 * @param p Particle.
 * @param metal Metal specie.
 * @param fluxes Fluxes accross the interface.
 * @param dt Time step for the flux exchange.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_fluxes_right(struct part *restrict p, const int metal,
                                   const double fluxes[4], const float dt) {
  p->chemistry_data.metal_mass_riemann[metal] += fluxes[0] * dt;
  p->chemistry_data.flux_riemann[metal][0] += fluxes[1] * dt;
  p->chemistry_data.flux_riemann[metal][1] += fluxes[2] * dt;
  p->chemistry_data.flux_riemann[metal][2] += fluxes[3] * dt;
}

/**
 * @brief Time integrate the source term in the flux relaxation equation
 *
 * The flux relaxation equation is :
 *       dflux/dt + K/tau * grad Z = - flux/tau.
 *
 * The - flux/tau is the source term.
 *
 * @param p Particle.
 * @param metal Metal specie.
 * @param dt Time step for the flux integration (in physical units).
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_integrate_flux_source_term(
    struct part *restrict p, const int metal, const float dt,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {

  struct chemistry_part_data *chd = &p->chemistry_data;
  const double tau = chd->tau;
  const double exp_decay = exp(-dt / tau);
  const double one_minus_exp_decay = - expm1(-dt/tau);
  const double flux_current[3] = {chd->flux[metal][0], chd->flux[metal][1],
				  chd->flux[metal][2]};
  double flux_parabolic[3];
  chemistry_get_physical_parabolic_flux(p, metal, flux_parabolic, chem_data,
					cosmo);

  for (int i = 0; i < 3; i++) {
    if (tau != 0.0) {
      /* Note that parabolic flux already includes the minus term. */
      chd->flux[metal][i] =
	flux_current[i] * exp_decay + flux_parabolic[i] * one_minus_exp_decay;
    } else {
      /* The asymptotic solution is that Flux(t+Delta t) = Flux_parabolic(t) */
      chd->flux[metal][i] = flux_parabolic[i];
    }
  }
}

/**
 * @brief Compute the metal mass flux for the Riemann problem with the given
 * left and right state, and interface normal, surface area and velocity.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi The #part pi.
 * @param pj The #part pj.
 * @param metal Metal specie.
 * @param UL Left diffusion state variables (in physical units).
 * @param UR Right diffusion state variables (in physical units).
 * @param WL Left state variables (in physical units).
 * @param WR Right state variables (in physical units).
 * @param n_unit Unit vector of the interface.
 * @param Anorm Surface area of the interface (in physical units).
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 * @param fluxes (return) The resulting flux at the interface (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_compute_flux(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const int metal, const double UL[4],
    const double UR[4], const float WL[5], const float WR[5],
    const float n_unit[3], const float Anorm,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double fluxes[4]) {

  /* Predict the diffusion driver q at the cell interface */
  double qL, qR;
  if (chem_data->diffusion_mode == isotropic_constant) {
    /* For constant isotropic case, U = q = rho*Z.
       This is already predicted at the cell interface, nothing else to do. */
    qL = UL[0];
    qR = UR[0];
  } else {
    /* In these cases, U = rho*Z, q = Z */
    qL = chemistry_get_metal_mass_fraction(pi, metal);
    qR = chemistry_get_metal_mass_fraction(pj, metal);

    chemistry_gradients_predict_Z(pi, pj, metal, dx, cosmo, &qL, &qR);
  }

  /* (flux[3], q / tau * K) */
  double hyper_flux_L[4][3];
  chemistry_get_hyperbolic_flux(pi, metal, UL, qL, chem_data, cosmo, hyper_flux_L);
  double hyper_flux_R[4][3];
  chemistry_get_hyperbolic_flux(pj, metal, UR, qR, chem_data, cosmo, hyper_flux_R);

  chemistry_check_unphysical_hyperbolic_flux(hyper_flux_L);
  chemistry_check_unphysical_hyperbolic_flux(hyper_flux_R);

  chemistry_riemann_solve_for_flux(dx, pi, pj, UL, UR, WL, WR, hyper_flux_L,
				   hyper_flux_R, Anorm, n_unit, metal,
				   chem_data, cosmo, fluxes);

  /* Anorm is already in physical units here. */
  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H */
