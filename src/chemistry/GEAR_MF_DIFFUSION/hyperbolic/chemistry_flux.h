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

#include "../chemistry_unphysical.h"
#include "../chemistry_struct.h"
#include "../chemistry_riemann_HLL.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/hyperbolic/chemistry_flux.h
 * @brief Header dealing with hyperbolic diffusion fluxes.
 */

/**
 * @brief Reset the metal mass fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_reset_fluxes(struct part* restrict p) {
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
chemistry_part_update_fluxes_left(struct part* restrict p,
				  const int metal,
				  const double fluxes[4],
				  const float dt) {
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
chemistry_part_update_fluxes_right(struct part* restrict p,
				   const int metal,
				   const double fluxes[4],
				   const float dt) {
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
    struct part* restrict p, const int metal, const float dt,
    const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo) {

  struct chemistry_part_data* chd = &p->chemistry_data;
  const double tau = chd->tau;
  const double exp_decay = exp(-dt / tau);
  const double one_minus_exp_decay = exp_decay + 1.0;

  const double flux_current[3] = {chd->flux[metal][0], chd->flux[metal][1],
				  chd->flux[metal][2]};
  double flux_parabolic[3];
  chemistry_get_physical_parabolic_flux(p, metal, flux_parabolic, chem_data,
					cosmo);

  for (int i = 0; i < 3; i++) {
    /* Note that parabolic flux already includes the minus term. */
    chd->flux[metal][i] = flux_current[i] * exp_decay + flux_parabolic[i]*one_minus_exp_decay;
  }
}

/**
 * @brief Compute the flux of the hyperbolic conservation law for a given
 * state U.
 *
 * @param p Particle.
 * @param U the state (metal density, diffusion flux) to use
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 * @param hypflux (return) resulting flux F(U) of the hyperbolic conservation
 * law (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_get_hyperbolic_flux(
    const struct part* restrict p, const int metal, const double U[4],
    const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo, double hypflux[4][3]) {

  /* Flux part (first row) */
  hypflux[0][0] = U[1];
  hypflux[0][1] = U[2];
  hypflux[0][2] = U[3];

  const double tau = p->chemistry_data.tau;
  double K[3][3];
  chemistry_get_physical_matrix_K(p, chem_data, cosmo, K);

  /* Get the right diffusion driver */
  double q;
  if (chem_data->diffusion_mode == isotropic_constant) {
    q = chemistry_get_physical_metal_density(p, metal, cosmo);
  } else {
    q = chemistry_get_metal_mass_fraction(p, metal);
  }

  const double multiplier = q / tau;

  /* The matrix part: q / tau * K */
  hypflux[1][0] = K[0][0] * multiplier;
  hypflux[1][1] = K[0][1] * multiplier;
  hypflux[1][2] = K[0][2] * multiplier;
  hypflux[2][0] = K[1][0] * multiplier;
  hypflux[2][1] = K[1][1] * multiplier;
  hypflux[2][2] = K[1][2] * multiplier;
  hypflux[3][0] = K[2][0] * multiplier;
  hypflux[3][1] = K[2][1] * multiplier;
  hypflux[3][2] = K[2][2] * multiplier;
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
 * @param fluxes (return) The resulting flux at the interface (of size 1).
 */
__attribute__((always_inline)) INLINE static void chemistry_compute_flux(
    const float dx[3], const struct part* restrict pi,
    const struct part* restrict pj, const int metal,
    const double UL[4], const double UR[4],
    const float WL[5], const float WR[5], const float n_unit[3],
    const float Anorm, const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo, double fluxes[4]) {

  /* (flux[3], q / tau * K) */
  double hyper_flux_L[4][3];
  chemistry_get_hyperbolic_flux(pi, metal, UL, chem_data, cosmo, hyper_flux_L);
  double hyper_flux_R[4][3];
  chemistry_get_hyperbolic_flux(pj, metal, UR, chem_data, cosmo, hyper_flux_R);

/* #ifdef SWIFT_RT_DEBUG_CHECKS */
  /* Check that the fluxes are meaningful */
  chemistry_check_unphysical_hyperbolic_flux(hyper_flux_L);
  chemistry_check_unphysical_hyperbolic_flux(hyper_flux_R);
/* #endif */

  chemistry_riemann_solve_for_flux(dx, pi, pj, UL, UR, WL, WR, hyper_flux_L,
				   hyper_flux_R, Anorm, n_unit, metal, chem_data,
				   cosmo, fluxes);

  /* Anorm is already in physical units here. */
  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H */
