/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H

#include "chemistry_getters.h"
#include "chemistry_properties.h"

/* Import the right file */
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
#include "hyperbolic/chemistry_flux.h"
#else
#include "parabolic/chemistry_flux.h"
#endif

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry_flux.h
 * @brief Main header dealing with fluxes.
 *
 * */

/**
 * @brief Get the metal mass fluxes for the given particle.
 *
 * Notice that this function returns the solution to the Riemann problem. Hence
 * the flux is 1D.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @return flux Fluxes for the particle (array of size 1).
 */
__attribute__((always_inline)) INLINE double chemistry_get_metal_mass_fluxes(
    const struct part *restrict p, int metal) {
  return p->chemistry_data.flux.metal_mass[metal];
}

/**
 * @brief Reset the metal mass fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_reset_mass_fluxes(struct part *restrict p) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    p->chemistry_data.flux.metal_mass[i] = 0.0;
  }
}

/**
 * @brief Limit the metal mass flux to avoid negative metal masses.
 *
 * The flux limiter MUST be symmetric under pi <--> pj.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param F_diff (return) Array to write diffusion flux component into.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_limit_metal_mass_flux(const struct part *restrict pi,
                                const struct part *restrict pj, const int metal,
                                double fluxes[4], const float dt,
                                const int interaction_mode) {

  /* Convert the raw riemann mass derivative to mass */
  double metal_mass_interface = fluxes[0] * dt;
  if (metal_mass_interface == 0.0) return;

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get some convenient variables */
  const double mZi_0 = chi->metal_mass[metal];
  const double mZj_0 = chj->metal_mass[metal];
  const double upwind_mass = (metal_mass_interface > 0.0) ? mZi_0 : mZj_0;

  /* Do not allow to remove more mass, since you don't have it... */
  if (upwind_mass <= 0.0) {
    fluxes[0] = 0.0;
    fluxes[1] = 0.0;
    fluxes[2] = 0.0;
    fluxes[3] = 0.0;
    return;
  }

  const int pi_is_active = pi->chemistry_data.flux.dt > 0.f;
  const int pj_is_active = pj->chemistry_data.flux.dt > 0.f;

  /* The Noise Gate: If the flux is smaller than machine epsilon relative to the
   * mass, it's just noise. Clipping it here can prevent the 'ratchet' effect.
   */
  const double eps = 1e-15;
  if (fabs(metal_mass_interface) < upwind_mass * eps) {
    /* Set to 0 to kill noise */
    fluxes[0] = 0.0;
    fluxes[1] = 0.0;
    fluxes[2] = 0.0;
    fluxes[3] = 0.0;
    return;
  }

  /* Limit mass exchange to not overshoot too much. Smooth Rational Limiter :
     It behaves like a hard cut at 0.5 * mass, but follows a curve:
		     factor = 1 / (1 + |flux_mass| / source_mass) */
  const double safety_scale = 0.5;
  const double x = fabs(metal_mass_interface) / (upwind_mass * safety_scale);
  const double factor = 1.0 / (1.0 + x);
  const double flux_init = fluxes[0];
  fluxes[0] *= factor;
  fluxes[1] *= factor;
  fluxes[2] *= factor;
  fluxes[3] *= factor;

  if (GEAR_FVPM_DIFFUSION_FLUX_LIMITER_VERBOSITY > 0 && factor < 1e-1) {
    message(
	    "[%lld, %lld] Flux limiting, flux = %e, final_flux = %e, factor = %e,"
	    " mZi_r = %e, mZj_r = %e, upwind_mass = %e, mZi = %e, mZj = %e | mode = %d, i_active = %d, j_active = %d",
	    pi->id, pj->id, flux_init*dt, fluxes[0]*dt, factor, mZi_0, mZj_0, upwind_mass,
	    chi->metal_mass[metal], chj->metal_mass[metal], interaction_mode, pi_is_active, pj_is_active);
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H  */
