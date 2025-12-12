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
  return p->chemistry_data.metal_mass_riemann[metal];
}

/**
 * @brief Reset the metal mass fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_reset_mass_fluxes(struct part *restrict p) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    p->chemistry_data.metal_mass_riemann[i] = 0.0;
  }
}

/**
 * @brief Limit the metal mass flux to avoid negative metal masses.
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
                                double fluxes[4], const float dt) {
#ifdef GEAR_FVPM_DIFFUSION_FLUX_LIMITER_AGGRESSIVE_RESCALING
  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Convert the raw riemann mass derivative to mass */
  double metal_mass_interface = fluxes[0] * dt;

  /* Use the updated metal masses to ensure that the final result won't be
   * negative */
  const double mZi = chi->metal_mass[metal] + chi->metal_mass_riemann[metal];
  const double mZj = chj->metal_mass[metal] + chj->metal_mass_riemann[metal];

  /* This one seemed to work for a certain time */
  const double upwind_mass = (metal_mass_interface > 0.0) ? mZi : mZj;

  /* Choose upwind mass to determine a stability bound on the maximum allowed */
  /* mass exchange, (we do this to prevent negative masses under all */
  /* circumstances) */
  const double max_mass = 0.9 * fabs(upwind_mass);
  if (fabs(metal_mass_interface) > 0.0 &&
      fabs(metal_mass_interface) > max_mass) {
    const double factor = max_mass / fabs(metal_mass_interface);
    const double flux_init = fluxes[0];
    fluxes[0] *= factor;
    fluxes[1] *= factor;
    fluxes[2] *= factor;
    fluxes[3] *= factor;
    if (GEAR_FVPM_DIFFUSION_FLUX_LIMITER_VERBOSITY > 1) {
      message(
          "[%lld, %lld] Flux limiting, flux = %e, final_flux = %e, factor = %e,"
          " mZi_r = %e, mZj_r = %e, upwind_mass = %e, mZi = %e, mZj = %e",
          pi->id, pj->id, flux_init, fluxes[0], factor, mZi, mZj, upwind_mass,
          chi->metal_mass[metal], chj->metal_mass[metal]);
    }
  }
#endif /* GEAR_FVPM_DIFFUSION_FLUX_LIMITER_AGGRESSIVE */
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H  */
