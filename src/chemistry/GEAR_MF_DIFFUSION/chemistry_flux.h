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
#ifndef SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_FLUX_H
#define SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_FLUX_H

#include "chemistry_getters.h"
#include "chemistry_properties.h"

/* Import the right file */
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
#include "hyperbolic/chemistry_flux.h"
#else
#include "parabolic/chemistry_flux.h"
#endif

/**
 * @file src/chemistry/GEAR_FVPM_DIFFUSION/chemistry_flux.h
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

  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  const double mZ_source = (metal_mass_interface > 0.0) ? mZi_0 : mZj_0;
  const double mZ_sink = (metal_mass_interface > 0.0) ? mZj_0 : mZi_0;
  const double mass_sink = (metal_mass_interface > 0.0) ? mj : mi;

  /* Source constaint: Do not allow to remove more mass, since you don't have
   * it... */
  if (mZ_source <= 0.0) {
    for (int k = 0; k < 4; k++) fluxes[k] = 0.0;
    return;
  }

  /* --- Capacity Constraint ---
   * Ensures the sink particle does not exceed Z = 1.0. This prevents
   * unphysical enrichment where metal mass exceeds gas mass, which would break
   * secondary chemistry/cooling solvers.
   */
  double capacity = mass_sink - mZ_sink;
  if (capacity <= 0.0) {
    for (int k = 0; k < 4; k++) fluxes[k] = 0.0;
    return;
  }

  /* --- Noise Gate ---
   * Truncates fluxes that are below machine epsilon relative to the source,
   * since it's just noise. Clipping it here can prevent the 'ratchet' effect.
   */
  const double eps = GEAR_FVPM_DIFF_NOISE_GATE;
  if (fabs(metal_mass_interface) < mZ_source * eps) {
    /* Set to 0 to kill noise */
    for (int k = 0; k < 4; k++) fluxes[k] = 0.0;
    return;
  }

  /* --- Stability & Startup Logic ---
   * We compute an 'effective' limit for the flux to prevent overshoots.
   *
   * Stability: Limits the flux to a fraction of the sink's current mass
   * to prevent the 'neighbor flooding' effect (oscillations).
   *
   * Startup: If the sink is pristine, the relative limit would be be zero,
   * stalling the PDE front. We use a 'startup_fraction' of the source mass to
   * kick-start diffusion into the vacuum.
   */
  const double safety_scale = GEAR_FVPM_DIFF_LIMITER_SAFETY;
  const double relative_change_limit = GEAR_FVPM_DIFF_LIMITER_SINK_STABILITY;
  const double startup_fraction = GEAR_FVPM_DIFF_LIMITER_STARTUP;
  double effective_limit_mass = min(mZ_source, capacity);

  /* Select the tighter of the source/capacity limit and the sink stability
   * limit */
  double sink_stability_limit =
      max(mZ_sink * relative_change_limit, mZ_source * startup_fraction);
  effective_limit_mass =
      min(effective_limit_mass, sink_stability_limit / safety_scale);

  /* --- Smooth Rational Limiter ---
   * Softly attenuates the flux using a rational factor: 1 / (1 + x).
   * This provides a C0-continuous transition toward the limit rather than
   * a hard discontinuous clip, which improves numerical convergence.
   */
  const double x = fabs(metal_mass_interface) /
                   (effective_limit_mass * safety_scale + 1e-40);
  const double factor = 1.0 / (1.0 + x);
  const double flux_init = fluxes[0];
  for (int k = 0; k < 4; k++) fluxes[k] *= factor;

  if (GEAR_FVPM_DIFF_FLUX_LIMITER_VERBOSITY > 0 && factor < 1e-1) {
    const int pi_is_active = pi->chemistry_data.flux.dt > 0.f;
    const int pj_is_active = pj->chemistry_data.flux.dt > 0.f;
    message(
        "[%lld, %lld | %d] Flux limiting, flux = %e, final_flux = %e, factor = "
        "%e,"
        " | mZi_0 = %e, mZj_0 = %e, kappa = %e, kappa = %e | mZ_source = %e, "
        " mZ_sink = %e, sink_capacity = %e | interaction_mode = %d | "
        " i_active = %d, j_active = %d",
        pi->id, pj->id, metal, flux_init * dt, fluxes[0] * dt, factor, mZi_0,
        mZj_0, chi->kappa, chj->kappa, mZ_source, mZ_sink, capacity,
        interaction_mode, pi_is_active, pj_is_active);
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_FLUX_H  */
