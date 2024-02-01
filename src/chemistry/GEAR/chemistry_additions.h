/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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

#ifndef SWIFT_CHEMISTRY_GEAR_ADDITIONS_H
#define SWIFT_CHEMISTRY_GEAR_ADDITIONS_H

/**
 * @brief Resets the metal mass fluxes for schemes that use them.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void chemistry_reset_mass_fluxes(
    struct part* restrict p) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.metal_mass_fluxes[i] = 0.f;
  }
}

/**
 * @brief Extra operations done during the kick. This needs to be
 * done before the particle mass is updated in the hydro_kick_extra.
 *
 * @param p Particle to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 */
__attribute__((always_inline)) INLINE static void chemistry_kick_extra(
    struct part* p, float dt_therm, float dt_grav, float dt_hydro,
    float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {
  /* For hydro schemes that exchange mass fluxes between the particles,
   * we want to advect the metals. */
  if (p->flux.dt > 0.) {

    /* Check for vacuum? */
    if (p->conserved.mass + p->flux.mass <= 0.) {
      for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
        p->chemistry_data.metal_mass[i] = 0.;
        p->chemistry_data.smoothed_metal_mass_fraction[i] = 0.;
      }
      chemistry_reset_mass_fluxes(p);
      /* Nothing left to do */
      return;
    }

    /* apply the metal mass fluxes and reset them */
    const double* metal_fluxes = p->chemistry_data.metal_mass_fluxes;
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      p->chemistry_data.metal_mass[i] =
          fmax(p->chemistry_data.metal_mass[i] + metal_fluxes[i], 0.);
    }

#ifdef SWIFT_DEBUG_CHECKS
    double sum = 0.;
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT - 1; i++) {
      sum += p->chemistry_data.metal_mass[i];
    }
    const double total_metal_mass =
        p->chemistry_data.metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT - 1];
    if (sum > total_metal_mass)
      error(
          "Sum of element-wise metal masses grew larger than total metal "
          "mass!");
    if (total_metal_mass > p->conserved.mass + p->flux.mass)
      error("Total metal mass grew larger than the particle mass!");
#endif
    chemistry_reset_mass_fluxes(p);
  }
}

/**
 * @brief update metal mass fluxes between two interacting particles during
 * hydro_iact_(non)sym(...) calls.
 *
 * Metals are advected. I.e. a particle loses metals according to its own
 * metal mass fractions and gains mass according to the neighboring particle's
 * mass fractions.
 *
 * @param pi first interacting particle
 * @param pj second interacting particle
 * @param mass_flux the mass flux between these two particles.
 * @param flux_dt the time-step over which the fluxes are exchanged
 * @param mode 0: non-symmetric interaction, update i only. 1: symmetric
 * interaction.
 **/
__attribute__((always_inline)) INLINE static void runner_iact_chemistry_fluxes(
    struct part* restrict pi, struct part* restrict pj, float mass_flux,
    float flux_dt, int mode) {

  const double mass_flux_integrated = mass_flux * flux_dt;
  const float mi = pi->conserved.mass;
  const float mi_inv = mi > 0 ? 1.f / mi : 0.f;
  const float mj = pj->conserved.mass;
  const float mj_inv = mj > 0 ? 1.f / mj : 0.f;

  /* Convention: a positive mass flux means that pi is losing said mass and pj
   * is gaining it. */
  if (mass_flux > 0.f) {
    /* pi is losing mass */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      pi->chemistry_data.metal_mass_fluxes[i] -=
          mass_flux_integrated * pi->chemistry_data.metal_mass[i] * mi_inv;
    }
  } else {
    /* pi is gaining mass: */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      pi->chemistry_data.metal_mass_fluxes[i] -=
          mass_flux_integrated * pj->chemistry_data.metal_mass[i] * mj_inv;
    }
  }

  /* update pj as well, even if it is inactive (flux.dt < 0.) */
  if (mode == 1 || pj->flux.dt < 0.f) {
    if (mass_flux > 0.f) {
      /* pj is gaining mass */
      for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
        pj->chemistry_data.metal_mass_fluxes[i] +=
            mass_flux_integrated * pi->chemistry_data.metal_mass[i] * mi_inv;
      }
    } else {
      /* pj is losing mass */
      for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
        pj->chemistry_data.metal_mass_fluxes[i] +=
            mass_flux_integrated * pj->chemistry_data.metal_mass[i] * mj_inv;
      }
    }
  }
}

#endif  // SWIFT_CHEMISTRY_GEAR_ADDITIONS_H
