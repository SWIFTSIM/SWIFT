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

#ifndef SWIFT_CHEMISTRY_EAGLE_ADDITIONS_H
#define SWIFT_CHEMISTRY_EAGLE_ADDITIONS_H

/**
 * @brief Resets the metal mass fluxes for schemes that use them.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void chemistry_reset_mass_fluxes(
    struct part* restrict p) {
  for (int i = 0; i < chemistry_element_count; i++) {
    p->chemistry_data.metal_mass_fluxes[i] = 0.f;
  }
  p->chemistry_data.metal_mass_flux_total = 0.f;
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

    /* update the metal mass fractions */

    /* First compute the current metal masses */
    const float current_mass_total = p->conserved.mass;
    const float current_metal_mass_total =
        current_mass_total * p->chemistry_data.metal_mass_fraction_total;
    float current_metal_masses[chemistry_element_count];
    for (int i = 0; i < chemistry_element_count; i++) {
      current_metal_masses[i] =
          current_mass_total * p->chemistry_data.metal_mass_fraction[i];
    }

    /* Add the mass fluxes */
    const float new_mass_total = current_mass_total + p->flux.mass;

    /* Check for vacuum? */
    if (new_mass_total <= 0.) {
      for (int i = 0; i < chemistry_element_count; i++) {
        p->chemistry_data.metal_mass_fraction[i] = 0.f;
        p->chemistry_data.smoothed_metal_mass_fraction[i] = 0.f;
      }
      p->chemistry_data.metal_mass_fraction_total = 0.f;
      chemistry_reset_mass_fluxes(p);
      /* Nothing left to do */
      return;
    }

    const float new_metal_mass_total =
        current_metal_mass_total + p->chemistry_data.metal_mass_flux_total;
    float new_metal_masses[chemistry_element_count];
    for (int i = 0; i < chemistry_element_count; i++) {
      new_metal_masses[i] = fmaxf(
          current_metal_masses[i] + p->chemistry_data.metal_mass_fluxes[i],
          0.f);
    }

    /* Finally update the metal mass ratios and reset the fluxes */
    const float new_mass_total_inv = 1.f / new_mass_total;
    p->chemistry_data.metal_mass_fraction_total =
        new_metal_mass_total * new_mass_total_inv;
    for (int i = 0; i < chemistry_element_count; i++) {
      p->chemistry_data.metal_mass_fraction[i] =
          new_metal_masses[i] * new_mass_total_inv;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (p->chemistry_data.metal_mass_fraction_total > 1.)
      error("Total metal mass fraction grew larger than 1!");
    float sum = 0.f;
    for (int i = 0; i < chemistry_element_count; i++) {
      sum += p->chemistry_data.metal_mass_fraction[i];
    }
    sum -= p->chemistry_data.metal_mass_fraction[chemistry_element_H];
    sum -= p->chemistry_data.metal_mass_fraction[chemistry_element_He];
    if (sum > p->chemistry_data.metal_mass_fraction_total)
      error(
          "Sum of element-wise metal mass fractions grew larger than total "
          "metal mass fraction!");
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

  const float mass_flux_integrated = mass_flux * flux_dt;

  /* Convention: a positive mass flux means that pi is losing said mass and pj
   * is gaining it. */
  if (mass_flux > 0.f) {
    /* pi is losing mass */
    pi->chemistry_data.metal_mass_flux_total -=
        mass_flux_integrated * pi->chemistry_data.metal_mass_fraction_total;
    for (int i = 0; i < chemistry_element_count; i++) {
      pi->chemistry_data.metal_mass_fluxes[i] -=
          mass_flux_integrated * pi->chemistry_data.metal_mass_fraction[i];
    }
  } else {
    /* pi is gaining mass: */
    pi->chemistry_data.metal_mass_flux_total -=
        mass_flux_integrated * pj->chemistry_data.metal_mass_fraction_total;
    for (int i = 0; i < chemistry_element_count; i++) {
      pi->chemistry_data.metal_mass_fluxes[i] -=
          mass_flux_integrated * pj->chemistry_data.metal_mass_fraction[i];
    }
  }

  /* update pj as well, even if it is inactive (flux.dt < 0.) */
  if (mode == 1 || pj->flux.dt < 0.f) {
    if (mass_flux > 0.f) {
      /* pj is gaining mass */
      pj->chemistry_data.metal_mass_flux_total +=
          mass_flux_integrated * pi->chemistry_data.metal_mass_fraction_total;
      for (int i = 0; i < chemistry_element_count; i++) {
        pj->chemistry_data.metal_mass_fluxes[i] +=
            mass_flux_integrated * pi->chemistry_data.metal_mass_fraction[i];
      }
    } else {
      /* pj is losing mass */
      pj->chemistry_data.metal_mass_flux_total +=
          mass_flux_integrated * pj->chemistry_data.metal_mass_fraction_total;
      for (int i = 0; i < chemistry_element_count; i++) {
        pj->chemistry_data.metal_mass_fluxes[i] +=
            mass_flux_integrated * pj->chemistry_data.metal_mass_fraction[i];
      }
    }
  }
}

#endif  // SWIFT_CHEMISTRY_EAGLE_ADDITIONS_H
