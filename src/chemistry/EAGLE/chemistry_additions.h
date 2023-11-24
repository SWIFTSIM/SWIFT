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

#ifndef SWIFT_EAGLE_CHEMISTRY_ADDITIONS_H
#define SWIFT_EAGLE_CHEMISTRY_ADDITIONS_H

/**
 * @brief update metal mass fluxes between two interacting particles during
 * hydro_iact_(non)sym(...) calls.
 *
 * @param pi first interacting particle
 * @param pj second interacting particle
 * @param mass_flux the mass flux between these two particles.
 * @param flux_dt the time-step over which the fluxes are exchanged
 * @param mode 0: non-symmetric interaction, update i only. 1: symmetric
 * interaction.
 **/
__attribute__((always_inline)) INLINE static void runner_iact_chemistry_fluxes(
    struct part *restrict pi, struct part *restrict pj, float mass_flux,
    float flux_dt, int mode) {
#ifdef HYDRO_DOES_MASS_FLUX
  /* Metals are advected. I.e. a particle looses metals according to its own
   * metal mass fractions and gains mass according to the neighboring particle's
   * mass fractions. */

  const float mass_flux_integrated = mass_flux * flux_dt;

  /* Convention: a positive mass flux means that pi is loosing said mass and pj
   * is gaining it. */
  if (mass_flux > 0.f) {
    /* pi is loosing mass */
    for (int i = 0; i < chemistry_element_count; i++) {
      pi->chemistry_data.metal_mass_flux_total -=
          mass_flux_integrated * pi->chemistry_data.metal_mass_fraction_total;
      pi->chemistry_data.metal_mass_fluxes[i] -=
          mass_flux_integrated * pi->chemistry_data.metal_mass_fraction[i];
    }
  } else {
    /* pi is gaining mass: */
    for (int i = 0; i < chemistry_element_count; i++) {
      pi->chemistry_data.metal_mass_flux_total -=
          mass_flux_integrated * pj->chemistry_data.metal_mass_fraction_total;
      pi->chemistry_data.metal_mass_fluxes[i] -=
          mass_flux_integrated * pj->chemistry_data.metal_mass_fraction[i];
    }
  }

  /* update pj as well, even if it is inactive (flux.dt < 0.) */
  if (mode == 1 || pj->flux.dt < 0.f) {
    if (mass_flux > 0.f) {
      /* pj is gaining mass */
      for (int i = 0; i < chemistry_element_count; i++) {
        pj->chemistry_data.metal_mass_flux_total +=
            mass_flux_integrated * pi->chemistry_data.metal_mass_fraction_total;
        pj->chemistry_data.metal_mass_fluxes[i] +=
            mass_flux_integrated * pi->chemistry_data.metal_mass_fraction[i];
      }
    } else {
      /* pj is loosing mass */
      for (int i = 0; i < chemistry_element_count; i++) {
        pj->chemistry_data.metal_mass_flux_total +=
            mass_flux_integrated * pj->chemistry_data.metal_mass_fraction_total;
        pj->chemistry_data.metal_mass_fluxes[i] +=
            mass_flux_integrated * pj->chemistry_data.metal_mass_fraction[i];
      }
    }
  }
#endif
}

#endif  // SWIFT_EAGLE_CHEMISTRY_ADDITIONS_H
