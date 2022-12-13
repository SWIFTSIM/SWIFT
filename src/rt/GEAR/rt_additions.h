/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GEAR_ADDITIONS_H
#define SWIFT_RT_GEAR_ADDITIONS_H

/**
 * @file src/rt/GEAR/rt_additions.h
 * @brief additional functions required for files outside of the RT modules
 * moved to separate file to avoid circular inclusions
 * */

/**
 * @brief update mass fluxes between two interacting particles during
 * hydro_iact_(non)sym(...) calls.
 *
 * @param pi first interacting particle
 * @param pj second interacting particle
 * @param mass_flux the mass flux between these two particles
 * @param mode 0: non-symmetric interaction, update i only. 1: symmetric
 * interaction.
 **/
__attribute__((always_inline)) INLINE static void rt_part_update_mass_fluxes(
    struct part* restrict pi, struct part* restrict pj, float mass_flux,
    int mode) {

  /* If a particle is losing mass, then it loses mass according to
   * its own mass fractions. If it's gaining mass, it's gaining mass
   * according to the interacting particle's mass fractions. */

  /* get the time step for the flux exchange. This is always the smallest time
     step among the two particles */
  const float mindt =
      (pj->flux.dt > 0.0f) ? fminf(pi->flux.dt, pj->flux.dt) : pi->flux.dt;

  const float mdt = mass_flux * mindt;

  /* Convention: negative flux for "left" particle pi */
  if (mass_flux < 0.f) {
    /* "left" particle is gaining mass */
    pi->rt_data.mass_flux.HI -= pj->rt_data.tchem.mass_fraction_HI * mdt;
    pi->rt_data.mass_flux.HII -= pj->rt_data.tchem.mass_fraction_HII * mdt;
    pi->rt_data.mass_flux.HeI -= pj->rt_data.tchem.mass_fraction_HeI * mdt;
    pi->rt_data.mass_flux.HeII -= pj->rt_data.tchem.mass_fraction_HeII * mdt;
    pi->rt_data.mass_flux.HeIII -= pj->rt_data.tchem.mass_fraction_HeIII * mdt;
  } else {
    /* "left" particle is losing mass */
    pi->rt_data.mass_flux.HI -= pi->rt_data.tchem.mass_fraction_HI * mdt;
    pi->rt_data.mass_flux.HII -= pi->rt_data.tchem.mass_fraction_HII * mdt;
    pi->rt_data.mass_flux.HeI -= pi->rt_data.tchem.mass_fraction_HeI * mdt;
    pi->rt_data.mass_flux.HeII -= pi->rt_data.tchem.mass_fraction_HeII * mdt;
    pi->rt_data.mass_flux.HeIII -= pi->rt_data.tchem.mass_fraction_HeIII * mdt;
  }

  if (mode == 1 || (pj->flux.dt < 0.f)) {
    /* Update right particle as well, even if it's inactive */

    if (mass_flux > 0.f) {
      /* "right" particle is gaining mass */
      pj->rt_data.mass_flux.HI += pi->rt_data.tchem.mass_fraction_HI * mdt;
      pj->rt_data.mass_flux.HII += pi->rt_data.tchem.mass_fraction_HII * mdt;
      pj->rt_data.mass_flux.HeI += pi->rt_data.tchem.mass_fraction_HeI * mdt;
      pj->rt_data.mass_flux.HeII += pi->rt_data.tchem.mass_fraction_HeII * mdt;
      pj->rt_data.mass_flux.HeIII +=
          pi->rt_data.tchem.mass_fraction_HeIII * mdt;
    } else {
      /* "right" particle is losing mass */
      pj->rt_data.mass_flux.HI += pj->rt_data.tchem.mass_fraction_HI * mdt;
      pj->rt_data.mass_flux.HII += pj->rt_data.tchem.mass_fraction_HII * mdt;
      pj->rt_data.mass_flux.HeI += pj->rt_data.tchem.mass_fraction_HeI * mdt;
      pj->rt_data.mass_flux.HeII += pj->rt_data.tchem.mass_fraction_HeII * mdt;
      pj->rt_data.mass_flux.HeIII +=
          pj->rt_data.tchem.mass_fraction_HeIII * mdt;
    }
  }
}

#endif /* SWIFT_RT_GEAR_ADDITIONS_H */
