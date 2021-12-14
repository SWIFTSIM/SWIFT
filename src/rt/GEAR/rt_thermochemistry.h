/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_THERMOCHEMISTRY_GEAR_H
#define SWIFT_RT_THERMOCHEMISTRY_GEAR_H

#include "rt_ionization_equilibrium.h"

/**
 * @file src/rt/GEAR/rt_thermochemistry.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * thermochemistry related functions.
 */

/**
 * @brief initialize particle quantities relevant for the thermochemistry.
 *
 * @param p part to work with
 * @param rt_props rt_properties struct
 * @param phys_const physical constants struct
 * @param us unit system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void rt_tchem_first_init_part(
    struct part* restrict p, const struct rt_props* rt_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo) {

  if (rt_props->set_equilibrium_initial_ionization_mass_fractions) {
    float XHI, XHII, XHeI, XHeII, XHeIII;
    rt_ion_equil_get_mass_fractions(&XHI, &XHII, &XHeI, &XHeII, &XHeIII, p,
                                    rt_props, phys_const, us, cosmo);
    p->rt_data.tchem.mass_fraction_HI = XHI;
    p->rt_data.tchem.mass_fraction_HII = XHII;
    p->rt_data.tchem.mass_fraction_HeI = XHeI;
    p->rt_data.tchem.mass_fraction_HeII = XHeII;
    p->rt_data.tchem.mass_fraction_HeIII = XHeIII;
  } else if (rt_props->set_initial_ionization_mass_fractions) {
    p->rt_data.tchem.mass_fraction_HI = rt_props->mass_fraction_HI_init;
    p->rt_data.tchem.mass_fraction_HII = rt_props->mass_fraction_HII_init;
    p->rt_data.tchem.mass_fraction_HeI = rt_props->mass_fraction_HeI_init;
    p->rt_data.tchem.mass_fraction_HeII = rt_props->mass_fraction_HeII_init;
    p->rt_data.tchem.mass_fraction_HeIII = rt_props->mass_fraction_HeIII_init;
  }
}

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 */
__attribute__((always_inline)) INLINE static void rt_do_thermochemistry(
    struct part* restrict p, const struct rt_props* rt_props) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (!p->rt_data.debug_injection_done)
    error("Trying to do thermochemistry when injection step hasn't been done");
  if (!p->rt_data.debug_gradients_done)
    error("Trying to do thermochemistry when gradient step hasn't been done");
  if (!p->rt_data.debug_transport_done)
    error("Trying to do thermochemistry when transport step hasn't been done");

  p->rt_data.debug_thermochem_done += 1;
#endif

  if (rt_props->skip_thermochemistry) return;
}

#endif /* SWIFT_RT_THERMOCHEMISTRY_GEAR_H */
