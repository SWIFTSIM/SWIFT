/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GEAR_SPECIES_H
#define SWIFT_RT_GEAR_SPECIES_H

#define SWIFT_RT_GRACKLE_DUST

/**
 * @file src/rt/GEAR/rt_species.h
 * @brief header file concerning (ionizing) species.
 **/

enum rt_ionizing_species {
  rt_ionizing_species_HI = 0,
  rt_ionizing_species_HeI,
  rt_ionizing_species_HeII,
  rt_ionizing_species_count
};

/**
 * @file src/rt/GEAR/rt_species.h
 * @brief header file concerning species.
 **/

enum rt_species {
#if GEARRT_GRACKLE_MODE >= 1
  rt_species_HI = 0,
  rt_species_HII,
  rt_species_HeI,
  rt_species_HeII,
  rt_species_HeIII,
  rt_species_e,
#endif

#if GEARRT_GRACKLE_MODE >= 2
  rt_species_HM,
  rt_species_H2I,
  rt_species_H2II,
  #ifdef SWIFT_RT_GRACKLE_DUST
    rt_species_dust,
    rt_species_SNe_ThisTimeStep,
    rt_species_isrf_habing,
    rt_species_He_gas,
    rt_species_C_gas,
    rt_species_N_gas,
    rt_species_O_gas,
    rt_species_Ne_gas,
    rt_species_Mg_gas,
    rt_species_Si_gas,
    rt_species_S_gas,
    rt_species_Ca_gas,
    rt_species_Fe_gas,
    rt_species_He_dust,
    rt_species_C_dust,
    rt_species_N_dust,
    rt_species_O_dust,
    rt_species_Ne_dust,
    rt_species_Mg_dust,
    rt_species_Si_dust,
    rt_species_S_dust,
    rt_species_Ca_dust,
    rt_species_Fe_dust,
  #endif
#endif
  rt_species_count
};

/**
 * Get the ionizing energy in erg for all ionizing species.
 *
 * @param E_ion (return) the ionizing energies.
 **/
__attribute__((always_inline)) INLINE static void
rt_species_get_ionizing_energy(double E_ion[rt_ionizing_species_count]) {

  E_ion[rt_ionizing_species_HI] = 2.179e-11;   /* 13.60 eV in erg */
  E_ion[rt_ionizing_species_HeI] = 3.940e-11;  /* 24.59 eV in erg */
  E_ion[rt_ionizing_species_HeII] = 8.719e-11; /* 54.42 eV in erg */
}

#endif /* SWIFT_RT_GEAR_SPECIES_H */
