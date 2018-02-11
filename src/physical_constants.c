/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "physical_constants.h"

/* Local headers. */
#include "error.h"
#include "physical_constants_cgs.h"
#include "restart.h"

/**
 * @brief Converts physical constants to the internal unit system
 *
 * @param us The current internal system of units.
 * @param internal_const The physical constants to initialize.
 */
void phys_const_init(struct unit_system *us,
                     struct phys_const *internal_const) {

  /* Units are declared as {U_M, U_L, U_t, U_I, U_T} */

  const float dimension_G[5] = {-1, 3, -2, 0, 0};
  internal_const->const_newton_G =
      const_newton_G_cgs / units_general_cgs_conversion_factor(us, dimension_G);

  const float dimension_c[5] = {0, 1, -1, 0, 0};
  internal_const->const_speed_light_c =
      const_speed_light_c_cgs /
      units_general_cgs_conversion_factor(us, dimension_c);

  const float dimension_h[5] = {1, -2, -1, 0, 0};
  internal_const->const_planck_h =
      const_planck_h_cgs / units_general_cgs_conversion_factor(us, dimension_h);
  internal_const->const_planck_hbar =
      const_planck_hbar_cgs /
      units_general_cgs_conversion_factor(us, dimension_h);

  const float dimension_k[5] = {1, 2, -2, 0, -1};
  internal_const->const_boltzmann_k =
      const_boltzmann_k_cgs /
      units_general_cgs_conversion_factor(us, dimension_k);

  const float dimension_thomson[5] = {0, 2, 0, 0, 0};
  internal_const->const_thomson_cross_section =
      const_thomson_cross_section_cgs /
      units_general_cgs_conversion_factor(us, dimension_thomson);

  const float dimension_ev[5] = {1, 2, -2, 0, 0};
  internal_const->const_electron_volt =
      const_electron_volt_cgs /
      units_general_cgs_conversion_factor(us, dimension_ev);

  const float dimension_charge[5] = {0, 0, -1, 1, 0};
  internal_const->const_electron_charge =
      const_electron_charge_cgs /
      units_general_cgs_conversion_factor(us, dimension_charge);

  const float dimension_mass[5] = {1, 0, 0, 0, 0};
  internal_const->const_electron_mass =
      const_electron_mass_cgs /
      units_general_cgs_conversion_factor(us, dimension_mass);
  internal_const->const_proton_mass =
      const_proton_mass_cgs /
      units_general_cgs_conversion_factor(us, dimension_mass);
  internal_const->const_solar_mass =
      const_solar_mass_cgs /
      units_general_cgs_conversion_factor(us, dimension_mass);
  internal_const->const_earth_mass =
      const_earth_mass_cgs /
      units_general_cgs_conversion_factor(us, dimension_mass);

  const float dimension_time[5] = {0, 0, 1, 0, 0};
  internal_const->const_year =
      const_year_cgs / units_general_cgs_conversion_factor(us, dimension_time);

  const float dimension_length[5] = {0, 1, 0, 0, 0};
  internal_const->const_astronomical_unit =
      const_astronomical_unit_cgs /
      units_general_cgs_conversion_factor(us, dimension_length);
  internal_const->const_parsec =
      const_parsec_cgs /
      units_general_cgs_conversion_factor(us, dimension_length);
  internal_const->const_light_year =
      const_light_year_cgs /
      units_general_cgs_conversion_factor(us, dimension_length);
}

void phys_const_print(struct phys_const *internal_const) {

  message("%25s = %e", "Gravitational constant",
          internal_const->const_newton_G);
  message("%25s = %e", "Speed of light", internal_const->const_speed_light_c);
  message("%25s = %e", "Planck constant", internal_const->const_planck_h);
  message("%25s = %e", "Boltzmann constant", internal_const->const_boltzmann_k);
  message("%25s = %e", "Thomson cross-section",
          internal_const->const_thomson_cross_section);
  message("%25s = %e", "Electron-Volt", internal_const->const_electron_volt);
  message("%25s = %e", "Year", internal_const->const_year);
  message("%25s = %e", "Astronomical Unit",
          internal_const->const_astronomical_unit);
  message("%25s = %e", "Parsec", internal_const->const_parsec);
  message("%25s = %e", "Solar mass", internal_const->const_solar_mass);
}

/**
 * @brief Write a phys_const struct to the given FILE as a stream of bytes.
 *
 * @param internal_const the struct
 * @param stream the file stream
 */
void phys_const_struct_dump(const struct phys_const *internal_const,
                            FILE *stream) {
  restart_write_blocks((void *)internal_const, sizeof(struct phys_const), 1,
                       stream, "physconst", "phys_const params");
}

/**
 * @brief Restore a phys_const struct from the given FILE as a stream of
 * bytes.
 *
 * @param internal_const the struct
 * @param stream the file stream
 */
void phys_const_struct_restore(const struct phys_const *internal_const,
                               FILE *stream) {
  restart_read_blocks((void *)internal_const, sizeof(struct phys_const), 1,
                      stream, NULL, "phys_const params");
}
