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
 * Some constants can be overwritten by the YAML file values. If the
 * param argument is NULL, no overwriting is done.
 *
 * @param us The current internal system of units.
 * @param params The parsed parameter file.
 * @param internal_const The physical constants to initialize.
 */
void phys_const_init(const struct unit_system *us, struct swift_params *params,
                     struct phys_const *internal_const) {

  /* Units are declared as {U_M, U_L, U_t, U_I, U_T} */

  const float dimension_G[5] = {-1, 3, -2, 0, 0}; /* [g^-1 cm^3 s^-2] */
  internal_const->const_newton_G =
      const_newton_G_cgs / units_general_cgs_conversion_factor(us, dimension_G);

  /* Overwrite G if present in the file */
  if (params != NULL) {
    internal_const->const_newton_G = parser_get_opt_param_double(
        params, "PhysicalConstants:G", internal_const->const_newton_G);
  }

  const float dimension_c[5] = {0, 1, -1, 0, 0}; /* [cm s^-1] */
  internal_const->const_speed_light_c =
      const_speed_light_c_cgs /
      units_general_cgs_conversion_factor(us, dimension_c);

  const float dimension_h[5] = {1, 2, -1, 0, 0}; /* [g cm^2 s^-1] */
  internal_const->const_planck_h =
      const_planck_h_cgs / units_general_cgs_conversion_factor(us, dimension_h);
  internal_const->const_planck_hbar =
      const_planck_hbar_cgs /
      units_general_cgs_conversion_factor(us, dimension_h);

  const float dimension_k[5] = {1, 2, -2, 0, -1}; /* [g cm^2 s^-2 K^-1] */
  internal_const->const_boltzmann_k =
      const_boltzmann_k_cgs /
      units_general_cgs_conversion_factor(us, dimension_k);

  const float dimension_Na[5] = {0, 0, 0, 0, 0}; /* [ - ] */
  internal_const->const_avogadro_number =
      const_avogadro_number_cgs /
      units_general_cgs_conversion_factor(us, dimension_Na);

  const float dimension_thomson[5] = {0, 2, 0, 0, 0}; /* [cm^2] */
  internal_const->const_thomson_cross_section =
      const_thomson_cross_section_cgs /
      units_general_cgs_conversion_factor(us, dimension_thomson);

  const float dimension_stefan[5] = {1, 0, -3, 0, -4}; /* [g s^-3 K^-4] */
  internal_const->const_stefan_boltzmann =
      const_stefan_boltzmann_cgs /
      units_general_cgs_conversion_factor(us, dimension_stefan);

  const float dimension_ev[5] = {1, 2, -2, 0, 0}; /* [g cm^2 s^-2] */
  internal_const->const_electron_volt =
      const_electron_volt_cgs /
      units_general_cgs_conversion_factor(us, dimension_ev);

  const float dimension_charge[5] = {0, 0, 1, 1, 0}; /* [A s] */
  internal_const->const_electron_charge =
      const_electron_charge_cgs /
      units_general_cgs_conversion_factor(us, dimension_charge);

  const float dimension_mass[5] = {1, 0, 0, 0, 0}; /* [g] */
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

  const float dimension_time[5] = {0, 0, 1, 0, 0}; /* [s] */
  internal_const->const_year =
      const_year_cgs / units_general_cgs_conversion_factor(us, dimension_time);

  const float dimension_length[5] = {0, 1, 0, 0, 0}; /* [cm] */
  internal_const->const_astronomical_unit =
      const_astronomical_unit_cgs /
      units_general_cgs_conversion_factor(us, dimension_length);
  internal_const->const_parsec =
      const_parsec_cgs /
      units_general_cgs_conversion_factor(us, dimension_length);
  internal_const->const_light_year =
      const_light_year_cgs /
      units_general_cgs_conversion_factor(us, dimension_length);

  const float dimension_temperature[5] = {0, 0, 0, 0, 1}; /* [K] */
  internal_const->const_T_CMB_0 =
      const_T_CMB_0_cgs /
      units_general_cgs_conversion_factor(us, dimension_temperature);

  const float dimension_Yp[5] = {0, 0, 0, 0, 0}; /* [ - ] */
  internal_const->const_primordial_He_fraction =
      const_primordial_He_fraction_cgs /
      units_general_cgs_conversion_factor(us, dimension_Yp);

  const float dimension_reduced_hubble[5] = {0, 0, -1, 0, 0}; /* [s^-1] */
  internal_const->const_reduced_hubble =
      const_reduced_hubble_cgs /
      units_general_cgs_conversion_factor(us, dimension_reduced_hubble);
}

/**
 * @brief Print the value of the physical constants to stdout.
 *
 * @param internal_const The constants in the internal unit system.
 */
void phys_const_print(const struct phys_const *internal_const) {

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
  message("%25s = %e", "km/s/Mpc", internal_const->const_reduced_hubble);
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
