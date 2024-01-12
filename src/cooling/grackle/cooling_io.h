
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_GRACKLE_IO_H
#define SWIFT_COOLING_GRACKLE_IO_H

/* Local includes */
#include "cooling_properties.h"
#include "cooling_struct.h"
#include "io_properties.h"
#include "physical_constants.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of cooling  to the file
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param cooling The #cooling_function_data
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, hid_t h_grp_columns,
    const struct cooling_function_data* cooling) {

#if COOLING_GRACKLE_MODE == 0
  io_write_attribute_s(h_grp, "Cooling Model", "Grackle");
#elif COOLING_GRACKLE_MODE == 1
  io_write_attribute_s(h_grp, "Cooling Model", "Grackle1");
#elif COOLING_GRACKLE_MODE == 2
  io_write_attribute_s(h_grp, "Cooling Model", "Grackle2");
#elif COOLING_GRACKLE_MODE == 3
  io_write_attribute_s(h_grp, "Cooling Model", "Grackle3");
#else
  error("This function should be called only with one of the Grackle cooling.");
#endif
}
#endif

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {

  int num = 0;

#if COOLING_GRACKLE_MODE >= 1
  /* List what we want to write */
  list[0] =
      io_make_output_field("HI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HI_frac, "HI mass fraction");

  list[1] =
      io_make_output_field("HII", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HII_frac, "HII mass fraction");

  list[2] =
      io_make_output_field("HeI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HeI_frac, "HeI mass fraction");

  list[3] =
      io_make_output_field("HeII", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HeII_frac, "HeII mass fraction");

  list[4] =
      io_make_output_field("HeIII", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HeIII_frac, "HeIII mass fraction");

  list[5] =
      io_make_output_field("e", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.e_frac, "free electron mass fraction");

  num += 6;
#endif

#if COOLING_GRACKLE_MODE >= 2
  list[6] =
      io_make_output_field("HM", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HM_frac, "H- mass fraction");

  list[7] =
      io_make_output_field("H2I", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.H2I_frac, "H2I mass fraction");

  list[8] =
      io_make_output_field("H2II", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.H2II_frac, "H2II mass fraction");

  num += 3;
#endif

#if COOLING_GRACKLE_MODE >= 3
  list[9] =
      io_make_output_field("DI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.DI_frac, "DI mass fraction");

  list[10] =
      io_make_output_field("DII", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.DII_frac, "DII mass fraction");

  list[11] =
      io_make_output_field("HDI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HDI_frac, "HDI mass fraction");
  num += 3;
#endif

  return num;
}

/**
 * @brief Parser the parameter file and initialize the #cooling_function_data
 *
 * @param parameter_file The parser parameter file
 * @param cooling The cooling properties to initialize
 * @param phys_const The #phys_const.
 */
__attribute__((always_inline)) INLINE static void cooling_read_parameters(
    struct swift_params* parameter_file, struct cooling_function_data* cooling,
    const struct phys_const* phys_const) {

  parser_get_param_string(parameter_file, "GrackleCooling:cloudy_table",
                          cooling->cloudy_table);

  cooling->primordial_chemistry = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:primordial_chemistry",
      COOLING_GRACKLE_MODE);

  if (cooling->primordial_chemistry < 0)
    error("Primordial chemistry cannot be below 0");

  if (cooling->primordial_chemistry > COOLING_GRACKLE_MODE)
    error("Cannot run primordial chemistry %i when compiled with %i",
          cooling->primordial_chemistry, COOLING_GRACKLE_MODE);

  cooling->H2_three_body_rate = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:H2_three_body_rate", 0);

  cooling->H2_cie_cooling = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:H2_cie_cooling", 0);

  cooling->H2_on_dust =
      parser_get_opt_param_int(parameter_file, "GrackleCooling:H2_on_dust", 0);

  cooling->local_dust_to_gas_ratio = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:local_dust_to_gas_ratio", -1);

  cooling->cmb_temperature_floor = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:cmb_temperature_floor", 1);

  cooling->with_uv_background =
      parser_get_param_int(parameter_file, "GrackleCooling:with_UV_background");

  cooling->redshift =
      parser_get_param_double(parameter_file, "GrackleCooling:redshift");

  cooling->with_metal_cooling =
      parser_get_param_int(parameter_file, "GrackleCooling:with_metal_cooling");

  cooling->use_radiative_transfer = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:use_radiative_transfer", 0);

  cooling->RT_heating_rate = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:RT_heating_rate_cgs", 0);

  cooling->RT_HI_ionization_rate = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:RT_HI_ionization_rate_cgs", 0);

  cooling->RT_HeI_ionization_rate = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:RT_HeI_ionization_rate_cgs", 0);

  cooling->RT_HeII_ionization_rate = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:RT_HeII_ionization_rate_cgs", 0);

  cooling->RT_H2_dissociation_rate = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:RT_H2_dissociation_rate_cgs", 0);

  cooling->volumetric_heating_rates = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:volumetric_heating_rates_cgs", 0);

  cooling->specific_heating_rates = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:specific_heating_rates_cgs", 0);

  cooling->HydrogenFractionByMass = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:HydrogenFractionByMass", 0.76);

  /* Self shielding */
  cooling->self_shielding_method = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:self_shielding_method", 0);

  if (cooling->self_shielding_method == -1) {
    cooling->self_shielding_threshold = parser_get_param_float(
        parameter_file, "GrackleCooling:self_shielding_threshold_atom_per_cm3");
  }

  cooling->HydrogenFractionByMass = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:HydrogenFractionByMass", 0.76);

  /* Initial step convergence */
  cooling->max_step = parser_get_opt_param_int(
      parameter_file, "GrackleCooling:max_steps", 10000);

  cooling->convergence_limit = parser_get_opt_param_double(
      parameter_file, "GrackleCooling:convergence_limit", 1e-2);

  /* Thermal time */
  cooling->thermal_time = parser_get_param_double(
      parameter_file, "GrackleCooling:thermal_time_myr");
  cooling->thermal_time *= phys_const->const_year * 1e6;
}

#endif /* SWIFT_COOLING_GRACKLE_IO_H */
