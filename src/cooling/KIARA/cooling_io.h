
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
#ifndef SWIFT_COOLING_KIARA_IO_H
#define SWIFT_COOLING_KIARA_IO_H

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
  io_write_attribute_s(h_grp, "Cooling Model", "Grackle0");
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

INLINE static void convert_part_HI_mass(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HI_frac;
}

INLINE static void convert_part_H2_mass(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  float H2_frac = 0.;
  //const float X_H = chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
#if COOLING_GRACKLE_MODE >= 2
  H2_frac = xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac;
  *ret = hydro_get_mass(p) * H2_frac;
#else
  if ( p->sf_data.SFR > 0 ) H2_frac = 1. - xp->cooling_data.HI_frac;
  *ret = hydro_get_mass(p) * H2_frac;
#endif
}

INLINE static void convert_part_HeII_mass(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HeII_frac;
}

INLINE static void convert_part_e_density(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, float* ret) {

  *ret = (float)xp->cooling_data.e_frac;
}

#if COOLING_GRACKLE_MODE >= 2
INLINE static void convert_part_dust_mass(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, float* ret) {

  *ret = (float)xp->cooling_data.dust_mass;
}

INLINE static void convert_part_dust_temperature(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, float* ret) {

  *ret = (float)xp->cooling_data.dust_temperature;
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
  list[num] = io_make_output_field_convert_part(
      "AtomicHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HI_mass, "Atomic hydrogen masses.");
  num ++;

  list[num] =
      io_make_output_field_convert_part(
      "MolecularHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_H2_mass, "Molecular hydrogen masses.");
  num ++;

  list[num] =
      io_make_output_field_convert_part(
      "HeIIMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HeII_mass, "HeII masses.");
  num ++;

  list[num] = io_make_output_field_convert_part(
      "ElectronNumberDensities", FLOAT, 1, UNIT_CONV_NUMBER_DENSITY, -3.f,
      parts, xparts, convert_part_e_density,
      "Electron number densities");
  num ++;


  /*
  list[0] =
      io_make_output_field("HI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           //cooling_data.HI_frac, "HI mass fraction");

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
#endif*/

#if COOLING_GRACKLE_MODE >= 2
  list[num] =
      io_make_output_field( "SubgridTemperatures", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      			   cooling_data.subgrid_temp, "Temperature of subgrid ISM in K");
  num ++;

  list[num] =
      io_make_output_field( "SubgridDensities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f, parts,
      			   cooling_data.subgrid_dens, "Comoving mass density of subgrid ISM");
  num ++;

  list[num] =
      io_make_output_field_convert_part("DustMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
                           convert_part_dust_mass, "Total mass in dust");
  num ++;

  list[num] =
      io_make_output_field_convert_part("DustTemperatures", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
                           convert_part_dust_temperature, "Dust temperature in subgrid dust model, in K");
  num ++;

  list[num] =
      io_make_output_field( "G0", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      			   chemistry_data.G0, "Interstellar radiation field in Habing units");
  num ++;
#endif

/*#if COOLING_GRACKLE_MODE >= 3
  list += num;

  list[0] =
      io_make_output_field("DI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.DI_frac, "DI mass fraction");

  list[1] =
      io_make_output_field("DII", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.DII_frac, "DII mass fraction");

  list[2] =
      io_make_output_field("HDI", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           cooling_data.HDI_frac, "HDI mass fraction");

  num += 3;
  */
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

  parser_get_param_string(parameter_file, "SIMBACooling:cloudy_table",
                          cooling->cloudy_table);

  cooling->with_uv_background =
      parser_get_param_int(parameter_file, "SIMBACooling:with_UV_background");

  cooling->redshift =
      parser_get_param_double(parameter_file, "SIMBACooling:redshift");

  cooling->with_metal_cooling =
      parser_get_param_int(parameter_file, "SIMBACooling:with_metal_cooling");

  cooling->provide_volumetric_heating_rates = parser_get_opt_param_int(
      parameter_file, "SIMBACooling:provide_volumetric_heating_rates", -1);

  cooling->provide_specific_heating_rates = parser_get_opt_param_int(
      parameter_file, "SIMBACooling:provide_specific_heating_rates", 1);

  /* Self shielding */
  cooling->self_shielding_method = parser_get_opt_param_int(
      parameter_file, "SIMBACooling:self_shielding_method", 3);

  /* Initial step convergence */
  cooling->max_step =
      parser_get_opt_param_int(parameter_file, "SIMBACooling:max_steps", 10000);

  cooling->convergence_limit = parser_get_opt_param_double(
      parameter_file, "SIMBACooling:convergence_limit", 1e-2);

  cooling->thermal_time =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:thermal_time_myr", 0.);
  cooling->thermal_time *= phys_const->const_year * 1e6;

  /* flag to turn on dust evolution option, only works for GRACKLE_CHEMISTRY>=2 (KIARA) */
  cooling->use_grackle_dust_evol =
      parser_get_opt_param_int(parameter_file, "SIMBACooling:use_grackle_dust_evol", 1);
#if COOLING_GRACKLE_MODE <= 1
  message("WARNING: Dust evol not implemented in SIMBA; use KIARA instead.");
  cooling->use_grackle_dust_evol = 0;
#endif

  /* These are dust parameters for KIARA's dust model (MODE>=2); irrelevant otherwise */
  cooling->dust_destruction_eff =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:dust_destruction_eff", 0.3);

  cooling->dust_sne_coeff =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:dust_sne_coeff", 1.0);

  cooling->dust_sne_shockspeed =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:dust_sne_shockspeed", 100.0);

  cooling->dust_grainsize =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:dust_grainsize", 0.1);

  cooling->dust_growth_densref =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:dust_growth_densref", 2.3e-20);

  cooling->dust_growth_tauref =
      parser_get_opt_param_double(parameter_file, "SIMBACooling:dust_growth_tauref", 1.0);

}

#endif /* SWIFT_COOLING_KIARA_IO_H */
