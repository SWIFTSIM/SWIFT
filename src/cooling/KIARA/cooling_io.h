
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
#include "cooling.h"
#include "engine.h"
#include "io_properties.h"
#include "cooling_properties.h"

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
    const struct cooling_function_data *cooling) {

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

INLINE static void convert_part_HI_mass(const struct engine *e,
                                        const struct part *p,
                                        const struct xpart *xp, float *ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HI_frac;
}

INLINE static void convert_part_H2_mass(const struct engine *e,
                                        const struct part *p,
                                        const struct xpart *xp, float *ret) {

  float H2_frac = 0.;
  // const float X_H =
  // chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
#if COOLING_GRACKLE_MODE >= 2
  H2_frac = xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac;
  *ret = hydro_get_mass(p) * p->cooling_data.subgrid_fcold * H2_frac;
#else
  if (p->sf_data.SFR > 0) H2_frac = 1. - xp->cooling_data.HI_frac;
  *ret = hydro_get_mass(p) * H2_frac;
#endif
}

INLINE static void convert_part_HII_mass(const struct engine *e,
                                         const struct part *p,
                                         const struct xpart *xp, float *ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HII_frac;
}

INLINE static void convert_part_HeI_mass(const struct engine *e,
                                         const struct part *p,
                                         const struct xpart *xp, float *ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HeI_frac;
}

INLINE static void convert_part_HeII_mass(const struct engine *e,
                                          const struct part *p,
                                          const struct xpart *xp, float *ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HeII_frac;
}

INLINE static void convert_part_HeIII_mass(const struct engine *e,
                                           const struct part *p,
                                           const struct xpart *xp, float *ret) {

  *ret = hydro_get_mass(p) * xp->cooling_data.HeIII_frac;
}

INLINE static void convert_part_e_density(const struct engine *e,
                                          const struct part *p,
                                          const struct xpart *xp, float *ret) {

  *ret = (float)xp->cooling_data.e_frac;
}

INLINE static void convert_part_T(const struct engine *e, const struct part *p,
                                  const struct xpart *xp, float *ret) {

  const float u = hydro_get_physical_internal_energy(p, xp, e->cosmology);
  const float ne = xp->cooling_data.e_frac;
  *ret = cooling_convert_u_to_temp(u, ne, e->cooling_func, p);
}

#ifdef RT_NONE
INLINE static void convert_mass_fractions(const struct engine *engine,
                                          const struct part *part,
                                          const struct xpart *xpart,
                                          float *ret) {

  ret[0] = (float)xpart->cooling_data.HI_frac;
  ret[1] = (float)xpart->cooling_data.HII_frac;
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
    const struct part *parts, const struct xpart *xparts,
    struct io_props *list) {

  int num = 0;

#if COOLING_GRACKLE_MODE >= 1
  /* List what we want to write */
  list[num] = io_make_output_field_convert_part(
      "AtomicHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HI_mass, "Atomic hydrogen (HI) masses.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "IonizedHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HII_mass, "Ionized hydrogen (HII) masses.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "MolecularHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_H2_mass, "Molecular hydrogen (H2) masses.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "HeIMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HeII_mass, "HeI masses.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "HeIIMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HeII_mass, "HeII masses.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "HeIIIMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HeIII_mass, "HeIII masses.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "ElectronNumberDensities", FLOAT, 1, UNIT_CONV_NO_UNITS, -3.f,
      parts, xparts, convert_part_e_density, "Electron number densities"
      "in units of the hydrogen density.");
  num++;

  list[num] = io_make_output_field_convert_part(
      "Temperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts, xparts,
      convert_part_T, "Temperatures of the overall gas particles.");
  num++;

#if COOLING_GRACKLE_MODE >= 2
  list[num] = io_make_output_field(
      "SubgridTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts,
      cooling_data.subgrid_temp, "Temperatures of the cold phase"
      				 "of the subgrid ISM gas particles.");
  num++;

  list[num] =
      io_make_output_field("SubgridDensities", FLOAT, 1, UNIT_CONV_DENSITY,
                           -3.f, parts, cooling_data.subgrid_dens,
                           "Mass densities in physical units of the "
			   "subgrid ISM gas particles.");
  num++;

  list[num] = io_make_output_field(
      "SubgridColdISMFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      cooling_data.subgrid_fcold,
      "Fraction of gas particle masses in cold component of subgrid ISM.");
  num++;

  list[num] =
      io_make_output_field("DustMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
                           cooling_data.dust_mass, "Total masses in dust.");
  num++;

  list[num] = io_make_output_field("DustMassFractions", FLOAT,
                                   chemistry_element_count, UNIT_CONV_NO_UNITS,
                                   0.f, parts, cooling_data.dust_mass_fraction,
                                   "Fractions of the particles' masses that "
                                   "are in dust for a given element.");
  num++;

  list[num] =
      io_make_output_field("DustTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE,
                           0.f, parts, cooling_data.dust_temperature,
                           "Dust temperatures in subgrid ISM dust model.");
  num++;

  list[num] = io_make_output_field(
      "CoolingTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, parts,
      cooling_data.mixing_layer_cool_time,
      "Cooling times for the gas particle. If it's currently a firehose wind"
      "particle (decoupling_delay_time>0), this is the mixing layer cooling time.");
  num++;
#endif
#endif

#ifdef RT_NONE
  list[num] = io_make_output_field_convert_part(
      "IonMassFractions", FLOAT, 2, UNIT_CONV_NO_UNITS, 0, parts, xparts,
      convert_mass_fractions, "Mass fractions of all constituent species.");
  num++;
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
    struct swift_params *parameter_file, struct cooling_function_data *cooling,
    const struct phys_const *phys_const, const struct unit_system *us) {

  parser_get_param_string(parameter_file, "KIARACooling:cloudy_table",
                          cooling->cloudy_table);

  cooling->with_uv_background =
      parser_get_param_int(parameter_file, "KIARACooling:with_UV_background");

  cooling->redshift =
      parser_get_param_double(parameter_file, "KIARACooling:redshift");

  cooling->with_metal_cooling =
      parser_get_param_int(parameter_file, "KIARACooling:with_metal_cooling");

  cooling->provide_volumetric_heating_rates = parser_get_opt_param_int(
      parameter_file, "KIARACooling:provide_volumetric_heating_rates", -1);

  cooling->provide_specific_heating_rates = parser_get_opt_param_int(
      parameter_file, "KIARACooling:provide_specific_heating_rates", 1);

  /* Use lookup tables when outside ISM */
  cooling->use_tables_outside_ism = parser_get_opt_param_int(
      parameter_file, "KIARACooling:use_tables_outside_ism", 0);

  /* Self shielding */
  cooling->self_shielding_method = parser_get_opt_param_int(
      parameter_file, "KIARACooling:self_shielding_method", 3);

  /* What to do with adiabatic du/dt when in ISM mode */
  cooling->ism_adiabatic_heating_method = parser_get_opt_param_int(
      parameter_file, "KIARACooling:ism_adiabatic_heating_method", 1);

  /* Initial step convergence */
  cooling->max_step = parser_get_opt_param_int(
      parameter_file, "KIARACooling:grackle_max_steps", 500);

  cooling->timestep_accuracy = parser_get_opt_param_double(
      parameter_file, "KIARACooling:timestep_accuracy", 0.2);

  cooling->grackle_damping_interval = parser_get_opt_param_double(
      parameter_file, "KIARACooling:grackle_damping_interval", 5);

  cooling->thermal_time = parser_get_opt_param_double(
      parameter_file, "KIARACooling:thermal_time_myr", 0.);
  cooling->thermal_time *= phys_const->const_year * 1e6;

  /* flag to turn on dust evolution option, only works for GRACKLE_CHEMISTRY>=2
   * (KIARA) */
  cooling->use_grackle_dust_evol = parser_get_opt_param_int(
      parameter_file, "KIARACooling:use_grackle_dust_evol", 1);
#if COOLING_GRACKLE_MODE <= 1
  message("WARNING: Dust evol not implemented in SIMBA; use KIARA instead.");
  cooling->use_grackle_dust_evol = 0;
#endif

  /* These are dust parameters for KIARA's dust model (MODE>=2); irrelevant
   * otherwise */
  cooling->dust_destruction_eff = parser_get_opt_param_double(
      parameter_file, "KIARACooling:dust_destruction_eff", 0.3);

  cooling->dust_sne_coeff = parser_get_opt_param_double(
      parameter_file, "KIARACooling:dust_sne_coeff", 1.0);

  cooling->dust_sne_shockspeed = parser_get_opt_param_double(
      parameter_file, "KIARACooling:dust_sne_shockspeed", 100.0);

  cooling->dust_grainsize = parser_get_opt_param_double(
      parameter_file, "KIARACooling:dust_grainsize", 0.1);

  cooling->dust_growth_densref = parser_get_opt_param_double(
      parameter_file, "KIARACooling:dust_growth_densref", 2.3e-20);

  cooling->dust_growth_tauref = parser_get_opt_param_double(
      parameter_file, "KIARACooling:dust_growth_tauref", 1.0);

  cooling->cold_ISM_frac = parser_get_opt_param_double(
      parameter_file, "KIARACooling:cold_ISM_frac", 1.0);

  cooling->G0_computation_method = parser_get_opt_param_int(
      parameter_file, "KIARACooling:G0_computation_method", 3);

  cooling->G0_multiplier = parser_get_opt_param_double(
      parameter_file, "KIARACooling:G0_multiplier", 1.0);

  cooling->max_subgrid_density = parser_get_opt_param_double(
      parameter_file, "KIARACooling:max_subgrid_density_g_p_cm3", FLT_MAX);
  /* convert to internal units */
  cooling->max_subgrid_density /=
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  cooling->subgrid_threshold_n_H_inv = parser_get_opt_param_double(
      parameter_file, "KIARACooling:subgrid_threshold_n_H_cgs", 0.13);
  /* convert to internal units, take inverse to save compute time */
  cooling->subgrid_threshold_n_H_inv /=
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  cooling->subgrid_threshold_n_H_inv = 1.f / cooling->subgrid_threshold_n_H_inv;

  cooling->subgrid_threshold_T = parser_get_opt_param_double(
      parameter_file, "KIARACooling:subgrid_threshold_T_K", 1.e4);
  /* convert to internal units */
  cooling->subgrid_threshold_T /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  cooling->subgrid_warm_ism_EOS = parser_get_opt_param_double(
      parameter_file, "KIARACooling:subgrid_warm_ism_EOS", 0.f);

  cooling->entropy_floor_margin = parser_get_opt_param_double(
      parameter_file, "KIARACooling:entropy_floor_margin_dex", 1.0);
  cooling->entropy_floor_margin = pow(10.f, cooling->entropy_floor_margin);

  cooling->self_enrichment_metallicity = parser_get_opt_param_double(
      parameter_file, "KIARACooling:self_enrichment_metallicity", 0.f);

  cooling->do_cooling_in_rt = parser_get_opt_param_int(
      parameter_file, "KIARACooling:do_cooling_in_rt", 0);
}

#endif /* SWIFT_COOLING_KIARA_IO_H */
