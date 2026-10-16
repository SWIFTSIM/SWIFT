
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
#ifndef SWIFT_COOLING_GEAJAY_IO_H
#define SWIFT_COOLING_GEAJAY_IO_H

/* Local includes */
#include "adiabatic_index.h"
#include "cooling.h"
#include "cooling_properties.h"
#include "cooling_struct.h"
#include "engine.h"
#include "io_properties.h"
#include "physical_constants.h"
#include "units.h"

#include <float.h>
#include <geajay_message.h>
#include <geajay_options.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

#if COOLING_GEAJAY_MODE == 0
  io_write_attribute_s(h_grp, "Cooling Model", "geajay");
#elif COOLING_GEAJAY_MODE == 1
  io_write_attribute_s(h_grp, "Cooling Model", "geajay1");
#elif COOLING_GEAJAY_MODE == 2
  io_write_attribute_s(h_grp, "Cooling Model", "geajay2");
#elif COOLING_GEAJAY_MODE == 3
  io_write_attribute_s(h_grp, "Cooling Model", "geajay3");
#else
  error("This function should be called only with one of the geajay cooling.");
#endif

  /* Names of the components of CoolingRates */
  const int name_length = 32;
  const int num_channels = cooling_rate_channel_count();
  char names[GEAJAY_RATE_NCHANNELS_COSMIC_RAY][32];
  for (int k = 0; k < num_channels; k++)
    snprintf(names[k], name_length, "%s", cooling_rate_channel_name(k));
  hsize_t dims[1] = {num_channels};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, name_length);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "CoolingRates", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names[0]);
  H5Dclose(dset);
  H5Tclose(type);
  H5Sclose(space);
}
#endif

/**
 * @brief Cooling time of a particle for the snapshot.
 *
 * @param e The #engine.
 * @param p The #part.
 * @param xp The #xpart.
 * @param ret (return) The physical cooling time.
 */
INLINE static void convert_part_cooling_time(const struct engine *e,
                                             const struct part *p,
                                             const struct xpart *xp,
                                             double *ret) {
  *ret =
      cooling_time(e->physical_constants, e->internal_units,
                   e->hydro_properties, e->cosmology, e->cooling_func, p, xp);
}

/**
 * @brief Energy rates of a particle for the snapshot.
 *
 * @param e The #engine.
 * @param p The #part.
 * @param xp The #xpart.
 * @param ret (return) The physical rate of each channel.
 */
INLINE static void convert_part_cooling_rates(const struct engine *e,
                                              const struct part *p,
                                              const struct xpart *xp,
                                              double *ret) {
  cooling_energy_rates(e->physical_constants, e->internal_units,
                       e->hydro_properties, e->cosmology, e->cooling_func, p,
                       xp, ret);
}

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

#if COOLING_GEAJAY_MODE >= 1
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

#if COOLING_GEAJAY_MODE >= 2
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

#if COOLING_GEAJAY_MODE >= 3
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

  list[num] = io_make_physical_output_field_convert_part(
      "CoolingTimes", DOUBLE, 1, UNIT_CONV_TIME, 0.f, parts, xparts,
      /*can convert to comoving=*/1, convert_part_cooling_time,
      "Cooling times u / (du/dt) of the gas at its stored state, not a step "
      "average: u is the energy of the last kick, not the drifted "
      "InternalEnergies. Negative for net cooling, +inf at a net rate of "
      "exactly 0, NaN if geajay refused the particle; do not store with a "
      "lossy filter. Cosmological runs: the redshift is the one of the end of "
      "the step that contains the snapshot.");
  num++;

  list[num] = io_make_physical_output_field_convert_part(
      "CoolingRates", DOUBLE, cooling_rate_channel_count(),
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, -3.f * hydro_gamma_minus_one,
      parts, xparts,
      /*can convert to comoving=*/1, convert_part_cooling_rates,
      "Energy rates du/dt of the gas at the state of CoolingTimes, not a step "
      "average, one per process (names in SubgridScheme/NamedColumns/"
      "CoolingRates); "
      "positive heats. Their sum is the net rate of CoolingTimes up to double "
      "rounding, without the H2 formation heating that the solver adds in a "
      "step. NaN if geajay "
      "refused the particle. Cosmological runs: the redshift is the one of the "
      "end of the step that contains the snapshot.");
  num++;

#if COOLING_GEAJAY_MODE >= 1
  if (cooling_with_cosmic_rays()) {
    list[num] = io_make_output_field(
        "CosmicRayValues", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
        cooling_data.cosmic_ray_value,
        "Cosmic-ray value of the cooling, proper: the primary ionisation rate "
        "per H in s^-1 or the cosmic-ray energy density in eV cm^-3, as "
        "GeajayCooling:cosmic_ray_input says (not read with the input "
        "constant).");
    num++;

    list[num] = io_make_output_field(
        "CosmicRayEnergyLossRates", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
        cooling_data.cosmic_ray_energy_loss_rate,
        "Coulomb and hadronic energy loss rate of the cosmic rays in the last "
        "solve, eV cm^-3 s^-1, proper, averaged over the step; 0 unless the "
        "input is the energy density or when no solve ran, NaN after a failed "
        "solve.");
    num++;
  }
#endif

  return num;
}

/**
 * @brief Read a geajay choice parameter given by the name of its value or by
 * its integer.
 *
 * @param parameter_file The parser parameter file.
 * @param key The parameter, e.g. "GeajayCooling:h2_cie_cooling".
 * @param option The geajay option, e.g. "h2_cie_cooling".
 * @param def The default value.
 * @param swift_name A name handled by SWIFT, or NULL.
 * @param swift_value The value of swift_name.
 * @return The value.
 */
INLINE static int cooling_read_choice(struct swift_params *parameter_file,
                                      const char *key, const char *option,
                                      int def, const char *swift_name,
                                      int swift_value) {
  char text[PARSER_MAX_LINE_SIZE], def_text[32];
  snprintf(def_text, sizeof(def_text), "%d", def);
  parser_get_opt_param_string(parameter_file, key, text, def_text);

  char *end;
  const long number = strtol(text, &end, 10);
  if (end != text && *end == '\0') return (int)number;

  if (swift_name != NULL && strcmp(text, swift_name) == 0) return swift_value;

  int value;
  if (geajay_choice_from_string(option, text, &value) == GEAJAY_FAIL)
    error("%s: %s%s%s%s", key, geajay_last_error_message(),
          swift_name != NULL ? ", or \"" : "",
          swift_name != NULL ? swift_name : "", swift_name != NULL ? "\"" : "");
  return value;
}

/**
 * @brief Read the cosmic-ray options of geajay and the initial value of the
 * particles.
 *
 * Only the options present in the parameter file are passed to geajay, which
 * refuses a value that the chosen input does not read.
 *
 * @param parameter_file The parser parameter file.
 * @param cooling The cooling properties to initialize.
 */
INLINE static void cooling_read_cosmic_ray_parameters(
    struct swift_params *parameter_file,
    struct cooling_function_data *cooling) {

  static const char *options[] = {"cosmic_rays",
                                  "cosmic_ray_input",
                                  "cosmic_ray_rate_per_energy_density",
                                  "cosmic_ray_constant_rate",
                                  "cosmic_ray_heating",
                                  "cosmic_ray_heat_per_ionisation",
                                  "cosmic_ray_coulomb_loss_coefficient",
                                  "cosmic_ray_hadronic_loss_coefficient"};
  const int num_options = sizeof(options) / sizeof(options[0]);

  geajay_cosmic_ray_settings *settings = &cooling->cosmic_ray_settings;
  geajay_cosmic_ray_settings_init(settings);

  char key[PARSER_MAX_LINE_SIZE], text[PARSER_MAX_LINE_SIZE];
  for (int i = 0; i < num_options; i++) {
    snprintf(key, sizeof(key), "GeajayCooling:%s", options[i]);
    if (!parser_does_param_exist(parameter_file, key)) continue;
    parser_get_param_string(parameter_file, key, text);
    if (geajay_cosmic_ray_settings_set(settings, options[i], text) ==
        GEAJAY_FAIL)
      error("%s: geajay refused the value '%s': %s", key, text,
            geajay_last_error_message());
  }

  /* The initial value is SWIFT's: read only by a per-particle input */
  const char *initial_key = "GeajayCooling:cosmic_ray_initial_value";
  const int has_initial = parser_does_param_exist(parameter_file, initial_key);
  const int per_particle = settings->network != GEAJAY_COSMIC_RAYS_OFF &&
                           settings->input != GEAJAY_COSMIC_RAY_INPUT_CONSTANT;
  cooling->cosmic_ray_initial_value = 0.f;
  if (per_particle) {
    if (!has_initial)
      error("%s is required with the cosmic-ray input %s.", initial_key,
            geajay_cosmic_ray_choice_name("cosmic_ray_input", settings->input));
    const double value = parser_get_param_double(parameter_file, initial_key);
    if (!(value >= 0.) || value > FLT_MAX)
      error("%s must be >= 0 and fit in a float, got %e.", initial_key, value);
    cooling->cosmic_ray_initial_value = (float)value;
  } else if (has_initial) {
    error(
        "%s is read only with cosmic rays on and a per-particle input "
        "(ionisation_rate or energy_density).",
        initial_key);
  }

  if (settings->network == GEAJAY_COSMIC_RAYS_OFF) return;

  if (geajay_cosmic_ray_settings_check(settings) == GEAJAY_FAIL)
    error("geajay refused the cosmic-ray options: %s",
          geajay_last_error_message());

  if (cooling->primordial_chemistry < 1)
    error(
        "GeajayCooling:cosmic_rays needs a species network (H_He or above): "
        "GeajayCooling:primordial_chemistry_network is %d.",
        cooling->primordial_chemistry);
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

  parser_get_param_string(parameter_file, "GeajayCooling:cloudy_table_file",
                          cooling->cloudy_table);

  cooling->primordial_chemistry = cooling_read_choice(
      parameter_file, "GeajayCooling:primordial_chemistry_network",
      "species_network", COOLING_GEAJAY_MODE, NULL, 0);

  if (cooling->primordial_chemistry < 0)
    error("Primordial chemistry cannot be below 0");

  if (cooling->primordial_chemistry > COOLING_GEAJAY_MODE)
    error("Cannot run primordial chemistry %i when compiled with %i",
          cooling->primordial_chemistry, COOLING_GEAJAY_MODE);

  cooling->H2_three_body_rate =
      cooling_read_choice(parameter_file, "GeajayCooling:h2_three_body_fit",
                          "h2_three_body_fit", 0, NULL, 0);

  cooling->H2_cie_cooling =
      cooling_read_choice(parameter_file, "GeajayCooling:h2_cie_cooling",
                          "h2_cie_cooling", 0, NULL, 0);

  cooling->H2_on_dust = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:h2_formation_on_dust", 0);

  cooling->H2_photodissociation_heating = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:h2_photodissociation_heating", 0);
#if COOLING_GEAJAY_MODE < 2
  if (cooling->H2_photodissociation_heating)
    error(
        "GeajayCooling:h2_photodissociation_heating needs "
        "COOLING_GEAJAY_MODE >= 2 (H2 species tracked); this SWIFT was "
        "compiled with mode %d.",
        COOLING_GEAJAY_MODE);
#endif

  cooling->local_dust_to_gas_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:reference_dust_to_gas_ratio", -1);

  cooling->cmb_temperature_floor = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:metal_cooling_cmb_floor", 1);

  cooling->agora_cmb_temperature_floor = (char)parser_get_opt_param_int(
      parameter_file, "GeajayCooling:agora_cmb_temperature_floor", 0);

  cooling->with_uv_background =
      parser_get_param_int(parameter_file, "GeajayCooling:with_UV_background");

  cooling->redshift =
      parser_get_param_double(parameter_file, "GeajayCooling:redshift");

  cooling->with_metal_cooling =
      parser_get_param_int(parameter_file, "GeajayCooling:with_metal_cooling");

  cooling->use_radiative_transfer = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:use_radiative_transfer_fields", 0);

  cooling->RT_heating_rate = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:RT_heating_rate_cgs", 0);

  cooling->RT_HI_ionization_rate = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:RT_HI_ionization_rate_cgs", 0);

  cooling->RT_HeI_ionization_rate = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:RT_HeI_ionization_rate_cgs", 0);

  cooling->RT_HeII_ionization_rate = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:RT_HeII_ionization_rate_cgs", 0);

  cooling->RT_H2_dissociation_rate = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:RT_H2_dissociation_rate_cgs", 0);

  /* Debug/test-only: see cooling_function_data.disable_cooling_for_
     debugging's own doxygen. Never set in a production run. */
  cooling->disable_cooling_for_debugging = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:disable_cooling_for_debugging", 0);

  cooling->volumetric_heating_rates = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:volumetric_heating_rates_cgs", 0);

  cooling->specific_heating_rates = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:specific_heating_rates_cgs", 0);

  cooling->HydrogenFractionByMass = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:hydrogen_mass_fraction", -1.);

  cooling->initial_nHII_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nHII_to_nH_ratio", -1);

  cooling->initial_nHeI_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nHeI_to_nH_ratio", -1);

  cooling->initial_nHeII_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nHeII_to_nH_ratio", -1);

  cooling->initial_nHeIII_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nHeIII_to_nH_ratio", -1);

  cooling->initial_nDI_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nDI_to_nH_ratio", -1);

  cooling->initial_nDII_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nDII_to_nH_ratio", -1);

  cooling->initial_nHM_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nHM_to_nH_ratio", -1);

  cooling->initial_nH2I_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nH2I_to_nH_ratio", -1);

  cooling->initial_nH2II_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nH2II_to_nH_ratio", -1);

  cooling->initial_nHDI_to_nH_ratio = parser_get_opt_param_double(
      parameter_file, "GeajayCooling:initial_nHDI_to_nH_ratio", -1);

  /* Self shielding */
  cooling->self_shielding_method =
      cooling_read_choice(parameter_file, "GeajayCooling:uv_self_shielding",
                          "uv_self_shielding", 0, "density_threshold", -1);

  if (cooling->self_shielding_method == -1) {
    cooling->self_shielding_threshold = parser_get_param_float(
        parameter_file, "GeajayCooling:self_shielding_threshold_atom_per_cm3");
  }

  cooling->H2_self_shielding = cooling_read_choice(
      parameter_file, "GeajayCooling:h2_shielding_length_method",
      "h2_shielding_length_method", 0, NULL, 0);

  /* Mode 1 differences the density against the six neighbouring grid
     cells, which do not exist when geajay is called on one particle. */
  if (cooling->H2_self_shielding == 1)
    error(
        "GeajayCooling:h2_shielding_length_method = 1 (Sobolev) is not "
        "supported: its "
        "Sobolev-like length reads neighbouring grid cells, and SWIFT calls "
        "geajay one particle at a time. Use 2 or given_length (kernel "
        "support radius), or 3 or Jeans (local Jeans length).");
  if (cooling->H2_self_shielding < 0 || cooling->H2_self_shielding > 3)
    error("GeajayCooling:h2_shielding_length_method must be 0, 2 or 3, got %d.",
          cooling->H2_self_shielding);

  char H2_path[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(parameter_file,
                              "GeajayCooling:H2_self_shielding_path", H2_path,
                              "kernel_radius");
  if (strcmp(H2_path, "kernel_radius") == 0)
    cooling->H2_self_shielding_path_in_kernel_radii = 1.0f;
  else if (strcmp(H2_path, "kernel_diameter") == 0)
    cooling->H2_self_shielding_path_in_kernel_radii = 2.0f;
  else
    error(
        "GeajayCooling:H2_self_shielding_path must be kernel_radius or "
        "kernel_diameter, got '%s'.",
        H2_path);

  /* Sub-cycle iteration limit */
  cooling->max_step = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:max_substeps", 10000);
  if (cooling->max_step < 1)
    error("GeajayCooling:max_substeps must be >= 1, got %d.",
          cooling->max_step);

  /* Retries of a failed solve, each one halving the sub-step. Default 2:
   * two retries give 4x the iteration budget, which covered every stiff
   * cell measured (14670-25574 iterations against the 10000 default),
   * so stiff cells keep being integrated instead of skipping their
   * cooling. */
  cooling->subcycle_on_failure = parser_get_opt_param_int(
      parameter_file, "GeajayCooling:subcycle_on_failure", 2);
  if (cooling->subcycle_on_failure < 0 || cooling->subcycle_on_failure > 20)
    error("GeajayCooling:subcycle_on_failure must be in [0, 20], got %d.",
          cooling->subcycle_on_failure);

  /* Thermal time */
  cooling->thermal_time =
      parser_get_param_double(parameter_file, "GeajayCooling:thermal_time_myr");
  cooling->thermal_time *= phys_const->const_year * 1e6;

  char table_edge[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(parameter_file,
                              "GeajayCooling:table_density_edge", table_edge,
                              "extrapolate");
  if (geajay_table_edge_from_string(table_edge, &cooling->table_density_edge) ==
      GEAJAY_FAIL)
    error(
        "GeajayCooling:table_density_edge must be extrapolate or edge, got "
        "'%s'.",
        table_edge);

  cooling->proton_mass_cgs = phys_const->const_proton_mass *
                             units_cgs_conversion_factor(us, UNIT_CONV_MASS);

  cooling_read_cosmic_ray_parameters(parameter_file, cooling);
}

#endif /* SWIFT_COOLING_GEAJAY_IO_H */
