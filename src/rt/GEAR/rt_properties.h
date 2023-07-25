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
#ifndef SWIFT_RT_PROPERTIES_GEAR_H
#define SWIFT_RT_PROPERTIES_GEAR_H

#include "rt_grackle_utils.h"
#include "rt_interaction_rates.h"
#include "rt_parameters.h"
#include "rt_stellar_emission_model.h"

#include <string.h>

/**
 * @file src/rt/GEAR/rt_properties.h
 * @brief Main header file for the 'GEAR' radiative transfer scheme
 * properties.
 */

#define RT_IMPLEMENTATION "GEAR M1closure"

#if defined(RT_RIEMANN_SOLVER_GLF)
#define RT_RIEMANN_SOLVER_NAME "GLF Riemann Solver"
#elif defined(RT_RIEMANN_SOLVER_HLL)
#define RT_RIEMANN_SOLVER_NAME "HLL Riemann Solver"
#else
#error "No valid choice of RT Riemann solver has been selected"
#endif

/**
 * @brief Properties of the 'GEAR' radiative transfer model
 */
struct rt_props {

  /* Which stellar emission model to use */
  enum rt_stellar_emission_models stellar_emission_model;

  /* (Lower) frequency bin edges for photon groups */
  float photon_groups[RT_NGROUPS];

  /* Global constant stellar emission rates */
  double stellar_const_emission_rates[RT_NGROUPS];

  /* CFL condition */
  float CFL_condition;

  /* Factor to limit cooling time by */
  float f_limit_cooling_time;

  /* do we set initial ionization mass fractions manually? */
  int set_initial_ionization_mass_fractions;
  int set_equilibrium_initial_ionization_mass_fractions;

  /* initial mass fractions for ionization species */
  /* the following are required for manually setting exact values */
  float mass_fraction_HI_init;
  float mass_fraction_HII_init;
  float mass_fraction_HeI_init;
  float mass_fraction_HeII_init;
  float mass_fraction_HeIII_init;
  /* float number_density_electrons_init; [> todo: do we need this? <] */

  /* Hydrogen and Helium mass fractions of the non-metal portion of the gas */
  float hydrogen_mass_fraction;
  float helium_mass_fraction;

  /* Skip thermochemistry? For testing/debugging only! */
  int skip_thermochemistry;

  /* Re-do thermochemistry recursively if difference in internal energy is too
   * big? */
  int max_tchem_recursion;

  /* Optionally restrict maximal timestep for stars */
  float stars_max_timestep;

  /* Which stellar spectrum type to use? */
  int stellar_spectrum_type;
  /* If constant: get max frequency */
  double const_stellar_spectrum_max_frequency;
  /* If blackbody: get temperature */
  double stellar_spectrum_blackbody_T;

  /* Storage for integrated photoionization cross sections */
  /* Note: they are always in cgs. */
  double** energy_weighted_cross_sections;
  double** number_weighted_cross_sections;
  /* Mean photon energy in frequency bin for user provided spectrum. In erg.*/
  double average_photon_energy[RT_NGROUPS];
  /* Integral over photon numbers of user provided spectrum. */
  double photon_number_integral[RT_NGROUPS];

  /* Grackle Stuff */
  /* ------------- */

  /*! grackle unit system */
  code_units grackle_units;

  /*! grackle chemistry data */
  chemistry_data grackle_chemistry_data;

  /*! grackle chemistry data storage
   * (needed for local function calls) */
  chemistry_data_storage grackle_chemistry_rates;

  /*! use case B recombination? */
  int case_B_recombination;

  /*! make grackle talkative? */
  int grackle_verbose;

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* radiation emitted by stars this step. This is not really a property,
   * but a placeholder to sum up a global variable. It's being reset
   * every timestep. */
  unsigned long long debug_radiation_emitted_this_step;

  /* total radiation emitted by stars. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_emitted_tot;

  /* radiation absorbed by gas this step. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_absorbed_this_step;

  /* total radiation absorbed by gas. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_absorbed_tot;

  /* Max number of subcycles per hydro step */
  int debug_max_nr_subcycles;
#endif
};

/* Some declarations to avoid cyclical inclusions. */
/* ----------------------------------------------- */
/* Keep the declarations for *after* the definition of rt_props struct */

/**
 * @brief allocate and pre-compute the averaged cross sections
 * for each photon group and ionizing species.
 * Declare this here to avoid cyclical inclusions.
 *
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param us internal units struct
 **/
void rt_cross_sections_init(struct rt_props* restrict rt_props,
                            const struct phys_const* restrict phys_const,
                            const struct unit_system* restrict us);

/* Now for the good stuff                */
/* ------------------------------------- */

/**
 * @brief Print the RT model.
 *
 * @param rtp The #rt_props
 */
__attribute__((always_inline)) INLINE static void rt_props_print(
    const struct rt_props* rtp) {

  /* Only the master print */
  if (engine_rank != 0) return;

  message("Radiative transfer scheme: '%s'", RT_IMPLEMENTATION);
  message("RT Riemann Solver used: '%s'", RT_RIEMANN_SOLVER_NAME);
  char messagestring[200] = "Using photon frequency bins: [ ";
  char freqstring[20];
  for (int g = 0; g < RT_NGROUPS; g++) {
    sprintf(freqstring, "%.3g ", rtp->photon_groups[g]);
    strcat(messagestring, freqstring);
  }
  strcat(messagestring, "]");
  message("%s", messagestring);

  if (rtp->stellar_emission_model == rt_stellar_emission_model_const) {
    strcpy(messagestring, "Using constant stellar emission rates: [ ");
    for (int g = 0; g < RT_NGROUPS; g++) {
      sprintf(freqstring, "%.3g ", rtp->stellar_const_emission_rates[g]);
      strcat(messagestring, freqstring);
    }
    strcat(messagestring, "]");
    message("%s", messagestring);
  } else if (rtp->stellar_emission_model ==
             rt_stellar_emission_model_IlievTest) {
    message("Using Iliev+06 Test 4 stellar emission model.");
  } else {
    error("Unknown stellar emission model %d", rtp->stellar_emission_model);
  }

  if (rtp->set_equilibrium_initial_ionization_mass_fractions)
    message(
        "Setting initial ionization mass fractions "
        "assuming ionization equilibrium");
  if (rtp->set_initial_ionization_mass_fractions)
    message(
        "Using initial ionization mass fractions specified in parameter file");
  if (rtp->skip_thermochemistry)
    message("WARNING: Thermochemistry will be skipped.");
}

/**
 * @brief Initialize the global properties of the RT scheme.
 *
 * @param rtp The #rt_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void rt_props_init(
    struct rt_props* rtp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    struct cosmology* cosmo) {

  /* Make sure we reset debugging counters correctly after
   * zeroth step. */

  /* Read in photon frequency group properties */
  /* ----------------------------------------- */
  if (RT_NGROUPS <= 0) {
    error(
        "You need to run GEAR-RT with at least 1 photon group, "
        "you have %d",
        RT_NGROUPS);
  } else {
    /* !! Keep the frequencies in Hz for now. !! */
    parser_get_param_float_array(params, "GEARRT:photon_groups_Hz", RT_NGROUPS,
                                 rtp->photon_groups);
  }

  /* Sanity check: photon group edges must be in increasing order. */
  for (int g = 0; g < RT_NGROUPS - 1; g++) {
    if (rtp->photon_groups[g + 1] <= rtp->photon_groups[g])
      error(
          "Photon frequency bin edges need to be in increasing order. "
          "Found index %d %.3g <= %.3g",
          g + 1, rtp->photon_groups[g + 1], rtp->photon_groups[g]);
  }

  /* Get stellar emission rate model related parameters */
  /* -------------------------------------------------- */

  /* First initialize everything */
  for (int g = 0; g < RT_NGROUPS; g++) {
    rtp->stellar_const_emission_rates[g] = 0.;
  }

  rtp->stellar_emission_model = rt_stellar_emission_model_none;

  char stellar_model_str[80];
  parser_get_param_string(params, "GEARRT:stellar_luminosity_model",
                          stellar_model_str);

  if (strcmp(stellar_model_str, "const") == 0) {
    rtp->stellar_emission_model = rt_stellar_emission_model_const;
  } else if (strcmp(stellar_model_str, "IlievTest4") == 0) {
    rtp->stellar_emission_model = rt_stellar_emission_model_IlievTest;
  } else {
    error("Unknown stellar luminosity model '%s'", stellar_model_str);
  }

  if (rtp->stellar_emission_model == rt_stellar_emission_model_const) {
    /* Read the luminosities from the parameter file */
    double emission_rates[RT_NGROUPS];
    parser_get_param_double_array(params,
                                  "GEARRT:const_stellar_luminosities_LSol",
                                  RT_NGROUPS, emission_rates);
    const double unit_power = units_cgs_conversion_factor(us, UNIT_CONV_POWER);
    const double unit_power_inv = 1. / unit_power;
    for (int g = 0; g < RT_NGROUPS; g++) {
      rtp->stellar_const_emission_rates[g] = emission_rates[g] * unit_power_inv;
    }
  } else if (rtp->stellar_emission_model ==
             rt_stellar_emission_model_IlievTest) {
    /* Nothing to do here */
  } else {
    error("Unknown stellar emission model %d", rtp->stellar_emission_model);
  }

  /* get reduced speed of light factor */
  /* --------------------------------- */
  const float f_r = parser_get_param_float(params, "GEARRT:f_reduce_c");
  if (f_r <= 0.f)
    error("Invalid speed of light reduction factor: %.3e <= 0.", f_r);
  rt_params.reduced_speed_of_light = phys_const->const_speed_light_c * f_r;
  rt_params.reduced_speed_of_light_inverse =
      1.f / rt_params.reduced_speed_of_light;

  /* get CFL condition */
  /* ----------------- */
  const float CFL = parser_get_param_float(params, "GEARRT:CFL_condition");
  if (CFL <= 0.f) error("Invalid CFL number: %.3e <= 0.", CFL);
  rtp->CFL_condition = CFL;

  const float f_limit_cooling_time = parser_get_opt_param_float(
      params, "GEARRT:f_limit_cooling_time", /*default=*/0.6);
  if (f_limit_cooling_time < 0.f)
    error("Invalid cooling time reduction factor: %.3e < 0.",
          f_limit_cooling_time);
  else if (f_limit_cooling_time == 0.f)
    message("Warning: Computation of cooling time will be skipped");
  rtp->f_limit_cooling_time = f_limit_cooling_time;

  /* Get thermochemistry set-up */
  /* -------------------------- */
  rtp->hydrogen_mass_fraction =
      parser_get_param_float(params, "GEARRT:hydrogen_mass_fraction");
  rtp->helium_mass_fraction = 1.f - rtp->hydrogen_mass_fraction;
  if (rtp->hydrogen_mass_fraction <= 0.f || rtp->hydrogen_mass_fraction > 1.f)
    error("Invalid hydrogen mass fraction: %g", rtp->hydrogen_mass_fraction);

  /* Are we manually overwriting initial mass fractions of H and He? */
  rtp->set_initial_ionization_mass_fractions = parser_get_opt_param_int(
      params, "GEARRT:set_initial_ionization_mass_fractions",
      /* default = */ 0);
  if (rtp->set_initial_ionization_mass_fractions) {
    /* Read in mass fractions */
    rtp->mass_fraction_HI_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HI");
    rtp->mass_fraction_HII_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HII");
    rtp->mass_fraction_HeI_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HeI");
    rtp->mass_fraction_HeII_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HeII");
    rtp->mass_fraction_HeIII_init =
        parser_get_param_float(params, "GEARRT:mass_fraction_HeIII");

    /* Temporary check neglecting metals. Make sure we sum up to 1. */
    const float h_sum =
        rtp->mass_fraction_HI_init + rtp->mass_fraction_HII_init;
    if (fabsf(h_sum - rtp->hydrogen_mass_fraction) > 1e-4)
      error(
          "Inconsistent H mass fractions: XH_tot %.6g != XHI %.6g + XHII %.6g",
          rtp->hydrogen_mass_fraction, rtp->mass_fraction_HI_init,
          rtp->mass_fraction_HII_init);

    const float he_sum = rtp->mass_fraction_HeI_init +
                         rtp->mass_fraction_HeII_init +
                         rtp->mass_fraction_HeIII_init;
    if (fabsf(he_sum - rtp->helium_mass_fraction) > 1e-4)
      error(
          "Inconsistent He mass fractions: XHe_tot %.6g != XHeI %.6g + XHeII "
          "%.6g + XHeIII %.6g",
          rtp->helium_mass_fraction, rtp->mass_fraction_HeI_init,
          rtp->mass_fraction_HeII_init, rtp->mass_fraction_HeIII_init);

    const float mass_fraction_sum = h_sum + he_sum;
    if (fabsf(mass_fraction_sum - 1.f) > 1e-5)
      error("Constituent species mass fraction sums up to %.6f, I expect 1.0",
            mass_fraction_sum);
  } else {
    /* Initialize properties to deliberately bogus values */
    rtp->mass_fraction_HI_init = -1.f;
    rtp->mass_fraction_HII_init = -1.f;
    rtp->mass_fraction_HeI_init = -1.f;
    rtp->mass_fraction_HeII_init = -1.f;
    rtp->mass_fraction_HeIII_init = -1.f;
  }

  /* Are we setting up initial mass fractions in equilibrium? */
  rtp->set_equilibrium_initial_ionization_mass_fractions =
      parser_get_opt_param_int(
          params, "GEARRT:set_equilibrium_initial_ionization_mass_fractions",
          /* default = */ 0);

  if (rtp->set_equilibrium_initial_ionization_mass_fractions &&
      rtp->set_initial_ionization_mass_fractions)
    error(
        "Can't use equilibrium initial ionization mass fractions "
        "simultaneously with manually set mass fractions. Pick one.");

  /* Are we skipping thermochemistry? */
  rtp->skip_thermochemistry = parser_get_opt_param_int(
      params, "GEARRT:skip_thermochemistry", /* default = */ 0);

  /* Are we re-doing thermochemistry? */
  rtp->max_tchem_recursion = parser_get_opt_param_int(
      params, "GEARRT:max_tchem_recursion", /* default = */ 0);

  /* Stellar Spectra */
  /* --------------- */

  /* Initialize conditional parameters to bogus values */
  rtp->const_stellar_spectrum_max_frequency = -1.;
  rtp->stellar_spectrum_blackbody_T = -1.;

  rtp->stellar_spectrum_type =
      parser_get_param_int(params, "GEARRT:stellar_spectrum_type");
  if (rtp->stellar_spectrum_type == 0) {
    /* Constant spectrum: Read additional parameter */
    /* TODO: also translate back to internal units at later. For now, keep it in
     * Hz */
    rtp->const_stellar_spectrum_max_frequency = parser_get_param_float(
        params, "GEARRT:stellar_spectrum_const_max_frequency_Hz");
  } else if (rtp->stellar_spectrum_type == 1) {
    /* Blackbody spectrum: Read additional parameter */
    rtp->stellar_spectrum_blackbody_T = parser_get_param_float(
        params, "GEARRT:stellar_spectrum_blackbody_temperature_K");
    rtp->stellar_spectrum_blackbody_T /=
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  } else {
    error("Selected unknown stellar spectrum type %d",
          rtp->stellar_spectrum_type);
  }

  /* Maximal Star Timestep? */
  /* ---------------------- */
  rtp->stars_max_timestep = parser_get_opt_param_float(
      params, "GEARRT:stars_max_timestep", /*default=*/FLT_MAX);
  /* Turn off if negative value given */
  if (rtp->stars_max_timestep < 0.f) rtp->stars_max_timestep = FLT_MAX;
  /* Better safe than sorry */
  if (rtp->stars_max_timestep == 0.f)
    error("You are restricting star time step to 0. That's a no-no.");

#ifdef SWIFT_RT_DEBUG_CHECKS
  rtp->debug_radiation_emitted_tot = 0ULL;
  rtp->debug_radiation_emitted_this_step = 0ULL;

  rtp->debug_radiation_absorbed_tot = 0ULL;
  rtp->debug_radiation_absorbed_this_step = 0ULL;

  /* Don't make it an optional parameter here so we crash
   * if I forgot to provide it */
  rtp->debug_max_nr_subcycles =
      parser_get_param_int(params, "TimeIntegration:max_nr_rt_subcycles");
#endif

  /* Grackle setup */
  /* ------------- */
  rtp->grackle_verbose =
      parser_get_opt_param_int(params, "GEARRT:grackle_verbose", /*default=*/0);
  rtp->case_B_recombination = parser_get_opt_param_int(
      params, "GEARRT:case_B_recombination", /*default=*/1);
  rt_init_grackle(&rtp->grackle_units, &rtp->grackle_chemistry_data,
                  &rtp->grackle_chemistry_rates, rtp->hydrogen_mass_fraction,
                  rtp->grackle_verbose, rtp->case_B_recombination, us);

  /* Pre-compute interaction rates/cross sections */
  /* -------------------------------------------- */

  rtp->energy_weighted_cross_sections = NULL;
  rtp->number_weighted_cross_sections = NULL;
  rt_cross_sections_init(rtp, phys_const, us);

  /* Finishers */
  /* --------- */
}

/**
 * @brief Write an RT properties struct to the given FILE as a
 * stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_dump(
    const struct rt_props* props, FILE* stream) {

  restart_write_blocks((void*)props, sizeof(struct rt_props), 1, stream,
                       "RT props", "RT properties struct");
  /* The RT parameters, in particular the reduced speed of light, are
   * not defined at compile time. So we need to read them in again. */
  restart_write_blocks(&rt_params, sizeof(struct rt_parameters), 1, stream,
                       "RT global parameters", "RT global parameters struct");
}

/**
 * @brief Restore an RT properties struct from the given FILE as
 * a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 */
__attribute__((always_inline)) INLINE static void rt_struct_restore(
    struct rt_props* props, FILE* stream, const struct phys_const* phys_const,
    const struct unit_system* us) {

  restart_read_blocks((void*)props, sizeof(struct rt_props), 1, stream, NULL,
                      "RT properties struct");
  /* Set up stuff that needs array allocation */
  rt_init_grackle(&props->grackle_units, &props->grackle_chemistry_data,
                  &props->grackle_chemistry_rates,
                  props->hydrogen_mass_fraction, props->grackle_verbose,
                  props->case_B_recombination, us);

  props->energy_weighted_cross_sections = NULL;
  props->number_weighted_cross_sections = NULL;
  rt_cross_sections_init(props, phys_const, us);

  /* The RT parameters, in particular the reduced speed of light, are
   * not defined at compile time. So we need to write them down. */
  restart_read_blocks(&rt_params, sizeof(struct rt_parameters), 1, stream, NULL,
                      "RT global parameters struct");
}

#endif /* SWIFT_RT_PROPERTIES_GEAR_H */
