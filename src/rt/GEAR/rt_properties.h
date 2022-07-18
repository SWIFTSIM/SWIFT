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

#include "hydro.h" /* Need hydro_gamma */
#include "rt_parameters.h"

#include <grackle.h>
#include <string.h>

/**
 * @file src/rt/GEAR/rt_properties.h
 * @brief Main header file for the 'GEAR' radiative transfer scheme
 * properties.
 */

#define RT_IMPLEMENTATION "GEAR M1closure"

/**
 * @brief Properties of the 'GEAR' radiative transfer model
 */
struct rt_props {

  /* Are we using constant stellar emission rates? */
  int use_const_emission_rates;

  /* Frequency bin edges for photon groups
   * Includes 0 as leftmost edge, doesn't include infinity as
   * rightmost bin edge*/
  float photon_groups[RT_NGROUPS];

  /* Global constant stellar emission rates */
  double stellar_const_emission_rates[RT_NGROUPS];

  /* CFL condition */
  float CFL_condition;

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
  float number_density_electrons_init; /* todo: do we need this? */

  /* Hydrogen and Helium mass fractions of the non-metal portion of the gas */
  float hydrogen_mass_fraction;
  float helium_mass_fraction;

  /* Skip thermochemistry? For testing/debugging only! */
  int skip_thermochemistry;

  /* Optionally restrict maximal timestep for stars */
  float stars_max_timestep;

  /* Which stellar spectrum type to use? */
  int stellar_spectrum_type;
  /* If constant: get max frequency */
  double const_stellar_spectrum_max_frequency;
  /* If blackbody: get temperature */
  double stellar_spectrum_blackbody_T;

  /* Storage for integrated photoionization cross sections */
  double** energy_weighted_cross_sections;
  double** number_weighted_cross_sections;

  /* Grackle Stuff */
  /* ------------- */

  /*! grackle unit system */
  code_units grackle_units;

  /*! grackle chemistry data */
  chemistry_data grackle_chemistry_data;

  /* TODO: cleanup later with all other grackle stuff */
  /*! grackle chemistry data storage
   * (needed for local function calls) */
  /* chemistry_data_storage* grackle_chemistry_rates; */

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

  /* Total radiation energy in the gas. It's being reset every step. */
  float debug_total_radiation_conserved_energy[RT_NGROUPS];
  float debug_total_star_emitted_energy[RT_NGROUPS];

  /* Files to write energy budget to after every step */
  FILE* conserved_energy_filep;
  FILE* star_emitted_energy_filep;
#endif
};

/* Declare this here to avoid cyclical inclusions.
 * It needs to be after definition of rt_props, but
 * before the call to rt_props_init.*/
void rt_interaction_rates_init(struct rt_props* restrict rt_props,
                               const struct phys_const* restrict phys_const,
                               const struct unit_system* restrict us);

/**
 * @brief open up files to write some debugging check outputs.
 * This function is temporary for development, and shouldn't stay
 * long.
 *
 * @param rtp #rt_props struct
 * @param mode open files with this mode. "w" for new file, "a" for append.
 **/
#ifdef SWIFT_RT_DEBUG_CHECKS
static void rt_props_open_debugging_files(struct rt_props* rtp,
                                          const char* mode) {
#ifdef WITH_MPI
  return;
#endif

  rtp->conserved_energy_filep = fopen("RT_conserved_energy_budget.txt", mode);
  if (rtp->conserved_energy_filep == NULL)
    error("Couldn't open RT conserved energy budget file to write in");
  rtp->star_emitted_energy_filep = fopen("RT_star_injected_energy.txt", mode);
  if (rtp->star_emitted_energy_filep == NULL)
    error("Couldn't open RT star energy budget file to write in");

  if (strcmp(mode, "w") == 0 && rtp->use_const_emission_rates) {
    /* If we're starting a new file, dump the header first */
    FILE* files[2] = {rtp->conserved_energy_filep,
                      rtp->star_emitted_energy_filep};
    for (int f = 0; f < 2; f++) {
      fprintf(files[f], "# Emission rates: ");
      const double solar_luminosity = 3.826e33; /* erg/s */
      for (int g = 0; g < RT_NGROUPS; g++) {
        fprintf(files[f], "%12.6e ",
                rtp->stellar_const_emission_rates[g] * solar_luminosity);
      }
      fprintf(files[f], "\n");
    }
  }
};
#endif

/**
 * @brief initialize grackle during rt_props_init
 *
 * @param rtp #rt_props struct
 * @param us #unit_system struct
 **/
__attribute__((always_inline)) INLINE static void rt_props_init_grackle(
    struct rt_props* rtp, const struct unit_system* us) {

  /* TODO: cleanup later with all other grackle stuff */
  /* #ifdef SWIFT_RT_DEBUG_CHECKS */
  /*   grackle_verbose = 1; */
  /* #endif */

  /* Initialize units */
  /* ---------------- */
  /* we assume all quantities to be physical, not comoving */
  rtp->grackle_units.a_units = 1.0;
  rtp->grackle_units.a_value = 1.0;
  rtp->grackle_units.comoving_coordinates = 0;
  rtp->grackle_units.density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  rtp->grackle_units.length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  rtp->grackle_units.time_units =
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  rtp->grackle_units.velocity_units =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  /* Chemistry Parameters */
  /* -------------------- */
  /* More details on https://grackle.readthedocs.io/en/latest/Parameters.html */

  /* TODO: cleanup later with all other grackle stuff */
  /* rtp->grackle_chemistry_data = _set_default_chemistry_parameters(); */
  if (set_default_chemistry_parameters(&rtp->grackle_chemistry_data) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }
  /* chemistry on */
  rtp->grackle_chemistry_data.use_grackle = 1;
  /* cooling on */
  /* NOTE: without cooling on, it also won't heat... */
  rtp->grackle_chemistry_data.with_radiative_cooling = 1;
  /* 6 species atomic H and He. */
  rtp->grackle_chemistry_data.primordial_chemistry = 1;
  /* dust processes */
  rtp->grackle_chemistry_data.dust_chemistry = 0;
  /* H2 formation on dust */
  rtp->grackle_chemistry_data.h2_on_dust = 0;
  /* metal cooling (uses Cloudy) off (for now) */
  rtp->grackle_chemistry_data.metal_cooling = 0;
  /* no cooling below CMB temperature */
  rtp->grackle_chemistry_data.cmb_temperature_floor = 1;
  /* UV background off */
  rtp->grackle_chemistry_data.UVbackground = 0;
  /* data file - currently not used. */
  rtp->grackle_chemistry_data.grackle_data_file = "";
  /* adiabatic index */
  rtp->grackle_chemistry_data.Gamma = hydro_gamma;
  /* we'll provide grackle with ionization and heating rates from RT */
  rtp->grackle_chemistry_data.use_radiative_transfer = 1;
  /* fraction by mass of Hydrogen in the metal-free portion of the gas */
  rtp->grackle_chemistry_data.HydrogenFractionByMass =
      rtp->hydrogen_mass_fraction;

  /* TODO: cleanup later with all other grackle stuff */
  /* Initialize the chemistry_data_storage object to be
   * able to use local functions */
  /* chemistry_data_storage* chem_data_storage = */
  /*     malloc(sizeof(chemistry_data_storage)); */
  /* rtp->grackle_chemistry_rates = chem_data_storage; */
  /* if (!_initialize_chemistry_data(&rtp->grackle_chemistry_data, */
  /*                                 rtp->grackle_chemistry_rates, */
  /*                                 &rtp->grackle_units)) { */
  if (initialize_chemistry_data(&rtp->grackle_units) == 0) {
    error("Error in _initialize_chemistry_data");
  }
}

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
  char messagestring[200] = "Using photon frequency bins: [0.";
  char freqstring[20];
  for (int g = 1; g < RT_NGROUPS; g++) {
    sprintf(freqstring, ", %.3g", rtp->photon_groups[g]);
    strcat(messagestring, freqstring);
  }
  strcat(messagestring, "]");
  message("%s", messagestring);

  if (rtp->use_const_emission_rates) {
    strcpy(messagestring, "Using constant stellar emission rates: [ ");
    for (int g = 0; g < RT_NGROUPS; g++) {
      sprintf(freqstring, "%.3g ", rtp->stellar_const_emission_rates[g]);
      strcat(messagestring, freqstring);
    }
    strcat(messagestring, "]");
    message("%s", messagestring);
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
  } else if (RT_NGROUPS == 1) {
    rtp->photon_groups[0] = 0.f;
  } else {
    float frequencies[RT_NGROUPS - 1];
    /* !! Keep the frequencies in Hz for now. !! */
    parser_get_param_float_array(params, "GEARRT:photon_groups_Hz",
                                 RT_NGROUPS - 1, frequencies);
    for (int g = 0; g < RT_NGROUPS - 1; g++) {
      rtp->photon_groups[g + 1] = frequencies[g];
    }
    rtp->photon_groups[0] = 0.f;
  }

  /* Are we using constant emission rates? */
  /* ------------------------------------- */
  rtp->use_const_emission_rates = parser_get_opt_param_int(
      params, "GEARRT:use_const_emission_rates", /* default = */ 0);

  if (rtp->use_const_emission_rates) {
    double emission_rates[RT_NGROUPS];
    parser_get_param_double_array(params, "GEARRT:star_emission_rates_LSol",
                                  RT_NGROUPS, emission_rates);
    const double unit_power = units_cgs_conversion_factor(us, UNIT_CONV_POWER);
    const double unit_power_inv = 1. / unit_power;
    for (int g = 0; g < RT_NGROUPS; g++) {
      rtp->stellar_const_emission_rates[g] = emission_rates[g] * unit_power_inv;
    }
  } else {
    /* kill the run for now */
    error("GEAR-RT can't run without constant stellar emission rates for now.");
  }

  /* get reduced speed of light factor */
  /* --------------------------------- */
  const float f_r = parser_get_param_float(params, "GEARRT:f_reduce_c");
  rt_params.reduced_speed_of_light = phys_const->const_speed_light_c * f_r;
  rt_params.reduced_speed_of_light_inverse =
      1.f / rt_params.reduced_speed_of_light;

  /* get CFL condition */
  /* ----------------- */
  const float CFL = parser_get_param_float(params, "GEARRT:CFL_condition");
  rtp->CFL_condition = CFL;

  /* Get thermochemistry set-up */
  /* -------------------------- */
  rtp->hydrogen_mass_fraction =
      parser_get_param_float(params, "GEARRT:hydrogen_mass_fraction");
  rtp->helium_mass_fraction = 1.f - rtp->hydrogen_mass_fraction;

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

  for (int g = 0; g < RT_NGROUPS; g++)
    rtp->debug_total_star_emitted_energy[g] = 0.f;

  rt_props_open_debugging_files(rtp, "w");

#endif

  /* Grackle setup */
  /* ------------- */
  rt_props_init_grackle(rtp, us);

  /* Pre-compute interaction rates/cross sections */
  /* -------------------------------------------- */
  rt_interaction_rates_init(rtp, phys_const, us);

  /* Finishers */
  /* --------- */

  /* After initialisation, print params to screen */
  rt_props_print(rtp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Radiative transfer initialized");
  }
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
 */
__attribute__((always_inline)) INLINE static void rt_struct_restore(
    struct rt_props* props, FILE* stream) {

  restart_read_blocks((void*)props, sizeof(struct rt_props), 1, stream, NULL,
                      "RT properties struct");
  /* The RT parameters, in particular the reduced speed of light, are
   * not defined at compile time. So we need to write them down. */
  restart_read_blocks(&rt_params, sizeof(struct rt_parameters), 1, stream, NULL,
                      "RT global parameters struct");
#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Reset the file pointers for temporary stats */
  rt_props_open_debugging_files(props, "a");
#endif
}

#endif /* SWIFT_RT_PROPERTIES_GEAR_H */
