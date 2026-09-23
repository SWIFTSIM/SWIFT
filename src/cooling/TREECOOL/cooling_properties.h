/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_PROPERTIES_TREECOOL_H
#define SWIFT_COOLING_PROPERTIES_TREECOOL_H

/**
 * @file src/cooling/TREECOOL/cooling_properties.h
 * @brief Structures related to the TREECOOL cooling function (Katz,
 * Weinberg & Hernquist 1996).
 */

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "parser.h"

/*! Number of entries in the tables of rate coefficients. */
#define treecool_cooling_N_temperature 2001

/*! Maximal number of redshift entries in the TREECOOL table. */
#define treecool_cooling_max_N_redshifts 512

/*! Number of columns (excluding the redshift one) in the TREECOOL table. */
#define treecool_cooling_N_columns 6

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /* ------------------------------------------------------------------ */
  /* Tables of rate coefficients, tabulated in log10(T / K).            */
  /* All the coefficients are in physical cgs units.                    */
  /* ------------------------------------------------------------------ */

  /*! Collisional excitation cooling of HI [erg * cm^3 * s^-1] */
  double table_Beta_H0[treecool_cooling_N_temperature];

  /*! Collisional excitation cooling of HeII [erg * cm^3 * s^-1] */
  double table_Beta_Hep[treecool_cooling_N_temperature];

  /*! Free-free (Bremsstrahlung) cooling [erg * cm^3 * s^-1] */
  double table_Beta_ff[treecool_cooling_N_temperature];

  /*! Recombination rate of HII [cm^3 * s^-1] */
  double table_Alpha_Hp[treecool_cooling_N_temperature];

  /*! Recombination rate of HeII [cm^3 * s^-1] */
  double table_Alpha_Hep[treecool_cooling_N_temperature];

  /*! Recombination rate of HeIII [cm^3 * s^-1] */
  double table_Alpha_Hepp[treecool_cooling_N_temperature];

  /*! Dielectronic recombination rate of HeII [cm^3 * s^-1] */
  double table_Alpha_d[treecool_cooling_N_temperature];

  /*! Collisional ionization rate of HI [cm^3 * s^-1] */
  double table_Gamma_eH0[treecool_cooling_N_temperature];

  /*! Collisional ionization rate of HeI [cm^3 * s^-1] */
  double table_Gamma_eHe0[treecool_cooling_N_temperature];

  /*! Collisional ionization rate of HeII [cm^3 * s^-1] */
  double table_Gamma_eHep[treecool_cooling_N_temperature];

  /*! Lowest log10(T / K) in the tables of rate coefficients */
  double log10_T_min;

  /*! Highest log10(T / K) in the tables of rate coefficients */
  double log10_T_max;

  /*! Spacing of the tables of rate coefficients in log10(T / K) */
  double delta_log10_T;

  /*! Inverse of the spacing of the tables in log10(T / K) */
  double inv_delta_log10_T;

  /* ------------------------------------------------------------------ */
  /* TREECOOL table giving the UV background as a function of redshift. */
  /* ------------------------------------------------------------------ */

  /*! Path to the TREECOOL file */
  char TREECOOL_file[PARSER_MAX_LINE_SIZE];

  /*! Number of redshift entries read from the TREECOOL file */
  int N_redshifts;

  /*! log10(1 + z) of the entries in the TREECOOL table */
  double TREECOOL_log10_1_plus_z[treecool_cooling_max_N_redshifts];

  /*! Photo-ionization rate of HI in the TREECOOL table [s^-1] */
  double TREECOOL_gamma_H0[treecool_cooling_max_N_redshifts];

  /*! Photo-ionization rate of HeI in the TREECOOL table [s^-1] */
  double TREECOOL_gamma_He0[treecool_cooling_max_N_redshifts];

  /*! Photo-ionization rate of HeII in the TREECOOL table [s^-1] */
  double TREECOOL_gamma_Hep[treecool_cooling_max_N_redshifts];

  /*! Photo-heating rate of HI in the TREECOOL table [erg * s^-1] */
  double TREECOOL_epsilon_H0[treecool_cooling_max_N_redshifts];

  /*! Photo-heating rate of HeI in the TREECOOL table [erg * s^-1] */
  double TREECOOL_epsilon_He0[treecool_cooling_max_N_redshifts];

  /*! Photo-heating rate of HeII in the TREECOOL table [erg * s^-1] */
  double TREECOOL_epsilon_Hep[treecool_cooling_max_N_redshifts];

  /* ------------------------------------------------------------------ */
  /* UV background interpolated to the current redshift.                */
  /* These are updated once per time-step by cooling_update().          */
  /* ------------------------------------------------------------------ */

  /*! Is the UV background switched on at the current redshift? */
  int UV_background_on;

  /*! Photo-ionization rate of HI at the current redshift [s^-1] */
  double gamma_H0;

  /*! Photo-ionization rate of HeI at the current redshift [s^-1] */
  double gamma_He0;

  /*! Photo-ionization rate of HeII at the current redshift [s^-1] */
  double gamma_Hep;

  /*! Photo-heating rate of HI at the current redshift [erg * s^-1] */
  double epsilon_H0;

  /*! Photo-heating rate of HeI at the current redshift [erg * s^-1] */
  double epsilon_He0;

  /*! Photo-heating rate of HeII at the current redshift [erg * s^-1] */
  double epsilon_Hep;

  /*! Temperature of the CMB at the current redshift [K] */
  double T_CMB;

  /*! (1 + z)^4 at the current redshift */
  double one_plus_z_to_the_4;

  /* ------------------------------------------------------------------ */
  /* Composition of the (primordial) gas.                               */
  /* ------------------------------------------------------------------ */

  /*! Hydrogen mass fraction */
  double X_H;

  /*! Helium mass fraction */
  double Y_He;

  /*! Ratio of the Helium to Hydrogen number densities (n_He / n_H) */
  double y_He;

  /* ------------------------------------------------------------------ */
  /* Model parameters.                                                  */
  /* ------------------------------------------------------------------ */

  /*! Switch to the rapid cooling regime above this value of dt / t_cool.
   * A negative value means that the slow-cooling regime is always used. */
  float rapid_cooling_threshold;

  /*! Redshift above which the UV background is switched off */
  float UV_background_start_redshift;

  /*! Include the inverse Compton cooling off the CMB? */
  int with_Compton_cooling;

  /* ------------------------------------------------------------------ */
  /* Conversion factors and constants.                                  */
  /* ------------------------------------------------------------------ */

  /*! Conversion factor from internal units to cgs for internal energy */
  double internal_energy_to_cgs;

  /*! Conversion factor from cgs to internal units for internal energy */
  double internal_energy_from_cgs;

  /*! Conversion factor from internal units to cgs for density */
  double density_to_cgs;

  /*! Conversion factor from internal units to cgs for time */
  double time_to_cgs;

  /*! Proton mass in cgs units [g] */
  double proton_mass_cgs;

  /*! Inverse of the proton mass in cgs units [g^-1] */
  double inv_proton_mass_cgs;

  /*! Boltzmann constant in cgs units [erg * K^-1] */
  double boltzmann_k_cgs;

  /*! Temperature of the CMB at redshift zero [K] */
  double T_CMB_0;

  /*! Minimal internal energy per unit mass in physical cgs units
   * [erg * g^-1] */
  double u_min_cgs;
};

#endif /* SWIFT_COOLING_PROPERTIES_TREECOOL_H */
