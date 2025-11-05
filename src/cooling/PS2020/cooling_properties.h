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
#ifndef SWIFT_COOLING_PROPERTIES_PS2020_H
#define SWIFT_COOLING_PROPERTIES_PS2020_H

#define colibre_table_path_name_length 500

/**
 * @brief struct containing cooling tables
 */
struct cooling_tables {

  /* array of all mean particle masses mu (temperature) */
  float *Tmu;

  /* array of all mean particle masses mu (internal energy) */
  float *Umu;

  /* array of all cooling processes (temperature) */
  float *Tcooling;

  /* array of all cooling processes (internal energy) */
  float *Ucooling;

  /* array of all heating processes (temperature) */
  float *Theating;

  /* array of all heating processes (internal energy) */
  float *Uheating;

  /* array of all electron abundances (temperature) */
  float *Telectron_fraction;

  /* array of all electron abundances (internal energy) */
  float *Uelectron_fraction;

  /* array to get T from U */
  float *T_from_U;

  /* array to get U from T */
  float *U_from_T;

  /* array of equilibrium temperatures */
  float *logTeq;

  /* array of mean particle masses at equilibrium temperatures */
  float *meanpartmass_Teq;

  /* array of pressures at equilibrium temperatures */
  float *logPeq;

  /* array of hydrogen fractions at equilibrium temperature */
  float *logHfracs_Teq;

  /* array of all hydrogen fractions */
  float *logHfracs_all;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling tables */
  struct cooling_tables table;

  /*! Redshift bins */
  float *Redshifts;

  /*! Hydrogen number density bins */
  float *nH;

  /*! Temperature bins */
  float *Temp;

  /*! Metallicity bins */
  float *Metallicity;

  /*! Internal energy bins */
  float *Therm;

  /*! Abundance ratios for each metallicity bin and for each included element */
  float *LogAbundances;
  float *Abundances;
  float *Abundances_inv;

  /*! Atomic masses for all included elements */
  float *atomicmass;
  float *atomicmass_inv;

  /*! Mass fractions of all included elements */
  float *LogMassFractions;
  float *MassFractions;

  /*! Index for solar metallicity in the metallicity dimension */
  int indxZsol;

  /*! Solar metallicity (metal mass fraction) */
  float *Zsol;

  /*! Inverse of solar metallicity (metal mass fraction) */
  float *Zsol_inv;

  /*! Filepath to the directory containing the HDF5 cooling tables */
  char cooling_table_path[colibre_table_path_name_length];

  /* Distance from EOS to use thermal equilibrium temperature for subgrid props
   */
  float dlogT_EOS;

  /*! Redshift of H reionization */
  float H_reion_z;

  /*! H reionization energy in CGS units */
  float H_reion_heat_cgs;

  /*! Have we already done H reionization? */
  int H_reion_done;

  /*! Ca over Si abundance divided by the solar ratio for these elements */
  float Ca_over_Si_ratio_in_solar;

  /*! S over Si abundance divided by the solar ratio for these elements */
  float S_over_Si_ratio_in_solar;

  /*! Redshift of He reionization */
  float He_reion_z_centre;

  /*! Spread of the He reionization */
  float He_reion_z_sigma;

  /*! He reionization energy in CGS units */
  float He_reion_heat_cgs;

  /*! Internal energy conversion from internal units to CGS (for quick access)
   */
  double internal_energy_to_cgs;

  /*! Internal energy conversion from CGS to internal units (for quick access)
   */
  double internal_energy_from_cgs;

  /*! Pressure conversion from internal units to CGS (for quick access) */
  double pressure_to_cgs;

  /*! Number density conversion from internal units to CGS (for quick access) */
  double number_density_to_cgs;

  /*! Number density conversion from CGS to internal units (for quick access) */
  double number_density_from_cgs;

  /*! Density conversion from internal units to CGS (for quick access) */
  double density_to_cgs;

  /*! Density conversion from CGS to internal units (for quick access) */
  double density_from_cgs;

  /*! Inverse of proton mass in cgs (for quick access) */
  double inv_proton_mass_cgs;

  /*! Proton mass in cgs (for quick access) */
  double proton_mass_cgs;

  /*! Logarithm base 10 of the Boltzmann constant in CGS (for quick access) */
  double log10_kB_cgs;

  /*! Temperatur of the CMB at present day (for quick access) */
  double T_CMB_0;

  /*! Compton rate in cgs units */
  double compton_rate_cgs;

  /*! sigma_T * k_B / (m_e * c^2) in internal units */
  double y_compton_factor;

  /*! Minimal temperature allowed for the gas particles */
  double Tmin;

  /*! Minimal internal energy in cgs allowed for the gas particles */
  double umin_cgs;

  /*! Threshold to switch between rapid and slow cooling regimes. */
  double rapid_cooling_threshold;
};

/**
 * @brief Subgrid properties to calculate
 */
enum cooling_subgrid_properties {
  cooling_compute_subgrid_density,
  cooling_compute_subgrid_temperature,
  cooling_compute_subgrid_HI_fraction,
  cooling_compute_subgrid_HII_fraction,
  cooling_compute_subgrid_H2_fraction
};

#endif /* SWIFT_COOLING_PROPERTIES_PS2020_H */
