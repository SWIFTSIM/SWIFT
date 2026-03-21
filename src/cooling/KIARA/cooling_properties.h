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
#ifndef SWIFT_COOLING_PROPERTIES_KIARA_H
#define SWIFT_COOLING_PROPERTIES_KIARA_H

/* include grackle */
#include <grackle.h>

/**
 * @file src/cooling/grackle/cooling_properties.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Filename of the Cloudy Table */
  char cloudy_table[200];

  /*! Enable/Disable UV backgroud */
  int with_uv_background;

  /*! Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /*! unit system */
  code_units units;

  /*! Temperature of the CMB at present day (for quick access) */
  double T_CMB_0;

  /*! k_Boltz/m_p plus conversion factor for converting u<->T */
  double temp_to_u_factor;

  /*! Convert time to Myr */
  double time_to_Myr;

  /*! conversion unit factor for rate of change of thermal energy */
  double dudt_units;

  /*! grackle chemistry data */
  chemistry_data chemistry;

  /*! Enable/Disable metal cooling */
  int with_metal_cooling;

  /*! User provide volumetric heating rates */
  int provide_volumetric_heating_rates;

  /*! User provide specific heating rates */
  int provide_specific_heating_rates;

  /*! Self shielding method (1 -> 3 for grackle's ones, 0 for none and -1 for
   * GEAR) */
  int self_shielding_method;

  /*! What to do with adiabatic du/dt when in ISM mode: 0=no adiabatic heating,
   * 1=evolve in grackle, 2=use to evaporate cold ism */
  int ism_adiabatic_heating_method;

  /*! number of step max for first init */
  int max_step;

  /*! max fractional change in quantities in single grackle substep */
  double timestep_accuracy;

  /*! parameter to control how fast grackle damps oscillatory behaviour
   * (lower=more aggressive) */
  int grackle_damping_interval;

  /*! Duration for switching off cooling after an event (e.g. supernovae) */
  double thermal_time;

  /*! track dust growth and destruction (only available in KIARA) */
  int use_grackle_dust_evol;

  /*! track H2 formation; this is set within the code based on selection options
   */
  int use_grackle_h2_form;

  /*! G0 conversion factors, scales to MW value based on local/global galaxy
   * props */
  double G0_factor1;
  double G0_factor2;
  double G0_factorSNe;

  /*! Dust parameters; see sample yml file */
  double dust_destruction_eff;
  double dust_sne_coeff;
  double dust_sne_shockspeed;
  double dust_grainsize;
  double dust_growth_densref;
  double dust_growth_tauref;

  /*! For dust model, need self-enrichment up to a small metallicity to
   * kick-start dust */
  double self_enrichment_metallicity;

  /*! For subgrid model (eg KIARA) need a subgrid ISM fraction */
  double cold_ISM_frac;

  /*! For Grackle subgrid model, choose way to determine G0: 1=Local SFR density; 
   * 2=Global sSFR; 3=2 if sSFR != 0, else 1; -3: vice versa */
  int G0_computation_method;

  /*! For Grackle subgrid model, arbitrary multiplier for G0 */
  double G0_multiplier;

  /*! For Grackle subgrid model, set max density to avoid pointlessly
   * over-iterating in Grackle */
  double max_subgrid_density;

  /*! For Grackle subgrid model, inverse of threshold nH above which multi-phase
   * ISM model kicks in */
  double subgrid_threshold_n_H_inv;

  /*! For Grackle subgrid model, temperature at threshold nH */
  double subgrid_threshold_T;

  /*! For Grackle subgrid model, Power-law eqn of state for warm ISM component
   * above threshold n_H */
  double subgrid_warm_ism_EOS;

  /*! For Grackle subgrid model, factor above entropy floor allowed to be in
   * subgrid mode */
  double entropy_floor_margin;

  /*! Option to use Cloudy lookup tables when outside ISM */
  int use_tables_outside_ism;

  /*! When using radiative transfer, set this on if you want thermochemistry 
   * done within RT modules (eg in KIARART) */
  int do_cooling_in_rt;
};

#endif /* SWIFT_COOLING_PROPERTIES_KIARA_H */
