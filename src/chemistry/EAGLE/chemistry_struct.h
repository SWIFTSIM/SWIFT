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
#ifndef SWIFT_CHEMISTRY_STRUCT_EAGLE_H
#define SWIFT_CHEMISTRY_STRUCT_EAGLE_H

#include <config.h>

/**
 * @brief The individual elements traced in the EAGLE model.
 */
enum chemistry_element {
  chemistry_element_H = 0,
  chemistry_element_He,
  chemistry_element_C,
  chemistry_element_N,
  chemistry_element_O,
  chemistry_element_Ne,
  chemistry_element_Mg,
  chemistry_element_Si,
  chemistry_element_Fe,
  chemistry_element_count
};

/**
 * @brief Global chemical abundance information in the EAGLE model.
 */
struct chemistry_global_data {

  /*! Fraction of the particle mass in given elements at the start of the run */
  float initial_metal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals at the start of the run */
  float initial_metal_mass_fraction_total;
};

/**
 * @brief Chemical abundances traced by the #part in the EAGLE model.
 */
struct chemistry_part_data {

  /*! Fraction of the particle mass in a given element */
  float metal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals */
  float metal_mass_fraction_total;

#ifdef HYDRO_DOES_MASS_FLUX
  /*! Mass fluxes of the metals in a given element */
  float metal_mass_fluxes[chemistry_element_count];

  /*! Mass flux in *all* metals */
  float metal_mass_flux_total;
#endif

  /*! Smoothed fraction of the particle mass in a given element */
  float smoothed_metal_mass_fraction[chemistry_element_count];

  /*! Smoothed fraction of the particle mass in *all* metals */
  float smoothed_metal_mass_fraction_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Fraction of total gas mass in metals coming from SNIa */
  float metal_mass_fraction_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Fraction of total gas mass in metals coming from AGB */
  float metal_mass_fraction_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Fraction of total gas mass in metals coming from SNII */
  float metal_mass_fraction_from_SNII;

  /*! Fraction of total gas mass in Iron coming from SNIa */
  float iron_mass_fraction_from_SNIa;

  /*! Smoothed fraction of total gas mass in Iron coming from SNIa */
  float smoothed_iron_mass_fraction_from_SNIa;
};

#define chemistry_spart_data chemistry_part_data

/**
 * @brief Chemical abundances traced by the #bpart in the EAGLE model.
 */
struct chemistry_bpart_data {

  /*! Mass in a given element */
  float metal_mass[chemistry_element_count];

  /*! Mass in *all* metals */
  float metal_mass_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Metal mass coming from SNIa */
  float metal_mass_from_SNIa;

  /*! Metal mass coming from AGB */
  float metal_mass_from_AGB;

  /*! Metal mass coming from SNII */
  float metal_mass_from_SNII;

  /*! Iron mass coming from SNIa */
  float iron_mass_from_SNIa;

  /*! Metallicity of converted part. */
  float formation_metallicity;

  /*! Smoothed metallicity of converted part. */
  float smoothed_formation_metallicity;
};

/**
 * @brief Chemical abundances traced by the #sink in the EAGLE model.
 *
 * Nothing here.
 */
struct chemistry_sink_data {};

#endif /* SWIFT_CHEMISTRY_STRUCT_EAGLE_H */
