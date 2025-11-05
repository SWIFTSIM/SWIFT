/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PS2020_COOL_TABLES_H
#define SWIFT_PS2020_COOL_TABLES_H

/**
 * @file src/cooling/PS2020/cooling_tables.h
 * @brief PS2020 cooling function
 */

/* Config parameters. */
#include <config.h>

/*! Number of different bins along the temperature axis of the tables */
#define colibre_cooling_N_temperature 86

/*! Number of different bins along the redshift axis of the tables */
#define colibre_cooling_N_redshifts 46

/*! Number of different bins along the density axis of the tables */
#define colibre_cooling_N_density 71

/*! Number of different bins along the metallicity axis of the tables */
#define colibre_cooling_N_metallicity 11

/*! Number of different bins along the internal energy axis of the tables */
#define colibre_cooling_N_internalenergy 191

/*! Number of different cooling channels in the tables */
#define colibre_cooling_N_cooltypes 22

/*! Number of different heating channels in the tables */
#define colibre_cooling_N_heattypes 24

/*! Number of different electron fractions (each element - other atoms
 *  + tot prim + tot metal + tot)  in the tables */
#define colibre_cooling_N_electrontypes 14

/*! Number of different elements in the tables */
#define colibre_cooling_N_elementtypes 12

/**
 * @brief Elements present in the tables
 */
enum colibre_cooling_element {
  element_H,
  element_He,
  element_C,
  element_N,
  element_O,
  element_Ne,
  element_Mg,
  element_Si,
  element_S,
  element_Ca,
  element_Fe,
  element_OA
};

/**
 * @brief Hydrogen species
 */
enum colibre_hydrogen_species { neutral = 0, ionized = 1, molecular = 2 };

/**
 * @brief Cooling channels beyond the metal lines
 */
enum colibre_cooling_channels {
  cooltype_H2 = element_OA + 1,
  cooltype_molecules,
  cooltype_HD,
  cooltype_NetFFH,
  cooltype_NetFFM,
  cooltype_eeBrems,
  cooltype_Compton,
  cooltype_Dust
};

/**
 * @brief Heating channels beyond the metal lines
 */
enum colibre_heating_channels {
  heattype_H2 = element_OA + 1,
  heattype_COdiss,
  heattype_CosmicRay,
  heattype_UTA,
  heattype_line,
  heattype_Hlin,
  heattype_ChaT,
  heattype_HFF,
  heattype_Compton,
  heattype_Dust
};

/* Pre-declaration */
struct cooling_function_data;

void get_cooling_redshifts(struct cooling_function_data *cooling);
void read_cooling_header(struct cooling_function_data *cooling);
void read_cooling_tables(struct cooling_function_data *cooling);

#endif
