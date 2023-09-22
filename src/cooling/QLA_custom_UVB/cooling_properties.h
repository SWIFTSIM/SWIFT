/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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
#ifndef SWIFT_COOLING_PROPERTIES_CONST_LAMBDA_H
#define SWIFT_COOLING_PROPERTIES_CONST_LAMBDA_H

/**
 * @file src/cooling/const_lambda/cooling_properties.h
 * @brief Structures related to the "constant lambda" cooling function.
 *
 * This model assumes a constant cooling rate Lambda irrespective of redshift
 * or density.
 */

#define SMALLNUM 1.0e-60
#define COOLLIM 0.1
#define HEATLIM 20.0
#define eV_to_K 11606.0
#define eV_to_erg 1.60184e-12
#define MAX_TABLESIZE 500 /* Max # of lines in TREECOOL */

/* data for gas state */
typedef struct {
  double ne, necgs, nHcgs;
  double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
  double gJH0ne, gJHe0ne, gJHepne;
  double nH0, nHp, nHep, nHe0, nHepp;
  double XH, yhelium;
  double mhboltz;
  double ethmin; /* minimum internal energy for neutral gas */
  double mu;
} GasState;

/* tabulated rates */
typedef struct {
  double BetaH0, BetaHep, Betaff;
  double AlphaHp, AlphaHep, Alphad, AlphaHepp;
  double GammaeH0, GammaeHe0, GammaeHep;
} RateTable;

/* photo-ionization/heating rate table */
typedef struct {
  float variable;       /* logz for UVB */
  float gH0, gHe, gHep; /* photo-ionization rates */
  float eH0, eHe, eHep; /* photo-heating rates */
} PhotoTable;

/* current interpolated photo-ionization/heating rates */
typedef struct {
  char J_UV;
  double gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
} PhotoCurrent;

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling rate / nH^2 in physical cgs units [erg * s^-1 * cm^3] */
  double lambda_nH2_cgs;

  /*! Conversion factor from internal units to cgs for density */
  double conv_factor_density_to_cgs;

  /*! Conversion factor from internal units from cgs for internal energy
   * derivative */
  double conv_factor_energy_rate_from_cgs;

  /*! Conversion factor from internal units from cgs for internal energy
   * derivative */
  double conv_factor_energy_to_cgs;

  /*! Inverse of the proton mass in cgs units [g^-1] */
  double proton_mass_cgs_inv;

  /*! Constant multiplication factor for time-step criterion */
  float cooling_tstep_mult;

  /*! Use rapid cooling? */
  int rapid_cooling;

  GasState gs;      /*!< gas state */
  RateTable *RateT; /*!< tabulated rates */
  PhotoTable
      *PhotoTUVB;  /*!< photo-ionization/heating rate table for UV background */
  PhotoCurrent pc; /*!< current interpolated photo rates */
  int NheattabUVB; /*!< length of UVB photo table */
};

#endif /* SWIFT_COOLING_PROPERTIES_CONST_LAMBDA_H */
