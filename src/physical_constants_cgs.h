/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PHYSICAL_CONSTANTS_CGS_H
#define SWIFT_PHYSICAL_CONSTANTS_CGS_H

/**
 * @file physical_constants_cgs.h
 * @brief Physical constants in the CGS unit system.
 *
 * The constants declared in this file should _NOT_ be used directly
 * Users should use the converted values in the phys_const structure
 * where all the constants are defined in the system of units specified
 * in the parameter file.
 *
 * Of special interest if the value of G which is 0.03% different from the
 * default value adopted by Gadget-2 (G = 6.672e-8).
 *
 * All values are taken from C. Patrignani et al. (Particle Data Group), Chin.
 * Phys. C, 40, 100001 (2016) and 2017 update.
 * http://pdg.lbl.gov/2017/reviews/rpp2017-rev-phys-constants.pdf
 * http://pdg.lbl.gov/2017/reviews/rpp2017-rev-astrophysical-constants.pdf
 */

/*! Newton's gravitation constant [g^-1 cm^3 s^-2] */
const double const_newton_G_cgs = 6.67408e-8;

/*! Speed of light in vacuum [cm s^-1] */
const double const_speed_light_c_cgs = 2.99792458e10;

/*! Planck's constant [g cm^2 s^-1] */
const double const_planck_h_cgs = 6.626070040e-27;

/*! Planck's reduced constant [g cm^2 s^-1] */
const double const_planck_hbar_cgs = 1.054571800e-27;

/*! Boltzmann's constant [g cm^2 s^-2 K^-1] */
const double const_boltzmann_k_cgs = 1.38064852e-16;

/*! Avogadro number [-] */
const double const_avogadro_number_cgs = 6.022140857e23;

/*! Thomson cross-section [cm^2] */
const double const_thomson_cross_section_cgs = 6.6524587158e-25;

/*! Stefan-Boltzmann constant [g s^-3 K^-4] */
const double const_stefan_boltzmann_cgs = 5.670367e-5;

/*! Elementary charge [A s] */
const double const_electron_charge_cgs = 1.6021766208e-19;

/*! Electron-Volt [g cm^2 s^-2] */
const double const_electron_volt_cgs = 1.6021766208e-12;

/*! Mass of the electron [g] */
const double const_electron_mass_cgs = 9.10938356e-28;

/*! Mass of the proton [g] */
const double const_proton_mass_cgs = 1.672621898e-24;

/*! Tropical year [s] */
const double const_year_cgs = 3.15569252e7;

/*! Astronomical unit [cm] */
const double const_astronomical_unit_cgs = 1.49597870700e13;

/*! Parsec [cm] */
const double const_parsec_cgs = 3.08567758149e18;

/*! Light-year [cm] */
const double const_light_year_cgs = 9.46053e17;

/*! Mass of the Sun [g] */
const double const_solar_mass_cgs = 1.98848e33;

/*! Mass of the Earth [g] */
const double const_earth_mass_cgs = 5.9724e27;

/*! Temperature of the CMB at present day [K] */
const double const_T_CMB_0_cgs = 2.7255;

/*! Primordial Helium fraction [-] */
const double const_primordial_He_fraction_cgs = 0.245;

/*! Reduced Hubble constant units (i.e. H_0 / h == 100 km / s / Mpc in CGS)
 * [s^-1] */
const double const_reduced_hubble_cgs = 3.2407792894458e-18;

#endif /* SWIFT_PHYSICAL_CONSTANTS_CGS_H */
