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

/* The constants declared in this file should _NOT_ be used directly */
/* Users should use the converted values in the phys_const structure */
/* where all the constants are defined in the system of units specified */
/* in the parameter file. */

/* All values are taken from K.A. Olive et al. (Particle Data Group), Chin. */
/* Phys. C, 38, 090001 (2014) and 2015 update. */
/* http://pdg.lbl.gov/2015/reviews/rpp2015-rev-phys-constants.pdf */
/* http://pdg.lbl.gov/2015/reviews/rpp2015-rev-astrophysical-constants.pdf */

/* Newton's gravitation constant */
const double const_newton_G_cgs = 6.67408e-8; /* g^-1 cm^3 s^-2 */

/* Speed of light in vacuum */
const double const_speed_light_c_cgs = 2.99792458e10; /* cm s^-1 */

/* Planck's constant */
const double const_planck_h_cgs = 6.626070040e-27; /* g cm^-2 s^-1 */

/* Planck's reduced constant */
const double const_planck_hbar_cgs = 1.054571800e-27; /* g cm^-2 s^-1 */

/* Boltzmann's constant */
const double const_boltzmann_k_cgs = 1.38064852e-16; /* g cm^2 s^-2 K^-1 */

/* Thomson cross-section */
const double const_thomson_cross_section_cgs = 6.6524587158e-25; /* cm^2 */

/* Elementary charge */
const double const_electron_charge_cgs = 1.6021766208e-19; /* A s^-1 */

/* Electron-Volt */
const double const_electron_volt_cgs = 1.6021766208e-12; /* g cm^2 s^-2 */

/* Mass of the electron */
const double const_electron_mass_cgs = 9.10938356e-28; /* g */

/* Mass of the proton */
const double const_proton_mass_cgs = 1.672621898e-24; /* g */

/* Tropical year */
const double const_year_cgs = 3.15569252e7; /* s */

/* Astronomical unit */
const double const_astronomical_unit_cgs = 1.49597870700e13; /* cm */

/* Parsec */
const double const_parsec_cgs = 3.08567758149e18; /* cm */

/* Light-year */
const double const_light_year_cgs = 9.46053e17; /* cm */

/* Mass of the Sun */
const double const_solar_mass_cgs = 1.9885e33; /* g */

/* Mass of the Earth */
const double const_earth_mass_cgs = 5.9726e27; /* g */

#endif /* SWIFT_PHYSICAL_CONSTANTS_CGS_H */
