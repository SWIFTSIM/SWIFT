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

/* physical constants in cgs */
#define NEWTON_GRAVITY_CGS 6.672e-8f
#define SOLAR_MASS_IN_CGS 1.989e33f
#define PARSEC_IN_CGS 3.086e18f
#define PROTON_MASS_IN_CGS 1.6726231e24f
#define YEAR_IN_CGS 3.154e+7f

/* Hydrodynamical constants. */
#define const_hydro_gamma (5.0f / 3.0f)

#endif /* SWIFT_PHYSICAL_CONSTANTS_CGS_H */
