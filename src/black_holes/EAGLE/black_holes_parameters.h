/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_BLACK_HOLES_PARAMETERS_H
#define SWIFT_EAGLE_BLACK_HOLES_PARAMETERS_H

/* Configuration file */
#include "config.h"

/**
 * @file EAGLE/black_holes_parameters.h
 * @brief Parameters of the EAGLE black holes
 *        model that need to be defined at compile time.
 */

/*! Maximal distance for merging particles in units of the (spline not Plummer)
 *  softening length. */
#define const_max_merging_distance_ratio 3.f

/*! Maximal distance for repositioning particles in units of the (spline not
 * Plummer) softening length. */
#define const_max_repositioning_distance_ratio 3.f

#endif /* SWIFT_EAGLE_BLACK_HOLES_PARAMETERS_H */
