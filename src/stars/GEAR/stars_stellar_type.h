/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_STARS_GEAR_STELLAR_TYPE_H
#define SWIFT_STARS_GEAR_STELLAR_TYPE_H

/**
 * @file src/stars/GEAR/stars_stellar_type.h
 * @brief header file concerning the stellar type of the star particle.
 **/

/**
 * @brief The stellar type.
 *
 * Star particles can represent a single star ("single_star"), a stellar
 * population from a continuous IMF or a stellar population from a whole IMF.
 */
enum stellar_type {
  single_star = 0,                /* particle representing a single star */
  star_population_continuous_IMF, /* particle representing a population of the
                                     continuous part of the IMF */
  star_population, /* particle representing a population with the whole IMF */
  stellar_type_count
};

#endif /* SWIFT_STARS_GEAR_STELLAR_TYPE_H */
