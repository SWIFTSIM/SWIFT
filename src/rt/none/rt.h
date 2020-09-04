/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_NONE_H
#define SWIFT_RT_NONE_H

/**
 * @file src/rt/none/rt.h
 * @brief Main header file for no radiative transfer scheme.
 */

/**
 * @brief First initialisation of the RT extra hydro partilce data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_xpart(struct xpart* restrict xp) {}

/**
 * @brief Initialisation of the RT extra hydro partilce data.
 */
__attribute__((always_inline)) INLINE static void rt_init_xpart(struct xpart* restrict xp) {
}


/**
 * @brief First initialisation of the RT extra star partilce data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart( struct spart* restrict sp) {
}


/**
 * @brief First initialisation of the RT extra star partilce data.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart( struct spart* restrict sp) {
}

#endif /* SWIFT_RT_NONE_H */
