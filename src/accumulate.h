/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ACCUMULATE_H
#define SWIFT_ACCUMULATE_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "atomic.h"
#include "minmax.h"

/**
 * @file accumulate.h
 * @brief Defines a series of functions used to update the fields of a particle
 * atomically. These functions should be used in any task that can run in
 * parallel to another one.
 */

/**
 * @brief Add x to the value stored at the location address (int version)
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to add to *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_add_i(
    atomic_int *const address, const int x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address += x;
#else
  atomic_add(address, x);
#endif
}

/**
 * @brief Add x to the value stored at the location address (long long version)
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to add to *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_add_ll(
    atomic_llong *const address, const long long x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address += x;
#else
  atomic_add(address, x);
#endif
}

/**
 * @brief Add x to the value stored at the location address (float version)
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to add to *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_add_f(
    atomic_float *const address, const float x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address += x;
#else
  atomic_add_f(address, x);
#endif
}

/**
 * @brief Add x to the value stored at the location address (double version)
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to add to *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_add_d(
    atomic_double *const address, const double x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address += x;
#else
  atomic_add_d(address, x);
#endif
}

/**
 * @brief Add 1 to the value stored at the location address (int version)
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 */
__attribute__((always_inline)) INLINE static void accumulate_inc_i(
    atomic_int *const address) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  (*address)++;
#else
  atomic_inc(address);
#endif
}

/**
 * @brief Add 1 to the value stored at the location address (long long version)
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 */
__attribute__((always_inline)) INLINE static void accumulate_inc_ll(
    atomic_llong *const address) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  (*address)++;
#else
  atomic_inc(address);
#endif
}

/**
 * @brief Compute the max of x and the value storedd at the location address
 * and store the value at the address (int8_t version).
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to max against *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_max_c(
    atomic_int8 *const address, const int8_t x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address = max(*address, x);
#else
  atomic_max_c(address, x);
#endif
}

/**
 * @brief Compute the max of x and the value storedd at the location address
 * and store the value at the address (int version).
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to max against *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_max_i(
    atomic_int *const address, const int x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address = max(*address, x);
#else
  atomic_max(address, x);
#endif
}

/**
 * @brief Compute the max of x and the value storedd at the location address
 * and store the value at the address (float version).
 *
 * When SWIFT_TASKS_WITHOUT_ATOMICS is *not* defined this function uses an
 * atomic operation.
 *
 * @param address The address to update.
 * @param x The value to max against *address.
 */
__attribute__((always_inline)) INLINE static void accumulate_max_f(
    atomic_float *const address, const float x) {

#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
  *address = max(*address, x);
#else
  atomic_max_f(address, x);
#endif
}

#endif /* SWIFT_ACCUMULATE_H */
