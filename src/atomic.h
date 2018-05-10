/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_ATOMIC_H
#define SWIFT_ATOMIC_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "inline.h"
#include "minmax.h"

#define atomic_add(v, i) __sync_fetch_and_add(v, i)
#define atomic_sub(v, i) __sync_fetch_and_sub(v, i)
#define atomic_or(v, i) __sync_fetch_and_or(v, i)
#define atomic_inc(v) atomic_add(v, 1)
#define atomic_dec(v) atomic_sub(v, 1)
#define atomic_cas(v, o, n) __sync_val_compare_and_swap(v, o, n)
#define atomic_swap(v, n) __sync_lock_test_and_set(v, n)

/**
 * @brief Atomic min operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_min_f(
    volatile float* address, float y) {

  int* int_ptr = (int*)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type assumed, old, new;
  old.as_float = *address;

  do {
    assumed.as_int = old.as_int;
    new.as_float = min(old.as_float, y);
    old.as_int = atomic_cas(int_ptr, assumed.as_int, new.as_int);
  } while (assumed.as_int != old.as_int);
}

/**
 * @brief Atomic max operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_max_f(
    volatile float* address, float y) {

  int* int_ptr = (int*)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type assumed, old, new;
  old.as_float = *address;

  do {
    assumed.as_int = old.as_int;
    new.as_float = max(old.as_float, y);
    old.as_int = atomic_cas(int_ptr, assumed.as_int, new.as_int);
  } while (assumed.as_int != old.as_int);
}

/**
 * @brief Atomic add operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_add_f(
    volatile float* address, float y) {

  int* int_ptr = (int*)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type assumed, old, new;
  old.as_float = *address;

  do {
    assumed.as_int = old.as_int;
    new.as_float = old.as_float + y;
    old.as_int = atomic_cas(int_ptr, assumed.as_int, new.as_int);
  } while (assumed.as_int != old.as_int);
}

#endif /* SWIFT_ATOMIC_H */
