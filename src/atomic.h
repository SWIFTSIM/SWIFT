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
#define atomic_and(v, i) __sync_fetch_and_and(v, i)
#define atomic_inc(v) atomic_add(v, 1)
#define atomic_dec(v) atomic_sub(v, 1)
#define atomic_cas(v, o, n) __sync_val_compare_and_swap(v, o, n)
#define atomic_swap(v, n) __sync_lock_test_and_set(v, n)
#define atomic_read(v) __sync_val_compare_and_swap(v, 0, 0)
//#define atomic_write(v, i) __sync_bool_compare_and_swap(v, *v, i)

__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_write(volatile int *const address, const int y) {

  const int old_val = atomic_read(address);
  __sync_bool_compare_and_swap(address, old_val, y);
}

__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_write_u(volatile unsigned int *const address, const unsigned int y) {

  const unsigned int old_val = atomic_read(address);
  __sync_bool_compare_and_swap(address, old_val, y);
}

__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_write_f(volatile float *const address, const float y) {
  int *const address_int = (int *const)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type yy;
  yy.as_float = y;

  atomic_write(address_int, yy.as_int);
}

__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_write_c(volatile char *const address, const char y) {
  int *const address_int = (int *const)address;

  typedef union {
    char as_char;
    int as_int;
  } cast_type;

  cast_type yy;
  yy.as_char = y;

  atomic_write(address_int, yy.as_int);
}

__attribute__((always_inline, no_sanitize("thread"))) INLINE static float
atomic_read_f(volatile float *const address) {

  int *const address_int = (int *const)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type ret;
  ret.as_int = atomic_read(address_int);
  return ret.as_float;
}

__attribute__((always_inline, no_sanitize("thread"))) INLINE static unsigned int
atomic_read_u(volatile unsigned int *const address) {

  int *const address_int = (int *const)address;

  typedef union {
    unsigned int as_uint;
    int as_int;
  } cast_type;

  cast_type ret;
  ret.as_int = atomic_read(address_int);
  return ret.as_uint;
}

/**
 * @brief Atomic min operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point min that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_min_f(volatile float *const address, const float y) {

  int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = atomic_read_f(address);

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = min(old_val.as_float, y);
    old_val.as_int = atomic_cas(int_ptr, test_val.as_int, new_val.as_int);
  } while (test_val.as_int != old_val.as_int);
}

/**
 * @brief Atomic min operation on ints.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_min(
    volatile int *address, int y) {

  int *int_ptr = (int *)address;

  int test_val, old_val, new_val;
  old_val = *address;

  do {
    test_val = old_val;
    new_val = min(old_val, y);
    old_val = atomic_cas(int_ptr, test_val, new_val);
  } while (test_val != old_val);
}

/**
 * @brief Atomic min operation on doubles.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point min that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_min_d(
    volatile double *const address, const double y) {

  long long *const long_long_ptr = (long long *)address;

  typedef union {
    double as_double;
    long long as_long_long;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_double = *address;

  do {
    test_val.as_long_long = old_val.as_long_long;
    new_val.as_double = min(old_val.as_double, y);
    old_val.as_long_long =
        atomic_cas(long_long_ptr, test_val.as_long_long, new_val.as_long_long);
  } while (test_val.as_long_long != old_val.as_long_long);
}

/**
 * @brief Atomic max operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point max that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_max_f(volatile float *const address, const float y) {

  int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = atomic_read_f(address);

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = max(old_val.as_float, y);
    old_val.as_int = atomic_cas(int_ptr, test_val.as_int, new_val.as_int);
  } while (test_val.as_int != old_val.as_int);
}

/**
 * @brief Atomic max operation on doubles.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point max that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_max_d(
    volatile double *const address, const double y) {

  long long *const long_long_ptr = (long long *)address;

  typedef union {
    double as_double;
    long long as_long_long;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_double = *address;

  do {
    test_val.as_long_long = old_val.as_long_long;
    new_val.as_double = max(old_val.as_double, y);
    old_val.as_long_long =
        atomic_cas(long_long_ptr, test_val.as_long_long, new_val.as_long_long);
  } while (test_val.as_long_long != old_val.as_long_long);
}

/**
 * @brief Atomic add operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_add_f(volatile float *const address, const float y) {

  int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = atomic_read_f(address);

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = old_val.as_float + y;
    old_val.as_int = atomic_cas(int_ptr, test_val.as_int, new_val.as_int);
  } while (test_val.as_int != old_val.as_int);
}

/**
 * @brief Atomic add operation on doubles.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline, no_sanitize("thread"))) INLINE static void
atomic_add_d(volatile double *const address, const double y) {

  long long *const long_long_ptr = (long long *)address;

  typedef union {
    double as_double;
    long long as_long_long;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_double = *address;

  do {
    test_val.as_long_long = old_val.as_long_long;
    new_val.as_double = old_val.as_double + y;
    old_val.as_long_long =
        atomic_cas(long_long_ptr, test_val.as_long_long, new_val.as_long_long);
  } while (test_val.as_long_long != old_val.as_long_long);
}

#endif /* SWIFT_ATOMIC_H */
