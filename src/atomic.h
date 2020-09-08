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

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "inline.h"
#include "minmax.h"

#include <stdatomic.h>

#define atomic_add(v, i) atomic_fetch_add(v, i)
#define atomic_and(v, i) atomic_fetch_and(v, i)
#define atomic_sub(v, i) atomic_fetch_sub(v, i)
#define atomic_or(v, i) atomic_fetch_or(v, i)
#define atomic_inc(v) atomic_fetch_add(v, 1)
#define atomic_dec(v) atomic_fetch_sub(v, 1)
#define atomic_cas_fast(o, e, d) atomic_compare_exchange_weak(o, e, d)
#define atomic_write(o, v) atomic_store(o, v)
#define atomic_write_u(o, v) atomic_store(o, v)
#define atomic_write_c(o, v) atomic_store(o, v)
#define atomic_write_f(o, v) atomic_write(o, v)
#define atomic_load_f(v) atomic_load(v)
#define atomic_load_u(v) atomic_load(v)
typedef _Atomic(float) atomic_float;
typedef _Atomic(double) atomic_double;
typedef _Atomic(uint32_t) atomic_uint32;
typedef _Atomic(uint16_t) atomic_uint16;
typedef _Atomic(int8_t) atomic_int8;
//typedef _Atomic(int16_t) atomic_short;
#define atomic_cas(obj,expected,desired) _Generic((desired), \
  int: bool_atomic_compare_and_swap_i, \
  long long: bool_atomic_compare_and_swap_ll, \
  float: bool_atomic_compare_and_swap_f, \
  char: bool_atomic_compare_and_swap_c, \
  default: bool_atomic_compare_and_swap_ll \
  )(obj,expected,desired)
#define atomic_swap(v, n) atomic_exchange(v,n)

__attribute__((always_inline)) INLINE _Bool bool_atomic_compare_and_swap_c( atomic_char *obj, char *expected, char desired){
  char preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong(obj, expected, desired);
  *expected = preexpected;
  return retval;
}

__attribute__((always_inline)) INLINE _Bool bool_atomic_compare_and_swap_i( atomic_int *obj, int *expected, int desired){
  int preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong(obj, expected, desired);
  *expected = preexpected;
  return retval;
}

__attribute__((always_inline)) INLINE _Bool bool_atomic_compare_and_swap_ll( atomic_llong *obj, long long *expected, long long desired){
  long long preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong(obj, expected, desired);
  *expected = preexpected;
  return retval;
}

__attribute__((always_inline)) INLINE _Bool bool_atomic_compare_and_swap_f( atomic_float *obj, float *expected, float desired){
  float preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong( obj, expected, desired);
  *expected = preexpected;
  return retval;
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
__attribute__((always_inline)) INLINE static void atomic_min_f( atomic_float *obj, float const y){
  float test_val, new_val;

  do{
    test_val = atomic_load(obj);
    new_val = fminf(test_val, y);
  }while(!atomic_cas_fast( obj, &test_val, new_val));
}



__attribute__((always_inline)) INLINE static void atomic_max(atomic_int *obj, const int y) {
  int test_val, new_val;

  do{
    test_val = atomic_load(obj);
    new_val = max(test_val, y);
  }while(!atomic_cas_fast(obj, &test_val, new_val));
}

__attribute__((always_inline)) INLINE static void atomic_max_c(atomic_int8 *obj, const int8_t y) {
  int8_t test_val, new_val;

  do{
    test_val = atomic_load(obj);
    new_val = max(test_val, y);
  }while(!atomic_cas_fast(obj, &test_val, new_val));
}

__attribute__((always_inline)) INLINE static void atomic_max_f( atomic_float *obj, float const y){
  float test_val, new_val;

  do{
    test_val = atomic_load(obj);
    new_val = fmaxf(test_val, y);
  }while(!atomic_cas_fast(obj, &test_val, new_val));

}

__attribute__((always_inline)) INLINE static void atomic_max_d( atomic_double *obj, double const y){
  double test_val, new_val;

  do{
    test_val = atomic_load(obj);
    new_val = fmax(test_val, y);
  }while(!atomic_cas_fast(obj, &test_val, new_val));
}

__attribute__((always_inline)) INLINE static void atomic_add_f(
    atomic_float *const address, const float y) {
  atomic_int *const int_ptr = (atomic_int *) address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, new_val;

  do{
    test_val.as_int = atomic_load(address);
    new_val.as_float = test_val.as_float + y;
  } while(!atomic_cas(int_ptr, &test_val.as_int, new_val.as_int));
}

/**
 * @brief Atomic add operation on doubles.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the double add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_add_d(
    atomic_double *const address, const double y) {
  atomic_int *const int_ptr = (atomic_int *) address;

  typedef union {
    double as_double;
    int as_int;
  } cast_type;

  cast_type test_val, new_val;

  do{
    test_val.as_int = atomic_load(address);
    new_val.as_double = test_val.as_double + y;
  } while(!atomic_cas(int_ptr, &test_val.as_int, new_val.as_int));
}

#endif /* SWIFT_ATOMIC_H */
