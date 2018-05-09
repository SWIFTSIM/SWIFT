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
 * @param Atomic min operation on floats.
 */
__attribute__((always_inline)) INLINE void atomic_min_f(float const* x,
                                                        float y) {
  int done = 0;
  while (!done) {
    const float val = *x;
    done = __sync_bool_compare_and_swap((int*)x, val, min(val, y));
  }
}

/**
 * @param Atomic max operation on floats.
 */
__attribute__((always_inline)) INLINE void atomic_max_f(float const* x,
                                                        float y) {
  int done = 0;
  while (!done) {
    const float val = *x;
    done = __sync_bool_compare_and_swap((int*)x, val, max(val, y));
  }
}

/**
 * @param Atomic add operation on floats.
 */
__attribute__((always_inline)) INLINE void atomic_add_f(float const* x,
                                                        float y) {
  int done = 0;
  while (!done) {
    const float val = *x;
    done = __sync_bool_compare_and_swap((int*)x, val, val + y);
  }
}

#endif /* SWIFT_ATOMIC_H */
