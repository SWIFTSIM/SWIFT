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
#ifndef SWIFT_LOCK_H
#define SWIFT_LOCK_H

/* Some standard headers. */
#include <pthread.h>

/* Includes. */
#include "atomic.h"

#ifdef PTHREAD_SPINLOCK
#include <pthread.h>
#define swift_lock_type pthread_spinlock_t
#define lock_init(l) (pthread_spin_init(l, PTHREAD_PROCESS_PRIVATE) != 0)
#define lock_destroy(l) (pthread_spin_destroy(l) != 0)
#define lock_lock(l) (pthread_spin_lock(l) != 0)
#define lock_trylock(l) (pthread_spin_lock(l) != 0)
#define lock_unlock(l) (pthread_spin_unlock(l) != 0)
#define lock_unlock_blind(l) pthread_spin_unlock(l)
#define lock_static_initializer ((pthread_spinlock_t)0)

#elif defined(PTHREAD_LOCK)
#include <pthread.h>
#define swift_lock_type pthread_mutex_t
#define lock_init(l) (pthread_mutex_init(l, NULL) != 0)
#define lock_destroy(l) (pthread_mutex_destroy(l) != 0)
#define lock_lock(l) (pthread_mutex_lock(l) != 0)
#define lock_trylock(l) (pthread_mutex_trylock(l) != 0)
#define lock_unlock(l) (pthread_mutex_unlock(l) != 0)
#define lock_unlock_blind(l) pthread_mutex_unlock(l)
#define lock_static_initializer PTHREAD_MUTEX_INITIALIZER

#else
#define swift_lock_type volatile int
#define lock_init(l) (*(l) = 0)
#define lock_destroy(l) 0
INLINE static int lock_lock(volatile int *l) {
  while (atomic_cas(l, 0, 1) != 0)
    ;
  return 0;
}
#define lock_trylock(l) ((*(l)) ? 1 : atomic_cas(l, 0, 1))
#define lock_unlock(l) (atomic_cas(l, 1, 0) != 1)
#define lock_unlock_blind(l) atomic_cas(l, 1, 0)
#define lock_static_initializer 0
#endif

#endif /* SWIFT_LOCK_H */
