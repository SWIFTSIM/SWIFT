/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_BARRIER_H
#define SWIFT_BARRIER_H

/**
 * @file barrier.h
 * @brief Define the thread barriers if the POSIX implementation on this system
 * does not.
 *
 * The pthread barriers are only an option of the POSIX norm and they are not
 * necessarily implemented. One example is OSX where all the rest of POSIX
 * exists but not the barriers.
 * We implement them here in a simple way to allow for SWIFT to run on such
 * systems but this may lead to poorer performance.
 *
 * Note that we only define the three functions we need. This is a textbook
 * implementation of a barrier that uses the common POSIX features (mutex,
 * conditions and broadcasts).
 *
 * If the pthread barriers exist (Linux systems), we default to them.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers */
#include <pthread.h>

/* Does this POSIX implementation provide barriers? */
#ifdef HAVE_PTHREAD_BARRIERS

#define swift_barrier_t pthread_barrier_t
#define swift_barrier_wait pthread_barrier_wait
#define swift_barrier_init pthread_barrier_init
#define swift_barrier_destroy pthread_barrier_destroy

#else

/* Local headers */
#include "error.h"
#include "inline.h"

/**
 * @brief An ersatz of POSIX barriers to be used on systems that don't provide
 * the good ones.
 */
typedef struct {

  /*! Barrier mutex */
  pthread_mutex_t mutex;

  /*! Condition to open the barrier */
  pthread_cond_t condition;

  /*! Total number of threads */
  int limit;

  /*! Number of threads that reached the barrier */
  int count;

} swift_barrier_t;

/**
 * @brief Initialise a barrier object.
 *
 * @param barrier The #swift_barrier_t to initialise
 * @param unused Unused parameter (NULL) as we don't support barrier attributes.
 * @param count The number of threads that will wait at the barrier.
 */
static INLINE int swift_barrier_init(swift_barrier_t *barrier, void *unused,
                                     unsigned int count) {
  /* Initialise the mutex */
  if (pthread_mutex_init(&barrier->mutex, 0) != 0)
    error("Error initializing the barrier mutex");

  /* Initialise the condition */
  if (pthread_cond_init(&barrier->condition, 0) != 0)
    error("Error initializing the barrier condition");

  barrier->limit = count;
  barrier->count = 0;

  /* All is good */
  return 0;
}

/**
 * @brief Make a set of threads wait at the barrier
 *
 * Note that once all threads have reached the barrier, we also
 * reset the barrier to state where it is ready to be re-used
 * without calling swift_barrier_init.
 *
 * @param barrier The (initialised) #swift_barrier_t to wait at.
 */
static INLINE int swift_barrier_wait(swift_barrier_t *barrier) {

  /* Start by locking the barrier */
  pthread_mutex_lock(&barrier->mutex);

  /* One more thread has gone home*/
  barrier->count++;

  /* Are threads still running? */
  if (barrier->count < barrier->limit) {

    /* We need to make the thread wait until everyone is back */
    pthread_cond_wait(&barrier->condition, &(barrier->mutex));

    /* Release the mutex */
    pthread_mutex_unlock(&barrier->mutex);

    /* Say that this was not the last thread */
    return 0;

  } else { /* Everybody is home */

    /* Open the barrier (i.e. release the threads blocked in the while loop) */
    pthread_cond_broadcast(&barrier->condition);

    /* Re-initialize the barrier */
    barrier->count = 0;

    /* Release the mutex */
    pthread_mutex_unlock(&barrier->mutex);

    /* Say that we are all done */
    return 1;
  }
}

/**
 * @brief Destroy a barrier object
 *
 * Note that if destroy is called before a barrier is open, we return
 * an error message and do not attempt to wait for the barrier to open
 * before destroying it.
 *
 * @param barrier The #swift_barrier_t object to destroy.
 */
static INLINE int swift_barrier_destroy(swift_barrier_t *barrier) {

  /* Destroy the pthread things */
  pthread_cond_destroy(&barrier->condition);
  pthread_mutex_destroy(&barrier->mutex);

  /* All is good */
  return 0;
}

#endif /* HAVE_PTHREAD_BARRIERS */

#endif /* SWIFT_BARRIER_H */
