//
// Created by yuyttenh on 24/03/22.
//

/**
 * @file queues.h
 *
 * @brief Generates code for a int LIFO queue and an int3 FIFO queue
 */

#ifndef SWIFTSIM_SHADOWSWIFT_QUEUES_H
#define SWIFTSIM_SHADOWSWIFT_QUEUES_H

#include "error.h"

/**
 * @brief A simple int2 tuple.
 */
typedef struct int2 {
  int _0;
  int _1;
} int2;

/**
 * @brief A simple int3 tuple.
 */
typedef struct int3 {
  int _0;
  int _1;
  int _2;
} int3;

#define QUEUE_SAFETY_CHECKS

#define QUEUE_TYPE int
#include "queues/generic_lifo_queue.h"
#undef QUEUE_TYPE

#define QUEUE_TYPE int2
#include "queues/generic_lifo_queue.h"
#undef QUEUE_TYPE

#define QUEUE_TYPE int3
#include "queues/generic_fifo_queue.h"
#undef QUEUE_TYPE

#define QUEUE_TYPE int
#include "queues/generic_fifo_queue.h"
#undef QUEUE_TYPE

#endif  // SWIFTSIM_SHADOWSWIFT_QUEUES_H
