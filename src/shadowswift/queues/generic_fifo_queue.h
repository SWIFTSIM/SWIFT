//
// Created by yuyttenh on 24/03/22.
//

/**
 * @file generic_fifo_queue.h
 *
 * @brief Generates code for a FIFO queue with values of the type QUEUE_TYPE,
 * which has to be defined when importing this file.
 */

/* QUEUE_TYPE and optionally QUEUE_NAME must be defined before importing this
 * header! */

#define PASTE(x, y) x##_##y

#ifndef QUEUE_NAME
#define _NAME(f) PASTE(f, fifo_queue)
#define QUEUE_NAME _NAME(QUEUE_TYPE)
#endif

#define _INIT(f) PASTE(f, init)
#define QUEUE_INIT _INIT(QUEUE_NAME)

#define _DESTROY(f) PASTE(f, destroy)
#define QUEUE_DESTROY _DESTROY(QUEUE_NAME)

#define _RESET(f) PASTE(f, reset)
#define QUEUE_RESET _RESET(QUEUE_NAME)

#define _IS_EMPTY(f) PASTE(f, is_empty)
#define QUEUE_IS_EMPTY _IS_EMPTY(QUEUE_NAME)

#define _PUSH(f) PASTE(f, push)
#define QUEUE_PUSH _PUSH(QUEUE_NAME)

#define _POP(f) PASTE(f, pop)
#define QUEUE_POP _POP(QUEUE_NAME)

/**@brief
 * Generic definition of a FIFO queue
 */
struct QUEUE_NAME {
  QUEUE_TYPE *values;
  int size;
  int start;
  int end;
};

inline static void QUEUE_INIT(struct QUEUE_NAME *q, int size) {
  q->values = (QUEUE_TYPE *)malloc(size * sizeof(QUEUE_TYPE));
  q->size = size;
  q->start = 0;
  q->end = 0;
}

inline static void QUEUE_DESTROY(struct QUEUE_NAME *q) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(q->values != NULL);
#endif
  free(q->values);
  q->values = NULL;
  q->start = -1;
  q->end = -1;
}

inline static void QUEUE_RESET(struct QUEUE_NAME *q) {
  q->start = 0;
  q->end = 0;
}

inline static int QUEUE_IS_EMPTY(struct QUEUE_NAME *q) {
  return q->start == q->end;
}

inline static void QUEUE_PUSH(struct QUEUE_NAME *q, QUEUE_TYPE value) {
#ifdef QUEUE_SAFETY_CHECKS
  assert(q->size > 0);
#endif
  if (q->size == q->end) {
    q->size <<= 1;
    q->values = realloc(q->values, q->size * sizeof(QUEUE_TYPE));
  }
  q->values[q->end] = value;
  q->end++;
}

inline static QUEUE_TYPE QUEUE_POP(struct QUEUE_NAME *q) {
#ifdef QUEUE_SAFETY_CHECKS
  if (QUEUE_IS_EMPTY(q)) {
    fprintf(stderr, "Trying to pop from empty queue!");
    abort();
  }
#endif
  QUEUE_TYPE v = q->values[q->start];
  q->start++;
  return v;
}

/* Undefine QUEUE_NAME and QUEUE_TYPE */
#undef QUEUE_NAME
#undef QUEUE_TYPE
#undef _NAME
