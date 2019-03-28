/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_MEMUSE_H
#define SWIFT_MEMUSE_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include <stdlib.h>

/* API. */
void memuse_use(long *size, long *resident, long *shared, long *text,
                long *data, long *library, long *dirty);
const char *memuse_process(int inmb);

#ifdef SWIFT_MEMUSE_REPORTS
void memuse_log_dump(const char *filename);
void memuse_log_dump_error(int rank);
void memuse_log_allocation(const char *label, void *ptr, int allocated,
                           size_t size);
#else

/* No-op when not reporting. */
#define memuse_log_allocation(label, ptr, allocated, size)
#endif

/**
 * @brief allocate aligned memory. The use and results are the same as the
 *        posix_memalign function. This function should be used for any
 *        significant allocations and consistently labelled.
 *
 * @param label a symbolic label for the memory, i.e. "parts".
 * @param memptr pointer to the allocated memory.
 * @param alignment alignment boundary.
 * @param size the quantity of bytes to allocate.
 * @result zero on success, otherwise an error code.
 */
__attribute__((always_inline)) inline int swift_memalign(const char *label,
                                                         void **memptr,
                                                         size_t alignment,
                                                         size_t size) {
  int result = posix_memalign(memptr, alignment, size);
#ifdef SWIFT_MEMUSE_REPORTS
  if (result == 0) {
    memuse_log_allocation(label, *memptr, 1, size);
  } else {
    /* Failed allocations are interesting as well. */
    memuse_log_allocation(label, NULL, -1, size);
  }
#endif
  return result;
}

/**
 * @brief allocate memory. The use and results are the same as the
 *        malloc function. This function should be used for any
 *        _significant_ allocations and consistently labelled.
 *        Do not use this function for small or high frequency
 *        allocations in production code.
 *
 * @param label a symbolic label for the memory, i.e. "parts".
 * @param size the quantity of bytes to allocate.
 * @result pointer to the allocated memory or NULL on failure.
 */
__attribute__((always_inline)) inline void *swift_malloc(const char *label,
                                                         size_t size) {
  void *memptr = malloc(size);
#ifdef SWIFT_MEMUSE_REPORTS
  if (memptr != NULL) {
    memuse_log_allocation(label, memptr, 1, size);
  } else {
    /* Failed allocations are interesting as well. */
    memuse_log_allocation(label, NULL, -1, size);
  }
#endif
  return memptr;
}

/**
 * @brief allocate zeroed memory. The use and results are the same as the
 *        calloc function. This function should be used for any
 *        _significant_ allocations and consistently labelled.
 *        Do not use this function for small or high frequency
 *        allocations in production code.
 *
 * @param label a symbolic label for the memory, i.e. "parts".
 * @param nmemb number of element to allocate.
 * @param size the size of each element in bytes.
 * @result pointer to the allocated memory or NULL on failure.
 */
__attribute__((always_inline)) inline void *swift_calloc(const char *label,
                                                         size_t nmemb,
                                                         size_t size) {
  void *memptr = calloc(nmemb, size);
#ifdef SWIFT_MEMUSE_REPORTS
  if (memptr != NULL) {
    memuse_log_allocation(label, memptr, 1, size * nmemb);
  } else {
    /* Failed allocations are interesting as well. */
    memuse_log_allocation(label, NULL, -1, size * nmemb);
  }
#endif
  return memptr;
}

/**
 * @brief free aligned memory. The use and results are the same as the
 *        free function. The label should match a prior call to swift_memalign
 *        or swift_malloc.
 *
 * @param label a symbolic label for the memory, i.e. "parts".
 * @param ptr pointer to the allocated memory.
 */
__attribute__((always_inline)) inline void swift_free(const char *label,
                                                      void *ptr) {
  free(ptr);
#ifdef SWIFT_MEMUSE_REPORTS
  memuse_log_allocation(label, ptr, 0, 0);
#endif
  return;
}

#endif /* SWIFT_MEMUSE_H */
