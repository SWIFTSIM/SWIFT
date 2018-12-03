#ifndef __LOGGER_TOOLS_H__
#define __LOGGER_TOOLS_H__

#include "../config.h"

#ifdef HAVE_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRING_SIZE 200

struct header;

#define error(s, ...)							\
  ({									\
    fflush(stdout);							\
    fprintf(stderr, "%s:%s():%i: " s "\n",				\
            __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__);		\
    abort();								\
  })

#define message(s, ...)                                             \
  ({                                                                \
    printf("%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__, \
           ##__VA_ARGS__);                                          \
  })

int tools_get_next_chunk(const struct header *h, void *map, size_t *offset,
                         int fd);
int _tools_get_next_chunk_backward(const struct header *h, void *map,
                                   size_t *offset, int fd);
int _tools_get_next_chunk_forward(const struct header *h, void *map,
                                  size_t *offset);
void tools_reverse_offset(const struct header *h, void *map, size_t *offset);
void tools_check_offset(const struct header *h, void *map, size_t *offset);

#endif  //__LOGGER_TOOLS_H__
