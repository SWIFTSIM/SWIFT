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

#ifdef HAVE_PYTHON
/* Set the error message for python and return the error code
 * WARNING for python, you need to return NULL and not the error code
 */
#define error(err, s, ...)                  \
  ({                                        \
    error_no_return(err, s, ##__VA_ARGS__); \
    return err;                             \
  })

/* Same as error but does not return the error code.
 * Should be used for python functions */
#define error_no_return(err, s, ...)                                      \
  ({                                                                      \
    char error_msg[STRING_SIZE];                                          \
    sprintf(error_msg, "%s:%s():%i: " s ": %s\n", __FILE__, __FUNCTION__, \
            __LINE__, ##__VA_ARGS__, strerror(err));                      \
    PyErr_SetString(PyExc_RuntimeError, error_msg);                       \
  })

#else

#define error(err, s, ...)                  \
  ({                                        \
    error_no_return(err, s, ##__VA_ARGS__); \
    exit(1);                                \
  })

#define error_no_return(err, s, ...)                                      \
  ({                                                                      \
    char error_msg[STRING_SIZE];                                          \
    sprintf(error_msg, "%s:%s():%i: " s ": %s\n", __FILE__, __FUNCTION__, \
            __LINE__, ##__VA_ARGS__, strerror(err));                      \
  })

#endif

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
int tools_reverse_offset(const struct header *h, void *map, size_t *offset);
int tools_check_offset(const struct header *h, void *map, size_t *offset);

#endif  //__LOGGER_TOOLS_H__
