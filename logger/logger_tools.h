#ifndef __TOOLS_H__
#define __TOOLS_H__

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>

#define STRING_SIZE 200

struct header;

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

#define message(s, ...)                                             \
  ({                                                                \
    printf("%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__, \
           ##__VA_ARGS__);                                          \
  })

/**
 * @brief get the offset of the next corresponding chunk
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of the next chunk
 * @param fd file id
 *
 * @return error code, -1 if no next chunk
 */
int tools_get_next_chunk(const struct header *h, void *map, size_t *offset,
                         int fd);

/**
 * @brief internal function of #tools_get_next_chunk. Should not be used (very
 * slow)
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of the next chunk
 * @param fd file id
 *
 * @return error code, -1 if no next chunk
 */
int _tools_get_next_chunk_backward(const struct header *h, void *map,
                                   size_t *offset, int fd);

/**
 * @brief internal function of #tools_get_next_chunk. Should not be used outside
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of the next chunk
 *
 * @return error code, -1 if no next chunk
 */
int _tools_get_next_chunk_forward(const struct header *h, void *map,
                                  size_t *offset);

/**
 * @brief switch side offset
 *
 * From current chunk, switch side of the offset of the previous one.
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of next chunk
 *
 * @return error code
 */
int tools_reverse_offset(const struct header *h, void *map, size_t *offset);

/**
 * @brief debugging function checking the offset of a chunk
 *
 * Compare the mask with the one pointed by the header.
 * if the chunk is a particle, check the id too.
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: chunk offset, Out: offset after the chunk
 *
 * @return error code
 */
int tools_check_offset(const struct header *h, void *map, size_t *offset);

#endif  //__TOOLS_H__
