#ifndef __LOGGER_HEADER_H__
#define __LOGGER_HEADER_H__

#include "tools.h"

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>

#define LOGGER_VERSION_SIZE 20
#define LOGGER_NBER_SIZE 4
#define LOGGER_OFFSET_SIZE 7
#define LOGGER_MASK_SIZE 1

struct header {
  /* Logger version */
  char version[STRING_SIZE];

  /* offset of the first header */
  size_t offset_first;

  /* Number of bytes for names */
  size_t name;
  
  /* number of masks */
  size_t nber_mask;

  /* list of masks */
  size_t *masks;

  /* list of mask name */
  char **masks_name;

  /* size of data in a mask */
  size_t *masks_size;

  /* offset direction */
  int forward_offset;
};

/**
 * @brief print a header struct
 */
void header_print(const struct header *h);

/**
 * @brief free allocated memory
 */
void header_free(struct header *h);

/**
 * @brief check if field is present in header
 *
 * @param field name of the requested field
 * @param ind (return value) indice of the requested field
 */
int header_is_present_and_get_index(const struct header *h, const char *field,
                                    size_t *ind);

/**
 * @brief check if field is present in header
 *
 * @param field name of the requested field
 */
int header_is_present(const struct header *h, const char *field);

/**
 * @brief read the logger header
 *
 * @param h out: header
 * @param map file mapping
 */
int header_read(struct header *h, void *map);

/**
 * @brief count number of bits in a given mask
 *
 * @param h #header file structure
 * @param mask mask to compute
 *
 * @return number of bits in mask
 */
size_t header_get_mask_size(const struct header *h, const size_t mask);

/**
 * @brief Inverse the offset direction
 *
 * @param h #header file structure
 * @param map file mapping
 *
 * @return error code
 */
int header_change_offset_direction(struct header *h, void *map);

#endif  // __LOGGER_HEADER_H__
