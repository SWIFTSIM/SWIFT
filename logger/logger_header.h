#ifndef __LOGGER_HEADER_H__
#define __LOGGER_HEADER_H__

#include "logger_tools.h"

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

void header_print(const struct header *h);
void header_free(struct header *h);
int header_field_is_present(const struct header *h, const char *field,
			    size_t *ind);
void header_read(struct header *h, void *map);
size_t header_get_mask_size(const struct header *h, const size_t mask);
void header_change_offset_direction(struct header *h, void *map);

#endif  // __LOGGER_HEADER_H__
