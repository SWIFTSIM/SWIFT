#include "logger_header.h"

#include "logger_io.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief print a header struct
 *
 * @param h The #header
 */
void header_print(const struct header *h) {
#ifdef SWIFT_DEBUG_CHECKS
  printf("Debug checks enabled\n");
#endif
  printf("Version:          %s\n", h->version);
  printf("First Offset:     %lu\n", h->offset_first);
  char direction[20];
  if (h->forward_offset)
    strcpy(direction, "Forward");
  else
    strcpy(direction, "Backward");
  printf("Offset direction: %s\n", direction);
  printf("Number masks:     %lu\n", h->nber_mask);

  for (size_t i = 0; i < h->nber_mask; i++) {
    printf("\tMask:  %s\n", h->masks_name[i]);
    printf("\tValue: %lu\n", h->masks[i]);
    printf("\tSize:  %lu\n", h->masks_size[i]);
    printf("\n");
  }
};

/**
 * @brief free allocated memory
 *
 * @param h The #header
 */
void header_free(struct header *h) {
  for (size_t i = 0; i < h->nber_mask; i++) {
    free(h->masks_name[i]);
  }
  free(h->masks_name);
  free(h->masks);
  free(h->masks_size);
};

/**
 * @brief check if field is present in header
 *
 * @param h The #header
 * @param field name of the requested field
 * @param ind (return value) indice of the requested field
 */
int header_field_is_present(const struct header *h, const char *field,
			    size_t *ind) {
  for (size_t i = 0; i < h->nber_mask; i++) {
    if (strcmp(h->masks_name[i], field) == 0) {
      if (ind != NULL) {
        *ind = i;
      }
      return 1;
    }
  }

  return 0;
};

/**
 * @brief Inverse the offset direction
 *
 * @param h #header file structure
 * @param map file mapping
 *
 */
void header_change_offset_direction(struct header *h, void *map) {
  h->forward_offset = !h->forward_offset;
  size_t offset = LOGGER_VERSION_SIZE;

  io_write_data(map, LOGGER_NBER_SIZE, &h->forward_offset, &offset);
}

/**
 * @brief read the logger header
 *
 * @param h out: header
 * @param map file mapping
 */
void header_read(struct header *h, void *map) {
  size_t offset = 0;

  /* read version */
  io_read_data(map, LOGGER_VERSION_SIZE, &h->version, &offset);

  /* read offset direction */
  h->forward_offset = 0;
  io_read_data(map, LOGGER_NBER_SIZE, &h->forward_offset, &offset);

  if (h->forward_offset != 0 && h->forward_offset != 1)
    error("Non boolean value for the offset direction (%i)",
          h->forward_offset);

  /* read offset to first data */
  h->offset_first = 0;
  io_read_data(map, LOGGER_OFFSET_SIZE, &h->offset_first, &offset);

  /* read name size */
  h->name = 0;
  io_read_data(map, LOGGER_NBER_SIZE, &h->name, &offset);

  /* check if value defined in this file is large enough */
  if (STRING_SIZE < h->name) {
    error("Name too large in dump file");
  }

  /* read number of masks */
  h->nber_mask = 0;
  io_read_data(map, LOGGER_NBER_SIZE, &h->nber_mask, &offset);

  /* allocate memory */
  h->masks = malloc(sizeof(size_t) * h->nber_mask);
  h->masks_name = malloc(sizeof(void *) * h->nber_mask);
  h->masks_size = malloc(sizeof(size_t) * h->nber_mask);

  /* loop over all masks */
  for (size_t i = 0; i < h->nber_mask; i++) {
    /* read mask name */
    h->masks_name[i] = malloc(h->name);
    io_read_data(map, h->name, h->masks_name[i], &offset);

    /* get mask value */
    h->masks[i] = 1 << i;

    /* read mask data size */
    h->masks_size[i] = 0;
    io_read_data(map, LOGGER_NBER_SIZE, &h->masks_size[i], &offset);
  }

  if (offset != h->offset_first) {
#ifdef SWIFT_DEBUG_CHECKS
    header_print(h);
#endif
    error("Wrong header size (in header %li, current %li)",
          h->offset_first, offset);
  }

};

/**
 * @brief count number of bits in a given mask
 *
 * @param h #header file structure
 * @param mask mask to compute
 *
 * @return number of bits in mask
 */
size_t header_get_mask_size(const struct header *h, const size_t mask) {
  size_t count = 0;
  for (size_t i = 0; i < h->nber_mask; i++) {
    if (mask & h->masks[i]) {
      count += h->masks_size[i];
    }
  }
  return count;
}
