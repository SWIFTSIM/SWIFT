#include "logger_tools.h"
#include "header.h"
#include "io.h"

#include "particle.h"

#include <stdio.h>

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
                         int fd) {
  if (h->forward_offset)
    return _tools_get_next_chunk_forward(h, map, offset);
  else
    return _tools_get_next_chunk_backward(h, map, offset, fd);
}

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
                                  size_t *offset) {
  size_t diff_offset = 0;

  /* read offset */
  int test = io_read_mask(h, map, offset, NULL, &diff_offset);
  if (test != 0) return test;

  if (diff_offset == 0) return -1;

  /* set absolute offset */
  *offset += diff_offset - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;
  return 0;
}

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
                                   size_t *offset, int fd) {
#ifndef SWIFT_DEBUG_CHECKS
  error(ENOTSUP, "Should not be used, method too slow");
#endif
  size_t current_offset = *offset;
  size_t chunk_header = LOGGER_MASK_SIZE + LOGGER_OFFSET_SIZE;

  size_t file_size = 0;
  int test = io_get_file_size(fd, &file_size);

  if (test != 0) return test;

  while (current_offset < file_size) {
    size_t mask = 0;
    size_t prev_offset;
    test = io_read_mask(h, map, &current_offset, &mask, &prev_offset);
    if (test != 0) return test;

    prev_offset = current_offset - prev_offset - chunk_header;
    if (*offset == prev_offset) {
      *offset = current_offset - chunk_header;
      return 0;
    }

    current_offset += header_get_mask_size(h, mask);
  }

  return -1;
}

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
int tools_reverse_offset(const struct header *h, void *map, size_t *offset) {
  size_t mask = 0;
  size_t prev_offset = 0;
  const size_t cur_offset = *offset;

  /* read mask + offset */
  int error_code = io_read_mask(h, map, offset, &mask, &prev_offset);
  if (error_code != 0) return error_code;

  /* write offset of zero (in case it is the last chunk) */
  const size_t zero = 0;
  *offset -= LOGGER_OFFSET_SIZE;
  error_code = io_write_data(map, LOGGER_OFFSET_SIZE, &zero, offset);
  if (error_code != 0) error(errno, "Unable to write data while reversing");

  /* set offset after current chunk */
  *offset += header_get_mask_size(h, mask);

  /* first chunks do not have a previous partner */
  if (prev_offset == cur_offset) return 0;

  if (prev_offset > cur_offset)
    error(EIO, "Unexpected offset, header %lu, current %lu", prev_offset,
          cur_offset);

  /* modify previous offset */
  size_t abs_prev_offset = cur_offset - prev_offset + LOGGER_MASK_SIZE;
  error_code =
      io_write_data(map, LOGGER_OFFSET_SIZE, &prev_offset, &abs_prev_offset);
  if (error_code != 0) error(errno, "Unable to write offset");

#ifdef SWIFT_DEBUG_CHECKS
  size_t prev_mask = 0;
  abs_prev_offset -= LOGGER_MASK_SIZE + LOGGER_OFFSET_SIZE;
  error_code = io_read_mask(h, map, &abs_prev_offset, &prev_mask, NULL);

  if (prev_mask != mask)
    error(EIO, "Unexpected mask: %lu (at %lu), got %lu (at %lu)", mask, *offset,
          prev_mask, prev_offset);

#endif  // SWIFT_DEBUG_CHECKS

  return 0;
}

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
int tools_check_offset(const struct header *h, void *map, size_t *offset) {
#ifndef SWIFT_DEBUG_CHECKS
  error(ENOTSUP, "Should not check in non debug mode");
#endif

  size_t tmp = *offset;

  size_t mask;
  size_t pointed_offset;

  /* read mask + offset */
  int error_code = io_read_mask(h, map, offset, &mask, &pointed_offset);
  if (error_code != 0) return error_code;

  /* get absolute offset */
  if (h->forward_offset)
    pointed_offset += tmp;
  else {
    if (tmp < pointed_offset)
      error(EIO, "Offset too large (%lu) at %lu with mask %lu", pointed_offset,
            tmp, mask);
    pointed_offset = tmp - pointed_offset;
  }

  /* set offset after current chunk */
  *offset += header_get_mask_size(h, mask);

  if (pointed_offset == tmp || pointed_offset == 0) return 0;

  /* read mask of the pointed chunk */
  size_t pointed_mask = 0;
  error_code = io_read_mask(h, map, &pointed_offset, &pointed_mask, NULL);

  /* check masks */
  if (pointed_mask != mask)
    error(EIO, "Error in the offset (mask %lu != %lu) at %lu and %lu", mask,
          pointed_mask,
          *offset - header_get_mask_size(h, mask) - LOGGER_MASK_SIZE -
              LOGGER_OFFSET_SIZE,
          pointed_offset - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE);

  if (pointed_mask == 128) return 0;

  struct particle part;
  error_code = particle_read(&part, h, map, &tmp, 0, reader_const, NULL);

  size_t id = part.id;
  tmp = pointed_offset - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;
  error_code = particle_read(&part, h, map, &tmp, 0, reader_const, NULL);

  if (id != part.id)
    error(EIO, "Offset wrong, id incorrect (%lu != %lu) at %lu", id, part.id,
          tmp);

  return 0;
}
