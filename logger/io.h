#ifndef __SWIFT_LOGGER_IO_H__
#define __SWIFT_LOGGER_IO_H__

#include "header.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief get the size of a file
 *
 * @param fd file id
 * @param size out: file size
 *
 * @return error code
 */
int io_get_file_size(int fd, size_t *size);

/**
 * @brief Open a file and map it
 *
 * @param filename file to read
 * @param fd out: file id
 * @param map out: file mapping
 *
 * @return error code
 */
int io_open_file(char *filename, int *fd, void **map);

/**
 * @brief Close a file and unmap it
 *
 * @param fd file id
 * @param map file mapping
 *
 * @return error code
 */
int io_close_file(int *fd, void **map);

/**
 * @brief read a single value in a file
 *
 * @param map file mapping
 * @param size size of the chunk to read
 * @param p pointer where to store the data
 * @param offset In: position to read, Out: shifted by size
 *
 * @return error code
 */
int io_read_data(void *map, const size_t size, void *p, size_t *offset);

/**
 * @brief write a single value in a file
 *
 * @param map file mapping
 * @param size size of the chunk to write
 * @param p pointer to the data
 * @param offset In: position to write, Out: shifted by size
 *
 * @return error code
 */
int io_write_data(void *map, const size_t size, const void *p, size_t *offset);

/**
 * @brief read a maks with its offset
 *
 * @param h #header file structure
 * @param map file mapping
 * @param offset In: position in the file, Out: shifted by the mask + offset
 * size
 * @param mask mask read
 * @param diff_offset offset difference to previous/next corresponding chunk
 *
 * @return error code
 */
int io_read_mask(const struct header *h, void *map, size_t *offset,
                 size_t *mask, size_t *diff_offset);

#endif  // __SWIFT_LOGGER_IO_H__
