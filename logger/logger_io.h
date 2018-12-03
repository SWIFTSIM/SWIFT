#ifndef __SWIFT_LOGGER_IO_H__
#define __SWIFT_LOGGER_IO_H__

#include "logger_header.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>

void io_get_file_size(int fd, size_t *size);
void io_open_file(char *filename, int *fd, void **map);
void io_close_file(int *fd, void **map);
void io_read_data(void *map, const size_t size, void *p, size_t *offset);
void io_write_data(void *map, const size_t size, const void *p, size_t *offset);
void io_read_mask(const struct header *h, void *map, size_t *offset,
                 size_t *mask, size_t *diff_offset);

#endif  // __SWIFT_LOGGER_IO_H__
