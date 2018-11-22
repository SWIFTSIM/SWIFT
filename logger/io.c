#include "io.h"
#include "header.h"
#include "logger_tools.h"

/* file size */
#include <sys/stat.h>
/* mapping */
#include <sys/mman.h>
/* open */
#include <fcntl.h>

int io_get_file_size(int fd, size_t *size) {
  struct stat s;
  int status = fstat(fd, &s);
  if (status != 0) error(errno, "Unable to get file size");
  *size = s.st_size;

  return 0;
}

int io_open_file(char *filename, int *fd, void **map) {
  /* open file */
  *fd = open(filename, O_RDWR);
  if (*fd == -1) error(errno, "Unable to open file %s", filename);

  /* get file size */
  size_t size = 0;
  int status = io_get_file_size(*fd, &size);
  if (status != 0) return status;

  /* map memory */
  *map = mmap(NULL, size, PROT_WRITE | PROT_READ, MAP_SHARED, *fd, 0);
  if (map == MAP_FAILED)
    error(errno, "Failed to allocate map of size %zi bytes.", size);

  return 0;
}

int io_close_file(int *fd, void **map) {
  /* get file size */
  size_t size = 0;
  int status = io_get_file_size(*fd, &size);
  if (status != 0) return status;

  /* unmap */
  if (munmap(*map, size) != 0) error(errno, "Unable to unmap the file");

  close(*fd);

  return 0;
}

int io_read_mask(const struct header *h, void *map, size_t *offset,
                 size_t *mask, size_t *diff_offset) {
  /* read mask */
  if (mask) {
    *mask = 0;
    memcpy(mask, map + *offset, LOGGER_MASK_SIZE);
  }
  *offset += LOGGER_MASK_SIZE;

  /* read offset */
  if (diff_offset) {
    *diff_offset = 0;
    memcpy(diff_offset, map + *offset, LOGGER_OFFSET_SIZE);
  }
  *offset += LOGGER_OFFSET_SIZE;

  return 0;
}

int io_read_data(void *map, const size_t size, void *p, size_t *offset) {
  memcpy(p, map + *offset, size);
  *offset += size;

  return 0;
};

int io_write_data(void *map, const size_t size, const void *p, size_t *offset) {
  memcpy(map + *offset, p, size);
  *offset += size;

  return 0;
};
