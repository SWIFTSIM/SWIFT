#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "logger_header.h"
#include "logger_io.h"
#include "logger_tools.h"

/**
 * @brief get the size of a file
 *
 * @param fd file id
 * @param size out: file size
 *
 * @return error code
 */
int io_get_file_size(int fd, size_t *size) {
  struct stat s;
  int status = fstat(fd, &s);
  if (status != 0) error(errno, "Unable to get file size");
  *size = s.st_size;

  return 0;
}

/**
 * @brief Open a file and map it
 *
 * @param filename file to read
 * @param fd out: file id
 * @param map out: file mapping
 *
 * @return error code
 */
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

/**
 * @brief Close a file and unmap it
 *
 * @param fd file id
 * @param map file mapping
 *
 * @return error code
 */
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
int io_read_data(void *map, const size_t size, void *p, size_t *offset) {
  memcpy(p, map + *offset, size);
  *offset += size;

  return 0;
};

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
int io_write_data(void *map, const size_t size, const void *p, size_t *offset) {
  memcpy(map + *offset, p, size);
  *offset += size;

  return 0;
};
