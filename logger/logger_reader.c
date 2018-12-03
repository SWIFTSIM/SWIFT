#include "logger_header.h"
#include "logger_io.h"

/**
 * @brief Reverse offset in dump file
 *
 * @param filename string filename of the dump file
 * @param verbose Verbose level
 */
void reverse_offset(char *filename, int verbose) {
  struct header h;

  /* open file */
  int fd;
  void *map;
  io_open_file(filename, &fd, &map);
  
  /* read header */
  header_read(&h, map);

  if (verbose > 0) {
    header_print(&h);
  }

  /* check offset direction */
  if (h.forward_offset) {
    error("Offset are already reversed");
  }

  /* compute file size */
  size_t sz;
  io_get_file_size(fd, &sz);
  
  size_t offset;

#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  if (verbose > 0) {
    printf("Check offsets...\n");
  }
  offset = h.offset_first;
  while (offset < sz) {
    tools_check_offset(&h, map, &offset);
  }
  if (verbose > 0) {
    printf("Check done\n");
  }
#endif

  /* reverse header offset */
  header_change_offset_direction(&h, map);

  offset = h.offset_first;

  /* reverse chunks */
  if (verbose > 0) {
    printf("Reversing offsets...\n");
  }
  while (offset < sz) {
    tools_reverse_offset(&h, map, &offset);
  }
  if (verbose > 0) {
    printf("Reversing done\n");
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  if (verbose > 0) {
    printf("Check offsets...\n");
  }
  offset = h.offset_first;
  while (offset < sz) {
    tools_check_offset(&h, map, &offset);
  }
  if (verbose > 0) {
    printf("Check done\n");
  }
#endif

  /* free internal variables */
  header_free(&h);

  io_close_file(&fd, &map);

}
