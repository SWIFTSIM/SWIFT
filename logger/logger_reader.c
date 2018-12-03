#include "logger_header.h"
#include "logger_io.h"

/**
 * @brief Reverse offset in dump file
 *
 * @param filename string filename of the dump file
 */
void reverse_offset(char *filename) {
  struct header h;

  /* open file */
  int fd;
  void *map;
  io_open_file(filename, &fd, &map);
  
  /* read header */
  header_read(&h, map);

  header_print(&h);

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
  printf("Check offsets...\n");
  offset = h.offset_first;
  while (offset < sz) {
    tools_check_offset(&h, map, &offset);
  }
  printf("Check done\n");
#endif

  /* reverse header offset */
  header_change_offset_direction(&h, map);

  offset = h.offset_first;

  /* reverse chunks */
  printf("Reversing offsets...\n");
  while (offset < sz) {
    tools_reverse_offset(&h, map, &offset);
  }
  printf("Reversing done\n");

#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  printf("Check offsets...\n");
  offset = h.offset_first;
  while (offset < sz) {
    tools_check_offset(&h, map, &offset);
  }
  printf("Check done\n");
#endif

  /* free internal variables */
  header_free(&h);

  io_close_file(&fd, &map);

}
