#include "timeline.h"
#include "io.h"

/**
 * @brief convert an integer time to a real time
 *
 * @param ti integer time to convert
 * @param timeBase time base of the integer time
 *
 * @return converted time
 */
double time_convert_to_double(const integertime_t ti, const double timeBase) {
  return ti * timeBase;  // should add timebegin
}

/**
 * @brief convert an integer time to a real time
 *
 * @param ti integer time to convert
 * @param timeBase time base of the integer time
 *
 * @return converted time
 */
integertime_t time_convert_to_integer(const double ti, const double timeBase) {
  return ti / timeBase;  // should add timebegin
}

/**
 * @brief read a time stamp
 *
 * @param timestamp timestamp read
 * @param time time read
 * @param h #header file structure
 * @param map file mapping
 * @param offset In: position in the file, Out: shifted by the timestamp
 */
int time_read(integertime_t *timestamp, double *time, const struct header *h,
              void *map, size_t *offset) {
  int error_code = 0;
  size_t mask = 0;
  size_t prev_offset = 0;
  *timestamp = 0;
  *time = 0;

  /* read chunck header */
  error_code = io_read_mask(h, map, offset, &mask, &prev_offset);
  if (error_code != 0) return error_code;

#ifdef SWIFT_DEBUG_CHECKS
  size_t ind = 0;

  /* check if timestamp is present */
  if (!header_is_present_and_get_index(h, "timestamp", &ind))
    error(EIO, "Header does not contain a timestamp");

  /* check if timestamp */
  if (h->masks[ind] != mask) error(EIO, "Not a timestamp");
#endif

  /* read data */
  // TODO
  error_code =
      io_read_data(map, sizeof(unsigned long long int), timestamp, offset);
  error_code = io_read_data(map, sizeof(double), time, offset);

  return error_code;
}

/**
 * @brief get offset of first timestamp
 *
 * @param h file #header
 * @param map file mapping
 * @param offset out: offset of first timestamp
 *
 * @return error code
 */
int time_first_timestamp(const struct header *h, void *map, size_t *offset) {
  *offset = h->offset_first;
  int test;

  size_t i;

  if (!header_is_present_and_get_index(h, "timestamp", &i))
    error(EIO, "Time stamp not present in header");

  size_t tmp = *offset;
  size_t mask = 0;
  test = io_read_mask(h, map, offset, &mask, NULL);
  if (test != 0) return test;

  if (mask != h->masks[i]) error(EIO, "Dump should begin by timestep");

  *offset = tmp;
  return 0;
}

/**
 * @brief Initialize a time array
 *
 * @param t #time_array to initialize
 * @param h #header file structure
 * @param map file mapping
 * @param fd file id
 */
int time_array_init(struct time_array *t, const struct header *h, void *map,
                    int fd) {

  t->next = NULL;
  t->prev = NULL;

  /* get first time stamp */
  size_t offset = 0;
  time_first_timestamp(h, map, &offset);

  integertime_t timestamp = 0;
  double time = 0;

  /* get file size */
  size_t file_size;
  int test = io_get_file_size(fd, &file_size);

  while (offset < file_size) {

    /* read time */
    t->offset = offset;
    size_t tmp_offset = offset;
    time_read(&timestamp, &time, h, map, &tmp_offset);
    t->timestamp = timestamp;
    t->time = time;

    /* get next chunk */
    test = tools_get_next_chunk(h, map, &offset, fd);
    if (test == -1)
      break;
    else if (test != 0)
      return test;

    /* allocate next time_array */
    struct time_array *tmp = malloc(sizeof(struct time_array));
    tmp->prev = t;
    tmp->next = NULL;
    t->next = tmp;
    t = tmp;
  }

  /* free unused time_array */
  struct time_array *tmp = t->prev;
  tmp->next = NULL;
  free(t);

  return 0;
}

/**
 * @brief access the time of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return integer time of the chunk
 */
integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset) {
  const struct time_array *tmp = time_array_get_time_array(t, offset);
  return tmp->timestamp;
}

/**
 * @brief access the time of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return time of the chunk
 */
double time_array_get_time(struct time_array *t, const size_t offset) {
  const struct time_array *tmp = time_array_get_time_array(t, offset);
  return tmp->time;
}

/**
 * @brief access the #time_array of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return pointer to the requested #time_array
 */
struct time_array *time_array_get_time_array(struct time_array *t,
                                             const size_t offset) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!t) return 0;
#endif
  while (t->next) {
    if (t->offset > offset) break;

    t = t->next;
  }

  if (t->prev == NULL) return t;

  struct time_array *tmp = t->prev;
  return tmp;
}

/**
 * @brief free memory of a #time_array
 *
 * @param t #time_array to free
 */
void time_array_free(struct time_array *t) {
  struct time_array *tmp;
  while (t->next) {
    tmp = t->next;
    free(t);
    t = tmp;
  }
}

/**
 * @brief print a #time_array
 *
 * @param t #time_array to print
 */
void time_array_print(const struct time_array *t) {
  size_t threshold = 4;

  size_t n = time_array_count(t);
  size_t i = 0;
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lli (%g)", n, t->timestamp, t->time);

  while (t->next) {
    i += 1;
    t = t->next;
    if (i < threshold || i > up_threshold)
      printf(", %lli (%g)", t->timestamp, t->time);

    if (i == threshold) printf(", ...");
  }

  printf("]\n");
}

/**
 * @brief print a #time_array (offset)
 *
 * @param t #time_array to print
 */
void time_array_print_offset(const struct time_array *t) {
  size_t threshold = 4;

  size_t n = time_array_count(t);
  size_t i = 0;
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lu", n, t->offset);

  while (t->next) {
    i += 1;
    t = t->next;
    if (i < threshold || i > up_threshold) printf(", %lu", t->offset);

    if (i == threshold) printf(", ...");
  }

  printf("]\n");
}

/**
 * @brief count number of element in #time_array
 *
 * @param t #time_array to count
 *
 * @return number of element
 */
size_t time_array_count(const struct time_array *t) {
  size_t count = 1;
  while (t->next) {
    t = t->next;
    count += 1;
  }

  return count;
}
