#ifndef __TIMELINE_H__
#define __TIMELINE_H__

#include "header.h"
#include "tools.h"

typedef char timebin_t;
typedef long long integertime_t;

struct time_array {
  void *next;
  void *prev;
  integertime_t timestamp;
  double time;
  size_t offset;
};

/**
 * @brief convert an integer time to a real time
 *
 * @param ti integer time to convert
 * @param timeBase time base of the integer time
 *
 * @return converted time
 */
double time_convert_to_double(const integertime_t ti, const double timeBase);

/**
 * @brief convert an integer time to a real time
 *
 * @param ti integer time to convert
 * @param timeBase time base of the integer time
 *
 * @return converted time
 */
integertime_t time_convert_to_integer(const double ti, const double timeBase);

/**
 * @brief read a time stamp
 *
 * @param timestamp timestamp read
 * @param time time read
 * @param h #header file structure
 * @param map file mapping
 * @param offset In: position in the file, Out: shifted by the timestamp
 */
int time_read(integertime_t *timestamp, double *time, const struct header *h, void *map,
              size_t *offset);

/**
 * @brief Initialize a time array
 *
 * @param t #time_array to initialize
 * @param h #header file structure
 * @param map file mapping
 * @param fd file id
 */
int time_array_init(struct time_array *t, const struct header *h, void *map,
                    int fd);

/**
 * @brief access the time of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return integer time of the chunk
 */
integertime_t time_array_get_integertime(struct time_array *t, const size_t offset);

/**
 * @brief access the time of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return time of the chunk
 */
double time_array_get_time(struct time_array *t, const size_t offset);

/**
 * @brief access the #time_array of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return pointer to the requested #time_array
 */
struct time_array *time_array_get_time_array(struct time_array *t,
                                             const size_t offset);

/**
 * @brief free memory of a #time_array
 *
 * @param t #time_array to free
 */
void time_array_free(struct time_array *t);

/**
 * @brief print a #time_array
 *
 * @param t #time_array to print
 */
void time_array_print(const struct time_array *t);

/**
 * @brief print a #time_array (offset)
 *
 * @param t #time_array to print
 */
void time_array_print_offset(const struct time_array *t);

/**
 * @brief count number of element in #time_array
 *
 * @param t #time_array to count
 *
 * @return number of element
 */
size_t time_array_count(const struct time_array *t);

/**
 * @brief get offset of first timestamp
 *
 * @param h file #header
 * @param map file mapping
 * @param offset out: offset of first timestamp
 *
 * @return error code
 */
int time_first_timestamp(const struct header *h, void *map, size_t *offset);

#endif  // __TIMELINE_H__
