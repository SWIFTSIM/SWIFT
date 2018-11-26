#ifndef __TIMELINE_H__
#define __TIMELINE_H__

#include "header.h"
#include "logger_tools.h"

typedef char timebin_t;
typedef long long integertime_t;

struct time_array {
  /* Pointer to next element */
  void *next;

  /* Pointer to prev element */
  void *prev;

  /* Integertime of this timestamp */
  integertime_t timestamp;

  /* Double time of this timestamp */
  double time;

  /* Offset in the file of this timestamp */
  size_t offset;
};

double time_convert_to_double(const integertime_t ti, const double timeBase);
integertime_t time_convert_to_integer(const double ti, const double timeBase);
int time_read(integertime_t *timestamp, double *time, const struct header *h,
              void *map, size_t *offset);
int time_array_init(struct time_array *t, const struct header *h, void *map,
                    int fd);
integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset);
double time_array_get_time(struct time_array *t, const size_t offset);
struct time_array *time_array_get_time_array(struct time_array *t,
                                             const size_t offset);
void time_array_free(struct time_array *t);
void time_array_print(const struct time_array *t);
void time_array_print_offset(const struct time_array *t);
size_t time_array_count(const struct time_array *t);
int time_first_timestamp(const struct header *h, void *map, size_t *offset);

#endif  // __TIMELINE_H__
