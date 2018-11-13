#if defined(HAVE_SETAFFINITY) && defined(SWIFT_DEBUG_TASKS) && defined( WITH_PERF )

#ifndef SWIFT_PERF_H
#define SWIFT_PERF_H
#define _GNU_SOURCE
#include <asm/unistd.h>

static long swift_perf_event_open(struct perf_event_attr *hw_event, pid_t pit,
                                  int cpu, int group_fd, unsigned_long flags){
   return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);

}

#undef _GNU_SOURCE

#endif
#endif
