/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Peter W. Draper (p.w.draper@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/**
 *  @file memuse.c
 *  @brief file of routines to report about memory use in SWIFT.
 *  Note reports are in KB.
 */

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

/* Local defines. */
#include "memuse.h"

/* Local includes. */
#include "atomic.h"
#include "clocks.h"
#include "engine.h"

#ifdef SWIFT_MEMUSE_REPORTS

/* Also recorded in logger. */
extern int engine_rank;
extern int engine_cstep;

/* Entry for logger of memory allocations and deallocations in a step. */
#define MEMUSE_MAXLAB 64
struct memuse_log_entry {

  /* Rank in action. */
  int rank;

  /* Step of action. */
  int step;

  /* Whether allocated or deallocated. */
  int allocated;

  /* Memory allocated in bytes. */
  size_t size;

  /* Address of memory. */
  void *ptr;

  /* Time of this action. */
  ticks tic;

  /* Label associated with the memory. */
  char label[MEMUSE_MAXLAB + 1];
};

/* The log of allocations and frees. */
static struct memuse_log_entry *memuse_log = NULL;
static size_t log_size = 0;
static size_t log_count = 0;
static size_t log_done = 0;

#define MEMUSE_INITLOG 1000000
static void log_reallocate(size_t ind) {

  if (ind == 0) {

    /* Need to perform initialization. Be generous. */
    if ((memuse_log = (struct memuse_log_entry *)
         malloc(sizeof(struct memuse_log_entry) * MEMUSE_INITLOG)) == NULL)
      error("Failed to allocate memuse log.");

    /* Last action. */
    log_size = MEMUSE_INITLOG;

  } else {
    struct memuse_log_entry *new_log;
    if ((new_log = (struct memuse_log_entry *)
         malloc(sizeof(struct memuse_log_entry) * log_size * 2)) == NULL)
      error("Failed to re-allocate memuse log.");

    /* Wait for all writes to the old buffer to complete. */
    while (log_done < log_size);

    /* Copy to new buffer. */
    memcpy(new_log, memuse_log, sizeof(struct memuse_log_entry) * log_count);
    free(memuse_log);
    memuse_log = new_log;

    /* Last action. */
    log_size *= 2;
  }
}

static void log_allocation(const char *label, void *ptr, int allocated,
                           size_t size) {
  size_t ind = atomic_inc(&log_count);
  if (ind == log_size) log_reallocate(ind);

  /* Other threads wait for space. */
  while (ind > log_size);

  /* Record the log. */
  memuse_log[ind].rank = engine_rank;
  memuse_log[ind].step = engine_cstep;
  memuse_log[ind].allocated = allocated;
  memuse_log[ind].size = size;
  memuse_log[ind].ptr = ptr;
  strncpy(memuse_log[ind].label, label, MEMUSE_MAXLAB);
  memuse_log[ind].label[MEMUSE_MAXLAB] = '\0';
  memuse_log[ind].tic = getticks();
  atomic_inc(&log_done);
}

/**
 * @brief dump the log to a file and reset, if anything to dump.
 */
void memuse_log_dump(const char *filename) {

  /* Skip if nothing allocated this step. */
  if (log_count == 0) return;

  /* Open the output file. */
  FILE *fd;
  if ((fd = fopen(filename, "w")) == NULL)
    error("Failed to create memuse log file '%s'.", filename);

  /* Write a header. */
  fprintf(fd, "# Current use: %s\n", memuse_process(1));
  fprintf(fd, "# cpufreq: %lld\n", clocks_get_cpufreq());
  fprintf(fd, "# tic adr rank step allocated label size\n");

  for (size_t k = 0; k < log_count; k++) {
    fprintf(fd, "%lld %p %d %d %d %s %zd\n", memuse_log[k].tic,
            memuse_log[k].ptr, memuse_log[k].rank, memuse_log[k].step,
            memuse_log[k].allocated, memuse_log[k].label, memuse_log[k].size);
  }

  /* Clear the log. */
  log_count = 0;

  /* Close the file. */
  fflush(fd);
  fclose(fd);
}

#endif /* SWIFT_MEMUSE_REPORTS */

/**
 * @brief allocate aligned memory. The use and results are the same as the
 *        posix_memalign function.
 *
 * @param label a symbolic label for the memory, i.e. "parts".
 * @param memptr pointer to the allocated memory.
 * @param alignment alignment boundary.
 * @param size the quantity of bytes to allocate.
 * @result zero on success, otherwise an error code.
 */
int swift_memalign(const char *label, void **memptr, size_t alignment,
                   size_t size) {
  int result = posix_memalign(memptr, alignment, size);

#ifdef SWIFT_MEMUSE_REPORTS
  if (result == 0) log_allocation(label, *memptr, 1, size);
#endif
  return result;
}

/**
 * @brief free aligned memory. The use and results are the same as the
 *        free function.
 *
 * @param label a symbolic label for the memory, i.e. "parts", should match
 *              call used to allocate the memory.
 * @param ptr pointer to the allocated memory.
 */
void swift_free(const char *label, void *ptr) {
  free(ptr);

#ifdef SWIFT_MEMUSE_REPORTS
  log_allocation(label, ptr, 0, 0);
#endif
  return;
}

/**
 * @brief parse the process /proc/self/statm file to get the process
 *        memory use (in KB). Top field in ().
 *
 * @param size     total virtual memory (VIRT)
 * @param resident resident non-swapped memory (RES)
 * @param share    shared (mmap'd) memory  (SHR)
 * @param trs      text (exe) resident set (CODE)
 * @param lrs      library resident set
 * @param drs      data+stack resident set (DATA)
 * @param dt       dirty pages (nDRT)
 */
void memuse_use(long *size, long *resident, long *share, long *trs, long *lrs,
                long *drs, long *dt) {

  /* Open the file. */
  FILE *file = fopen("/proc/self/statm", "r");
  if (file != NULL) {
    int nscan = fscanf(file, "%ld %ld %ld %ld %ld %ld %ld", size, resident,
                       share, trs, lrs, drs, dt);

    if (nscan == 7) {
      /* Convert pages into bytes. Usually 4096, but could be 512 on some
       * systems so take care in conversion to KB. */
      long sz = sysconf(_SC_PAGESIZE);
      *size *= sz;
      *resident *= sz;
      *share *= sz;
      *trs *= sz;
      *lrs *= sz;
      *drs *= sz;
      *dt *= sz;

      *size /= 1024;
      *resident /= 1024;
      *share /= 1024;
      *trs /= 1024;
      *lrs /= 1024;
      *drs /= 1024;
      *dt /= 1024;
    } else {
      error("Failed to read sufficient fields from /proc/self/statm");
    }
    fclose(file);
  } else {
    error("Failed to open /proc/self/statm");
  }
}

/**
 * @brief Return a string with the current memory use of the process described.
 *
 * Not thread safe.
 *
 * @param inmb if true then report in MB, not KB.
 *
 * @result the memory use of the process, note make a copy if not used
 * immediately.
 */
const char *memuse_process(int inmb) {
  static char buffer[256];
  long size;
  long resident;
  long share;
  long trs;
  long lrs;
  long drs;
  long dt;
  memuse_use(&size, &resident, &share, &trs, &lrs, &drs, &dt);

  if (inmb) {
      snprintf(buffer, 256,
               "VIRT = %f SHR = %f CODE = %f DATA = %f "
               "RES = %f (MB)",
               size/1024.0, share/1024.0, trs/1024.0, drs/1024.0, resident/1024.0);
  } else {
      snprintf(buffer, 256,
               "VIRT = %ld SHR = %ld CODE = %ld DATA = %ld "
               "RES = %ld (KB)",
               size, share, trs, drs, resident);
  }
  return buffer;
}
