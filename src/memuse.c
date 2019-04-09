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
extern int engine_current_step;

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

  /* Relative time of this action. */
  ticks dtic;

  /* Label associated with the memory. */
  char label[MEMUSE_MAXLAB + 1];
};

/* The log of allocations and frees. */
static struct memuse_log_entry *memuse_log = NULL;
static volatile size_t memuse_log_size = 0;
static volatile size_t memuse_log_count = 0;
static volatile size_t memuse_log_done = 0;

#define MEMUSE_INITLOG 1000000
static void memuse_log_reallocate(size_t ind) {

  if (ind == 0) {

    /* Need to perform initialization. Be generous. */
    if ((memuse_log = (struct memuse_log_entry *)malloc(
             sizeof(struct memuse_log_entry) * MEMUSE_INITLOG)) == NULL)
      error("Failed to allocate memuse log.");

    /* Last action. */
    memuse_log_size = MEMUSE_INITLOG;

  } else {
    struct memuse_log_entry *new_log;
    if ((new_log = (struct memuse_log_entry *)malloc(
             sizeof(struct memuse_log_entry) * memuse_log_size * 2)) == NULL)
      error("Failed to re-allocate memuse log.");

    /* Wait for all writes to the old buffer to complete. */
    while (memuse_log_done < memuse_log_size)
      ;

    /* Copy to new buffer. */
    memcpy(new_log, memuse_log,
           sizeof(struct memuse_log_entry) * memuse_log_count);
    free(memuse_log);
    memuse_log = new_log;

    /* Last action. */
    memuse_log_size *= 2;
  }
}

/**
 * @brief Log an allocation or deallocation of memory.
 *
 * @param label the label associated with the memory.
 * @param ptr the memory pointer.
 * @param allocated whether this is an allocation or deallocation.
 * @param size the size in byte of memory allocated, set to 0 when
 *             deallocating.
 */
void memuse_log_allocation(const char *label, void *ptr, int allocated,
                           size_t size) {
  size_t ind = atomic_inc(&memuse_log_count);

  /* If we are at the current size we need more space. */
  if (ind == memuse_log_size) memuse_log_reallocate(ind);

  /* Other threads wait for space. */
  while (ind > memuse_log_size)
    ;

  /* Record the log. */
  memuse_log[ind].rank = engine_rank;
  memuse_log[ind].step = engine_current_step;
  memuse_log[ind].allocated = allocated;
  memuse_log[ind].size = size;
  memuse_log[ind].ptr = ptr;
  strncpy(memuse_log[ind].label, label, MEMUSE_MAXLAB);
  memuse_log[ind].label[MEMUSE_MAXLAB] = '\0';
  memuse_log[ind].dtic = getticks() - clocks_start_ticks;
  atomic_inc(&memuse_log_done);
}

/**
 * @brief dump the log to a file and reset, if anything to dump.
 *
 * @param filename name of file for log dump.
 */
void memuse_log_dump(const char *filename) {

  /* Skip if nothing allocated this step. */
  if (memuse_log_count == 0) return;

  /* Open the output file. */
  FILE *fd;
  if ((fd = fopen(filename, "w")) == NULL)
    error("Failed to create memuse log file '%s'.", filename);

  /* Write a header. */
  fprintf(fd, "# Current use: %s\n", memuse_process(1));
  fprintf(fd, "# cpufreq: %lld\n", clocks_get_cpufreq());
  fprintf(fd, "# dtic adr rank step allocated label size\n");

  for (size_t k = 0; k < memuse_log_count; k++) {
    fprintf(fd, "%lld %p %d %d %d %s %zd\n", memuse_log[k].dtic,
            memuse_log[k].ptr, memuse_log[k].rank, memuse_log[k].step,
            memuse_log[k].allocated, memuse_log[k].label, memuse_log[k].size);
  }

  /* Clear the log. */
  memuse_log_count = 0;

  /* Close the file. */
  fflush(fd);
  fclose(fd);
}

/**
 * @brief dump the log for using the given rank to generate a standard
 *        name for the output. Used when exiting in error.
 *
 * @param rank the rank exiting in error.
 */
void memuse_log_dump_error(int rank) {
  char filename[60];
  sprintf(filename, "memuse-error-report-rank%d.txt", rank);
  memuse_log_dump(filename);
}

#endif /* SWIFT_MEMUSE_REPORTS */

/**
 * @brief parse the process /proc/self/statm file to get the process
 *        memory use (in KB). Top field in ().
 *
 * @param size     total virtual memory (VIRT/VmSize)
 * @param resident resident non-swapped memory (RES/VmRSS)
 * @param shared   shared (mmap'd) memory  (SHR, RssFile+RssShmem)
 * @param text     text (exe) resident set (CODE, note also includes data
 *                 segment, so is considered broken for Linux)
 * @param data     data+stack resident set (DATA, note also includes library,
 *                 so is considered broken for Linux)
 * @param library  library resident set (0 for Linux)
 * @param dirty    dirty pages (nDRT = 0 for Linux)
 */
void memuse_use(long *size, long *resident, long *shared, long *text,
                long *data, long *library, long *dirty) {

  /* Open the file. */
  FILE *file = fopen("/proc/self/statm", "r");
  if (file != NULL) {
    int nscan = fscanf(file, "%ld %ld %ld %ld %ld %ld %ld", size, resident,
                       shared, text, library, data, dirty);

    if (nscan == 7) {
      /* Convert pages into bytes. Usually 4096, but could be 512 on some
       * systems so take care in conversion to KB. */
      long sz = sysconf(_SC_PAGESIZE);
      *size *= sz;
      *resident *= sz;
      *shared *= sz;
      *text *= sz;
      *library *= sz;
      *data *= sz;
      *dirty *= sz;

      *size /= 1024;
      *resident /= 1024;
      *shared /= 1024;
      *text /= 1024;
      *library /= 1024;
      *data /= 1024;
      *dirty /= 1024;
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
  long shared;
  long text;
  long library;
  long data;
  long dirty;
  memuse_use(&size, &resident, &shared, &text, &data, &library, &dirty);

  if (inmb) {
    snprintf(buffer, 256,
             "VIRT = %.3f SHR = %.3f CODE = %.3f DATA = %.3f "
             "RES = %.3f (MB)",
             size / 1024.0, shared / 1024.0, text / 1024.0, data / 1024.0,
             resident / 1024.0);
  } else {
    snprintf(buffer, 256,
             "VIRT = %ld SHR = %ld CODE = %ld DATA = %ld "
             "RES = %ld (KB)",
             size, shared, text, data, resident);
  }
  return buffer;
}
