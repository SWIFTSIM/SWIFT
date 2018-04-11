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
#include <sys/types.h>
#include <unistd.h>

/* Local defines. */
#include "memuse.h"

/* Local includes. */
#include "clocks.h"
#include "engine.h"

#ifdef SWIFT_MEMUSE_REPORTS

/* For output messages. */
extern int engine_rank;
extern int engine_cstep;

/**
 *  @brief Report a memory allocation or use in bytes.
 *
 *  @param what a name for the report, "parts", "gparts" etc.
 *  @param file caller code file name
 *  @param function caller code function.
 *  @param line caller code file line no.
 *  @param bytes the number of bytes that have been allocated
 */
void memuse_report__(const char *what, const char *file, const char *function,
                     int line, size_t bytes) {
  char buffer[32];
  sprintf(buffer, "%zd", bytes / 1024);
  memuse_report_str__(what, file, function, line, buffer);
}

/**
 *  @brief Report a memory allocation or use formatted description.
 *
 *  @param what a name for the report, "parts", "gparts" etc.
 *  @param file caller code file name
 *  @param function caller code function.
 *  @param line caller code file line no.
 *  @param description the report, values should be in KB, the
 *         result of memuse_process() is suitable.
 */
void memuse_report_str__(const char *what, const char *file,
                         const char *function, int line,
                         const char *description) {
#ifdef WITH_MPI
  printf("[%04i] %s :memuse: %i:%s %s:%s:%i %s\n", engine_rank,
         clocks_get_timesincestart(), engine_cstep, what, file, function, line,
         description);
#else
  printf("%s :memuse: %i:%s %s:%s:%i %s\n", clocks_get_timesincestart(),
         engine_cstep, what, file, function, line, description);
#endif
}

#endif /* SWIFT_MEMUSE_REPORTS */

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
 * @result the memory use of the process, note make a copy if not used
 * immediately.
 */
const char *memuse_process() {
  static char buffer[256];
  long size;
  long resident;
  long share;
  long trs;
  long lrs;
  long drs;
  long dt;
  memuse_use(&size, &resident, &share, &trs, &lrs, &drs, &dt);

  snprintf(buffer, 256,
           "VIRT = %ld SHR = %ld CODE = %ld DATA = %ld "
           "RES = %ld",
           size, share, trs, drs, resident);
  return buffer;
}
