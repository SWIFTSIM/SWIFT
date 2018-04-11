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
#ifndef SWIFT_MEMUSE_H
#define SWIFT_MEMUSE_H

/* Config parameters. */
#include "../config.h"

/* Public API. */
void memuse_use(long *size, long *resident, long *share, long *trs, long *lrs,
                long *drs, long *dt);
const char *memuse_process();

/* Reports are a no-op unless wanted. */
#ifdef SWIFT_MEMUSE_REPORTS
void memuse_report__(const char *what, const char *file, const char *function,
                     int line, size_t bytes);
void memuse_report_str__(const char *what, const char *file,
                         const char *function, int line,
                         const char *description);

#define memuse_report(what, size) \
  memuse_report__(what, __FILE__, __FUNCTION__, __LINE__, size)
#define memuse_report_str(what, description) \
  memuse_report_str__(what, __FILE__, __FUNCTION__, __LINE__, description)
#else
#define memuse_report(what, size)
#define memuse_report_str(what, description)
#endif

#endif /* SWIFT_MEMUSE_H */
