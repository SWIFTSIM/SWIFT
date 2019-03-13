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
int swift_memalign(const char *label, void **memptr, size_t alignment,
                   size_t size);
void swift_free(const char *label, void *ptr);

void memuse_use(long *size, long *resident, long *share, long *trs, long *lrs,
                long *drs, long *dt);
const char *memuse_process(void);

#ifdef SWIFT_MEMUSE_REPORTS
void memuse_log_dump(const char *filename);
#endif

#endif /* SWIFT_MEMUSE_H */
