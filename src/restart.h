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
#ifndef SWIFT_RESTART_H
#define SWIFT_RESTART_H

#include <stdio.h>

struct engine;

void restart_write(struct engine *e, const char *filename);
void restart_read(struct engine *e, const char *filename);

char **restart_locate(const char *dir, const char *basename, int *nfiles);
void restart_locate_free(int nfiles, char **files);
int restart_genname(const char *dir, const char *basename, int nodeID,
                    char *name, int size);

void restart_read_blocks(void *ptr, size_t size, size_t nblocks, FILE *stream,
                         char *label, const char *errstr);
void restart_write_blocks(void *ptr, size_t size, size_t nblocks, FILE *stream,
                          const char *label, const char *errstr);

int restart_stop_now(const char *dir, int cleanup);

void restart_save_previous(const char *filename);
void restart_remove_previous(const char *filename);

void restart_resubmit(const char *command);

#endif /* SWIFT_RESTART_H */
