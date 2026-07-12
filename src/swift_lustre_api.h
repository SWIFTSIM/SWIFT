/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_LUSTRE_API_H
#define SWIFT_LUSTRE_API_H

/* For size_t and FILE. */
#include <stdlib.h>

/* Structure to store information about an OST. */
struct swift_ost_info {
  int index;   /* OST index */
  size_t size; /* Size in bytes */
  size_t used; /* Used in bytes */
};

/* Structure to store a scan of all the OSTs for a mount point. */
struct swift_ost_store {
  struct swift_ost_info *infos;
  int count;     /* Count of active OSTs */
  int fullcount; /* Count of OSTs available (only different when culled) */
  int size;      /* Space available for storing OST infos */
};

/* Public functions. */

/* OST scanning and selection. */
void swift_ost_select(struct swift_ost_store *ost_infos, const char *path,
                      int minfree, int writetest, int verbose);
int swift_ost_scan(const char *path, struct swift_ost_store *ost_infos);
void swift_ost_cull(struct swift_ost_store *ost_infos, int minfree);
void swift_ost_remove(struct swift_ost_store *ost_infos, int index);
int swift_ost_next(struct swift_ost_store *ost_infos, int *arrayindex,
                   int count);

/* OST store. */
void swift_ost_store_init(struct swift_ost_store *ost_infos);
void swift_ost_store_alloc(struct swift_ost_store *ost_infos, int size);
void swift_ost_store_copy(struct swift_ost_store *ost_infos_src,
                          struct swift_ost_store *ost_infos_dst);
void swift_ost_store_free(struct swift_ost_store *ost_infos);
void swift_ost_store_print(struct swift_ost_store *ost_infos, int verbose);
void swift_ost_store_write(FILE *file, struct swift_ost_store *ost_infos);

/* File striping. */
int swift_create_striped_file(const char *filename, int offset, int count,
                              int *usedoffset);

#endif /* SWIFT_LUSTRE_API_H */
