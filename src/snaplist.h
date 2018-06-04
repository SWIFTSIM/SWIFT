/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_SNAPLIST_H
#define SWIFT_SNAPLIST_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cosmology.h"

enum SNAPLIST_TYPE {
  SNAPLIST_AGE,
  SNAPLIST_REDSHIFT,
  SNAPLIST_SCALE_FACTOR,
};

struct snaplist {
  double *times;
  size_t size;
  size_t max_size;
};


void snaplist_read_file(struct snaplist *snaplist, const char* filename, struct cosmology *cosmo,
			size_t max_size);
void snaplist_print(const struct snaplist *snaplist);
void snaplist_clean(struct snaplist *snaplist);
void snaplist_struct_dump(struct snaplist *list, FILE *stream);
void snaplist_struct_restore(struct snaplist * list, FILE *stream);

#endif // SWIFT_SNAPLIST_H
