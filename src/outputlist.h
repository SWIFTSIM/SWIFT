/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausamman (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_OUTPUTLIST_H
#define SWIFT_OUTPUTLIST_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cosmology.h"

enum OUTPUTLIST_TYPE {
  OUTPUTLIST_AGE,
  OUTPUTLIST_REDSHIFT,
  OUTPUTLIST_SCALE_FACTOR,
};

struct outputlist {
  double *times;
  size_t size;
};

void outputlist_read_file(struct outputlist *outputlist, const char *filename,
                          struct cosmology *cosmo);
void outputlist_print(const struct outputlist *outputlist);
void outputlist_clean(struct outputlist *outputlist);
void outputlist_struct_dump(struct outputlist *list, FILE *stream);
void outputlist_struct_restore(struct outputlist *list, FILE *stream);

#endif /* SWIFT_OUTPUTLIST_H */
