/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#ifndef SWIFT_HDF5_OBJECT_TO_BLOB_H
#define SWIFT_HDF5_OBJECT_TO_BLOB_H

/* Config parameters. */
#include "../config.h"

#ifdef HAVE_HDF5

#include <hdf5.h>

void hdf5_object_to_blob(hid_t group_id, char *name, size_t *len, void **data);
void blob_to_hdf5_object(size_t len, void *data, hid_t dest_id, char *name);

#endif /* HAVE_HDF5 */

#endif /* SWIFT_HDF5_OBJECT_TO_BLOB_H */
