/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FOF_CATALOGUE_IO_H
#define SWIFT_FOF_CATALOGUE_IO_H

/* Config parameters. */
#include <config.h>

#ifdef WITH_FOF

void write_fof_hdf5_catalogue(const struct fof_props *props,
                              long long num_groups, const struct engine *e);

#endif /* WITH_FOF */

#endif /* SWIFT_FOF_CATALOGUE_IO_H */
