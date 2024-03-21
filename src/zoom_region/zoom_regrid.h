/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will J. Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_ZOOM_REGRID_H
#define SWIFT_ZOOM_REGRID_H

/* Local includes */
#include "space.h"

int zoom_need_regrid(const struct space *s, const int new_cdim[3]);
void zoom_prepare_cells(struct space *s, const int zoom_cdim[3], int verbose);
void zoom_allocate_cells(struct space *s);

#endif /* SWIFT_ZOOM_REGRID_H */
