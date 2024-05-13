/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
 *               2024 Will J. Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_ZOOM_H
#define SWIFT_ZOOM_H

/* Local includes */
#include "timeline.h"

/* Avoid cyclic inclusions. */
struct swift_params;
struct space;
struct cell;
struct engine;

/* Zoom region and cell grid initialisation functions. */
void zoom_props_init(struct swift_params *params, struct space *s,
                     const int verbose);
void zoom_region_init(struct space *s, const int verbose);

/* Construct top level cells with a zoom region. */
void zoom_construct_tl_cells(struct space *s, const integertime_t ti_current,
                             int verbose);

/* Linking zoom cells to void leaves. */
void zoom_link_void_leaves(struct space *s, struct cell *c);

/* Space regridding functions. */
int zoom_need_regrid(const struct space *s, const int new_cdim[3]);
void zoom_prepare_cells(struct space *s, const int zoom_cdim[3], int verbose);
void zoom_allocate_cells(struct space *s);

/* Void cell tree construction functions. */
void zoom_void_space_split(struct space *s, int verbose);

/* Task creation functions. */
void zoom_engine_make_self_gravity_tasks(struct space *s, struct engine *e);
void zoom_engine_make_void_gravity_tasks(struct space *s, struct engine *e);

/* Zoom specific unskip. */
void zoom_unskip_void_tasks(struct scheduler *s);

#endif /* SWIFT_ZOOM_H */