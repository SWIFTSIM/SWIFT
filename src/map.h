/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Includes. */
#ifndef SWIFT_MAP_H
#define SWIFT_MAP_H

#include "cell.h"
#include "part.h"

void map_cells_plot(struct cell *c, void *data);
void map_check(struct part *p, struct cell *c, void *data);
void map_cellcheck(struct cell *c, void *data);
void map_maxdepth(struct cell *c, void *data);
void map_count(struct part *p, struct cell *c, void *data);
void map_wcount_min(struct part *p, struct cell *c, void *data);
void map_wcount_max(struct part *p, struct cell *c, void *data);
void map_h_min(struct part *p, struct cell *c, void *data);
void map_h_max(struct part *p, struct cell *c, void *data);
void map_stars_h_max(struct spart *p, struct cell *c, void *data);
void map_icount(struct part *p, struct cell *c, void *data);
void map_dump(struct part *p, struct cell *c, void *data);

#endif /* SWIFT_MAP_H */
