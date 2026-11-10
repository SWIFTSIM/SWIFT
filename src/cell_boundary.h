/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 SWIFT developers.
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
#ifndef SWIFT_CELL_BOUNDARY_H
#define SWIFT_CELL_BOUNDARY_H

struct cell;
struct engine;

/** Is the particle outside the box? The upper face retains its old convention. */
static inline int cell_boundary_outside(const double x[3],
                                        const double dim[3]) {
  return x[0] < 0. || x[0] > dim[0] || x[1] < 0. || x[1] > dim[1] ||
         x[2] < 0. || x[2] > dim[2];
}

void cell_boundary_part(struct cell *c, const struct engine *e);
void cell_boundary_gpart(struct cell *c, const struct engine *e);
void cell_boundary_spart(struct cell *c, const struct engine *e);
void cell_boundary_bpart(struct cell *c, const struct engine *e);
void cell_boundary_sink(struct cell *c, const struct engine *e);

#endif /* SWIFT_CELL_BOUNDARY_H */
