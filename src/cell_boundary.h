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

/**
 * @brief Has the particle left the box through an open face?
 *
 * Reflective walls never count: particles crossing them are mirrored back.
 * The upper face retains its old convention.
 *
 * @param x The particle position.
 * @param dim The box dimensions.
 * @param walls Bit mask of reflective faces (see #boundary_wall).
 */
static inline int cell_boundary_outside(const double x[3], const double dim[3],
                                        const unsigned int walls) {
  for (int k = 0; k < 3; k++) {
    if (x[k] < 0. && !(walls & (1u << (2 * k)))) return 1;
    if (x[k] > dim[k] && !(walls & (1u << (2 * k + 1)))) return 1;
  }
  return 0;
}

/**
 * @brief Mirror a position back across the reflective walls it crossed.
 *
 * @param x The particle position (updated).
 * @param dim The box dimensions.
 * @param walls Bit mask of reflective faces (see #boundary_wall).
 * @param shift (return) The displacement applied along each axis.
 *
 * @return 1 if the particle was reflected, 0 otherwise.
 */
static inline int cell_boundary_reflect(double x[3], const double dim[3],
                                        const unsigned int walls,
                                        double shift[3]) {
  int reflected = 0;
  for (int k = 0; k < 3; k++) {
    shift[k] = 0.;
    if (x[k] < 0. && (walls & (1u << (2 * k))))
      shift[k] = -2. * x[k];
    else if (x[k] > dim[k] && (walls & (1u << (2 * k + 1))))
      shift[k] = 2. * (dim[k] - x[k]);
    else
      continue;
    x[k] += shift[k];
    reflected = 1;
  }
  return reflected;
}

void cell_boundary_part(struct cell *c, const struct engine *e);
void cell_boundary_gpart(struct cell *c, const struct engine *e);
void cell_boundary_spart(struct cell *c, const struct engine *e);
void cell_boundary_bpart(struct cell *c, const struct engine *e);
void cell_boundary_sink(struct cell *c, const struct engine *e);

#endif /* SWIFT_CELL_BOUNDARY_H */
