/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_STARS_DEBUG_H
#define SWIFT_EAGLE_STARS_DEBUG_H

__attribute__((always_inline)) INLINE static void stars_debug_particle(
    const struct spart* p) {
  printf(
      "x=[%.3e,%.3e,%.3e], "
      "v_full=[%.3e,%.3e,%.3e] p->mass=%.3e \n t_begin=%d, t_end=%d\n",
      p->x[0], p->x[1], p->x[2], p->v_full[0], p->v_full[1], p->v_full[2],
      p->mass, p->ti_begin, p->ti_end);
}

#endif /* SWIFT_EAGLE_STARS_DEBUG_H */
