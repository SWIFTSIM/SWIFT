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

__attribute__((always_inline))
    INLINE static void hydro_debug_particle(struct part* p, struct xpart* xp) {
  printf(
      "x=[%.16e,%.16e,%.16e], "
      "v=[%.3e,%.3e,%.3e], a=[%.3e,%.3e,%.3e], volume=%.3e\n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], p->a_hydro[0],
      p->a_hydro[1], p->a_hydro[2], p->geometry.volume);
}
