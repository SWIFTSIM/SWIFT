/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DEFAULT_GRAVITY_DEBUG_H
#define SWIFT_DEFAULT_GRAVITY_DEBUG_H

__attribute__((always_inline)) INLINE static void gravity_debug_particle(
    const struct gpart* p) {
  printf(
      "mass=%.3e time_bin=%d\n"
      "x=[%.5e,%.5e,%.5e], v_full=[%.5e,%.5e,%.5e], a=[%.5e,%.5e,%.5e]\n",
      p->mass, p->time_bin, p->x[0], p->x[1], p->x[2], p->v_full[0],
      p->v_full[1], p->v_full[2], p->a_grav[0], p->a_grav[1], p->a_grav[2]);
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  printf("pot=%e\n", p->potential);
#endif
#ifdef SWIFT_DEBUG_CHECKS
  printf("num_interacted=%lld ti_drift=%lld ti_kick=%lld\n", p->num_interacted,
         p->ti_drift, p->ti_kick);
#endif
}

#endif /* SWIFT_DEFAULT_GRAVITY_DEBUG_H */
