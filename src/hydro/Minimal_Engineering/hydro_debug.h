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
#ifndef SWIFT_MINIMAL_ENGINEERING_HYDRO_DEBUG_H
#define SWIFT_MINIMAL_ENGINEERING_HYDRO_DEBUG_H

/**
 * @file Minimal/hydro_debug.h
 * @brief Minimal weakly compressible implementation of SPH (Debugging routines)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term without switches is implemented. Based on Moris et al 1997
 */

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  printf(
      "\n "
      "x=[%.6g, %.6g, %.6g], v=[%.3g, %.3g, %.3g], \n "
      "v_full=[%.3g, %.3g, %.3g], a=[%.3g, %.3g, %.3g], \n "
      "m=%.3g, u=%.3g, P=%.3g, \n "
      "v_sig=%.3g, h=%.3g, wcount=%.3g, rho=%.3g, \n "
      "time_bin=%d wakeup=%d \n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], xp->v_full[0],
      xp->v_full[1], xp->v_full[2], p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->mass, p->u, hydro_get_comoving_pressure(p),
      p->force.v_sig, p->h,
      p->density.wcount, p->rho, p->time_bin,
      p->limiter_data.wakeup);
}

#endif /* SWIFT_MINIMAL_HYDRO_DEBUG_H */
