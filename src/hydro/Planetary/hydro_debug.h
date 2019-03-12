/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_HYDRO_DEBUG_H
#define SWIFT_PLANETARY_HYDRO_DEBUG_H

/**
 * @file Planetary/hydro_debug.h
 * @brief Minimal conservative implementation of SPH (Debugging routines)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term without switches is implemented. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of
 * Price, D., Journal of Computational Physics, 2012, Volume 231, Issue 3,
 * pp. 759-794.
 */

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  printf(
      "\n "
      "x=[%.6g, %.6g, %.6g], v=[%.3g, %.3g, %.3g], \n "
      "v_full=[%.3g, %.3g, %.3g], a=[%.3g, %.3g, %.3g], \n "
      "m=%.3g, u=%.3g, du/dt=%.3g, P=%.3g, c_s=%.3g, \n "
      "v_sig=%.3g, h=%.3g, dh/dt=%.3g, wcount=%.3g, rho=%.3g, \n "
      "dh_drho=%.3g, time_bin=%d, wakeup=%d mat_id=%d \n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], xp->v_full[0],
      xp->v_full[1], xp->v_full[2], p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->mass, p->u, p->u_dt, hydro_get_comoving_pressure(p),
      p->force.soundspeed, p->force.v_sig, p->h, p->force.h_dt,
      p->density.wcount, p->rho, p->density.rho_dh, p->time_bin, p->wakeup,
      p->mat_id);
}

#endif /* SWIFT_PLANETARY_HYDRO_DEBUG_H */
