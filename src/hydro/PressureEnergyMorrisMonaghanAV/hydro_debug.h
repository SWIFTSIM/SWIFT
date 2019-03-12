/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk) &
 *                    Josh Borrow (joshua.borrow@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENERGY_MORRIS_HYDRO_DEBUG_H
#define SWIFT_PRESSURE_ENERGY_MORRIS_HYDRO_DEBUG_H
/**
 * @file PressureEnergy/hydro_debug.h
 * @brief P-U conservative implementation of SPH (Debugging routines)
 *
 * The thermal variable is the internal energy (u). A simple variable
 * viscosity term (Morris & Monaghan 1997) with a Balsara switch is
 * implemented.
 */

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  printf(
      "x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e],v_full=[%.3e,%.3e,%.3e] \n a=[%.3e,%.3e,%.3e], "
      "u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e\n"
      "h=%.3e, dh/dt=%.3e wcount=%d, m=%.3e, dh_drho=%.3e, rho=%.3e, \n"
      "p_dh=%.3e, p_bar=%.3e \n"
      "time_bin=%d, wakeup=%d alpha=%.3e\n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], xp->v_full[0],
      xp->v_full[1], xp->v_full[2], p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->u, p->u_dt, p->force.v_sig, hydro_get_comoving_pressure(p), p->h,
      p->force.h_dt, (int)p->density.wcount, p->mass, p->density.rho_dh, p->rho,
      p->density.pressure_bar_dh, p->pressure_bar, p->time_bin, p->wakeup,
      p->alpha);
}

#endif /* SWIFT_PRESSURE_ENERGY_MORRIS_HYDRO_DEBUG_H */
