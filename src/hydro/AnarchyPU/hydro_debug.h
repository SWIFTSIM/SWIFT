/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl) &
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
#ifndef SWIFT_ANARCHY_PU_HYDRO_DEBUG_H
#define SWIFT_ANARCHY_PU_HYDRO_DEBUG_H
/**
 * @file AnarchyPU/hydro_debug.h
 * @brief P-U conservative implementation of SPH,
 *        with added ANARCHY physics (Cullen & Denhen 2011 AV,
 *        Price 2008 thermal diffusion (Debugging routines).
 */

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  warning("[PID%lld] part:", p->id);
  warning(
      "[PID%lld] x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e]",
      p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2]);
  warning(
      "[PID%lld] a=[%.3e,%.3e,%.3e], "
      "u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e",
      p->id, p->a_hydro[0], p->a_hydro[1], p->a_hydro[2], p->u, p->u_dt,
      p->viscosity.v_sig, hydro_get_comoving_pressure(p));
  warning(
      "[PID%lld] h=%.3e, dh/dt=%.3e wcount=%d, m=%.3e, dh_drho=%.3e, rho=%.3e",
      p->id, p->h, p->force.h_dt, (int)p->density.wcount, p->mass,
      p->density.rho_dh, p->rho);
  warning("[PID%lld] p_dh=%.3e, p_bar=%.3e, alpha=%.3e", p->id,
          p->density.pressure_bar_dh, p->pressure_bar, p->viscosity.alpha);
  warning("[PID%lld] time_bin=%d", p->id, p->time_bin);
  if (xp != NULL) {
    warning("[PID%lld] xpart:", p->id);
    warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]", p->id, xp->v_full[0],
            xp->v_full[1], xp->v_full[2]);
  }
}

#endif /* SWIFT_ANARCHY_PU_HYDRO_DEBUG_H */
