/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
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

#ifndef SWIFT_GASOLINE_HYDRO_DEBUG_H
#define SWIFT_GASOLINE_HYDRO_DEBUG_H

/**
 * @file Gasoline/hydro_debug.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added Gasoline physics (Wadsley+ 2017) (Debugging routines)
 */

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  warning("[PID%lld] part:", p->id);
  warning("[PID%lld] x=[%.3e,%.3e,%.3e]", p->id, p->x[0], p->x[1], p->x[2]);
  warning("[PID%lld] v=[%.3e,%.3e,%.3e]", p->id, p->v[0], p->v[1], p->v[2]);
  warning("[PID%lld] a=[%.3e,%.3e,%.3e]", p->id, p->a_hydro[0], p->a_hydro[1],
          p->a_hydro[2]);
  warning("[PID%lld] u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e", p->id, p->u,
          p->u_dt, p->viscosity.v_sig, hydro_get_comoving_pressure(p));
  warning("[PID%lld] h=%.3e, dh/dt=%.3e wcount=%d, m=%.3e, dh_drho=%.3e", p->id,
          p->h, p->force.h_dt, (int)p->density.wcount, p->mass,
          p->density.rho_dh);
  warning(
      "[PID%lld] alpha=%.3e, time_bin=%d, rho=%.3e, velocity_gradient=["
      "[%.3e,%.3e,%.3e],"
      "[%.3e,%.3e,%.3e],"
      "[%.3e,%.3e,%.3e]]",
      p->id, p->viscosity.alpha, p->time_bin, p->rho,
      p->viscosity.velocity_gradient[0][0],
      p->viscosity.velocity_gradient[0][1],
      p->viscosity.velocity_gradient[0][2],
      p->viscosity.velocity_gradient[1][0],
      p->viscosity.velocity_gradient[1][1],
      p->viscosity.velocity_gradient[1][2],
      p->viscosity.velocity_gradient[2][0],
      p->viscosity.velocity_gradient[2][1],
      p->viscosity.velocity_gradient[2][2]);
  warning("[PID%lld] smooth_pressure_gradient=[%.3e,%.3e,%.3e]", p->id,
          p->smooth_pressure_gradient[0], p->smooth_pressure_gradient[1],
          p->smooth_pressure_gradient[2]);
  warning("[PID%lld] weighted_wcount=%.3e", p->id, p->weighted_wcount);

  if (xp != NULL) {
    warning("[PID%lld] xpart:", p->id);
    warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]", p->id, xp->v_full[0],
            xp->v_full[1], xp->v_full[2]);
  }
}

#endif /* SWIFT_GASOLINE_HYDRO_DEBUG_H */
