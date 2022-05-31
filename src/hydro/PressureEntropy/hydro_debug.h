/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PRESSURE_ENTROPY_HYDRO_DEBUG_H
#define SWIFT_PRESSURE_ENTROPY_HYDRO_DEBUG_H

/**
 * @file PressureEntropy/hydro_debug.h
 * @brief Pressure-Entropy implementation of SPH (Debugging routines)
 *
 * The thermal variable is the entropy (S) and the entropy is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows eqautions (19), (21) and (22) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */
__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  warning("[PID%lld] part:", p->id);
  warning(
      "[PID%lld] x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e], a=[%.3e,%.3e,%.3e], "
      "h=%.3e, wcount=%.3f, wcount_dh=%.3e, m=%.3e, dh_drho=%.3e, rho=%.3e, "
      "rho_bar=%.3e, P=%.3e, dP_dh=%.3e, P_over_rho2=%.3e, S=%.3e, S^1/g=%.3e, "
      "dS/dt=%.3e, c=%.3e v_sig=%e dh/dt=%.3e time_bin=%d wakeup=%d",
      p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
      p->a_hydro[0], p->a_hydro[1], p->a_hydro[2], p->h, p->density.wcount,
      p->density.wcount_dh, p->mass, p->density.rho_dh, p->rho, p->rho_bar,
      hydro_get_comoving_pressure(p), p->density.pressure_dh,
      p->force.P_over_rho2, p->entropy, p->entropy_one_over_gamma,
      p->entropy_dt, p->force.soundspeed, p->force.v_sig, p->force.h_dt,
      p->time_bin, p->limiter_data.wakeup);
  if (xp != NULL) {
    warning("[PID%lld] xpart:", p->id);
    warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]", p->id, xp->v_full[0],
            xp->v_full[1], xp->v_full[2]);
  }
}

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_DEBUG_H */
