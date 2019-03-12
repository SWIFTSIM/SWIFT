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
#ifndef SWIFT_GIZMO_MFM_HYDRO_DEBUG_H
#define SWIFT_GIZMO_MFM_HYDRO_DEBUG_H

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  printf(
      "x=[%.16e,%.16e,%.16e], "
      "v=[%.3e,%.3e,%.3e], "
      "a=[%.3e,%.3e,%.3e], "
      "h=%.3e, "
      "time_bin=%d, "
      "wakeup=%d, "
      "rho=%.3e, "
      "P=%.3e, "
      "gradients={"
      "rho=[%.3e,%.3e,%.3e], "
      "v=[[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e]], "
      "P=[%.3e,%.3e,%.3e]}, "
      "limiter={"
      "rho=[%.3e,%.3e], "
      "v=[[%.3e,%.3e],[%.3e,%.3e],[%.3e,%.3e]], "
      "P=[%.3e,%.3e], "
      "maxr=%.3e}, "
      "conserved={"
      "momentum=[%.3e,%.3e,%.3e], "
      "mass=%.3e, "
      "energy=%.3e}, "
      "geometry={"
      "volume=%.3e, "
      "matrix_E=[[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e]]}, "
      "timestepvars={"
      "vmax=%.3e},"
      "density={"
      "wcount_dh=%.3e, "
      "wcount=%.3e}\n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], p->a_hydro[0],
      p->a_hydro[1], p->a_hydro[2], p->h, p->time_bin, p->wakeup, p->rho, p->P,
      p->gradients.rho[0], p->gradients.rho[1], p->gradients.rho[2],
      p->gradients.v[0][0], p->gradients.v[0][1], p->gradients.v[0][2],
      p->gradients.v[1][0], p->gradients.v[1][1], p->gradients.v[1][2],
      p->gradients.v[2][0], p->gradients.v[2][1], p->gradients.v[2][2],
      p->gradients.P[0], p->gradients.P[1], p->gradients.P[2],
      p->limiter.rho[0], p->limiter.rho[1], p->limiter.v[0][0],
      p->limiter.v[0][1], p->limiter.v[1][0], p->limiter.v[1][1],
      p->limiter.v[2][0], p->limiter.v[2][1], p->limiter.P[0], p->limiter.P[1],
      p->limiter.maxr, p->conserved.momentum[0], p->conserved.momentum[1],
      p->conserved.momentum[2], p->conserved.mass, p->conserved.energy,
      p->geometry.volume, p->geometry.matrix_E[0][0],
      p->geometry.matrix_E[0][1], p->geometry.matrix_E[0][2],
      p->geometry.matrix_E[1][0], p->geometry.matrix_E[1][1],
      p->geometry.matrix_E[1][2], p->geometry.matrix_E[2][0],
      p->geometry.matrix_E[2][1], p->geometry.matrix_E[2][2],
      p->timestepvars.vmax, p->density.wcount_dh, p->density.wcount);
}

#endif /* SWIFT_GIZMO_MFM_HYDRO_DEBUG_H */
