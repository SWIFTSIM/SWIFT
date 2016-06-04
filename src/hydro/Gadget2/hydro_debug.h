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

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    struct part* p, struct xpart* xp) {
  printf(
      "x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e],v_full=[%.3e,%.3e,%.3e] \n a=[%.3e,%.3e,%.3e],\n "
      "h=%.3e, "
      "wcount=%d, wcount_dh=%.3e, m=%.3e, dh_drho=%.3e, rho=%.3e, P=%.3e, "
      "S=%.3e, "
      "dS/dt=%.3e, c=%.3e\n"
      "divV=%.3e, curlV=%.3e, rotV=[%.3e,%.3e,%.3e]  \n "
      "v_sig=%e dh/dt=%.3e t_begin=%d, t_end=%d\n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], xp->v_full[0],
      xp->v_full[1], xp->v_full[2], p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->h, (int)p->density.wcount, p->density.wcount_dh, p->mass, p->rho_dh,
      p->rho, p->force.pressure, p->entropy, p->entropy_dt, p->force.soundspeed,
      p->div_v, p->force.curl_v, p->density.rot_v[0], p->density.rot_v[1],
      p->density.rot_v[2], p->force.v_sig, p->h_dt, p->ti_begin, p->ti_end);
}
