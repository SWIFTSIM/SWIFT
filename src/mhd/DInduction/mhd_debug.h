/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DI_MHD_DEBUG_H
#define SWIFT_DI_MHD_DEBUG_H

/**
 * @brief Print out the mhd fields of a particle.
 *
 * Function used for debugging purposes.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_debug_particle(
    const struct part *p, const struct xpart *xp) {
  printf(
      "Bfld=[%.3e,%.3e,%.3e], "
      "Bpred=[%.3e,%.3e,%.3e], "
      "dBdt=[%.3e,%.3e,%.3e], \n"
      "divB=%.3e, Q1=%.3e, Q0=%.3e, phi=%.3e\n",
      xp->mhd_data.Bfld_full[0], xp->mhd_data.Bfld_full[1], xp->mhd_data.Bfld_full[2],
      p->mhd_data.BPred[0], p->mhd_data.BPred[1], p->mhd_data.BPred[2],
      p->mhd_data.dBdt[0], p->mhd_data.dBdt[1], p->mhd_data.dBdt[2],
      p->mhd_data.divB, p->mhd_data.Q1, p->mhd_data.Q0, p->mhd_data.phi);
}
#endif /* SWIFT_DI_MHD_DEBUG_H */
