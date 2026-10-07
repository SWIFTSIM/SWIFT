/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_STRUCT_TREECOOL_H
#define SWIFT_COOLING_STRUCT_TREECOOL_H

/**
 * @file src/cooling/TREECOOL/cooling_struct.h
 * @brief Particle data related to the TREECOOL cooling function (Katz,
 * Weinberg & Hernquist 1996).
 */

/**
 * @brief Properties of the cooling stored in the #part data.
 */
struct cooling_part_data {};

/**
 * @brief Properties of the cooling stored in the #xpart data.
 */
struct cooling_xpart_data {

  /*! Net energy radiated away by this particle since the start of the run
   * (negative if photo-heating dominated) */
  float radiated_energy;

  /*! Electron number density in units of the Hydrogen number density.
   *
   * This is the ionization equilibrium at the energy reached at the end of the
   * last cooling step of this particle. It is used as the starting guess of
   * the iterative solvers, written to the snapshots and returned by
   * cooling_get_electron_density(). It carries no information that the model
   * could not recompute from scratch. */
  float electron_fraction;
};

#endif /* SWIFT_COOLING_STRUCT_TREECOOL_H */
