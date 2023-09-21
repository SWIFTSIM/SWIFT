/*******************************************************************************
* This file is part of SWIFT.
* Copyright (c) 2021 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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
#ifndef SWIFT_CHEMISTRY_ADDITIONS_H
#define SWIFT_CHEMISTRY_ADDITIONS_H

/**
 * @file Contains some forward declarations of chemistry functions used outside
 * of the chemistry code files, to avoid circular inclusions/implicit
 * declarations.
 **/

__attribute__((always_inline)) INLINE static void runner_iact_chemistry_fluxes(
    struct part *restrict pi, struct part *restrict pj, float mass_flux,
    float flux_dt, int mode);

#endif  // SWIFT_CHEMISTRY_ADDITIONS_H
