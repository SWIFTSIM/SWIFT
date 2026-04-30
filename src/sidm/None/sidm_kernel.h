/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_NONE_SIDM_KERNEL_H
#define SWIFT_NONE_SIDM_KERNEL_H

/**
 * @brief Kernel overlap assuming top-hat kernel. 
 *
 * @param r Comoving distance between the two particles.
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 */
INLINE static float sidm_kernel_overlap_tophat(const float r, const float hi,
                                               const float hj) {}


#endif /* SWIFT_NONE_SIDM_KERNEL_H */
