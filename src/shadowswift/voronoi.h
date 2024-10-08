/*******************************************************************************
* This file is part of SWIFT.
* Copyright (c) 2024 Matthieu Schaller (schaller@strw.leidenuniv.nl)
*                             Yolan Uyttenhove (Yolan.Uyttenhove@UGent.be)
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

#ifndef SWIFTSIM_SHADOWSWIFT_VORONOI_H
#define SWIFTSIM_SHADOWSWIFT_VORONOI_H

struct voronoi {
  int pair_count[27];
};

struct voronoi_pair {};

__attribute__((always_inline)) INLINE static void voronoi_destroy(struct voronoi* v) {}

#endif  // SWIFTSIM_SHADOWSWIFT_VORONOI_H
