/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STRUCT_SPHM1RT_H
#define SWIFT_RT_STRUCT_SPHM1RT_H

/**
 * @file src/rt/SPHM1RT/rt_struct.h
 * @brief Main header file for no radiative transfer struct.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/* Additional RT data in hydro particle struct */
struct rt_part_data {};

/* Additional RT data in star particle struct */
struct rt_spart_data {};

#endif /* SWIFT_RT_STRUCT_SPHM1RT_H */
