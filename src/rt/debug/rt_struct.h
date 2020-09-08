/*******************************************************************************
 * This file is part of SWIFT.
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
#ifndef SWIFT_RT_STRUCT_DEBUG_H
#define SWIFT_RT_STRUCT_DEBUG_H

/**
 * @file src/rt/debug/rt_struct.h
 * @brief Main header file for the debug radiative transfer struct.
 */

struct rt_xpart_data {
  int iact_stars; /* how many stars this particle interacted with */
  int calls_tot;  /* total number of calls to this particle during entire run */
  int calls_per_step; /* calls per time step to this particle */
  int calls_self;
  int calls_pair;
};

struct rt_spart_data {
  int iact_hydro; /* how many hydro particles this particle interacted with */
  int calls_tot;  /* total number of calls to this particle during entire run */
  int calls_per_step; /* calls per time step to this particle */
  int calls_self;
  int calls_pair;
};

#endif /* SWIFT_RT_STRUCT_DEBUG_H */
