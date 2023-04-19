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

/* A large number of macros are defined in runner_doiact_hydro.h
 * based on the values of some macros defined in runner_main.c
 * For each loop definition (density, force, etc.) the values are different
 * so we need to undefine everything before redefining them.
 * This is done here in a single place to avoid code duplication and missing
 * undefs if they were scattered all over the place */
#undef IACT
#undef IACT_NONSYM
#undef IACT_MHD
#undef IACT_NONSYM_MHD
#undef IACT_STARS
#undef IACT_BH_GAS
#undef IACT_BH_BH
#undef GET_MU0
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* these are defined in runner_doiact_functions_hydro.h at every #include */
#undef PART_IS_ACTIVE
#undef CELL_IS_ACTIVE
#undef CELL_ARE_PART_DRIFTED
#undef DO_DRIFT_DEBUG_CHECKS
