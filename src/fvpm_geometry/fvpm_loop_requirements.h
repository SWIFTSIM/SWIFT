/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_FVPM_LOOP_REQUIREMENTS_H
#define SWIFT_FVPM_LOOP_REQUIREMENTS_H

/**
 * @file fvpm_loop_requirements.h
 * @brief What the users of the FVPM geometry need from the neighbour loops.
 *
 * The hydro schemes declare their own needs in part.h and the RT schemes in
 * rt.h. This header only depends on config.h, so that every file that
 * includes part.h sees the same MPI switches.
 *
 * chemistry_flux_updates_inactive_neighbours: does the chemistry flux
 * exchange of the force and the chemistry flux prep loops also write an
 * inactive pj? The FVPM chemistry updates both particles of a pair, so the
 * rank of an inactive particle must run its pairs with foreign active
 * particles too. The hydro calls of these loops are not changed by it.
 */

/* Config parameters. */
#include <config.h>

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
#define chemistry_flux_updates_inactive_neighbours 1

/* Only over MPI: the switch also changes the unskip and the task graph of a
   single rank. */
#if defined(WITH_MPI) && !defined(MPI_SYMMETRIC_FORCE_INTERACTION)
#define MPI_SYMMETRIC_FORCE_INTERACTION
#endif

#else
#define chemistry_flux_updates_inactive_neighbours 0
#endif

#endif /* SWIFT_FVPM_LOOP_REQUIREMENTS_H */
