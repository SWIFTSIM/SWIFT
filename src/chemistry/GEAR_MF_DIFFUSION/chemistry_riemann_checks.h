/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_CHECKS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_CHECKS_H

#include "error.h"

/* Import the right chemistry_part_data definition */
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
#include "hyperbolic/chemistry_riemann_checks.h"
#else
#include "parabolic/chemistry_riemann_checks.h"
#endif

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_CHECKS_H */
