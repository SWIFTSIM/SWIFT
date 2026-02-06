/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MHD_H
#define SWIFT_MHD_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "part.h"

/* Import the right functions */
#if defined(NONE_MHD)
#include "./mhd/None/mhd.h"
#include "./mhd/None/mhd_iact.h"
#define MHD_IMPLEMENTATION "No MHD scheme"
#elif defined(DIRECT_INDUCTION_MHD)
#include "./mhd/DirectInduction/mhd.h"
#define EXTRA_HYDRO_LOOP
#include "./mhd/DirectInduction/mhd_iact.h"
#define MHD_IMPLEMENTATION "MHD scheme using direct induction"
#elif defined(DIRECT_INDUCTION_FEDE_MHD)
#include "./mhd/DInduction/mhd.h"
#include "./mhd/DInduction/mhd_iact.h"
#define MHD_IMPLEMENTATION "MHD scheme using direct induction BASE"
#define EXTRA_HYDRO_LOOP
#elif defined(VECTOR_POTENTIAL_MHD)
#include "./mhd/VPotential/mhd.h"
#include "./mhd/VPotential/mhd_iact.h"
#define MHD_IMPLEMENTATION "Vector potentials"
#define EXTRA_HYDRO_LOOP
#else
#error "Invalid choice of MHD variant"
#endif

#endif /* SWIFT_MHD_H */
