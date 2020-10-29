/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_HYDRO_LOGGER_H
#define SWIFT_HYDRO_LOGGER_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "align.h"
#include "logger.h"
#include "part_type.h"
#include "timeline.h"

/* Import the right function */
#if defined(MINIMAL_SPH)
#error TODO
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_logger.h"
#elif defined(HOPKINS_PE_SPH)
#error TODO
#elif defined(HOPKINS_PU_SPH)
#error TODO
#elif defined(HOPKINS_PU_SPH_MONAGHAN)
#error TODO
#elif defined(PHANTOM_SPH)
#error TODO
#elif defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
#error TODO
#elif defined(SHADOWFAX_SPH)
#error TODO
#elif defined(PLANETARY_SPH)
#error TODO
#elif defined(SPHENIX_SPH)
#error TODO
#elif defined(ANARCHY_PU_SPH)
#error TODO
#else
#error "Invalid choice of SPH variant"
#endif

#endif /* SWIFT_HYDRO_LOGGER_H */
