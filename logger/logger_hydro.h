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
#ifndef SWIFT_LOGGER_HYDRO_H
#define SWIFT_LOGGER_HYDRO_H

/* Config parameters. */
#include "../config.h"

/* Import the right functions */
#if defined(MINIMAL_SPH)
#error TODO
#include "./hydro/Minimal/logger_hydro.h"
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/logger_hydro.h"
#elif defined(HOPKINS_PE_SPH)
#error TODO
#include "./hydro/PressureEntropy/logger_hydro.h"
#elif defined(HOPKINS_PU_SPH)
#error TODO
#include "./hydro/PressureEnergy/logger_hydro.h"
#elif defined(HOPKINS_PU_SPH_MONAGHAN)
#error TODO
#include "./hydro/PressureEnergyMorrisMonaghanAV/logger_hydro.h"
#elif defined(PHANTOM_SPH)
#error TODO
#include "./hydro/Phantom/logger_hydro.h"
#elif defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
#error TODO
#include "./hydro/Gizmo/logger_hydro.h"
#elif defined(SHADOWFAX_SPH)
#error TODO
#include "./hydro/Shadowswift/logger_hydro.h"
#elif defined(PLANETARY_SPH)
#error TODO
#include "./hydro/Planetary/logger_hydro.h"
#elif defined(SPHENIX_SPH)
#include "./hydro/SPHENIX/logger_hydro.h"
#elif defined(ANARCHY_PU_SPH)
#error TODO
#include "./hydro/AnarchyPU/logger_hydro.h"
#else
#error "Invalid choice of SPH variant"
#endif

#endif /* SWIFT_HYDRO_H */
