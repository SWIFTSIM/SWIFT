/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_HYDRO_H
#define SWIFT_HYDRO_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "hydro_properties.h"
#include "kernel_hydro.h"
#include "part.h"

/* Import the right functions */
#if defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro.h"
#include "./hydro/Minimal/hydro_iact.h"
#define SPH_IMPLEMENTATION "Minimal version of SPH (e.g. Price 2010)"
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro.h"
#include "./hydro/Gadget2/hydro_iact.h"
#define SPH_IMPLEMENTATION "Gadget-2 version of SPH (Springel 2005)"
#elif defined(HOPKINS_PE_SPH)
#include "./hydro/PressureEntropy/hydro.h"
#include "./hydro/PressureEntropy/hydro_iact.h"
#define SPH_IMPLEMENTATION "Pressure-Entropy SPH (Hopkins 2013)"
#elif defined(HOPKINS_PU_SPH)
#include "./hydro/PressureEnergy/hydro.h"
#include "./hydro/PressureEnergy/hydro_iact.h"
#define SPH_IMPLEMENTATION "Pressure-Energy SPH (Hopkins 2013)"
#elif defined(HOPKINS_PU_SPH_MONAGHAN)
#include "./hydro/PressureEnergyMorrisMonaghanAV/hydro.h"
#include "./hydro/PressureEnergyMorrisMonaghanAV/hydro_iact.h"
#define SPH_IMPLEMENTATION                                                \
  "Pressure-Energy SPH (Hopkins 2013) with a Morris and Monaghan (1997) " \
  "variable artificial viscosity."
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro.h"
#include "./hydro/Default/hydro_iact.h"
#define SPH_IMPLEMENTATION "Default version of SPH"
#elif defined(GIZMO_MFV_SPH)
#include "./hydro/GizmoMFV/hydro.h"
#include "./hydro/GizmoMFV/hydro_iact.h"
#define SPH_IMPLEMENTATION "GIZMO MFV (Hopkins 2015)"
#elif defined(GIZMO_MFM_SPH)
#include "./hydro/GizmoMFM/hydro.h"
#include "./hydro/GizmoMFM/hydro_iact.h"
#define SPH_IMPLEMENTATION "GIZMO MFM (Hopkins 2015)"
#elif defined(SHADOWFAX_SPH)
#include "./hydro/Shadowswift/hydro.h"
#include "./hydro/Shadowswift/hydro_iact.h"
#define SPH_IMPLEMENTATION \
  "Shadowfax moving mesh (Vandenbroucke and De Rijcke 2016)"
#elif defined(PLANETARY_SPH)
#include "./hydro/Planetary/hydro.h"
#include "./hydro/Planetary/hydro_iact.h"
#define SPH_IMPLEMENTATION "Minimal version of SPH with multiple materials"
#elif defined(ANARCHY_DU_SPH)
#include "./hydro/AnarchyDU/hydro.h"
#include "./hydro/AnarchyDU/hydro_iact.h"
#define SPH_IMPLEMENTATION "ANARCHY (Density-Energy) SPH (Borrow+ in prep)"
#elif defined(ANARCHY_PU_SPH)
#include "./hydro/AnarchyPU/hydro.h"
#include "./hydro/AnarchyPU/hydro_iact.h"
#define SPH_IMPLEMENTATION \
  "ANARCHY (Pressure-Energy) SPH (Dalla Vecchia+ in prep)"
#else
#error "Invalid choice of SPH variant"
#endif

#endif /* SWIFT_HYDRO_H */
