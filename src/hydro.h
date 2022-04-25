/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#if defined(NONE_SPH)
#include "./hydro/None/hydro.h"
#include "./hydro/None/hydro_iact.h"
#define SPH_IMPLEMENTATION "No hydro scheme"
#elif defined(MINIMAL_SPH)
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
#elif defined(PHANTOM_SPH)
#include "./hydro/Phantom/hydro.h"
#include "./hydro/Phantom/hydro_iact.h"
#define SPH_IMPLEMENTATION "PHANTOM SPH reference implementation (Price 2018)"
#elif defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
#include "./hydro/Gizmo/hydro.h"
#include "./hydro/Gizmo/hydro_iact.h"
#elif defined(SHADOWFAX_SPH)
#include "./hydro/Shadowswift/hydro.h"
#include "./hydro/Shadowswift/hydro_iact.h"
#define SPH_IMPLEMENTATION \
  "Shadowfax moving mesh (Vandenbroucke and De Rijcke 2016)"
#elif defined(PLANETARY_SPH)
#include "./hydro/Planetary/hydro.h"
#include "./hydro/Planetary/hydro_iact.h"
#define SPH_IMPLEMENTATION "Minimal version of SPH with multiple materials"
#elif defined(SPHENIX_SPH)
#include "./hydro/SPHENIX/hydro.h"
#include "./hydro/SPHENIX/hydro_iact.h"
#define SPH_IMPLEMENTATION "SPHENIX (Borrow+ 2020)"
#elif defined(GASOLINE_SPH)
#include "./hydro/Gasoline/hydro.h"
#include "./hydro/Gasoline/hydro_iact.h"
#define SPH_IMPLEMENTATION "Gasoline-2 (Wadsley+ 2017)"
#elif defined(ANARCHY_PU_SPH)
#include "./hydro/AnarchyPU/hydro.h"
#include "./hydro/AnarchyPU/hydro_iact.h"
#define SPH_IMPLEMENTATION \
  "ANARCHY (Pressure-Energy) SPH (Dalla Vecchia+ in prep)"
#else
#error "Invalid choice of SPH variant"
#endif

/* Check whether this scheme implements the density checks */
#ifdef SWIFT_HYDRO_DENSITY_CHECKS
#if !defined(SPHENIX_SPH) && !defined(PLANETARY_SPH)
#error \
    "Can only use the hydro brute-force density checks with the SPHENIX or PLANETARY hydro schemes."
#endif
#endif

struct engine;
struct space;

void hydro_exact_density_compute(struct space *s, const struct engine *e,
                                 const int check_force);
void hydro_exact_density_check(struct space *s, const struct engine *e,
                               const float rel_tol, const int check_force);

#endif /* SWIFT_HYDRO_H */
