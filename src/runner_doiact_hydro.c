/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "chemistry.h"
#include "engine.h"
#include "mhd.h"
#include "pressure_floor_iact.h"
#include "rt.h"
#include "runner.h"
#include "runner_doiact_hydro_vec.h"
#include "runner_doiact_sinks.h"
#include "sink_iact.h"
#include "sink_properties.h"
#include "space_getsid.h"
#include "star_formation_iact.h"
#include "timers.h"
#include "timestep_limiter_iact.h"

/* Import the density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_functions_hydro.h"
#include "runner_doiact_undef.h"

/* Import the gradient loop functions (if required). */
#ifdef EXTRA_HYDRO_LOOP
#define FUNCTION gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_GRADIENT
#include "runner_doiact_functions_hydro.h"
#include "runner_doiact_undef.h"
#endif

/* Import the force loop functions. */
#define FUNCTION force
#define FUNCTION_TASK_LOOP TASK_LOOP_FORCE
#include "runner_doiact_functions_hydro.h"
#include "runner_doiact_undef.h"

/* Import the RT gradient loop functions */
#define FUNCTION rt_gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_GRADIENT
#include "runner_doiact_functions_hydro.h"
#include "runner_doiact_undef.h"

/* Import the RT transport (force) loop functions. */
#define FUNCTION rt_transport
#define FUNCTION_TASK_LOOP TASK_LOOP_RT_TRANSPORT
#include "runner_doiact_functions_hydro.h"
#include "runner_doiact_undef.h"

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

#if defined(WITH_VECTORIZATION) && defined(GADGET2_SPH)
#error \
    "The GEAR FVPM chemistry FCT prep loop is not implemented for the vectorized GADGET2 force loop."
#endif

/* No-op hydro/MHD kernels for the chemistry FCT prep loop: the loop template
   is instantiated with the force-loop code paths (so the prep pass visits the
   exact same pair evaluations as the flux exchange), but only the chemistry
   prep interaction must do any work. */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_fct_prep(const float r2, const float dx[3],
                               const float hi, const float hj,
                               struct part *restrict pi,
                               struct part *restrict pj, const float a,
                               const float H) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_chemistry_fct_prep(const float r2, const float dx[3],
                                      const float hi, const float hj,
                                      struct part *restrict pi,
                                      struct part *restrict pj, const float a,
                                      const float H) {}

__attribute__((always_inline)) INLINE static void
runner_iact_mhd_chemistry_fct_prep(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   struct part *restrict pj, const double mu_0,
                                   const float a, const float H) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mhd_chemistry_fct_prep(const float r2, const float dx[3],
                                          const float hi, const float hj,
                                          struct part *restrict pi,
                                          struct part *restrict pj,
                                          const double mu_0, const float a,
                                          const float H) {}

/* Route the force-loop chemistry diffusion calls of the template to the FCT
   prep pass (accumulate donor outflow sums, apply nothing). The timebin
   updates also run here; they are idempotent min() updates, so repeating them
   in the force loop is harmless. runner_iact_(nonsym_)diffusion are empty for
   this scheme (see chemistry_iact.h) so this redirect is a no-op in itself;
   it is kept so any private out-of-tree scheme sharing this template still
   gets its own diffusion hooks routed to the prep pass unchanged. The new
   flux-exchange extension gets the identical redirect, alongside. */
#define runner_iact_diffusion runner_iact_diffusion_fct_prep
#define runner_iact_nonsym_diffusion runner_iact_nonsym_diffusion_fct_prep
#define runner_iact_chemistry_flux_exchange \
  runner_iact_chemistry_flux_exchange_fct_prep

/* Import the chemistry FCT prep loop functions. */
#define FUNCTION chemistry_fct_prep
#define FUNCTION_TASK_LOOP TASK_LOOP_FORCE
#include "runner_doiact_functions_hydro.h"
#include "runner_doiact_undef.h"

#undef runner_iact_diffusion
#undef runner_iact_nonsym_diffusion
#undef runner_iact_chemistry_flux_exchange

#endif /* GEAR FVPM diffusion chemistry */
