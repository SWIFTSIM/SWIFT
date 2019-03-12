/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (ptcedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_CONST_H
#define SWIFT_CONST_H

/* SPH Viscosity constants. */
/* Cosmology default beta=3.0. Planetary default beta=4.0
 * Alpha can be set in the parameter file.
 * Beta is defined as in e.g. Price (2010) Eqn (103) */
#define const_viscosity_beta 3.0f

/* SPH Thermal conductivity constants. */
#define const_conductivity_alpha \
  1.f /* Value taken from (Price,2008), not used in legacy gadget mode */

/* Time integration constants. */
#define const_max_u_change 0.1f

/* Time-step limiter maximal difference in signal velocity */
#define const_limiter_max_v_sig_ratio 4.1f

/* Type of gradients to use (GIZMO_SPH only) */
/* If no option is chosen, no gradients are used (first order scheme) */
//#define GRADIENTS_SPH
#define GRADIENTS_GIZMO

/* Types of slope limiter to use (GIZMO_SPH only) */
/* Different slope limiters can be combined */
#define SLOPE_LIMITER_PER_FACE
#define SLOPE_LIMITER_CELL_WIDE

/* Types of flux limiter to use (GIZMO_SPH only) */
#define GIZMO_FLUX_LIMITER

/* Options to control the movement of particles for GIZMO_SPH. */
/* This option disables particle movement */
//#define GIZMO_FIX_PARTICLES
/* Try to keep cells regular by adding a correction velocity. */
//#define GIZMO_STEER_MOTION
/* Use the total energy instead of the thermal energy as conserved variable. */
//#define GIZMO_TOTAL_ENERGY

/* Options to control handling of unphysical values (GIZMO_SPH only). */
/* In GIZMO, mass and energy (and hence density and pressure) can in principle
   become negative, which will cause unwanted behaviour that can make the code
   crash.
   If no options are selected below, we assume (and pray) that this will not
   happen, and add no restrictions to how these variables are treated. */
/* Check for unphysical values and crash if they occur. */
//#define GIZMO_UNPHYSICAL_ERROR
/* Check for unphysical values and reset them to safe values. */
#define GIZMO_UNPHYSICAL_RESCUE
/* Show a warning message if an unphysical value was reset (only works if
   GIZMO_UNPHYSICAL_RESCUE is also selected). */
#ifdef SWIFT_DEBUG_CHECKS
#define GIZMO_UNPHYSICAL_WARNING
#endif

/* Parameters that control how GIZMO handles pathological particle
   configurations. */
/* Show a warning message if a pathological configuration has been detected. */
//#define GIZMO_PATHOLOGICAL_WARNING
/* Crash if a pathological configuration has been detected. */
//#define GIZMO_PATHOLOGICAL_ERROR
/* Maximum allowed gradient matrix condition number. If the condition number of
   the gradient matrix (defined in equation C1 in Hopkins, 2015) is larger than
   this value, we artificially increase the number of neighbours to get a more
   homogeneous sampling. */
#define const_gizmo_max_condition_number 100.0f
/* Correction factor applied to the particle wcount to force more neighbours if
   the condition number is too large. */
#define const_gizmo_w_correction_factor 0.9f
/* Lower limit on the wcount correction factor. If the condition number is still
   too high after this wcount correction has been applied, we give up on the
   gradient matrix and use SPH gradients instead. */
#define const_gizmo_min_wcorr 0.5f

/* Types of gradients to use for SHADOWFAX_SPH */
/* If no option is chosen, no gradients are used (first order scheme) */
#define SHADOWFAX_GRADIENTS

/* SHADOWFAX_SPH slope limiters */
#define SHADOWFAX_SLOPE_LIMITER_PER_FACE
#define SHADOWFAX_SLOPE_LIMITER_CELL_WIDE

/* Options to control SHADOWFAX_SPH */
/* This option disables cell movement */
//#define SHADOWFAX_FIX_CELLS
/* This option enables cell steering, i.e. trying to keep the cells regular by
   adding a correction to the cell velocities.*/
#define SHADOWFAX_STEER_CELL_MOTION
/* This option evolves the total energy instead of the thermal energy */
//#define SHADOWFAX_TOTAL_ENERGY

/* Source terms */
#define SOURCETERMS_NONE
//#define SOURCETERMS_SN_FEEDBACK

#endif /* SWIFT_CONST_H */
