/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (ptcedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* I/O Constant; this determines the relative tolerance between the value of
 * redshift read from the snapshot, and the value from the parameter file. This
 * current value asserts that they must match within 0.1%. */
#define io_redshift_tolerance 1e-3f

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

/* Options controlling ShadowSWIFT */
/* Option to enable gradients for ShadowSWIFT */
/* If disabled, no gradients are used (first order scheme) */
#define SHADOWSWIFT_MESHLESS_GRADIENTS
#ifndef SHADOWSWIFT_MESHLESS_GRADIENTS
#define SHADOWSWIFT_GRADIENTS
#endif
/* Always activate the slope limiters if we use gradients (the scheme becomes
 * unstable otherwise) */
#ifdef SHADOWSWIFT_GRADIENTS
/*! @brief Option controlling which type of gradient calculation is used */
#define SHADOWSWIFT_GRADIENTS_WLS
/* Always activate the slope limiters if we use gradients (the scheme becomes
 * unstable otherwise) */
#define SHADOWSWIFT_SLOPE_LIMITER_PER_FACE
#define SHADOWSWIFT_SLOPE_LIMITER_CELL_WIDE
#endif
#ifdef SHADOWSWIFT_MESHLESS_GRADIENTS
#define SHADOWSWIFT_SLOPE_LIMITER_PER_FACE
#define SHADOWSWIFT_SLOPE_LIMITER_MESHLESS
#endif
/* Option controlling output of grids */
/*! @brief Option to enable time extrapolation */
#define SHADOWSWIFT_EXTRAPOLATE_TIME

/*! @brief Option controlling output of grids */
//#define SHADOWSWIFT_OUTPUT_GRIDS

/* Options controlling acceleration strategies*/
/*! @brief Option to enable the hilbert order insertion during the grid construction */
#define SHADOWSWIFT_HILBERT_ORDERING
/*! @brief Option to enable the bvh acceleration structure for neighbour searching */
#define SHADOWSWIFT_BVH

/* Options controlling particle movement */
/*! @brief This option disables cell movement */
//#define SHADOWSWIFT_FIX_PARTICLES
/*! @brief This option enables cell steering, i.e. trying to keep the cells regular by
 * adding a correction to the cell velocities.*/
#ifndef SHADOWSWIFT_FIX_PARTICLES
#define SHADOWSWIFT_STEER_MOTION
#endif

/*! @brief This option enables boundary conditions for non-periodic ShadowSWIFT runs */
#define VACUUM_BC 0
#define REFLECTIVE_BC 1
#define OPEN_BC 2
#define INFLOW_BC 3
#define RADIAL_INFLOW_BC 4
#define SHADOWSWIFT_BC REFLECTIVE_BC

/* Options controlling behaviour of the code when unphysical situations are
 * encountered */
/*! @brief This option tries to recover from unphysical situations */
#define SHADOWSWIFT_UNPHYSICAL_RESCUE
#ifdef SHADOWSWIFT_UNPHYSICAL_RESCUE
/*! @brief Show a warning message if an unphysical value was reset */
#define SHADOWSWIFT_UNPHYSICAL_WARNING
#else
/*! @brief This option halts the execution in the case of unphysical conditions */
#define SHADOWSWIFT_UNPHYSICAL_ERROR
#endif

/* Source terms */
#define SOURCETERMS_NONE
//#define SOURCETERMS_SN_FEEDBACK

/* GRACKLE doesn't really like exact zeroes, so use something
 * comparatively small instead. */
#define RT_GEAR_TINY_MASS_FRACTION 1.e-20

#endif /* SWIFT_CONST_H */
