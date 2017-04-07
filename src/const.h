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
#define const_viscosity_alpha 0.8f
#define const_viscosity_alpha_min \
  0.1f /* Values taken from (Price,2004), not used in legacy gadget mode */
#define const_viscosity_alpha_max \
  2.0f /* Values taken from (Price,2004), not used in legacy gadget mode */
#define const_viscosity_length \
  0.1f /* Values taken from (Price,2004), not used in legacy gadget mode */

/* SPH Thermal conductivity constants. */
#define const_conductivity_alpha \
  1.f /* Value taken from (Price,2008), not used in legacy gadget mode */

/* Time integration constants. */
#define const_max_u_change 0.1f

/* Thermal energy per unit mass used as a constant for the isothermal EoS */
#define const_isothermal_internal_energy 20.2615290634f

/* Type of gradients to use (GIZMO_SPH only) */
/* If no option is chosen, no gradients are used (first order scheme) */
//#define GRADIENTS_SPH
#define GRADIENTS_GIZMO

/* Types of slope limiter to use (GIZMO_SPH only) */
/* Different slope limiters can be combined */
#define SLOPE_LIMITER_PER_FACE
#define SLOPE_LIMITER_CELL_WIDE

/* Options to control the movement of particles for GIZMO_SPH. */
/* This option disables particle movement */
//#define GIZMO_FIX_PARTICLES
//#define GIZMO_TOTAL_ENERGY

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
