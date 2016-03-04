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

#include  "physical_constants.h"

/* Hydrodynamical constants. */
#define const_hydro_gamma (5.0f / 3.0f)

/* SPH Viscosity constants. */
#define const_viscosity_alpha \
  0.8f /* Used in the legacy gadget-2 SPH mode only */
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
#define const_cfl 0.1f
#define const_ln_max_h_change                                           \
  0.231111721f /* Particle can't change volume by more than a factor of \
                  2=1.26^3 over one time step */
#define const_max_u_change 0.1f

/* Neighbour search constants. */
#define const_eta_kernel \
  1.2349f /* Corresponds to 48 ngbs with the cubic spline kernel */
#define const_delta_nwneigh 0.1f
#define const_smoothing_max_iter 30
#define CUBIC_SPLINE_KERNEL

/* Gravity stuff. */
#define const_theta_max                                   \
  0.57735f /* Opening criteria, which is the ratio of the \
              cell distance over the cell width. */

#define const_epsilon 0.0014f         /* Gravity blending distance. */
#define const_iepsilon 714.285714286f /* Inverse gravity blending distance. */
#define const_iepsilon2 (const_iepsilon* const_iepsilon)
#define const_iepsilon3 (const_iepsilon2* const_iepsilon)
#define const_iepsilon4 (const_iepsilon2* const_iepsilon2)
#define const_iepsilon5 (const_iepsilon3* const_iepsilon2)
#define const_iepsilon6 (const_iepsilon3* const_iepsilon3)

/* SPH variant to use */
//#define MINIMAL_SPH
//#define GADGET2_SPH
//#define DEFAULT_SPH
#define NO_SPH

/* Gravity properties */
#define GRAVITY
/* valid choices DEFAULT_GRAVITY || EXTERNAL_POTENTIAL */
#define EXTERNAL_POTENTIAL
//#define DEFAULT_GRAVITY


/* System of units */
#define const_unit_length_in_cgs       (1000 *  PARSEC_IN_CGS )  /* kpc */
#define const_unit_mass_in_cgs         (SOLAR_MASS_IN_CGS)       /* solar mass */
#define const_unit_velocity_in_cgs     (1e5)                     /* km/s */ 

/* Derived constants */
#define const_unit_time_in_cgs (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
#define const_G (NEWTON_GRAVITY_CGS * const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/ (const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs))

/* External Potential Constants */
#define External_Potential_X  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Y  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Z  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Mass (1e10 * SOLAR_MASS_IN_CGS / const_unit_mass_in_cgs)

#endif /* SWIFT_CONST_H */
