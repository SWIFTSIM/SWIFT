/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ADAPTIVE_SOFTENING_H
#define SWIFT_ADAPTIVE_SOFTENING_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "error.h"
#include "gravity_properties.h"
#include "inline.h"
#include "kernel_hydro.h"

#ifdef ADAPTIVE_SOFTENING

/* Force the gravity tasks to take place after the density loop */
#define gravity_after_hydro_density 1

/* Verify the correct hydro kernel is being used. */
#ifndef WENDLAND_C2_KERNEL
#error \
    "The adaptive softening terms are only defined for the Wendland-C2 kernel used in the gravity scheme."
#endif

/* Verify we are using the allowed hydro schemes */
#if defined(GIZMO_MFM_SPH) || defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
#error "Adaptive softening only implemented for the SPH schemes."
#endif

/* Verify we are using the correct gravity scheme */
#ifndef MULTI_SOFTENING_GRAVITY
#error \
    "Adaptive softening only implemented for the Multi-softening gravity scheme."
#endif

/**
 * @ifdef Update the gravity softening based on the gas' smoothing length.
 *
 * @param gp the #gpart to update.
 * @param p The #part.
 * @param grav_props The properties of the gravity scheme.
 */
INLINE static void gravity_update_softening(
    struct gpart* gp, const struct part* p,
    const struct gravity_props* grav_props) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gp != p->gpart) error("Unlinked part and gpart!");
#endif

  const float new_softening = p->h * kernel_gamma;

  /* Update the softening but respect limits */
  if (new_softening > grav_props->max_adaptive_softening)
    gp->epsilon = grav_props->max_adaptive_softening;
  else if (new_softening < grav_props->min_adaptive_softening)
    gp->epsilon = grav_props->min_adaptive_softening;
  else
    gp->epsilon = new_softening;
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * @param p The particle to act upon
 */
INLINE static void adaptive_softening_init_part(struct part* p) {

  p->adaptive_softening_data.zeta = 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Finish the equation 26 of Price & Monaghan 2007, MNRAS, 374, 4
 * by adding the pre-factors. Also adds the G/2 factor of eq. 27.
 *
 * @param p The particle to act upon.
 * @param props The properties of the gravity scheme.
 */
INLINE static void adaptive_softening_end_density(
    struct part* p, const struct gravity_props* props) {

  /* Finish calculation of the adaptive softening prefactor (zeta)
   * Pre-factor in eq. 26 of Price & Monaghan 2007 */
  const float rho_inv = 1.f / p->rho;
  const float h_drho = -p->h * rho_inv * hydro_dimension_inv;
  p->adaptive_softening_data.zeta *= h_drho;

  /* Multiply in the Gravitational constant
   * See the prefactor in eq. 27 of Price & Monaghan 2007 */
  p->adaptive_softening_data.zeta *= 0.5f * props->G_Newton;
}

#else

/* Gravity tasks can take place at the same time as hydro */
#define gravity_after_hydro_density 0

/**
 * @ifdef Update the gravity softening based on the gas' smoothing length.
 *
 * Softening is fixed: Nothing to do here.
 *
 * @param gp the #gpart to update.
 * @param p The #part.
 * @param grav_props The properties of the gravity scheme.
 */
INLINE static void gravity_update_softening(
    struct gpart* gp, const struct part* p,
    const struct gravity_props* grav_props) {}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 */
INLINE static void adaptive_softening_init_part(struct part* p) {}

/**
 * @brief Finishes the density calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 * @param props The properties of the gravity scheme.
 */
INLINE static void adaptive_softening_end_density(
    struct part* p, const struct gravity_props* props) {}

#endif /* ADAPTIVE_SOFTENING */

#endif /* SWIFT_ADAPTIVE_SOFTENING_H */
