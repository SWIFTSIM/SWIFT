/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MULTIPOLE_ACCEPT_H
#define SWIFT_MULTIPOLE_ACCEPT_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "binomial.h"
#include "gravity_properties.h"
#include "integer_power.h"
#include "kernel_long_gravity.h"
#include "minmax.h"
#include "multipole_struct.h"

/**
 * @brief Compute the inverse of the force estimator entering the MAC
 *
 * Note that in the unsofted case, the first condition is naturally
 * never reached (as H == 0). In the non-periodic (non-truncated) case
 * the second condition is never reached (as r_s == inf, r_s_inv == 0).
 *
 * @param H The spline softening length.
 * @param r_s_inv The inverse of the scale of the gravity mesh.
 * @param r2 The square of the distance between the multipoles.
 */
__attribute__((const)) INLINE static float gravity_f_MAC_inverse(
    const float H, const float r_s_inv, const float r2) {

  if (r2 < (25.f / 81.f) * H * H) {

    /* Below softening radius */
    return (25.f / 81.f) * H * H;

  } else if (r_s_inv * r_s_inv * r2 > (25.f / 9.f)) {

    /* Above truncation radius */
    return (9.f / 25.f) * r_s_inv * r_s_inv * r2 * r2;

  } else {

    /* Normal Newtonian case */
    return r2;
  }
}

/**
 * @brief Checks whether The multipole in B can be used to update the field
 * tensor in A.
 *
 * We use the MAC of Dehnen 2014 eq. 16.
 *
 * Note: this is *not* symmetric in A<->B unless the purely geometric criterion
 * is used.
 *
 * @param props The properties of the gravity scheme.
 * @param A The gravity tensors that we want to update (sink).
 * @param B The gravity tensors that act as a source.
 * @param r2 The square of the distance between the centres of mass of A and B.
 * @param use_rebuild_sizes Are we considering the sizes at the last tree-build
 * (1) or current sizes (0)?
 * @param periodic Are we using periodic BCs?
 */
__attribute__((nonnull, pure)) INLINE static int gravity_M2L_accept(
    const struct gravity_props *props, const struct gravity_tensors *restrict A,
    const struct gravity_tensors *restrict B, const float r2,
    const int use_rebuild_sizes, const int periodic) {

  /* Order of the expansion */
  const int p = SELF_GRAVITY_MULTIPOLE_ORDER;

  /* Sizes of the multipoles */
  const float rho_A = use_rebuild_sizes ? A->r_max_rebuild : A->r_max;
  const float rho_B = use_rebuild_sizes ? B->r_max_rebuild : B->r_max;

  /* Get the softening */
  const float max_softening =
      max(A->m_pole.max_softening, B->m_pole.max_softening);

  /* Compute the error estimator (without the 1/M_B term that cancels out) */
  float E_BA_term = 0.f;
  for (int n = 0; n <= p; ++n) {
    E_BA_term +=
        binomial(p, n) * B->m_pole.power[n] * integer_powf(rho_A, p - n);
  }
  E_BA_term *= 8.f;
  if (rho_A + rho_B > 0.f) {
    E_BA_term *= max(rho_A, rho_B);
    E_BA_term /= (rho_A + rho_B);
  }

  /* Compute r^p */
#if SELF_GRAVITY_MULTIPOLE_ORDER % 2 == 1
  const float r_to_p = integer_powf(sqrtf(r2), p);
#else
  const float r_to_p = integer_powf(r2, (p / 2));
#endif

  float f_MAC_inv;
  if (props->consider_truncation_in_MAC) {
    f_MAC_inv = gravity_f_MAC_inverse(max_softening, props->r_s_inv, r2);
  } else {
    f_MAC_inv = r2;
  }

  /* Get the mimimal acceleration in A */
  const float min_a_grav = A->m_pole.min_old_a_grav_norm;

  /* Get the relative tolerance */
  const float eps = props->adaptive_tolerance;

  /* Get the basic geometric critical angle */
  const float theta_crit = props->theta_crit;
  const float theta_crit2 = theta_crit * theta_crit;

  /* Get the sum of the multipole sizes */
  const float rho_sum = rho_A + rho_B;

  if (props->use_advanced_MAC) {

#ifdef SWIFT_DEBUG_CHECKS
    if (min_a_grav == 0.) error("Acceleration is 0");
#endif

    /* Test the different conditions */

    /* Condition 1: We are in the converging part of the Taylor expansion */
    const int cond_1 = rho_sum * rho_sum < r2;

    /* Condition 2: We are not below softening */
    const int cond_2 =
        props->use_tree_below_softening || max_softening * max_softening < r2;

    /* Condition 3: The contribution is accurate enough
     * (E_BA * (1 / r^(p)) * ((1 / r^2) * W) < eps * a_min) */
    const int cond_3 = E_BA_term < eps * min_a_grav * r_to_p * f_MAC_inv;

    return cond_1 && cond_2 && cond_3;

  } else {

    /* Condition 1: We are obeying the purely geometric criterion */
    const int cond_1 = rho_sum * rho_sum < theta_crit2 * r2;

    /* Condition 2: We are not below softening */
    const int cond_2 =
        props->use_tree_below_softening || max_softening * max_softening < r2;

    return cond_1 && cond_2;
  }
}

/**
 * @brief Checks whether The multipole in B can be used to update the field
 * tensor in A and whether the multipole in A can be used to update the field
 * tensor in B.
 *
 * We use the MAC of Dehnen 2014 eq. 16.
 *
 * @param props The properties of the gravity scheme.
 * @param A The first set of multipole and gravity tensors.
 * @param B The second set of multipole and gravity tensors.
 * @param r2 The square of the distance between the centres of mass of A and B.
 * @param use_rebuild_sizes Are we considering the sizes at the last tree-build
 * (1) or current sizes (0)?
 * @param periodic Are we using periodic BCs?
 */
__attribute__((nonnull, pure)) INLINE static int gravity_M2L_accept_symmetric(
    const struct gravity_props *props, const struct gravity_tensors *restrict A,
    const struct gravity_tensors *restrict B, const float r2,
    const int use_rebuild_sizes, const int periodic) {

  return gravity_M2L_accept(props, A, B, r2, use_rebuild_sizes, periodic) &&
         gravity_M2L_accept(props, B, A, r2, use_rebuild_sizes, periodic);
}

/**
 * @brief Checks whether The multipole in B can be used to update the particle
 * pa
 *
 * We use the MAC of Dehnen 2014 eq. 16.
 *
 * @param props The properties of the gravity scheme.
 * @param pa The particle we want to compute forces for (sink)
 * @param B The gravity tensors that act as a source.
 * @param r2 The square of the distance between pa and the centres of mass of B.
 * @param periodic Are we using periodic BCs?
 */
__attribute__((nonnull, pure)) INLINE static int gravity_M2P_accept(
    const struct gravity_props *props, const struct gpart *pa,
    const struct gravity_tensors *B, const float r2, const int periodic) {

  /* Order of the expansion */
  const int p = SELF_GRAVITY_MULTIPOLE_ORDER;

  /* Sizes of the multipoles */
  const float rho_B = B->r_max;

  /* Get the maximal softening */
  const float max_softening =
      max(B->m_pole.max_softening, gravity_get_softening(pa, props));

#ifdef SWIFT_DEBUG_CHECKS
  if (rho_B == 0.) error("Size of multipole B is 0!");
#endif

  /* Compute the error estimator (without the 1/M_B term that cancels out) */
  const float E_BA_term = 8.f * B->m_pole.power[p];

  /* Compute r^p */
#if SELF_GRAVITY_MULTIPOLE_ORDER % 2 == 1
  const float r_to_p = integer_powf(sqrtf(r2), p);
#else
  const float r_to_p = integer_powf(r2, (p / 2));
#endif

  float f_MAC_inv;
  if (props->consider_truncation_in_MAC) {
    f_MAC_inv = gravity_f_MAC_inverse(max_softening, props->r_s_inv, r2);
  } else {
    f_MAC_inv = r2;
  }

  /* Get the estimate of the acceleration */
  const float old_a_grav = pa->old_a_grav_norm;

  /* Get the relative tolerance */
  const float eps = props->adaptive_tolerance;

  /* Get the basic geometric critical angle */
  const float theta_crit = props->theta_crit;
  const float theta_crit2 = theta_crit * theta_crit;

  if (props->use_advanced_MAC) {

#ifdef SWIFT_DEBUG_CHECKS
    if (old_a_grav == 0.) error("Acceleration is 0");
#endif

    /* Test the different conditions */

    /* Condition 1: We are in the converging part of the Taylor expansion */
    const int cond_1 = rho_B * rho_B < r2;

    /* Condition 2: We are not below softening */
    const int cond_2 =
        props->use_tree_below_softening || max_softening * max_softening < r2;

    /* Condition 3: The contribution is accurate enough
     * (E_BA * (1 / r^(p)) * ((1 / r^2) * W) < eps * a_min) */
    const int cond_3 = E_BA_term < eps * old_a_grav * r_to_p * f_MAC_inv;

    return cond_1 && cond_2 && cond_3;

  } else {

    /* Condition 1: We are obeying the purely geometric criterion */
    const int cond_1 = rho_B * rho_B < theta_crit2 * r2;

    /* Condition 2: We are not below softening */
    const int cond_2 =
        props->use_tree_below_softening || max_softening * max_softening < r2;

    return cond_1 && cond_2;
  }
}

#endif /* SWIFT_MULTIPOLE_ACCEPT_H */
