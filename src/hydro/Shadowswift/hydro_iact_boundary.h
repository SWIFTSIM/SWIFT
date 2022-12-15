/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_IACT_BOUNDARY_H
#define SWIFT_SHADOWSWIFT_HYDRO_IACT_BOUNDARY_H

#include "hydro_flux.h"
#include "hydro_getters.h"
#include "hydro_gradients.h"
#include "hydro_part.h"
#include "hydro_setters.h"

__attribute__((always_inline)) INLINE static void
runner_iact_boundary_set_primitives(struct part *p_boundary,
                                    const struct part *p,
                                    const struct hydro_space *hs) {
#if SHADOWSWIFT_BC == VACUUM_BC
  p_boundary->rho = 0.;
  p_boundary->v[0] = 0.;
  p_boundary->v[1] = 0.;
  p_boundary->v[2] = 0.;
  p_boundary->P = 0.;
#elif SHADOWSWIFT_BC == OPEN_BC
  /* Nothing to do here, primitives are already equal since p_boundary is
   * created from a copy of p */
#ifdef SWIFT_DEBUG_CHECKS
  assert(p_boundary->rho == p->rho);
  assert(p_boundary->v[0] == p->v[0]);
  assert(p_boundary->v[1] == p->v[1]);
  assert(p_boundary->v[2] == p->v[2]);
  assert(p_boundary->P == p->P);
#endif
#elif SHADOWSWIFT_BC == REFLECTIVE_BC
  /* We just need to reflect the velocity component perpendicular to the
   * boundary */
#ifdef SWIFT_DEBUG_CHECKS
  assert(p_boundary->rho == p->rho);
  assert(p_boundary->P == p->P);
#endif
  for (int i = 0; i < 3; i++) {
    if (p_boundary->x[i] != p->x[i]) p_boundary->v[i] = -p_boundary->v[i];
  }
#elif SHADOWSWIFT_BC == INFLOW_BC
  if (p->x[0] != p_boundary->x[0] &&
      hs->velocity * (p->x[0] - p_boundary->x[0]) > 0.) {
    /* Only enforce inflow boundary conditions on the side from which the gas
     * is flowing in. */
    p_boundary->rho = hs->density;
    p_boundary->v[0] = hs->velocity;
    p_boundary->v[1] = 0.f;
    p_boundary->v[2] = 0.f;
    p_boundary->P = hs->pressure;
  } else {
    /* Open boundary conditions */
#ifdef SWIFT_DEBUG_CHECKS
    assert(p_boundary->rho == p->rho);
    assert(p_boundary->v[0] == p->v[0]);
    assert(p_boundary->v[1] == p->v[1]);
    assert(p_boundary->v[2] == p->v[2]);
    assert(p_boundary->P == p->P);
#endif
  }
#elif SHADOWSWIFT_BC == RADIAL_INFLOW_BC
  /* Get the direction to the center of the simulation */
  double dx[3] = {
      hs->center[0] - p_boundary->x[0],
      hs->center[1] - p_boundary->x[1],
      hs->center[2] - p_boundary->x[2],
  };
  double norm = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  p_boundary->rho = hs->density;
  p_boundary->v[0] = hs->velocity * dx[0] / norm;
  p_boundary->v[1] = hs->velocity * dx[1] / norm;
  p_boundary->v[2] = hs->velocity * dx[2] / norm;
  p_boundary->P = hs->pressure;
#else
#error "Unknown boundary condition!"
#endif
}

/**
 * @brief Update the slope estimates of particles pi and pj.
 *
 * @param p The "real" particle. Must always be active.
 * @param p_boundary The imaginary boundary particle.
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Whether to do the slope estimate for both particles.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_boundary_slope_estimate(struct part *p, struct part *p_boundary,
                                    double const *centroid, float surface_area,
                                    const struct hydro_space *hs) {
  const double shift[3] = {0., 0., 0.};
  runner_iact_boundary_set_primitives(p_boundary, p, hs);
  runner_iact_slope_estimate(p, p_boundary, centroid, surface_area, shift, 0);
}

/**
 * @brief Collect info necessary for limiting the gradient estimates.
 *
 * @param p The "real" particle. Must always be active.
 * @param p_boundary The imaginary boundary particle.
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Whether to do the slope limiting for both particles.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_boundary_slope_limiter(struct part *p, struct part *p_boundary,
                                   const double *centroid, float surface_area,
                                   const struct hydro_space *hs) {

  const double shift[3] = {0., 0., 0.};
  runner_iact_boundary_set_primitives(p_boundary, p, hs);
  runner_iact_slope_limiter(p, p_boundary, centroid, surface_area, shift, 0);
}

__attribute__((always_inline)) INLINE static float fsgnf(float x) {
  if (x > 0.f)
    return 1.f;
  else if (x < 0.f)
    return -1.f;
  else
    return 0.f;
}

__attribute__((always_inline)) INLINE static float
riemann_solve_reflective_for_pressure(const float rho, const float v,
                                      const float P) {

  /* calculate sound speeds */
  float a = sqrtf(hydro_gamma * P / rho);

  /* Check for vacuum */
  if (!rho || (hydro_two_over_gamma_minus_one * a <= -v)) {
    float Whalf[5];
    float WL[5] = {rho, v, 0.f, 0.f, P};
    float WR[5] = {rho, -v, 0.f, 0.f, P};
    float n_unit[3] = {1.f, 0.f, 0.f};
    riemann_solve_vacuum(WL, WR, v, -v, a, a, Whalf, n_unit);
    return Whalf[4];
  }

  /* Solve eq. 4.5 from Toro for pressure */
  if (v > 0.f) {
    /* Shock waves */
    float A = hydro_two_over_gamma_plus_one / rho;
    float B = hydro_gamma_minus_one_over_gamma_plus_one * P;
    return (v * sqrtf(4.f * P * A + 4.f * A * B + v * v) + 2.f * P * A +
            v * v) /
           (2.f * A);
  } else {
    /* Rarefaction waves */
    float res = 0.5f * v * hydro_gamma_minus_one / a + 1.f;
    return pow_two_gamma_over_gamma_minus_one(res) * P;
  }
}

/**
 * @brief The flux calculation between particle i and j
 *
 * This method calculates the surface area of the interface between particle i
 * and particle j, as well as the interface position and velocity. These are
 * then used to reconstruct and predict the primitive variables, which are then
 * fed to a Riemann solver that calculates a flux. This flux is used to update
 * the conserved variables of both particles.
 *
 * This method also calculates the maximal velocity used to calculate the time
 * step.
 *
 * @param p The "real" particle. Must always be active.
 * @param p_boundary The imaginary boundary particle.
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Unused, flux exchange must be manifestly symmetric.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_boundary_flux_exchange(struct part *p, struct part *p_boundary,
                                   double const *centroid, float surface_area,
                                   const struct hydro_space *hs) {

#if SHADOWSWIFT_BC == VACUUM_BC || SHADOWSWIFT_BC == INFLOW_BC || \
    SHADOWSWIFT_BC == RADIAL_INFLOW_BC
  const double shift[3] = {0., 0., 0.};
  runner_iact_boundary_set_primitives(p_boundary, p, hs);
  /* Set gradients of boundary particle to 0. */
  hydro_gradients_init(p_boundary);
  runner_iact_flux_exchange(p, p_boundary, centroid, surface_area, shift, 0);
#else
  /* Get flux vector */
  /* Vector from pj to pi */
  float dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = (float)(p->x[k] - p_boundary->x[k]);
  }
  /* Normal vector at interface (pointing to pj) */
  float n_unit[3] = {fsgnf(-dx[0]), fsgnf(-dx[1]), fsgnf(-dx[2])};
  /* Get primitive variables of pi. */
  float WL[6];
  hydro_part_get_primitive_variables(p, WL);
  /* Calculate flux: */
  float totflux[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
#if SHADOWSWIFT_BC == OPEN_BC
  /* Open BC: Whalf = (rho_l, v_l, P_l), interface velocity = 0.0, so we can
   * calculate the flux exactly. */
  riemann_flux_from_half_state(WL, &WL[1], n_unit, totflux);
#elif SHADOWSWIFT_BC == REFLECTIVE_BC
  /* Reflective BC: Solve riemann problem exactly for pressure, all other terms
   * in the flux are 0 because they involve the velocity (which is 0
   * perpendicular to the interface). */
  /* calculate velocity for 1D riemann problem */
  float v = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  float P = riemann_solve_reflective_for_pressure(p->rho, v, p->P);
  totflux[1] = n_unit[0] * P;
  totflux[2] = n_unit[1] * P;
  totflux[3] = n_unit[2] * P;
#endif
  /* Take area and time integrals */
  for (int i = 0; i < 5; i++) totflux[i] *= p->flux.dt * surface_area;
  /* Calculate entropy flux */
  totflux[5] = totflux[0] * p->A;
  hydro_part_update_fluxes_left(p, totflux, dx);
#endif
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_IACT_BOUNDARY_H */
