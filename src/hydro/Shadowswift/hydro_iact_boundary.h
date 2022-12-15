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

__attribute__((always_inline)) INLINE static void riemann_solve_reflective(
    const float *WL, const float *WR, const float *n_unit, float *Whalf) {
  /* calculate velocities in interface frame */
  float vL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  float vR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* calculate sound speeds */
  float aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  float aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  if (riemann_is_vacuum(WL, WR, vL, vR, aL, aR)) {
    riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);
    return;
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
#if SHADOWSWIFT_BC == OPEN_BC
  /* Open BC: Whalf = (rho_l, v_l, P_l) */
  float rho_inv = p->rho > 0.f ? 1.f / p->rho : 0.f;
  float totflux[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
  float v_lab;
  if (p->x[0] != p_boundary->x[0]) {
    v_lab = p->v[0];
  } else if (p->x[1] != p_boundary->x[1]) {
    v_lab = p->v[1];
  } else {
    v_lab = p->v[2];
  }
  totflux[0] = p->rho * v_lab;
  float momemtum_flux = p->rho * v_lab * v_lab + P;
  totflux[1] = fsgn(dx[0]) * momentum_flux;
  totflux[2] = fsgn(dx[1]) * momentum_flux;
  totflux[3] = fsgn(dx[2]) * momentum_flux;
  totflux[4] = p->rho * (0.5f * v_lab * v_lab +
                         hydro_one_over_gamma_minus_one * p->P * rho_inv) +
               p->P * v_lab;
#elif SHADOWSWIFT_BC == REFLECTIVE_BC
  /* Reflective BC: Solve riemann problem exactly */
  float WL[5], WR[5], Whalf[3];
  float n_unit[3] = {fsgnf(dx[0]), fsgnf(dx[1]), fsgnf(dx[2])};
  hydro_part_get_primitive_variables(p, WL);
  hydro_part_get_primitive_variables(p_boundary, WR);

  riemann_solve_reflective(WL, WR, n_unit, Whalf);

#ifdef SWIFT_DEBUG_CHECKS
  /* TODO check that correct terms are 0 */
#endif

  float totflux[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
  if (n_unit[0] != 0.f) {

  } else if (n_unit[1] != 0.f) {

  } else {

  }
#endif
  for (int i = 0; i < 5; i++) totflux[i] *= p->flux.dt * surface_area;
  hydro_part_update_fluxes_left(p, totflux, dx);
#endif
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_IACT_BOUNDARY_H */
