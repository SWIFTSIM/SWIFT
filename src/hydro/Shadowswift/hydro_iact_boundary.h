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
#include "riemann/riemann_common.h"

/** @brief Set the primitive quantities of p_boundary to be the ones of p in
 * such a way that the face between p_boundary and p acts as a reflective
 * boundary.
 *
 * Concretely, this function mirrors the density and pressure and reflects the
 * velocity in the rest frame of the face around the face.
 *
 * @param p_boundary The fictive boundary particle
 * @param p The real particle to reflect.
 */
__attribute__((always_inline)) INLINE static void runner_reflect_primitives(
    struct part *p_boundary, const struct part *p, const double *centroid) {
  /* Copy rho and P */
  p_boundary->rho = p->rho;
  p_boundary->P = p->P;

  /* Calculate reflected velocity */
  float dx[3] = {p_boundary->x[0] - p->x[0], p_boundary->x[1] - p->x[1],
                 p_boundary->x[2] - p->x[2]};
  double midpoint[3] = {0.5 * (p->x[0] + p_boundary->x[0]),
                        0.5 * (p->x[1] + p_boundary->x[1]),
                        0.5 * (p->x[2] + p_boundary->x[2])};
  float vij[3];
  hydro_get_interface_velocity(p->v_full, p_boundary->v_full, dx, midpoint,
                               centroid, vij);
  /* Particle velocity in interface frame*/
  float v[3] = {p->v[0] - vij[0], p->v[1] - vij[1], p->v[2] - vij[2]};
  float r_inv = 1.f / sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  float v_dot_n = (v[0] * dx[0] + v[1] * dx[1] + v[2] * dx[2]) * r_inv;
  /* Reflected velocity in lab frame */
  for (int i = 0; i < 3; i++) {
    p_boundary->v[i] = v[i] - 2 * v_dot_n * dx[i] * r_inv + vij[i];
  }
}

/** @brief Set the primitive quantities of an imaginary boundary particle to
 * respect the proper boundary conditions of the simulation
 *
 * @param p_boundary The boundary #part whose primitives get updated.
 * @param p The real #part subjected to boundary conditions
 * @param centroid The centroid of the face between p and p_boundary
 * @param hs The #hydro_space.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_boundary_set_primitives(struct part *p_boundary,
                                    const struct part *p,
                                    const double *centroid,
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
  runner_reflect_primitives(p_boundary, p, centroid);
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
  runner_iact_boundary_set_primitives(p_boundary, p, centroid, hs);
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
  runner_iact_boundary_set_primitives(p_boundary, p, centroid, hs);
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

/** @brief Solve the riemann problem for a reflective boundary exactly.
 *
 * We do this in the frame of reference of the boundary and assume that in
 * that frame (rho_l, v_l, P_l) = (rho_r, -v_r, P_r).
 *
 * @param W Left primitives
 * @param n_unit Normal vector to the interface
 * @param Whalf (return) The primitives at the interface.
 */
__attribute__((always_inline)) INLINE static void
riemann_solve_reflective_boundary(const float *W, const float *n_unit,
                                  float *Whalf) {
  /* Extract local variables for 1D riemann problem */
  float rho = W[0];
  float v = W[1] * n_unit[0] + W[2] * n_unit[1] + W[3] * n_unit[2];
  float P = W[4];

  /* calculate sound speed */
  float a = sqrtf(hydro_gamma * P / rho);

  /* Check for vacuum */
  if (!rho || (hydro_two_over_gamma_minus_one * a <= -v)) {
    riemann_solve_vacuum(W, W, v, -v, a, a, Whalf, n_unit);
    return;
  }
  float P_half, rho_half;
  /* Solve eq. 4.5 from Toro for pressure */
  if (v > 0.f) {
    /* Shock waves */
    float A = hydro_two_over_gamma_plus_one / rho;
    float B = hydro_gamma_minus_one_over_gamma_plus_one * P;
    P_half =
        (v * sqrtf(4.f * P * A + 4.f * A * B + v * v) + 2.f * P * A + v * v) /
        (2.f * A);
    float P_frac = P_half / P;
    rho_half = rho * (hydro_gamma_minus_one_over_gamma_plus_one + P_frac) /
               (hydro_gamma_minus_one_over_gamma_plus_one * P_frac + 1.f);
  } else {
    /* Rarefaction waves */
    float res = 0.5f * v * hydro_gamma_minus_one / a + 1.f;
    P_half = pow_two_gamma_over_gamma_minus_one(res) * P;
    rho_half = rho * pow_one_over_gamma(P_half / P);
  }

  Whalf[0] = rho_half;
  Whalf[1] = W[1] - v * n_unit[0];
  Whalf[2] = W[2] - v * n_unit[1];
  Whalf[3] = W[3] - v * n_unit[2];
  Whalf[4] = P_half;
}

/** @brief For a reflective boundary, we are only interested in the pressure at
 * the boundary since all other terms in the flux will be zero
 * (due to the fact that the velocity at the boundary will be zero).
 * The riemann problem can easily be solved exactly for this case.
 *
 * If the boundary is moving, we solve the reflective riemann problem in the
 * rest frame of the interface.
 *
 * @param p The #part subjected to reflective boundary conditions
 * @param p_boundary The image of p reflected across the boundary, we only use
 * the position of this particle.
 * @param surface_area The surface area of the reflective boundary
 */
__attribute__((always_inline)) INLINE static void
runner_iact_boundary_reflective_flux_exchange(struct part *p,
                                              struct part *p_boundary,
                                              float surface_area,
                                              const double *centroid) {
  /* Vector from pj to pi */
  float dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = (float)(p->x[k] - p_boundary->x[k]);
  }
  float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  float r_inv = 1.f / r;
  /* Normal vector at interface (pointing to pj) */
  float n_unit[3] = {-dx[0] * r_inv, -dx[1] * r_inv, -dx[2] * r_inv};
  /* Interface velocity */
  float vij[3];
  double midpoint[3] = {0.5 * (p->x[0] + p_boundary->x[0]),
                        0.5 * (p->x[1] + p_boundary->x[1]),
                        0.5 * (p->x[2] + p_boundary->x[2])};
  hydro_get_interface_velocity(p->v_full, p_boundary->v_full, dx, midpoint,
                               centroid, vij);
  /* Get primitive variables of pi. */
  float WL[6];
  hydro_part_get_primitive_variables(p, WL);
  /* Boost to rest frame of interface */
  WL[1] -= vij[0];
  WL[2] -= vij[1];
  WL[3] -= vij[2];

  /* Reflective BC: Solve riemann problem exactly (in interface frame) */
  float Whalf[5];
  riemann_solve_reflective_boundary(WL, n_unit, Whalf);
  /* Calculate flux: */
  float totflux[6];
  riemann_flux_from_half_state(Whalf, vij, n_unit, totflux);
  /* Take area and time integrals */
  for (int i = 0; i < 5; i++) totflux[i] *= p->flux.dt * surface_area;
  /* Calculate entropy flux */
  totflux[5] = totflux[0] * p->A;

  /* Check output */
  if (totflux[0] != totflux[0]) error("NaN Mass flux!");
  if (totflux[1] != totflux[1]) error("NaN Velocity flux!");
  if (totflux[2] != totflux[2]) error("NaN Velocity flux!");
  if (totflux[3] != totflux[3]) error("NaN Velocity flux!");
  if (totflux[4] != totflux[4]) error("NaN Energy flux!");
  hydro_part_update_fluxes_left(p, totflux, dx);
}

/**
 * @brief The flux calculation between a "real" particle i and a boundary
 * particle j.
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
#elif SHADOWSWIFT_BC == OPEN_BC
  /* Open BC: Whalf = (rho_l, v_l, P_l), interface velocity = 0.0, so we can
   * calculate the flux exactly. */

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
  /* Get flux vector */
  float vij[3] = {0.f, 0.f, 0.f};
  riemann_flux_from_half_state(WL, vij, n_unit, totflux);
  /* Take area and time integrals */
  for (int i = 0; i < 5; i++) totflux[i] *= p->flux.dt * surface_area;
  /* Calculate entropy flux */
  totflux[5] = totflux[0] * p->A;
  hydro_part_update_fluxes_left(p, totflux, dx);
#elif SHADOWSWIFT_BC == REFLECTIVE_BC
  runner_iact_boundary_reflective_flux_exchange(p, p_boundary, surface_area,
                                                centroid);
#endif
}

/** @brief Do the interaction between an interior and exterior boundary particle
 * if necessary.
 *
 * This function applies to boundary particles inside the simulation (obstacles)
 * and not the global boundary conditions of the simulation. Faces between an
 * interior (only connected to other boundary particles) and exterior boundary
 * particle are treated as reflective. This function assumes that boundary
 * particles are either stationary or move together as a rigid body.
 *
 * @param part_left The left particle to consider
 * @param part_right The right particle to consider
 * @param left_is_active Whether to consider the left particle (we only treat
 * local active boundary particles).
 * @param right_is_active Whether to consider the right particle.
 * @param centroid The centroid of the face between both particles
 * @param surface_area The surface area of the face between both particles.
 * @param shift The shift to apply to the right particle to bring it into the
 * frame of the left
 * @return 1 When we found a pair of boundary particles, 0 otherwise. */
__attribute__((always_inline)) INLINE static int
runner_doiact_boundary_particle(struct part *part_left, struct part *part_right,
                                int left_is_active, int right_is_active,
                                const double *centroid, double surface_area,
                                const double *shift) {
#ifdef SWIFT_BOUNDARY_PARTICLES
  /* Interaction between an interior and exterior boundary particle? */
  if (part_left->id < SWIFT_BOUNDARY_PARTICLES &&
      part_left->id >= space_boundary_parts_interior &&
      part_right->id < space_boundary_parts_interior) {
    /* Anything to do here? */
    if (!left_is_active) return 1;

    /* Create placeholder particle to apply boundary conditions on */
    struct part p_boundary = *part_right;
#if FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE
    runner_iact_boundary_reflective_flux_exchange(
        part_left, &p_boundary, surface_area, centroid);
#else
    runner_reflect_primitives(&p_boundary, part_left, centroid);
    IACT(part_left, &p_boundary, centroid, surface_area, shift, 0);
#endif
    /* Nothing left to do for this pair*/
    return 1;
  } else if (part_right->id < SWIFT_BOUNDARY_PARTICLES &&
             part_right->id >= space_boundary_parts_interior &&
             part_left->id < space_boundary_parts_interior) {
    /* Anything to do here? */
    if (!right_is_active) return 1;

    /* Create placeholder particle to apply boundary conditions on */
    struct part p_boundary = *part_left;
#if FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE
    runner_iact_boundary_reflective_flux_exchange(part_right, &p_boundary,
                                                  surface_area, centroid);
#else
    /* The pair needs to be flipped around */
    /* The midpoint from the reference frame of the right particle */
    double midpoint[3] = {centroid[0] - shift[0], centroid[1] - shift[1],
                          centroid[2] - shift[2]};
    /* The reversed shift */
    double r_shift[3] = {-shift[0], -shift[1], -shift[2]};
    runner_reflect_primitives(&p_boundary, part_right, centroid);
    IACT(part_right, &p_boundary, midpoint, surface_area, r_shift, 0);
#endif
    /* Nothing left to do for this pair*/
    return 1;
  } else if (part_left->id < space_boundary_parts_interior &&
             part_right->id < space_boundary_parts_interior) {
    /* No flux exchange between two interior boundary particles */
    return 1;
  }
#endif
  /* No interaction between boundary particles took place */
  return 0;
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_IACT_BOUNDARY_H */
