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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_IACT_H
#define SWIFT_SHADOWSWIFT_HYDRO_IACT_H

#include "hydro_flux.h"
#include "hydro_getters.h"
#include "hydro_gradients.h"
#include "hydro_setters.h"
#include "rt_additions.h"

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
 * @param pi Particle i (the "left" particle). This particle must always be
 * active.
 * @param pj Particle j (the "right" particle).
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 */
__attribute__((always_inline)) INLINE static void runner_iact_flux(
    struct part *pi, struct part *pj, double const *centroid,
    float surface_area, const double *shift) {

  /* Initialize local variables */
  /* Vector from pj to pi */
  float dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = pi->x[k] - pj->x[k] - shift[k];
  }
  const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = (float)sqrt(r2);

  /* Midpoint between pj and pi */
  double midpoint[3];
  for (int k = 0; k < 3; k++) {
    midpoint[k] = 0.5 * (pi->x[k] + pj->x[k] + shift[k]);
  }

  /* Primitive quantities */
  float Wi[5], Wj[5];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* calculate the maximal signal velocity */
  double vmax = 0.0f;
  if (Wi[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(pi->rho, pi->P);
  }
  if (Wj[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(pj->rho, pj->P);
  }

  double dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
               (pi->v[2] - pj->v[2]) * dx[2];
  /* Velocity on the axis linking the particles */
  /* This velocity will be the same as dvdr for MFM, so hopefully this gets
     optimised out. */
  double dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
                   (Wi[3] - Wj[3]) * dx[2];
  /* We only care about this velocity for particles moving towards each others
   */
  dvdotdx = min3(dvdr, dvdotdx, 0.f);

  /* Get the signal velocity */
  vmax -= dvdotdx / r;

  /* Store the signal velocity */
  pi->timestepvars.vmax = (float)fmax(pi->timestepvars.vmax, vmax);
  pj->timestepvars.vmax = (float)fmax(pj->timestepvars.vmax, vmax);

  /* particle velocities */
  double vi[3], vj[3];
  for (int k = 0; k < 3; k++) {
    vi[k] = pi->v[k];
    vj[k] = pj->v[k];
  }

  /* Compute interface velocity, see Springel 2010 (33) */
  float vij[3];
  double fac = ((vj[0] - vi[0]) * (centroid[0] - midpoint[0]) +
                (vj[1] - vi[1]) * (centroid[1] - midpoint[1]) +
                (vj[2] - vi[2]) * (centroid[2] - midpoint[2])) /
               r2;
  vij[0] = 0.5f * (vi[0] + vj[0]) + fac * dx[0];
  vij[1] = 0.5f * (vi[1] + vj[1]) + fac * dx[1];
  vij[2] = 0.5f * (vi[2] + vj[2]) + fac * dx[2];
#if defined(SWIFT_DEBUG_CHECKS) && defined(SHADOWFAX_FIX_CELLS)
  assert(vij[0] == 0.f && vij[1] == 0.f && vij[2] == 0.);
#endif

  float xij_i[3];
  for (int k = 0; k < 3; k++) {
    xij_i[k] = centroid[k] - pi->x[k];
  }
  hydro_gradients_predict(pi, pj, pi->h, pj->h, dx, r, xij_i, /*min_dt,*/ Wi, Wj);

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  /* get the time step for the flux exchange. This is always the smallest time
     step among the two particles */
  const float min_dt = (pj->flux.dt > 0.f)
                           ? fminf(pi->flux.dt, pj->flux.dt)
                           : pi->flux.dt;

#ifdef SWIFT_DEBUG_CHECKS
  assert(pi->flux.dt >= 0);
  assert(min_dt >= 0);
#endif
  if (pj->rho == 0 && pi->rho != 0 && min_dt > 0) {
    pi->fluid_v[0] = pi->fluid_v[0];
  }

  float totflux[5];

  /* compute the normal vector of the interface */
  float n_unit[3];
  for (int k = 0; k < 3; ++k) {
    n_unit[k] = (float)(-dx[k] / r);
  }

  hydro_compute_flux(Wi, Wj, n_unit, vij, surface_area, min_dt, totflux);

  hydro_part_update_fluxes_left(pi, totflux, dx);

  if (pj->flux.dt < 0.0) {
    hydro_part_update_fluxes_right(pj, totflux, dx);
  }
}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

/**
* @brief Calculate the gradient interaction between particle i and particle j
*
* This method wraps around hydro_gradients_collect, which can be an empty
* method, in which case no gradients are used.
*
* @param r2 Comoving squared distance between particle i and particle j.
* @param dx Comoving distance vector between the particles (dx = pi->x -
* pj->x).
* @param hi Comoving smoothing-length of particle i.
* @param hj Comoving smoothing-length of particle j.
* @param pi Particle i.
* @param pj Particle j.
* @param a Current scale factor.
* @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
  hydro_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
* @brief Calculate the gradient interaction between particle i and particle j:
* non-symmetric version
*
* This method wraps around hydro_gradients_nonsym_collect, which can be an
* empty method, in which case no gradients are used.
*
* @param r2 Comoving squared distance between particle i and particle j.
* @param dx Comoving distance vector between the particles (dx = pi->x -
* pj->x).
* @param hi Comoving smoothing-length of particle i.
* @param hj Comoving smoothing-length of particle j.
* @param pi Particle i.
* @param pj Particle j.
* @param a Current scale factor.
* @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
  hydro_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_IACT_H */
