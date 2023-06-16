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
#include "hydro_part.h"
#include "hydro_setters.h"

/**
 * @brief Update the slope estimates of particles pi and pj.
 *
 * @param pi Particle i (the "left" particle). This particle must always be
 * active.
 * @param pj Particle j (the "right" particle).
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Whether to do the slope estimate for both particles.
 */
__attribute__((always_inline)) INLINE static void runner_iact_slope_estimate(
    struct part *pi, struct part *pj, double const *centroid,
    float surface_area, const double *shift, int symmetric) {
#ifdef SHADOWSWIFT_GRADIENTS
  if (!surface_area) {
    /* particle is not a cell neighbour: do nothing */
    return;
  }

  /* Initialize local variables */
  const double dx[3] = {pi->x[0] - pj->x[0] - shift[0],
                        pi->x[1] - pj->x[1] - shift[1],
                        pi->x[2] - pj->x[2] - shift[2]};
  const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* c is supposed to be the vector pointing from the midpoint of pi and pj to
     the centroid of the face between pi and pj.
     The coordinates of the centroid of the face of the voronoi cell of particle
     pi are given in the case of periodic boundary conditions. */
  const double c[3] = {centroid[0] - 0.5 * (pi->x[0] + pj->x[0] + shift[0]),
                       centroid[1] - 0.5 * (pi->x[1] + pj->x[1] + shift[1]),
                       centroid[2] - 0.5 * (pi->x[2] + pj->x[2] + shift[2])};

  /* Update gradient estimate pi */
  double r = sqrt(r2);
  hydro_slope_estimate_collect(pi, pj, c, dx, r, surface_area);

  /* Also update gradient estimate pj? */
  if (symmetric) {
    double mindx[3];
    mindx[0] = -dx[0];
    mindx[1] = -dx[1];
    mindx[2] = -dx[2];
    hydro_slope_estimate_collect(pj, pi, c, mindx, r, surface_area);
  }
#endif
}

/**
 * @brief Collect info necessary for limiting the gradient estimates.
 *
 * @param pi Particle i (the "left" particle). This particle must always be
 * active.
 * @param pj Particle j (the "right" particle).
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Whether to do the slope limiting for both particles.
 */
__attribute__((always_inline)) INLINE static void runner_iact_slope_limiter(
    struct part *pi, struct part *pj, const double *centroid,
    float surface_area, const double *shift, int symmetric) {
#ifdef SHADOWSWIFT_SLOPE_LIMITER_CELL_WIDE
  float f_ij[3] = {centroid[0] - pi->x[0], centroid[1] - pi->x[1],
                   centroid[2] - pi->x[2]};
  hydro_slope_limit_cell_collect(pi, pj, f_ij);

  /* Also treat pj? */
  if (symmetric) {
    float f_ji[3] = {centroid[0] - pj->x[0] - shift[0],
                     centroid[1] - pj->x[1] - shift[1],
                     centroid[2] - pj->x[2] - shift[2]};
    hydro_slope_limit_cell_collect(pj, pi, f_ji);
  }
#endif
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
 * @param pi Particle i (the "left" particle). This particle must always be
 * active.
 * @param pj Particle j (the "right" particle).
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Unused, flux exchange must be manifestly symmetric.
 */
__attribute__((always_inline)) INLINE static void runner_iact_flux_exchange(
    struct part *pi, struct part *pj, double const *centroid,
    float surface_area, const double *shift, int symmetric) {

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
  float Wi[6], Wj[6];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* calculate the maximal signal velocity */
  double vmax = 0.0f;
  if (Wi[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(Wi[0], Wi[4]);
  }
  if (Wj[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(Wj[0], Wj[4]);
  }

  /* Velocity on the axis linking the particles */
  double dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
                   (Wi[3] - Wj[3]) * dx[2];
  /* We only care about this velocity for particles moving towards each others
   */
  dvdotdx = fmin(dvdotdx, 0.f);

  /* Get the signal velocity */
  vmax -= dvdotdx / r;

  /* Store the signal velocity */
  pi->timestepvars.vmax = (float)fmax(pi->timestepvars.vmax, vmax);
  if (pj->flux.dt >= 0)
    pj->timestepvars.vmax = (float)fmax(pj->timestepvars.vmax, vmax);

  /* Store the maximal kinetic energy of the neighbours */
  float mi = pi->conserved.mass;
  float mj = pj->conserved.mass;
  float mi_inv = mi > 0.f ? 1.f / mi : 0.f;
  float mj_inv = mj > 0.f ? 1.f / mj : 0.f;
  float v_rel[3] = {
      mj_inv * pj->conserved.momentum[0] - mi_inv * pi->conserved.momentum[0],
      mj_inv * pj->conserved.momentum[1] - mi_inv * pi->conserved.momentum[1],
      mj_inv * pj->conserved.momentum[2] - mi_inv * pi->conserved.momentum[2]};
  float Ekin_j =
      0.5f * mj *
      (v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]);
  pi->timestepvars.Ekin = fmaxf(pi->timestepvars.Ekin, Ekin_j);
  if (pj->flux.dt >= 0) {
    float Ekin_i =
        0.5f * mi *
        (v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]);
    pj->timestepvars.Ekin = fmaxf(pi->timestepvars.Ekin, Ekin_i);
  }

  /* Compute interface velocity */
  float vij[3];
  hydro_get_interface_velocity(pi->v_full, pj->v_full, dx, midpoint, centroid, vij);

  /* get the time step for the flux exchange. This is always the smallest time
     step among the two particles */
  const float min_dt =
      (pj->flux.dt > 0.f) ? fminf(pi->flux.dt, pj->flux.dt) : pi->flux.dt;

  float xij_i[3];
  for (int k = 0; k < 3; k++) {
    xij_i[k] = centroid[k] - pi->x[k];
  }
  hydro_gradients_predict(pi, pj, dx, r, xij_i, min_dt, Wi, Wj);

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

#ifdef SWIFT_DEBUG_CHECKS
  assert(pi->flux.dt >= 0);
  assert(min_dt >= 0);
#endif

  float totflux[6];

  /* compute the normal vector of the interface */
  float n_unit[3];
  for (int k = 0; k < 3; ++k) {
    n_unit[k] = (float)(-dx[k] / r);
  }

  hydro_compute_flux(Wi, Wj, n_unit, vij, surface_area, min_dt, totflux);

  hydro_part_update_fluxes_left(pi, totflux, dx);
  /* We always update the fluxes for the right particle as well, to make
   * flux exchange manifestly symmetric. */
  hydro_part_update_fluxes_right(pj, totflux, dx);
}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
  error("Shouldn't call this function!");
}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
  error("Shouldn't call this function!");
}

/**
 * @brief Only used with meshless gradients scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
#ifdef SHADOWSWIFT_MESHLESS_GRADIENTS
  hydro_gradients_collect(r2, dx, hi, hj, pi, pj);
#else
  error("Shouldn't call this function!");
#endif
}

/**
 * @brief Only used with meshless gradients scheme.
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
#ifdef SHADOWSWIFT_MESHLESS_GRADIENTS
  /* Only do the gradient interaction if pi also falls within the search
   * radius of pj */
  if (r2 < hj * hj * kernel_gamma2)
    hydro_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
#else
  error("Shouldn't call this function!");
#endif
}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
  error("Shouldn't call this function!");
}

/**
 * @brief Not used in the ShadowSWIFT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {
  error("Shouldn't call this function!");
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_IACT_H */
