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
  hydro_gradients_collect(pi, pj, c, dx, r, surface_area);

  /* Also update gradient estimate pj? */
  if (symmetric) {
    double mindx[3];
    mindx[0] = -dx[0];
    mindx[1] = -dx[1];
    mindx[2] = -dx[2];
    hydro_gradients_collect(pj, pi, c, mindx, r, surface_area);
  }
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
  pj->timestepvars.vmax = (float)fmax(pj->timestepvars.vmax, vmax);

  /* particle velocities */
  double vi[3], vj[3];
  for (int k = 0; k < 3; k++) {
    vi[k] = pi->v_full[k];
    vj[k] = pj->v_full[k];
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
#if defined(SWIFT_DEBUG_CHECKS) && defined(SHADOWSWIFT_FIX_PARTICLES)
  assert(vij[0] == 0.f && vij[1] == 0.f && vij[2] == 0.);
#endif

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

  float totflux[5];

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
 * @brief Update the do apoptosis flags of a pair of neighbouring particles
 *
 * This interaction ensures that we do not remove two neighbouring particles at
 * the same time, as this would lead to inconsistencies
 *
 * @param pi Particle i (the "left" particle). In asymmetric mode this is the
 * active particle
 * @param pj Particle j (the "right" particle).
 * @param centroid Centroid of the face between pi and pj.
 * @param surface_area Surface area of the face.
 * @param shift Shift to apply to the coordinates of pj.
 * @param symmetric Whether both particles are active
 */
__attribute__((always_inline)) INLINE static void runner_iact_apoptosis(
    struct part *pi, struct part *pj, double const *centroid,
    float surface_area, const double *shift, int symmetric) {

  /* If only the left particle is active, it is always allowed to be removed */
  if (!symmetric) return;

  /* Is there a conflict? */
  if (pi->derefinement.do_apoptosis && pj->derefinement.do_apoptosis) {
    /* remove the particle that was flagged at the earliest time first */
    if (pi->derefinement.ti_apoptosis < pj->derefinement.ti_apoptosis) {
      pj->derefinement.do_apoptosis = 0;
    } else if (pj->derefinement.ti_apoptosis < pi->derefinement.ti_apoptosis) {
      pi->derefinement.do_apoptosis = 0;
    } else {
      /* Particles were flagged at the same time. Remove the smallest one
       * first... */
      if (pi->geometry.volume < pj->geometry.volume) {
        pj->derefinement.do_apoptosis = 0;
      } else if (pj->geometry.volume < pi->geometry.volume) {
        pi->derefinement.do_apoptosis = 0;
      } else {
        /* Particles also have same volume, remove neither... */
        pi->derefinement.do_apoptosis = 0;
        pj->derefinement.do_apoptosis = 0;
      }
    }
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
 * @brief Not used in ShadowSWIFT.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

/**
 * @brief Not used in ShadowSWIFT
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
    const float H) {}

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
