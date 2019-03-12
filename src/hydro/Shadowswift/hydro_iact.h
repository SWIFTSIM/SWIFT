/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#include "adiabatic_index.h"
#include "hydro_gradients.h"
#include "riemann.h"
#include "voronoi_algorithm.h"

/**
 * @brief Calculate the Voronoi cell by interacting particle pi and pj
 *
 * This method wraps around voronoi_cell_interact().
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float mindx[3];

  voronoi_cell_interact(&pi->cell, dx, pj->id);
  mindx[0] = -dx[0];
  mindx[1] = -dx[1];
  mindx[2] = -dx[2];
  voronoi_cell_interact(&pj->cell, mindx, pi->id);
}

/**
 * @brief Calculate the Voronoi cell by interacting particle pi with pj
 *
 * This method wraps around voronoi_cell_interact().
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  voronoi_cell_interact(&pi->cell, dx, pj->id);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  hydro_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  hydro_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Common part of the flux calculation between particle i and j
 *
 * Since the only difference between the symmetric and non-symmetric version
 * of the flux calculation  is in the update of the conserved variables at the
 * very end (which is not done for particle j if mode is 0 and particle j is
 * active), both runner_iact_force and runner_iact_nonsym_force call this
 * method, with an appropriate mode.
 *
 * This method retrieves the oriented surface area and face midpoint for the
 * Voronoi face between pi and pj (if it exists). It uses the midpoint position
 * to reconstruct the primitive quantities (if gradients are used) at the face
 * and then uses the face quantities to estimate a flux through the face using
 * a Riemann solver.
 *
 * This method also calculates the maximal velocity used to calculate the time
 * step.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_fluxes_common(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, int mode, float a, float H) {

  float r = sqrtf(r2);
  int k;
  float A;
  float xij_i[3];
  float vmax, dvdotdx;
  float vi[3], vj[3], vij[3];
  float Wi[5], Wj[5];
  float n_unit[3];

  A = voronoi_get_face(&pi->cell, pj->id, xij_i);
  if (A == 0.0f) {
    /* this neighbour does not share a face with the cell, return */
    return;
  }

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    vi[k] = pi->force.v_full[k]; /* particle velocities */
    vj[k] = pj->force.v_full[k];
  }
  Wi[0] = pi->primitives.rho;
  Wi[1] = pi->primitives.v[0];
  Wi[2] = pi->primitives.v[1];
  Wi[3] = pi->primitives.v[2];
  Wi[4] = pi->primitives.P;
  Wj[0] = pj->primitives.rho;
  Wj[1] = pj->primitives.v[0];
  Wj[2] = pj->primitives.v[1];
  Wj[3] = pj->primitives.v[2];
  Wj[4] = pj->primitives.P;

  /* calculate the maximal signal velocity */
  vmax = 0.0f;
  if (Wi[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(Wi[0], Wi[4]);
  }

  if (Wj[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(Wj[0], Wj[4]);
  }

  dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
            (Wi[3] - Wj[3]) * dx[2];
  if (dvdotdx > 0.) {
    vmax -= dvdotdx / r;
  }

  pi->timestepvars.vmax = fmaxf(pi->timestepvars.vmax, vmax);
  if (mode == 1) {
    pj->timestepvars.vmax = fmaxf(pj->timestepvars.vmax, vmax);
  }

  /* compute the normal vector of the interface */
  for (k = 0; k < 3; ++k) {
    n_unit[k] = -dx[k] / r;
  }

  /* Compute interface velocity */
  float fac = (vi[0] - vj[0]) * (xij_i[0] + 0.5f * dx[0]) +
              (vi[1] - vj[1]) * (xij_i[1] + 0.5f * dx[1]) +
              (vi[2] - vj[2]) * (xij_i[2] + 0.5f * dx[2]);
  fac /= r;
  vij[0] = 0.5f * (vi[0] + vj[0]) - fac * dx[0];
  vij[1] = 0.5f * (vi[1] + vj[1]) - fac * dx[1];
  vij[2] = 0.5f * (vi[2] + vj[2]) - fac * dx[2];

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  hydro_gradients_predict(pi, pj, hi, hj, dx, r, xij_i, Wi, Wj);

  /* we don't need to rotate, we can use the unit vector in the Riemann problem
   * itself (see GIZMO) */

  if (Wi[0] < 0.0f || Wj[0] < 0.0f || Wi[4] < 0.0f || Wj[4] < 0.0f) {
    printf("WL: %g %g %g %g %g\n", pi->primitives.rho, pi->primitives.v[0],
           pi->primitives.v[1], pi->primitives.v[2], pi->primitives.P);
#ifdef USE_GRADIENTS
    printf("dWL: %g %g %g %g %g\n", dWi[0], dWi[1], dWi[2], dWi[3], dWi[4]);
#endif
    printf("gradWL[0]: %g %g %g\n", pi->primitives.gradients.rho[0],
           pi->primitives.gradients.rho[1], pi->primitives.gradients.rho[2]);
    printf("gradWL[1]: %g %g %g\n", pi->primitives.gradients.v[0][0],
           pi->primitives.gradients.v[0][1], pi->primitives.gradients.v[0][2]);
    printf("gradWL[2]: %g %g %g\n", pi->primitives.gradients.v[1][0],
           pi->primitives.gradients.v[1][1], pi->primitives.gradients.v[1][2]);
    printf("gradWL[3]: %g %g %g\n", pi->primitives.gradients.v[2][0],
           pi->primitives.gradients.v[2][1], pi->primitives.gradients.v[2][2]);
    printf("gradWL[4]: %g %g %g\n", pi->primitives.gradients.P[0],
           pi->primitives.gradients.P[1], pi->primitives.gradients.P[2]);
    printf("WL': %g %g %g %g %g\n", Wi[0], Wi[1], Wi[2], Wi[3], Wi[4]);
    printf("WR: %g %g %g %g %g\n", pj->primitives.rho, pj->primitives.v[0],
           pj->primitives.v[1], pj->primitives.v[2], pj->primitives.P);
#ifdef USE_GRADIENTS
    printf("dWR: %g %g %g %g %g\n", dWj[0], dWj[1], dWj[2], dWj[3], dWj[4]);
#endif
    printf("gradWR[0]: %g %g %g\n", pj->primitives.gradients.rho[0],
           pj->primitives.gradients.rho[1], pj->primitives.gradients.rho[2]);
    printf("gradWR[1]: %g %g %g\n", pj->primitives.gradients.v[0][0],
           pj->primitives.gradients.v[0][1], pj->primitives.gradients.v[0][2]);
    printf("gradWR[2]: %g %g %g\n", pj->primitives.gradients.v[1][0],
           pj->primitives.gradients.v[1][1], pj->primitives.gradients.v[1][2]);
    printf("gradWR[3]: %g %g %g\n", pj->primitives.gradients.v[2][0],
           pj->primitives.gradients.v[2][1], pj->primitives.gradients.v[2][2]);
    printf("gradWR[4]: %g %g %g\n", pj->primitives.gradients.P[0],
           pj->primitives.gradients.P[1], pj->primitives.gradients.P[2]);
    printf("WR': %g %g %g %g %g\n", Wj[0], Wj[1], Wj[2], Wj[3], Wj[4]);
    error("Negative density or pressure!\n");
  }

  float totflux[5];
  riemann_solve_for_flux(Wi, Wj, n_unit, vij, totflux);

  /* Update conserved variables */
  /* eqn. (16) */
  pi->conserved.flux.mass -= A * totflux[0];
  pi->conserved.flux.momentum[0] -= A * totflux[1];
  pi->conserved.flux.momentum[1] -= A * totflux[2];
  pi->conserved.flux.momentum[2] -= A * totflux[3];
  pi->conserved.flux.energy -= A * totflux[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
  float ekin = 0.5f * (pi->primitives.v[0] * pi->primitives.v[0] +
                       pi->primitives.v[1] * pi->primitives.v[1] +
                       pi->primitives.v[2] * pi->primitives.v[2]);
  pi->conserved.flux.energy += A * totflux[1] * pi->primitives.v[0];
  pi->conserved.flux.energy += A * totflux[2] * pi->primitives.v[1];
  pi->conserved.flux.energy += A * totflux[3] * pi->primitives.v[2];
  pi->conserved.flux.energy -= A * totflux[0] * ekin;
#endif

  /* here is how it works:
     Mode will only be 1 if both particles are ACTIVE and they are in the same
     cell. In this case, this method IS the flux calculation for particle j, and
     we HAVE TO UPDATE it.
     Mode 0 can mean several things: it can mean that particle j is INACTIVE, in
     which case we NEED TO UPDATE it, since otherwise the flux is lost from the
     system and the conserved variable is not conserved.
     It can also mean that particle j sits in another cell and is ACTIVE. In
     this case, the flux exchange for particle j is done TWICE and we SHOULD NOT
     UPDATE particle j.
     ==> we update particle j if (MODE IS 1) OR (j IS INACTIVE)
  */
  if (mode == 1 || pj->force.active == 0) {
    pj->conserved.flux.mass += A * totflux[0];
    pj->conserved.flux.momentum[0] += A * totflux[1];
    pj->conserved.flux.momentum[1] += A * totflux[2];
    pj->conserved.flux.momentum[2] += A * totflux[3];
    pj->conserved.flux.energy += A * totflux[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
    ekin = 0.5f * (pj->primitives.v[0] * pj->primitives.v[0] +
                   pj->primitives.v[1] * pj->primitives.v[1] +
                   pj->primitives.v[2] * pj->primitives.v[2]);
    pj->conserved.flux.energy -= A * totflux[1] * pj->primitives.v[0];
    pj->conserved.flux.energy -= A * totflux[2] * pj->primitives.v[1];
    pj->conserved.flux.energy -= A * totflux[3] * pj->primitives.v[2];
    pj->conserved.flux.energy += A * totflux[0] * ekin;
#endif
  }
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_fluxes_common with mode 1.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 1, a, H);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * This method calls runner_iact_fluxes_common with mode 0.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 0, a, H);
}

/**
 * @brief Timestep limiter loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_limiter(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  /* Nothing to do here if both particles are active */
}

/**
 * @brief Timestep limiter loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_limiter(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  /* Wake up the neighbour? */
  if (pi->timestepvars.vmax >
      const_limiter_max_v_sig_ratio * pj->timestepvars.vmax) {

    pj->wakeup = time_bin_awake;
  }
}
