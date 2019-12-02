/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_H
#define SWIFT_GIZMO_HYDRO_H

//#define GIZMO_LLOYD_ITERATION

#include "hydro_getters.h"
#include "hydro_setters.h"
#include "hydro_space.h"

#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro_velocities.h"
#endif

#include <float.h>

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param hydro_properties Pointer to the hydro parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct hydro_props* restrict hydro_properties,
    const struct cosmology* restrict cosmo) {

  const float CFL_condition = hydro_properties->CFL_condition;

#ifdef GIZMO_LLOYD_ITERATION
  return CFL_condition;
#endif

  float W[5];
  hydro_part_get_primitive_variables(p, W);

  /* v_full is the actual velocity of the particle, v is its
     hydrodynamical velocity. The time step depends on the relative difference
     of the two. */
  float vrel[3];
  vrel[0] = W[1] - xp->v_full[0];
  vrel[1] = W[2] - xp->v_full[1];
  vrel[2] = W[3] - xp->v_full[2];
  float vmax =
      sqrtf(vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]) +
      sqrtf(hydro_gamma * W[4] / W[0]);
  vmax = max(vmax, p->timestepvars.vmax);

  const float psize = cosmo->a * cosmo->a *
                      powf(p->geometry.volume / hydro_dimension_unit_sphere,
                           hydro_dimension_inv);
  float dt = FLT_MAX;
  if (vmax > 0.0f) {
    dt = psize / vmax;
  }
  return CFL_condition * dt;
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * This method is no longer used, as Gizmo is now unaware of the actual particle
 * time step.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part* p, float dt) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dt == 0.0f) {
    error("Zero time step assigned to particle!");
  }

  if (dt != dt) {
    error("NaN time step assigned to particle!");
  }
#endif
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * In this case, we copy the particle velocities into the corresponding
 * primitive variable field. We do this because the particle velocities in GIZMO
 * can be independent of the actual fluid velocity. The latter is stored as a
 * primitive variable and integrated using the linear momentum, a conserved
 * variable.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    struct part* p, struct xpart* xp) {

  float W[5], Q[5];

  W[0] = 0.0f;
  W[1] = p->v[0];
  W[2] = p->v[1];
  W[3] = p->v[2];
  W[4] = 0.0f;

  Q[0] = p->conserved.mass;
  Q[1] = Q[0] * W[1];
  Q[2] = Q[0] * W[2];
  Q[3] = Q[0] * W[3];
#if defined(EOS_ISOTHERMAL_GAS)
  Q[4] = Q[0] * gas_internal_energy_from_entropy(0.0f, 0.0f);
#else
  Q[4] = p->conserved.energy * Q[0];
#endif

#ifdef GIZMO_TOTAL_ENERGY
  Q[4] += 0.5f * (Q[1] * W[1] + Q[2] * W[2] + Q[3] * W[3]);
#endif

#ifdef GIZMO_LLOYD_ITERATION
  /* overwrite all variables to make sure they have safe values */
  W[0] = 1.0f;
  W[1] = 0.0f;
  W[2] = 0.0f;
  W[3] = 0.0f;
  W[4] = 1.0f;

  Q[0] = 1.0f;
  Q[1] = 0.0f;
  Q[2] = 0.0f;
  Q[3] = 0.0f;
  Q[4] = 1.0f;

  p->v[0] = 0.0f;
  p->v[1] = 0.0f;
  p->v[2] = 0.0f;
#endif

  p->time_bin = 0;

  hydro_part_set_primitive_variables(p, W);
  hydro_part_set_conserved_variables(p, Q);

#if defined(GIZMO_MFV_SPH)
  /* initialize the particle velocity based on the primitive fluid velocity */
  hydro_velocities_init(p, xp);
#elif defined(GIZMO_MFM_SPH)
  /* initialize the particle velocity based on the primitive fluid velocity */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
#endif

  /* ignore accelerations present in the initial condition */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* we cannot initialize wcorr in init_part, as init_part gets called every
     time the density loop is repeated, and the whole point of storing wcorr
     is to have a way of remembering that we need more neighbours for this
     particle */
  hydro_part_set_wcorr(p, 1.0f);
}

/**
 * @brief Prepares a particle for the volume calculation.
 *
 * Simply makes sure all necessary variables are initialized to zero.
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part* p, const struct hydro_space* hs) {

  p->density.wcount = 0.0f;
  p->density.wcount_dh = 0.0f;
  p->geometry.volume = 0.0f;
  p->geometry.matrix_E[0][0] = 0.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 0.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 0.0f;
  p->geometry.centroid[0] = 0.0f;
  p->geometry.centroid[1] = 0.0f;
  p->geometry.centroid[2] = 0.0f;
}

#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro.h"
#define SPH_IMPLEMENTATION "GIZMO MFV (Hopkins 2015)"
#elif defined(GIZMO_MFM_SPH)
#include "MFM/hydro.h"
#define SPH_IMPLEMENTATION "GIZMO MFM (Hopkins 2015)"
#endif

#endif /* SWIFT_GIZMO_HYDRO_H */
