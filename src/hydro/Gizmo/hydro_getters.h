/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_GETTERS_H
#define SWIFT_GIZMO_HYDRO_GETTERS_H

#include "cosmology.h"
#include "equation_of_state.h"

/**
 * @brief Get a 5-element state vector W containing the primitive hydrodynamic
 * variables.
 *
 * @param p Particle.
 * @param W Pointer to the array in which the result needs to be stored (of size
 * 5 or more).
 */
__attribute__((always_inline)) INLINE static void
hydro_part_get_primitive_variables(const struct part* restrict p, float* W) {

  W[0] = p->rho;
  W[1] = p->fluid_v[0];
  W[2] = p->fluid_v[1];
  W[3] = p->fluid_v[2];
  W[4] = p->P;
}

/**
 * @brief Get the gradients of the primitive variables for the given particle.
 *
 * @param p Particle.
 * @param drho Density gradient (of size 3 or more).
 * @param ddvx x velocity gradient (of size 3 or more).
 * @param ddvy y velocity gradient (of size 3 or more).
 * @param ddvz z velocity gradient (of size 3 or more).
 * @param dP Pressure gradient (of size 3 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_part_get_gradients(
    const struct part* restrict p, float* drho, float* dvx, float* dvy,
    float* dvz, float* dP) {

  drho[0] = p->gradients.rho[0];
  drho[1] = p->gradients.rho[1];
  drho[2] = p->gradients.rho[2];

  dvx[0] = p->gradients.v[0][0];
  dvx[1] = p->gradients.v[0][1];
  dvx[2] = p->gradients.v[0][2];
  dvy[0] = p->gradients.v[1][0];
  dvy[1] = p->gradients.v[1][1];
  dvy[2] = p->gradients.v[1][2];
  dvz[0] = p->gradients.v[2][0];
  dvz[1] = p->gradients.v[2][1];
  dvz[2] = p->gradients.v[2][2];

  dP[0] = p->gradients.P[0];
  dP[1] = p->gradients.P[1];
  dP[2] = p->gradients.P[2];
}

/**
 * @brief Get the slope limiter variables for the given particle.
 *
 * @param p Particle.
 * @param rholim Minimum and maximum density of neighbours (of size 2 or more).
 * @param vxlim Minimum and maximum x velocity of neighbours (of size 2 or
 * more).
 * @param vylim Minimum and maximum y velocity of neighbours (of size 2 or
 * more).
 * @param vzlim Minimum and maximum z velocity of neighbours (of size 2 or
 * more).
 * @param Plim Minimum and maximum pressure of neighbours (of size 2 or more).
 * @param rmax Maximum distance of any neighbour (of size 1 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_part_get_slope_limiter(
    const struct part* restrict p, float* rholim, float* vxlim, float* vylim,
    float* vzlim, float* Plim, float* rmax) {

  rholim[0] = p->limiter.rho[0];
  rholim[1] = p->limiter.rho[1];

  vxlim[0] = p->limiter.v[0][0];
  vxlim[1] = p->limiter.v[0][1];
  vylim[0] = p->limiter.v[1][0];
  vylim[1] = p->limiter.v[1][1];
  vzlim[0] = p->limiter.v[2][0];
  vzlim[1] = p->limiter.v[2][1];

  Plim[0] = p->limiter.P[0];
  Plim[1] = p->limiter.P[1];

  rmax[0] = p->limiter.maxr;
}

/**
 * @brief Returns the comoving internal energy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part* restrict p,
                                   const struct xpart* restrict xp) {

  if (p->rho > 0.0f)
    return gas_internal_energy_from_pressure(p->rho, p->P);
  else
    return 0.;
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy(const struct part* restrict p,
                                   const struct xpart* restrict xp,
                                   const struct cosmology* cosmo) {

  return cosmo->a_factor_internal_energy *
         hydro_get_comoving_internal_energy(p, xp);
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_internal_energy(const struct part* restrict p,
                                           const struct cosmology* cosmo) {

  return hydro_get_physical_internal_energy(p, /*xp=*/NULL, cosmo);
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_internal_energy(const struct part* restrict p) {

  return hydro_get_comoving_internal_energy(p, NULL);
}

/**
 * @brief Returns the comoving entropy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part* restrict p, const struct xpart* restrict xp) {

  if (p->rho > 0.0f) {
    return gas_entropy_from_pressure(p->rho, p->P);
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct cosmology* cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return hydro_get_comoving_entropy(p, NULL);
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(const struct part* restrict p) {

  return hydro_get_comoving_entropy(p, NULL);
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_entropy(const struct part* restrict p,
                                   const struct cosmology* cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return hydro_get_comoving_entropy(p, NULL);
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(const struct part* restrict p) {

  if (p->rho > 0.0f)
    return gas_soundspeed_from_pressure(p->rho, p->P);
  else
    return 0.;
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(const struct part* restrict p,
                              const struct cosmology* cosmo) {

  return cosmo->a_factor_sound_speed * hydro_get_comoving_soundspeed(p);
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    const struct part* restrict p) {

  return p->P;
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_pressure(
    const struct part* restrict p, const struct cosmology* cosmo) {

  return cosmo->a_factor_pressure * p->P;
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    const struct part* restrict p) {

  return p->conserved.mass;
}

/**
 * @brief Returns the velocities drifted to the current time of a particle.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle.
 * @param dt_kick_hydro The time (for hydro accelerations) since the last kick.
 * @param dt_kick_grav The time (for gravity accelerations) since the last kick.
 * @param v (return) The velocities at the current time.
 */
__attribute__((always_inline)) INLINE static void hydro_get_drifted_velocities(
    const struct part* restrict p, const struct xpart* xp, float dt_kick_hydro,
    float dt_kick_grav, float v[3]) {

  if (p->conserved.mass > 0.) {
    const float m_inv = 1.0f / p->conserved.mass;
    v[0] = p->fluid_v[0] + p->flux.momentum[0] * dt_kick_hydro * m_inv;
    v[1] = p->fluid_v[1] + p->flux.momentum[1] * dt_kick_hydro * m_inv;
    v[2] = p->fluid_v[2] + p->flux.momentum[2] * dt_kick_hydro * m_inv;
  } else {
    v[0] = p->fluid_v[0];
    v[1] = p->fluid_v[1];
    v[2] = p->fluid_v[2];
  }

  // MATTHIEU: Bert is this correct? Also, we need to add the mesh kick!
  if (p->gpart) {
    v[0] += p->gpart->a_grav[0] * dt_kick_grav;
    v[1] += p->gpart->a_grav[1] * dt_kick_grav;
    v[2] += p->gpart->a_grav[2] * dt_kick_grav;
  }
}

/**
 * @brief Compute the fluid velocity in the reference frame co-moving with the
 * particle.
 *
 * @param p The #part
 * @param v_rel (return) The relative fluid velocity.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_get_relative_fluid_velocity(const struct part* p, float* v_rel) {
  v_rel[0] = p->fluid_v[0] - p->v[0];
  v_rel[1] = p->fluid_v[1] - p->v[1];
  v_rel[2] = p->fluid_v[2] - p->v[2];
}

/**
 * @brief Returns the time derivative of co-moving internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(const struct part* restrict p) {

  float W[5];
  hydro_part_get_primitive_variables(p, W);

  float v_rel[3];
  hydro_part_get_relative_fluid_velocity(p, v_rel);

  if (W[0] <= 0.0f) {
    return 0.0f;
  }

  const float rho_inv = 1.f / W[0];

  float gradrho[3], gradvx[3], gradvy[3], gradvz[3], gradP[3];
  hydro_part_get_gradients(p, gradrho, gradvx, gradvy, gradvz, gradP);

  const float divv = gradvx[0] + gradvy[1] + gradvz[2];

  float gradu[3] = {0.f, 0.f, 0.f};
  for (int i = 0; i < 3; i++) {
    gradu[i] = hydro_one_over_gamma_minus_one * rho_inv *
               (gradP[i] - rho_inv * W[4] * gradrho[i]);
  }

  const float du_dt =
      -(v_rel[0] * gradu[0] + v_rel[1] * gradu[1] + v_rel[2] * gradu[2]) -
      rho_inv * W[4] * divv;

  return du_dt;
}

/**
 * @brief Returns the time derivative of physical internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy_dt(const struct part* restrict p,
                                      const struct cosmology* cosmo) {

  return hydro_get_comoving_internal_energy_dt(p) *
         cosmo->a_factor_internal_energy;
}

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
hydro_part_geometry_well_behaved(const struct part* restrict p) {

  return p->geometry.wcorr > const_gizmo_min_wcorr;
}

#endif /* SWIFT_GIZMO_HYDRO_GETTERS_H */
