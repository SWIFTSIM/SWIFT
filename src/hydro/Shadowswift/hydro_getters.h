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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_GETTERS_H
#define SWIFT_SHADOWSWIFT_HYDRO_GETTERS_H

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
hydro_part_get_primitive_variables(const struct part* restrict p, float W[6]) {

  W[0] = p->rho;
  W[1] = p->v[0];
  W[2] = p->v[1];
  W[3] = p->v[2];
  W[4] = p->P;
  W[5] = p->A;
}

/**
 * @brief Get a 5-element state vector Q containing the conserved hydrodynamic
 * variables.
 *
 * @param p Particle.
 * @param Q Pointer to the array in which the result needs to be stored (of size
 * 5 or more).
 */
__attribute__((always_inline)) INLINE static void
hydro_part_get_conserved_variables(const struct part* restrict p, float* Q) {

  Q[0] = p->conserved.mass;
  Q[1] = p->conserved.momentum[0];
  Q[2] = p->conserved.momentum[1];
  Q[3] = p->conserved.momentum[2];
  Q[4] = p->conserved.energy;
  Q[5] = p->conserved.entropy;
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
    float* dvz, float* dP, float* dA) {

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

  dA[0] = p->gradients.A[0];
  dA[1] = p->gradients.A[1];
  dA[2] = p->gradients.A[2];
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_density(
    const struct part* restrict p) {

  return p->rho;
}

/**
 * @brief Returns the physical density of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    const struct part* restrict p, const struct cosmology* cosmo) {

  return cosmo->a3_inv * p->rho;
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
    return 0.f;
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

  return hydro_get_comoving_internal_energy(p, /*xp*/ NULL);
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
    return 0.f;
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
    return 0.f;
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
    v[0] = p->v[0] + p->flux.momentum[0] * dt_kick_hydro * m_inv;
    v[1] = p->v[1] + p->flux.momentum[1] * dt_kick_hydro * m_inv;
    v[2] = p->v[2] + p->flux.momentum[2] * dt_kick_hydro * m_inv;
  } else {
    v[0] = p->v[0];
    v[1] = p->v[1];
    v[2] = p->v[2];
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
  v_rel[0] = p->v[0] - p->v_full[0];
  v_rel[1] = p->v[1] - p->v_full[1];
  v_rel[2] = p->v[2] - p->v_full[2];
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

  float W[6];
  hydro_part_get_primitive_variables(p, W);

  if (W[0] <= 0.0f) {
    return 0.0f;
  }

  const float rho_inv = 1.f / W[0];

  float gradrho[3], gradvx[3], gradvy[3], gradvz[3], gradP[3], gradA[3];
  hydro_part_get_gradients(p, gradrho, gradvx, gradvy, gradvz, gradP, gradA);

  float v_rel[3];
  hydro_part_get_relative_fluid_velocity(p, v_rel);
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
 * @brief Returns the comoving particle size (~radius).
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_psize(
    const struct part* restrict p) {
  return pow_inv_dimension(p->geometry.volume / hydro_dimension_unit_sphere);
}

/**
 * @brief Returns the physical particle size (~radius).
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_psize(
    const struct part* restrict p, const struct cosmology* cosmo) {
  return cosmo->a * hydro_get_comoving_psize(p);
}

/**
 * @brief Returns the velocity of the interface (at its centroid) between two
 * particles.
 *
 * @param vi The velocity of the first #part (Particle velocity, not fluid
 * velocity!)
 * @param vj The velocity of the second #part
 * @param dx Vector pointing from pj to pi.
 * @param midpoint The midpoint between the two particles.
 * @param centroid The centroid of the interface
 * @param vij (return) The interface velocity
 */
__attribute__((always_inline)) INLINE static void hydro_get_interface_velocity(
    const float* vi, const float* vj, const float* dx, const double* midpoint,
    const double* centroid, float* vij) {
  /* Interface velocity, see Springel 2010 (33) */
  float fac = (float)(((vj[0] - vi[0]) * (centroid[0] - midpoint[0]) +
                       (vj[1] - vi[1]) * (centroid[1] - midpoint[1]) +
                       (vj[2] - vi[2]) * (centroid[2] - midpoint[2])) /
                      (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]));
  vij[0] = 0.5f * (vi[0] + vj[0]) + fac * dx[0];
  vij[1] = 0.5f * (vi[1] + vj[1]) + fac * dx[1];
  vij[2] = 0.5f * (vi[2] + vj[2]) + fac * dx[2];
#if defined(SWIFT_DEBUG_CHECKS) && defined(SHADOWSWIFT_FIX_PARTICLES)
  assert(vij[0] == 0.f && vij[1] == 0.f && vij[2] == 0.);
#endif
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_GETTERS_H */
