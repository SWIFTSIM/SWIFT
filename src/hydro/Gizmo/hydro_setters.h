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
#ifndef SWIFT_GIZMO_HYDRO_SETTERS_H
#define SWIFT_GIZMO_HYDRO_SETTERS_H

#include "pressure_floor.h"

/**
 * @brief Set the primitive variables for the given particle to the given
 * values.
 *
 * @param p Particle.
 * @param W Primitive variables.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_set_primitive_variables(struct part* restrict p, const float* W) {

  p->rho = W[0];
  p->fluid_v[0] = W[1];
  p->fluid_v[1] = W[2];
  p->fluid_v[2] = W[3];
  p->P = W[4];
}

/**
 * @brief Set the conserved variables for the given particle to the given
 * values.
 *
 * @param p Particle.
 * @param Q Conserved variables.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_set_conserved_variables(struct part* restrict p, const float* Q) {

  p->conserved.mass = Q[0];
  p->conserved.momentum[0] = Q[1];
  p->conserved.momentum[1] = Q[2];
  p->conserved.momentum[2] = Q[3];
  p->conserved.energy = Q[4];
}

/**
 * @brief Set the gradients for the given particle to zero.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_gradients(
    struct part* restrict p) {

  p->gradients.rho[0] = 0.0f;
  p->gradients.rho[1] = 0.0f;
  p->gradients.rho[2] = 0.0f;

  p->gradients.v[0][0] = 0.0f;
  p->gradients.v[0][1] = 0.0f;
  p->gradients.v[0][2] = 0.0f;
  p->gradients.v[1][0] = 0.0f;
  p->gradients.v[1][1] = 0.0f;
  p->gradients.v[1][2] = 0.0f;
  p->gradients.v[2][0] = 0.0f;
  p->gradients.v[2][1] = 0.0f;
  p->gradients.v[2][2] = 0.0f;

  p->gradients.P[0] = 0.0f;
  p->gradients.P[1] = 0.0f;
  p->gradients.P[2] = 0.0f;
}

/**
 * @brief Set the gradients for the given particle to the given values.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_set_gradients(
    struct part* restrict p, const float* gradrho, const float* gradvx,
    const float* gradvy, const float* gradvz, const float* gradP) {

  p->gradients.rho[0] = gradrho[0];
  p->gradients.rho[1] = gradrho[1];
  p->gradients.rho[2] = gradrho[2];

  p->gradients.v[0][0] = gradvx[0];
  p->gradients.v[0][1] = gradvx[1];
  p->gradients.v[0][2] = gradvx[2];
  p->gradients.v[1][0] = gradvy[0];
  p->gradients.v[1][1] = gradvy[1];
  p->gradients.v[1][2] = gradvy[2];
  p->gradients.v[2][0] = gradvz[0];
  p->gradients.v[2][1] = gradvz[1];
  p->gradients.v[2][2] = gradvz[2];

  p->gradients.P[0] = gradP[0];
  p->gradients.P[1] = gradP[1];
  p->gradients.P[2] = gradP[2];
}

/**
 * @brief Update the gradients for the given particle with the given
 * contributions.
 *
 * @param p Particle.
 * @param drho Density gradient contribution.
 * @param dvx x velocity gradient contribution.
 * @param dvy y velocity gradient contribution.
 * @param dvz z velocity gradient contribution.
 * @param dP Pressure gradient contribution.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_gradients(
    struct part* restrict p, const float* drho, const float* dvx,
    const float* dvy, const float* dvz, const float* dP) {

  p->gradients.rho[0] += drho[0];
  p->gradients.rho[1] += drho[1];
  p->gradients.rho[2] += drho[2];

  p->gradients.v[0][0] += dvx[0];
  p->gradients.v[0][1] += dvx[1];
  p->gradients.v[0][2] += dvx[2];
  p->gradients.v[1][0] += dvy[0];
  p->gradients.v[1][1] += dvy[1];
  p->gradients.v[1][2] += dvy[2];
  p->gradients.v[2][0] += dvz[0];
  p->gradients.v[2][1] += dvz[1];
  p->gradients.v[2][2] += dvz[2];

  p->gradients.P[0] += dP[0];
  p->gradients.P[1] += dP[1];
  p->gradients.P[2] += dP[2];
}

/**
 * @brief Normalise the gradients for the given particle with the given
 * normalisation factor.
 *
 * @param p Particle.
 * @param norm Normalisation factor.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_normalise_gradients(struct part* restrict p, const float norm) {

  p->gradients.rho[0] *= norm;
  p->gradients.rho[1] *= norm;
  p->gradients.rho[2] *= norm;

  p->gradients.v[0][0] *= norm;
  p->gradients.v[0][1] *= norm;
  p->gradients.v[0][2] *= norm;
  p->gradients.v[1][0] *= norm;
  p->gradients.v[1][1] *= norm;
  p->gradients.v[1][2] *= norm;
  p->gradients.v[2][0] *= norm;
  p->gradients.v[2][1] *= norm;
  p->gradients.v[2][2] *= norm;

  p->gradients.P[0] *= norm;
  p->gradients.P[1] *= norm;
  p->gradients.P[2] *= norm;
}

/**
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(
    struct part* restrict p, float m) {

  p->conserved.mass = m;
}

/**
 * @brief Sets the time derivative of the co-moving internal energy of a
 * particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the comoving internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(struct part* restrict p,
                                      const float du_dt) {

  const float old_du_dt = hydro_get_comoving_internal_energy_dt(p);

  p->flux.energy += p->conserved.mass * (du_dt - old_du_dt) * p->flux.dt;
}

/**
 * @brief Sets the time derivative of the physical internal energy of a particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The time derivative of the physical internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(struct part* restrict p,
                                      const struct cosmology* restrict cosmo,
                                      const float du_dt) {

  hydro_set_comoving_internal_energy_dt(
      p, du_dt / cosmo->a_factor_internal_energy);
}

/**
 * @brief Sets the comoving internal energy of a particle
 *
 * @param p The particle of interest.
 * @param u The comoving internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy(struct part* p, const float u) {

  const float mass = p->conserved.mass;
  if (mass <= 0.0f) {
    return;
  }

  const float Etherm = mass * u;

#ifdef GIZMO_TOTAL_ENERGY
  const float Ekin = 0.5f *
                     (p->conserved.momentum[0] * p->conserved.momentum[0] +
                      p->conserved.momentum[1] * p->conserved.momentum[1] +
                      p->conserved.momentum[2] * p->conserved.momentum[2]) /
                     mass;

  const float Etot = Ekin + Etherm;
  p->conserved.energy = Etot;
#else
  p->conserved.energy = Etherm;
#endif

  p->P = gas_pressure_from_internal_energy(p->rho, u);
}

/**
 * @brief Sets the physical entropy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param entropy The physical entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_physical_entropy(
    struct part* p, struct xpart* xp, const struct cosmology* cosmo,
    const float entropy) {

  const float u = gas_internal_energy_from_entropy(p->rho, entropy);
  hydro_set_comoving_internal_energy(p, u);
}

/**
 * @brief Sets the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy(struct part* p, struct xpart* xp,
                                   const struct cosmology* cosmo,
                                   const float u) {

  hydro_set_comoving_internal_energy(p, u / cosmo->a_factor_internal_energy);
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(
    struct part* p, const struct cosmology* cosmo,
    const struct pressure_floor_props* pressure_floor, const float u) {

  hydro_set_comoving_internal_energy(p, u / cosmo->a_factor_internal_energy);
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part* restrict p, float alpha) {
  /* Purposefully left empty */
}

/**
 * @brief Update the value of the viscosity alpha to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(struct part* restrict p) {
  /* Purposefully left empty */
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
 * @brief Modifies the thermal state of a particle to the imposed internal
 * energy
 *
 * This overrides the current state of the particle but does *not* change its
 * time-derivatives
 *
 * @param p The particle
 * @param u The new internal energy
 */
__attribute__((always_inline)) INLINE static void hydro_set_internal_energy(
    struct part* restrict p, float u) {

  /* conserved.energy is NOT the specific energy (u), but the total thermal
     energy (u*m) */
  p->conserved.energy = u * p->conserved.mass;
#ifdef GIZMO_TOTAL_ENERGY
  /* add the kinetic energy */
  p->conserved.energy += 0.5f * p->conserved.mass *
                         (p->conserved.momentum[0] * p->fluid_v[0] +
                          p->conserved.momentum[1] * p->fluid_v[1] +
                          p->conserved.momentum[2] * p->fluid_v[2]);
#endif
  p->P = hydro_gamma_minus_one * p->rho * u;
}

/**
 * @brief Modifies the thermal state of a particle to the imposed entropy
 *
 * This overrides the current state of the particle but does *not* change its
 * time-derivatives
 *
 * @param p The particle
 * @param S The new entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_entropy(
    struct part* restrict p, float S) {

  p->conserved.energy = S * pow_gamma_minus_one(p->rho) *
                        hydro_one_over_gamma_minus_one * p->conserved.mass;
#ifdef GIZMO_TOTAL_ENERGY
  /* add the kinetic energy */
  p->conserved.energy += 0.5f * p->conserved.mass *
                         (p->conserved.momentum[0] * p->fluid_v[0] +
                          p->conserved.momentum[1] * p->fluid_v[1] +
                          p->conserved.momentum[2] * p->fluid_v[2]);
#endif
  p->P = S * pow_gamma(p->rho);
}

/**
 * @brief Overwrite the initial internal energy of a particle.
 *
 * Note that in the cases where the thermodynamic variable is not
 * internal energy but gets converted later, we must overwrite that
 * field. The conversion to the actual variable happens later after
 * the initial fake time-step.
 *
 * @param p The #part to write to.
 * @param u_init The new initial internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_init_internal_energy(struct part* p, float u_init) {

  /* We store the initial energy per unit mass in the energy
   * variable as the conversion to energy will be done later,
   * in hydro_first_init_part(). */
  p->conserved.energy = u_init;
}

/**
 * @brief Set the *particle* velocity of a particle.
 *
 * This must be the velocity at which the particle is drifted during its
 * timestep (i.e. xp.v_full).
 *
 * @param p The #part to write to.
 * @param v The new particle velocity.
 */
__attribute__((always_inline)) INLINE static void hydro_set_particle_velocity(
    struct part* p, float* v) {
  p->v[0] = v[0];
  p->v[1] = v[1];
  p->v[2] = v[2];
}

#endif /* SWIFT_GIZMO_HYDRO_SETTERS_H */
