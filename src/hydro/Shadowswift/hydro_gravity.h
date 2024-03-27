//
// Created by yuyttenh on 30/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_GRAVITY_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_GRAVITY_H

/**
 * @brief Add the gravitational contribution to the fluid velocity drift.
 *
 * @param fluid_v Fluid velocity.
 * @param v (Undrifted) particle velocity.
 * @param v_full (Drifted) particle velocity.
 */
__attribute__((always_inline)) INLINE static void
hydro_gravity_extra_velocity_drift(struct part* p) {
  /* This is no longer used */
}

/**
 * @brief Get the term required to update the energy due to the change in
 * gravitational energy.
 *
 * @param dt_kick_corr Time step for the potential energy correction.
 * @param dt_grav Time step for the (optional) kinetic energy correction.
 * @param p Particle.
 * @param momentum Momentum of the particle, explicitly requested so that it is
 * clear from the code that the momentum needs to be updated after the call to
 * this function.
 * @param a_grav Gravitational acceleration.
 * @return Term used to update the energy variable.
 */
__attribute__((always_inline)) INLINE static float
hydro_gravity_energy_update_term(const float dt_kick_corr,
                                 const struct part* restrict p,
                                 const float* momentum, const float* a_grav,
                                 const float* grav_kick_factor) {

  float dE =
      -0.5f * dt_kick_corr *
      (p->gravity.mflux[0] * a_grav[0] + p->gravity.mflux[1] * a_grav[1] +
       p->gravity.mflux[2] * a_grav[2]);
  dE += momentum[0] * grav_kick_factor[0] + momentum[1] * grav_kick_factor[1] +
        momentum[2] * grav_kick_factor[2];
  return dE;
}

/**
 * @brief Get the term required to update the mass due to the mass flux.
 *
 * @param mass_flux Mass flux rate.
 * @param dt Time step (in comoving units).
 * @return Mass flux update term.
 */
__attribute__((always_inline)) INLINE static float
hydro_gravity_mass_update_term(const float mass_flux, const float dt) {
  return mass_flux * dt;
}

/**
 * @brief Update the mass of the gpart associated with the given particle after
 * the mass has been updated with the hydrodynamical mass flux.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_gravity_update_gpart_mass(struct part* restrict p) {

  if (p->gpart) {
    /* Make sure the gpart knows the mass has changed. */
    p->gpart->mass = p->conserved.mass;
  }
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_GRAVITY_H
