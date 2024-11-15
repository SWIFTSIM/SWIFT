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
 * @param dt_kick_corr Time step for the kinetic energy correction due to mass
 * fluxes.
 * @param p Particle.
 * @param a_grav Gravitational acceleration.
 * @param a_grav_prev Gravitational acceleration at the previous full timestep.
 * @param grav_kick Gravitational kick vector also used to update momentum:
 * $dt (m^n a^n + m^(n+1) a^(n+1)$.
 * @return Term used to update the energy variable.
 */
__attribute__((always_inline)) INLINE static float
hydro_gravity_energy_update_term(const float dt_kick_corr,
                                 const struct part* restrict p,
                                 const float* a_grav, const float* a_grav_prev,
                                 const float* grav_kick) {

  float dE =
      -0.5f *
      (p->gravity.mflux[0] *
           (dt_kick_corr * a_grav[0] + p->gravity.dt_corr * a_grav_prev[0]) +
       p->gravity.mflux[1] *
           (dt_kick_corr * a_grav[1] + p->gravity.dt_corr * a_grav_prev[1]) +
       p->gravity.mflux[2] *
           (dt_kick_corr * a_grav[2] + p->gravity.dt_corr * a_grav_prev[2]));

  dE += p->v_full[0] * grav_kick[0] + p->v_full[1] * grav_kick[1] +
        p->v_full[2] * grav_kick[2];

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
 * @brief Applies the gravitational work term at the face between pi and pj to
 * both particles.
 */
__attribute__((always_inline)) INLINE static void
hydro_grav_work_from_half_state(struct part* pi, struct part* pj,
                                const double* shift, const float* Whalf,
                                const float* vij, const double* cij,
                                const float* n_unit, const float area,
                                const float dt) {
  /* Deboost velocity at interface to lab frame */
  const float v_half_lab[3] = {
      vij[0] + Whalf[1],
      vij[1] + Whalf[2],
      vij[2] + Whalf[3],
  };
  const float ri[3] = {
      pi->geometry.centroid[0] + pi->x[0],
      pi->geometry.centroid[1] + pi->x[1],
      pi->geometry.centroid[2] + pi->x[2],
  };
  const float rj[3] = {
      pj->geometry.centroid[0] + pj->x[0] + shift[0],
      pj->geometry.centroid[1] + pj->x[1] + shift[1],
      pj->geometry.centroid[2] + pj->x[2] + shift[2],
  };

  float v_dot_c_i = 0.f;
  float v_dot_c_j = 0.f;
  for (int i = 0; i < 3; i++) {
    v_dot_c_i += (v_half_lab[i] - pi->v_full[i]) * (cij[i] - ri[i]);
    v_dot_c_j += (v_half_lab[i] - pj->v_full[i]) * (cij[i] - rj[i]);
  }
  for (int i = 0; i < 3; i++) {
    pi->gravity.mflux[i] += Whalf[0] * v_dot_c_i * area * n_unit[i];
    pj->gravity.mflux[i] -= Whalf[0] * v_dot_c_j * area * n_unit[i];
  }
}

/**
 * @brief Applies the gravitational work term at the face between pi and pj to
 * both particles.
 *
 * NOTE: This is only an approximation to the (more) exact gravitational work
 * term computed by #hydro_grav_work_from_half_state().
 */
__attribute__((always_inline)) INLINE static void
hydro_grav_work_from_mass_flux(struct part* pi, struct part* pj, float* dx,
                               float mass_flux, const float dt) {
  for (int i = 0; i < 3; i++) {
    pi->gravity.mflux[i] -= 0.5f * mass_flux * dx[i];
    pj->gravity.mflux[i] -= 0.5f * mass_flux * dx[i];
  }
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
