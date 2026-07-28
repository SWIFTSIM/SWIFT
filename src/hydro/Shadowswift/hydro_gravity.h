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
 * @param v_full Particle velocity for this timestep
 * @return Term used to update the energy variable.
 */
__attribute__((always_inline)) INLINE static float
hydro_gravity_energy_update_term(const float dt_kick_corr1,
                                 const float dt_kick_corr2,
                                 const float* a_grav1, const float* a_grav2,
                                 const float* mflux1, const float* mflux2,
                                 const float* v_full1, const float* v_full2,
                                 const float m1dt1, const float m2dt2,
                                 const float* grav_kick) {

  /* Developers note: this is follows the springel 2010 Eq 94 and all following
   * equations that are based from it. In some cases, you will see a 1/2
   * prefactor for the contributions. This is a result of averaging over the
   * two timestep contributions. However, if dt_kick_corr is the actual
   * half-timestep, this eliminates the need for an average, and the sum of the
   * two time integrations is one full timestep */

  /* Solve as in Springel 2010 Eq 94. Carry out as in Hopkins 2015 Eq H2.
   * This is a slightly modified method of Springel, but we do both kicks
   * at the same time, but as an average of n and n+1.
   *
   * These terms are not particularly efficiently calculated, but seeing as
   * this term alone has been the cause of an enormous amount of headache,
   * instructive and readable code has been prioritised.
   */

  /* Start with momentum kicks */
  float dE_momentum = v_full1[0] * grav_kick[0] + v_full1[1] * grav_kick[1] +
        v_full1[2] * grav_kick[2];

  float dE_momentum1 = m1dt1 *
    (a_grav1[0] * v_full1[0] +
      a_grav1[1] * v_full1[1] +
        a_grav1[2] * v_full1[2]);

  float dE_momentum2 = m2dt2 *
    (a_grav2[0] * v_full2[0] +
      a_grav2[1] * v_full2[1] +
        a_grav2[2] * v_full2[2]);

  /* Contribution from timestep n */
  float grav_work1 = dt_kick_corr1 *
    (a_grav1[0] * mflux1[0] +
      a_grav1[1] * mflux1[1] +
        a_grav1[2] * mflux1[2]);

  /* Contribution from timestep n+1 */
  float grav_work2 = dt_kick_corr2 *
    (a_grav2[0] * mflux2[0] +
      a_grav2[1] * mflux2[1] +
        a_grav2[2] * mflux2[2]);

  const float dE2 = dE_momentum1 + dE_momentum2 + grav_work1 + grav_work2;

  return dE2;
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
 *
 * @param pi, pj The particles
 * @param shift Shift to appie to pj's coordinates, if wrapping
 * around simulation box
 * @param Whalf The state of primitive vectors at the interface (output from
 * riemann solver, in reference frame of interface)
 * @param vij Interface velocity
 * @param cij Interface centroid
 * @param n_unit Interface normal. Note that this is also the direction from
 * pi->x to pj->x (voronoi faces are perpendicular to the corresponding Delaunay
 * edges).
 * @param area Interface area
 * @param dt the timestep over which fluxes are exchanged. (currently unused)
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
    v_dot_c_i += (v_half_lab[i] - pi->v_part_full[i]) * (cij[i] - ri[i]);
    v_dot_c_j += (v_half_lab[i] - pj->v_part_full[i]) * (cij[i] - rj[i]);
  }
  for (int i = 0; i < 3; i++) {

    /* Ammended: Integrate the mass exchange! */
    pi->gravity.mflux[i] += Whalf[0] * v_dot_c_i * area * n_unit[i] * dt;
    pj->gravity.mflux[i] -= Whalf[0] * v_dot_c_j * area * n_unit[i] * dt;

    /* Care for directions here.. Not fully convinced of everything */


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

  /* Now defunct */
}

/**
 * @brief Applies the gravitational work term at the face between pi and pj to
 * both particles.
 *
 * NOTE: This is only an approximation to the (more) exact gravitational work
 * term computed by #hydro_grav_work_from_half_state().
 *
 * This operates under a different assumption than
 * hydro_grav_work_from_mass_flux(): This simplifies many terms, this is shown
 * in the Thesis eq 233. It assumes the mass flux AND the distance travelled
 * in this work term.
 *
 * The first approximation concerns mass flux. I assume this assumption is fine.
 * The second is that the work done is carried out over a distance that is
 * half the distance between two centroids, ie that the face centroid is
 * halfway between centroid j and centroid i. This is not correct.
 *
 * In this approximation we assume that the mass exchanged for this work is
 * the mass flux, but that the distance over which we exert the work is
 * from the centroid of the cell to the face of the cell.
 *
 */
__attribute__((always_inline)) INLINE static void
hydro_grav_work_from_mass_flux_half_approximated(struct part* pi,
  struct part* pj, float mass_flux,
  const double* shift, const double* cij, const float dt) {


  /* Ideally apply to centroid, since the potential is evaluated there
   * There is potential to test this with generator positions, but since grav
   * is evaluated through centroid, it would be better to do this.
   */
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

  for (int i = 0; i < 3; i++) {
    /* So: Following the Hopkins H2 prescription, it should be
     * Mwork,i = (ri - cij)*dm_i
     *         = -(ri - cij)dm_ij
     *
     * Mwork, j = (rj - cij)*dm_j
     *          = (rj - cij)*dm_ij
     *
     * where dm_i is the mass flux into i and dm_j follows
     * dm_ij is the mass_flux found from the solver. Recall that we will do
     * pi->flux.mass -= mass_flux, hence the source of minus term.
     *
     * of course, mass_flux here is a rate, we multiply by dt to get an actual
     * exchange of mass
     *
     * Recall: These are "corrective terms" to energy flux, the signs are
     * + / - but the result is the same effect since the vectors are opposite.
     * If mass leaves i into j, and acceleration is into j, this will
     * reduce the energy loss of i (term is + in i) but increase the energy
     * gain inside j (also + in j). Counter-intuitive to think, but this is then
     * basically a source/sink term of energy.
     *
     */
    pi->gravity.mflux[i] -= dt * mass_flux * (ri[i] - cij[i]);
    pj->gravity.mflux[i] += dt * mass_flux * (rj[i] - cij[i]);
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

/**
 * @brief Update xp->mflux to be the mflux vector at end of kick2. Used in next
 * timestep as a mflux at timestep n value.
 *
 * Also sets the previous dt, used to rescale the mflux accordingly.
 *
 * @param p Particle.
 * @param xp The extra part of p
 */
__attribute__((always_inline)) INLINE static void
hydro_gravity_xp_mflux(struct part* p, struct xpart* xp) {

  for (int i = 0; i < 3; i++) {
    xp->mflux[i] = p->gravity.mflux[i];
    p->flux.dt_previous = p->flux.dt;
  }
}


#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_GRAVITY_H
