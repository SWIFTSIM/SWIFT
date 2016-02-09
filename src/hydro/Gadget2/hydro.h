/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    struct part* p, struct xpart* xp) {

  /* Acceleration */
  float ac = sqrtf(p->a[0] * p->a[0] + p->a[1] * p->a[1] + p->a[2] * p->a[2]);
  ac = fmaxf(ac, 1e-30);

  const float dt_accel = sqrtf(2.f);  // MATTHIEU

  /* CFL condition */
  const float dt_cfl = 2.f * const_cfl * p->h / p->v_sig;

  return fminf(dt_cfl, dt_accel);
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the variaous density tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_init_part(struct part* p) {
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->rho_dh = 0.f;
  p->div_v = 0.f;
  p->curl_v = 0.f;
  p->rot_v[0] = 0.f;
  p->rot_v[1] = 0.f;
  p->rot_v[2] = 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 *
 * @param p The particle to act upon
 * @param time The current time
 */
__attribute__((always_inline))
    INLINE static void hydro_end_density(struct part* p, float time) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;
  const float ih2 = ih * ih;
  const float ih4 = ih2 * ih2;

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->rho_dh -= 3.0f * p->mass * kernel_root * kernel_igamma;
  p->density.wcount += kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= ih * ih2;
  p->rho_dh *= ih4;
  p->density.wcount *= (4.0f / 3.0 * M_PI * kernel_gamma3);
  p->density.wcount_dh *= ih * (4.0f / 3.0 * M_PI * kernel_gamma3);

  /* Compute the derivative term */
  p->rho_dh = 1.f / (1.f + 0.33333333f * p->h * p->rho_dh / p->rho);

  /* Finish calculation of the velocity curl */
  p->rot_v[0] *= ih4;
  p->rot_v[1] *= ih4;
  p->rot_v[2] *= ih4;

  /* And compute the norm of the curl */
  p->curl_v = sqrtf(p->rot_v[0] * p->rot_v[0] + p->rot_v[1] * p->rot_v[1] +
                    p->rot_v[2] * p->rot_v[2]) /
              p->rho;

  /* Finish calculation of the velocity divergence */
  p->div_v *= ih4 / p->rho;

  /* Compute the pressure */
  const float dt = time - 0.5f * (p->t_begin + p->t_end);
  p->pressure =
      (p->entropy + p->entropy_dt * dt) * powf(p->rho, const_hydro_gamma);
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * Computes viscosity term, conduction term and smoothing length gradient terms.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_prepare_force(struct part* p, struct xpart* xp) {

  /* Some smoothing length multiples. */
  /* const float h = p->h; */
  /* const float ih = 1.0f / h; */
  /* const float ih2 = ih * ih; */
  /* const float ih4 = ih2 * ih2; */

  /* /\* Pre-compute some stuff for the balsara switch. *\/ */
  /* const float normDiv_v = fabs(p->density.div_v / p->rho * ih4); */
  /* const float normCurl_v = sqrtf(p->density.curl_v[0] * p->density.curl_v[0]
   * + */
  /* 				 p->density.curl_v[1] * p->density.curl_v[1] +
   */
  /* 				 p->density.curl_v[2] * p->density.curl_v[2]) /
   */
  /*   p->rho * ih4; */

  /* /\* Compute this particle's sound speed. *\/ */
  /* const float u = p->u; */
  /* const float fc = p->force.c =  */
  /*   sqrtf(const_hydro_gamma * (const_hydro_gamma - 1.0f) * u); */

  /* /\* Compute the P/Omega/rho2. *\/ */
  /* xp->omega = 1.0f + 0.3333333333f * h * p->rho_dh / p->rho; */
  /* p->force.POrho2 = u * (const_hydro_gamma - 1.0f) / (p->rho * xp->omega); */

  /* /\* Balsara switch *\/ */
  /* p->force.balsara = */
  /*   normDiv_v / (normDiv_v + normCurl_v + 0.0001f * fc * ih); */
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking place in the variaous force tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_reset_acceleration(struct part* p) {

  /* Reset the acceleration. */
  p->a[0] = 0.0f;
  p->a[1] = 0.0f;
  p->a[2] = 0.0f;

  /* Reset the time derivatives. */
  p->entropy_dt = 0.0f;

  /* Reset maximal signal velocity */
  p->v_sig = 0.0f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
 * @param t0 The time at the start of the drift
 * @param t1 The time at the end of the drift
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, float t0, float t1) {

  // p->rho *= expf(-p->div_v * dt);
  // p->h *= expf(0.33333333f * p->div_v * dt)

  const float dt_entr = t1 - 0.5f * (p->t_begin + p->t_end);
  p->pressure =
      (p->entropy + p->entropy_dt * dt_entr) * powf(p->rho, const_hydro_gamma);
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the forces and accelerationsby the appropiate constants
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_end_force(struct part* p) {
  
  p->entropy_dt *=
      (const_hydro_gamma - 1.f) * powf(p->rho, -(const_hydro_gamma - 1.f));
}

/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_kick_extra(struct part* p, float dt) {

  /* Do not decrease the entropy (temperature) by more than a factor of 2*/
  const float entropy_change = p->entropy_dt * dt;
  if (entropy_change > -0.5f * p->entropy)
    p->entropy += entropy_change;
  else
    p->entropy *= 0.5f;

  /* Do not 'overcool' when timestep increases */
  if (p->entropy + 0.5f * p->entropy_dt * dt < 0.5f * p->entropy)
    p->entropy_dt = -0.5f * p->entropy / dt;
}

/**
 * @brief Converts hydro quantity of a particle
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_convert_quantities(struct part* p) {

  p->entropy = (const_hydro_gamma - 1.f) * p->entropy *
               powf(p->rho, -(const_hydro_gamma - 1.f));
}
