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

  /* CFL condition */
  const float dt_cfl = 2.f * const_cfl * kernel_gamma * p->h / p->force.v_sig;

  /* Limit change in u */
  const float dt_u_change = (p->force.u_dt != 0.0f)
                          ? fabsf(const_max_u_change * p->u / p->force.u_dt)
                          : FLT_MAX;

  return fminf(dt_cfl, dt_u_change);
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
  p->density.div_v = 0.f;
  p->density.curl_v[0] = 0.f;
  p->density.curl_v[1] = 0.f;
  p->density.curl_v[2] = 0.f;
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

  /* Final operation on the density. */
  p->rho = ih * ih2 * (p->rho + p->mass * kernel_root);
  p->rho_dh = (p->rho_dh - 3.0f * p->mass * kernel_root) * ih4;
  p->density.wcount =
      (p->density.wcount + kernel_root) * (4.0f / 3.0 * M_PI * kernel_gamma3);
  p->density.wcount_dh =
      p->density.wcount_dh * ih * (4.0f / 3.0 * M_PI * kernel_gamma3);
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * Computes viscosity term, conduction term and smoothing length gradient terms.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param time The current time
 */
__attribute__((always_inline))
INLINE static void hydro_prepare_force(struct part* p, struct xpart* xp, float time) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;
  const float ih2 = ih * ih;
  const float ih4 = ih2 * ih2;

  /* Pre-compute some stuff for the balsara switch. */
  const float normDiv_v = fabs(p->density.div_v / p->rho * ih4);
  const float normCurl_v = sqrtf(p->density.curl_v[0] * p->density.curl_v[0] +
                                 p->density.curl_v[1] * p->density.curl_v[1] +
                                 p->density.curl_v[2] * p->density.curl_v[2]) /
                           p->rho * ih4;

  /* Compute this particle's sound speed. */
  const float u = p->u;
  const float fc = p->force.c =
      sqrtf(const_hydro_gamma * (const_hydro_gamma - 1.0f) * u);

  /* Compute the P/Omega/rho2. */
  xp->omega = 1.0f + 0.3333333333f * h * p->rho_dh / p->rho;
  p->force.POrho2 = u * (const_hydro_gamma - 1.0f) / (p->rho * xp->omega);

  /* Balsara switch */
  p->force.balsara = normDiv_v / (normDiv_v + normCurl_v + 0.0001f * fc * ih);

  /* Viscosity parameter decay time */
  const float tau = h / (2.f * const_viscosity_length * p->force.c);

  /* Viscosity source term */
  const float S = fmaxf(-normDiv_v, 0.f);

  /* Compute the particle's viscosity parameter time derivative */
  const float alpha_dot = (const_viscosity_alpha_min - p->alpha) / tau +
                          (const_viscosity_alpha_max - p->alpha) * S;

  /* Update particle's viscosity paramter */
  p->alpha += alpha_dot * (p->t_end - p->t_begin);
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
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->force.u_dt = 0.0f;
  p->h_dt = 0.0f;
  p->force.v_sig = 0.0f;
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
  float u, w;

  const float dt = t1 - t0;

  /* Predict internal energy */
  w = p->force.u_dt / p->u * dt;
  if (fabsf(w) < 0.01f) /* 1st order expansion of exp(w) */
    u = p->u *=
        1.0f + w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w)));
  else
    u = p->u *= expf(w);

  /* Predict gradient term */
  p->force.POrho2 = u * (const_hydro_gamma - 1.0f) / (p->rho * xp->omega);
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the forces and accelerationsby the appropiate constants
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_end_force(struct part* p) {}



/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
   INLINE static void hydro_kick_extra(struct part* p, float dt) {}

/**
 * @brief Converts hydro quantity of a particle
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_convert_quantities(struct part* p) {}
