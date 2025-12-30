/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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

/* Include header */
#include "feedback.h"

/* Local includes */
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part *p, struct xpart *xp,
                          const struct engine *e) {

  /* Did the particle receive a supernovae */
  if (xp->feedback_data.delta_mass == 0) return;

  const struct cosmology *cosmo = e->cosmology;
  const struct pressure_floor_props *pressure_floor = e->pressure_floor_props;

  /* Turn off the cooling */
  xp->cooling_data.time_last_event = e->time;

  /* Update mass */
  const float dm = xp->feedback_data.delta_mass;
  const float old_mass = hydro_get_mass(p);
  const float new_mass = old_mass + dm;

  if (xp->feedback_data.delta_mass < 0.) {
    error("Delta mass smaller than 0");
  }

  /* Update the mass of p, as well as its gpart's friend */
  hydro_set_mass(p, new_mass);

  /* Update the density */
  p->rho *= new_mass / old_mass;

  /* Update internal energy */
  const float u =
      hydro_get_physical_internal_energy(p, xp, cosmo) * old_mass / new_mass;
  const float u_new = u + xp->feedback_data.delta_u;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

  /* Compute correction therm to account for multiple feedback events. This
     terms allows to recover energy conservation. If there is only one
     feedback that affected p and xp, f_corr = 1. */
  if (e->feedback_props->enable_multiple_SN_momentum_correction_factor &&
      xp->feedback_data.number_SN > 1) {
    const double f_corr =
        feedback_compute_momentum_correction_factor_for_multiple_sn_events(
            p, xp, old_mass, new_mass);

    /* Update the xpart accumulated dp */
    xp->feedback_data.delta_p[0] *= f_corr;
    xp->feedback_data.delta_p[1] *= f_corr;
    xp->feedback_data.delta_p[2] *= f_corr;
  }

  /* Update the velocities */
  for (int i = 0; i < 3; i++) {
    const float dv = xp->feedback_data.delta_p[i] / new_mass;

    xp->v_full[i] += dv;
    p->v[i] += dv;

    /* Reset the values */
    xp->feedback_data.delta_p[i] = 0.0;
  }

  /* Reset the values */
  xp->feedback_data.delta_u = 0.0;
  xp->feedback_data.delta_E_kin = 0.0;
  xp->feedback_data.delta_mass = 0.0;
  xp->feedback_data.number_SN = 0.0;
}

/**
 * @brief Finishes the #part density calculation.
 *
 * Save the #part density.wcount for use in feedback loops.
 *
 * @param p The particle to act upon
 * @param xp The extra particle to act upon
 */
__attribute__((always_inline)) INLINE void feedback_end_density(
    struct part *p, struct xpart *xp) {
  p->feedback_data.density.wcount = p->density.wcount;
}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here in the GEAR model.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
void feedback_reset_part(struct part *p, struct xpart *xp) {}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart *sp, const struct engine *e) {
  /* The particle is inactive if its birth_scale_factor or birth_time is
   * negative */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return 0;

  return sp->feedback_data.will_do_feedback;
}

/**
 * @brief Should this particle inject anything as supernovae feedback?
 *
 * This function is used in feedback_iact.h to determine if we have SN feedback
 * to inject in the gas.
 *
 * @param sp The #spart.
 */
int feedback_should_inject_SN_feedback(const struct spart *sp) {
  /* Do we have supernovae? */
  return sp->feedback_data.energy_ejected != 0;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart *sp) {
  sp->feedback_data.enrichment_weight = 0.f;

  sp->feedback_data.f_sum_minus_term[0] = 0.0;
  sp->feedback_data.f_sum_minus_term[1] = 0.0;
  sp->feedback_data.f_sum_minus_term[2] = 0.0;

  sp->feedback_data.f_sum_plus_term[0] = 0.0;
  sp->feedback_data.f_sum_plus_term[1] = 0.0;
  sp->feedback_data.f_sum_plus_term[2] = 0.0;

  sp->feedback_data.accumulator.E_total = 0.0;
  sp->feedback_data.accumulator.beta_1 = 0.0;
  sp->feedback_data.accumulator.beta_2 = 0.0;
  sp->feedback_data.weighted_gas_density = 0.0;
  sp->feedback_data.weighted_gas_metallicity = 0.0;

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  sp->feedback_data.fluxes_conservation_check.delta_m = 0.0;
  sp->feedback_data.fluxes_conservation_check.delta_p_norm = 0.0;

  sp->feedback_data.fluxes_conservation_check.delta_p[0] = 0.0;
  sp->feedback_data.fluxes_conservation_check.delta_p[1] = 0.0;
  sp->feedback_data.fluxes_conservation_check.delta_p[2] = 0.0;
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
 */
void feedback_reset_feedback(struct spart *sp,
                             const struct feedback_props *feedback_props) {
  sp->feedback_data.energy_ejected = 0.0;
  sp->feedback_data.enrichment_weight = 0.f;

  sp->feedback_data.f_sum_minus_term[0] = 0.0;
  sp->feedback_data.f_sum_minus_term[1] = 0.0;
  sp->feedback_data.f_sum_minus_term[2] = 0.0;

  sp->feedback_data.f_sum_plus_term[0] = 0.0;
  sp->feedback_data.f_sum_plus_term[1] = 0.0;
  sp->feedback_data.f_sum_plus_term[2] = 0.0;

  sp->feedback_data.accumulator.E_total = 0.0;
  sp->feedback_data.accumulator.beta_1 = 0.0;
  sp->feedback_data.accumulator.beta_2 = 0.0;
  sp->feedback_data.weighted_gas_density = 0.0;
  sp->feedback_data.weighted_gas_metallicity = 0.0;

#ifdef SWIFT_FEEDBACK_DEBUG_CHECKS
  sp->feedback_data.fluxes_conservation_check.delta_m = 0.0;
  sp->feedback_data.fluxes_conservation_check.delta_p_norm = 0.0;

  sp->feedback_data.fluxes_conservation_check.delta_p[0] = 0.0;
  sp->feedback_data.fluxes_conservation_check.delta_p[1] = 0.0;
  sp->feedback_data.fluxes_conservation_check.delta_p[2] = 0.0;
#endif /* SWIFT_FEEDBACK_DEBUG_CHECKS */
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
void feedback_prepare_spart(struct spart *sp,
                            const struct feedback_props *feedback_props) {}

/**
 * @brief Prepare a #spart for the feedback task.
 *
 * This is called in the stars ghost task.
 *
 * In here, we only need to add the missing coefficients.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
void feedback_prepare_feedback(struct spart *restrict sp,
                               const struct feedback_props *feedback_props,
                               const struct cosmology *cosmo,
                               const struct unit_system *us,
                               const struct phys_const *phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology) {}

/**
 * @brief Compute the scalar weight for the feedback. This scalar weight is
 * used to compute the vector weight.
 *
 * This function needs to be called in loop 1.
 *
 * Note: i = star, j = gas.
 * Note 2: This is scale-factor free --> no conversion needed from/to comoving.
 *
 * Reference: https://arxiv.org/abs/1707.07010
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle.
 * @param dx_ij_plus (return) Projection vector plus. Pointer to array of
 * size 3.
 * @param dx_ij_minus (return) Projection vector minus. Pointer to array of
 * size 3.
 * @param scalar_weigth_j (return) Scalar weight.
 */
__attribute__((always_inline)) INLINE void feedback_compute_scalar_weight(
    const float r2, const float *dx, const float hi, const float hj,
    const struct spart *restrict si, const struct part *restrict pj,
    double dx_ij_plus[3], double dx_ij_minus[3], double *scalar_weight_j) {

  const float r = sqrtf(r2);

  /* Kernel derivatives evaluation */
  const float u_ij = r / hi;
  const float u_jj = r / hj;
  float dW_ij_dr_j, dW_jj_dr_j, dummy_W;

  kernel_deval(u_ij, &dummy_W, &dW_ij_dr_j);
  kernel_deval(u_jj, &dummy_W, &dW_jj_dr_j);

  dW_ij_dr_j *= pow_dimension_plus_one(1.0 / hi); /* 1/h_i^(d+1) */
  dW_jj_dr_j *= pow_dimension_plus_one(1.0 / hj); /* 1/h_j^(d+1) */

  /* Ensure they are positive (norm of the gradient) */
  dW_ij_dr_j = fabsf(dW_ij_dr_j);
  dW_jj_dr_j = fabsf(dW_jj_dr_j);

  /* Compute the projection vectors (scale-factors cancel out) */
  dx_ij_plus[0] = max(dx[0], 0.0) / r;
  dx_ij_plus[1] = max(dx[1], 0.0) / r;
  dx_ij_plus[2] = max(dx[2], 0.0) / r;

  dx_ij_minus[0] = min(dx[0], 0.0) / r;
  dx_ij_minus[1] = min(dx[1], 0.0) / r;
  dx_ij_minus[2] = min(dx[2], 0.0) / r;

  /* This is simply dx/r, i.e the unit vector of dx (scale-factors cancel out)
   */
  const double dx_ij_hat[3] = {(dx_ij_plus[0] + dx_ij_minus[0]),
                               (dx_ij_plus[1] + dx_ij_minus[1]),
                               (dx_ij_plus[2] + dx_ij_minus[2])};

  /* The star wcount has been computed in the feedback density loop. It need to
   * be multiplied by 1/h^d. The gas wcount cannot be retrived from here. So we
   * approximate it using rho/mass. */
  /* Note: Gizmo does not compute the matrix E for stars, as we can understand
     from the paper. Moreover, Gizmo uses wcount = rho/mass for the gas and not
     the MFV/M volume. */
  const double n_bar_i = si->density.wcount;
  const double n_bar_j = pj->feedback_data.density.wcount;
  const double n_bar_i_2_inv = 1.0 / (n_bar_i * n_bar_i);
  const double n_bar_j_2_inv = 1.0 / (n_bar_j * n_bar_j);

  /* Compute the face orientation vector */
  const double A_j[3] = {
      (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j) * dx_ij_hat[0],
      (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j) * dx_ij_hat[1],
      (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j) * dx_ij_hat[2]};

  /* Prepare the computation of the scalar weight.
     Notice that the scale factors cancel out between number_1 (a^2) and
     number_2 (a^2). */
  const double number_1 =
      A_j[0] * dx_ij_hat[0] + A_j[1] * dx_ij_hat[1] + A_j[2] * dx_ij_hat[2];
  const double number_2 = M_PI * r2;
  const double denom = sqrt(1 + number_1 / number_2);

  /* Finally compute the scalar weight = solid angle fraction the gas j occupy
   * around the star. */
  *scalar_weight_j = 0.5 * (1.0 - 1.0 / denom);
}

/**
 * @brief Compute the non-normalized vector weight for the feedback.
 *
 * This function needs to be called after loop 1 (in loop 2 and final feedback
 * computations), i.e. it needs the accumulation of scalar weights to properly
 * compute the vector weights.
 *
 * Note: This is scale-factor free --> no conversion needed from/to comoving.
 *
 * Reference: https://arxiv.org/abs/1707.07010
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle.
 * @param f_plus_i (return) Vector factor f_minus. Pointer to array of size 3.
 * @param f_minus_i (return) Vector factor f_minus. Pointer to array of size 3.
 * @param w_j (return) Non-noralized vector weight. Pointer to array of size 3.
 */
__attribute__((always_inline)) INLINE void
feedback_compute_vector_weight_non_normalized(
    const float r2, const float *dx, const float hi, const float hj,
    const struct spart *restrict si, const struct part *restrict pj,
    double f_plus_i[3], double f_minus_i[3], double w_j[3]) {
  double dx_ij_plus[3], dx_ij_minus[3], scalar_weight_j;
  feedback_compute_scalar_weight(r2, dx, hi, hj, si, pj, dx_ij_plus,
                                 dx_ij_minus, &scalar_weight_j);

  /* Now, that we have accumulated the sums, we can compute the f_plus and
     f_minus */
  double value_plus[3] = {1 + (si->feedback_data.f_sum_minus_term[0] *
                               si->feedback_data.f_sum_minus_term[0]) /
                                  (si->feedback_data.f_sum_plus_term[0] *
                                   si->feedback_data.f_sum_plus_term[0]),
                          1 + (si->feedback_data.f_sum_minus_term[1] *
                               si->feedback_data.f_sum_minus_term[1]) /
                                  (si->feedback_data.f_sum_plus_term[1] *
                                   si->feedback_data.f_sum_plus_term[1]),
                          1 + (si->feedback_data.f_sum_minus_term[2] *
                               si->feedback_data.f_sum_minus_term[2]) /
                                  (si->feedback_data.f_sum_plus_term[2] *
                                   si->feedback_data.f_sum_plus_term[2])};

  /* In rare cases, f_sum_plus_term can have a component that is 0. Since
     division by 0 will give inf and further operations will give NaNs, give a
     large number to represent inf. */
  if (isinf(value_plus[0])) {
    value_plus[0] = FLT_MAX;
  }
  if (isinf(value_plus[1])) {
    value_plus[1] = FLT_MAX;
  }
  if (isinf(value_plus[2])) {
    value_plus[2] = FLT_MAX;
  }

  /* Compute the vector factor f_plus */
  f_plus_i[0] = sqrt(0.5 * value_plus[0]);
  f_plus_i[1] = sqrt(0.5 * value_plus[1]);
  f_plus_i[2] = sqrt(0.5 * value_plus[2]);

  double value_minus[3] = {1 + (si->feedback_data.f_sum_plus_term[0] *
                                si->feedback_data.f_sum_plus_term[0]) /
                                   (si->feedback_data.f_sum_minus_term[0] *
                                    si->feedback_data.f_sum_minus_term[0]),
                           1 + (si->feedback_data.f_sum_plus_term[1] *
                                si->feedback_data.f_sum_plus_term[1]) /
                                   (si->feedback_data.f_sum_minus_term[1] *
                                    si->feedback_data.f_sum_minus_term[1]),
                           1 + (si->feedback_data.f_sum_plus_term[2] *
                                si->feedback_data.f_sum_plus_term[2]) /
                                   (si->feedback_data.f_sum_minus_term[2] *
                                    si->feedback_data.f_sum_minus_term[2])};

  /* In rare cases, f_sum_minus_term can have a component that is 0. Since
     division by 0 will give inf and further operations will give NaNs, give a
     large number to represent inf. */
  if (isinf(value_minus[0])) {
    value_minus[0] = FLT_MAX;
  }
  if (isinf(value_minus[1])) {
    value_minus[1] = FLT_MAX;
  }
  if (isinf(value_minus[2])) {
    value_minus[2] = FLT_MAX;
  }

  /* Compute the vector factor f_minus */
  f_minus_i[0] = sqrt(0.5 * value_minus[0]);
  f_minus_i[1] = sqrt(0.5 * value_minus[1]);
  f_minus_i[2] = sqrt(0.5 * value_minus[2]);

  /* Now compute the vector weight (non-normalized) */
  w_j[0] = scalar_weight_j *
           (dx_ij_plus[0] * f_plus_i[0] + dx_ij_minus[0] * f_minus_i[0]);
  w_j[1] = scalar_weight_j *
           (dx_ij_plus[1] * f_plus_i[1] + dx_ij_minus[1] * f_minus_i[1]);
  w_j[2] = scalar_weight_j *
           (dx_ij_plus[2] * f_plus_i[2] + dx_ij_minus[2] * f_minus_i[2]);
}

/**
 * @brief Compute the non-normalized vector weight for the feedback.
 *
 * This function need to be called after loop 2, i.e. it needs the accumulation
 * of f_plus and f_minus numerator and denominator to properly compute the
 * noramlized vector weights.
 *
 * This is scale-factor free --> no conversion needed from/to comoving.
 *
 * Reference: https://arxiv.org/abs/1707.07010
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle.
 * @param f_plus_i (return) Vector factor f_minus. Pointer to array of size 3.
 * @param f_minus_i (return) Vector factor f_minus. Pointer to array of size 3.
 * @param w_j_bar (return) Normalized vector weight. Pointer to array of size 3.
 */
__attribute__((always_inline)) INLINE void
feedback_compute_vector_weight_normalized(const float r2, const float *dx,
                                          const float hi, const float hj,
                                          const struct spart *restrict si,
                                          const struct part *restrict pj,
                                          double w_j_bar[3]) {
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj,
                                                f_plus_i, f_minus_i, w_j);

  /* The normalized vector weight */
  w_j_bar[0] = w_j[0] / si->feedback_data.enrichment_weight;
  w_j_bar[1] = w_j[1] / si->feedback_data.enrichment_weight;
  w_j_bar[2] = w_j[2] / si->feedback_data.enrichment_weight;
}

/**
 * @brief Compute the physical terminal momentum of a SN explosion. This is the
 *  momentum the blastwave can give to the gas after the energy-conserving
 * phase.
 *
 * This function is used if we do not resolve the enery-conserving phase.
 *
 * Note: This function compute the terminal momentum in the same way as in
 * Fire-3 (https://arxiv.org/abs/2404.16987).
 *
 * @param sp Star particle undergoing SN explosion.
 * @param p Gas particle receiving the terminal momentum.
 * @param xp The #xpart.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE double
feedback_get_physical_SN_terminal_momentum(const struct spart *restrict sp,
                                           const struct part *restrict p,
                                           const struct xpart *restrict xp,
                                           const struct phys_const *phys_const,
                                           const struct unit_system *us,
                                           const struct cosmology *cosmo) {

  /* Terminal momentum 0 (in internal units). Note the 1e-5 term since we want
     it in km and not cm. */
  const double p_terminal_0 =
      2.5e5 * phys_const->const_solar_mass * 1e-5 *
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  /* In erg */
  const double E_ej = sp->feedback_data.energy_ejected *
                      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  const double ten_to_51 = 1e51;

  /* Get velocity factor. Currently, this is = 1. See the attached paper. */
  const double velocity_factor = 1;

  /* Get metallicity factor */
  const double Z_mean = sp->feedback_data.weighted_gas_metallicity;
  const double Z_sun = 0.0134;
  double metallicity_factor = 0.0;

  if (Z_mean / Z_sun < 0.01) {
    metallicity_factor = 2;
  } else if ((0.01 <= Z_mean / Z_sun) && (Z_mean / Z_sun <= 1)) {
    metallicity_factor = pow(Z_mean / Z_sun, -0.18);
  } else /* Z/Z_sun > 1 */ {
    metallicity_factor = pow(Z_mean / Z_sun, -0.14);
  }

  /* Get number density factor in cgs */
  const double m_p_cgs = phys_const->const_proton_mass *
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double density_mean =
      sp->feedback_data.weighted_gas_density * cosmo->a3_inv *
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) / m_p_cgs;

  double density_factor = 0.0;

  if (density_mean < 0.001) {
    density_factor = 2.63;
  } else /* >= 0.001 */ {
    density_factor = pow(density_mean, -0.143);
  }

  /* This is in internal units */
  double p_terminal = p_terminal_0 * E_ej / ten_to_51 * density_factor *
                      metallicity_factor * velocity_factor;
  return p_terminal;
}

/**
 * @brief Compute the physical supernova cooling radius.
 *
 * This function is used if we do not resolve the enery-conserving phase.
 *
 * Note: This functions is used by the mechanical feedback mode 1.
 *
 * @param sp Star particle undergoing SN explosion.
 * @param p Gas particle receiving the terminal momentum.
 * @param xp The #xpart.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE float
feedback_get_physical_SN_cooling_radius(const struct spart *restrict sp,
                                        float p_SN_initial, float p_terminal,
                                        const struct cosmology *cosmo) {

  const float m_ej = sp->feedback_data.mass_ejected;

  /* Convert to physical units */
  const float mean_density =
      sp->feedback_data.weighted_gas_density * cosmo->a3_inv;

  /* Compute the cooling radius */
  const float second_part =
      p_terminal * p_terminal / (p_SN_initial * p_SN_initial) - 1;
  const float r_cool =
      pow(3.0 * m_ej * second_part / (4.0 * M_PI * mean_density), 1.0 / 3.0);

  return r_cool;
}

/**
 * @brief Compute a correction factor to ensure energy conservation when
 * multiple feedback (supernovae) events affect the same #part in a given
 * timestep.
 *
 * Note: This function is called in feedback_update_part().
 *
 * This is scale-factor free --> no conversion needed from/to comoving.
 *
 * Reference: https://arxiv.org/abs/2203.00040
 *
 * @param p The #part to correct.
 * @param xp The #xpart.
 * @param old_mass The mass before feeback events.
 * @param new_mass The mass after feeback events.
 */
__attribute__((always_inline)) INLINE double
feedback_compute_momentum_correction_factor_for_multiple_sn_events(
    struct part *p, struct xpart *xp, float old_mass, float new_mass) {
  const float dm = xp->feedback_data.delta_mass;

  const double p_old[3] = {old_mass * xp->v_full[0], old_mass * xp->v_full[1],
                           old_mass * xp->v_full[2]};
  const double p_old_norm_2 =
      p_old[0] * p_old[0] + p_old[1] * p_old[1] + p_old[2] * p_old[2];
  const double p_tilde_norm_2 =
      p_old_norm_2 * dm / old_mass +
      2 * (old_mass + dm) * xp->feedback_data.delta_E_kin;

  const double dp[3] = {xp->feedback_data.delta_p[0],
                        xp->feedback_data.delta_p[1],
                        xp->feedback_data.delta_p[2]};
  const double dp_norm_2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  const double p_old_times_dp =
      p_old[0] * dp[0] + p_old[1] * dp[1] + p_old[2] * dp[2];

  /* Finally compute the corrector factor */
  const double sqrt_argument =
      fabs(p_old_times_dp * p_old_times_dp - p_tilde_norm_2 * dp_norm_2);
  const double f_corr =
      (-p_old_times_dp + sqrt(sqrt_argument)) / p_tilde_norm_2;

  return f_corr;
}
