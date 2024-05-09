/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#include "stellar_evolution.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part* p, struct xpart* xp,
			  const struct engine* e) {

  /* Did the particle receive a supernovae */
  if (xp->feedback_data.delta_mass == 0) return;

  const struct cosmology* cosmo = e->cosmology;
  const struct pressure_floor_props* pressure_floor = e->pressure_floor_props;

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
  p->gpart->mass = p->mass ;

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
  if (e->feedback_props->enable_multiple_SN_momentum_correction_factor
      && xp->feedback_data.number_SN > 1) {
    const double p_old[3] = {old_mass*xp->v_full[0], old_mass*xp->v_full[1], old_mass*xp->v_full[2]};
    const double p_old_norm_2 =  p_old[0]*p_old[0] + p_old[1]*p_old[1] + p_old[2]*p_old[2];
    const double p_tilde_norm_2 = p_old_norm_2*dm/old_mass + 2*(old_mass + dm)*xp->feedback_data.delta_E_kin;

    const double dp[3] = {xp->feedback_data.delta_p[0], xp->feedback_data.delta_p[1], xp->feedback_data.delta_p[2]};
    const double dp_norm_2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
    const double p_old_times_dp = p_old[0]*dp[0] + p_old[1]*dp[1] + p_old[2]*dp[2];

    /* Finally compute the corrector factor */
    const double sqrt_argument = fabs(p_old_times_dp*p_old_times_dp - p_tilde_norm_2*dp_norm_2);
    const double f_corr = (- p_old_times_dp + sqrt(sqrt_argument))/p_tilde_norm_2;

    message("f_corr = %e", f_corr);

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
    xp->feedback_data.delta_p[i] = 0;
  }

  /* Reset the values */
  xp->feedback_data.delta_u = 0.;
  xp->feedback_data.delta_E_kin = 0.;
  xp->feedback_data.delta_mass = 0;
  /* xp->feedback_data.number_SN = 0; */
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
void feedback_reset_part(struct part* p, struct xpart* xp) {}

/**
 * @brief Compute the times for the stellar model.
 *
 * This function assumed to be called in the time step task.
 *
 * @param sp The #spart to act upon
 * @param with_cosmology Are we running with the cosmological expansion?
 * @param cosmo The current cosmological model.
 * @param star_age_beg_of_step (output) Age of the star at the beginning of the
 * step.
 * @param dt_enrichment (output) Time step for the stellar evolution.
 * @param ti_begin_star (output) Integer time at the beginning of the time step.
 * @param ti_current The current time (in integer)
 * @param time_base The time base.
 * @param time The current time (in double)
 */
void compute_time(struct spart* sp, const int with_cosmology,
                  const struct cosmology* cosmo, double* star_age_beg_of_step,
                  double* dt_enrichment, integertime_t* ti_begin_star,
                  const integertime_t ti_current, const double time_base,
                  const double time) {
  const integertime_t ti_step = get_integer_timestep(sp->time_bin);
  *ti_begin_star = get_integer_time_begin(ti_current, sp->time_bin);

  /* Get particle time-step */
  double dt_star;
  if (with_cosmology) {
    dt_star = cosmology_get_delta_time(cosmo, *ti_begin_star,
                                       *ti_begin_star + ti_step);
  } else {
    dt_star = get_timestep(sp->time_bin, time_base);
  }

  /* Calculate age of the star at current time */
  double star_age_end_of_step;
  if (with_cosmology) {
    if (cosmo->a > (double)sp->birth_scale_factor)
      star_age_end_of_step = cosmology_get_delta_time_from_scale_factors(
          cosmo, (double)sp->birth_scale_factor, cosmo->a);
    else
      star_age_end_of_step = 0.;
  } else {
    star_age_end_of_step = max(time - (double)sp->birth_time, 0.);
  }

  /* Get the length of the enrichment time-step */
  *dt_enrichment = feedback_get_enrichment_timestep(sp, with_cosmology, cosmo,
                                                    time, dt_star);

  *star_age_beg_of_step = star_age_end_of_step - *dt_enrichment;
}

/**
 * @brief Will this star particle want to do feedback during the next time-step?
 *
 * This is called in the time step task.
 *
 * In GEAR, we compute the full stellar evolution here.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param ti_current The current time (in integer)
 * @param time_base The time base.
 * @param time The physical time in internal units.
 */
void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base) {

  /* Compute the times */
  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin = 0;
  compute_time(sp, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin, ti_current, time_base, time);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;
  sp->feedback_data.will_do_feedback = 0;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
  if (star_age_beg_step + dt_enrichment < 0) {
    error("Negative age for a star");
  }
#endif

  /* Ensure that the age is positive (rounding errors) */
  const double star_age_beg_step_safe =
      star_age_beg_step < 0 ? 0 : star_age_beg_step;

  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metal =
      chemistry_get_star_total_iron_mass_fraction_for_feedback(sp);
  const float threshold = feedback_props->metallicity_max_first_stars;

  const struct stellar_model* model =
      metal < threshold ? &feedback_props->stellar_model_first_stars
                        : &feedback_props->stellar_model;

  /* Compute the stellar evolution including SNe energy */
  stellar_evolution_evolve_spart(sp, model, cosmo, us, phys_const, ti_begin,
                                 star_age_beg_step_safe, dt_enrichment);

  /* Set the particle as doing some feedback */
  sp->feedback_data.will_do_feedback = sp->feedback_data.energy_ejected != 0.;
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart* sp, const struct engine* e) {
  return sp->feedback_data.will_do_feedback;
}

/**
 * @brief Returns the length of time since the particle last did
 * enrichment/feedback.
 *
 * @param sp The #spart.
 * @param with_cosmology Are we running with cosmological time integration on?
 * @param cosmo The cosmological model.
 * @param time The current time (since the Big Bang / start of the run) in
 * internal units.
 * @param dt_star the length of this particle's time-step in internal units.
 * @return The length of the enrichment step in internal units.
 */
double feedback_get_enrichment_timestep(const struct spart* sp,
                                        const int with_cosmology,
                                        const struct cosmology* cosmo,
                                        const double time,
                                        const double dt_star) {
  return dt_star;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart* sp) {
  sp->feedback_data.enrichment_weight = 0.f;

  sp->feedback_data.f_plus_num[0] = 0.0;
  sp->feedback_data.f_plus_num[1] = 0.0;
  sp->feedback_data.f_plus_num[2] = 0.0;

  sp->feedback_data.f_plus_denom[0] = 0.0;
  sp->feedback_data.f_plus_denom[1] = 0.0;
  sp->feedback_data.f_plus_denom[2] = 0.0;

  sp->feedback_data.f_minus_num[0] = 0.0;
  sp->feedback_data.f_minus_num[1] = 0.0;
  sp->feedback_data.f_minus_num[2] = 0.0;

  sp->feedback_data.f_minus_denom[0] = 0.0;
  sp->feedback_data.f_minus_denom[1] = 0.0;
  sp->feedback_data.f_minus_denom[2] = 0.0;

  sp->feedback_data.E_total_accumulator = 0;
  sp->feedback_data.beta_1_accumulator = 0;
  sp->feedback_data.beta_2_accumulator = 0;
  sp->feedback_data.sum_gas_density = 0;
  sp->feedback_data.sum_gas_metallicity = 0;
  sp->feedback_data.density_wcount = 0;

  sp->feedback_data.delta_m_check = 0.0;
  sp->feedback_data.delta_p_norm_check = 0.0;

  sp->feedback_data.delta_p_check[0] = 0 ;
  sp->feedback_data.delta_p_check[1] = 0 ;
  sp->feedback_data.delta_p_check[2] = 0 ;
}

/**
 * @brief Reset the feedback field when the spart is not
 * in a correct state for feeedback_will_do_feedback.
 *
 * This function is called in the timestep task.
 */
void feedback_init_after_star_formation(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Give appropriate values to the star by setting everything to 0 */
  feedback_init_spart(sp);

  sp->feedback_data.energy_ejected = 0.f;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
 */
void feedback_reset_feedback(struct spart* sp,
                             const struct feedback_props* feedback_props) {
  sp->feedback_data.energy_ejected = 0;
  sp->feedback_data.enrichment_weight = 0.f;

  sp->feedback_data.f_plus_num[0] = 0.0;
  sp->feedback_data.f_plus_num[1] = 0.0;
  sp->feedback_data.f_plus_num[2] = 0.0;

  sp->feedback_data.f_plus_denom[0] = 0.0;
  sp->feedback_data.f_plus_denom[1] = 0.0;
  sp->feedback_data.f_plus_denom[2] = 0.0;

  sp->feedback_data.f_minus_num[0] = 0.0;
  sp->feedback_data.f_minus_num[1] = 0.0;
  sp->feedback_data.f_minus_num[2] = 0.0;

  sp->feedback_data.f_minus_denom[0] = 0.0;
  sp->feedback_data.f_minus_denom[1] = 0.0;
  sp->feedback_data.f_minus_denom[2] = 0.0;

  sp->feedback_data.E_total_accumulator = 0;
  sp->feedback_data.beta_1_accumulator = 0;
  sp->feedback_data.beta_2_accumulator = 0;
  sp->feedback_data.sum_gas_density = 0;
  sp->feedback_data.sum_gas_metallicity = 0;
  sp->feedback_data.density_wcount = 0;

  sp->feedback_data.delta_m_check = 0.0;
  sp->feedback_data.delta_p_norm_check = 0.0;

  sp->feedback_data.delta_p_check[0] = 0 ;
  sp->feedback_data.delta_p_check[1] = 0 ;
  sp->feedback_data.delta_p_check[2] = 0 ;
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
void feedback_first_init_spart(struct spart* sp,
                               const struct feedback_props* feedback_props) {
  /* Initialize the feedback struct for the first time */
  feedback_init_spart(sp);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;
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
void feedback_prepare_spart(struct spart* sp,
                            const struct feedback_props* feedback_props) {}

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
void feedback_prepare_feedback(struct spart* restrict sp,
                               const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo,
                               const struct unit_system* us,
                               const struct phys_const* phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology) {  }

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  restart_write_blocks((void*)feedback, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");

  stellar_evolution_dump(&feedback->stellar_model, stream);
  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_dump(&feedback->stellar_model_first_stars, stream);
  }
}

/**
 * @brief Restore a feedback struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {

  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  stellar_evolution_restore(&feedback->stellar_model, stream);

  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_restore(&feedback->stellar_model_first_stars, stream);
  }
}

/**
 * @brief Clean the allocated memory.
 *
 * @param feedback the #feedback_props.
 */
void feedback_clean(struct feedback_props* feedback) {

  stellar_evolution_clean(&feedback->stellar_model);
  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_clean(&feedback->stellar_model_first_stars);
  }
}


/**
 * @brief Compute the scalar weight for the feedback. This scalar weight is
 * used to compute the vector weight. 
 *
 * This function need to be called in loop 1.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle.
 * @param dx_ij_plus (return) Projection vector plus. Pointer to array of size 3.
 * @param dx_ij_minus (return) Projection vector minus. Pointer to array of size 3.
 */
__attribute__((always_inline)) INLINE
double feedback_compute_scalar_weight(const float r2, const float *dx,
				      const float hi, const float hj,
				      const struct spart *restrict si,
				      const struct part *restrict pj,
				      double* dx_ij_plus,
				      double* dx_ij_minus) {

  const float r = sqrtf(r2);

  /* Kernel derivatives evaluation */
  const float u_ij = r/hi;
  const float u_jj = r/hj;
  float dW_ij_dr_j, dW_jj_dr_j, dummy_W;

  kernel_deval(u_ij, &dummy_W, &dW_ij_dr_j);
  kernel_deval(u_jj, &dummy_W, &dW_jj_dr_j);

  dW_ij_dr_j *= pow_dimension_plus_one(1.0/hi); /* 1/h_i^(d+1) */
  dW_jj_dr_j *= pow_dimension_plus_one(1.0/hj); /* 1/h_j^(d+1) */

  /* Ensure they are positive (norm of the gradient) */
  dW_ij_dr_j = fabsf(dW_ij_dr_j); 
  dW_jj_dr_j = fabsf(dW_jj_dr_j);

  /* Compute the projection vectors */
  dx_ij_plus[0] = max(dx[0], 0.0)/r;
  dx_ij_plus[1] = max(dx[1], 0.0)/r;
  dx_ij_plus[2]	= max(dx[2], 0.0)/r;

  dx_ij_minus[0] = min(dx[0], 0.0)/r;
  dx_ij_minus[1] = min(dx[1], 0.0)/r;
  dx_ij_minus[2] = min(dx[2], 0.0)/r;

  /* This is simply dx/r, i.e the unit vector of dx */
  const double dx_ij_hat[3] = {(dx_ij_plus[0] + dx_ij_minus[0]),
			       (dx_ij_plus[1] + dx_ij_minus[1]),
			       (dx_ij_plus[2] + dx_ij_minus[2])};


  /* The star wcount has been computed in the feedback density loop. It need to be
   * multiplied by 1/h^d.
   * The gas wcount cannot be retrived from here. So we approximate it using
   * rho/mass. */
  double n_bar_i = si->feedback_data.density_wcount *  pow_dimension(1.0/hi);
  double n_bar_j = pj->rho/pj->mass;
  double n_bar_i_2_inv = 1.0/(n_bar_i*n_bar_i);
  double n_bar_j_2_inv = 1.0/(n_bar_j*n_bar_j);

  /* Compute the face orientation vector */
  double A_j[3] = {(n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[0],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[1],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[2]};

  /* Prepare the computation of the scalar weight */
  double number_1 = A_j[0]*dx_ij_hat[0] +  A_j[1]*dx_ij_hat[1] + A_j[2]*dx_ij_hat[2];
  double number_2 = M_PI * r2;  
  double denom = sqrt(1 + number_1/number_2);

  /* Finally compute the scalar weight = solid angle fraction the gas j occupy around the star. */
  double scalar_weight_j = 0.5*(1.0 - 1.0/denom) ;

  return scalar_weight_j;
}


/**
 * @brief Compute the non-normalized vector weight for the feedback.
 *
 * This function need to be called after loop 1 (in loop 2 and final feedback
 * computations), i.e. it needs the accumulation of scalar weights to properly
 * compute the vector weights. 
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
__attribute__((always_inline)) INLINE
void feedback_compute_vector_weight_non_normalized(const float r2, const float *dx,
						     const float hi, const float hj,
						     const struct spart *restrict si,
						     const struct part *restrict pj,
						     double* f_plus_i,
						     double* f_minus_i,
						     double* w_j) {
  double dx_ij_plus[3];
  double dx_ij_minus[3];
  double scalar_weight_j = feedback_compute_scalar_weight(r2, dx, hi, hj, si, pj,
							  dx_ij_plus, dx_ij_minus);

  /* Now, that we have accumulated the sums, we can compute the f_plus and
     f_minus */
  const double value_plus[3] = {1 + (si->feedback_data.f_plus_num[0]*si->feedback_data.f_plus_num[0])/(si->feedback_data.f_plus_denom[0]*si->feedback_data.f_plus_denom[0]),
		     1 + (si->feedback_data.f_plus_num[1]*si->feedback_data.f_plus_num[1])/(si->feedback_data.f_plus_denom[1]*si->feedback_data.f_plus_denom[1]),
		     1 + (si->feedback_data.f_plus_num[2]*si->feedback_data.f_plus_num[2])/(si->feedback_data.f_plus_denom[2]*si->feedback_data.f_plus_denom[2])};

  /* Compute the vector factor f_plus */
  f_plus_i[0] = sqrt(0.5*value_plus[0]);
  f_plus_i[1] = sqrt(0.5*value_plus[1]);
  f_plus_i[2] = sqrt(0.5*value_plus[2]);

  const double value_minus[3] = {1 + (si->feedback_data.f_minus_num[0]*si->feedback_data.f_minus_num[0])/(si->feedback_data.f_minus_denom[0]*si->feedback_data.f_minus_denom[0]),
		     1 + (si->feedback_data.f_minus_num[1]*si->feedback_data.f_minus_num[1])/(si->feedback_data.f_minus_denom[1]*si->feedback_data.f_minus_denom[1]),
		     1 + (si->feedback_data.f_minus_num[2]*si->feedback_data.f_minus_num[2])/(si->feedback_data.f_minus_denom[2]*si->feedback_data.f_minus_denom[2])};

  /* Compute the vector factor f_minus */
  f_minus_i[0] = sqrt(0.5*value_minus[0]);
  f_minus_i[1] = sqrt(0.5*value_minus[1]);
  f_minus_i[2] = sqrt(0.5*value_minus[2]);

  /* Now compute the vector weight (non-normalized) */
  w_j[0] = scalar_weight_j*(dx_ij_plus[0]*f_plus_i[0] + dx_ij_minus[0]*f_minus_i[0]);
  w_j[1] = scalar_weight_j*(dx_ij_plus[1]*f_plus_i[1] + dx_ij_minus[1]*f_minus_i[1]);
  w_j[2] = scalar_weight_j*(dx_ij_plus[2]*f_plus_i[2] + dx_ij_minus[2]*f_minus_i[2]);
}


/**
 * @brief Compute the non-normalized vector weight for the feedback.
 *
 * This function need to be called after loop 2, i.e. it needs the accumulation
 * of f_plus and f_minus numerator and denominator to properly compute the
 * noramlized vector weights. 
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
__attribute__((always_inline)) INLINE
void feedback_compute_vector_weight_normalized(const float r2, const float *dx,
						     const float hi, const float hj,
						     const struct spart *restrict si,
						     const struct part *restrict pj,
						     double* w_j_bar) {
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj, f_plus_i, f_minus_i, w_j);

  /* The normalized vector weight */
  w_j_bar[0] = w_j[0]/si->feedback_data.enrichment_weight;
  w_j_bar[1] = w_j[1]/si->feedback_data.enrichment_weight;
  w_j_bar[2] = w_j[2]/si->feedback_data.enrichment_weight;
}

/**
 * @brief Compute the terminal momentum of a SN explosion. This is the momentum
 * the blastwave can give to the gas after the energy-conserving phase.
 *
 * This function is used if we do not resolve the enery-conserving phase.
 *
 * Note: This function compute the terminal momentum in the same way as in
 * Fire-3 (PAPERS).
 *
 * @param sp Star particle undergoing SN explosion.
 * @param p Gas particle receiving the terminal momentum.
 * @param xp The #xpart.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
double feedback_get_SN_terminal_momentum(const struct spart* restrict sp,
					 const struct part* restrict p,
					 const struct xpart* restrict xp,
					 const struct phys_const* phys_const,
					 const struct unit_system* us) {

  /* Terminal momentum 0 (in internal units) */
  const double p_terminal_0 = 2.5e5*phys_const->const_solar_mass*1e-5*units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  /* In erg */
  const double E_ej = sp->feedback_data.energy_ejected*units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  const double ten_to_51 = 1e51;

  /* Get velocity factor. Currently, this is = 1. See PAPER 2024. */
  const double velocity_factor = 1;

  /* Compute the number of neighbours. Its need to multiply wcount by 1/h^d */
  const double n_neighbours = sp->feedback_data.density_wcount * pow_dimension(1.0/sp->h);

  /* Get metallicity factor */
  const double Z_mean = sp->feedback_data.sum_gas_metallicity/(n_neighbours);
  const double Z_sun = 0.0134; /* Find the exact value somewhere */
  double metallicity_factor = 0.0;

  if (Z_mean/Z_sun < 0.01) {
    metallicity_factor = 2;
  } else if ((0.01 <= Z_mean/Z_sun) && (Z_mean/Z_sun<= 1)) {
    metallicity_factor = pow(Z_mean/Z_sun, -0.18);
  } else /* Z/Z_sun > 1 */ {
    metallicity_factor = pow(Z_mean/Z_sun, -0.14);
  }

  /* Get number density factor in cgs */
  const double m_p_cgs = phys_const->const_proton_mass * units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  double density_mean = sp->feedback_data.sum_gas_density / n_neighbours;
  density_mean = density_mean*units_cgs_conversion_factor(us, UNIT_CONV_DENSITY)/m_p_cgs;

  double density_factor = 0.0;

  if (density_mean < 0.001) {
    density_factor = 2.63;
  } else /* >= 0.001 */ {
    density_factor = pow(density_mean, -0.143);
  }

  /* This is in internal units */
  double p_terminal = p_terminal_0*E_ej/ten_to_51*density_factor*metallicity_factor*velocity_factor;
  return p_terminal;
}
