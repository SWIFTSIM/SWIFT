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
#ifndef SWIFT_GEAR_MECHANICAL_FEEDBACK_IACT_H
#define SWIFT_GEAR_MECHANICAL_FEEDBACK_IACT_H

/* Local includes */
#include "feedback.h"
#include "hydro.h"
#include "mechanical_feedback_iact.h"
#include "random.h"
#include "timestep_sync_part.h"

#include <math.h>

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * In GEAR, this function does nothing. What we need is the
 * star->density.wcount computed in runner_iact_nonsym_stars_density().
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *si, const struct part *pj,
                                    const struct xpart *xpj,
                                    const struct cosmology *cosmo,
                                    const struct feedback_props *fb_props,
                                    const integertime_t ti_current) {

  const float r_max_2 = fb_props->r_max * fb_props->r_max;

  /* If the particle is farther than the maximal radius, it does not receive
     feedback. Hence, do not count it. */
  if (r2 > r_max_2) {
    return;
  }

  /* Do we have supernovae? */
  if (!feedback_should_inject_SN_feedback(si)) {
    return;
  }
}

/**
 * @brief Prepare the feedback by computing the required quantities (loop 1).
 * Used for updating properties of gas particles required for the feedback.
 *
 * Note: Not used in GEAR_mechanical.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const struct feedback_props *fb_props,
                                  const integertime_t ti_current) {

  const float r_max_2 = fb_props->r_max * fb_props->r_max;

  /* If the particle is farther than the maximal radius, it does not receive
     feedback */
  if (r2 > r_max_2) {
    return;
  }

  /* Do we have supernovae? */
  if (!feedback_should_inject_SN_feedback(si)) {
    return;
  }
}

/**
 * @brief Prepare the feedback by computing the required quantities (loop 2).
 * Used for updating properties of star particles required for the feedback.
 *
 * In GEAR, compute the first quantities required for the enrichment vector
 * weights.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const struct feedback_props *fb_props,
                                  const integertime_t ti_current) {

  const float r_max_2 = fb_props->r_max * fb_props->r_max;

  /* If the particle is farther than the maximal radius, it does not receive
     feedback */
  if (r2 > r_max_2) {
    return;
  }

  /* Do we have supernovae? */
  if (!feedback_should_inject_SN_feedback(si)) {
    return;
  }

  /* Accumulate the sum in the numerator and denominator of f_plus and f_minus
   */
  double dx_ij_plus[3], dx_ij_minus[3], scalar_weight_j;
  feedback_compute_scalar_weight(r2, dx, hi, hj, si, pj, dx_ij_plus,
                                 dx_ij_minus, &scalar_weight_j);

  si->feedback_data.f_sum_plus_term[0] += scalar_weight_j * fabs(dx_ij_plus[0]);
  si->feedback_data.f_sum_plus_term[1] += scalar_weight_j * fabs(dx_ij_plus[1]);
  si->feedback_data.f_sum_plus_term[2] += scalar_weight_j * fabs(dx_ij_plus[2]);

  si->feedback_data.f_sum_minus_term[0] +=
      scalar_weight_j * fabs(dx_ij_minus[0]);
  si->feedback_data.f_sum_minus_term[1] +=
      scalar_weight_j * fabs(dx_ij_minus[1]);
  si->feedback_data.f_sum_minus_term[2] +=
      scalar_weight_j * fabs(dx_ij_minus[2]);
}

/**
 * @brief Prepare the feedback by computing the required quantities (loop 2).
 * Used for updating properties of star particles required for the feedback.
 *
 * In GEAR, compute the second quantities required for the enrichment vector
 * weights.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep3(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const struct feedback_props *fb_props,
                                  const integertime_t ti_current) {

  const float r_max_2 = fb_props->r_max * fb_props->r_max;

  /* If the particle is farther than the maximal radius, it does not receive
     feedback */
  if (r2 > r_max_2) {
    return;
  }

  /* Do we have supernovae? */
  if (!feedback_should_inject_SN_feedback(si)) {
    return;
  }

  /* Now we can compute f_plus and f_minus for the star */
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj,
                                                f_plus_i, f_minus_i, w_j);

  /* Accumulate w_j norm for later */
  const double w_j_norm_2 = w_j[0] * w_j[0] + w_j[1] * w_j[1] + w_j[2] * w_j[2];
  si->feedback_data.enrichment_weight += sqrt(w_j_norm_2);
}

#if FEEDBACK_GEAR_MECHANICAL_MODE == 2
/**
 * @brief Prepare the feedback by computing the required quantities (loop 3).
 *
 * Only for mechanical feedback 2. Accumulate values to compute quantities that
 * take into account the star-gas motion.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep4(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const struct feedback_props *fb_props,
                                  const integertime_t ti_current) {

  const float r_max_2 = fb_props->r_max * fb_props->r_max;

  /* If the particle is farther than the maximal radius, it does not receive
     feedback */
  if (r2 > r_max_2) {
    return;
  }

  /* Do we have supernovae? */
  if (!feedback_should_inject_SN_feedback(si)) {
    return;
  }

  /* Compute w_j_bar. */
  double w_j_bar[3];
  feedback_compute_vector_weight_normalized(r2, dx, hi, hj, si, pj, w_j_bar);
  const double w_j_bar_norm_2 = w_j_bar[0] * w_j_bar[0] +
                                w_j_bar[1] * w_j_bar[1] +
                                w_j_bar[2] * w_j_bar[2];
  const double w_j_bar_norm = sqrt(w_j_bar_norm_2);

  /* If p does not contribute, skip the computations to avoid NaN */
  if (w_j_bar_norm == 0) {
    return;
  }

  /* Get some properties for our computations */
  const float mj = hydro_get_mass(pj);
  const float m_ej = si->feedback_data.mass_ejected;
  const double dm = max(w_j_bar_norm * m_ej, FLT_MIN);

  /* Accumulate (pay attention to the conversions to physical units) */
  const float v_ij[3] = {pj->v[0] - si->v[0], pj->v[1] - si->v[1],
                         pj->v[2] - si->v[2]};

  /* Calculate the velocity with the Hubble flow */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const float a2H = a * a * H;
  const float v_ij_plus_H_flow[3] = {
      a2H * dx[0] + v_ij[0], a2H * dx[1] + v_ij[1], a2H * dx[2] + v_ij[2]};

  /* Compute the _physical_ relative velocity between the particles */
  const float v_ij_p[3] = {v_ij_plus_H_flow[0] * cosmo->a_inv,
                           v_ij_plus_H_flow[1] * cosmo->a_inv,
                           v_ij_plus_H_flow[2] * cosmo->a_inv};

  const float v_ij_p_norm_2 =
      v_ij_p[0] * v_ij_p[0] + v_ij_p[1] * v_ij_p[1] + v_ij_p[2] * v_ij_p[2];

  /* w_j_bar_hat refers to w_j_bar/|w_j_bar| */
  const double v_ij_p_times_w_j_bar_hat =
      (v_ij_p[0] * w_j_bar[0] + v_ij_p[1] * w_j_bar[1] +
       v_ij_p[2] * w_j_bar[2]) /
      w_j_bar_norm;
  const double w_prime_ij = w_j_bar_norm / (1 + dm / mj);

  /* Notice that we will multiply by 0.5*m_ej later on */
  si->feedback_data.accumulator.E_total += w_prime_ij * v_ij_p_norm_2;

  /* Notice that we need the small epsilon (total available kinetic energy) to
     finish the computation of this. The small epsilon is determined by E_tot */
  si->feedback_data.accumulator.beta_1 += w_prime_ij * v_ij_p_times_w_j_bar_hat;

  /* Notice that we will multiply by m_ej later on */
  si->feedback_data.accumulator.beta_2 += w_prime_ij * w_j_bar_norm / mj;

  /* Compute the comoving weigthed average of the gas properties around the star
     with our isotropic weighting scheme. */
  si->feedback_data.weighted_gas_density += w_j_bar_norm * pj->rho;
  si->feedback_data.weighted_gas_metallicity +=
      w_j_bar_norm * chemistry_get_total_metal_mass_fraction_for_feedback(pj);
}

#endif /*  FEEDBACK_GEAR_MECHANICAL_MODE == 2 */

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param fb_props Properties of the feedback scheme.
 * @param constants The physical constants (in internal units).
 * @param us The internal system of units.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, const integertime_t ti_current) {

  const float r_max_2 = fb_props->r_max * fb_props->r_max;

  /* If the particle is farther than the maximal radius, it does not receive
     feedback */
  if (r2 > r_max_2) {
#ifdef SWIFT_DEBUG_CHECKS
    warning(
        "Particle %lld has a distance (%e) bigger than r_max = %e. It will not "
        "receive feedback!",
        pj->id, sqrt(r2), fb_props->r_max);
#endif
    return;
  }

  /* Compute the w_j_bar. */
  double w_j_bar[3];
  feedback_compute_vector_weight_normalized(r2, dx, hi, hj, si, pj, w_j_bar);
  const double w_j_bar_norm_2 = w_j_bar[0] * w_j_bar[0] +
                                w_j_bar[1] * w_j_bar[1] +
                                w_j_bar[2] * w_j_bar[2];
  const double w_j_bar_norm = sqrt(w_j_bar_norm_2);

  /* If the particle does not contribute, skip the computations. This
   * avoids 1./0. */
  if (w_j_bar_norm == 0) {
    return;
  }

  /* Do we have supernovae? */
  if (feedback_should_inject_SN_feedback(si)) {

    /****************************************************************************
     * Compute the common properties for both modes
     ****************************************************************************/
    /* Here just get the feedback properties we want to distribute (in physical
       units) */
    const float E_ej = si->feedback_data.energy_ejected;
    const float mj = hydro_get_mass(pj);
    const float m_ej = si->feedback_data.mass_ejected;

    /* Distribute mass... (the max avoids to have dm=0 and 1/0 divisions) */
    const double dm = max(w_j_bar_norm * m_ej, FLT_MIN);
    const double new_mass = mj + dm;
    xpj->feedback_data.delta_mass += dm;

    /* ... metals */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      pj->chemistry_data.metal_mass[i] +=
          w_j_bar_norm * si->feedback_data.metal_mass_ejected[i];
    }

    /* Calculate the velocity with the Hubble flow */
    const float a = cosmo->a;
    const float H = cosmo->H;
    const float a2H = a * a * H;
    const float vi_plus_H_flow[3] = {a2H * si->x[0] + si->v[0],
                                     a2H * si->x[1] + si->v[1],
                                     a2H * si->x[2] + si->v[2]};
    const float vj_plus_H_flow[3] = {a2H * pj->x[0] + xpj->v_full[0],
                                     a2H * pj->x[1] + xpj->v_full[1],
                                     a2H * pj->x[2] + xpj->v_full[2]};

    /* Compute the _physical_ relative velocity between the particles */
    const float v_i_p[3] = {vi_plus_H_flow[0] * cosmo->a_inv,
                            vi_plus_H_flow[1] * cosmo->a_inv,
                            vi_plus_H_flow[2] * cosmo->a_inv};

    const float v_j_p[3] = {vj_plus_H_flow[0] * cosmo->a_inv,
                            vj_plus_H_flow[1] * cosmo->a_inv,
                            vj_plus_H_flow[2] * cosmo->a_inv};

    /****************************************************************************
     * Now we treat the fluxes distribution differently for each mode
     ****************************************************************************/
    double dU = 0.0;
    double dKE = 0.0;
    double dp_prime[3] = {0.0, 0.0, 0.0};

    runner_iact_nonsym_mechanical_feedback_apply(
        r2, si, pj, xpj, w_j_bar, w_j_bar_norm, v_i_p, v_j_p, E_ej, m_ej, mj,
        dm, new_mass, cosmo, fb_props, phys_const, us, &dU, &dKE, dp_prime);

    /* Now we can give momentum, thermal and kinetic energy to the xpart.
       Note: Do not give momentum for the isotropy check test. Momentum pushes
       particles too efficiently and then the python face area computations are
       not exacly the same as SWIFT. */
#if !defined(SWIFT_TEST_FEEDBACK_ISOTROPY_CHECK)
    /* Convert to comoving units */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.delta_p[i] += dp_prime[i] * cosmo->a;
    }
#endif /* !defined SWIFT_TEST_FEEDBACK_ISOTROPY_CHECK */

    /* Note: This is physical internal energy. See feedback_update_part(). */
    xpj->feedback_data.delta_u += dU / new_mass;
    xpj->feedback_data.delta_E_kin += dKE;
    xpj->feedback_data.number_SN += 1;

    /* Only use this for non-cosmological simulations. In cosmological
       simulations, this suppresses the momentum effects (to be investigated
       why). */
    if (cosmo->a == 1.0 && cosmo->a_inv == 1.0 && cosmo->z == 0.0) {
      /* Update the signal velocity of gas particles that receive a kick. From
         Chaikin et al. (2023) (also implemented in EAGLE_kinetic) */
      const double dp_prime_norm_2 = dp_prime[0] * dp_prime[0] +
                                     dp_prime[1] * dp_prime[1] +
                                     dp_prime[2] * dp_prime[2];
      const float dp_prime_norm = sqrt(dp_prime_norm_2);
      const float dv_phys = dp_prime_norm / new_mass;
      hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, dv_phys);
    }
  }

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}

#endif /* SWIFT_GEAR_MECHANICAL_FEEDBACK_IACT_H */
