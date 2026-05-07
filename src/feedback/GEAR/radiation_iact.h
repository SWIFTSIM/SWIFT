/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_IACT_GEAR_H
#define SWIFT_RADIATION_IACT_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_iact.h
 * @brief Subgrid radiation feedback for GEAR. This file contains the generic
 * functions to be called in feedback_iact.h or
 * feedback_prepare_feedback(). The radiation model is split into this
 * functions so that they can be called by the mechanical feedback without code
 * duplication.
 */

#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "feedback_properties.h"
#include "radiation.h"
#include "random.h"

/**
 * @brief Radiation density interaction between two particles (non-symmetric).
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
radiation_iact_nonsym_feedback_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, const struct part *pj, const struct xpart *xpj,
    const struct cosmology *cosmo, const struct feedback_props *fb_props,
    const struct hydro_props *hydro_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current) {

  /* Exit if no radiation policies are enabled */
  if (fb_props->radiation_policy == radiation_policy_none) {
    return;
  }

  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(ui, &wi, &wi_dx);

  /* Unit vector pointing to pj */
  float dx_unit[3];
  for (int k = 0; k < 3; ++k) {
    dx_unit[k] = dx[k] / r;
  }

  /* Gradient of the kernel */
  float gradW[3];
  for (int k = 0; k < 3; ++k) {
    gradW[k] = wi_dx * dx_unit[k];
  }

  /* Gradient of the density */
  for (int k = 0; k < 3; ++k) {
    si->feedback_data.grad_rho_star[k] += mj * gradW[k];
  }

  /* Metallicity at the star location */
  si->feedback_data.Z_star += chemistry_get_total_metal_mass_fraction_for_feedback(pj) * wi;
}

/**
 * @brief Prepare a #spart for the radiation feedback task. Here we perform the
 * photoionization of HII regions.
 *
 * This is called in the feedback_prepare_feedback(), which is called in the
 * stars ghost task.
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
__attribute__((always_inline)) INLINE static void
feedback_prepare_radiation_feedback(
    struct spart *restrict sp, const struct feedback_props *feedback_props,
    const struct cosmology *cosmo, const struct unit_system *us,
    const struct phys_const *phys_const, const double star_age_beg_step,
    const double dt, const double time, const integertime_t ti_begin,
    const int with_cosmology) {
  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */

  sp->feedback_data.grad_rho_star[0] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[1] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[2] *= hi_inv_dim_plus_one;

  sp->feedback_data.Z_star *= hi_inv / sp->feedback_data.enrichment_weight;

  /* const float Sigma_gas =
   * radiation_get_comoving_gas_column_density_at_star(sp); */
  /* const float kappa_IR = radiation_get_physical_IR_opacity(sp, us,
   * phys_const, cosmo); */
  /* const float tau_IR = radiation_get_physical_IR_optical_depth(sp, us,
   * phys_const, cosmo); */

  /* message( */
  /*     "rho_star = %e, Grad rho = (%e %e %e), Z = %e, Sigma_gas = %e, kappa_IR
   * " */
  /*     "= %e, tau_IR = %e", */
  /*     sp->feedback_data.rho_star, sp->feedback_data.grad_rho_star[0], */
  /*     sp->feedback_data.grad_rho_star[1], sp->feedback_data.grad_rho_star[2],
   */
  /*     sp->feedback_data.Z_star, Sigma_gas, kappa_IR, tau_IR); */
}

/**
 * @brief Radiation feedback interaction between two particles
 * (non-symmetric). Used for updating properties of gas particles neighbouring
 * a star particle.
 *
 * Here we tag particles within HII regions and apply radiation pressure.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
radiation_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current, const double time_base) {

  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);

  /* Get the kernel for hi. */
  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;

  /* Compute inverse enrichment weight */
  const double si_inv_weight = si->feedback_data.enrichment_weight == 0
                                   ? 0.
                                   : 1. / si->feedback_data.enrichment_weight;
  const double weight = mj * wi * si_inv_weight;

  /* Compute radiation pressure */
  if (si->feedback_data.radiation.L_bol != 0.0) {
    const float Delta_t = get_timestep(si->time_bin, time_base);
    const float p_rad = radiation_get_star_physical_radiation_pressure(
        si, Delta_t, phys_const, us, cosmo);
    const float delta_p_rad = weight * p_rad;

    /* Add the radiation pressure radially outwards from the star. Notice the
       conversion to comoving units. */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.radiation.delta_p[i] -=
          delta_p_rad * dx[i] / r * cosmo->a;
    }
  }

  /*
     5. Transport the emergent FUV radiation. And then compute the
     photohelectric heating. We assume that the effect is only local and so we
     do not transport radiation.
  */
}

/**
 * @brief Update the properties of the particle due to radiation feedback.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_update_part_radiation(struct part *p, struct xpart *xp,
                               const struct engine *e,
                               const float initial_mass) {

  /* Here, wo only update radiation pressure. The gas cooling state is updated
     before cooling */
  if (xp->feedback_data.hit_by_radiation) {
    for (int i = 0; i < 3; i++) {
      /* We use the initial mass of the gas, i.e. before any winds or SN */
      const float dv = xp->feedback_data.radiation.delta_p[i] / initial_mass;
      xp->v_full[i] += dv;
      p->v[i] += dv;

      /* Reset */
      xp->feedback_data.radiation.delta_p[i] = 0;
    }
    xp->feedback_data.hit_by_radiation = 0;
  }
}

#endif /* SWIFT_RADIATION_IACT_GEAR_H */
