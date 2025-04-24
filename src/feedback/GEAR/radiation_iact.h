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

#include "feedback.h"
#include "error.h"
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

  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(ui, &wi, &wi_dx);

  /* Gather data to compute the column density with the sobolev approximation.
   */
  si->feedback_data.rho_star += mj * wi;

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
  si->feedback_data.Z_star +=
    pj->chemistry_data.metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT - 1] * wi;

  /* Gather neighbours data for HII ionization */
  if (!radiation_is_part_ionized(phys_const, hydro_props, us, cosmo, cooling,
                                 pj, xpj)) {
    /* If a particle is already ionized, it won't be able to ionize again so do
       not gather its data. */
    const double Delta_dot_N_ion = radiation_get_part_rate_to_fully_ionize(
									   phys_const, hydro_props, us, cosmo, cooling, pj, xpj);

    /* Compute the size of the array that we want to sort. If the current
     * function is called for the first time (at this time-step for this star),
     * then si->num_ngbs = 1 and there is nothing to sort. Note that the
     * maximum size of the sorted array cannot be larger then the maximum
     * number of rays. */
    const int arr_size =
      min(si->feedback_data.num_ngbs, GEAR_STROMGREN_NUMBER_NEIGHBOURS);

    /* Minimise separation between the gas particles and the BH. The rays
     * structs with smaller ids in the ray array will refer to the particles
     * with smaller distances to the BH. */
    stromgren_sort_distance(r, si->feedback_data.radiation.stromgren_sphere,
                            arr_size, Delta_dot_N_ion);
  }
}

/**
 * @brief Prepare a #spart for the radiation feedback task. Here we perform the
 * photoionization of HII refions.
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
__attribute__((always_inline)) INLINE static
void feedback_prepare_radiation_feedback(struct spart* restrict sp,
                               const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo,
                               const struct unit_system* us,
                               const struct phys_const* phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology) {
   /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */

  sp->feedback_data.rho_star *= hi_inv;
  sp->feedback_data.grad_rho_star[0] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[1] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[2] *= hi_inv_dim_plus_one;

  sp->feedback_data.Z_star *= hi_inv / sp->feedback_data.rho_star;

  const float Sigma_gas = radiation_get_comoving_gas_column_density_at_star(sp);
  const float kappa_IR = radiation_get_physical_IR_opacity(sp, us, phys_const, cosmo);
  const float tau_IR = radiation_get_physical_IR_optical_depth(sp, us, phys_const, cosmo);

  message(
      "rho_star = %e, Grad rho = (%e %e %e), Z = %e, Sigma_gas = %e, kappa_IR "
      "= %e, tau_IR = %e",
      sp->feedback_data.rho_star, sp->feedback_data.grad_rho_star[0],
      sp->feedback_data.grad_rho_star[1], sp->feedback_data.grad_rho_star[2],
      sp->feedback_data.Z_star, Sigma_gas, kappa_IR, tau_IR);

  /*----------------------------------------*/
  /* Do the HII ionization */

  if (feedback_props->do_photoionization) {
    const struct stromgren_shell_data* stromgren =
        sp->feedback_data.radiation.stromgren_sphere;
    const int num_ngb =
        min(sp->feedback_data.num_ngbs, GEAR_STROMGREN_NUMBER_NEIGHBOURS);

    /* Loop over the sorted gas neighbours */
    for (int i = 0; i < num_ngb; i++) {
      /* This is recomputed at each iteration */
      const double dot_N_ion = radiation_get_star_ionization_rate(sp);

      if (dot_N_ion <= 0.0) {
        sp->feedback_data.radiation.dot_N_ion = 0.0;
        break;
      }

      const double Delta_dot_N_ion = stromgren[i].Delta_N_dot;

      if (Delta_dot_N_ion <= dot_N_ion) {
        /* We can fully ionize this particle */
        /* Update the Stromgren sphere radius */
        sp->feedback_data.radiation.R_stromgren = stromgren[i].distance;

        /* Consume the photons */
        radiation_consume_ionizing_photons(sp, Delta_dot_N_ion);
      } else {
        /* If we cannot fully ionize, compute a probability to determine if we
           fully ionize pj or not and draw the random number.  */
        const float proba = dot_N_ion / Delta_dot_N_ion;
        const float random_number =
            random_unit_interval(sp->id, ti_begin, random_number_HII_regions);

        /* If we are lucky or we are the first particle, do the ionization */
        if (random_number <= proba || i == 0) {
          /* Update the Stromgren sphere radius */
          sp->feedback_data.radiation.R_stromgren = stromgren[i].distance;
        }

        /* Consume the photons in all cases */
        radiation_consume_ionizing_photons(sp, Delta_dot_N_ion);
      }
    }
  }
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

  /* 3. Photoionization - HII region:
     From step 1 we know N_dot_ion = L_ion / (h nu_ion) . We will use a
     simple stromgren sphere approximation. For each particle:
     a) Test if the particle is already ionized : T > 10^4 or particle was
     flagged to be in an ionized region.
     b) If it is not ionized, compute the ioninzing rate needed to fully
     ionize:
     \Delta N_dot_j = N(H)_j beta n_e_j
     N(H)_j = X_H m_part / (mu m_proton) (the number of H atoms)
     with beta = 3e-13 cm^3 / s is the recombination coefficient, n_e_j the
     electron number density assuming full ionization, X_H is the hydrogen
     mass fraction and mu the molecular weight.
     c) If \Delta N_dot_j <= N_dot_ion:
     tag the particle as being in a HII region
     consume the photons: N_ion -= Delta N_dot_j

     else:
     determine randomly if the particle is ionized by computing the
     proba p = N_dot_io / \Delta N_dot_j
     If rand_number <= proba :
     tag the particle as being in a HII region
     consume the photons: N_ion -= Delta N_dot_j

     d) For the particles tagged as ionized:
     set the temperature (internal energy) to the
     min(current temperature + heat added from the energy of the ionisation,
     equilibrium HII region tem from collisional cooling)
     Set the incident rad and FUV flux to the stromgren value --> to
     compute the inonizing

     Concretely,
     u_new = min(u + delta U, U_collisional),
     delta U = N_H * E_ion / m_gas,
     E_ion = 13.6 eV = 2.18e-11 erg
     Gamma = \Delta N_dot_j / N_H.

     4. Radiation pressure:
     p_rad_tot = Delta t / c * Sum_nu L_abs_nu

     Which bands? Onlyt the IR I guess.

     5. Transport the emergent FUV radiation. And then compute the
     photohelectric heating. We assume that the effect is only local and so we
     do not transport radiation.
  */

  /* Photoionization */
  if (fb_props->do_photoionization) {
    const float R_stromgren = si->feedback_data.radiation.R_stromgren;
    if (r <= R_stromgren) {
      /* message("Found particle to ionize ! r= %e, id = %lld", r, pj->id); */
      /* Tag the particle */
      radiation_tag_part_as_ionized(pj, xpj);
    }
  }

  /* Compute radiation pressure */
  if (fb_props->radiation_pressure_efficiency != 0) {
    const float Delta_t = get_timestep(si->time_bin, time_base);
    const float p_rad =
        fb_props->radiation_pressure_efficiency *
      radiation_get_star_physical_radiation_pressure(si, Delta_t, us, phys_const, cosmo);
    const float delta_p_rad = weight * p_rad;

    /* Add the radiation pressure radially outwards from the star. Notice the
       conversion to comoving units. */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.radiation.delta_p[i] -= delta_p_rad * dx[i] / r * cosmo->a;
    }
  }
}

#endif /* SWIFT_RADIATION_IACT_GEAR_H */
