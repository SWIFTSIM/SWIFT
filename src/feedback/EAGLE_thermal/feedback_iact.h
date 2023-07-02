/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_FEEDBACK_IACT_THERMAL_H
#define SWIFT_EAGLE_FEEDBACK_IACT_THERMAL_H

/* Local includes */
#include "feedback.h"
#include "random.h"
#include "rays.h"
#include "timestep_sync_part.h"
#include "tracers.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
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

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* We found a neighbour! */
  si->feedback_data.to_collect.ngb_N++;

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.to_collect.ngb_mass += mj;

  /* Update counter of total (integer) neighbours */
  si->feedback_data.to_collect.num_ngbs++;

  /* Contribution to the star's surrounding gas density */
  si->feedback_data.to_collect.ngb_rho += mj * wi;

  const float Zj = chemistry_get_total_metal_mass_fraction_for_feedback(pj);

  /* Contribution to the star's surrounding metallicity (metal mass fraction */
  si->feedback_data.to_collect.ngb_Z += mj * Zj * wi;

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */
  const float rho = hydro_get_comoving_density(pj);
  if (rho != 0.f)
    si->feedback_data.to_collect.enrichment_weight_inv += wi / rho;

  /* Choose SNII feedback model */
  switch (fb_props->feedback_model) {
    case SNII_isotropic_model: {

      /* Compute arc lengths in stellar isotropic feedback and collect
       * relevant data for later use in the feedback_apply loop */

      /* Loop over rays */
      for (int i = 0; i < eagle_SNII_feedback_num_of_rays; i++) {

        /* We generate two random numbers that we use
         * to randomly select the direction of the ith ray */

        /* Two random numbers in [0, 1[ */
        const double rand_theta_SNII = random_unit_interval_part_ID_and_index(
            si->id, i, ti_current,
            random_number_isotropic_SNII_feedback_ray_theta);
        const double rand_phi_SNII = random_unit_interval_part_ID_and_index(
            si->id, i, ti_current,
            random_number_isotropic_SNII_feedback_ray_phi);

        /* Compute arclength */
        ray_minimise_arclength(dx, r, si->feedback_data.SNII_rays + i,
                               /*ray_type=*/ray_feedback_thermal, pj->id,
                               rand_theta_SNII, rand_phi_SNII, mj,
                               /*ray_ext=*/NULL, /*v=*/NULL);
      }
      break;
    }
    case SNII_minimum_distance_model: {
      /* Compute the size of the array that we want to sort. If the current
       * function is called for the first time (at this time-step for this
       * star), then bi->num_ngbs = 1 and there is nothing to sort. Note that
       * the maximum size of the sorted array cannot be larger then the maximum
       * number of rays. */
      const int arr_size = min(si->feedback_data.to_collect.ngb_N,
                               eagle_SNII_feedback_num_of_rays);

      /* Minimise separation between the gas particles and the star. The rays
       * structs with smaller ids in the ray array will refer to the particles
       * with smaller distances to the star. */
      ray_minimise_distance(r, si->feedback_data.SNII_rays, arr_size, pj->id,
                            mj);
      break;
    }
    case SNII_minimum_density_model: {
      /* Compute the size of the array that we want to sort. If the current
       * function is called for the first time (at this time-step for this
       * star), then bi->num_ngbs = 1 and there is nothing to sort. Note that
       * the maximum size of the sorted array cannot be larger then the maximum
       * number of rays. */
      const int arr_size = min(si->feedback_data.to_collect.ngb_N,
                               eagle_SNII_feedback_num_of_rays);

      /* Minimise separation between the gas particles and the star. The rays
       * structs with smaller ids in the ray array will refer to the particles
       * with smaller distances to the star. */
      ray_minimise_distance(rho, si->feedback_data.SNII_rays, arr_size, pj->id,
                            mj);
      break;
    }
    case SNII_random_ngb_model: {
      /* Compute the size of the array that we want to sort. If the current
       * function is called for the first time (at this time-step for this
       * star), then bi->num_ngbs = 1 and there is nothing to sort. Note that
       * the maximum size of the sorted array cannot be larger then the maximum
       * number of rays. */
      const int arr_size = min(si->feedback_data.to_collect.ngb_N,
                               eagle_SNII_feedback_num_of_rays);

      /* To mimic a random draw among all the particles in the kernel, we
       * draw random distances in [0,1) and then pick the particle(s) with
       * the smallest of these 'fake' distances */
      const float dist = random_unit_interval_two_IDs(
          si->id, pj->id, ti_current, random_number_stellar_feedback_1);

      /* Minimise separation between the gas particles and the BH. The rays
       * structs with smaller ids in the ray array will refer to the particles
       * with smaller 'fake' distances to the BH. */
      ray_minimise_distance(dist, si->feedback_data.SNII_rays, arr_size, pj->id,
                            mj);
      break;
    }
  }
}

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
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const integertime_t ti_current) {

#ifdef SWIFT_DEBUG_CHECKS
  if (si->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback from a star that should not");
#endif

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Gas particle density */
  const float rho_j = hydro_get_comoving_density(pj);

  /* Compute weighting for distributing feedback quantities */
  float Omega_frac;
  if (rho_j != 0.f) {
    Omega_frac = si->feedback_data.to_distribute.enrichment_weight * wi / rho_j;
  } else {
    Omega_frac = 0.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (Omega_frac < 0. || Omega_frac > 1.01)
    error(
        "Invalid fraction of material to distribute for star ID=%lld "
        "Omega_frac=%e count since last enrich=%d",
        si->id, Omega_frac, si->count_since_last_enrichment);
#endif

  /* Update particle mass */
  const double current_mass = hydro_get_mass(pj);
  const double delta_mass = si->feedback_data.to_distribute.mass * Omega_frac;
  const double new_mass = current_mass + delta_mass;

  hydro_set_mass(pj, new_mass);

  /* Inverse of the new mass */
  const double new_mass_inv = 1. / new_mass;

  /* Update total metallicity */
  const double current_metal_mass_total =
      pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const double delta_metal_mass_total =
      si->feedback_data.to_distribute.total_metal_mass * Omega_frac;
  const double new_metal_mass_total =
      current_metal_mass_total + delta_metal_mass_total;

  pj->chemistry_data.metal_mass_fraction_total =
      new_metal_mass_total * new_mass_inv;

  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const double delta_metal_mass =
        si->feedback_data.to_distribute.metal_mass[elem] * Omega_frac;
    const double new_metal_mass = current_metal_mass + delta_metal_mass;

    pj->chemistry_data.metal_mass_fraction[elem] =
        new_metal_mass * new_mass_inv;
  }

  /* Update iron mass fraction from SNIa  */
  const double current_iron_from_SNIa_mass =
      pj->chemistry_data.iron_mass_fraction_from_SNIa * current_mass;
  const double delta_iron_from_SNIa_mass =
      si->feedback_data.to_distribute.Fe_mass_from_SNIa * Omega_frac;
  const double new_iron_from_SNIa_mass =
      current_iron_from_SNIa_mass + delta_iron_from_SNIa_mass;

  pj->chemistry_data.iron_mass_fraction_from_SNIa =
      new_iron_from_SNIa_mass * new_mass_inv;

  /* Update mass from SNIa  */
  const double delta_mass_from_SNIa =
      si->feedback_data.to_distribute.mass_from_SNIa * Omega_frac;

  pj->chemistry_data.mass_from_SNIa += delta_mass_from_SNIa;

  /* Update metal mass fraction from SNIa */
  const double current_metal_mass_from_SNIa =
      pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  const double delta_metal_mass_from_SNIa =
      si->feedback_data.to_distribute.metal_mass_from_SNIa * Omega_frac;
  const double new_metal_mass_from_SNIa =
      current_metal_mass_from_SNIa + delta_metal_mass_from_SNIa;

  pj->chemistry_data.metal_mass_fraction_from_SNIa =
      new_metal_mass_from_SNIa * new_mass_inv;

  /* Update mass from SNII  */
  const double delta_mass_from_SNII =
      si->feedback_data.to_distribute.mass_from_SNII * Omega_frac;

  pj->chemistry_data.mass_from_SNII += delta_mass_from_SNII;

  /* Update metal mass fraction from SNII */
  const double current_metal_mass_from_SNII =
      pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  const double delta_metal_mass_from_SNII =
      si->feedback_data.to_distribute.metal_mass_from_SNII * Omega_frac;
  const double new_metal_mass_from_SNII =
      current_metal_mass_from_SNII + delta_metal_mass_from_SNII;

  pj->chemistry_data.metal_mass_fraction_from_SNII =
      new_metal_mass_from_SNII * new_mass_inv;

  /* Update mass from AGB  */
  const double delta_mass_from_AGB =
      si->feedback_data.to_distribute.mass_from_AGB * Omega_frac;

  pj->chemistry_data.mass_from_AGB += delta_mass_from_AGB;

  /* Update metal mass fraction from AGB */
  const double current_metal_mass_from_AGB =
      pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  const double delta_metal_mass_from_AGB =
      si->feedback_data.to_distribute.metal_mass_from_AGB * Omega_frac;
  const double new_metal_mass_from_AGB =
      current_metal_mass_from_AGB + delta_metal_mass_from_AGB;

  pj->chemistry_data.metal_mass_fraction_from_AGB =
      new_metal_mass_from_AGB * new_mass_inv;

  /* Compute the current kinetic energy */
  const double current_v2 = xpj->v_full[0] * xpj->v_full[0] +
                            xpj->v_full[1] * xpj->v_full[1] +
                            xpj->v_full[2] * xpj->v_full[2];
  const double current_kinetic_energy_gas =
      0.5 * cosmo->a2_inv * current_mass * current_v2;

  /* Compute the current thermal energy */
  const double current_thermal_energy =
      current_mass * hydro_get_physical_internal_energy(pj, xpj, cosmo);

  /* Apply conservation of momentum */

  /* Update velocity following change in gas mass */
  xpj->v_full[0] *= current_mass * new_mass_inv;
  xpj->v_full[1] *= current_mass * new_mass_inv;
  xpj->v_full[2] *= current_mass * new_mass_inv;

  /* Update velocity following addition of mass with different momentum */
  xpj->v_full[0] += delta_mass * new_mass_inv * si->v[0];
  xpj->v_full[1] += delta_mass * new_mass_inv * si->v[1];
  xpj->v_full[2] += delta_mass * new_mass_inv * si->v[2];

  /* Compute the new kinetic energy */
  const double new_v2 = xpj->v_full[0] * xpj->v_full[0] +
                        xpj->v_full[1] * xpj->v_full[1] +
                        xpj->v_full[2] * xpj->v_full[2];
  const double new_kinetic_energy_gas = 0.5 * cosmo->a2_inv * new_mass * new_v2;

  /* Energy injected
   * (thermal SNIa + kinetic energy of ejecta + kinetic energy of star) */
  const double injected_energy =
      si->feedback_data.to_distribute.energy * Omega_frac;

  /* Apply energy conservation to recover the new thermal energy of the gas
   *
   * Note: in some specific cases the new_thermal_energy could be lower
   * than the current_thermal_energy, this is mainly the case if the change
   * in mass is relatively small and the velocity vectors between both the
   * gas particle and the star particle have a small angle. */
  double new_thermal_energy = current_kinetic_energy_gas +
                              current_thermal_energy + injected_energy -
                              new_kinetic_energy_gas;

  /* In rare configurations the new thermal energy could become negative.
   * We must prevent that even if that implies a slight violation of the
   * conservation of total energy.
   * The minimum energy (in units of energy not energy per mass) is
   * the total particle mass (including the mass to distribute) at the
   * minimal internal energy per unit mass */
  const double min_u = hydro_props->minimal_internal_energy * new_mass;

  new_thermal_energy = max(new_thermal_energy, min_u);

  /* Convert this to a specific thermal energy */
  const double u_new_enrich = new_thermal_energy * new_mass_inv;

  /* Do the energy injection. */
  hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new_enrich);
  hydro_set_drifted_physical_internal_energy(pj, cosmo, /*pfloor=*/NULL,
                                             u_new_enrich);

  /* Finally, SNII stochastic feedback */

  /* Get the total number of SNII thermal energy injections per stellar
   * particle at this time-step */
  const int N_of_SNII_thermal_energy_inj =
      si->feedback_data.to_distribute.SNII_num_of_thermal_energy_inj;

  /* Are we doing some SNII feedback? */
  if (N_of_SNII_thermal_energy_inj > 0) {

    int N_of_SNII_energy_inj_received_by_gas = 0;

    /* Find out how many rays this gas particle has received. */
    for (int i = 0; i < N_of_SNII_thermal_energy_inj; i++) {
      if (pj->id == si->feedback_data.SNII_rays[i].id_min_length)
        N_of_SNII_energy_inj_received_by_gas++;
    }

    /* If the number of SNII energy injections > 0, do SNII feedback */
    if (N_of_SNII_energy_inj_received_by_gas > 0) {

      /* Compute new energy of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u = si->feedback_data.to_distribute.SNII_delta_u;
      const double u_new =
          u_init + delta_u * (float)N_of_SNII_energy_inj_received_by_gas;

      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, /*pfloor=*/NULL,
                                                 u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Mark this particle has having been heated by supernova feedback */
      tracers_after_feedback(xpj);

      /* message( */
      /*     "We did some heating! id %llu star id %llu probability %.5e " */
      /*     "random_num %.5e du %.5e du/ini %.5e", */
      /*     pj->id, si->id, 0., 0., delta_u, delta_u / u_init); */

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }
}

#endif /* SWIFT_EAGLE_FEEDBACK_IACT_THERMAL_H */
