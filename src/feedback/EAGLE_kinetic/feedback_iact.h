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
#ifndef SWIFT_EAGLE_FEEDBACK_KINETIC_IACT_H
#define SWIFT_EAGLE_FEEDBACK_KINETIC_IACT_H

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

  /* Compute arc lengths in stellar isotropic feedback and collect
   * relevant data for later use in the feedback_apply loop */
  for (int i = 0; i < eagle_SNII_feedback_num_of_rays; i++) {

    /* We generate two random numbers that we use
     * to randomly select the direction of the ith ray
     * and the associated mirror ray in SNII feedback */

    /* Two random numbers in [0, 1[ */
    const double rand_theta_SNII = random_unit_interval_part_ID_and_index(
        si->id, i, ti_current, random_number_isotropic_SNII_feedback_ray_theta);
    const double rand_phi_SNII = random_unit_interval_part_ID_and_index(
        si->id, i, ti_current, random_number_isotropic_SNII_feedback_ray_phi);

    /* Compute arclength for the true particle (SNII kinetic feedback) */
    ray_minimise_arclength(dx, r, si->feedback_data.SNII_rays_true + i,
                           ray_feedback_kinetic_true, pj->id, rand_theta_SNII,
                           rand_phi_SNII, mj,
                           si->feedback_data.SNII_rays_ext_true + i, pj->v);
    /* Compute arclength for the mirror particle (SNII kinetic feedback) */
    ray_minimise_arclength(dx, r, si->feedback_data.SNII_rays_mirr + i,
                           ray_feedback_kinetic_mirr, pj->id, rand_theta_SNII,
                           rand_phi_SNII, mj,
                           si->feedback_data.SNII_rays_ext_mirr + i, pj->v);
  }
}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  const struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {

  /* Get the the number of SNII kinetic energy injections per stellar
   * particle at this time-step */
  const int N_of_SNII_kinetic_events =
      si->feedback_data.to_distribute.SNII_num_of_kinetic_energy_inj;

  /* Loop over the SNII kick events. In each event, two gas
   * particles are kicked in exactly the opposite directions. */
  for (int i = 0; i < N_of_SNII_kinetic_events; i++) {

    /* Find the particle that is closest to the ith ray OR the ith mirror ray
     */
    if (pj->id == si->feedback_data.SNII_rays_true[i].id_min_length ||
        pj->id == si->feedback_data.SNII_rays_mirr[i].id_min_length) {

      /* If this spart has the largest id among all sparts that want to kick
       * this gas particle in this time-step, then the gas particle will save
       * the id of this spart. */
      if (pj->feedback_data.SNII_star_largest_id < si->id) {
        /* Update the largest stellar id carried by the gas particle */
        pj->feedback_data.SNII_star_largest_id = si->id;
      }
    }
  }
}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {

  /* Get the the number of SNII kinetic energy injections per stellar
   * particle at this time-step */
  const int N_of_SNII_kinetic_events =
      si->feedback_data.to_distribute.SNII_num_of_kinetic_energy_inj;

  for (int i = 0; i < N_of_SNII_kinetic_events; i++) {

    /* Find the particle that is closest to the ith ray */
    if (pj->id == si->feedback_data.SNII_rays_true[i].id_min_length) {

      /* Does this gas particle want to be kicked by this stellar particle
       * via ray i? If so, store this information in the ith ray extra struct */
      if (pj->feedback_data.SNII_star_largest_id == si->id) {

        si->feedback_data.SNII_rays_ext_true[i].status =
            ray_feedback_kick_allowed;

        /* If we are using maximum_number_of_rays > 1, then for a given spart,
         * as soon as we have found the first ray that points at this gas part,
         * we stop. Otherwise, the same spart might kick the same gas part
         * twice in the same time-step (or even more times). */
        break;
      }
    }
    /* Same as above but for the mirror ith ray */
    else if (pj->id == si->feedback_data.SNII_rays_mirr[i].id_min_length) {
      if (pj->feedback_data.SNII_star_largest_id == si->id) {
        si->feedback_data.SNII_rays_ext_mirr[i].status =
            ray_feedback_kick_allowed;
        break;
      }
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

  /* SNII stochastic kinetic feedback begins.
   *
   * To conserve linear momentum, it is done before the particle velocity
   * is recomputed due to the change in particle mass */

  /* Get the the number of SNII kinetic energy injections from this star
   * particle at this time-step */
  const int N_of_SNII_kinetic_events =
      si->feedback_data.to_distribute.SNII_num_of_kinetic_energy_inj;

  double E_kinetic_unused = 0.0;

  /* Are we doing some SNII kinetic feedback? */
  if (N_of_SNII_kinetic_events > 0) {

    /* Loop over the number of SNII kinetic events. In each event, two
     * particles are kicked in exactly the opposite directions. */
    for (int i = 0; i < N_of_SNII_kinetic_events; i++) {

      /* Find whether we are are looking at the particle that is closest to the
       * ith ray OR the ith mirror ray */
      if (pj->id == si->feedback_data.SNII_rays_true[i].id_min_length ||
          pj->id == si->feedback_data.SNII_rays_mirr[i].id_min_length) {

        /* Get the SNII feedback kinetic energy per pair.
         * We thus divide the  total kinetic energy we have from the star
         * particle by the number of events (in each event, two particles are
         * kicked) */
        const double energy_per_pair =
            si->feedback_data.to_distribute.SNII_E_kinetic /
            N_of_SNII_kinetic_events;

        /* Are we kicking or heating? If at least one gas part in the pair does
         * not want to be kicked by this spart, then we heat */
        if (si->feedback_data.SNII_rays_ext_true[i].status ==
                ray_feedback_kick_allowed &&
            si->feedback_data.SNII_rays_ext_mirr[i].status ==
                ray_feedback_kick_allowed) {

          /* Which particles have we caught: the original or the mirror one? */
          const ray_feedback_type ray_type =
              (pj->id == si->feedback_data.SNII_rays_mirr[i].id_min_length)
                  ? ray_feedback_kinetic_mirr
                  : ray_feedback_kinetic_true;

          /* Two random numbers in [0, 1[
           * Note: this are the same numbers we drew in the density loop! */
          const double rand_theta = random_unit_interval_part_ID_and_index(
              si->id, i, ti_current,
              random_number_isotropic_SNII_feedback_ray_theta);
          const double rand_phi = random_unit_interval_part_ID_and_index(
              si->id, i, ti_current,
              random_number_isotropic_SNII_feedback_ray_phi);

          /* Initialise the kick velocity vector and its modulus */
          float v_kick[3] = {0.f, 0.f, 0.f};
          float v_kick_abs = 0.f;

          /* Get the mass of the gas particles *before* any enrichment mass
           * (from this star or another) was added */
          const double mass_true = si->feedback_data.SNII_rays_true[i].mass;
          const double mass_mirr = si->feedback_data.SNII_rays_mirr[i].mass;

          if (mass_true > 0.0 && mass_mirr > 0.0 && energy_per_pair > 0.0) {

            /* Mass = 0 means that the ray does not point to any gas particle.
             * We need to check it for both the original and mirror ray.
             * We also make sure the energy is positive to avoid the possibility
             * of division by zero in the function below */

            /* Compute the physical kick velocity in internal units */
            ray_kinetic_feedback_compute_kick_velocity(
                v_kick, &v_kick_abs, si->feedback_data.SNII_rays_ext_true + i,
                si->feedback_data.SNII_rays_ext_mirr + i, ray_type,
                energy_per_pair, cosmo, current_mass, si->v, rand_theta,
                rand_phi, mass_true, mass_mirr);
          }

          /* Do the kicks by updating the particle velocity.
           *
           * Note that xpj->v_full = a^2 * dx/dt, with x the comoving
           * coordinate. Therefore, a physical kick, dv, gets translated into a
           * code velocity kick, a * dv */
          xpj->v_full[0] += v_kick[0] * cosmo->a;
          xpj->v_full[1] += v_kick[1] * cosmo->a;
          xpj->v_full[2] += v_kick[2] * cosmo->a;

          /* Update the signal velocity of the particle based on the velocity
           * kick
           */
          hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, v_kick_abs);

          /* Synchronize the particle on the timeline */
          timestep_sync_part(pj);

        } else {

          /* In the absence of a kick event, store the unused kinetic energy
           * for heating */
          E_kinetic_unused = 0.5 * energy_per_pair;
        }
      }
    }
  }

  /* Now account in the fully energy and momentum conserving way for the
   * change in gas particle mass, energy and momentum due to AGB feedback
   * energy and stellar ejecta (with the mass contributed at this time-step
   * by all available feedback channels) moving at the star's velocity */

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
   * which may include extra energy from the failed kinetic injection attempt.
   *
   * Note: in some specific cases the new_thermal_energy could be lower
   * than the current_thermal_energy, this is mainly the case if the change
   * in mass is relatively small and the velocity vectors between both the
   * gas particle and the star particle have a small angle. */
  double new_thermal_energy = current_kinetic_energy_gas +
                              current_thermal_energy + injected_energy +
                              E_kinetic_unused - new_kinetic_energy_gas;

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

  /* Synchronize the particle on the timeline if we've got some extra
   * thermal energy from the SNII kicks that have not occured */
  if (E_kinetic_unused) {
    timestep_sync_part(pj);
  }
}

#endif /* SWIFT_EAGLE_FEEDBACK_KINETIC_IACT_H */
