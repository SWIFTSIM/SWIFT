/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_KIARA_FEEDBACK_IACT_H
#define SWIFT_KIARA_FEEDBACK_IACT_H

/* Local includes */
#include "random.h"
#include "rays.h"
#include "timestep_sync_part.h"
#include "tools.h"
#include "tracers.h"

/**
 * @brief Compute the mean DM velocity around a star. (non-symmetric).
 *
 * @param si First sparticle.
 * @param gj Second particle (not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_dm_vel_sum(struct spart *si, const struct gpart *gj,
                                       int *dm_ngb_N,
                                       float dm_mean_velocity[3]) {}

/**
 * @brief Compute the DM velocity dispersion around a star. (non-symmetric).
 *
 * @param si First sparticle.
 * @param gj Second particle.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_dm_vel_disp(struct spart *si,
                                        const struct gpart *gj,
                                        const float dm_mean_velocity[3]) {}

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

  /* Ignore wind in density computation */
  /*if (pj->feedback_data.decoupling_delay_time > 0.f) return;*/

  const float rho = hydro_get_comoving_density(pj);
  if (rho <= 0.f) return;

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
  si->feedback_data.ngb_N++;

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.ngb_mass += mj;

  /* Update counter of total (integer) neighbours */
  si->feedback_data.num_ngbs++;

  /* Contribution to the star's surrounding gas density */
  si->feedback_data.ngb_rho += mj * wi;

  const float Zj = chemistry_get_total_metal_mass_fraction_for_feedback(pj);

  /* Contribution to the star's surrounding metallicity (metal mass fraction */
  si->feedback_data.ngb_Z += mj * Zj * wi;

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */
  si->feedback_data.enrichment_weight_inv += wi / rho;

}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  const struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {}

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
    const struct feedback_props *fb_props, 
    const integertime_t ti_current) {

  /* Ignore decoupled particles */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  /* Gas particle density */
  const float rho_j = hydro_get_comoving_density(pj);
  if (rho_j <= 0.f) return;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Compute weighting for distributing feedback quantities */
  const float Omega_frac = si->feedback_data.enrichment_weight * wi / rho_j;

  /* Never apply feedback if Omega_frac is bigger than or equal to unity */
  if (Omega_frac > 1.0) {
    error("Problem with neighbors: Omega_frac=%g wi=%g rho_j=%g",
            Omega_frac, wi, rho_j);
  }

#ifdef SIMBA_DEBUG_CHECKS
  if (Omega_frac < 0. || Omega_frac > 1.01)
    warning(
        "Invalid fraction of material to distribute for star ID=%lld "
        "Omega_frac=%e count since last enrich=%d",
        si->id, Omega_frac, si->count_since_last_enrichment);
#endif

  /* Update particle mass */
  const double current_mass = hydro_get_mass(pj);
  const double delta_mass = si->feedback_data.mass * Omega_frac;
  const double new_mass = current_mass + delta_mass;

  hydro_set_mass(pj, new_mass);

  /* Inverse of the new mass */
  const double new_mass_inv = 1. / new_mass;

  /* Update total metallicity */
  const double current_metal_mass_total =
      pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const double delta_metal_mass_total =
      si->feedback_data.total_metal_mass * Omega_frac;
  const double new_metal_mass_total =
      current_metal_mass_total + delta_metal_mass_total;

  pj->chemistry_data.metal_mass_fraction_total =
      new_metal_mass_total * new_mass_inv;

  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const double delta_metal_mass =
        si->feedback_data.metal_mass[elem] * Omega_frac;
    const double new_metal_mass = current_metal_mass + delta_metal_mass;

    pj->chemistry_data.metal_mass_fraction[elem] =
        new_metal_mass * new_mass_inv;
  }

  /* Compute the current thermal energy */
  const double current_thermal_energy =
      current_mass * hydro_get_physical_internal_energy(pj, xpj, cosmo);

  /* Energy injected */
  const double injected_energy =
      si->feedback_data.energy * Omega_frac;

  double new_thermal_energy = current_thermal_energy + injected_energy;

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
  hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new_enrich);

#if COOLING_GRACKLE_MODE >= 2
  if (0) {
    /* NOT USED: Compute G0 contribution from star to the gas particle in Habing units of 
     * 1.6e-3 erg/s/cm^2. Note that this value includes the 4*pi geometric factor */
    const float length_to_physical_cm = cosmo->a * fb_props->length_to_kpc * 3.08567758e21f;
    /* Compute a softened distance from star to gas particle */
    const double r2_in_cm = (r2 + 0.01*hi*hi) * length_to_physical_cm * length_to_physical_cm;
    const double r_in_cm = sqrt(r2_in_cm);
  
    /* Compute self-shielding from H2, from Schauer et al. 2015 eq 8,9 */
    /* H attenuation factor */
    const double NH_cgs = hydro_get_physical_density(pj, cosmo) * fb_props->rho_to_n_cgs * r_in_cm;
    const double xH = NH_cgs / 2.85e23;
    const double fH_shield = pow(1.f+xH,-1.62) * exp(-0.149*xH);
    /* H2 attenuation factor */
    const double NH2_cgs = pj->sf_data.H2_fraction * NH_cgs;
    const double DH2_cgs = 1.e-5 * sqrt(2.*1.38e-16*cooling_get_subgrid_temperature(pj, xpj) / 3.346e-24);
    const double xH2 = NH2_cgs / 8.465e13;
    const double fH2_shield = 0.9379/pow(1.f+xH2/DH2_cgs,1.879) + 0.03465/pow(1.f+xH2,0.473) * exp(-2.293e-4*sqrt(1+xH2));
    //message("G0 shield: r=%g xH2=%g xH=%g fH2=%g fH=%g\n",r_in_cm/3.086e21,xH2,xH,fH2_shield,fH_shield);
  
    if (si->feedback_data.lum_habing > -10.) {  
      pj->chemistry_data.G0 += fH2_shield * fH_shield * pow(10.,si->feedback_data.lum_habing) / (1.6e-3 * r2_in_cm);
    }
  }

  /* Compute kernel-smoothed contribution to number of SNe going off this timestep */
  pj->feedback_data.SNe_ThisTimeStep += si->feedback_data.SNe_ThisTimeStep * Omega_frac;

  /* Spread dust ejecta to gas */
  xpj->cooling_data.dust_mass = 0.f;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_dust_mass =
        xpj->cooling_data.dust_mass_fraction[elem] * current_mass;
    const double delta_dust_mass =
        si->feedback_data.delta_dust_mass[elem] * Omega_frac;

    xpj->cooling_data.dust_mass_fraction[elem] =
        (current_dust_mass + delta_dust_mass) * new_mass_inv;
    /* Sum up each element to get total dust mass */
    if (elem > chemistry_element_C) xpj->cooling_data.dust_mass += current_dust_mass + delta_dust_mass;
  }
  if (xpj->cooling_data.dust_mass > pj->mass) {
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      message("DUST EXCEEDS MASS elem=%d md=%g\n",elem, xpj->cooling_data.dust_mass_fraction[elem]);
    }
    error("DUST EXCEEDS MASS mgas=%g  mdust=%g\n",pj->mass, xpj->cooling_data.dust_mass);
  }

#endif

}

#endif /* SWIFT_KIARA_FEEDBACK_IACT_H */
