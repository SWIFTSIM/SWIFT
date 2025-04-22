/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_FEEDBACK_IACT_H
#define SWIFT_GEAR_FEEDBACK_IACT_H

/* Local includes */
#include "error.h"
#include "feedback.h"
#include "hydro.h"
#include "random.h"
#include "timestep_sync_part.h"
#include "radiation.h"
#include "random.h"

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
				    const struct hydro_props *hydro_props,
				    const struct phys_const* phys_const,
				    const struct unit_system* us,
				    const struct cooling_function_data* cooling,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(ui, &wi, &wi_dx);

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  /* The normalization by 1 / h^d is done in feedback.h */
  si->feedback_data.enrichment_weight += mj * wi;

  /* Contribution to the number of neighbours */
  si->feedback_data.num_ngbs += 1;

  /* Radiation */
  /* Compute the column density with the sobolev approx: here we need to compute rho
     and | grad rho | at the star location using the gas particles. */
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
  si->feedback_data.Z_star += pj->chemistry_data.metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT-1] * wi;

  /* Gather neighbours data for HII ionization */
  if (!radiation_is_part_ionized(phys_const, hydro_props, us, cosmo, cooling, pj, xpj)) {
    /* If a particle is already ionized, it won't be able to ionize again so do
       not gather its data. */
    const double Delta_dot_N_ion = radiation_get_part_rate_to_fully_ionize(phys_const, hydro_props, us, cosmo, cooling, pj, xpj);

    /* Compute the size of the array that we want to sort. If the current
     * function is called for the first time (at this time-step for this star),
     * then si->num_ngbs = 1 and there is nothing to sort. Note that the
     * maximum size of the sorted array cannot be larger then the maximum
     * number of rays. */
    const int arr_size = min(si->feedback_data.num_ngbs, GEAR_STROMGREN_NUMBER_NEIGHBOURS);

    /* Minimise separation between the gas particles and the BH. The rays
     * structs with smaller ids in the ray array will refer to the particles
     * with smaller distances to the BH. */
    stromgren_sort_distance(r, si->feedback_data.radiation.stromgren_sphere, arr_size, Delta_dot_N_ion);
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
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props,  const struct phys_const* phys_const,
    const struct unit_system* us, const struct cooling_function_data* cooling,
    const integertime_t ti_current, const double time_base) {

  const double e_sn = si->feedback_data.energy_ejected;

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

  /* Mass received */
  const double m_ej = si->feedback_data.mass_ejected;
  const double weight = mj * wi * si_inv_weight;
  const double dm = m_ej * weight;
  const double new_mass = mj + dm;

  if (e_sn != 0.0) {
    /* Energy received */
    const double du = e_sn * weight / new_mass;

    xpj->feedback_data.delta_mass += dm;
    xpj->feedback_data.delta_u += du;

    /* Compute momentum received. */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.delta_p[i] += dm * (si->v[i] - xpj->v_full[i]);
    }

    /* Add the metals */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      pj->chemistry_data.metal_mass[i] +=
        weight * si->feedback_data.metal_mass_ejected[i];
    }
  }

  /* TODO: Distribute pre-SN */
  /* 1. Here we know
     - Get Column density Sigma: eq E3
     - Get tau_nu = kappa * Simga (eq E1)
     - The emitted luminosity in each band */

  /* 2. Local extinction around the star: L_abs_nu = (1 - exp(-tau_nu)) L_nu
     L_emitted_nu = exp(-tau_nu) L_nu
     For the IR band, we have:
     L_IR = Sum_{nu=FUR, UV, Opt} L_abs_nu

     For the ioninzing band, we do not do the local extinction. We treat it
     separately.

     3. Photoionization - HII region:
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

  /* Here we must choose the model according to the metallicity to distinguish
     pop III fomr pop II */
  /* const struct radiation* radiation = &fb_props->stellar_model.radiation; */

  /* 3. Photoionization */
  const float R_stromgren = si->feedback_data.radiation.R_stromgren;
  if (r <= R_stromgren) {
    message("Found particle to ionize ! r= %e, id = %lld", r, pj->id);
    /* Tag the particle */
    radiation_tag_part_as_ionized(pj, xpj);
  }

  /* 4. Compute radiation pressure */
  /* TODO: DO we want to compute it here or at the same locations than SN
     feedback? */
  if (fb_props->radiation_pressure_efficiency != 0) {
    const float Delta_t = get_timestep(si->time_bin, time_base);
    const float p_rad = fb_props->radiation_pressure_efficiency * radiation_get_star_radiation_pressure(si, Delta_t, us, phys_const);
    const float delta_p_rad = weight * p_rad;

    /* Add the radiation pressure radially outwards from the star */
    for (int i = 0; i < 3; i++) {
      xpj->feedback_data.radiation.delta_p[i] -= delta_p_rad * dx[i] / r;
    }
  }

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}

#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
