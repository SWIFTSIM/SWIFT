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
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  /* The normalization by 1 / h^d is done in feedback.h */
  si->feedback_data.enrichment_weight += mj * wi;



  /* Radiation */
  /* Compute the culumn density with the sobolev approx: here we need to compute rho
     and | grad rho | at the star location using the gas particles. */
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
    const integertime_t ti_current) {

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
     Gamma = \Delta N_dot_j / N_H

     Idea: In the original algorithm, Hopkins sort the particles but here
     we do not. Maybe add a proba such that the closest particles have a
     higer chances of being ionized.

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
  /* NOTE: Probalby create a fct to store the computations so that we only have
     to call the functions for the mech fbk */
  /* Is pj already ionized ? If yes, there is nothing to do here. */

  if (!radiation_is_part_ionized(phys_const, hydro_props, us, cosmo, cooling, pj, xpj)) {
    const double dot_N_ion = radiation_get_star_ionisation_rate(si);
    const double Delta_dot_N_ion = radiation_get_part_rate_to_fully_ionize(phys_const, hydro_props, us, cosmo, cooling, pj, xpj);

    /* message("dot_N_ion = %e, Delta_N_dot_ion = %e", dot_N_ion, Delta_dot_N_ion); */

    /* Compute a probability to determine if we fully ionize pj or not and
       draw the random number. */
    /* Note: In the current version we distribute in a weighted manner the
       ionisation */
    const float proba = weight * dot_N_ion / Delta_dot_N_ion;
    const float random_number = random_unit_interval(si->id, ti_current, random_number_HII_regions);
    /* For later when we compute the stromgren radius on the fly
      const int do_ionization = (dot_N_ion <= Delta_dot_N_ion) ? 1 :
      (random_number <= proba);
    */
    const int do_ionization = (random_number <= proba);

    /* message("proba = %e, do_ionization = %d", proba, do_ionization); */

  if (do_ionization) {
    /* Tag the particle */
    radiation_tag_part_as_ionized(pj, xpj);

    /* For on-the-fly stromgren sphere */
    /* radiation_consume_ionizing_photons(si, Delta_dot_N_ion); */
  }
}

  /* 4. Compute radiation pressure */
  /* const float p_rad = radiation_compute_radiation_pressure(sj); */
  /* const float delta_p_rad = weight * p_rad; */

  /* Add the radiation pressure radially outwards from the star */
  

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);
}

#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
