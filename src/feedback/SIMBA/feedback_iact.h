/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SIMBA_FEEDBACK_IACT_H
#define SWIFT_SIMBA_FEEDBACK_IACT_H

/* Local includes */
#include "random.h"
#include "active.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xp Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float *dx,
                                    const float hi, const float hj,
                                    struct spart *restrict si,
                                    const struct part *restrict pj,
                                    const struct xpart *restrict xp,
                                    const struct cosmology *restrict cosmo,
                                    const integertime_t ti_current) {}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * In SIMBA model kick particle to model injection of energy due to SN, add
 * rest of energy as heating. ALEXEI: update comment
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xp Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  const struct spart *restrict si,
                                  struct part *restrict pj,
                                  struct xpart *restrict xp,
                                  const struct cosmology *restrict cosmo,
                                  const integertime_t ti_current) {
  
  if (ti_current == 0 || part_is_decoupled(pj)) return;

  /* Get the probability of doing feedback */
  //const float mass_frac = si->feedback_data.to_distribute.wind_mass/pj->mass;
  //const double prob_kick = 1. - exp(-mass_frac);
  const double prob_kick = 1./50.; // ALEXEI: placeholder
  const double prob_heat = 0.3; // ALEXEI: placeholder

  /* First we kick a particle */
  /* Draw a random number (Note mixing both IDs) */
  const double rand_kick = random_unit_interval(si->id + pj->id, ti_current, random_number_stellar_feedback);

  // ALEXEI: there might be some correlation here. maybe come up with a different seed (or think whether needed at all)
  const double rand_heat = random_unit_interval(si->id * pj->id, ti_current, random_number_stellar_feedback);

  /* Are we lucky? */
  if (rand_kick < prob_kick) {
    /* kick particle */
    float v_new[3];
    // ALEXEI: temporary simple definition. Change to be consistent with sfr_eff.c: 1565
    for (int i = 0; i < 3; i++) v_new[i] = si->feedback_data.to_distribute.v_kick; 
    hydro_set_velocity(pj,v_new);

    /* Heat particle */
    // probability GALSF_SUBGRID HOTWIND = 0.3 in SIMBA
    // Make sure star particle doesn't heat multiple times, i.e. it only launches eta*sm
    if (rand_heat > prob_heat) {
      const float u_init = hydro_get_physical_internal_energy(pj, xp, cosmo);
      const float u_new = u_init + si->feedback_data.to_distribute.delta_u;
      hydro_set_physical_internal_energy(pj, xp, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);
    }

    /* Set delaytime before which the particle cannot interact */
    pj->delay_time = si->feedback_data.to_distribute.simba_delay_time;
    pj->time_bin = time_bin_decoupled;
    message("spart id %llu decoupled particle %llu delay_time %.5e rand %.5e prob %.5e", si->id, pj->id, pj->delay_time, rand_kick, prob_kick);
  }

}

#endif /* SWIFT_SIMBA_FEEDBACK_IACT_H */
