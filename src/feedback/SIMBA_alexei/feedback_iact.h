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
runner_iact_nonsym_feedback_density(
    const float r2, const float *dx, const float hi, const float hj,
    struct spart *restrict si, const struct part *restrict pj,
    const struct xpart *restrict xp, const struct cosmology *restrict cosmo,
    const struct feedback_props *fb_props, const integertime_t ti_current) {}

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
runner_iact_nonsym_feedback_apply(
    const float r2, const float *dx, const float hi, const float hj,
    const struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xp, const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const integertime_t ti_current) {

  /* Get the probability of doing feedback */
  // Compute mass loading which will determine probability
  const float prob = sqrtf(pj->potential) /
                     50.;  // ALEXEI: just set to random constant for now

  if (prob > 0) {
    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(si->id + pj->id, ti_current,
                                            random_number_stellar_feedback);

    /* Are we lucky? */
    if (rand < prob) {
      /* kick particle */
      for (int i = 0; i < 3; i++)
        pj->v[i] += si->feedback_data.to_distribute.v_kick;

      /* Heat particle */
      const float u_init = hydro_get_physical_internal_energy(pj, xp, cosmo);
      const float u_new = u_init + si->feedback_data.to_distribute.delta_u;
      hydro_set_physical_internal_energy(pj, xp, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      /* Set delaytime before which the particle cannot interact */
      pj->delay_time = si->feedback_data.to_distribute.simba_delay_time;
    }
  }
}

#endif /* SWIFT_SIMBA_FEEDBACK_IACT_H */
