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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_SETTERS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_SETTERS_H

#include "chemistry_getters.h"
#include "chemistry_struct.h"
#include "hydro.h"
#include "kernel_hydro.h"

/**
 * @brief Set the gradients for the given particle to zero.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_reset_gradients(struct part *restrict p) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chd->gradients.Z[i][0] = 0.0f;
    chd->gradients.Z[i][1] = 0.0f;
    chd->gradients.Z[i][2] = 0.0f;
  }

  chd->gradients.v[0][0] = 0.0f;
  chd->gradients.v[0][1] = 0.0f;
  chd->gradients.v[0][2] = 0.0f;
  chd->gradients.v[1][0] = 0.0f;
  chd->gradients.v[1][1] = 0.0f;
  chd->gradients.v[1][2] = 0.0f;
  chd->gradients.v[2][0] = 0.0f;
  chd->gradients.v[2][1] = 0.0f;
  chd->gradients.v[2][2] = 0.0f;

  chd->filtered.grad_v_tilde[0][0] = 0.0f;
  chd->filtered.grad_v_tilde[0][1] = 0.0f;
  chd->filtered.grad_v_tilde[0][2] = 0.0f;
  chd->filtered.grad_v_tilde[1][0] = 0.0f;
  chd->filtered.grad_v_tilde[1][1] = 0.0f;
  chd->filtered.grad_v_tilde[1][2] = 0.0f;
  chd->filtered.grad_v_tilde[2][0] = 0.0f;
  chd->filtered.grad_v_tilde[2][1] = 0.0f;
  chd->filtered.grad_v_tilde[2][2] = 0.0f;
}

/**
 * @brief Set the gradients for the given particle to the given values.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param gradF Metal mass fraction gradient (of size 3) to set.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_set_metal_mass_fraction_gradients(struct part *restrict p,
                                                 int metal,
                                                 const double gradF[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients.Z[metal][0] = gradF[0];
  chd->gradients.Z[metal][1] = gradF[1];
  chd->gradients.Z[metal][2] = gradF[2];
}

/**
 * @brief Update the diffusion gradients for the given particle with the
 * given contributions.
 *
 * @param p Particle.
 * @param metal metal specie index to update (0 <= metal <
 * GEAR_CHEMISTRY_ELEMENT_COUNT).
 * @param dF Metal mass fraction gradient (of size 3).
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_metal_mass_fraction_gradients(struct part *restrict p,
                                                    int metal, double dF[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients.Z[metal][0] += dF[0];
  chd->gradients.Z[metal][1] += dF[1];
  chd->gradients.Z[metal][2] += dF[2];
}

/**
 * @brief Update the velocity gradients for the given particle with the
 * given contributions.
 *
 * @param p Particle
 * @param dvx x Velocity gradient contribution.
 * @param dvy y Velocity gradient contribution.
 * @param dvz z Velocity gradient contribution.
 * @param dvx_tilde x Velocity tilde gradient contribution.
 * @param dvy_tilde y Velocity tilde gradient contribution.
 * @param dvz_tilde z Velocity tilde gradient contribution.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_hydro_gradients(struct part *restrict p, float dvx[3],
                                      float dvy[3], float dvz[3],
                                      float dvx_tilde[3], float dvy_tilde[3],
                                      float dvz_tilde[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients.v[0][0] += dvx[0];
  chd->gradients.v[0][1] += dvx[1];
  chd->gradients.v[0][2] += dvx[2];
  chd->gradients.v[1][0] += dvy[0];
  chd->gradients.v[1][1] += dvy[1];
  chd->gradients.v[1][2] += dvy[2];
  chd->gradients.v[2][0] += dvz[0];
  chd->gradients.v[2][1] += dvz[1];
  chd->gradients.v[2][2] += dvz[2];

  chd->filtered.grad_v_tilde[0][0] += dvx_tilde[0];
  chd->filtered.grad_v_tilde[0][1] += dvx_tilde[1];
  chd->filtered.grad_v_tilde[0][2] += dvx_tilde[2];
  chd->filtered.grad_v_tilde[1][0] += dvy_tilde[0];
  chd->filtered.grad_v_tilde[1][1] += dvy_tilde[1];
  chd->filtered.grad_v_tilde[1][2] += dvy_tilde[2];
  chd->filtered.grad_v_tilde[2][0] += dvz_tilde[0];
  chd->filtered.grad_v_tilde[2][1] += dvz_tilde[1];
  chd->filtered.grad_v_tilde[2][2] += dvz_tilde[2];
}

/**
 * @brief Normalise the gradients for the given particle with the given
 * normalisation factor.
 *
 * @param p Particle.
 * @param norm Normalisation factor.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_normalise_gradients(struct part *restrict p, const float norm) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chd->gradients.Z[i][0] *= norm;
    chd->gradients.Z[i][1] *= norm;
    chd->gradients.Z[i][2] *= norm;
  }

  chd->gradients.v[0][0] *= norm;
  chd->gradients.v[0][1] *= norm;
  chd->gradients.v[0][2] *= norm;
  chd->gradients.v[1][0] *= norm;
  chd->gradients.v[1][1] *= norm;
  chd->gradients.v[1][2] *= norm;
  chd->gradients.v[2][0] *= norm;
  chd->gradients.v[2][1] *= norm;
  chd->gradients.v[2][2] *= norm;

  chd->filtered.grad_v_tilde[0][0] *= norm;
  chd->filtered.grad_v_tilde[0][1] *= norm;
  chd->filtered.grad_v_tilde[0][2] *= norm;
  chd->filtered.grad_v_tilde[1][0] *= norm;
  chd->filtered.grad_v_tilde[1][1] *= norm;
  chd->filtered.grad_v_tilde[1][2] *= norm;
  chd->filtered.grad_v_tilde[2][0] *= norm;
  chd->filtered.grad_v_tilde[2][1] *= norm;
  chd->filtered.grad_v_tilde[2][2] *= norm;
}

/**
 * @brief Compute and set the ejected metal yields from supernovae events (SNII
 * and SNIa) for a star particle.
 *
 * This function calculates the total mass of metals ejected during supernova
 * feedback (Type II and Type Ia) for a given star particle. It combines the
 * yields from SNII and SNIa, accounts for unprocessed gas, and converts the
 * results into internal units.
 *
 * @param sp Pointer to the star particle structure (`struct spart`) where the
 * results will be stored.
 * @param m_snii Stellar mass involved per supernova II event.
 * @param m_non_processed Mass of unprocessed gas that retains the star's
 * initial metallicity.
 * @param number_snii Number of Type II supernovae events.
 * @param number_snia Number of Type Ia supernovae events.
 * @param snii_yields Array of metal yields per element for Type II supernovae.
 *                        The array size is `GEAR_CHEMISTRY_ELEMENT_COUNT`.
 * @param snia_yields Array of metal yields per element for Type Ia supernovae.
 *                        The array size is `GEAR_CHEMISTRY_ELEMENT_COUNT`.
 * @param phys_const Pointer to a structure containing physical constants.
 *
 * @note The resulting metal mass ejected per element is stored in:
 *       `sp->feedback_data.metal_mass_ejected[i]` for each element `i`.
 */
__attribute__((always_inline)) INLINE static void
chemistry_set_star_supernovae_ejected_yields(
    struct spart *restrict sp, const float mass_snii_event,
    const float m_non_processed, const int number_snii, const int number_snia,
    const float snii_yields[GEAR_CHEMISTRY_ELEMENT_COUNT],
    const float snia_yields[GEAR_CHEMISTRY_ELEMENT_COUNT],
    const struct phys_const *phys_const) {

  /* In MF diffusion, the last element correspond to the other untracked
     metals, not the sum of all metals. This ensure proper diffusion of the
     elements and consistency between the tracked elements and untracked ones */

  float snii_yields_new[GEAR_CHEMISTRY_ELEMENT_COUNT] = {0.f};
  float snia_yields_new[GEAR_CHEMISTRY_ELEMENT_COUNT] = {0.f};

  /* Get the sum of all explicitely tracked elements */
  float m_Z_tot_snii_tracked = 0.0;
  float m_Z_tot_snia_tracked = 0.0;
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT - 1; i++) {
    m_Z_tot_snii_tracked += snii_yields[i];
    m_Z_tot_snia_tracked += snia_yields[i];

    snii_yields_new[i] = snii_yields[i];
    snia_yields_new[i] = snia_yields[i];
  }

  const int last_elem = GEAR_CHEMISTRY_ELEMENT_COUNT - 1;
  snii_yields_new[last_elem] = snii_yields[last_elem] - m_Z_tot_snii_tracked;
  snia_yields_new[last_elem] =
      snia_yields_new[last_elem] - m_Z_tot_snia_tracked;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {

    /* Compute the mass fraction of metals */
    sp->feedback_data.metal_mass_ejected[i] =
        /* Supernovae II yields */
        snii_yields_new[i] +
        /* Gas contained in stars initial metallicity */
        chemistry_get_star_metal_mass_fraction_for_feedback(sp)[i] *
            m_non_processed;

    /* Convert it to total mass */
    sp->feedback_data.metal_mass_ejected[i] *= mass_snii_event * number_snii;

    /* Supernovae Ia yields */
    sp->feedback_data.metal_mass_ejected[i] += snia_yields_new[i] * number_snia;

    /* Convert everything in code units */
    sp->feedback_data.metal_mass_ejected[i] *= phys_const->const_solar_mass;
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_SETTERS_H */
