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
#ifndef SWIFT_EAGLE_CHEMISTRY_IACT_H
#define SWIFT_EAGLE_CHEMISTRY_IACT_H

/**
 * @file EAGLE/chemistry_iact.h
 * @brief Smooth metal interaction functions following the EAGLE model.
 */

/**
 * @brief do chemistry computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi;
  float wj;

  /* Get the masses. */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_eval(ui, &wi);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_eval(uj, &wj);

  /* Compute contribution to the smooth metallicity */
  for (int i = 0; i < chemistry_element_count; i++) {
    chi->smoothed_metal_mass_fraction[i] +=
        mj * chj->metal_mass_fraction[i] * wi;
    chj->smoothed_metal_mass_fraction[i] +=
        mi * chi->metal_mass_fraction[i] * wj;
  }

  // Smooth metal mass fraction of all metals
  chi->smoothed_metal_mass_fraction_total +=
      mj * chj->metal_mass_fraction_total * wi;
  chj->smoothed_metal_mass_fraction_total +=
      mi * chi->metal_mass_fraction_total * wj;

  // Smooth iron mass fraction from SNIa
  chi->smoothed_iron_mass_fraction_from_SNIa +=
      mj * chj->iron_mass_fraction_from_SNIa * wi;
  chj->smoothed_iron_mass_fraction_from_SNIa +=
      mi * chi->iron_mass_fraction_from_SNIa * wj;
}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_eval(ui, &wi);

  /* Compute contribution to the smooth metallicity */
  for (int i = 0; i < chemistry_element_count; i++) {
    chi->smoothed_metal_mass_fraction[i] +=
        mj * chj->metal_mass_fraction[i] * wi;
  }

  // Smooth metal mass fraction of all metals
  chi->smoothed_metal_mass_fraction_total +=
      mj * chj->metal_mass_fraction_total * wi;

  // Smooth iron mass fraction from SNIa
  chi->smoothed_iron_mass_fraction_from_SNIa +=
      mj * chj->iron_mass_fraction_from_SNIa * wi;
}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology) {}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology) {}

#endif /* SWIFT_EAGLE_CHEMISTRY_IACT_H */
