/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_GEAR_CHEMISTRY_IACT_H
#define SWIFT_GEAR_CHEMISTRY_IACT_H

/**
 * @file GEAR/chemistry_iact.h
 * @brief Smooth metal interaction functions following the GEAR version of smooth metalicity.
 *
 * The interactions computed here are the ones presented in Wiersma, Schaye et al. 2009
 */

#include "cache.h"
#include "minmax.h"
#include "chemistry_struct.h"

/**
 * @brief Smooth metalicity
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
const struct chemistry_data *chem_data) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;
  float wj, wj_dx;

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the smooth metallicity */
  for(int i=0; i < chemistry_element_count; i++) {
    chi->smoothed_metal_mass_fraction[i] += mj * chj->metal_mass_fraction[i] * wi;
    chj->smoothed_metal_mass_fraction[i] += mi * chi->metal_mass_fraction[i] * wj;
  }
}

/**
 * @brief Non symmetric smooth metalicity
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
float r2, float *dx, float hi, float hj, struct part *pi, const struct part *pj,
const struct chemistry_data *chem_data) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the smooth metallicity */
  for(int i=0; i < chemistry_element_count; i++) {
    chi->smoothed_metal_mass_fraction[i] += mj * chj->metal_mass_fraction[i] * wi;
  }

}



#endif /* SWIFT_GEAR_CHEMISTRY_IACT_H */
