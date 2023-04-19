/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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
#ifndef SWIFT_SPHM1RT_RT_SETTERS_H
#define SWIFT_SPHM1RT_RT_SETTERS_H

/**
 * @file src/rt/SPHM1RT/rt_setters.h
 * @brief Independent setter functions for the SPHM1RT scheme
 */

/**
 * @brief Sets the comoving radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p The particle of interest.
 * @param urad The comoving radiation energy per mass
 *
 */
__attribute__((always_inline)) INLINE static void
rt_set_comoving_urad_multifrequency(struct part* p,
                                    const float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].urad = urad[g];
  }
}

/**
 * @brief Sets the physical radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param urad The physical radiation energy per mass
 */
__attribute__((always_inline)) INLINE static void
rt_set_physical_urad_multifrequency(struct part* p,
                                    const struct cosmology* cosmo,
                                    const float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].urad = urad[g];
  }
}

/**
 * @brief Sets the comoving radiation flux per density of a particle
 * (note that the comoving and physical flux per density are the same in our
 * convention)
 *
 * @param p The particle of interest.
 * @param frad The comoving radiation flux
 */
__attribute__((always_inline)) INLINE static void
rt_set_comoving_frad_multifrequency(struct part* p, float frad[RT_NGROUPS][3]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].frad[0] = frad[g][0];
    p->rt_data.conserved[g].frad[1] = frad[g][1];
    p->rt_data.conserved[g].frad[2] = frad[g][2];
  }
}

/**
 * @brief Sets the physical radiation flux of a particle
 * (note that the comoving and physical flux are the same in our convention)
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param frad The comoving radiation flux
 */
__attribute__((always_inline)) INLINE static void
rt_set_physical_radiation_flux_multifrequency(struct part* p,
                                              const struct cosmology* cosmo,
                                              float frad[RT_NGROUPS][3]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].frad[0] = frad[g][0] * cosmo->a_inv;
    p->rt_data.conserved[g].frad[1] = frad[g][1] * cosmo->a_inv;
    p->rt_data.conserved[g].frad[2] = frad[g][2] * cosmo->a_inv;
  }
}

/**
 * @brief Sets the physical opacity (chi) of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param chi The physical opacity
 */
__attribute__((always_inline)) INLINE static void
rt_set_physical_radiation_opacity(struct part* p, const struct cosmology* cosmo,
                                  const float chi[RT_NGROUPS]) {

  /* avoid getting negative opacity */
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.params.chi[g] = max(chi[g], 0.f) * cosmo->a_inv * cosmo->a_inv;
  }
}

#endif /* SWIFT_SPHM1RT_RT_SETTERS_H */
