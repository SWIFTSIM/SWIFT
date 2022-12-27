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
#ifndef SWIFT_SPHM1RT_RT_GETTERS_H
#define SWIFT_SPHM1RT_RT_GETTERS_H

/**
 * @file src/rt/SPHM1RT/rt_getters.h
 * @brief Independent getter functions for the SPHM1RT scheme
 */

/**
 * @brief Returns the comoving radiation speed of a particle
 * @param p Pointer to the particle data.
 * @param a Current scale factor.
 * @return comoving reduced speed of light
 */
__attribute__((always_inline)) INLINE static float rt_get_comoving_cred(
    const struct part* restrict p, float a) {
  return p->rt_data.params.cred_phys / a;
}

/**
 * @brief Returns the physical radiation speed of a particle
 *
 * @param p Pointer to the particle data.
 * @param a Current scale factor.
 * @return physical reduced speed of light
 *
 */
__attribute__((always_inline)) INLINE static float rt_get_physical_cred(
    const struct part* restrict p, float a) {
  return p->rt_data.params.cred_phys;
}

/**
 * @brief Returns the comoving radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param urad The comoving radiation energy per mass.
 *
 */
__attribute__((always_inline)) INLINE static void
rt_get_comoving_urad_multifrequency(const struct part* restrict p,
                                    float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    urad[g] = p->rt_data.conserved[g].urad;
  }
}

/**
 * @brief Returns the physical radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure.
 * @param urad The physical radiation energy.
 *
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_urad_multifrequency(const struct part* restrict p,
                                    const struct cosmology* cosmo,
                                    float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    urad[g] = p->rt_data.conserved[g].urad;
  }
}

/**
 * @brief Returns the comoving radiation flux per gas density of a particle
 *
 * @param p Pointer to the particle data.
 * @param fradtemp The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
rt_get_comoving_frad_multifrequency(const struct part* restrict p,
                                    float fradtemp[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[g][0] = p->rt_data.conserved[g].frad[0];
    fradtemp[g][1] = p->rt_data.conserved[g].frad[1];
    fradtemp[g][2] = p->rt_data.conserved[g].frad[2];
  }
}

/**
 * @brief Returns the physical radiation flux per gas density of a particle
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure
 * @param fradtemp The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_frad_multifrequency(const struct part* restrict p,
                                    const struct cosmology* cosmo,
                                    float fradtemp[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[g][0] = p->rt_data.conserved[g].frad[0] * cosmo->a;
    fradtemp[g][1] = p->rt_data.conserved[g].frad[1] * cosmo->a;
    fradtemp[g][2] = p->rt_data.conserved[g].frad[2] * cosmo->a;
  }
}


/**
 * @brief Returns the physical radiation energy per mass of a particle
 * per time from propagation
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure.
 * @param uradrate The physical radiation energy injection rate.
 *
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_urad_propagation_rate(const struct part* restrict p,
                                    const struct cosmology* cosmo,
                                    float uradrate[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    uradrate[g] = p->rt_data.dconserved_dt[g].urad;
  }
}


/**
 * @brief Returns the physical radiation energy per mass of a particle
 * from injection
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure.
 * @param uradinj The physical radiation energy injected.
 *
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_urad_injection(const struct part* restrict p,
                                    const struct cosmology* cosmo,
                                    float uradinj[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    uradinj[g] = p->rt_data.dconserved_inj[g].urad;
  }
}

/**
 * @brief Returns the physical radiation flux per gas density
 * of a particle from propagation of photon
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure
 * @param fradinj The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_frad_propagation_rate(const struct part* restrict p,
                                    const struct cosmology* cosmo,
                                    float fradpro[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradpro[g][0] = p->rt_data.dconserved_dt[g].frad[0] * cosmo->a;
    fradpro[g][1] = p->rt_data.dconserved_dt[g].frad[1] * cosmo->a;
    fradpro[g][2] = p->rt_data.dconserved_dt[g].frad[2] * cosmo->a;
  }
}


/**
 * @brief Returns the physical radiation flux per gas density 
 * of a particle from injection of photon 
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure
 * @param fradinj The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_frad_injection(const struct part* restrict p,
                                    const struct cosmology* cosmo,
                                    float fradinj[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradinj[g][0] = p->rt_data.dconserved_inj[g].frad[0] * cosmo->a;
    fradinj[g][1] = p->rt_data.dconserved_inj[g].frad[1] * cosmo->a;
    fradinj[g][2] = p->rt_data.dconserved_inj[g].frad[2] * cosmo->a;
  }
}

#endif /* SWIFT_SPHM1RT_RT_GETTERS_H */
