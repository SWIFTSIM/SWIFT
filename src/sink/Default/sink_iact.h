/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_SINKS_IACT_H
#define SWIFT_DEFAULT_SINKS_IACT_H

/**
 * @brief do sink computation after the runner_iact_density (symmetric
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
 * @param cut_off_radius Sink cut off radius.
 */
__attribute__((always_inline)) INLINE static void runner_iact_sink(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

/**
 * @brief do sink computation after the runner_iact_density (non symmetric
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
 * @param cut_off_radius Sink cut off radius.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_sink(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First particle (sink).
 * @param pj Second particle (gas, not updated).
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sink *si, const struct part *pj, const int with_cosmology,
    const struct cosmology *cosmo, const struct gravity_props *grav_props,
    const struct sink_props *sink_props, const integertime_t ti_current,
    const double time) {

  float wi, wi_dx;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

#ifdef SWIFT_SINK_DENSITY_CHECKS
  si->rho_check += hydro_get_mass(pj) * wi;
  si->n_check += wi;
  si->N_check_density++;
#endif
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First particle (sink).
 * @param pj Second particle (gas, not updated).
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_density(
    const float r2, const float dx[3], const float ri, const float hj,
    struct sink *si, const struct part *pj, const int with_cosmology,
    const struct cosmology *cosmo, const struct gravity_props *grav_props,
    const struct sink_props *sink_props, const integertime_t ti_current,
    const double time) {

  /* Get r. */
  const float r = sqrtf(r2);

  if (r < ri) {
    /* Contribution to the number of neighbours in cutoff radius */
    si->num_ngbs++;
  }
}

/**
 * @brief Compute sink-sink swallow interaction (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing length of particle i.
 * @param hj Comoving smoothing length of particle j.
 * @param si First sink particle.
 * @param sj Second sink particle.
 * @param with_cosmology if we run with cosmology.
 * @param cosmo The cosmological parameters and properties.
 * @param grav_props The gravity scheme parameters and properties.
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_sink_swallow(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sink *restrict si, struct sink *restrict sj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct sink_props *sink_properties, const integertime_t ti_current,
    const double time) {}

/**
 * @brief Compute sink-gas swallow interaction (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 * @param with_cosmology if we run with cosmology.
 * @param cosmo The cosmological parameters and properties.
 * @param grav_props The gravity scheme parameters and properties.
 * @param sink_props the sink properties to use.
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sink *restrict si, struct part *restrict pj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct sink_props *sink_properties, const integertime_t ti_current,
    const double time) {}

#endif
