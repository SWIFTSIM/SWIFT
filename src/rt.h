/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_H
#define SWIFT_RT_H

/**
 * @file src/rt.h
 * @brief Branches between the different radiative transfer schemes. Also
 * contains some globally valid functions related to time bin data.
 */

/* Config parameters. */
#include <config.h>

/* Import the right RT definition */
#if defined(RT_NONE)
#include "./rt/none/rt.h"
#include "./rt/none/rt_iact.h"
#elif defined(RT_DEBUG)
#include "./rt/debug/rt.h"
#include "./rt/debug/rt_iact.h"
#elif defined(RT_GEAR)
#include "./rt/GEAR/rt.h"
#include "./rt/GEAR/rt_iact.h"
#elif defined(RT_SPHM1RT)
#include "./rt/SPHM1RT/rt.h"
#include "./rt/SPHM1RT/rt_iact.h"
#else
#error "Invalid choice of radiation scheme"
#endif

/**
 * @brief Prepare the rt *time step* quantities for a *hydro force* calculation.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void
rt_part_timestep_prepare_force(struct part *restrict p) {

  p->rt_time_data.min_ngb_time_bin = num_time_bins + 1;
}

/**
 * @brief Initialise RT time step data. This struct should be hidden from users,
 * so we do it for all schemes here.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void
rt_first_init_part_timestep_data(struct part *restrict p) {

  p->rt_time_data.time_bin = 0;
  rt_part_timestep_prepare_force(p);
}

/**
 * @brief Gathers neighbor data for RT during *hydro* force loops.
 * This is something all RT schemes should have in common, and something
 * users most likely never should be touching, so it's 'hidden' in the
 * top-level header file all schemes have in common.
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
__attribute__((always_inline)) INLINE static void runner_iact_rt_part_timebin(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

#ifndef RT_NONE
  /* Update the minimal time-bin */
  if (pj->rt_time_data.time_bin > 0) {
    pi->rt_time_data.min_ngb_time_bin =
        min(pi->rt_time_data.min_ngb_time_bin, pj->rt_time_data.time_bin);
  }

  if (pi->rt_time_data.time_bin > 0) {
    pj->rt_time_data.min_ngb_time_bin =
        min(pj->rt_time_data.min_ngb_time_bin, pi->rt_time_data.time_bin);
  }
#endif
}

/**
 * @brief Gathers neighbor data for RT during *hydro* force loops
 * (non-symmetric). This is something all RT schemes should have in common, and
 * something users most likely never should be touching, so it's 'hidden' in the
 * top-level header file all schemes have in common.
 *
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
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_part_timebin(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   const struct part *restrict pj,
                                   const float a, const float H) {

#ifndef RT_NONE
  /* Update the minimal time-bin */
  if (pj->rt_time_data.time_bin > 0) {
    pi->rt_time_data.min_ngb_time_bin =
        min(pi->rt_time_data.min_ngb_time_bin, pj->rt_time_data.time_bin);
  }
#endif
}

/**
 * @brief Prepare the rt *time step* quantities for a stars feedback calculation
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void
rt_spart_timestep_prepare_feedback(struct spart *restrict sp) {

  sp->min_rt_ngb_time_bin = num_time_bins + 1;
}

/**
 * @brief Initialise RT time step data.
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void
rt_first_init_spart_timestep_data(struct spart *restrict sp) {

  rt_spart_timestep_prepare_feedback(sp);
}

/**
 * @brief Gathers neighbor data for RT during *hydro* force loops
 * (non-symmetric). This is something all RT schemes should have in common, and
 * something users most likely never should be touching, so it's 'hidden' in the
 * top-level header file all schemes have in common.
 *
 *
 * @param r2 Comoving square distance between the two sparts.
 * @param dx Comoving vector separating sparts and part (sp - p).
 * @param hi Comoving smoothing-length of spart
 * @param hj Comoving smoothing-length of part
 * @param sp The spart.
 * @param p The part (not updated)
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_spart_timebin(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *restrict sp,
                                    const struct part *restrict p,
                                    const float a, const float H) {

#ifndef RT_NONE
  /* Update the minimal time-bin */
  if (p->rt_time_data.time_bin > 0) {
    sp->min_rt_ngb_time_bin =
        min(sp->min_rt_ngb_time_bin, p->rt_time_data.time_bin);
  }
#endif
}

#endif /* SWIFT_RT_H */
