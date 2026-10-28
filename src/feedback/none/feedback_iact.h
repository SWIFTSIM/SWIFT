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
#ifndef SWIFT_NONE_FEEDBACK_IACT_H
#define SWIFT_NONE_FEEDBACK_IACT_H

/* Local includes */
#include "feedback.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * Nothing to do here for the no-feedback model.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xp Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param hydro_props The properties of the hydro scheme.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The properties of the cooling scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, const struct part *pj, const struct xpart *xp,
    const struct cosmology *cosmo, const struct feedback_props *fb_props,
    const struct hydro_props *hydro_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current) {}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * Nothing to do here for the no-feedback model.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xp Extra particle data
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param fb_props Properties of the feedback scheme.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The properties of the cooling scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 * @param time_base The time base used to compute integer times.
 * @param with_cosmology Are we running with cosmology on?
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct spart *si, struct part *pj, struct xpart *xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const *phys_const,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const integertime_t ti_current, const double time_base,
    const int with_cosmology) {}

/**
 * @brief `grad(u)` interaction (symmetric): no-op, this
 * feedback model does not track the ISRF band fields.
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
__attribute__((always_inline)) INLINE static void runner_iact_isrf_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

/**
 * @brief `grad(u)` interaction (non-symmetric): no-op, this
 * feedback model does not track the ISRF band fields.
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
runner_iact_nonsym_isrf_gradient(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, const float a,
                                 const float H) {}

/**
 * @brief Force-loop ISRF packet hop between two particles (symmetric).
 * Nothing to do here.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param book_i_to_j Book the share of i into j.
 * @param book_j_to_i Book the share of j into i.
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_packet_hop(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const int book_i_to_j, const int book_j_to_i) {}

/**
 * @brief Force-loop ISRF packet hop from particle i into j (non-symmetric).
 * Nothing to do here.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Sender.
 * @param pj Receiver.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param book Book the share.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_packet_hop(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   struct part *restrict pj, const float a,
                                   const float H, const int book) {}

#endif /* SWIFT_NONE_FEEDBACK_IACT_H */
