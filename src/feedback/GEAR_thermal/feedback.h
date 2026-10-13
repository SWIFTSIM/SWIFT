/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_GEAR_H
#define SWIFT_FEEDBACK_GEAR_H

#include "../GEAR/feedback_common.h"
#include "../GEAR/stellar_evolution.h"
#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <string.h>
#include <strings.h>

void feedback_update_part(struct part *p, struct xpart *xp,
                          const struct engine *e);
void feedback_end_density(struct part *p, struct xpart *xp,
                          const struct engine *e);
void feedback_part_has_no_neighbours(struct part *p, struct xpart *xp,
                                     const struct engine *e);
void feedback_end_gradient(struct part *p, const struct engine *e);
void feedback_end_force(struct part *p, const struct engine *e);
float feedback_compute_part_timestep(const struct part *restrict p,
                                     const struct engine *e);
void feedback_reset_part(struct part *p, struct xpart *xp,
                         const struct engine *e);
void feedback_init_part(struct part *p, const struct engine *e);
void feedback_first_init_part(struct part *restrict p);
int feedback_is_active(const struct spart *sp, const struct engine *e);
int feedback_is_star_dead(const struct spart *sp);
int feedback_should_inject_SN_feedback(const struct spart *sp);
int feedback_should_inject_wind_feedback(const struct spart *sp);
int feedback_should_inject_SN_mass_only(const struct spart *sp);
int feedback_should_inject_feedback(const struct spart *sp);
float feedback_compute_momentum_correction_factor_for_multiple_sn_events(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo);
void feedback_accumulate_kinetic_energy_for_multiple_sn_events(
    struct xpart *xp, const float mj, const float new_mass,
    const float v_pec[3], const float v_hubble[3], const double dp[3],
    const double dp_ejecta[3]);
void feedback_accumulate_wind_for_multiple_sn_events(
    struct xpart *xp, const struct spart *si, const float dx[3], const float r2,
    const double weight, const float mj, const double dm_SW,
    const double new_mass, const struct cosmology *cosmo);
void feedback_accumulate_SN_for_multiple_sn_events(
    struct xpart *xp, const struct spart *si, const float mj,
    const double dm_SN, const double new_mass, const struct cosmology *cosmo,
    const int is_event);
float feedback_compute_residual_internal_energy_for_multiple_sn_events(
    const struct part *p, const struct xpart *xp, const struct cosmology *cosmo,
    const float old_mass, const float new_mass, const float f_corr);
void feedback_init_spart(struct spart *sp);
void feedback_reset_feedback(struct spart *sp,
                             const struct feedback_props *feedback_props);
void feedback_prepare_spart(struct spart *sp,
                            const struct feedback_props *feedback_props);
void feedback_prepare_feedback(struct spart *restrict sp,
                               const struct feedback_props *feedback_props,
                               const struct cosmology *cosmo,
                               const struct unit_system *us,
                               const struct phys_const *phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology);

/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The #feedback_props.
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour(struct feedback_props *feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "GEAR");
};

/**
 * @brief Pack a #part's HII claim for the gas owner.
 *
 * The pass stamp on @p p is cleared after the pack, so a later pack without a
 * new claim does not report it again. The entry is zeroed first so that its
 * padding is deterministic on the wire.
 *
 * @param p The #part to pack from.
 * @param data The destination entry.
 */
__attribute__((always_inline)) INLINE static void feedback_pack_hii_tag_report(
    struct part *restrict p, struct hii_tag_report *restrict data) {

  memset(data, 0, sizeof(struct hii_tag_report));
#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  data->is_ionized = p->feedback_data.is_ionized;
  data->star_id = p->feedback_data.star_id;
  data->end_time = p->feedback_data.end_time;
  data->r2 = p->feedback_data.r2;
  data->cost = p->feedback_data.cost;
  data->claimed_this_pass = p->feedback_data.claimed_this_pass;
  p->feedback_data.claimed_this_pass = 0;
#endif
}

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
/**
 * @brief Apply an incoming HII claim onto @p p.
 *
 * Only the claim fields move: the cooling caches belong to this rank's own
 * cooling, and claimed_this_pass to this rank's own next pack.
 *
 * @param p The #part to apply the claim to.
 * @param data The source entry.
 */
__attribute__((always_inline)) INLINE static void feedback_apply_hii_tag_claim(
    struct part *restrict p, const struct hii_tag_report *restrict data) {
  p->feedback_data.is_ionized = data->is_ionized;
  p->feedback_data.star_id = data->star_id;
  p->feedback_data.end_time = data->end_time;
  p->feedback_data.r2 = data->r2;
  p->feedback_data.cost = data->cost;
}
#endif

/**
 * @brief Unpack an HII claim and merge it with the claim @p p already holds.
 *
 * Only a stamped entry is merged. Against a claim of another star, the
 * smallest (r2, star_id) wins; the photons of the loser are lost and only
 * counted. A report from the star that already holds @p p is applied
 * directly and is not a collision.
 *
 * @param p The #part to unpack into.
 * @param data The source entry.
 * @param collision (return) Set to 1 if the entry met a claim of another
 * star, left untouched otherwise.
 * @param forfeited_cost (return) Photon cost of the losing claim, set only
 * when @p collision is set.
 */
__attribute__((always_inline)) INLINE static void
feedback_unpack_hii_tag_report(struct part *restrict p,
                               const struct hii_tag_report *restrict data,
                               char *collision, double *forfeited_cost) {

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  if (!data->claimed_this_pass) return;

  if (p->feedback_data.is_ionized &&
      data->star_id != p->feedback_data.star_id) {
    *collision = 1;
#ifdef SWIFT_DEBUG_CHECKS
    if (!data->is_ionized)
      error("Stamped HII report for part %lld carries no live claim.", p->id);
#endif
    const int incoming_wins = (data->r2 < p->feedback_data.r2) ||
                              (data->r2 == p->feedback_data.r2 &&
                               data->star_id < p->feedback_data.star_id);
    if (incoming_wins) {
      *forfeited_cost = p->feedback_data.cost;
      feedback_apply_hii_tag_claim(p, data);
    } else {
      *forfeited_cost = data->cost;
    }
    return;
  }

  feedback_apply_hii_tag_claim(p, data);
#endif
}

/**
 * @brief Pack a #part's post-cooling HII state for the ranks that hold a
 * foreign copy.
 *
 * The entry is zeroed first so that its padding is deterministic on the
 * wire.
 *
 * @param p The #part to pack from.
 * @param data The destination entry.
 */
__attribute__((always_inline)) INLINE static void
feedback_pack_hii_state_update(const struct part *restrict p,
                               struct hii_state_update *restrict data) {

  memset(data, 0, sizeof(struct hii_state_update));
#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  data->is_ionized = p->feedback_data.is_ionized;
  data->star_id = p->feedback_data.star_id;
  data->end_time = p->feedback_data.end_time;
  data->r2 = p->feedback_data.r2;
  data->cost = p->feedback_data.cost;
  data->T_eligibility = p->feedback_data.T_eligibility;
  data->mu_eligibility = p->feedback_data.mu_eligibility;
  data->neutral_H_frac = p->feedback_data.neutral_H_frac;
#endif
}

/**
 * @brief Unpack a #part's post-cooling HII state into its foreign copy.
 *
 * Only the HII fields are written. claimed_this_pass is this rank's own
 * bookkeeping and is kept: taking the sender's value would echo a phantom
 * claim back to the owner on the next report.
 *
 * @param p The #part to unpack into.
 * @param data The source entry.
 */
__attribute__((always_inline)) INLINE static void
feedback_unpack_hii_state_update(struct part *restrict p,
                                 const struct hii_state_update *restrict data) {

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)
  p->feedback_data.is_ionized = data->is_ionized;
  p->feedback_data.star_id = data->star_id;
  p->feedback_data.end_time = data->end_time;
  p->feedback_data.r2 = data->r2;
  p->feedback_data.cost = data->cost;
  p->feedback_data.T_eligibility = data->T_eligibility;
  p->feedback_data.mu_eligibility = data->mu_eligibility;
  p->feedback_data.neutral_H_frac = data->neutral_H_frac;
#endif
}

#endif /* SWIFT_FEEDBACK_GEAR_H */
