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
void feedback_end_density(struct part *p, struct xpart *xp);
void feedback_reset_part(struct part *p, struct xpart *xp);
int feedback_is_active(const struct spart *sp, const struct engine *e);
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
 * @brief Pack a #part's HII tag report-back entry (MPI plan S3.1).
 *
 * excess_photon_energy_HI and photoionization_rate_HI have no per-part
 * storage yet and are left at their placeholder values; r2, cost and
 * claimed_this_pass (S3.4) are read from this pass's claim, if any, and
 * the this-pass stamp on @p p is drained immediately after so a later
 * pack without a fresh claim does not re-report a stale one. The struct
 * is memset first so its tail padding is deterministic on the wire (MSAN
 * hygiene).
 *
 * @param p The #part to pack from.
 * @param data The destination entry.
 */
__attribute__((always_inline)) INLINE static void feedback_pack_hii_tag_report(
    struct part *restrict p, struct hii_tag_report *restrict data) {

  memset(data, 0, sizeof(struct hii_tag_report));
  data->tag = p->feedback_data;
#ifdef WITH_MPI
  data->r2 = p->feedback_data.r2;
  data->cost = p->feedback_data.cost;
  data->claimed_this_pass = p->feedback_data.claimed_this_pass;
  p->feedback_data.claimed_this_pass = 0;
#endif
}

/**
 * @brief Unpack a #part's HII tag report-back entry and merge it against
 * any existing claim (MPI plan S3.1, S3.4).
 *
 * Only a stamped entry (S3.1/F2) is merged; an unstamped entry's tag state
 * is stale ambient data riding along in the pack, not a claim to apply.
 * A stamped entry competing against an existing claim (@p p already
 * is_ionized) is resolved by the deterministic commutative rule: smallest
 * (r2, star_id) wins, both values already on hand on either side (S3.4,
 * no owner-side star lookup). The loser's photon spend is forfeit and
 * only reported through @p collision / @p forfeited_cost for the debug
 * counters; nothing here can recover it.
 *
 * @param p The #part to unpack into.
 * @param data The source entry.
 * @param collision (return) Set to 1 if this entry competed against an
 * existing claim, left untouched otherwise.
 * @param forfeited_cost (return) Photon cost of whichever claim lost, set
 * only when @p collision is set.
 */
__attribute__((always_inline)) INLINE static void
feedback_unpack_hii_tag_report(struct part *restrict p,
                               const struct hii_tag_report *restrict data,
                               char *collision, float *forfeited_cost) {

  if (!data->claimed_this_pass) return;

#ifdef WITH_MPI
  if (p->feedback_data.is_ionized) {
    *collision = 1;
    const int incoming_wins = (data->r2 < p->feedback_data.r2) ||
                              (data->r2 == p->feedback_data.r2 &&
                               data->tag.star_id < p->feedback_data.star_id);
    if (incoming_wins) {
      *forfeited_cost = p->feedback_data.cost;
      p->feedback_data = data->tag;
    } else {
      *forfeited_cost = data->cost;
    }
    return;
  }
#endif

  p->feedback_data = data->tag;
}

/**
 * @brief Pack a #part's post-cooling HII state update entry (MPI plan
 * S3.1b).
 *
 * The struct is memset first so its tail padding is deterministic on the
 * wire (MSAN hygiene).
 *
 * @param p The #part to pack from.
 * @param data The destination entry.
 */
__attribute__((always_inline)) INLINE static void
feedback_pack_hii_state_update(const struct part *restrict p,
                               struct hii_state_update *restrict data) {

  memset(data, 0, sizeof(struct hii_state_update));
#ifdef WITH_MPI
  data->T_eligibility = p->feedback_data.T_eligibility;
#endif
  data->neutral_H_frac = p->feedback_data.neutral_H_frac;
  data->tag = p->feedback_data;
}

/**
 * @brief Unpack a #part's post-cooling HII state update entry (MPI plan
 * S3.1b).
 *
 * tag is written first as the wholesale mirror, not the primary read path
 * (see the struct's doxygen); the two standalone fields are the
 * authoritative same-step values and must win if they ever diverge.
 *
 * @param p The #part to unpack into.
 * @param data The source entry.
 */
__attribute__((always_inline)) INLINE static void
feedback_unpack_hii_state_update(struct part *restrict p,
                                 const struct hii_state_update *restrict data) {

  p->feedback_data = data->tag;
#ifdef WITH_MPI
  p->feedback_data.T_eligibility = data->T_eligibility;
#endif
  p->feedback_data.neutral_H_frac = data->neutral_H_frac;
}

#endif /* SWIFT_FEEDBACK_GEAR_H */
