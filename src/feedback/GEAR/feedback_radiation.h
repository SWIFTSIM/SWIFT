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
#ifndef SWIFT_FEEDBACK_GEAR_FEEDBACK_RADIATION_H
#define SWIFT_FEEDBACK_GEAR_FEEDBACK_RADIATION_H

/**
 * @file src/feedback/GEAR/feedback_radiation.h
 * @brief Feedback hooks of the GEAR subgrid radiation that do not depend on
 * the feedback module: the snapshot ledgers of the packet transport, the MPI
 * exchange of the HII state and the engine hooks of the packet transport.
 */

#include "error.h"
#include "part.h"
#include "radiation_hii.h"
#include "radiation_packet.h"
#include "timestep_sync_part.h"

#include <string.h>

/**
 * @brief Write the packet ledgers of the gas, the stars and the sinks that
 * left the box and of the gas swallowed by black holes, summed over the ranks,
 * PE then LW, in snapshot units. Only with SWIFT_DEBUG_CHECKS_ISRF.
 *
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour_radiation(hid_t h_grp) {

#if defined(GEAR_SUBGRID_RADIATION_ISRF) && defined(SWIFT_DEBUG_CHECKS_ISRF)
  struct feedback_isrf_packet_removed_data l;
  double f;
  radiation_packet_left_box_total(&l, &f);
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    l.removed_energy[b] *= f;
    l.removed_photons[b] *= f;
    l.inherited_energy[b] *= f;
    l.inherited_photons[b] *= f;
  }
  io_write_attribute(h_grp, "PacketLeftBoxEscapedEnergies", DOUBLE,
                     l.removed_energy, ISRF_OPERATOR_COUNT);
  io_write_attribute(h_grp, "PacketLeftBoxEscapedPhotonNumbers", DOUBLE,
                     l.removed_photons, ISRF_OPERATOR_COUNT);
  io_write_attribute(h_grp, "PacketLeftBoxInheritedEnergies", DOUBLE,
                     l.inherited_energy, ISRF_OPERATOR_COUNT);
  io_write_attribute(h_grp, "PacketLeftBoxInheritedPhotonNumbers", DOUBLE,
                     l.inherited_photons, ISRF_OPERATOR_COUNT);

  struct feedback_isrf_packet_box_data x;
  radiation_packet_box_total(&x);
  double *v = (double *)&x;
  for (size_t k = 0; k < sizeof(x) / sizeof(double); k++) v[k] *= f;
  const struct feedback_isrf_packet_star_ledger_data *st = &x.stars_left_box;
  io_write_attribute(h_grp, "PacketStarsLeftBoxEmittedEnergies", DOUBLE,
                     st->emitted_energy, ISRF_OPERATOR_COUNT);
  io_write_attribute(h_grp, "PacketStarsLeftBoxEmittedPhotonNumbers", DOUBLE,
                     st->emitted_photons, ISRF_OPERATOR_COUNT);
  io_write_attribute(h_grp, "PacketStarsLeftBoxEscapedEnergies", DOUBLE,
                     st->escaped_energy, ISRF_OPERATOR_COUNT);
  io_write_attribute(h_grp, "PacketStarsLeftBoxEscapedPhotonNumbers", DOUBLE,
                     st->escaped_photons, ISRF_OPERATOR_COUNT);
  const struct {
    const char *name;
    const struct feedback_isrf_packet_removed_data *r;
  } groups[3] = {{"StarsLeftBox", &st->removed},
                 {"SinksLeftBox", &x.sinks_left_box},
                 {"SwallowedByBlackHoles", &x.gas_swallowed_by_black_holes}};
  for (int g = 0; g < 3; g++) {
    char name[64];
    snprintf(name, sizeof(name), "Packet%sRemovedEnergies", groups[g].name);
    io_write_attribute(h_grp, name, DOUBLE, groups[g].r->removed_energy,
                       ISRF_OPERATOR_COUNT);
    snprintf(name, sizeof(name), "Packet%sRemovedPhotonNumbers",
             groups[g].name);
    io_write_attribute(h_grp, name, DOUBLE, groups[g].r->removed_photons,
                       ISRF_OPERATOR_COUNT);
    snprintf(name, sizeof(name), "Packet%sInheritedEnergies", groups[g].name);
    io_write_attribute(h_grp, name, DOUBLE, groups[g].r->inherited_energy,
                       ISRF_OPERATOR_COUNT);
    snprintf(name, sizeof(name), "Packet%sInheritedPhotonNumbers",
             groups[g].name);
    io_write_attribute(h_grp, name, DOUBLE, groups[g].r->inherited_photons,
                       ISRF_OPERATOR_COUNT);
  }
#endif
}

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
 * A fresh claim (an untagged @p p, or a won collision) flags @p p for
 * time-step synchronization, as a claim on a local particle does; a renewal
 * by the holding star does not.
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
      timestep_sync_part(p);
    } else {
      *forfeited_cost = data->cost;
    }
    return;
  }

  const int fresh = !p->feedback_data.is_ionized && data->is_ionized;
  feedback_apply_hii_tag_claim(p, data);
  if (fresh) timestep_sync_part(p);
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

/**
 * @brief Start-up of the ISRF packet transport: sets the speed limit of the
 * packet light, with or without stars.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_startup(const struct engine *e) {
  radiation_packet_startup(e);
}

/**
 * @brief Follows the cell width after a regrid: the packet reach margin
 * depends on it.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_regrid(const struct engine *e) {
  radiation_packet_refresh_frame_width(e);
}

/**
 * @brief Prints the far light of the packet transport at the statistics
 * cadence. Collective.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_statistics(const struct engine *e) {
  radiation_packet_far_diagnostics(e);
}

/**
 * @brief Reports and resets the arrival bin counters of the step. Collective.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_step_diagnostics(const struct engine *e) {
  radiation_packet_step_diagnostics(e);
}

/**
 * @brief Checks the share of the packet light in the far arrival cells once
 * per step. Collective.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_step_report(const struct engine *e) {
  radiation_packet_far_share_report(e);
}

/**
 * @brief Prints the far share of the whole run. Collective.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_end_of_run(const struct engine *e) {
  radiation_packet_far_share_final(e);
}

/**
 * @brief Sums the ledger of the gas that left the box over the ranks before a
 * snapshot. Collective.
 *
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_radiation_transport_snapshot_reduce(const struct engine *e) {
  radiation_packet_reduce_left_box(e);
}

#endif /* SWIFT_FEEDBACK_GEAR_FEEDBACK_RADIATION_H */
