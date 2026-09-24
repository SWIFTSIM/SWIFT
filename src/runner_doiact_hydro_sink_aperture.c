/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "runner.h"
#include "sink.h"
#include "space_getsid.h"
#include "timers.h"

/* This object's header. */
#include "runner_doiact_hydro_sink_aperture.h"

/**
 * @brief Test one active gas particle against every (non-reserved) sink in
 * @p sinks, recording overlaps on the gas particle.
 *
 * @param e The #engine.
 * @param pi The active gas particle.
 * @param xpi Its #xpart.
 * @param sinks The sink array to test against.
 * @param scount The number of sinks in @p sinks.
 * @param with_cosmology Whether we are running with cosmology.
 * @param cosmo The #cosmology.
 * @param sink_props The #sink_props of this run.
 * @param r_acc_p Accretion radius of #pi if it forms a sink. Same for all
 * sinks, so the caller computes it once.
 */
static INLINE void runner_iact_hydro_sink_aperture_prep_sink_formation_sink(
    struct engine *e, struct part *restrict pi, struct xpart *restrict xpi,
    struct sink *restrict sinks, const int scount, const int with_cosmology,
    const struct cosmology *cosmo, const struct sink_props *sink_props,
    const float r_acc_p) {

  /* Most gas cannot form a sink. Skip the loop over the sinks for it. */
  if (!pi->sink_data.can_form_sink) return;

  /* Box size, or 0 if the box is not periodic */
  const struct space *s = e->s;
  const double dim[3] = {s->periodic ? s->dim[0] : 0.,
                         s->periodic ? s->dim[1] : 0.,
                         s->periodic ? s->dim[2] : 0.};

  for (int sjd = 0; sjd < scount; sjd++) {

    struct sink *restrict sj = &sinks[sjd];

    /* Ignore inhibited sinks and empty slots reserved for new sinks */
    if (sink_is_inhibited(sj, e) || sj->time_bin == time_bin_not_created)
      continue;

    sink_prepare_part_sink_formation_sink_criteria(
        e, pi, xpi, sj, with_cosmology, cosmo, sink_props, e->time, r_acc_p,
        dim);
  }
}

/**
 * @brief Compute cell self-interactions: every active gas particle in @p c
 * against every sink in @p c.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius (unused at the leaf: the physics
 * check in sink_prepare_part_sink_formation_sink_criteria() is itself the
 * only distance test needed; @p r_cut only drives DOSUB recursion).
 */
void runner_doself1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *c, const float r_cut) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct sink_props *sink_props = e->sink_properties;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->sinks.count == 0 || c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (!cell_are_part_drifted(c, e))
    error("Interacting undrifted cell (parts).");
  if (!cell_are_sink_drifted(c, e))
    error("Interacting undrifted cell (sinks).");
#endif

  const int count = c->hydro.count;
  const int scount = c->sinks.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct sink *restrict sinks = c->sinks.parts;

  for (int pid = 0; pid < count; pid++) {

    struct part *restrict pi = &parts[pid];

    if (part_is_inhibited(pi, e)) continue;
    if (!part_is_active(pi, e)) continue;

    struct xpart *restrict xpi = &xparts[pid];

    const float r_acc_p =
        (sink_props->use_fixed_r_cut ? sink_props->cut_off_radius
                                     : kernel_gamma * pi->h) *
        cosmo->a;

    runner_iact_hydro_sink_aperture_prep_sink_formation_sink(
        e, pi, xpi, sinks, scount, with_cosmology, cosmo, sink_props, r_acc_p);
  }

  TIMER_TOC(timer_doself_hydro_sink_aperture_prep_sink_formation_sink);
}

/**
 * @brief Non-symmetric half of the pair interaction: every active gas
 * particle in @p c_gas against every sink in @p c_sink.
 */
static void do_nonsym_pair1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, const struct cell *restrict c_gas,
    const struct cell *restrict c_sink) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct sink_props *sink_props = e->sink_properties;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  if (c_sink->sinks.count == 0 || c_gas->hydro.count == 0) return;
  if (c_gas->nodeID != e->nodeID) return;
  if (!cell_is_active_hydro(c_gas, e)) return;

  const int count = c_gas->hydro.count;
  const int scount = c_sink->sinks.count;
  struct part *restrict parts = c_gas->hydro.parts;
  struct xpart *restrict xparts = c_gas->hydro.xparts;
  struct sink *restrict sinks = c_sink->sinks.parts;

  for (int pid = 0; pid < count; pid++) {

    struct part *restrict pi = &parts[pid];

    if (part_is_inhibited(pi, e)) continue;
    if (!part_is_active(pi, e)) continue;

    struct xpart *restrict xpi = &xparts[pid];

    const float r_acc_p =
        (sink_props->use_fixed_r_cut ? sink_props->cut_off_radius
                                     : kernel_gamma * pi->h) *
        cosmo->a;

    runner_iact_hydro_sink_aperture_prep_sink_formation_sink(
        e, pi, xpi, sinks, scount, with_cosmology, cosmo, sink_props, r_acc_p);
  }
}

/**
 * @brief Compute the interactions between a cell pair: active gas in @p ci
 * against sinks in @p cj, and active gas in @p cj against sinks in @p ci.
 *
 * Naive (no sorted scan bounds): sink counts per cell are small, so the
 * cost is already linear in the active gas count.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void runner_dopair1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if ((ci->hydro.count == 0 || cj->sinks.count == 0) &&
      (cj->hydro.count == 0 || ci->sinks.count == 0))
    return;
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells (parts).");
  if (!cell_are_sink_drifted(ci, e) || !cell_are_sink_drifted(cj, e))
    error("Interacting undrifted cells (sinks).");
#endif

  do_nonsym_pair1_hydro_sink_aperture_prep_sink_formation_sink(r, ci, cj);
  do_nonsym_pair1_hydro_sink_aperture_prep_sink_formation_sink(r, cj, ci);

  TIMER_TOC(timer_dopair_hydro_sink_aperture_prep_sink_formation_sink);
}

/**
 * @brief Recursively compute self interactions for sub-cells.
 *
 * We stop splitting when 2 * r_cut >= 0.5 * dmin. The gas-gas loop uses
 * r_cut here. We need 2 * r_cut because the overlap test adds two radii.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius used by the gas-gas loop.
 * @param gettimer Whether to record a timer for this call.
 */
void runner_dosub_self1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *c, const float r_cut, const int gettimer) {

  TIMER_TIC;

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->sinks.count == 0 || !cell_is_active_hydro(c, e))
    return;

  if (!c->split || c->hydro.count < space_recurse_size_self_hydro ||
      2.f * r_cut >= 0.5f * c->dmin) {

    runner_doself1_hydro_sink_aperture_prep_sink_formation_sink(r, c, r_cut);

  } else {

    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        runner_dosub_self1_hydro_sink_aperture_prep_sink_formation_sink(
            r, c->progeny[k], r_cut, /*gettimer=*/0);

        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            runner_dosub_pair1_hydro_sink_aperture_prep_sink_formation_sink(
                r, c->progeny[k], c->progeny[j], r_cut, /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer)
    TIMER_TOC(timer_dosub_self_hydro_sink_aperture_prep_sink_formation_sink);
}

/**
 * @brief Recursively compute pair interactions for sub-cells.
 *
 * Recursion threshold is 2 * r_cut < 0.5 * dmin (see the self variant
 * above for why).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius used by the gas-gas loop.
 * @param gettimer Whether to record a timer for this call.
 */
void runner_dosub_pair1_hydro_sink_aperture_prep_sink_formation_sink(
    struct runner *r, struct cell *ci, struct cell *cj, const float r_cut,
    const int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Nothing to do if no active gas can see a sink in the other cell. Check
     both directions. */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;
  if ((ci->hydro.count == 0 || cj->sinks.count == 0) &&
      (cj->hydro.count == 0 || ci->sinks.count == 0))
    return;

  /* Get the pair direction and put the two cells in the standard order. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  if (!ci->split || ci->hydro.count < space_recurse_size_pair_hydro ||
      !cj->split || cj->hydro.count < space_recurse_size_pair_hydro ||
      2.f * r_cut >= 0.5f * ci->dmin) {

    runner_dopair1_hydro_sink_aperture_prep_sink_formation_sink(r, ci, cj);

  } else {

    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        runner_dosub_pair1_hydro_sink_aperture_prep_sink_formation_sink(
            r, ci->progeny[pid], cj->progeny[pjd], r_cut, /*gettimer=*/0);
      }
    }
  }

  if (gettimer)
    TIMER_TOC(timer_dosub_pair_hydro_sink_aperture_prep_sink_formation_sink);
}
