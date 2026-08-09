/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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

/* This object's header. */
#include "physical_constants.h"
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "sink.h"
#include "sink_iact.h"
#include "space_getsid.h"
#include "timers.h"

/**
 * @brief Process all the gas particles in a cell that have been flagged for
 * swallowing by a sink.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of sink for the particle with the corresponding
 * ID. If found, the sink swallows the gas particle and the gas particle is
 * removed. If the cell is local, we may be looking for a foreign sink, in which
 * case, we do not update the sink (that will be done on its node) but just
 * remove the gas particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;

  struct sink *sinks = s->sinks;
  const size_t nr_sink = s->nr_sinks;
#ifdef WITH_MPI
  struct sink *sinks_foreign = s->sinks_foreign;
  const size_t nr_sinks_foreign = s->nr_sinks_foreign;
#endif

  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;

  integertime_t ti_current = e->ti_current;
  integertime_t ti_beg_max = 0;

  /* Early abort?
   * (We only want cells for which we drifted the gas as these are
   * the only ones that could have gas particles that have been flagged
   * for swallowing) */
  if (c->hydro.count == 0 || c->hydro.ti_old_part != e->ti_current) {
    return;
  }

  /* Loop over the progeny ? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        runner_do_sinks_gas_swallow(r, cp, 0);

        /* Propagate the ti_beg_max from the leaves to the roots.
         * See bug fix below. */
        ti_beg_max = max(cp->hydro.ti_beg_max, ti_beg_max);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell
     * Note that the cell (and hence the parts) may be local or foreign. */
    const size_t nr_parts = c->hydro.count;
    for (size_t k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* Ignore inhibited particles (they have already been removed!) */
      if (part_is_inhibited(p, e)) continue;

      /* Get the ID of the sink that will swallow this part */
      const long long swallow_id = sink_get_part_swallow_id(&p->sink_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

#ifdef SWIFT_DEBUG_CHECKS
        if (p->ti_drift != e->ti_current)
          error("Trying to swallow an un-drifted particle.");
#endif

        /* ID of the sink swallowing this particle */
        const long long sink_id = swallow_id;

        /* Have we found this particle's sink already? */
        int found = 0;

        /* Let's look for the hungry sink in the local list */
        for (size_t i = 0; i < nr_sink; ++i) {

          /* Get a handle on the bpart. */
          struct sink *sp = &sinks[i];

          if (sp->id == sink_id) {

            /* Lock the space as we are going to work directly on the spart list
             */
            lock_lock(&s->lock);

            /* Swallow the gas particle (i.e. update the sink properties) */
            sink_swallow_part(sp, p, xp, e->cosmology);

            /* Release the space as we are done updating the spart */
            if (lock_unlock(&s->lock) != 0)
              error("Failed to unlock the space.");

#ifdef SWIFT_DEBUG_CHECKS
            message(
                "sink %lld (node %d) swallowing gas particle %lld (node %d, "
                "cellID=%lld, local case) at step %d",
                sp->id, e->nodeID, p->id, c->nodeID, c->cellID, e->step);
#endif

            /* If the gas particle is local, remove it */
            if (c->nodeID == e->nodeID) {

              lock_lock(&e->s->lock);

              /* Re-check that the particle has not been removed
               * by another thread before we do the deed. */
              if (!part_is_inhibited(p, e)) {

                /* Finally, remove the gas particle from the system
                 * Recall that the gpart associated with it is also removed
                 * at the same time. */
                cell_remove_part(e, c, p, xp);
              }

              if (lock_unlock(&e->s->lock) != 0)
                error("Failed to unlock the space!");
            }

            /* In any case, prevent the particle from being re-swallowed */
            sink_mark_part_as_swallowed(&p->sink_data);

            found = 1;
            break;
          }

        } /* Loop over local sinks */

#ifdef WITH_MPI
        /* We could also be in the case of a local gas particle being
         * swallowed by a foreign sink. In this case, we won't update the
         * sink but just remove the particle from the local list. */
        if (c->nodeID == e->nodeID && !found) {

          /* Let's look for the foreign hungry black hole */
          for (size_t i = 0; i < nr_sinks_foreign; ++i) {

            /* Get a handle on the sink. */
            struct sink *sink = &sinks_foreign[i];

            if (sink->id == sink_id) {
#ifdef SWIFT_DEBUG_CHECKS
              message(
                  "Sink %lld removing gas particle %lld (foreign sink case)",
                  sink->id, p->id);
#endif /* SWIFT_DEBUG_CHECKS */

              lock_lock(&e->s->lock);

              /* Re-check that the particle has not been removed
               * by another thread before we do the deed. */
              if (!part_is_inhibited(p, e)) {

                /* Finally, remove the gas particle from the system */
                cell_remove_part(e, c, p, xp);
              }

              if (lock_unlock(&e->s->lock) != 0)
                error("Failed to unlock the space!");

              found = 1;
              break;
            }
          } /* Loop over foreign sinks */
        } /* Is the cell local? */
#endif

        /* If we have a local particle, we must have found the sink in one
         * of our list of sinks. */
        if (c->nodeID == e->nodeID && !found) {
          error("Gas particle %lld could not find sink %lld to be swallowed",
                p->id, swallow_id);
        }
      } /* Part was flagged for swallowing */

      /* Bug fix : Change the hydro.ti_beg_max when a sink eats the last gas
       * particle possessing the ti_beg_max of the cell. We set hydro.ti_beg_max
       * to the max ti_beg of the remaining gas particle. Why this fix ?
       * Otherwise, we fail the check from cell_check_timesteps. This bug is
       * rare because it needs that the swallowed gas is the last part with the
       * ti_beg_max of the cell.
       * The same is not done for ti_end_min since it may inactivate cells that
       * need to perform sinks tasks.
       */

      if (part_is_inhibited(p, e)) continue;

      integertime_t ti_beg =
          get_integer_time_begin(ti_current + 1, p->time_bin);
      ti_beg_max = max(ti_beg, ti_beg_max);
    } /* Loop over the parts */
  } /* Cell is not split */

  /* Update ti_beg_max. See bug fix above. */
  if (ti_beg_max != c->hydro.ti_beg_max) {
    c->hydro.ti_beg_max = ti_beg_max;
  }
}

/**
 * @brief Processing of gas particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_self(struct runner *r, struct cell *c,
                                      int timer) {
#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != r->e->nodeID) error("Running self task on foreign node");
  if (!cell_is_active_sinks(c, r->e) && !cell_is_active_hydro(c, r->e))
    error("Running self task on inactive cell");
#endif

  runner_do_sinks_gas_swallow(r, c, timer);
}

/**
 * @brief Processing of gas particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_pair(struct runner *r, struct cell *ci,
                                      struct cell *cj, int timer) {
#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active sink */
  if (cell_is_active_sinks(cj, e)) runner_do_sinks_gas_swallow(r, ci, timer);
  if (cell_is_active_sinks(ci, e)) runner_do_sinks_gas_swallow(r, cj, timer);
}

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Count sinks in a cell's subtree flagged for a pending merger.
 *
 * A cell's sinks.parts/sinks.count already spans its whole subtree as a
 * contiguous sub-range of the space-wide sink array (see
 * cell_split.c:cell_split(), which sets each progeny's sinks.parts to a
 * sub-range of the parent's), so a flat loop over this cell alone is
 * enough -- no recursion into progeny is required.
 *
 * @param c The #cell (split or leaf).
 * @param e The #engine.
 * @return The number of non-inhibited sinks with a pending swallow_id.
 */
static int runner_sinks_count_pending_merges(struct cell *c,
                                             const struct engine *e) {
  int n_pending = 0;
  for (int k = 0; k < c->sinks.count; k++) {
    struct sink *sp = &c->sinks.parts[k];
    if (sink_is_inhibited(sp, e)) continue;
    if (sink_get_sink_swallow_id(&sp->merger_data) != -1) n_pending++;
  }
  return n_pending;
}
#endif

/**
 * @brief Process all the sink particles in a cell that have been flagged for
 * swallowing by a sink.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of sinks for the particle with the corresponding
 * ID. If found, the sink swallows the sink particle and the sink particle is
 * removed. If the cell is local, we may be looking for a foreign sink, in which
 * case, we do not update the sink (that will be done on its node) but just
 * remove the sink particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;

  struct sink *sinks = s->sinks;
  const size_t nr_sink = s->nr_sinks;
#ifdef WITH_MPI
  struct sink *sinks_foreign = s->sinks_foreign;
  const size_t nr_sinks_foreign = s->nr_sinks_foreign;
#endif

  struct sink *cell_sinks = c->sinks.parts;

  /* Early abort?
   * (We only want cells for which we drifted the sink as these are
   * the only ones that could have sink particles that have been flagged
   * for swallowing) */
  if (c->sinks.count == 0 || c->sinks.ti_old_part != e->ti_current) {
#ifdef SWIFT_DEBUG_CHECKS
    const int n_pending = runner_sinks_count_pending_merges(c, e);
    if (n_pending > 0)
      message(
          "DISPATCH PROBE A: runner_do_sinks_sink_swallow early-returns on "
          "cellID=%lld depth=%d nodeID=%d sinks.count=%d ti_old_part=%lld "
          "ti_current=%lld is_hydro_super=%d n_pending_marked=%d",
          c->cellID, c->depth, c->nodeID, c->sinks.count,
          (long long)c->sinks.ti_old_part, (long long)e->ti_current,
          c == c->hydro.super, n_pending);
#endif
    return;
  }

  /* Loop over the progeny ? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        runner_do_sinks_sink_swallow(r, cp, 0);
      }
    }
  } else {

    /* Loop over all the sinks particles in the cell
     * Note that the cell (and hence the sinks) may be local or foreign. */
    const size_t nr_cell_sinks = c->sinks.count;

    for (size_t k = 0; k < nr_cell_sinks; k++) {

      /* Get a handle on the part. */
      struct sink *const cell_sp = &cell_sinks[k];

      /* Ignore inhibited particles (they have already been removed!) */
      if (sink_is_inhibited(cell_sp, e)) continue;

      /* Get the ID of the sink that will swallow this sink */
      const long long swallow_id =
          sink_get_sink_swallow_id(&cell_sp->merger_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

#ifdef SWIFT_DEBUG_CHECKS
        if (cell_sp->ti_drift != e->ti_current)
          error("Trying to swallow an un-drifted particle.");
#endif

        /* ID of the sink swallowing this particle */
        const long long sink_id = swallow_id;

        /* Have we found this particle's sink already? */
        int found = 0;

        /* Let's look for the hungry sink in the local list */
        for (size_t i = 0; i < nr_sink; ++i) {

          /* Get a handle on the bpart. */
          struct sink *sp = &sinks[i];

          if (sp->id == sink_id) {

            /* Is the swallowing sink itself flagged for swallowing by
               another sink? */
            if (sink_get_sink_swallow_id(&sp->merger_data) != -1) {

              /* Pretend it was found and abort */
              sink_mark_sink_as_not_swallowed(&cell_sp->merger_data);
              found = 1;
              break;
            }

            /* Lock the space as we are going to work directly on the
             * space's bpart list */
            lock_lock(&s->lock);

            /* Swallow the sink particle (i.e. update the swallowing sink
             * properties with the properties of cell_sp) */
            sink_swallow_sink(sp, cell_sp, e->cosmology);

            /* Release the space as we are done updating the spart */
            if (lock_unlock(&s->lock) != 0)
              error("Failed to unlock the space.");

#ifdef SWIFT_DEBUG_CHECKS
            message(
                "sink %lld (node %d) swallowing sink particle %lld (node %d, "
                "cellID=%lld, local case) at step %d",
                sp->id, e->nodeID, cell_sp->id, c->nodeID, c->cellID, e->step);
#endif

            /* If the sink particle is local, remove it */
            if (c->nodeID == e->nodeID) {

              /* Finally, remove the sink particle from the system
               * Recall that the gpart associated with it is also removed
               * at the same time. */
              cell_remove_sink(e, c, cell_sp);
            }

            /* In any case, prevent the particle from being re-swallowed */
            sink_mark_sink_as_merged(&cell_sp->merger_data);

            found = 1;
            break;
          }

        } /* Loop over local sinks */

#ifdef WITH_MPI
        /* We could also be in the case of a local sink particle being
         * swallowed by a foreign sink. In this case, we won't update the
         * foreign sink but just remove the particle from the local list. */
        if (c->nodeID == e->nodeID && !found) {

          /* Let's look for the foreign hungry sink */
          for (size_t i = 0; i < nr_sinks_foreign; ++i) {

            /* Get a handle on the sink. */
            struct sink *sink = &sinks_foreign[i];

            if (sink->id == sink_id) {

              /* Is the swallowing sink itself flagged for swallowing by
                 another sink? */
              if (sink_get_sink_swallow_id(&sink->merger_data) != -1) {

                /* Pretend it was found and abort */
                sink_mark_sink_as_not_swallowed(&cell_sp->merger_data);
                found = 1;
                break;
              }

#ifdef SWIFT_DEBUG_CHECKS
              message(
                  "Sink %lld (foreign) swallowing sink particle %lld "
                  "(node %d, cellID=%lld, foreign sink case) at step %d",
                  sink->id, cell_sp->id, c->nodeID, c->cellID, e->step);
#endif /* SWIFT_DEBUG_CHECKS */

              /* Finally, remove the gas particle from the system */
              cell_remove_sink(e, c, cell_sp);

              found = 1;
              break;
            }
          } /* Loop over foreign sinks */
        } /* Is the cell local? */
#endif

        /* If we have a local particle, we must have found the sink in one
         * of our list of sinks. */
        if (c->nodeID == e->nodeID && !found) {
          error("sink particle %lld could not find sink %lld to be swallowed",
                cell_sp->id, swallow_id);
        }

      } /* Part was flagged for swallowing */
    } /* Loop over the parts */
  } /* Cell is not split */
}

/**
 * @brief Processing of sink particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow_self(struct runner *r, struct cell *c,
                                       int timer) {
#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  const int n_pending = runner_sinks_count_pending_merges(c, r->e);
  if (n_pending > 0)
    message(
        "DISPATCH PROBE B: runner_do_sinks_sink_swallow_self entered on "
        "ci_cellID=%lld ci_nodeID=%d depth=%d step=%d n_pending_marked=%d",
        c->cellID, c->nodeID, c->depth, r->e->step, n_pending);

  if (c->nodeID != r->e->nodeID) error("Running self task on foreign node");
  if (!cell_is_active_sinks(c, r->e) && !cell_is_active_hydro(c, r->e))
    error("Running self task on inactive cell");
#endif

  runner_do_sinks_sink_swallow(r, c, timer);
}

/**
 * @brief Processing of sink particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow_pair(struct runner *r, struct cell *ci,
                                       struct cell *cj, int timer) {
#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  const int ci_pending = runner_sinks_count_pending_merges(ci, e);
  const int cj_pending = runner_sinks_count_pending_merges(cj, e);
  if (ci_pending > 0 || cj_pending > 0)
    message(
        "DISPATCH PROBE B: runner_do_sinks_sink_swallow_pair entered on "
        "ci_cellID=%lld ci_nodeID=%d cj_cellID=%lld cj_nodeID=%d depth=%d "
        "step=%d ci_pending_marked=%d cj_pending_marked=%d",
        ci->cellID, ci->nodeID, cj->cellID, cj->nodeID, ci->depth, e->step,
        ci_pending, cj_pending);

  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");

  /* cell_unskip.c activates this task on cell_is_active_sinks(x) ||
   * cell_is_active_hydro(x), but we only process a side when the
   * neighbour is cell_is_active_sinks -- if a side is hydro-active but
   * not sinks-active, the task runs yet silently skips that side's
   * flagged sinks. Flag the exact scenario when it happens. */
  if ((ci->sinks.count > 0 || cj->sinks.count > 0) &&
      !cell_is_active_sinks(ci, e) && !cell_is_active_sinks(cj, e) &&
      (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e))) {
    message(
        "UNSKIP/RUNTIME MISMATCH: do_sink_swallow pair ci_cellID=%lld "
        "cj_cellID=%lld ran but neither side is sinks-active (only "
        "hydro-active) -- ci.active_hydro=%d cj.active_hydro=%d "
        "ci.sinks.count=%d cj.sinks.count=%d",
        ci->cellID, cj->cellID, cell_is_active_hydro(ci, e),
        cell_is_active_hydro(cj, e), ci->sinks.count, cj->sinks.count);
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Race probe: compare this foreign cell's cached activity read against its
   * read at the marking gate in DOPAIR1_BRANCH_SINKS. */
  const int cj_active_probe = cell_is_active_sinks(cj, e);
  const int ci_active_probe = cell_is_active_sinks(ci, e);
  if (cj->nodeID != e->nodeID && cj->sinks.count > 0)
    message(
        "RESOLVE_GATE foreign_cellID=%lld ti_end_min=%lld active=%d step=%d",
        cj->cellID, cj->sinks.ti_end_min, cj_active_probe, e->step);
  if (ci->nodeID != e->nodeID && ci->sinks.count > 0)
    message(
        "RESOLVE_GATE foreign_cellID=%lld ti_end_min=%lld active=%d step=%d",
        ci->cellID, ci->sinks.ti_end_min, ci_active_probe, e->step);
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active sink */
  if (cell_is_active_sinks(cj, e)) runner_do_sinks_sink_swallow(r, ci, timer);
  if (cell_is_active_sinks(ci, e)) runner_do_sinks_sink_swallow(r, cj, timer);
}

/**
 * @brief Compute the energies (kinetic, potential, etc ) of the gas particle
 * p and all quantities required for the formation of a sink.
 *
 * Note: This function iterates over gas particles and sink particles.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param p The #part.
 * @param xp The #xpart data of the particle p.
 */
void runner_do_prepare_part_sink_formation(struct runner *r, struct cell *c,
                                           struct part *restrict pi,
                                           struct xpart *restrict xpi) {
  struct engine *e = r->e;
  struct space *s = e->s;

  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const struct sink_props *sink_props = e->sink_properties;

  /* Loop over gas particles in this cell. Note that it means we are missing
   *gas particles in other cells.
   *
   * TODO (Darwin): This will be improved in the future with a proper self/pair
   *task search.
   */
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  /* Loop over all particles to find the neighbours within r_acc. Then,
     compute all quantities you need.  */
  for (int j = 0; j < count; j++) {

    /* Get a handle on the part */
    struct part *restrict pj = &parts[j];
    struct xpart *restrict xpj = &xparts[j];

    /* Ignore inhibited particles */
    if (part_is_inhibited(pj, e)) continue;

    /* Compute the quantities required to later decide to form a sink or not. */
    sink_prepare_part_sink_formation_gas_criteria(e, pi, xpi, pj, xpj, cosmo,
                                                  sink_props);
  } /* End of gas neighbour loop */

  /* Check that we are not forming a sink in the accretion radius of another
     one. The new sink may be swallowed by the older one.) */

  /* For the sinks, we can loop over all sinks in the space. This is an
     O(N_part_eligible*N_sink) search. We assume that N_sink < N_part, which
     make this brute force search feasible.
   *
   * TODO: In the future, we can optimise by adding a self/pair tasks */
  const int scount = s->nr_sinks;
  struct sink *restrict sinks = s->sinks;

  for (int j = 0; j < scount; j++) {

    /* Get a hold of the ith sinks in ci. */
    struct sink *restrict sj = &sinks[j];

    /* Ignore inhibited particles */
    if (sink_is_inhibited(sj, e)) continue;

    /* Compute the quantities required to later decide to form a sink or not. */
    sink_prepare_part_sink_formation_sink_criteria(
        e, pi, xpi, sj, with_cosmology, cosmo, sink_props, e->time);

  } /* End of sink neighbour loop */
}
