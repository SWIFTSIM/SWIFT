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
  error("MPI is not implemented yet for sink particles.");
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
        error("MPI is not implemented yet for sink particles.");
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
  error("MPI is not implemented yet for sink particles.");
#endif

  struct sink *cell_sinks = c->sinks.parts;

  /* Early abort?
   * (We only want cells for which we drifted the sink as these are
   * the only ones that could have sink particles that have been flagged
   * for swallowing) */
  if (c->sinks.count == 0 || c->sinks.ti_old_part != e->ti_current) {
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

            // message("sink %lld swallowing sink particle %lld", sp->id,
            // cell_sp->id);

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
        error("MPI is not implemented yet for sink particles.");
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

#ifdef SWIFT_DEBUG_CHECKS
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

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active sink */
  if (cell_is_active_sinks(cj, e)) runner_do_sinks_sink_swallow(r, ci, timer);
  if (cell_is_active_sinks(ci, e)) runner_do_sinks_sink_swallow(r, cj, timer);
}

/**
 * @brief Compute the energies (kinetic, potential, etc ) of the gas particle
 * #p and all quantities required for the formation of a sink.
 *
 * Note: This function iterates over gas particles and sink particles.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param p The #part.
 * @param xp The #xpart data of the particle #p.
 */
void runner_do_prepare_part_sink_formation(struct runner *r, struct cell *c,
                                           struct part *restrict p,
                                           struct xpart *restrict xp) {
  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct sink_props *sink_props = e->sink_properties;
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  /* Loop over all particles to find the neighbours within r_acc. Then,
     compute all quantities you need.  */
  for (int i = 0; i < count; i++) {

    /*Get a handle on the part */
    struct part *restrict pi = &parts[i];
    struct xpart *restrict xpi = &xparts[i];

    /* Compute the quantities required to later decide to form a sink or not. */
    sink_prepare_part_sink_formation_gas_criteria(e, p, xp, pi, xpi, cosmo,
                                                  sink_props);
  } /* End of gas neighbour loop */

  /* Shall we reset the values of the energies for the next timestep? No, it is
     done in cell_drift.c and space_init.c, for active particles. The
     potential is set in runner_others.c->runner_do_end_grav_force() */

  /* Check that we are not forming a sink in the accretion radius of another
     one. The new sink may be swallowed by the older one.) */
  const int scount = c->sinks.count;
  struct sink *restrict sinks = c->sinks.parts;

  for (int i = 0; i < scount; i++) {

    /* Get a hold of the ith sinks in ci. */
    struct sink *restrict si = &sinks[i];

    /* Compute the quantities required to later decide to form a sink or not. */
    sink_prepare_part_sink_formation_sink_criteria(e, p, xp, si, cosmo,
                                                   sink_props);

  } /* End of sink neighbour loop */
}
