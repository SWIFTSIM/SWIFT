/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include "black_holes.h"
#include "cell.h"
#include "engine.h"
#include "periodic.h"
#include "timers.h"

/**
 * @brief Process all the gas particles in a cell that have been flagged for
 * swallowing by a black hole.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of black holes for the particle with the corresponding
 * ID. If found, the BH swallows the gas particle and the gas particle is
 * removed. If the cell is local, we may be looking for a foreign BH, in which
 * case, we do not update the BH (that will be done on its node) but just remove
 * the gas particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_gas_swallow(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;
  const struct black_holes_props *props = e->black_holes_properties;
  const int use_nibbling = props->use_nibbling;

  struct bpart *bparts = s->bparts;
  const size_t nr_bpart = s->nr_bparts;
#ifdef WITH_MPI
  struct bpart *bparts_foreign = s->bparts_foreign;
  const size_t nr_bparts_foreign = s->nr_bparts_foreign;
#endif

  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;

  /* Nothing to do here if the cell is foreign and we are nibbling */
  if (c->nodeID != e->nodeID && use_nibbling) {
    return;
  }

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

        runner_do_gas_swallow(r, cp, 0);
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

      /* Update mass of associated gpart, to reflect potential changes from
       * nibbling. In this case, we are already done. */
      if (use_nibbling) {
        p->gpart->mass = hydro_get_mass(p);
        continue;
      }

      /* Get the ID of the black holes that will swallow this part */
      const long long swallow_id =
          black_holes_get_part_swallow_id(&p->black_holes_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

#ifdef SWIFT_DEBUG_CHECKS
        if (p->ti_drift != e->ti_current)
          error("Trying to swallow an un-drifted particle.");
#endif

        /* ID of the BH swallowing this particle */
        const long long BH_id = swallow_id;

        /* Have we found this particle's BH already? */
        int found = 0;

        /* Let's look for the hungry black hole in the local list */
        for (size_t i = 0; i < nr_bpart; ++i) {

          /* Get a handle on the bpart. */
          struct bpart *bp = &bparts[i];

          if (bp->id == BH_id) {

            /* Lock the space as we are going to work directly on the bpart list
             */
            lock_lock(&s->lock);

            /* Swallow the gas particle (i.e. update the BH properties) */
            black_holes_swallow_part(bp, p, xp, e->cosmology);

            /* Release the space as we are done updating the bpart */
            if (lock_unlock(&s->lock) != 0)
              error("Failed to unlock the space.");

            /* If the gas particle is local, remove it */
            if (c->nodeID == e->nodeID) {

              message("BH %lld removing gas particle %lld", bp->id, p->id);

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
            black_holes_mark_part_as_swallowed(&p->black_holes_data);

            found = 1;
            break;
          }

        } /* Loop over local BHs */

#ifdef WITH_MPI

        /* We could also be in the case of a local gas particle being
         * swallowed by a foreign BH. In this case, we won't update the
         * BH but just remove the particle from the local list. */
        if (c->nodeID == e->nodeID && !found) {

          /* Let's look for the foreign hungry black hole */
          for (size_t i = 0; i < nr_bparts_foreign; ++i) {

            /* Get a handle on the bpart. */
            struct bpart *bp = &bparts_foreign[i];

            if (bp->id == BH_id) {

              message("BH %lld removing gas particle %lld (foreign BH case)",
                      bp->id, p->id);

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
          } /* Loop over foreign BHs */
        } /* Is the cell local? */
#endif

        /* If we have a local particle, we must have found the BH in one
         * of our list of black holes. */
        if (c->nodeID == e->nodeID && !found) {
          error("Gas particle %lld could not find BH %lld to be swallowed",
                p->id, swallow_id);
        }
      } /* Part was flagged for swallowing */
    } /* Loop over the parts */
  } /* Cell is not split */
}

/**
 * @brief Processing of gas particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_gas_swallow_self(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != r->e->nodeID) error("Running self task on foreign node");
  if (!cell_is_active_black_holes(c, r->e))
    error("Running self task on inactive cell");
#endif

  runner_do_gas_swallow(r, c, timer);
}

/**
 * @brief Processing of gas particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_gas_swallow_pair(struct runner *r, struct cell *ci,
                                struct cell *cj, int timer) {

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active BH */
  if (cell_is_active_black_holes(cj, e)) runner_do_gas_swallow(r, ci, timer);
  if (cell_is_active_black_holes(ci, e)) runner_do_gas_swallow(r, cj, timer);
}

/**
 * @brief Process all the BH particles in a cell that have been flagged for
 * swallowing by a black hole.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of black holes for the particle with the corresponding
 * ID. If found, the BH swallows the BH particle and the BH particle is
 * removed. If the cell is local, we may be looking for a foreign BH, in which
 * case, we do not update the BH (that will be done on its node) but just remove
 * the BH particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_bh_swallow(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct black_holes_props *props = e->black_holes_properties;
  const int use_nibbling = props->use_nibbling;

  struct bpart *bparts = s->bparts;
  const size_t nr_bpart = s->nr_bparts;
#ifdef WITH_MPI
  struct bpart *bparts_foreign = s->bparts_foreign;
  const size_t nr_bparts_foreign = s->nr_bparts_foreign;
#endif

  struct bpart *cell_bparts = c->black_holes.parts;

  /* Early abort?
   * (We only want cells for which we drifted the BH as these are
   * the only ones that could have BH particles that have been flagged
   * for swallowing) */
  if (c->black_holes.count == 0 ||
      c->black_holes.ti_old_part != e->ti_current) {
    return;
  }

  /* Loop over the progeny ? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        runner_do_bh_swallow(r, cp, 0);
      }
    }
  } else {

    /* Loop over all the BH particles in the cell
     * Note that the cell (and hence the bparts) may be local or foreign. */
    const size_t nr_cell_bparts = c->black_holes.count;

    for (size_t k = 0; k < nr_cell_bparts; k++) {

      /* Get a handle on the part. */
      struct bpart *const cell_bp = &cell_bparts[k];

      /* Ignore inhibited particles (they have already been removed!) */
      if (bpart_is_inhibited(cell_bp, e)) continue;

      /* Update mass of associated gpart, to reflect potential changes from
       * nibbling. */
      if (use_nibbling && c->nodeID == e->nodeID) {
        cell_bp->gpart->mass = cell_bp->mass;
        cell_bp->gpart->v_full[0] = cell_bp->v[0];
        cell_bp->gpart->v_full[1] = cell_bp->v[1];
        cell_bp->gpart->v_full[2] = cell_bp->v[2];
      }

      /* Get the ID of the black holes that will swallow this bpart */
      const long long swallow_id =
          black_holes_get_bpart_swallow_id(&cell_bp->merger_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

#ifdef SWIFT_DEBUG_CHECKS
        if (cell_bp->ti_drift != e->ti_current)
          error("Trying to swallow an un-drifted particle.");
#endif

        /* ID of the BH swallowing this particle */
        const long long BH_id = swallow_id;

        /* Have we found this particle's BH already? */
        int found = 0;

        /* Let's look for the hungry black hole in the local list */
        for (size_t i = 0; i < nr_bpart; ++i) {

          /* Get a handle on the bpart. */
          struct bpart *bp = &bparts[i];

          if (bp->id == BH_id) {

            /* Is the swallowing BH itself flagged for swallowing by
               another BH? */
            if (black_holes_get_bpart_swallow_id(&bp->merger_data) != -1) {

              /* Pretend it was found and abort */
              black_holes_mark_bpart_as_not_swallowed(&cell_bp->merger_data);
              found = 1;
              break;
            }

            /* Lock the space as we are going to work directly on the
             * space's bpart list */
            lock_lock(&s->lock);

            /* Swallow the BH particle (i.e. update the swallowing BH
             * properties with the properties of cell_bp) */
            black_holes_swallow_bpart(bp, cell_bp, e->cosmology, e->time,
                                      with_cosmology, props,
                                      e->physical_constants);

            /* Release the space as we are done updating the bpart */
            if (lock_unlock(&s->lock) != 0)
              error("Failed to unlock the space.");

            message("BH %lld swallowing BH particle %lld", bp->id, cell_bp->id);

            /* If the BH particle is local, remove it */
            if (c->nodeID == e->nodeID) {

              message("BH %lld removing BH particle %lld", bp->id, cell_bp->id);

              /* Finally, remove the BH particle from the system
               * Recall that the gpart associated with it is also removed
               * at the same time. */
              cell_remove_bpart(e, c, cell_bp);
            }

            /* In any case, prevent the particle from being re-swallowed */
            black_holes_mark_bpart_as_merged(&cell_bp->merger_data);

            found = 1;
            break;
          }

        } /* Loop over local BHs */

#ifdef WITH_MPI

        /* We could also be in the case of a local BH particle being
         * swallowed by a foreign BH. In this case, we won't update the
         * foreign BH but just remove the particle from the local list. */
        if (c->nodeID == e->nodeID && !found) {

          /* Let's look for the foreign hungry black hole */
          for (size_t i = 0; i < nr_bparts_foreign; ++i) {

            /* Get a handle on the bpart. */
            struct bpart *bp = &bparts_foreign[i];

            if (bp->id == BH_id) {

              /* Is the swallowing BH itself flagged for swallowing by
                 another BH? */
              if (black_holes_get_bpart_swallow_id(&bp->merger_data) != -1) {

                /* Pretend it was found and abort */
                black_holes_mark_bpart_as_not_swallowed(&cell_bp->merger_data);
                found = 1;
                break;
              }

              message("BH %lld removing BH particle %lld (foreign BH case)",
                      bp->id, cell_bp->id);

              /* Finally, remove the gas particle from the system */
              cell_remove_bpart(e, c, cell_bp);

              found = 1;
              break;
            }
          } /* Loop over foreign BHs */
        } /* Is the cell local? */
#endif

        /* If we have a local particle, we must have found the BH in one
         * of our list of black holes. */
        if (c->nodeID == e->nodeID && !found) {
          error("BH particle %lld could not find BH %lld to be swallowed",
                cell_bp->id, swallow_id);
        }

      } /* Part was flagged for swallowing */
    } /* Loop over the parts */
  } /* Cell is not split */
}

/**
 * @brief Processing of bh particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_bh_swallow_self(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != r->e->nodeID) error("Running self task on foreign node");
  if (!cell_is_active_black_holes(c, r->e))
    error("Running self task on inactive cell");
#endif

  runner_do_bh_swallow(r, c, timer);
}

/**
 * @brief Processing of bh particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_bh_swallow_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer) {

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active BH */
  if (cell_is_active_black_holes(cj, e)) runner_do_bh_swallow(r, ci, timer);
  if (cell_is_active_black_holes(ci, e)) runner_do_bh_swallow(r, cj, timer);
}

/**
 * @brief Brute-force stellar accretion (TDE) for all active BHs in a cell.
 *
 * For each active BH, loops over all star particles in the simulation to:
 *   1. Compute the total stellar mass density within a 1 kpc physical aperture.
 *   2. Apply a TDE accretion rate (currently a placeholder constant) to the BH
 *      subgrid mass.
 *   3. Nibble the corresponding mass from the star particles within the aperture,
 *      proportional to each star's share of the total stellar mass.
 *
 * Note: this is an O(N_BH * N_stars) brute-force implementation. It is also
 * not thread-safe when multiple BH cells are processed concurrently.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_bh_stellar_accretion(struct runner *r, struct cell *c,
                                    int timer) {

  struct engine *e = r->e;
  const struct space *s = e->s;
  const struct cosmology *cosmo = e->cosmology;
  const struct unit_system *us = e->internal_units;
  const struct black_holes_props *props = e->black_holes_properties;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int periodic = s->periodic;

  /* Anything to do here? */
  if (c->black_holes.count == 0) return;
  if (!cell_is_active_black_holes(c, e)) return;

  /* Recurse over sub-cells. */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_bh_stellar_accretion(r, c->progeny[k], 0);
    return;
  }

  struct bpart *restrict bparts = c->black_holes.parts;
  const int count = c->black_holes.count;

  /* Physical aperture radius: 1 kpc in internal length units. */
  const double kpc_in_cm = 3.08567758e21;
  const double aperture_phys = kpc_in_cm / us->UnitLength_in_cgs;

  /* Comoving aperture threshold (positions are stored in comoving coords). */
  const double aperture_comoving = aperture_phys / cosmo->a;
  const double aperture_comoving2 = aperture_comoving * aperture_comoving;

  /* Physical volume of the aperture sphere. */
  const double aperture_volume =
      (4.0 / 3.0) * M_PI * aperture_phys * aperture_phys * aperture_phys;

  /* Radiative efficiency (same as gas accretion). */
  const double epsilon_r = props->epsilon_r;
  const double excess_fraction = 1.0 / (1.0 - epsilon_r);

  for (int i = 0; i < count; i++) {

    struct bpart *restrict bp = &bparts[i];
    if (!bpart_is_active(bp, e)) continue;

    /* Get this BH's timestep. */
    double dt;
    if (with_cosmology) {
      const integertime_t ti_step = get_integer_timestep(bp->time_bin);
      const integertime_t ti_begin =
          get_integer_time_begin(e->ti_current - 1, bp->time_bin);
      dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
    } else {
      dt = get_timestep(bp->time_bin, e->time_base);
    }

    /* --- First pass: accumulate total stellar mass within aperture --- */
    double M_total = 0.0;
    for (size_t j = 0; j < s->nr_sparts; j++) {
      const struct spart *sp = &s->sparts[j];
      if (spart_is_inhibited(sp, e)) continue;

      double dx = bp->x[0] - sp->x[0];
      double dy = bp->x[1] - sp->x[1];
      double dz = bp->x[2] - sp->x[2];

      if (periodic) {
        dx = nearest(dx, s->dim[0]);
        dy = nearest(dy, s->dim[1]);
        dz = nearest(dz, s->dim[2]);
      }

      if (dx * dx + dy * dy + dz * dz < aperture_comoving2)
        M_total += sp->mass;
    }

    /* Store physical stellar mass density on the BH. */
    bp->rho_stellar = (M_total > 0.0) ? (float)(M_total / aperture_volume) : 0.f;

    /* Compute TDE accretion and update BH subgrid mass. */
    const double tde_mass_deficit = black_holes_do_tde_accretion(bp, dt);

    message(
        "BH (ID %lld) z=%.4f  rho_stellar=%g (internal)  "
        "tde_mass_deficit=%g (internal)",
        bp->id, cosmo->z, (double)bp->rho_stellar, tde_mass_deficit);

    /* Nothing to nibble if no stars nearby or zero deficit. */
    if (M_total <= 0.0 || tde_mass_deficit <= 0.0) continue;

    /* --- Second pass: nibble mass from star particles within aperture --- */
    for (size_t j = 0; j < s->nr_sparts; j++) {
      struct spart *sp = &s->sparts[j];
      if (spart_is_inhibited(sp, e)) continue;

      double dx = bp->x[0] - sp->x[0];
      double dy = bp->x[1] - sp->x[1];
      double dz = bp->x[2] - sp->x[2];

      if (periodic) {
        dx = nearest(dx, s->dim[0]);
        dy = nearest(dy, s->dim[1]);
        dz = nearest(dz, s->dim[2]);
      }

      if (dx * dx + dy * dy + dz * dz >= aperture_comoving2) continue;

      /* Each star contributes proportionally to its share of M_total. */
      const double weight = (double)sp->mass / M_total;

      /* Mass removed from this star (includes the radiated fraction). */
      const float star_mass_loss =
          (float)(tde_mass_deficit * weight * excess_fraction);

      if (star_mass_loss >= sp->mass) {
        warning(
            "TDE nibbling would make star particle %lld mass negative "
            "(mass=%g, loss=%g). Skipping.",
            sp->id, (double)sp->mass, (double)star_mass_loss);
        continue;
      }

      sp->mass -= star_mass_loss;
      sp->gpart->mass -= star_mass_loss;
    }

    /* Add the net accreted mass (excluding radiation) to the BH dynamical
     * mass, consistent with how gas nibbling updates bp->mass. */
    bp->mass += (float)tde_mass_deficit;
    bp->gpart->mass += (float)tde_mass_deficit;
  }
}
