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
 *   2. Apply a TDE accretion rate to the BH subgrid mass and energy reservoir.
 *   3. Nibble the corresponding mass from the nearest star particle.
 *
 * Note: this is an O(N_BH * N_stars) brute-force implementation. It is also
 * not thread-safe when multiple BH cells are processed concurrently (mitigated
 * by a space-level lock around star mass writes).
 *
 * NOTE: star particles are not guaranteed to be fully drifted to the current
 * time at this point in the task graph. A proper fix requires a dedicated
 * BH-star pair task with drift_spart dependencies.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_bh_stellar_accretion(struct runner *r, struct cell *c,
                                    int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;
  const struct cosmology *cosmo = e->cosmology;
  const struct unit_system *us = e->internal_units;
  const struct black_holes_props *props = e->black_holes_properties;
  const struct phys_const *constants = e->physical_constants;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int periodic = s->periodic;

  /* Minimum star mass allowed after nibbling: 50% of mean baryon particle mass. */
  const double min_star_mass_for_nibbling =
      0.5 * s->initial_mean_mass_particles[swift_type_gas];

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

  /* Half-life for TDE mass loss: star loses 50% of available mass per Gyr. */
  const double gyr_in_cgs = 3.15576e16; /* 1 Gyr in seconds */
  const double t_half = gyr_in_cgs / us->UnitTime_in_cgs;
  const double ln2_over_t_half = log(2.0) / t_half;

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

    /* --- Single pass: accumulate total stellar mass and find nearest star --- */
    /* --- Also get mean velocity --- */
    double M_total = 0.0;
    double nearest_r2 = aperture_comoving2;
    struct spart *nearest_sp = NULL;
    int star_count = 0;
    double v_mean_x = 0.0;
    double v_mean_y = 0.0;
    double v_mean_z = 0.0;
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

      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < aperture_comoving2) {
        M_total += sp->mass;
        v_mean_x += sp->v[0] - bp->v_full[0];
			  v_mean_y += sp->v[1] - bp->v_full[1];
        v_mean_z += sp->v[2] - bp->v_full[2];
        star_count++;
      }  
      if (r2 < nearest_r2) {
        nearest_r2 = r2;
        nearest_sp = sp;
      }
    }
    /* finished pass, calculate mean velocity */
    v_mean_x /= star_count;
    v_mean_y /= star_count;
    v_mean_z /= star_count;

    /* --- Single pass 2: get density and velocity dispersion profiles by binning --- */
    /* Create log-spaced bins */
    int n_bins = 50;        // This is a variable to decide on
    double start = log10(dist_inner_internal);
    double stop = log10(dist_threshold_internal);
    double step = (stop - start) / n_bins;
    double bin_edges[n_bins + 1];
    double bin_centres[n_bins];
    for (int i = 0; i < n_bins + 1; i++){
        bin_edges[i] = pow(10, (start + step*i));
    }

    /* Initialize arrays */
    double mass_per_shell[n_bins];
    double density_per_shell[n_bins];
    double v_rel_x = 0.0;
    double v_rel_y = 0.0;
    double v_rel_z = 0.0;
    double count_per_shell[n_bins]; 
    double sig_per_shell_x[n_bins];
    double sig_per_shell_y[n_bins];
    double sig_per_shell_z[n_bins];
    double sig_per_shell[n_bins];
    double sig_total_x = 0.0;
    double sig_total_y = 0.0;
    double sig_total_z = 0.0;
    double sig_total = 0.0;
    for (int i = 0; i < n_bins; i++){
      bin_centres[i] = (bin_edges[i] + bin_edges[i+1]) / 2;
      count_per_shell[i] = 0.0;
      mass_per_shell[i] = 0.0;
      sig_per_shell_x[i] = 0.0;
      sig_per_shell_y[i] = 0.0;
      sig_per_shell_z[i] = 0.0;
      sig_per_shell[i] = 0.0;
    }

    /* Do second pass */
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

      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r1 = sqrt(r2)
      if (r2 < aperture_comoving2) {
			  v_rel_x = sp->v[0] - bp->v_full[0] - v_mean_x;
			  v_rel_y = sp->v[1] - bp->v_full[1] - v_mean_y;
			  v_rel_z = sp->v[2] - bp->v_full[2] - v_mean_z;

			  /* Binwise calculate velocity dispersion, mass-weighted */
        /* Bin stellar masses */
			  for (int j = 0; j < n_bins; j++){
			    if (r1 >= bin_edges[j] && r1 < bin_edges[j+1]){
            count_per_shell[j]++;
            mass_per_shell[j] += sp->mass;
            sig_per_shell_x[j] += v_rel_x * v_rel_x * sp->mass;
            sig_per_shell_y[j] += v_rel_y * v_rel_y * sp->mass;
            sig_per_shell_z[j] += v_rel_z * v_rel_z * sp->mass;
            sig_total_x += v_rel_x * v_rel_x * sp->mass;
            sig_total_y += v_rel_y * v_rel_y * sp->mass;
            sig_total_z += v_rel_z * v_rel_z * sp->mass;
			  	}
			  }
      }
    }

    /* With all particles binned, we loop over bins for density and final velocity dispersion */
    for (int i = 0; i < n_bins; i++){
      double shellvolume = (4.0/3.0) * M_PI * (pow(bin_edges[i+1], 3) - pow(bin_edges[i], 3));
    density_per_shell[i] = mass_per_shell[i] / shellvolume;

    if (count_per_shell[i] < 5){
      sig_per_shell[i] = NAN;
    } else{
      sig_per_shell_x[i] /= mass_per_shell[i];
      sig_per_shell_y[i] /= mass_per_shell[i];
      sig_per_shell_z[i] /= mass_per_shell[i];
      sig_per_shell[i] = sqrt((sig_per_shell_x[i] + sig_per_shell_y[i] + sig_per_shell_z[i]) / 3.0);
    }
    }

    /* Compute total velocity dispersion */
    sig_total_x /= M_total;
    sig_total_y /= M_total;
    sig_total_z /= M_total;
    sig_total = sqrt((sig_total_x + sig_total_y + sig_total_z) / 3.0);


    /* Store physical stellar mass density on the BH. */
    bp->rho_stellar = (M_total > 0.0) ? (float)(M_total / aperture_volume) : 0.f;

    /* Store stellar density profile, sigma profile and total sigma on the BH. */
    /* Can we do this? (store arrays) */

    /* Nothing to nibble if no star found. */
    if (nearest_sp == NULL) continue;

    /* Available mass on the nearest star above the nibbling floor. */
    const double available_mass =
        (double)nearest_sp->mass - min_star_mass_for_nibbling;
    if (available_mass <= 0.0) continue;

    /* Star loses 50% of available mass per Gyr (exponential decay,
     * half-life = 1 Gyr): star_mass_loss_rate = ln(2) / t_half * available_mass.
     * The BH gains less due to radiation: bh_accretion_rate = (1 - epsilon_r)
     * * star_mass_loss_rate. */
    const double star_mass_loss_rate = ln2_over_t_half * available_mass;

    /* Mass changes this timestep. */
    const double star_mass_loss = star_mass_loss_rate * dt;
    const double new_star_mass = (double)nearest_sp->mass - star_mass_loss;

    /* Check minimum mass before modifying any particles. */
    if (new_star_mass < min_star_mass_for_nibbling) {
      warning(
          "TDE nibbling would reduce star particle %lld (mass=%g) below "
          "minimum mass threshold. BH ID=%lld (mass=%g). Skipping.",
          nearest_sp->id, (double)nearest_sp->mass, bp->id, (double)bp->mass);
      continue;
    }

    /* Update BH subgrid mass and energy reservoir. */
    const double bh_mass_gain =
        black_holes_do_tde_accretion(bp, props, constants, star_mass_loss_rate, dt);

    message(
        "BH (ID %lld) z=%.4f  rho_stellar=%g (internal)  "
        "bh_mass_gain=%g (internal)",
        bp->id, cosmo->z, (double)bp->rho_stellar, bh_mass_gain);

    /* Lock the space to prevent concurrent writes from other BH cells
     * being processed simultaneously on different threads. */
    lock_lock(&s->lock);
    nearest_sp->mass = (float)new_star_mass;
    nearest_sp->gpart->mass = (float)new_star_mass;
    nearest_sp->mass_lost_to_tde += (float)star_mass_loss;
    if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space.");

    /* Update BH velocity to conserve momentum of the accreted mass,
     * mirroring the gas nibbling momentum update. */
    const double bp_mass_orig = (double)bp->mass;
    const double new_bp_mass = bp_mass_orig + bh_mass_gain;
    bp->v[0] = (float)((bp_mass_orig * bp->v[0] +
                        bh_mass_gain * nearest_sp->v[0]) / new_bp_mass);
    bp->v[1] = (float)((bp_mass_orig * bp->v[1] +
                        bh_mass_gain * nearest_sp->v[1]) / new_bp_mass);
    bp->v[2] = (float)((bp_mass_orig * bp->v[2] +
                        bh_mass_gain * nearest_sp->v[2]) / new_bp_mass);

    /* Add the net accreted mass (excluding radiation) to the BH dynamical
     * mass, consistent with how gas nibbling updates bp->mass. */
    bp->mass = (float)new_bp_mass;
    bp->gpart->mass = (float)new_bp_mass;
    bp->mass_gained_from_tde += (float)bh_mass_gain;
  }
}
