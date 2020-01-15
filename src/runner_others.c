/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "chemistry.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "logger.h"
#include "pressure_floor.h"
#include "space.h"
#include "star_formation.h"
#include "star_formation_logger.h"
#include "stars.h"
#include "timers.h"
#include "timestep_limiter.h"
#include "tracers.h"

/**
 * @brief Calculate gravity acceleration from external potential
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_grav_external(struct runner *r, struct cell *c, int timer) {

  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;
  const struct engine *e = r->e;
  const struct external_potential *potential = e->external_potential;
  const struct phys_const *constants = e->physical_constants;
  const double time = r->e->time;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_grav_external(r, c->progeny[k], 0);
  } else {

    /* Loop over the gparts in this cell. */
    for (int i = 0; i < gcount; i++) {

      /* Get a direct pointer on the part. */
      struct gpart *restrict gp = &gparts[i];

      /* Is this part within the time step? */
      if (gpart_is_active(gp, e)) {
        external_gravity_acceleration(time, potential, constants, gp);
      }
    }
  }

  if (timer) TIMER_TOC(timer_dograv_external);
}

/**
 * @brief Calculate gravity accelerations from the periodic mesh
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_grav_mesh(struct runner *r, struct cell *c, int timer) {

  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;
  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (!e->s->periodic) error("Calling mesh forces in non-periodic mode.");
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_grav_mesh(r, c->progeny[k], 0);
  } else {

    /* Get the forces from the gravity mesh */
    pm_mesh_interpolate_forces(e->mesh, e, gparts, gcount);
  }

  if (timer) TIMER_TOC(timer_dograv_mesh);
}

/**
 * @brief Calculate change in thermal state of particles induced
 * by radiative cooling and heating.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_cooling(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cooling_function_data *cooling_func = e->cooling_func;
  const struct phys_const *constants = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct entropy_floor_properties *entropy_floor_props = e->entropy_floor;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_cooling(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      struct part *restrict p = &parts[i];
      struct xpart *restrict xp = &xparts[i];

      if (part_is_active(p, e)) {

        double dt_cool, dt_therm;
        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(ti_current - 1, p->time_bin);

          dt_cool =
              cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
          dt_therm = cosmology_get_therm_kick_factor(e->cosmology, ti_begin,
                                                     ti_begin + ti_step);

        } else {
          dt_cool = get_timestep(p->time_bin, time_base);
          dt_therm = get_timestep(p->time_bin, time_base);
        }

        /* Let's cool ! */
        cooling_cool_part(constants, us, cosmo, hydro_props,
                          entropy_floor_props, cooling_func, p, xp, dt_cool,
                          dt_therm);
      }
    }
  }

  if (timer) TIMER_TOC(timer_do_cooling);
}

/**
 *
 */
void runner_do_star_formation(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct star_formation *sf_props = e->star_formation;
  const struct phys_const *phys_const = e->physical_constants;
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const struct hydro_props *restrict hydro_props = e->hydro_properties;
  const struct unit_system *restrict us = e->internal_units;
  struct cooling_function_data *restrict cooling = e->cooling_func;
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;
  const int current_stars_count = c->stars.count;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID)
    error("Running star formation task on a foreign node!");
#endif

  /* Anything to do here? */
  if (c->hydro.count == 0 || !cell_is_active_hydro(c, e)) {
    star_formation_logger_log_inactive_cell(&c->stars.sfh);
    return;
  }

  /* Reset the SFR */
  star_formation_logger_init(&c->stars.sfh);

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        /* Load the child cell */
        struct cell *restrict cp = c->progeny[k];

        /* Do the recursion */
        runner_do_star_formation(r, cp, 0);

        /* Update current cell using child cells */
        star_formation_logger_add(&c->stars.sfh, &cp->stars.sfh);
      }
  } else {

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Only work on active particles */
      if (part_is_active(p, e)) {

        /* Is this particle star forming? */
        if (star_formation_is_star_forming(p, xp, sf_props, phys_const, cosmo,
                                           hydro_props, us, cooling,
                                           entropy_floor)) {

          /* Time-step size for this particle */
          double dt_star;
          if (with_cosmology) {
            const integertime_t ti_step = get_integer_timestep(p->time_bin);
            const integertime_t ti_begin =
                get_integer_time_begin(ti_current - 1, p->time_bin);

            dt_star =
                cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);

          } else {
            dt_star = get_timestep(p->time_bin, time_base);
          }

          /* Compute the SF rate of the particle */
          star_formation_compute_SFR(p, xp, sf_props, phys_const, hydro_props,
                                     cosmo, dt_star);

          /* Add the SFR and SFR*dt to the SFH struct of this cell */
          star_formation_logger_log_active_part(p, xp, &c->stars.sfh, dt_star);

          /* Are we forming a star particle from this SF rate? */
          if (star_formation_should_convert_to_star(p, xp, sf_props, e,
                                                    dt_star)) {

#ifdef WITH_LOGGER
            /* Write the particle */
            /* Logs all the fields request by the user */
            // TODO select only the requested fields
            logger_log_part(e->logger, p,
                            logger_mask_data[logger_x].mask |
                                logger_mask_data[logger_v].mask |
                                logger_mask_data[logger_a].mask |
                                logger_mask_data[logger_u].mask |
                                logger_mask_data[logger_h].mask |
                                logger_mask_data[logger_rho].mask |
                                logger_mask_data[logger_consts].mask |
                                logger_mask_data[logger_special_flags].mask,
                            &xp->logger_data.last_offset,
                            /* special flags */ swift_type_stars);
#endif

            /* Convert the gas particle to a star particle */
            struct spart *sp = cell_convert_part_to_spart(e, c, p, xp);

            /* Did we get a star? (Or did we run out of spare ones?) */
            if (sp != NULL) {

              /* message("We formed a star id=%lld cellID=%d", sp->id,
               * c->cellID); */

              /* Copy the properties of the gas particle to the star particle */
              star_formation_copy_properties(p, xp, sp, e, sf_props, cosmo,
                                             with_cosmology, phys_const,
                                             hydro_props, us, cooling);

              /* Update the Star formation history */
              star_formation_logger_log_new_spart(sp, &c->stars.sfh);

#ifdef WITH_LOGGER
              /* Copy the properties back to the stellar particle */
              sp->logger_data = xp->logger_data;

              /* Write the s-particle */
              logger_log_spart(e->logger, sp,
                               logger_mask_data[logger_x].mask |
                                   logger_mask_data[logger_v].mask |
                                   logger_mask_data[logger_consts].mask,
                               &sp->logger_data.last_offset,
                               /* special flags */ 0);

              /* Set counter back to zero */
              sp->logger_data.steps_since_last_output = 0;
#endif
            }
          }

        } else { /* Are we not star-forming? */

          /* Update the particle to flag it as not star-forming */
          star_formation_update_part_not_SFR(p, xp, e, sf_props,
                                             with_cosmology);

        } /* Not Star-forming? */

      } else { /* is active? */

        /* Check if the particle is not inhibited */
        if (!part_is_inhibited(p, e)) {
          star_formation_logger_log_inactive_part(p, xp, &c->stars.sfh);
        }
      }
    } /* Loop over particles */
  }

  /* If we formed any stars, the star sorts are now invalid. We need to
   * re-compute them. */
  if (with_feedback && (c == c->top) &&
      (current_stars_count != c->stars.count)) {
    cell_set_star_resort_flag(c);
  }

  if (timer) TIMER_TOC(timer_do_star_formation);
}

/**
 * @brief End the hydro force calculation of all active particles in a cell
 * by multiplying the acccelerations by the relevant constants
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_end_hydro_force(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_end_hydro_force(r, c->progeny[k], 0);
  } else {

    const struct cosmology *cosmo = e->cosmology;
    const int count = c->hydro.count;
    struct part *restrict parts = c->hydro.parts;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];

      if (part_is_active(p, e)) {

        /* Finish the force loop */
        hydro_end_force(p, cosmo);
        timestep_limiter_end_force(p);
        chemistry_end_force(p, cosmo);

#ifdef SWIFT_BOUNDARY_PARTICLES

        /* Get the ID of the part */
        const long long id = p->id;

        /* Cancel hdyro forces of these particles */
        if (id < SWIFT_BOUNDARY_PARTICLES) {

          /* Don't move ! */
          hydro_reset_acceleration(p);

#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)

          /* Some values need to be reset in the Gizmo case. */
          hydro_prepare_force(p, &c->hydro.xparts[k], cosmo,
                              e->hydro_properties, 0);
#endif
        }
#endif
      }
    }
  }

  if (timer) TIMER_TOC(timer_end_hydro_force);
}

/**
 * @brief End the gravity force calculation of all active particles in a cell
 * by multiplying the acccelerations by the relevant constants
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_end_grav_force(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);
  const int with_black_holes = (e->policy & engine_policy_black_holes);

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_end_grav_force(r, c->progeny[k], 0);
  } else {

    const struct space *s = e->s;
    const int periodic = s->periodic;
    const float G_newton = e->physical_constants->const_newton_G;

    /* Potential normalisation in the case of periodic gravity */
    float potential_normalisation = 0.;
    if (periodic && with_self_gravity) {
      const double volume = s->dim[0] * s->dim[1] * s->dim[2];
      const double r_s = e->mesh->r_s;
      potential_normalisation = 4. * M_PI * e->total_mass * r_s * r_s / volume;
    }

    const int gcount = c->grav.count;
    struct gpart *restrict gparts = c->grav.parts;

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the gpart. */
      struct gpart *restrict gp = &gparts[k];

      if (gpart_is_active(gp, e)) {

        /* Finish the force calculation */
        gravity_end_force(gp, G_newton, potential_normalisation, periodic);

#ifdef SWIFT_MAKE_GRAVITY_GLASS

        /* Negate the gravity forces */
        gp->a_grav[0] *= -1.f;
        gp->a_grav[1] *= -1.f;
        gp->a_grav[2] *= -1.f;
#endif

#ifdef SWIFT_NO_GRAVITY_BELOW_ID

        /* Get the ID of the gpart */
        long long id = 0;
        if (gp->type == swift_type_gas)
          id = e->s->parts[-gp->id_or_neg_offset].id;
        else if (gp->type == swift_type_stars)
          id = e->s->sparts[-gp->id_or_neg_offset].id;
        else if (gp->type == swift_type_black_hole)
          error("Unexisting type");
        else
          id = gp->id_or_neg_offset;

        /* Cancel gravity forces of these particles */
        if (id < SWIFT_NO_GRAVITY_BELOW_ID) {

          /* Don't move ! */
          gp->a_grav[0] = 0.f;
          gp->a_grav[1] = 0.f;
          gp->a_grav[2] = 0.f;
        }
#endif

#ifdef SWIFT_DEBUG_CHECKS
        if ((e->policy & engine_policy_self_gravity) &&
            !(e->policy & engine_policy_black_holes)) {

          /* Let's add a self interaction to simplify the count */
          gp->num_interacted++;

          /* Check that this gpart has interacted with all the other
           * particles (via direct or multipoles) in the box */
          if (gp->num_interacted !=
              e->total_nr_gparts - e->count_inhibited_gparts) {

            /* Get the ID of the gpart */
            long long my_id = 0;
            if (gp->type == swift_type_gas)
              my_id = e->s->parts[-gp->id_or_neg_offset].id;
            else if (gp->type == swift_type_stars)
              my_id = e->s->sparts[-gp->id_or_neg_offset].id;
            else if (gp->type == swift_type_black_hole)
              error("Unexisting type");
            else
              my_id = gp->id_or_neg_offset;

            error(
                "g-particle (id=%lld, type=%s) did not interact "
                "gravitationally with all other gparts "
                "gp->num_interacted=%lld, total_gparts=%lld (local "
                "num_gparts=%zd inhibited_gparts=%lld)",
                my_id, part_type_names[gp->type], gp->num_interacted,
                e->total_nr_gparts, e->s->nr_gparts, e->count_inhibited_gparts);
          }
        }
#endif

        /* Deal with black holes' need of potentials */
        if (with_black_holes && gp->type == swift_type_black_hole) {
          const size_t offset = -gp->id_or_neg_offset;
          black_holes_store_potential_in_bpart(&s->bparts[offset], gp);
        }
        if (with_black_holes && gp->type == swift_type_gas) {
          const size_t offset = -gp->id_or_neg_offset;
          black_holes_store_potential_in_part(
              &s->parts[offset].black_holes_data, gp);
        }
      }
    }
  }
  if (timer) TIMER_TOC(timer_end_grav_force);
}

/**
 * @brief Write the required particles through the logger.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_logger(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_LOGGER
  TIMER_TIC;

  const struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e) &&
      !cell_is_active_stars(c, e))
    return;

  /* Recurse? Avoid spending too much time in useless cells. */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_logger(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be log */
      if (part_is_active(p, e)) {

        if (logger_should_write(&xp->logger_data, e->logger)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          logger_log_part(e->logger, p,
                          logger_mask_data[logger_x].mask |
                              logger_mask_data[logger_v].mask |
                              logger_mask_data[logger_a].mask |
                              logger_mask_data[logger_u].mask |
                              logger_mask_data[logger_h].mask |
                              logger_mask_data[logger_rho].mask |
                              logger_mask_data[logger_consts].mask,
                          &xp->logger_data.last_offset,
                          /* special flags */ 0);

          /* Set counter back to zero */
          xp->logger_data.steps_since_last_output = 0;
        } else
          /* Update counter */
          xp->logger_data.steps_since_last_output += 1;
      }
    }

    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* Write only the dark matter particles */
      if (gp->type != swift_type_dark_matter) continue;

      /* If particle needs to be log */
      if (gpart_is_active(gp, e)) {

        if (logger_should_write(&gp->logger_data, e->logger)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          logger_log_gpart(e->logger, gp,
                           logger_mask_data[logger_x].mask |
                               logger_mask_data[logger_v].mask |
                               logger_mask_data[logger_a].mask |
                               logger_mask_data[logger_consts].mask,
                           &gp->logger_data.last_offset,
                           /* Special flags */ 0);

          /* Set counter back to zero */
          gp->logger_data.steps_since_last_output = 0;
        } else
          /* Update counter */
          gp->logger_data.steps_since_last_output += 1;
      }
    }

    /* Loop over the sparts in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the part. */
      struct spart *restrict sp = &sparts[k];

      /* If particle needs to be log */
      if (spart_is_active(sp, e)) {

        if (logger_should_write(&sp->logger_data, e->logger)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          logger_log_spart(e->logger, sp,
                           logger_mask_data[logger_x].mask |
                               logger_mask_data[logger_v].mask |
                               logger_mask_data[logger_consts].mask,
                           &sp->logger_data.last_offset,
                           /* Special flags */ 0);

          /* Set counter back to zero */
          sp->logger_data.steps_since_last_output = 0;
        } else
          /* Update counter */
          sp->logger_data.steps_since_last_output += 1;
      }
    }
  }

  if (timer) TIMER_TOC(timer_logger);

#else
  error("Logger disabled, please enable it during configuration");
#endif
}

/**
 * @brief Recursively search for FOF groups in a single cell.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_fof_self(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_FOF

  TIMER_TIC;

  const struct engine *e = r->e;
  struct space *s = e->s;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const struct gpart *const gparts = s->gparts;
  const double search_r2 = e->fof_properties->l_x2;

  rec_fof_search_self(e->fof_properties, dim, search_r2, periodic, gparts, c);

  if (timer) TIMER_TOC(timer_fof_self);

#else
  error("SWIFT was not compiled with FOF enabled!");
#endif
}

/**
 * @brief Recursively search for FOF groups between a pair of cells.
 *
 * @param r runner task
 * @param ci cell i
 * @param cj cell j
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_fof_pair(struct runner *r, struct cell *ci, struct cell *cj,
                        int timer) {

#ifdef WITH_FOF

  TIMER_TIC;

  const struct engine *e = r->e;
  struct space *s = e->s;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const struct gpart *const gparts = s->gparts;
  const double search_r2 = e->fof_properties->l_x2;

  rec_fof_search_pair(e->fof_properties, dim, search_r2, periodic, gparts, ci,
                      cj);

  if (timer) TIMER_TOC(timer_fof_pair);
#else
  error("SWIFT was not compiled with FOF enabled!");
#endif
}
