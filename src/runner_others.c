/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

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
#include "csds.h"
#include "csds_io.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "fof.h"
#include "forcing.h"
#include "gravity.h"
#include "hydro.h"
#include "potential.h"
#include "pressure_floor.h"
#include "rt.h"
#include "runner_doiact_sinks.h"
#include "space.h"
#include "star_formation.h"
#include "star_formation_logger.h"
#include "stars.h"
#include "timers.h"
#include "timestep_limiter.h"
#include "tracers.h"
#include "black_holes.h"

#include "kernel_hydro.h"


// Linear search for parent in the cell ->inefficient but oh well
int find_parent_linear(const struct part *parts, int cnt, long long id) {
    for (int i = 0; i < cnt; i++) {
        if (parts[i].id == id)
            return i;
    }
    return -1; // not found
}

/* serialise the splits */
/**** split gas particles within BH hsml ****/
void runner_do_particle_split(struct runner *r, struct cell *c, int timer) {
  struct engine *e = r->e;
  //leaf data
  const int count_leaf = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID)
    error("Running gas splitting task on a foreign node!");
#endif

   /* Anything to do here? */
  if (count_leaf == 0 || !cell_is_active_hydro(c, e)) {
    return;
  }

  lock_lock(&c->top->stars.star_formation_lock);
  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {

	struct cell *restrict cp = c->progeny[k];
	runner_do_particle_split(r, cp, 0);
	
        /* Update the h_max */
        c->hydro.h_max = max(c->hydro.h_max, cp->hydro.h_max);
        c->hydro.h_max_active =
            max(c->hydro.h_max_active, cp->hydro.h_max_active);
      }
  } else {
    
    //hardcode splits for now
    // options are nsplit = 7, 13 (odd because we count the og central particle
    const int n_split = 7;
    
    // Now loop over gas particles at parent cell, this is broked but we get the full
    //set of parent splits in one go, probs best to save the full parent list and then split    
    int j = 0;

    //loop has to be long, idk why this works but it does...
    while (j < c->top->hydro.count){
      struct part *p  = &parts[j];
      struct xpart *xp = &xparts[j];

      //this only really needs to be calculated once
      const float child_mass = p->mass / n_split;  
      double new_h = p->h;
      
      if ((p->split_flag == 2) || (p->split_flag == 0)){
	j ++;
	continue;
      }
      
      // If gas particle is marked for splitting
      if (p->split_flag == 1) {
	//immediately mark as split
	p->split_flag = 2;                  	
	// Update the parent before any memmove
        p->mass = child_mass;
	
	/* force activation
        if (p->gpart){                
          p->gpart->mass = child_mass;
	  if (!gpart_is_active(p->gpart, e)){                                                                                                                           
             message("gparent is not active, ACTIVATE!");
             p->gpart->time_bin = e->min_active_bin;}   
	}                                                                                                                                                                               
	if (!part_is_active(p, e)){                                                                                                                                                                               
          message("parent is not active, ACTIVATE!");
          p->time_bin = e->min_active_bin;}
	*/
	//split into 3D shape
	message("------- splitting particles ------");
	double delta = p->h / 2;
	/*3D cross*/
	
	double offsets[6][3] = {
	  {+delta, 0.0, 0.0}, {-delta, 0.0, 0.0},
	  {0.0, +delta, 0.0}, {0.0, -delta, 0.0},
	  {0.0, 0.0, +delta}, {0.0, 0.0, -delta}
	};
		
	//check particle masses add up
#ifdef SWIFT_DEBUG_CHECKS
	double child_mass_sum = 0.0;
#endif

	// Spawn children
	// problem of looping whilst list shifting
	const long long parent_id = p->id;

	struct part parent = *p;    // copy all fields of parent 
	struct xpart parent_xp = *xp;
	/*
	message("Cell: %p, Parent before splits: id=%lld, mass=%e, pos=(%g,%g,%g), vel=(%g,%g,%g), a_hydro=(%g,%g,%g)",
		(void*)c, parent.id, parent.mass, parent.x[0], parent.x[1], parent.x[2],
		parent.v[0], parent.v[1], parent.v[2], parent.a_hydro[0], parent.a_hydro[1], parent.a_hydro[2]);
	*/
	//struct part **child_list = malloc((n_split-1) * sizeof(struct part *));
	//long long *child_ids = malloc((n_split-1) * sizeof(long long));
	
	for (int n = 0; n < (n_split-1); n++) { 
	  double pos_offsets[3] = {offsets[n][0], offsets[n][1], offsets[n][2]};
	  
	  struct part *child = cell_spawn_new_part_from_part(e, c, &parent, &parent_xp, child_mass, pos_offsets, new_h);
	  
	  //message("Cell: %p, child id: id=%lld has gpart id: %lld", (void*)c, child->id, child->gpart->id_or_neg_offset);
	  
	  //child_list[n] = child;
	  //child_ids[n] = child->id;
	  
	  int parent_index = find_parent_linear(c->top->hydro.parts, c->top->hydro.count, parent_id);
	  if (parent_index == -1) error("parent not found!");
	  
	  if (child == NULL)
	    error("Failed to spawn child particle in gas splitting.");
	  
#ifdef SWIFT_DEBUG_CHECKS
	  child_mass_sum += child->mass;
#endif
	  
	}

	//free(child_list);
	//free(child_ids);
	
	//check mass conservation
#ifdef SWIFT_DEBUG_CHECKS
	double total_mass = child_mass + child_mass_sum;
	double expected_mass = n_split * child_mass;
       
	if (fabs(total_mass - expected_mass) > 1e-8) // tighter tolerance
	  error("Mass conservation broken in splitting: parent+children = %e, expected %e",
		total_mass, n_split * child_mass);
#endif
	/*
	int parent_index = find_parent_linear(c->top->hydro.parts, c->top->hydro.count, parent_id); 
	message("Leaf Cell: %p,Parent ID after full split: %lld, mass=%e",                                                                                                                                
                c,c->top->hydro.parts[parent_index].id, c->top->hydro.parts[parent_index].mass);     
	*/
      }
      j ++;
    }    
  }   
  lock_unlock(&c->top->stars.star_formation_lock);
  /* If we formed any hydro parts, need to force flag for resort and drift*/
  if (c->hydro.count!= count_leaf) {
    //message("setting resort flag");
    //cell_set_hydro_resort_flag(c);
    cell_set_hydro_resort_flag(c->top);
    cell_set_flag(c->top, cell_flag_do_hydro_drift);
    cell_set_flag(c->top, cell_flag_do_grav_drift);
  }
}

extern const int sort_stack_size;

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

#ifdef SWIFT_DEBUG_CHECKS
      if (gp->time_bin == time_bin_not_created)
        error("Found an extra particle in external gravity.");
#endif

      /* Is this part within the time step? */
      if (gpart_is_active(gp, e)) {
        external_gravity_acceleration(time, potential, constants, gp);
      }
    }
  }

  if (timer) TIMER_TOC(timer_dograv_external);
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
  const struct pressure_floor_props *pressure_floor = e->pressure_floor_props;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;
  const double time = e->time;

  TIMER_TIC;

  /* Anything to do here? (i.e. does this cell need updating?) */
  if (!cell_is_active_hydro(c, e)) {
    return;
  }

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

      /* Anything to do here? (i.e. does this particle need updating?) */
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
                          entropy_floor_props, pressure_floor, cooling_func, p,
                          xp, dt_cool, dt_therm, time);
      }
    }
  }

  if (timer) TIMER_TOC(timer_do_cooling);
}

/**
 * @brief Spawns some stars from the sink particles.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_star_formation_sink(struct runner *r, struct cell *c,
                                   int timer) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct phys_const *phys_const = e->physical_constants;
  const struct sink_props *sink_props = e->sink_properties;
  const int count = c->sinks.count;
  struct sink *restrict sinks = c->sinks.parts;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const struct unit_system *restrict us = e->internal_units;
  const int current_stars_count = c->stars.count;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID)
    error("Running star formation task on a foreign node!");
#endif

  /* Anything to do here? */
  if (count == 0 || !cell_is_active_sinks(c, e)) {
    return;
  }

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        /* Load the child cell */
        struct cell *restrict cp = c->progeny[k];

        /* Do the recursion */
        runner_do_star_formation_sink(r, cp, 0);

        /* Update the h_max */
        c->stars.h_max = max(c->stars.h_max, cp->stars.h_max);
        c->stars.h_max_active =
            max(c->stars.h_max_active, cp->stars.h_max_active);
      }
  } else {

    /* Loop over the sink particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct sink *restrict s = &sinks[k];

      /* Only work on active particles */
      if (sink_is_active(s, e)) {

#ifdef WITH_CSDS
        error("TODO");
#endif

        /* Update the sink properties before spwaning stars */
        sink_update_sink_properties_before_star_formation(s, e, sink_props,
                                                          phys_const);

        /* Spawn as many stars as necessary
           - loop counter for the random seed.
           - Start by 1 as 0 is used at init (sink_copy_properties) */
        for (int star_counter = 1; sink_spawn_star(
                 s, e, sink_props, cosmo, with_cosmology, phys_const, us);
             star_counter++) {

          /* Create a new star with a mass s->target_mass */
          struct spart *sp = cell_spawn_new_spart_from_sink(e, c, s);
          if (sp == NULL)
            error("Run out of available star particles or gparts");

          /* Copy the properties to the star particle */
          sink_copy_properties_to_star(s, sp, e, sink_props, cosmo,
                                       with_cosmology, phys_const, us);

          /* Verify that we do not have too many stars in the leaf for
           * the sort task to be able to act. */
          if (c->stars.count > (1LL << sort_stack_size))
            error(
                "Too many stars in the cell tree leaf! The sorting task will "
                "not be able to perform its duties. Possible solutions: (1) "
                "The code need to be run with different star formation "
                "parameters to reduce the number of star particles created. OR "
                "(2) The size of the sorting stack must be increased in "
                "runner_sort.c.");

          /* Update the h_max */
          c->stars.h_max = max(c->stars.h_max, sp->h);
          c->stars.h_max_active = max(c->stars.h_max_active, sp->h);

          /* Update sink properties */
          sink_update_sink_properties_during_star_formation(
              s, sp, e, sink_props, phys_const, star_counter);
        } /* Loop over the stars to spawn */

        /* Update the sink after star formation */
        sink_update_sink_properties_after_star_formation(s, e, sink_props,
                                                         phys_const);
      } /* if sink_is_active */
    } /* Loop over the particles */
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
 * @brief Convert some hydro particles into stars depending on the star
 * formation model.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
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

        /* Update the h_max */
        c->stars.h_max = max(c->stars.h_max, cp->stars.h_max);
        c->stars.h_max_active =
            max(c->stars.h_max_active, cp->stars.h_max_active);

        /* Update the dx_max */
        if (star_formation_need_update_dx_max) {
          c->hydro.dx_max_part =
              max(cp->hydro.dx_max_part, c->hydro.dx_max_part);
          c->hydro.dx_max_sort =
              max(cp->hydro.dx_max_sort, c->hydro.dx_max_sort);
        }
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

            /* Convert the gas particle to a star particle */
            const int n_spart_spawn =
                star_formation_number_spart_to_spawn(p, xp, sf_props);
            const int n_spart_convert =
                star_formation_number_spart_to_convert(p, xp, sf_props);

#ifdef SWIFT_DEBUG_CHECKS
            if (n_spart_convert > 1 || n_spart_convert < 0)
              error("Invalid number of sparts to convert");
#endif

            int n_spart_to_create = n_spart_spawn + n_spart_convert;

            while (n_spart_to_create > 0) {

              struct spart *sp = NULL;
              int part_converted;

              /* Are we using a model that actually generates star particles? */
              if (swift_star_formation_model_creates_stars) {
		
                /* Check if we should create a new particle or transform one */
                if (n_spart_to_create == 1 && n_spart_convert == 1) {
                  /* Convert the gas particle to a star particle */
                  sp = cell_convert_part_to_spart(e, c, p, xp);
                  part_converted = 1;
#ifdef WITH_CSDS
                  /* Write the particle */
                  /* Logs all the fields request by the user */
                  // TODO select only the requested fields
                  csds_log_part(e->csds, p, xp, e, /* log_all */ 1,
                                csds_flag_change_type, swift_type_stars);
#endif
                } else {
                  /* Spawn a new spart (+ gpart) */
                  sp = cell_spawn_new_spart_from_part(e, c, p, xp);
                  part_converted = 0;
                }

              } else {

                /* We are in a model where spart don't exist
                 * --> convert the part to a DM gpart */
                cell_convert_part_to_gpart(e, c, p, xp);
                part_converted = 1;
              }

              /* Did we get a star? (Or did we run out of spare ones?) */
              if (sp != NULL) {

                /* Copy the properties of the gas particle to the star particle
                 */
                star_formation_copy_properties(
                    p, xp, sp, e, sf_props, cosmo, with_cosmology, phys_const,
                    hydro_props, us, cooling, part_converted);

                /* Update the Star formation history */
                star_formation_logger_log_new_spart(sp, &c->stars.sfh);

                /* Update the h_max */
                c->stars.h_max = max(c->stars.h_max, sp->h);
                c->stars.h_max_active = max(c->stars.h_max_active, sp->h);

                /* Update the displacement information */
                if (star_formation_need_update_dx_max) {
                  const float dx2_part = xp->x_diff[0] * xp->x_diff[0] +
                                         xp->x_diff[1] * xp->x_diff[1] +
                                         xp->x_diff[2] * xp->x_diff[2];
                  const float dx2_sort =
                      xp->x_diff_sort[0] * xp->x_diff_sort[0] +
                      xp->x_diff_sort[1] * xp->x_diff_sort[1] +
                      xp->x_diff_sort[2] * xp->x_diff_sort[2];

                  const float dx_part = sqrtf(dx2_part);
                  const float dx_sort = sqrtf(dx2_sort);

                  /* Note: no need to update quantities further up the tree as
                     this task is always called at the top-level */
                  c->hydro.dx_max_part = max(c->hydro.dx_max_part, dx_part);
                  c->hydro.dx_max_sort = max(c->hydro.dx_max_sort, dx_sort);
                }

#ifdef WITH_CSDS
                if (spawn_spart) {
                  /* Set to zero the csds data. */
                  csds_part_data_init(&sp->csds_data);
                } else {
                  /* Copy the properties back to the stellar particle */
                  sp->csds_data = xp->csds_data;
                }

                /* Write the s-particle */
                csds_log_spart(e->csds, sp, e, /* log_all */ 1,
                               csds_flag_create,
                               /* data */ 0);
#endif
              } else if (swift_star_formation_model_creates_stars) {

                /* Do something about the fact no star could be formed.
                   Note that in such cases a tree rebuild to create more free
                   slots has already been triggered by the function
                   cell_convert_part_to_spart() */
                star_formation_no_spart_available(e, p, xp);
              }

              /* We have spawned a particle, decrease the counter of particles
               * to create */
              n_spart_to_create--;

            } /* while n_spart_to_create > 0 */
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
 * @brief Creates sink particles.
 *
 * @param r runner task
 * @param c cell
 */
void runner_do_sink_formation(struct runner *r, struct cell *c) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct sink_props *sink_props = e->sink_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct hydro_props *restrict hydro_props = e->hydro_properties;
  const struct unit_system *restrict us = e->internal_units;
  struct cooling_function_data *restrict cooling = e->cooling_func;
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID)
    error("Running sink formation task on a foreign node!");
#endif

  /* Anything to do here? */
  if (c->hydro.count == 0 || !cell_is_active_hydro(c, e)) {
    return;
  }

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        /* Load the child cell */
        struct cell *restrict cp = c->progeny[k];

        /* Do the recursion */
        runner_do_sink_formation(r, cp);
      }
  } else {

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Only work on active particles */
      if (part_is_active(p, e)) {

        /* Loop over all particles to find the neighbours within r_acc. Then, */
        /* compute all quantities you need to decide to form a sink or not. */
        runner_do_prepare_part_sink_formation(r, c, p, xp);

        /* Is this particle star forming? */
        if (sink_is_forming(p, xp, sink_props, phys_const, cosmo, hydro_props,
                            us, cooling, entropy_floor)) {

          /* Time-step size for this particle */
          double dt_sink;
          if (with_cosmology) {
            const integertime_t ti_step = get_integer_timestep(p->time_bin);
            const integertime_t ti_begin =
                get_integer_time_begin(ti_current - 1, p->time_bin);

            dt_sink =
                cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
          } else {
            dt_sink = get_timestep(p->time_bin, time_base);
          }

          /* Are we forming a sink particle? */
          if (sink_should_convert_to_sink(p, xp, sink_props, e, dt_sink)) {

#ifdef WITH_CSDS
            error("TODO");
#endif

            /* Convert the gas particle to a sink particle */
            struct sink *sink = NULL;

            /* Convert the gas particle to a sink particle */
            sink = cell_convert_part_to_sink(e, c, p, xp);

            /* Did we get a sink? (Or did we run out of spare ones?) */
            if (sink != NULL) {

              /* Copy the properties of the gas particle to the star particle */
              sink_copy_properties(p, xp, sink, e, sink_props, cosmo,
                                   with_cosmology, phys_const, hydro_props, us,
                                   cooling);
            }
          }
        }
      }
    } /* Loop over particles */
  }
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
  const int with_cosmology = e->policy & engine_policy_cosmology;
 
  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_end_hydro_force(r, c->progeny[k], 0);
  } else {

    const struct cosmology *cosmo = e->cosmology;
    //const int count = c->hydro.count;
    struct part *restrict parts = c->hydro.parts;
    struct xpart *restrict xparts = c->hydro.xparts;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < c->hydro.count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      double dt = 0;
      if (part_is_active(p, e)) {

        if (with_cosmology) {
          /* Compute the time step. */
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(e->ti_current - 1, p->time_bin);

          dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
        } else {
          dt = get_timestep(p->time_bin, e->time_base);
        }

	

        /* Finish the force loop */
	//lily change hydro_end_force added xp,e and c
	//hydro_end_force((struct engine *)e, c, xp, p, cosmo, e->time);
	hydro_end_force(p, cosmo);
        mhd_end_force(p, cosmo);
        timestep_limiter_end_force(p);
        chemistry_end_force(p, cosmo, with_cosmology, e->time, dt);
	
        /* Apply the forcing terms (if any) */
        forcing_terms_apply(e->time, e->forcing_terms, e->s,
                            e->physical_constants, p, xp);

#ifdef SWIFT_BOUNDARY_PARTICLES

        /* Get the ID of the part */
        const long long id = p->id;

        /* Cancel hdyro forces of these particles */
        if (id < SWIFT_BOUNDARY_PARTICLES) {

          /* Don't move ! */
          hydro_reset_acceleration(p);
          mhd_reset_acceleration(p);

#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)

          /* Some values need to be reset in the Gizmo case. */
          hydro_prepare_force(p, &c->hydro.xparts[k], cosmo,
                              e->hydro_properties, e->pressure_floor_props,
                              /*dt_alpha=*/0, /*dt_therm=*/0);
          rt_prepare_force(p);
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
  const int with_sinks = (e->policy & engine_policy_sinks);

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
        gravity_end_force(gp, G_newton, potential_normalisation, periodic,
                          with_self_gravity);

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
        else if (gp->type == swift_type_sink)
          id = e->s->sinks[-gp->id_or_neg_offset].id;
        else if (gp->type == swift_type_black_hole)
          id = e->s->bparts[-gp->id_or_neg_offset].id;
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
            !(e->policy & engine_policy_black_holes) &&
            !(e->policy & engine_policy_star_formation) &&
            !(e->policy & engine_policy_sinks)) {

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
            else if (gp->type == swift_type_sink)
              my_id = e->s->sinks[-gp->id_or_neg_offset].id;
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

        /* Deal with sinks' need of potentials */
        if (with_sinks && gp->type == swift_type_gas) {
          const size_t offset = -gp->id_or_neg_offset;
          sink_store_potential_in_part(&s->parts[offset].sink_data, gp);
        }
      }
    }
  }
  if (timer) TIMER_TOC(timer_end_grav_force);
}

/**
 * @brief Write the required particles through the csds.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_csds(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_CSDS
  TIMER_TIC;

  const struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;

  if (c->black_holes.count != 0) {
    error("Black holes are not implemented in the csds.");
  }
  if (c->sinks.count != 0) {
    error("Sink particles are not implemented in the csds.");
  }

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e) &&
      !cell_is_active_stars(c, e))
    return;

  /* Recurse? Avoid spending too much time in useless cells. */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_csds(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be log */
      if (part_is_active(p, e)) {

        if (csds_should_write(&xp->csds_data, e->csds)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          csds_log_part(e->csds, p, xp, e, /* log_all_fields= */ 0,
                        csds_flag_none, /* flag_data= */ 0);
        } else
          /* Update counter */
          xp->csds_data.steps_since_last_output += 1;
      }
    }

    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* Write only the dark matter particles */
      if (gp->type != swift_type_dark_matter &&
          gp->type != swift_type_dark_matter_background)
        continue;

      /* If particle needs to be log */
      if (gpart_is_active(gp, e)) {

        if (csds_should_write(&gp->csds_data, e->csds)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          csds_log_gpart(e->csds, gp, e, /* log_all_fields= */ 0,
                         csds_flag_none, /* flag_data= */ 0);

        } else
          /* Update counter */
          gp->csds_data.steps_since_last_output += 1;
      }
    }

    /* Loop over the sparts in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the part. */
      struct spart *restrict sp = &sparts[k];

      /* If particle needs to be log */
      if (spart_is_active(sp, e)) {

        if (csds_should_write(&sp->csds_data, e->csds)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          csds_log_spart(e->csds, sp, e, /* Log_all_fields= */ 0,
                         csds_flag_none, /* flag_data= */ 0);
        } else
          /* Update counter */
          sp->csds_data.steps_since_last_output += 1;
      }
    }
  }

  if (timer) TIMER_TOC(timer_csds);

#else
  error("CSDS disabled, please enable it during configuration");
#endif
}

/**
 * @brief Recursively search for FOF groups in a single cell.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_fof_search_self(struct runner *r, struct cell *c, int timer) {

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
void runner_do_fof_search_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer) {

#ifdef WITH_FOF

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != cj->nodeID) error("Searching foreign cells!");
#endif

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

/**
 * @brief Recursively search for FOF groups in a single cell.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_fof_attach_self(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_FOF

  TIMER_TIC;

  const struct engine *e = r->e;
  struct space *s = e->s;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const struct gpart *const gparts = s->gparts;
  const double attach_r2 = e->fof_properties->l_x2;

  rec_fof_attach_self(e->fof_properties, dim, attach_r2, periodic, gparts,
                      s->nr_gparts, c);

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
void runner_do_fof_attach_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer) {

#ifdef WITH_FOF

  TIMER_TIC;

  const struct engine *e = r->e;
  struct space *s = e->s;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const struct gpart *const gparts = s->gparts;
  const double attach_r2 = e->fof_properties->l_x2;

  rec_fof_attach_pair(e->fof_properties, dim, attach_r2, periodic, gparts,
                      s->nr_gparts, ci, cj, e->nodeID == ci->nodeID,
                      e->nodeID == cj->nodeID);

  if (timer) TIMER_TOC(timer_fof_pair);
#else
  error("SWIFT was not compiled with FOF enabled!");
#endif
}

/**
 * @brief Finish up the transport step and do the thermochemistry
 *        for radiative transfer
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_rt_tchem(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const int count = c->hydro.count;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct rt_props *rt_props = e->rt_props;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct cosmology *cosmo = e->cosmology;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;

  /* Anything to do here? */
  if (count == 0) return;
  if (!cell_is_rt_active(c, e)) return;

  TIMER_TIC;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_rt_tchem(r, c->progeny[k], 0);
  } else {

    struct part *restrict parts = c->hydro.parts;
    struct xpart *restrict xparts = c->hydro.xparts;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Skip inhibited parts */
      if (part_is_inhibited(p, e)) continue;

      /* Skip inactive parts */
      if (!part_is_rt_active(p, e)) continue;

      /* Finish the force loop */
      const integertime_t ti_current_subcycle = e->ti_current_subcycle;
      const integertime_t ti_step =
          get_integer_timestep(p->rt_time_data.time_bin);
      const integertime_t ti_begin = get_integer_time_begin(
          ti_current_subcycle + 1, p->rt_time_data.time_bin);
      const integertime_t ti_end = ti_begin + ti_step;

      const double dt =
          rt_part_dt(ti_begin, ti_end, e->time_base, with_cosmology, cosmo);
#ifdef SWIFT_DEBUG_CHECKS
      if (ti_begin != ti_current_subcycle)
        error(
            "Particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
            "ti_step=%lld time_bin=%d wakeup=%d ti_current=%lld",
            ti_end, ti_begin, ti_step, p->time_bin, p->limiter_data.wakeup,
            ti_current_subcycle);
      if (dt < 0.)
        error("Got part with negative time-step: %lld, %.6g", p->id, dt);
#endif

      rt_finalise_transport(p, rt_props, dt, cosmo);

      /* And finally do thermochemistry */
      rt_tchem(p, xp, rt_props, cosmo, hydro_props, phys_const, us, dt);
    }
  }

  if (timer) TIMER_TOC(timer_do_rt_tchem);
}
