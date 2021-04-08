/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Local headers. */
#include "black_holes.h"
#include "cosmology.h"
#include "gravity.h"
#include "hydro.h"
#include "lightcone.h"
#include "lightcone_io.h"
#include "lightcone_replications.h"
#include "part.h"
#include "stars.h"
#include "timeline.h"

#ifndef SWIFT_LIGHTCONE_CROSSING_H
#define SWIFT_LIGHTCONE_CROSSING_H

/**
 * @brief Check if a particle crosses the lightcone during a drift.
 *
 * Here we don't assume anything about the particle type except
 * that it has a corresponding #gpart. x and v_full must be set
 * by the calling function because gp->x and gp->v_full may not
 * be the right values to use for some particle types.
 *
 * The particle type is checked if we decide to output the
 * particle, at which point we call a type specific function.
 *
 */
__attribute__((always_inline)) INLINE static void lightcone_check_particle_crosses(
     struct lightcone_props *props, struct replication_list *replication_list,
     const struct cosmology *c, const struct gpart *gp, const double *x,
     const float *v_full, const double dt_drift, const integertime_t ti_old,
     const integertime_t ti_current, const double cell_loc[3],
     void *extra1, void *extra2) {

  /* Are we making a lightcone and in the lightcone redshift range? */
  if(!props->enabled || replication_list->nrep==0)return;

  /* Are we using this particle type in the lightcone? */
  if(!props->use_type[gp->type])return;
  
  /* Unpack some variables we need */
  const double boxsize = props->boxsize;
  const double *observer_position = props->observer_position;
  const int nreps = replication_list->nrep;
  const struct replication *rep = replication_list->replication;

  /* Consistency check - are our limits on the drift endpoints good? */
  if(ti_old < props->ti_old || ti_current > props->ti_current)
    error("Particle drift is outside the range used to make replication list!");

  /* Determine expansion factor at start and end of the drift */
  const double a_start = c->a_begin * exp(ti_old     * c->time_base);
  const double a_end   = c->a_begin * exp(ti_current * c->time_base);

  /* Does this drift overlap the lightcone redshift range? If not, nothing to do. */
  if((a_start > props->a_at_z_min) || (a_end < props->a_at_z_max))return;

  /* Find comoving distance to these expansion factors */
  const double comoving_dist_start   = cosmology_get_comoving_distance(c, a_start);
  const double comoving_dist_2_start = comoving_dist_start*comoving_dist_start;
  const double comoving_dist_end     = cosmology_get_comoving_distance(c, a_end);
  const double comoving_dist_2_end   = comoving_dist_end*comoving_dist_end;

  /* Thickness of the 'shell' between the lightcone surfaces at start and end of drift.
     We use this as a limit on how far a particle can drift (i.e. assume v < c).*/
  const double boundary = comoving_dist_2_start - comoving_dist_2_end;

  /* Wrap particle starting coordinates to nearest it's parent cell */
  const double x_wrapped[3] = {box_wrap(x[0], cell_loc[0]-0.5*boxsize, cell_loc[0]+0.5*boxsize),
                               box_wrap(x[1], cell_loc[1]-0.5*boxsize, cell_loc[1]+0.5*boxsize),
                               box_wrap(x[2], cell_loc[2]-0.5*boxsize, cell_loc[2]+0.5*boxsize)};
  
  /* Get wrapped position relative to observer */
  const double x_wrapped_rel[3] = {x_wrapped[0] - observer_position[0],
                                   x_wrapped[1] - observer_position[1],
                                   x_wrapped[2] - observer_position[2]};

  /* Loop over periodic copies of the volume:
     
     Here we're looking for cases where a periodic copy of the particle
     is closer to the observer than the lightcone surface at the start
     of the drift, and further away than the lightcone surface at the
     end of the drift. I.e. the surface of the lightcone has swept over
     the particle as it contracts towards the observer.
   */
  for(int i=0; i<nreps; i+=1) {

    /* If all particles in this periodic replica are beyond the lightcone surface
       at the earlier time, then they already crossed the lightcone. Since the
       replications are in ascending order of rmin we don't need to check any
       more. */
    if(rep[i].rmin2 > comoving_dist_2_start)break;

    /* If all particles in this periodic replica start their drifts inside the
       lightcone surface, and are sufficiently far inside that their velocity
       can't cause them to cross the lightcone, then we don't need to consider
       this replication */
    if(rep[i].rmax2 + boundary < comoving_dist_2_end)continue;

    /* Get the coordinates of this periodic copy of the gpart relative to the observer */
    const double x_start[3] = {
      x_wrapped_rel[0] + rep[i].coord[0],
      x_wrapped_rel[1] + rep[i].coord[1],
      x_wrapped_rel[2] + rep[i].coord[2],
    };

    /* Get distance squared from the observer at start of drift */
    const double r2_start =
      x_start[0]*x_start[0]+
      x_start[1]*x_start[1]+
      x_start[2]*x_start[2];

    /* If particle is initially beyond the lightcone surface, it can't cross */
    if(r2_start > comoving_dist_2_start)continue;

    /* Get position of this periodic copy at the end of the drift */
    const double x_end[3] = {
      x_start[0] + dt_drift * v_full[0],
      x_start[1] + dt_drift * v_full[1],
      x_start[2] + dt_drift * v_full[2],
    };

    /* Get distance squared from the observer at end of drift */
    const double r2_end =
      x_end[0]*x_end[0]+
      x_end[1]*x_end[1]+
      x_end[2]*x_end[2];
    
    /* If particle is still within the lightcone surface at the end of the drift,
       it didn't cross*/
    if(r2_end < comoving_dist_2_end)continue;

    /* This periodic copy of the gpart crossed the lightcone during this drift.
       Now need to estimate when it crossed within the timestep.

       If r is the distance from the observer to this periodic copy of the particle,
       and it crosses after a fraction f of the drift:

       r_cross = r_start + (r_end - r_start) * f

       and if R is the comoving distance to the lightcone surface

       R_cross = R_start + (R_end - R_start) * f

       The particle crosses the lightcone when r_cross = R_cross, so

       r_start + (r_end - r_start) * f = R_start + (R_end - R_start) * f

       Solving for f:

       f = (r_start - R_start) / (R_end - R_start - r_end + r_start)

    */
    const double f = (sqrt(r2_start) - comoving_dist_start) /
      (comoving_dist_end-comoving_dist_start-sqrt(r2_end)+sqrt(r2_start));

    /* f should always be in the range 0-1 */
    const double eps = 1.0e-5;
    if((f < 0.0-eps) || (f > 1.0+eps)) error("Particle interpolated outside time step!");

    /* Compute position at crossing */
    const double x_cross[3] = {
      x_start[0] + dt_drift * f * v_full[0],
      x_start[1] + dt_drift * f * v_full[1],
      x_start[2] + dt_drift * f * v_full[2],
    };

    /* Check if particle is in our distance limits at crossing */
    const double r2_cross = (x_cross[0]*x_cross[0] +
                             x_cross[1]*x_cross[1] + 
                             x_cross[2]*x_cross[2]);
    if(r2_cross < props->r2_min || r2_cross > props->r2_max)continue;

    /* Compute expansion factor at crossing */
    const double a_cross = cosmology_scale_factor_at_comoving_distance(c, sqrt(r2_cross));

    /* Need to write out this particle */
    switch (gp->type) {

    case swift_type_gas: {

      const struct part *p = (struct part *) extra1;
      const struct xpart *xp = (struct xpart *) extra2;
      struct lightcone_gas_data data;
      lightcone_store_gas(gp, p, xp, a_cross, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_gas, &data);
 
    } break;

    case swift_type_stars: {

      const struct spart *sp = (struct spart *) extra1;
      struct lightcone_stars_data data;
      lightcone_store_stars(gp, sp, a_cross, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_stars, &data);

    } break;

    case swift_type_black_hole: {

      const struct bpart *bp = (struct bpart *) extra1;
      struct lightcone_black_hole_data data;
      lightcone_store_black_hole(gp, bp, a_cross, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_black_hole, &data);

    } break;

    case swift_type_dark_matter: {

      struct lightcone_dark_matter_data data;
      lightcone_store_dark_matter(gp, a_cross, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_dark_matter, &data);
 
    } break;

    case swift_type_dark_matter_background: {

      /* Assumed to have same properties as DM particles */
      struct lightcone_dark_matter_data data;
      lightcone_store_dark_matter(gp, a_cross, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_dark_matter_background, &data);
            
    } break;

    case swift_type_neutrino: {

      struct lightcone_neutrino_data data;
      lightcone_store_neutrino(gp, a_cross, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_neutrino, &data);

    } break;

    default:
      error("Particle type not supported in lightcones");
    }

  } /* Next periodic replication*/
}

#endif /* SWIFT_LIGHTCONE_CROSSING_H */
