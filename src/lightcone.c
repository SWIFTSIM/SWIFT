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

/* Some standard headers. */
#include <stdio.h>
#include <math.h>

/* This object's header. */
#include "lightcone.h"

/* Local headers */
#include "engine.h"
#include "error.h"
#include "cosmology.h"
#include "parser.h"
#include "periodic.h"
#include "periodic_replications.h"
#include "restart.h"
#include "space.h"
#include "timeline.h"

/**
 * @brief Dump lightcone_props struct to the output stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to write to.
 */
void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");
}


/**
 * @brief Restore lightcone_props struct from the output stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to read from.
 */
void lightcone_struct_restore(struct lightcone_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                      NULL, "lightcone_props");
}


/**
 * @brief Initialise the properties of the lightcone code.
 *
 * @param props the #lightcone_props structure to fill.
 * @param params the parameter file parser.
 */
void lightcone_init(struct lightcone_props *props,
                    const struct space *s,
                    struct swift_params *params) {
  
  /* Whether we generate lightcone output */
  props->enabled = 1;

  /* Redshift range for the lightcone */
  props->z_min = parser_get_param_double(params, "Lightcone:z_min");
  props->z_max = parser_get_param_double(params, "Lightcone:z_max");

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, "Lightcone:observer_position", 3,
                                props->observer_position);

  /* Get the size of the simulation box */
  props->boxsize = s->dim[0];
  if(s->dim[1] != s->dim[0] || s->dim[2] != s->dim[0])
    error("Lightcones require a cubic simulation box.");
}


/**
 * @brief Determine periodic copies of the simulation box which could
 * contribute to the lightcone.
 *
 * @param props the #lightcone_props structure
 * @param props the #cosmology structure
 * @param props the #space structure
 */
void lightcone_init_replication_list(struct lightcone_props *props,
                                     struct cosmology *cosmo,
                                     struct space *s) {
  /* 
     For now, we'll check all periodic replications between z_min and z_max.
     
     TODO:

     - on each timestep, get limits on times particles can be drifted between
       and regenerate the list of replications.
     - set lightcone_boundary in a more reasonable way
     - sort boxes by distance so we can exit sooner when checking for crossings?
  */

  /* Get the size of the simulation box */
  const double boxsize = props->boxsize;

  /* Convert redshift range to a distance range */
  double lightcone_rmin = cosmology_get_comoving_distance(cosmo, 1.0/(1.0+props->z_max));
  double lightcone_rmax = cosmology_get_comoving_distance(cosmo, 1.0/(1.0+props->z_min));;
  double lightcone_boundary = boxsize; /* Needs to be an upper limit on drift distance */
  if(lightcone_rmin > lightcone_rmax)
    error("Lightcone has rmin > rmax - check z_min and z_max parameters?");

  /* Determine periodic copies we need to search */
  replication_list_init(&props->replication_list, boxsize,
                        props->observer_position,
                        lightcone_rmin, lightcone_rmax,
                        lightcone_boundary);
}


/**
 * @brief Check if a gpart crosses the lightcone during a drift.
 *
 * @param e The #engine structure.
 * @param gp The #gpart to check.
 */
void lightcone_check_gpart_crosses(const struct engine *e, const struct gpart *gp,
                                   const double dt_drift, const integertime_t ti_old) {

  /* Unpack some variables we need */
  const struct lightcone_props *props = e->lightcone_properties;
  const double boxsize = props->boxsize;
  const double *observer_position = props->observer_position;
  const int nreps = props->replication_list.nrep;
  const struct replication *rep = props->replication_list.replication;
  const struct cosmology *c = e->cosmology;

  /* Determine expansion factor at start and end of the drift */
  const double a_start = c->a_begin * exp(ti_old * c->time_base);
  const double a_end   = c->a_begin * exp(ti_old * c->time_base + dt_drift);

  /* Find comoving distance to these expansion factors */
  const double comoving_dist_2_start = pow(cosmology_get_comoving_distance(c, a_start), 2.0);
  const double comoving_dist_2_end   = pow(cosmology_get_comoving_distance(c, a_end), 2.0);

  /* Wrap particle starting coordinates into the box. The particle can
     still wander out of the box as it is drifted, but we account for 
     that by including a boundary layer when we decide which periodic
     copies to search. */
  const double x_wrapped[3] = {box_wrap(gp->x[0], 0.0, boxsize),
                               box_wrap(gp->x[1], 0.0, boxsize),
                               box_wrap(gp->x[2], 0.0, boxsize)};
  
  /* Loop over periodic copies of the volume */
  for(int i=0; i<nreps; i+=1) {

    /* Get the coordinates of this periodic copy of the gpart relative to the observer */
    const double x_start[3] = {
      x_wrapped[0] + rep[i].coord[0]*boxsize - observer_position[0],
      x_wrapped[1] + rep[i].coord[1]*boxsize - observer_position[1],
      x_wrapped[2] + rep[i].coord[2]*boxsize - observer_position[2],
    };

    /* Get distance squared from the observer at start of drift */
    const double r2_start =
      x_start[0]*x_start[0]+
      x_start[1]*x_start[1]+
      x_start[2]*x_start[2];

    /* Get position of this periodic copy at the end of the drift */
    const double x_end[3] = {
      x_start[0] + dt_drift * gp->v_full[0],
      x_start[1] + dt_drift * gp->v_full[1],
      x_start[2] + dt_drift * gp->v_full[2],
    };

    /* Get distance squared from the observer at end of drift */
    const double r2_end =
      x_end[0]*x_end[0]+
      x_end[1]*x_end[1]+
      x_end[2]*x_end[2];
    
    if(r2_start < comoving_dist_2_start && r2_end > comoving_dist_2_end) {

      /* This periodic copy of the gpart crossed the lightcone during this drift */
      
      
      
    }

  }

}
