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
#include <string.h>
#include <hdf5.h>

/* This object's header. */
#include "lightcone.h"

/* Local headers */
#include "common_io.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "lightcone_io.h"
#include "lock.h"
#include "parser.h"
#include "particle_buffer.h"
#include "periodic.h"
#include "periodic_replications.h"
#include "restart.h"
#include "space.h"
#include "timeline.h"

/* Whether to dump the replication list */
//#define DUMP_REPLICATIONS
#ifdef DUMP_REPLICATIONS
static int output_nr = 0;
#endif

/* MPI rank for diagnostic messages */
extern int engine_rank;

/**
 * @brief Dump lightcone_props struct to the output stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to write to.
 */
void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  /* Don't dump the replication list - will regenerate it as needed */
  struct lightcone_props tmp = *props;
  tmp.replication_list.nrep = 0;
  tmp.replication_list.replication = NULL;
  tmp.have_replication_list = 0;

  /* Don't write out particle buffers - must flush before dumping restart. */
  memset(tmp.buffer, 0, sizeof(struct particle_buffer)*swift_type_count);

  restart_write_blocks((void *) &tmp, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");
}


/**
 * @brief Restore lightcone_props struct from the input stream.
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
 * If restarting, this is called after lightcone_struct_restore().
 *
 * @param props the #lightcone_props structure to fill.
 * @param params the parameter file parser.
 */
void lightcone_init(struct lightcone_props *props,
                    const struct space *s,
                    const struct cosmology *cosmo,
                    struct swift_params *params,
                    const int restart) {
  
  /* Whether we generate lightcone output */
  props->enabled = 1;

  /* Base name for output files */
  parser_get_param_string(params, "Lightcone:basename", props->basename);

  /* Redshift range for the lightcone */
  props->z_min = parser_get_param_double(params, "Lightcone:z_min");
  props->z_max = parser_get_param_double(params, "Lightcone:z_max");

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, "Lightcone:observer_position", 3,
                                props->observer_position);

  /* Write particles to disk if this many or more are in the buffer */
  props->max_particles_buffered = parser_get_opt_param_int(params, "Lightcone:max_particles_buffered", 100000);

  /* Whether we're doing a pencil beam */
  props->pencil_beam = parser_get_opt_param_int(params, "Lightcone:pencil_beam", 0);
  
  /* Direction of the pencil beam */
  for(int i=0; i<3; i+=1)
    props->view_vector[i] = 0.0;
  parser_get_opt_param_double_array(params, "Lightcone:view_vector", 3, props->view_vector);
  
  /* Radius of the pencil beam (radians) */
  props->view_radius = parser_get_opt_param_double(params, "Lightcone:view_radius", 0.0);

  if(props->pencil_beam) {

    /* Normalize the view vector */
    double dr2 = 0.0;
    for(int i=0; i<3; i+=1)
      dr2 += props->view_vector[i] * props->view_vector[i];
    if(dr2==0.0)error("Must specify non-zero Lightcone:view_vector if Lightcone:pencil_beam != 0");
    for(int i=0; i<3; i+=1)
      props->view_vector[i] /= sqrt(dr2);

    /* Sanity check on radius */
    if((props->view_radius <= 0.0) || (props->view_radius >= 0.5*M_PI))
      error("Must have 0 < Lightcone:view_radius < pi/2 radians");

  }

  /* Get the size of the simulation box */
  props->boxsize = s->dim[0];
  if(s->dim[1] != s->dim[0] || s->dim[2] != s->dim[0])
    error("Lightcones require a cubic simulation box.");

  /* Store expansion factors at z_min, z_max */
  props->a_at_z_min = 1.0/(1.0+props->z_min);
  props->a_at_z_max = 1.0/(1.0+props->z_max);

  /* Estimate number of particles which will be output.

     Assumptions:
     - flat cosmology (haven't implemented comoving volume calculation for non-flat)
     - uniform box
  */
  const long long nr_gparts = s->nr_gparts;
  long long total_nr_gparts;
#ifdef WITH_MPI
  MPI_Reduce(&nr_gparts, &total_nr_gparts, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  total_nr_gparts = nr_gparts;
#endif
  if(engine_rank==0) {
    const double lightcone_rmax = cosmology_get_comoving_distance(cosmo, props->a_at_z_max);
    const double lightcone_rmin = cosmology_get_comoving_distance(cosmo, props->a_at_z_min);
    const double volume = 4./3.*M_PI*(pow(lightcone_rmax, 3.)-pow(lightcone_rmin, 3.));
    const long long est_nr_output = total_nr_gparts / pow(props->boxsize, 3.0) * volume;
    message("gparts in lightcone (if uniform box+flat cosmology): %lld", est_nr_output);
  }

  /* Initially have no replication list */
  props->have_replication_list = 0;
  props->ti_old = 0;
  props->ti_current = 0;

  /* If we're not restarting, initialize various counters */
  if(!restart) {
    for(int i=0; i<swift_type_count; i+=1) {
      props->tot_num_particles_written[i] = 0;
      props->num_particles_written_to_file[i] = 0;
    }
    props->current_file = -1;
  }

  /* Always start a new file initially, whether starting or restarting */
  props->start_new_file = 1;

  /* Initialize particle output buffers */
  const size_t elements_per_block = 100000;

  particle_buffer_init(&props->buffer[swift_type_gas],
                       sizeof(struct lightcone_gas_data),
                       elements_per_block);

  particle_buffer_init(&props->buffer[swift_type_dark_matter],
                       sizeof(struct lightcone_dark_matter_data),
                       elements_per_block);

  particle_buffer_init(&props->buffer[swift_type_dark_matter_background],
                       sizeof(struct lightcone_dark_matter_background_data),
                       elements_per_block);

  particle_buffer_init(&props->buffer[swift_type_stars],
                       sizeof(struct lightcone_stars_data),
                       elements_per_block);

  particle_buffer_init(&props->buffer[swift_type_neutrino],
                       sizeof(struct lightcone_neutrino_data),
                       elements_per_block);

}


/**
 * @brief Flush any buffers which exceed the specified size.
 */
void lightcone_flush_buffers(struct lightcone_props *props,
                             int flush_all, int end_file) {

  /* Count how many types have data to write out */
  int types_to_flush = 0;

  /* Will flush any buffers with more particles than this */
  int max_to_buffer = props->max_particles_buffered;
  if(flush_all)max_to_buffer = 0;

  /* Loop over particle types we know how to write out */
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    switch(ptype) {
    case swift_type_gas:
    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
    case swift_type_stars:
    case swift_type_black_hole:
    case swift_type_neutrino:
      if(props->buffer[ptype].total_num_elements >= max_to_buffer &&
         props->buffer[ptype].total_num_elements > 0) {
        types_to_flush += 1;
      }
      break;
    default:
      /* Don't support this type, so do nothing */
      break;
    }
  }
  
  /* Check if there's anything to do */
  if(types_to_flush>0) {
    
    /* We have data to flush, so open or create the output file */
    hid_t file_id;
    char fname[FILENAME_BUFFER_SIZE];
    if(props->start_new_file) {

      /* Get the name of the next file */
      props->current_file += 1;
      sprintf(fname, "%s_%04d.%d.hdf5", props->basename, props->current_file, engine_rank);

      /* Create the file */
      file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if(file_id < 0)error("Unable to create new lightcone file: %s", fname);

      /* We have now written no particles to the current file */
      for(int ptype=0; ptype<swift_type_count; ptype+=1)
        props->num_particles_written_to_file[ptype] = 0;    

      /* We no longer need to create a new file */
      props->start_new_file = 0;

    } else {

      /* Re-open an existing file */
      sprintf(fname, "%s_%04d.%d.hdf5", props->basename, props->current_file, engine_rank);
      file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
      if(file_id < 0)error("Unable to open current lightcone file: %s", fname);

    }

    /* Loop over particle types */
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {

      const size_t num_to_write = props->buffer[ptype].total_num_elements;
      if(num_to_write >= max_to_buffer && num_to_write > 0) {
        switch(ptype) {
        case swift_type_gas:
          error("Lightcone output not implemented yet for gas");
          break;
        case swift_type_dark_matter:
          lightcone_write_dark_matter(props, file_id, ptype);
          particle_buffer_empty(&props->buffer[ptype]);
          props->num_particles_written_to_file[ptype] += num_to_write;
          break;
        case swift_type_dark_matter_background:
          error("Lightcone output not implemented yet for background DM");
          break;
        case swift_type_stars:
          error("Lightcone output not implemented yet for stars");
          break;
        case swift_type_black_hole:
          error("Lightcone output not implemented yet for black holes");
          break;
        case swift_type_neutrino:
          error("Lightcone output not implemented yet for neutrinos");
          break;
        default:
          /* Don't support this type, so do nothing */
          break;
        }
      }
    } /* Next particle type */

    /* We're done updating the output file */
    H5Fclose(file_id);
  }

  /* If we need to start a new file next time, record this */
  if(end_file)props->start_new_file = 1;

}


/**
 * @brief Determine periodic copies of the simulation box which could
 * contribute to the lightcone.
 *
 *                     \
 *           \          \
 *            |         |
 * Obs      A |    B    | C
 *            |         |
 *           /          /
 *          R1         /
 *                    R0
 *
 * Consider a single particle being drifted. Here R0 is the comoving
 * distance to the time the particle is drifted FROM. R1 is the comoving
 * distance to the time the particle is drifted TO on this step.
 *
 * Particles which are beyond the lightcone surface at the start of
 * their drift (C) cannot cross the lightcone on this step if v < c.
 * Particles between the lightcone surfaces at the start and end of
 * their drift (B) may cross the lightcone (and certainly will if they
 * have zero velocity).
 *
 * Particles just within the lightcone surface at the start of their
 * drift (A) may be able to cross the lightcone due to their velocity so
 * we need to allow a boundary layer on the inside edge of the shell.
 * If we assume v < c, then we can use a layer of thickness R0-R1.
 *
 * Here we compute the earliest and latest times particles may be drifted
 * between, find the corresponding comoving distances R0 and R1, reduce
 * the inner distance by R0-R1, and find all periodic copies of the 
 * simulation box which overlap this spherical shell.
 *
 * Later we use this list to know which periodic copies to check when
 * particles are drifted.
 *
 * @param props The #lightcone_props structure
 * @param cosmo The #cosmology structure
 * @param s The #space structure
 * @param time_min Start of the time step
 * @param time_max End of the time step
 */
void lightcone_init_replication_list(struct lightcone_props *props,
                                     const struct cosmology *cosmo,
                                     const integertime_t ti_old,
                                     const integertime_t ti_current,
                                     const double dt_max) {

  /* Deallocate the old list, if there is one */
  if(props->have_replication_list)replication_list_clean(&props->replication_list);

  /* Get the size of the simulation box */
  const double boxsize = props->boxsize;

  /* Get a lower limit on earliest time particle may be drifted from */
  float dt = cosmo->time_end - cosmo->time_begin;
  while (dt > dt_max) dt /= 2.f;
  timebin_t bin = get_time_bin(dt*cosmo->time_base_inv);
  integertime_t ti_lim = get_integer_time_begin(ti_old, bin);

  /* Get expansion factor at earliest and latest times particles might be drifted between */
  double a_current = cosmo->a_begin * exp(ti_current * cosmo->time_base);
  double a_old = cosmo->a_begin * exp(ti_lim * cosmo->time_base);
  if(a_old < cosmo->a_begin)a_old = cosmo->a_begin;

  /* Convert redshift range to a distance range */
  double lightcone_rmin = cosmology_get_comoving_distance(cosmo, a_current);
  double lightcone_rmax = cosmology_get_comoving_distance(cosmo, a_old);
  if(lightcone_rmin > lightcone_rmax)
    error("Lightcone has rmin > rmax");

  /* Allow inner boundary layer, assuming all particles have v < c.
     This is to account for particles moving during the time step. */
  double boundary = lightcone_rmax-lightcone_rmin;
  lightcone_rmin -= boundary;
  if(lightcone_rmin < 0)lightcone_rmin = 0;

  /* Determine periodic copies we need to search */
  replication_list_init(&props->replication_list, boxsize,
                        props->observer_position,
                        lightcone_rmin, lightcone_rmax,
                        props->pencil_beam, props->view_vector,
                        props->view_radius, boundary);

  /* Record that we made the list */
  props->have_replication_list = 1;

  /* Store times we used to make the list, for consistency check later */
  props->ti_old = ti_lim;
  props->ti_current = ti_current;

  /* Report the size of the list */
#ifdef DUMP_REPLICATIONS
  if(engine_rank==0) {
    message("no. of replications to check: %d", props->replication_list.nrep);
    message("shell to search inner radius=%e, outer radius=%e", lightcone_rmin,
            lightcone_rmax);
  }
#endif

  /* Write out the list, if required */
#ifdef DUMP_REPLICATIONS
  if(engine_rank==0) {
    char fname[500];
    sprintf(fname, "replication_list.%d.txt", output_nr);
    FILE *fd_rep = fopen(fname, "w");
    fprintf(fd_rep, "# Observer x, y, z\n");
    fprintf(fd_rep, "%e, %e, %e\n", props->observer_position[0],
            props->observer_position[1], props->observer_position[2]); 
    fprintf(fd_rep, "# Box size, inner radius, outer radius\n");
    fprintf(fd_rep, "%e, %e, %e\n", boxsize, lightcone_rmin-boundary, lightcone_rmax);
    fprintf(fd_rep, "# x, y, z, rmin2, rmax2\n");
    replication_list_write(&props->replication_list, fd_rep);
    fclose(fd_rep);
    output_nr += 1;
  }
#endif
}


/**
 * @brief Check if a gpart crosses the lightcone during a drift.
 *
 * Here we don't assume anything about the particle type except
 * that it has a corresponding #gpart. x and v_full must be set
 * by the calling function because gp->x and gp->v_full may not
 * be the right values to use for some particle types.
 *
 * The particle type is checked if we decide to output the
 * particle, at which point we call a type specific function.
 *
 * @param e The #engine structure.
 * @param gp The #gpart to check.
 */
void lightcone_check_particle_crosses(struct lightcone_props *props,
                                      const struct cosmology *c, const struct gpart *gp,
                                      const double *x, const float *v_full,
                                      const double dt_drift, const integertime_t ti_old,
                                      const integertime_t ti_current) {

  /* Are we making a lightcone? */
  if(!props->enabled)return;

  /* Unpack some variables we need */
  const double boxsize = props->boxsize;
  const double *observer_position = props->observer_position;
  const int nreps = props->replication_list.nrep;
  const struct replication *rep = props->replication_list.replication;

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

  /* Wrap particle starting coordinates into the box */
  const double x_wrapped[3] = {box_wrap(x[0], 0.0, boxsize),
                               box_wrap(x[1], 0.0, boxsize),
                               box_wrap(x[2], 0.0, boxsize)};
  
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

    /* Compute expansion factor at crossing (approximate, used for debugging) */    
    //const double ti_cross = ti_old + (ti_current-ti_old)*f;
    //const double a_cross = c->a_begin * exp(ti_cross * c->time_base);

    /* Compute position at crossing */
    const double x_cross[3] = {
      x_start[0] + dt_drift * f * v_full[0],
      x_start[1] + dt_drift * f * v_full[1],
      x_start[2] + dt_drift * f * v_full[2],
    };

    /* Check if the particle is in the field of view */
    if(props->pencil_beam) {
      
      /* Get distance to the particle at time of crossing */
      double r_cross = 0;
      for(i=0; i<3; i+=1)      
        r_cross += x_cross[i]*x_cross[i];
      r_cross = sqrt(r_cross);

      /* Find dot product of view vector and particle position at crossing */
      double dp = 0.0;
      for(i=0; i<3; i+=1)
        dp += x_cross[i]*props->view_vector[i];
      
      /* Find angle between line of sight and particle position
         (assuming view_vector is normalized) */
      double theta = acos(dp/r_cross);

      /* If particle is outside view radius at crossing, don't output it */
      if(theta > props->view_radius)continue;

    }

    /* Need to write out this particle */
    switch (gp->type) {

    case swift_type_gas: {
      /*
      const struct part *p = &parts[-gparts[i].id_or_neg_offset];
      const struct xpart *xp = &xparts[-gparts[i].id_or_neg_offset];
      */
      error("Gas particle lightcone output not implemented yet");
    } break;

    case swift_type_stars: {
      /*
      const struct spart *sp = &sparts[-gparts[i].id_or_neg_offset];
      */
      error("Star particle lightcone output not implemented yet");
    } break;

    case swift_type_black_hole: {
      /*
      const struct bpart *bp = &bparts[-gparts[i].id_or_neg_offset];
      */
      error("BH particle lightcone output not implemented yet");
    } break;

    case swift_type_dark_matter: {

      struct lightcone_dark_matter_data data;
      lightcone_store_dark_matter(gp, x_cross, &data);
      particle_buffer_append(props->buffer+swift_type_dark_matter, &data);
 
    } break;

    case swift_type_dark_matter_background: {
            
    } break;

    case swift_type_neutrino: {
      error("Neutrino particle lightcone output not implemented yet");
    } break;

    default:
      error("Particle type not implemented in lightcones.");
    }

  } /* Next periodic replication*/
}
