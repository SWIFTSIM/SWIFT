/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2022 Will Roper (w.roper@sussex.ac.uk)
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

/* Includes. */
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* This object's header. */
#include "engine.h"

/* Local includes. */
#include "cell.h"
#include "fof.h"
#include "fof_catalogue_io.h"
#include "halo_finder/halo.h"

/**
 * @brief Allocate the memory and initialise the arrays for a FOF calculation.
 *
 * @param s The #space to act on.
 * @param props The properties of the FOF structure.
 * @param props Are we running the halo finder?
 */
void halo_finder_allocate(const struct space *s,
                          struct fof_props *props) {

  const int verbose = s->e->verbose;
  const ticks total_tic = getticks();

  /* Safety check */
  if (!(s->e->policy & engine_policy_cosmology))
    error(
          "Attempting to run FoF on a simulation using cosmological "
          "information but cosmology was not initialised");
  
  /* Calculate the mean inter-particle separation as if we were in
     a scenario where the entire box was filled with high-resolution
     particles */
  const double Omega_cdm = s->e->cosmology->Omega_cdm;
  const double Omega_b = s->e->cosmology->Omega_b;
  const double Omega_m = Omega_cdm + Omega_b;
  const double critical_density_0 = s->e->cosmology->critical_density_0;
  double mean_matter_density;
  if (s->with_hydro)
    mean_matter_density = Omega_cdm * critical_density_0;
  else
    mean_matter_density = Omega_m * critical_density_0;


  /* Let's report some interesting fats */
  if (props->l_x_absolute != -1.) {

    if (s->e->nodeID == 0) 
      message("Host linking length is set to %e [internal units].",
              props->l_x_absolute);

    /* Define the subhalo linking length based on overdensity ratio. */
    if (props->find_subhalos) {
      props->sub_l_x = props->overdensity_ratio * props->l_x_absolute;
      props->sub_l_x2 = props->sub_l_x * props->sub_l_x;

      if (s->e->nodeID == 0)
        message("Subhalo linking length is set to %e [internal units].",
                props->sub_l_x);
    }

  } else {

    if (s->e->nodeID == 0)
      message(
          "Host linking length is set to %e [internal units] (%f of mean "
          "inter-DM-particle separation).",
          props->l_x, props->l_x_ratio);


    /* Define the subhalo linking length based on overdensity ratio. */
    if (props->find_subhalos) {
      props->sub_l_x_ratio = props->overdensity_ratio * props->l_x_ratio;
      props->sub_l_x = props->sub_l_x_ratio * props->l_x / props->l_x_ratio;
      props->sub_l_x2 = props->sub_l_x * props->sub_l_x;

      if (s->e->nodeID == 0)
      message(
          "Subhalo Linking length is set to %e [internal units] (%f of mean "
          "inter-DM-particle separation).",
          props->sub_l_x, props->sub_l_x_ratio);
    }
  }

  /* Now compute the velocity space linking lengths. */
  const double newt_G = s->e->physical_constants->const_newton_G;
  props->const_l_v = sqrt(newt_G / 2) *
    pow((4 * M_PI * props->host_overdensity * mean_matter_density) / 3,
        1.0 / 6.0);
  props->sub_const_l_v = sqrt(newt_G / 2) *
    pow((4 * M_PI * props->subhalo_overdensity * mean_matter_density) / 3,
        1.0 / 6.0);

  if (s->e->nodeID == 0)
      message(
          "Constant velocity Linking length (cosmology and halo independent) "
          "is set to %e [internal units].",
          props->const_l_v);

  if (s->e->nodeID == 0)
      message(
          "Constant velocity Linking length (cosmology and subhalo independent) "
          "is set to %e [internal units].",
          props->sub_const_l_v);

#ifdef WITH_MPI
  /* Check size of linking length against the top-level cell dimensions. */
  if (props->l_x2 > s->width[0] * s->width[0])
    error(
        "Linking length greater than the width of a top-level cell. Need to "
        "check more than one layer of top-level cells for links.");
#endif

  /* Allocate and initialise a group index array. */
  if (swift_memalign("fof_host_index", (void **)&props->host_index, 64,
                     s->nr_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle host indices for FOF search.");

    /* Allocate and initialise an array of halo objects. */
  if (swift_memalign("hosts", (void **)&props->hosts, halo_align,
                     s->nr_gparts * sizeof(struct halo)) != 0)
    error("Failed to allocate list of hosts for FOF search.");


  if (props->find_subhalos) {
    /* Allocate and initialise a group index array. */
    if (swift_memalign("fof_subhalo_index", (void **)&props->subhalo_index, 64,
                       s->nr_gparts * sizeof(size_t)) != 0)
      error("Failed to allocate list of particle subhalo indices for FOF search.");

    /* Allocate and initialise an array of halo objects. */
    if (swift_memalign("subhalos", (void **)&props->subhalos,
                       halo_align,
                       s->nr_gparts * sizeof(struct halo)) != 0)
      error("Failed to allocate list of subhalos for FOF search.");
    
  }

  ticks tic = getticks();

  /* Set initial host index */
  threadpool_map(&s->e->threadpool, fof_set_initial_group_index_mapper,
                 props->host_index, s->nr_gparts, sizeof(size_t),
                 threadpool_auto_chunk_size, props->host_index);

  if (verbose)
    message("Setting initial host index took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (props->find_subhalos) {
    tic = getticks();

    /* Set initial group index */
    threadpool_map(&s->e->threadpool, fof_set_initial_group_index_mapper,
                   props->subhalo_index, s->nr_gparts, sizeof(size_t),
                   threadpool_auto_chunk_size, props->subhalo_index);

    if (verbose)
      message("Setting initial subhalo index took: %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

#ifdef SWIFT_DEBUG_CHECKS
  ti_current = s->e->ti_current;
#endif

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - total_tic),
            clocks_getunit());
}

/**
 * @brief Run a halo finder search.
 *
 * @param e the engine
 * @param dump_results Are we writing group catalogues to output files?
 * @param dump_debug_results Are we writing a txt-file debug catalogue
 * (including BH seed info)?
 * @param seed_black_holes Are we seeding black holes?
 * @param foreign_buffers_allocated Are the foreign buffers currently
 * allocated?
 */
void engine_halo_finder(struct engine *e, const int dump_results,
                        const int dump_debug_results,
                        const int seed_black_holes,
                        const int foreign_buffers_allocated) {

#ifdef WITH_FOF

  struct fof_props *props = e->fof_properties;
  struct space *s = e->s;

  const ticks tic = getticks();

  /* Start by cleaning up the foreign buffers */
  if (foreign_buffers_allocated) {
#ifdef WITH_MPI
    space_free_foreign_parts(e->s, /*clear pointers=*/1);
#endif
  }

  /* Compute number of DM particles */
  const long long total_nr_baryons =
      e->total_nr_parts + e->total_nr_sparts + e->total_nr_bparts;
  const long long total_nr_dmparts =
      e->total_nr_gparts - e->total_nr_DM_background_gparts -
      e->total_nr_neutrino_gparts - total_nr_baryons;

  /* ---------------- First do the spatial FOF ---------------- */

  /* Initialise FOF parameters and allocate FOF arrays. */
  fof_allocate(e->s, total_nr_dmparts, e->fof_properties);

  /* Set current level. */
  props->current_level = fof_group;

  /* Make FOF tasks */
  engine_make_fof_tasks(e);

  /* and activate them. */
  engine_activate_fof_tasks(e);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

  /* Perform local FOF tasks. */
  engine_launch(e, "fof");

  /* Perform FOF search over foreign particles and
   * find groups which require black hole seeding.  */
  fof_search_tree(e->fof_properties, e->black_holes_properties,
                  e->physical_constants, e->cosmology, e->s, /*dump_results*/0,
                  /*dump_debug_results*/0, /*seed_black_holes*/0);

  /* Initialise halo finder parameters and allocate halo finder arrays. */
  halo_finder_allocate(e->s, props);

  /* ---------------- Run 6D host FOF ---------------- */

  /* Only bother if we have groups to refine. */
  /* TODO: We should write an empty file for consistency if
   * no groups are found. */
  if (props->num_groups > 0) {
    
    /* Set current level. */
    props->current_level = host_halo;

    /* Make the host halo tasks */
    engine_make_host_tasks(e);

    /* and activate them. */
    engine_activate_fof_tasks(e);

    /* Print the number of active tasks ? */
    if (e->verbose) engine_print_task_counts(e);

    /* Perform local host tasks. */
    engine_launch(e, "fof");

    /* Perform the host search to unify halos and compute properties. */
    halo_search_tree(props, e->physical_constants,  e->cosmology, e->s,
                     dump_results, dump_debug_results);

    /* Dump group data. */
    if (dump_results && !props->find_subhalos) {
#ifdef HAVE_HDF5
      write_fof_hdf5_catalogue(props, props->num_hosts, s->e,
                               /*is_halo_finder*/1);
#else
      error("Can't dump hdf5 catalogues with hdf5 switched off!");
#endif
    }
  }

  /* ---------------- Repeat for the subhalos ---------------- */

  /* Only search if subhalos are requested and there are hosts to search. */
  if (props->find_subhalos && props->num_hosts > 0) {

    /* Set current level to fof group. */
    props->current_level = sub_halo;
    
    /* Make the subhalo halo tasks */
    engine_make_subhalo_tasks(e);

    /* and activate them. */
    engine_activate_fof_tasks(e);

    /* Print the number of active tasks ? */
    if (e->verbose) engine_print_task_counts(e);

    /* Perform local host tasks. */
    engine_launch(e, "fof");

    /* Perform the subhalo search to unify halos and compute properties. */
    halo_search_tree(props, e->physical_constants,  e->cosmology, e->s,
                     dump_results, dump_debug_results);
    
    /* Dump group data. */
    if (dump_results) {
#ifdef HAVE_HDF5
      write_fof_hdf5_catalogue(props, props->num_subhalos, s->e,
                               /*is_halo_finder*/1);
#else
      error("Can't dump hdf5 catalogues with hdf5 switched off!");
#endif
    }
  }

  /* Restore the foreign buffers as they were*/
  if (foreign_buffers_allocated) {
#ifdef WITH_MPI
    engine_allocate_foreign_particles(e, /*fof=*/0);
#endif
  }

  /* Clean up arrays we don't want to carry around. */
  swift_free("fof_group_index", props->group_index);
  swift_free("fof_host_index", props->host_index);
  swift_free("fof_subhalo_index", props->subhalo_index);
  swift_free("groups", props->groups);
  swift_free("hosts", props->hosts);
  swift_free("subhalos", props->subhalos);
  props->group_index = NULL;
  props->host_index = NULL;
  props->subhalo_index = NULL;
  props->groups = NULL;
  props->hosts = NULL;
  props->subhalos = NULL;

  if (engine_rank == 0)
    message("Complete FOF search took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#else
  error("SWIFT was not compiled with FOF enabled!");
#endif
}
