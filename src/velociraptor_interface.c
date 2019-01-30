/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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
#include <unistd.h>

/* This object's header. */
#include "velociraptor_interface.h"

/* Local includes. */
#include "cooling.h"
#include "engine.h"
#include "hydro.h"
#include "swift_velociraptor_part.h"
#include "velociraptor_struct.h"

#ifdef HAVE_VELOCIRAPTOR

/**
 * @brief Structure for passing cosmological information to VELOCIraptor.
 */
struct cosmoinfo {

  /*! Current expansion factor of the Universe. (cosmology.a) */
  double atime;

  /*! Reduced Hubble constant (H0 / (100km/s/Mpc) (cosmology.h) */
  double littleh;

  /*! Matter density parameter (cosmology.Omega_m) */
  double Omega_m;

  /*! Radiation density parameter (cosmology.Omega_r) */
  double Omega_r;

  /*! Neutrino density parameter (0 in SWIFT) */
  double Omega_nu;

  /*! Neutrino density parameter (cosmology.Omega_k) */
  double Omega_k;

  /*! Baryon density parameter (cosmology.Omega_b) */
  double Omega_b;

  /*! Radiation constant density parameter (cosmology.Omega_lambda) */
  double Omega_Lambda;

  /*! Dark matter density parameter (cosmology.Omega_m - cosmology.Omega_b) */
  double Omega_cdm;

  /*! Dark-energy equation of state at the current time (cosmology.w)*/
  double w_de;
};

/**
 * @brief Structure for passing unit information to VELOCIraptor.
 */
struct unitinfo {

  /*! Length conversion factor to kpc. */
  double lengthtokpc;

  /*! Velocity conversion factor to km/s. */
  double velocitytokms;

  /*! Mass conversion factor to solar masses. */
  double masstosolarmass;

  /*! Potential conversion factor to (km/s)^2. */
  double energyperunitmass;

  /*! Newton's gravitationl constant (phys_const.const_newton_G)*/
  double gravity;

  /*! Hubble constant at the current redshift (cosmology.H) */
  double hubbleunit;
};

/**
 * @brief Structure to hold the location of a top-level cell.
 */
struct cell_loc {

  /*! Coordinates x,y,z */
  double loc[3];
};

/**
 * @brief Structure for passing simulation information to VELOCIraptor for a
 * given call.
 */
struct siminfo {

  /*! Size of periodic replications */
  double period;

  /*! Mass of the high-resolution DM particles in a zoom-in run. */
  double zoomhigresolutionmass;

  /*! Mean inter-particle separation of the DM particles */
  double interparticlespacing;

  /*! Spacial extent of the simulation volume */
  double spacedimension[3];

  /*! Number of top-level cells. */
  int numcells;

  /*! Locations of top-level cells. */
  struct cell_loc *cell_loc;

  /*! Top-level cell width. */
  double cellwidth[3];

  /*! Inverse of the top-level cell width. */
  double icellwidth[3];

  /*! Holds the node ID of each top-level cell. */
  int *cellnodeids;

  /*! Is this a cosmological simulation? */
  int icosmologicalsim;

  /*! Is this a zoom-in simulation? */
  int izoomsim;

  /*! Do we have DM particles? */
  int idarkmatter;

  /*! Do we have gas particles? */
  int igas;

  /*! Do we have star particles? */
  int istar;

  /*! Do we have BH particles? */
  int ibh;

  /*! Do we have other particles? */
  int iother;
};

/**
 * @brief Structure for group information back to swift
 */
struct groupinfo {

  /*! Index of a #gpart in the global array on this MPI rank */
  int index;

  /*! Group number of the #gpart. */
  long long groupID;
};

int InitVelociraptor(char *config_name, struct unitinfo unit_info,
                     struct siminfo sim_info, const int numthreads);

struct groupinfo *InvokeVelociraptor(
    const int snapnum, char *output_name, struct cosmoinfo cosmo_info,
    struct siminfo sim_info, const size_t num_gravity_parts,
    const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, const int *cell_node_ids,
    const int numthreads, const int return_group_flags,
    int *const num_in_groups);

#endif /* HAVE_VELOCIRAPTOR */

/**
 * @brief Temporary structure used for the data copy mapper.
 */
struct velociraptor_copy_data {
  const struct engine *e;
  struct swift_vel_part *swift_parts;
};

/**
 * @brief Mapper function to conver the #gpart into VELOCIraptor Particles.
 *
 * @param map_data The array of #gpart.
 * @param nr_gparts The number of #gpart.
 * @param extra_data Pointer to the #engine and to the array to fill.
 */
void velociraptor_convert_particles_mapper(void *map_data, int nr_gparts,
                                           void *extra_data) {

  /* Unpack the data */
  struct gpart *restrict gparts = (struct gpart *)map_data;
  struct velociraptor_copy_data *data =
      (struct velociraptor_copy_data *)extra_data;
  const struct engine *e = data->e;
  const struct space *s = e->s;
  struct swift_vel_part *swift_parts =
      data->swift_parts + (ptrdiff_t)(gparts - s->gparts);

  /* Handle on the other particle types */
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;
  const struct spart *sparts = s->sparts;

  /* Handle on the physics modules */
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cooling_function_data *cool_func = e->cooling_func;

  const float a_inv = e->cosmology->a_inv;

  /* Convert particle properties into VELOCIraptor units.
   * VELOCIraptor wants:
   * - Co-moving positions,
   * - Peculiar velocities,
   * - Co-moving potential,
   * - Physical internal energy (for the gas),
   * - Temperatures (for the gas).
   */
  for (int i = 0; i < nr_gparts; i++) {

    swift_parts[i].x[0] = gparts[i].x[0];
    swift_parts[i].x[1] = gparts[i].x[1];
    swift_parts[i].x[2] = gparts[i].x[2];

    swift_parts[i].v[0] = gparts[i].v_full[0] * a_inv;
    swift_parts[i].v[1] = gparts[i].v_full[1] * a_inv;
    swift_parts[i].v[2] = gparts[i].v_full[2] * a_inv;

    swift_parts[i].mass = gravity_get_mass(&gparts[i]);
    swift_parts[i].potential = gravity_get_comoving_potential(&gparts[i]);

    swift_parts[i].type = gparts[i].type;

    swift_parts[i].index = i;
#ifdef WITH_MPI
    swift_parts[i].task = e->nodeID;
#else
    swift_parts[i].task = 0;
#endif

    /* Set gas particle IDs from their hydro counterparts and set internal
     * energies. */
    switch (gparts[i].type) {

      case swift_type_gas: {
        const struct part *p = &parts[-gparts[i].id_or_neg_offset];
        const struct xpart *xp = &xparts[-gparts[i].id_or_neg_offset];

        swift_parts[i].id = parts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = hydro_get_drifted_physical_internal_energy(p, cosmo);
        swift_parts[i].T = cooling_get_temperature(phys_const, hydro_props, us,
                                                   cosmo, cool_func, p, xp);
      } break;

      case swift_type_stars:

        swift_parts[i].id = sparts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = 0.f;
        swift_parts[i].T = 0.f;
        break;

      case swift_type_dark_matter:

        swift_parts[i].id = gparts[i].id_or_neg_offset;
        swift_parts[i].u = 0.f;
        swift_parts[i].T = 0.f;
        break;

      default:
        error("Particle type not handled by VELOCIraptor.");
    }
  }
}

/**
 * @brief Initialise VELOCIraptor with configuration, units,
 * simulation info needed to run.
 *
 * @param e The #engine.
 */
void velociraptor_init(struct engine *e) {

#ifdef HAVE_VELOCIRAPTOR
  const ticks tic = getticks();

  /* Internal SWIFT units */
  const struct unit_system *swift_us = e->internal_units;

  /* CGS units and physical constants in CGS */
  struct unit_system cgs_us;
  units_init_cgs(&cgs_us);
  struct phys_const cgs_pc;
  phys_const_init(&cgs_us, /*params=*/NULL, &cgs_pc);

  /* Set unit conversions. */
  struct unitinfo unit_info;
  unit_info.lengthtokpc =
      units_cgs_conversion_factor(swift_us, UNIT_CONV_LENGTH) /
      (1000. * cgs_pc.const_parsec);
  unit_info.velocitytokms =
      units_cgs_conversion_factor(swift_us, UNIT_CONV_VELOCITY) / 1.0e5;
  unit_info.masstosolarmass =
      units_cgs_conversion_factor(swift_us, UNIT_CONV_MASS) /
      cgs_pc.const_solar_mass;
  unit_info.energyperunitmass =
      units_cgs_conversion_factor(swift_us, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
      (1.0e10);
  unit_info.gravity = e->physical_constants->const_newton_G;
  unit_info.hubbleunit = e->cosmology->H0 / e->cosmology->h;

  /* Gather some information about the simulation */
  struct siminfo sim_info;

  /* Are we running with cosmology? */
  if (e->policy & engine_policy_cosmology) {
    sim_info.icosmologicalsim = 1;
  } else {
    sim_info.icosmologicalsim = 0;
  }
  sim_info.izoomsim = 0;

  /* Tell VELOCIraptor what we have in the simulation */
  sim_info.idarkmatter = (e->total_nr_gparts - e->total_nr_parts > 0);
  sim_info.igas = (e->policy & engine_policy_hydro);
  sim_info.istar = (e->policy & engine_policy_stars);
  sim_info.ibh = 0;  // sim_info.ibh = (e->policy&engine_policy_bh);
  sim_info.iother = 0;

  /* Be nice, talk! */
  if (e->verbose) {
    message("VELOCIraptor conf: Length conversion factor: %e",
            unit_info.lengthtokpc);
    message("VELOCIraptor conf: Velocity conversion factor: %e",
            unit_info.velocitytokms);
    message("VELOCIraptor conf: Mass conversion factor: %e",
            unit_info.masstosolarmass);
    message("VELOCIraptor conf: Internal energy conversion factor: %e",
            unit_info.energyperunitmass);
    message("VELOCIraptor conf: G: %e", unit_info.gravity);
    message("VELOCIraptor conf: H0/h: %e", unit_info.hubbleunit);
    message("VELOCIraptor conf: Config file name: %s", e->stf_config_file_name);
    message("VELOCIraptor conf: Cosmological Simulation: %d",
            sim_info.icosmologicalsim);
  }

  /* Initialise VELOCIraptor. */
  if (InitVelociraptor(e->stf_config_file_name, unit_info, sim_info,
                       e->nr_threads) != 1)
    error("VELOCIraptor initialisation failed.");

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT not configure to run with VELOCIraptor.");
#endif /* HAVE_VELOCIRAPTOR */
}

/**
 * @brief Run VELOCIraptor with current particle data.
 *
 * @param e The #engine.
 * @param linked_with_snap Are we running at the same time as a snapshot dump?
 */
void velociraptor_invoke(struct engine *e, const int linked_with_snap) {

#ifdef HAVE_VELOCIRAPTOR

  /* Handle on the particles */
  const struct space *s = e->s;
  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_parts = s->nr_parts;
  const size_t nr_sparts = s->nr_sparts;
  const int nr_cells = s->nr_cells;
  const struct cell *cells_top = s->cells_top;

  /* Allow thread to run on any core for the duration of the call to
   * VELOCIraptor so that  when OpenMP threads are spawned
   * they can run on any core on the processor. */
  const int nr_cores = sysconf(_SC_NPROCESSORS_ONLN);
  pthread_t thread = pthread_self();

  /* Set affinity mask to include all cores on the CPU for VELOCIraptor. */
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  for (int j = 0; j < nr_cores; j++) CPU_SET(j, &cpuset);
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);

  /* Set cosmology information for this point in time */
  struct cosmoinfo cosmo_info;
  cosmo_info.atime = e->cosmology->a;
  cosmo_info.littleh = e->cosmology->h;
  cosmo_info.Omega_m = e->cosmology->Omega_m;
  cosmo_info.Omega_b = e->cosmology->Omega_b;
  cosmo_info.Omega_r = e->cosmology->Omega_r;
  cosmo_info.Omega_k = e->cosmology->Omega_k;
  cosmo_info.Omega_nu = 0.;
  cosmo_info.Omega_Lambda = e->cosmology->Omega_lambda;
  cosmo_info.Omega_cdm = e->cosmology->Omega_m - e->cosmology->Omega_b;
  cosmo_info.w_de = e->cosmology->w;

  /* Report the cosmo info we use */
  if (e->verbose) {
    message("VELOCIraptor conf: Scale factor: %e", cosmo_info.atime);
    message("VELOCIraptor conf: Little h: %e", cosmo_info.littleh);
    message("VELOCIraptor conf: Omega_m: %e", cosmo_info.Omega_m);
    message("VELOCIraptor conf: Omega_b: %e", cosmo_info.Omega_b);
    message("VELOCIraptor conf: Omega_Lambda: %e", cosmo_info.Omega_Lambda);
    message("VELOCIraptor conf: Omega_cdm: %e", cosmo_info.Omega_cdm);
    message("VELOCIraptor conf: w_de: %e", cosmo_info.w_de);
  }

  /* Update the simulation information */
  struct siminfo sim_info;

  /* Period of the box (Note we assume a cubic box!) */
  if (e->s->periodic) {
    sim_info.period = s->dim[0];
  } else {
    sim_info.period = 0.0;
  }

  /* Tell VELOCIraptor this is not a zoom-in simulation */
  sim_info.zoomhigresolutionmass = -1.0;

  /* Are we running with cosmology? */
  if (e->policy & engine_policy_cosmology) {
    sim_info.icosmologicalsim = 1;
    sim_info.izoomsim = 0;
    const size_t total_nr_baryons = e->total_nr_parts + e->total_nr_sparts;
    const size_t total_nr_dmparts = e->total_nr_gparts - total_nr_baryons;
    sim_info.interparticlespacing = sim_info.period / cbrt(total_nr_dmparts);
  } else {
    sim_info.icosmologicalsim = 0;
    sim_info.izoomsim = 0;
    sim_info.interparticlespacing = -1;
  }

  /* Set the spatial extent of the simulation volume */
  sim_info.spacedimension[0] = s->dim[0];
  sim_info.spacedimension[1] = s->dim[1];
  sim_info.spacedimension[2] = s->dim[2];

  /* Store number of top-level cells */
  sim_info.numcells = s->nr_cells;

  /* Size and inverse size of the top-level cells in VELOCIraptor units */
  sim_info.cellwidth[0] = s->cells_top[0].width[0];
  sim_info.cellwidth[1] = s->cells_top[0].width[1];
  sim_info.cellwidth[2] = s->cells_top[0].width[2];
  sim_info.icellwidth[0] = s->iwidth[0];
  sim_info.icellwidth[1] = s->iwidth[1];
  sim_info.icellwidth[2] = s->iwidth[2];

  ticks tic = getticks();

  /* Allocate and populate array of cell node IDs and positions. */
  int *cell_node_ids = NULL;
  if (posix_memalign((void **)&sim_info.cell_loc, SWIFT_STRUCT_ALIGNMENT,
                     s->nr_cells * sizeof(struct cell_loc)) != 0)
    error("Failed to allocate top-level cell locations for VELOCIraptor.");
  if (posix_memalign((void **)&cell_node_ids, SWIFT_STRUCT_ALIGNMENT,
                     nr_cells * sizeof(int)) != 0)
    error("Failed to allocate list of cells node IDs for VELOCIraptor.");

  for (int i = 0; i < s->nr_cells; i++) {
    cell_node_ids[i] = cells_top[i].nodeID;

    sim_info.cell_loc[i].loc[0] = cells_top[i].loc[0];
    sim_info.cell_loc[i].loc[1] = cells_top[i].loc[1];
    sim_info.cell_loc[i].loc[2] = cells_top[i].loc[2];
  }

  if (e->verbose) {
    message("VELOCIraptor conf: Space dimensions: (%e,%e,%e)",
            sim_info.spacedimension[0], sim_info.spacedimension[1],
            sim_info.spacedimension[2]);
    message("VELOCIraptor conf: No. of top-level cells: %d", sim_info.numcells);
    message(
        "VELOCIraptor conf: Top-level cell locations range: (%e,%e,%e) -> "
        "(%e,%e,%e)",
        sim_info.cell_loc[0].loc[0], sim_info.cell_loc[0].loc[1],
        sim_info.cell_loc[0].loc[2],
        sim_info.cell_loc[sim_info.numcells - 1].loc[0],
        sim_info.cell_loc[sim_info.numcells - 1].loc[1],
        sim_info.cell_loc[sim_info.numcells - 1].loc[2]);
  }

  /* Report timing */
  if (e->verbose)
    message("VR Collecting top-level cell info took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Mention the number of particles being sent */
  if (e->verbose)
    message(
        "VELOCIraptor conf: MPI rank %d sending %zu gparts to VELOCIraptor.",
        engine_rank, nr_gparts);

  /* Append base name with the current output number */
  char outputFileName[PARSER_MAX_LINE_SIZE + 128];

  /* What should the filename be? */
  if (linked_with_snap) {
    snprintf(outputFileName, PARSER_MAX_LINE_SIZE + 128,
             "stf_%s_%04i.VELOCIraptor", e->snapshot_base_name,
             e->snapshot_output_count);
  } else {
    snprintf(outputFileName, PARSER_MAX_LINE_SIZE + 128, "%s_%04i.VELOCIraptor",
             e->stf_base_name, e->stf_output_count);
  }

  /* What is the snapshot number? */
  int snapnum;
  if (linked_with_snap) {
    snapnum = e->snapshot_output_count;
  } else {
    snapnum = e->stf_output_count;
  }

  tic = getticks();

  /* Allocate and populate an array of swift_vel_parts to be passed to
   * VELOCIraptor. */
  struct swift_vel_part *swift_parts = NULL;
  if (posix_memalign((void **)&swift_parts, part_align,
                     nr_gparts * sizeof(struct swift_vel_part)) != 0)
    error("Failed to allocate array of particles for VELOCIraptor.");

  struct velociraptor_copy_data copy_data = {e, swift_parts};
  threadpool_map(&e->threadpool, velociraptor_convert_particles_mapper,
                 s->gparts, nr_gparts, sizeof(struct gpart), 0, &copy_data);

  /* Report timing */
  if (e->verbose)
    message("VR Collecting particle info took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Values returned by VELOCIRaptor */
  int num_gparts_in_groups = -1;
  struct groupinfo *group_info = NULL;

  /* Call VELOCIraptor. */
  group_info = (struct groupinfo *)InvokeVelociraptor(
      snapnum, outputFileName, cosmo_info, sim_info, nr_gparts, nr_parts,
      nr_sparts, swift_parts, cell_node_ids, e->nr_threads, linked_with_snap,
      &num_gparts_in_groups);

  /* Check that the ouput is valid */
  if (linked_with_snap && group_info == NULL && num_gparts_in_groups < 0) {
    error("Exiting. Call to VELOCIraptor failed on rank: %d.", e->nodeID);
  }
  if (!linked_with_snap && group_info != NULL) {
    error("VELOCIraptor returned an array whilst it should not have.");
  }

  /* Report timing */
  if (e->verbose)
    message("VR Invokation of velociraptor took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Assign the group IDs back to the gparts */
  if (linked_with_snap) {

    if (posix_memalign((void **)&s->gpart_group_data, part_align,
                       nr_gparts * sizeof(struct velociraptor_gpart_data)) != 0)
      error("Failed to allocate array of gpart data for VELOCIraptor i/o.");

    struct velociraptor_gpart_data *data = s->gpart_group_data;

    /* Zero the array (gparts not in groups have an ID of 0) */
    bzero(data, nr_gparts * sizeof(struct velociraptor_gpart_data));

    /* Copy the data at the right place */
    for (int i = 0; i < num_gparts_in_groups; i++) {
      data[group_info[i].index].groupID = group_info[i].groupID;
    }

    /* Report timing */
    if (e->verbose)
      message("VR Copying group information back took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());

    /* Free the array returned by VELOCIraptor */
    free(group_info);
  }

  /* Reset the pthread affinity mask after VELOCIraptor returns. */
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), engine_entry_affinity());

  /* Increase output counter (if not linked with snapshots) */
  if (!linked_with_snap) e->stf_output_count++;

#else
  error("SWIFT not configure to run with VELOCIraptor.");
#endif /* HAVE_VELOCIRAPTOR */
}
