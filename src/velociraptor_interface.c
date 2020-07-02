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
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/* This object's header. */
#include "velociraptor_interface.h"

/* Local includes. */
#include "cooling.h"
#include "engine.h"
#include "hydro.h"
#include "swift_velociraptor_part.h"
#include "threadpool.h"
#include "velociraptor_struct.h"

#ifdef HAVE_VELOCIRAPTOR

/**
 * @brief Structure for passing cosmological information to VELOCIraptor.
 *
 * This should match the structure cosmoinfo in the file src/swiftinterface.h
 * in the VELOCIraptor code.
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
 *
 * This should match the structure unitinfo in the file src/swiftinterface.h
 * in the VELOCIraptor code.
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
 *
 * This should match the structure siminfo in the file src/swiftinterface.h
 * in the VELOCIraptor code.
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

  /*! Number of top-level cells. */
  int numcellsperdim;

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

#ifdef HAVE_VELOCIRAPTOR_WITH_NOMASS
  /*! Mass of the DM particles */
  double mass_uniform_box;
#endif
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
  const struct bpart *bparts = s->bparts;

  /* Handle on the physics modules */
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cooling_function_data *cool_func = e->cooling_func;

  const float a_inv = e->cosmology->a_inv;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double pos_dithering[3] = {s->pos_dithering[0], s->pos_dithering[1],
                                   s->pos_dithering[2]};

  /* Convert particle properties into VELOCIraptor units.
   * VELOCIraptor wants:
   * - Un-dithered co-moving positions,
   * - Peculiar velocities,
   * - Co-moving potential,
   * - Physical internal energy (for the gas),
   * - Temperatures (for the gas).
   */
  for (int i = 0; i < nr_gparts; i++) {

    if (periodic) {
      swift_parts[i].x[0] =
          box_wrap(gparts[i].x[0] - pos_dithering[0], 0.0, dim[0]);
      swift_parts[i].x[1] =
          box_wrap(gparts[i].x[1] - pos_dithering[1], 0.0, dim[1]);
      swift_parts[i].x[2] =
          box_wrap(gparts[i].x[2] - pos_dithering[2], 0.0, dim[2]);
    } else {
      swift_parts[i].x[0] = gparts[i].x[0];
      swift_parts[i].x[1] = gparts[i].x[1];
      swift_parts[i].x[2] = gparts[i].x[2];
    }

    swift_parts[i].v[0] = gparts[i].v_full[0] * a_inv;
    swift_parts[i].v[1] = gparts[i].v_full[1] * a_inv;
    swift_parts[i].v[2] = gparts[i].v_full[2] * a_inv;

#ifndef HAVE_VELOCIRAPTOR_WITH_NOMASS
    swift_parts[i].mass = gravity_get_mass(&gparts[i]);
#endif

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

      case swift_type_black_hole:

        swift_parts[i].id = bparts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = 0.f;
        swift_parts[i].T = 0.f;
        break;

      case swift_type_dark_matter:

        swift_parts[i].id = gparts[i].id_or_neg_offset;
        swift_parts[i].u = 0.f;
        swift_parts[i].T = 0.f;
        break;

      case swift_type_dark_matter_background:

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

  /* Are we running a zoom? */
  if (e->s->with_DM_background) {
    sim_info.izoomsim = 1;
  } else {
    sim_info.izoomsim = 0;
  }

  /* Tell VELOCIraptor what we have in the simulation */
  sim_info.idarkmatter = (e->total_nr_gparts - e->total_nr_parts > 0);
  sim_info.igas = (e->policy & engine_policy_hydro);
  sim_info.istar = (e->policy & engine_policy_stars);
  sim_info.ibh = (e->policy & engine_policy_black_holes);
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
  error("SWIFT not configured to run with VELOCIraptor.");
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

  /* Start by freeing some of the unnecessary memory to give VR some breathing
     space */
#ifdef WITH_MPI
  space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);
#endif

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

    /* Are we running a zoom? */
    if (e->s->with_DM_background) {
      sim_info.izoomsim = 1;
    } else {
      sim_info.izoomsim = 0;
    }

    /* Collect the mass of the non-background gpart */
    double high_res_DM_mass = 0.;
    for (size_t i = 0; i < e->s->nr_gparts; ++i) {
      const struct gpart *gp = &e->s->gparts[i];
      if (gp->type == swift_type_dark_matter &&
          gp->time_bin != time_bin_inhibited &&
          gp->time_bin != time_bin_not_created) {
        high_res_DM_mass = gp->mass;
        break;
      }
    }

#ifdef WITH_MPI
    /* We need to all-reduce this in case one of the nodes had 0 DM particles.
     */
    MPI_Allreduce(MPI_IN_PLACE, &high_res_DM_mass, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
#endif

    const double Omega_m = e->cosmology->Omega_m;
    const double Omega_b = e->cosmology->Omega_b;
    const double critical_density_0 = e->cosmology->critical_density_0;

    /* Linking length based on the mean DM inter-particle separation
     * in the zoom region and assuming the mean density of the Universe
     * is used in the zoom region. */
    double mean_matter_density;
    if (s->with_hydro)
      mean_matter_density = (Omega_m - Omega_b) * critical_density_0;
    else
      mean_matter_density = Omega_m * critical_density_0;

    sim_info.interparticlespacing =
        cbrt(high_res_DM_mass / mean_matter_density);

  } else {
    sim_info.izoomsim = 0;
    sim_info.icosmologicalsim = 0;
    sim_info.interparticlespacing = -1.;
  }

#ifdef HAVE_VELOCIRAPTOR_WITH_NOMASS
  /* Assume all particles have the same mass */
  double DM_mass = 0.;
  for (size_t i = 0; i < e->s->nr_gparts; ++i) {
    const struct gpart *gp = &e->s->gparts[i];
    if (gp->time_bin != time_bin_inhibited &&
        gp->time_bin != time_bin_not_created) {
      DM_mass = gp->mass;
      break;
    }
  }
  sim_info.mass_uniform_box = DM_mass;
#endif

  /* Set the spatial extent of the simulation volume */
  sim_info.spacedimension[0] = s->dim[0];
  sim_info.spacedimension[1] = s->dim[1];
  sim_info.spacedimension[2] = s->dim[2];

  /* Store number of top-level cells */
  sim_info.numcells = s->nr_cells;
  sim_info.numcellsperdim = s->cdim[0]; /* We assume a cubic box! */
  if (s->cdim[0] != s->cdim[1] || s->cdim[0] != s->cdim[2])
    error("Trying to run VR on a non-cubic number of top-level cells");

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
  if (swift_memalign("VR.cell_loc", (void **)&sim_info.cell_loc,
                     SWIFT_STRUCT_ALIGNMENT,
                     s->nr_cells * sizeof(struct cell_loc)) != 0)
    error("Failed to allocate top-level cell locations for VELOCIraptor.");
  if (swift_memalign("VR.cell_nodeID", (void **)&cell_node_ids,
                     SWIFT_STRUCT_ALIGNMENT, nr_cells * sizeof(int)) != 0)
    error("Failed to allocate list of cells node IDs for VELOCIraptor.");

  for (int i = 0; i < s->nr_cells; i++) {
    cell_node_ids[i] = cells_top[i].nodeID;

    if (s->periodic) {
      sim_info.cell_loc[i].loc[0] =
          box_wrap(cells_top[i].loc[0] - s->pos_dithering[0], 0.0, s->dim[0]);
      sim_info.cell_loc[i].loc[1] =
          box_wrap(cells_top[i].loc[1] - s->pos_dithering[1], 0.0, s->dim[1]);
      sim_info.cell_loc[i].loc[2] =
          box_wrap(cells_top[i].loc[2] - s->pos_dithering[2], 0.0, s->dim[2]);
    } else {
      sim_info.cell_loc[i].loc[0] = cells_top[i].loc[0];
      sim_info.cell_loc[i].loc[1] = cells_top[i].loc[1];
      sim_info.cell_loc[i].loc[2] = cells_top[i].loc[2];
    }
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

  /* Generate directory name for this output - start with snapshot directory, if
   * specified */
  char outputDirName[FILENAME_BUFFER_SIZE] = "";
  if (strcmp(e->snapshot_subdir, engine_default_snapshot_subdir) != 0) {
    if (snprintf(outputDirName, FILENAME_BUFFER_SIZE, "%s/",
                 e->snapshot_subdir) >= FILENAME_BUFFER_SIZE) {
      error("FILENAME_BUFFER_SIZE is to small for snapshot directory name!");
    }
    if (engine_rank == 0) io_make_snapshot_subdir(e->snapshot_subdir);
#ifdef WITH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  /* Then create output-specific subdirectory if necessary */
  char subDirName[FILENAME_BUFFER_SIZE] = "";
  if (strcmp(e->stf_subdir_per_output, engine_default_stf_subdir_per_output) !=
      0) {
    if (snprintf(subDirName, FILENAME_BUFFER_SIZE, "%s%s_%04i/", outputDirName,
                 e->stf_subdir_per_output,
                 e->stf_output_count) >= FILENAME_BUFFER_SIZE) {
      error(
          "FILENAME_BUFFER_SIZE is to small for Velociraptor directory name!");
    }
    if (engine_rank == 0) io_make_snapshot_subdir(subDirName);
#ifdef WITH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  } else {
    /* Not making separate directories so subDirName=outputDirName */
    strncpy(subDirName, outputDirName, FILENAME_BUFFER_SIZE);
  }

  /* What should the filename be? */
  char outputFileName[FILENAME_BUFFER_SIZE];
  if (snprintf(outputFileName, FILENAME_BUFFER_SIZE, "%s%s_%04i.VELOCIraptor",
               subDirName, e->stf_base_name,
               e->stf_output_count) >= FILENAME_BUFFER_SIZE) {
    error("FILENAME_BUFFER_SIZE is too small for Velociraptor file name!");
  }

  tic = getticks();

  /* Allocate and populate an array of swift_vel_parts to be passed to
   * VELOCIraptor. */
  struct swift_vel_part *swift_parts = NULL;
  if (swift_memalign("VR.parts", (void **)&swift_parts, part_align,
                     nr_gparts * sizeof(struct swift_vel_part)) != 0)
    error("Failed to allocate array of particles for VELOCIraptor.");

  struct velociraptor_copy_data copy_data = {e, swift_parts};
  threadpool_map(&e->threadpool, velociraptor_convert_particles_mapper,
                 s->gparts, nr_gparts, sizeof(struct gpart),
                 threadpool_auto_chunk_size, &copy_data);

  /* Report timing */
  if (e->verbose)
    message("VR Collecting particle info took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Values returned by VELOCIRaptor */
  int num_gparts_in_groups = -1;
  struct groupinfo *group_info = NULL;

#ifdef SWIFT_MEMUSE_REPORTS
  char report_filename[60];
  sprintf(report_filename, "memuse-VR-report-rank%d-step%d.txt", e->nodeID,
          e->step);
  memuse_log_dump(report_filename);
#endif

  /* Call VELOCIraptor. */
  group_info = (struct groupinfo *)InvokeVelociraptor(
      e->stf_output_count, outputFileName, cosmo_info, sim_info, nr_gparts,
      nr_parts, nr_sparts, swift_parts, cell_node_ids, e->nr_threads,
      linked_with_snap, &num_gparts_in_groups);

  /* Report that the memory was freed */
  memuse_log_allocation("VR.cell_loc", sim_info.cell_loc, 0, 0);
  memuse_log_allocation("VR.cell_nodeID", cell_node_ids, 0, 0);
  memuse_log_allocation("VR.parts", swift_parts, 0, 0);

  /* Check that the ouput is valid */
  if (linked_with_snap && group_info == NULL && num_gparts_in_groups < 0) {
    error("Exiting. Call to VELOCIraptor failed on rank: %d.", e->nodeID);
  }
  if (!linked_with_snap && group_info != NULL) {
    error("VELOCIraptor returned an array whilst it should not have.");
  }

  /* Report timing */
  if (e->verbose)
    message("VR Invocation of velociraptor took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Assign the group IDs back to the gparts */
  if (linked_with_snap) {

    if (swift_memalign("VR.group_data", (void **)&s->gpart_group_data,
                       part_align,
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
    swift_free("VR.group_data", group_info);
  }

  /* Reset the pthread affinity mask after VELOCIraptor returns. */
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), engine_entry_affinity());

  /* Increase output counter (if not linked with snapshot) */
  if (!linked_with_snap) e->stf_output_count++;

  /* Record we have ran stf this timestep */
  e->stf_this_timestep = 1;

  /* Reallocate the memory that was freed earlier */
#ifdef WITH_MPI

  engine_allocate_foreign_particles(e);
#endif

#else
  error("SWIFT not configured to run with VELOCIraptor.");
#endif /* HAVE_VELOCIRAPTOR */
}
