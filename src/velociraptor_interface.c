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
#include "engine.h"
#include "hydro.h"
#include "swift_velociraptor_part.h"

#ifdef HAVE_VELOCIRAPTOR

/* Structure for passing cosmological information to VELOCIraptor. */
struct cosmoinfo {

  /*! Current expansion factor of the Universe. (cosmology.a) */
  double atime;

  /*! Reduced Hubble constant (H0 / (100km/s/Mpc) (cosmology.h) */
  double littleh;

  /*! Matter density parameter (cosmology.Omega_m) */
  double Omega_m;

  /*! Baryon density parameter (cosmology.Omega_b) */
  double Omega_b;

  /*! Radiation constant density parameter (cosmology.Omega_lambda) */
  double Omega_Lambda;

  /*! Dark matter density parameter (cosmology.Omega_m - cosmology.Omega_b) */
  double Omega_cdm;

  /*! Dark-energy equation of state at the current time (cosmology.w)*/
  double w_de;
};

/* Structure for passing unit information to VELOCIraptor. */
struct unitinfo {

  /* Length conversion factor to kpc. */
  double lengthtokpc;

  /* Velocity conversion factor to km/s. */
  double velocitytokms;

  /* Mass conversion factor to solar masses. */
  double masstosolarmass;

  /* Potential conversion factor. */
  double energyperunitmass;

  /*! Newton's gravitationl constant (phys_const.const_newton_G)*/
  double gravity;

  /*! Hubble constant at the current redshift (cosmology.H) */
  double hubbleunit;
};

/* Structure to hold the location of a top-level cell. */
struct cell_loc {

  /* Coordinates x,y,z */
  double loc[3];
};

/* Structure for passing simulation information to VELOCIraptor. */
struct siminfo {
  double period, zoomhigresolutionmass, interparticlespacing, spacedimension[3];

  /* Number of top-cells. */
  int numcells;

  /*! Locations of top-level cells. */
  struct cell_loc *cell_loc;

  /*! Top-level cell width. */
  double cellwidth[3];

  /*! Inverse of the top-level cell width. */
  double icellwidth[3];

  int icosmologicalsim;
};

/* VELOCIraptor interface. */
int InitVelociraptor(char *config_name, char *output_name,
                     struct cosmoinfo cosmo_info, struct unitinfo unit_info,
                     struct siminfo sim_info);

int InvokeVelociraptor(const size_t num_gravity_parts,
                       const size_t num_hydro_parts, const int snapnum,
                       struct swift_vel_part *swift_parts,
                       const int *cell_node_ids, char *output_name);

#endif /* HAVE_VELOCIRAPTOR */

/**
 * @brief Initialise VELOCIraptor with input and output file names along with
 * cosmological info needed to run.
 *
 * @param e The #engine.
 * @param linked_with_snap Are we running at the same time as a snapshot dump?
 */
void velociraptor_init(struct engine *e, const int linked_with_snap) {

#ifdef HAVE_VELOCIRAPTOR

  struct space *s = e->s;
  struct cosmoinfo cosmo_info;
  struct unitinfo unit_info;
  struct siminfo sim_info;
  const ticks tic = getticks();

  /* Set cosmological constants. */
  cosmo_info.atime = e->cosmology->a;
  cosmo_info.littleh = e->cosmology->h;
  cosmo_info.Omega_m = e->cosmology->Omega_m;
  cosmo_info.Omega_b = e->cosmology->Omega_b;
  cosmo_info.Omega_Lambda = e->cosmology->Omega_lambda;
  cosmo_info.Omega_cdm = e->cosmology->Omega_m - e->cosmology->Omega_b;
  cosmo_info.w_de = e->cosmology->w;

  message("Scale factor: %e", cosmo_info.atime);
  message("Little h: %e", cosmo_info.littleh);
  message("Omega_m: %e", cosmo_info.Omega_m);
  message("Omega_b: %e", cosmo_info.Omega_b);
  message("Omega_Lambda: %e", cosmo_info.Omega_Lambda);
  message("Omega_cdm: %e", cosmo_info.Omega_cdm);
  message("w_de: %e", cosmo_info.w_de);

  if (e->cosmology->w != -1.)
    error("w_de is not 1. It is: %lf", e->cosmology->w);

  /* Set unit conversions. */
  unit_info.lengthtokpc = 1.0;
  unit_info.velocitytokms = 1.0;
  unit_info.masstosolarmass = 1.0;
  unit_info.energyperunitmass = 1.0;
  unit_info.gravity = e->physical_constants->const_newton_G;
  unit_info.hubbleunit = e->cosmology->H0 / e->cosmology->h;

  message("Length conversion factor: %e", unit_info.lengthtokpc);
  message("Velocity conversion factor: %e", unit_info.velocitytokms);
  message("Mass conversion factor: %e", unit_info.masstosolarmass);
  message("Potential conversion factor: %e", unit_info.energyperunitmass);
  message("G: %e", unit_info.gravity);
  message("H: %e", unit_info.hubbleunit);

  /* TODO: Find the total number of DM particles when running with star
   * particles and BHs. */
  const int total_nr_dmparts = e->total_nr_gparts - e->total_nr_parts;

  /* Set periodicity information. */
  if (e->s->periodic) {
    /* Physical size of box in VELOCIraptor units (kpc). */
    sim_info.period = unit_info.lengthtokpc * s->dim[0];
  } else {
    sim_info.period = 0.0;
  }

  /* Are we running with cosmology? */
  if (e->policy & engine_policy_cosmology) {
    sim_info.icosmologicalsim = 1;
  } else {
    sim_info.icosmologicalsim = 0;
  }

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

  /* Gather the rest of the information */
  sim_info.zoomhigresolutionmass = -1.0; /* Placeholder. */
  sim_info.interparticlespacing = sim_info.period / cbrt(total_nr_dmparts);
  sim_info.spacedimension[0] = unit_info.lengthtokpc * s->dim[0];
  sim_info.spacedimension[1] = unit_info.lengthtokpc * s->dim[1];
  sim_info.spacedimension[2] = unit_info.lengthtokpc * s->dim[2];
  sim_info.numcells = s->nr_cells;

  sim_info.cellwidth[0] = unit_info.lengthtokpc * s->cells_top[0].width[0];
  sim_info.cellwidth[1] = unit_info.lengthtokpc * s->cells_top[0].width[1];
  sim_info.cellwidth[2] = unit_info.lengthtokpc * s->cells_top[0].width[2];

  sim_info.icellwidth[0] = s->iwidth[0] / unit_info.lengthtokpc;
  sim_info.icellwidth[1] = s->iwidth[1] / unit_info.lengthtokpc;
  sim_info.icellwidth[2] = s->iwidth[2] / unit_info.lengthtokpc;

  /* Only allocate cell location array on first call to velociraptor_init(). */
  if (e->cell_loc == NULL) {

    /* Allocate and populate top-level cell locations. */
    if (posix_memalign((void **)&(e->cell_loc), 32,
                       s->nr_cells * sizeof(struct cell_loc)) != 0)
      error("Failed to allocate top-level cell locations for VELOCIraptor.");

    for (int i = 0; i < s->nr_cells; i++) {
      e->cell_loc[i].loc[0] = unit_info.lengthtokpc * s->cells_top[i].loc[0];
      e->cell_loc[i].loc[1] = unit_info.lengthtokpc * s->cells_top[i].loc[1];
      e->cell_loc[i].loc[2] = unit_info.lengthtokpc * s->cells_top[i].loc[2];
    }
  }

  sim_info.cell_loc = e->cell_loc;

  message("Config file name: %s", e->stf_config_file_name);
  message("Period: %e", sim_info.period);
  message("Zoom high res mass: %e", sim_info.zoomhigresolutionmass);
  message("Inter-particle spacing: %e", sim_info.interparticlespacing);
  message("Cosmological: %d", sim_info.icosmologicalsim);
  message("Space dimensions: (%e,%e,%e)", sim_info.spacedimension[0],
          sim_info.spacedimension[1], sim_info.spacedimension[2]);
  message("No. of top-level cells: %d", sim_info.numcells);
  message("Top-level cell locations range: (%e,%e,%e) -> (%e,%e,%e)",
          sim_info.cell_loc[0].loc[0], sim_info.cell_loc[0].loc[1],
          sim_info.cell_loc[0].loc[2],
          sim_info.cell_loc[sim_info.numcells - 1].loc[0],
          sim_info.cell_loc[sim_info.numcells - 1].loc[1],
          sim_info.cell_loc[sim_info.numcells - 1].loc[2]);

  /* Initialise VELOCIraptor. */
  if (!InitVelociraptor(e->stf_config_file_name, outputFileName, cosmo_info,
                        unit_info, sim_info)) {
    error("Exiting. VELOCIraptor initialisation failed.");
  }

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

  const struct space *s = e->s;
  struct gpart *gparts = s->gparts;
  struct part *parts = s->parts;
  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_hydro_parts = s->nr_parts;
  const int nr_cells = s->nr_cells;

  const ticks tic = getticks();

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

  /* Allocate and populate array of cell node IDs. */
  int *cell_node_ids = NULL;
  if (posix_memalign((void **)&cell_node_ids, 32, nr_cells * sizeof(int)) != 0)
    error("Failed to allocate list of cells node IDs for VELOCIraptor.");

  for (int i = 0; i < nr_cells; i++) cell_node_ids[i] = s->cells_top[i].nodeID;

  /* Mention the number of particles being sent */
  if (e->verbose)
    message("MPI rank %d sending %zu gparts to VELOCIraptor.", engine_rank,
            nr_gparts);

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

  /* Allocate and populate an array of swift_vel_parts to be passed to
   * VELOCIraptor. */
  struct swift_vel_part *swift_parts = NULL;
  if (posix_memalign((void **)&swift_parts, part_align,
                     nr_gparts * sizeof(struct swift_vel_part)) != 0)
    error("Failed to allocate array of particles for VELOCIraptor.");

  bzero(swift_parts, nr_gparts * sizeof(struct swift_vel_part));

  const float a_inv = e->cosmology->a_inv;

  /* Convert particle properties into VELOCIraptor units */
  for (size_t i = 0; i < nr_gparts; i++) {

    swift_parts[i].x[0] = gparts[i].x[0];
    swift_parts[i].x[1] = gparts[i].x[1];
    swift_parts[i].x[2] = gparts[i].x[2];

    swift_parts[i].v[0] = gparts[i].v_full[0] * a_inv;
    swift_parts[i].v[1] = gparts[i].v_full[1] * a_inv;
    swift_parts[i].v[2] = gparts[i].v_full[2] * a_inv;

    swift_parts[i].mass = gravity_get_mass(&gparts[i]);
    swift_parts[i].potential = gravity_get_comoving_potential(&gparts[i]);

    swift_parts[i].type = gparts[i].type;

    /* Set gas particle IDs from their hydro counterparts and set internal
     * energies. */
    switch (gparts[i].type) {

      case swift_type_gas:

        swift_parts[i].id = parts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = hydro_get_drifted_physical_internal_energy(
            &parts[-gparts[i].id_or_neg_offset], e->cosmology);
        break;

      case swift_type_dark_matter:

        swift_parts[i].id = gparts[i].id_or_neg_offset;
        swift_parts[i].u = 0.f;
        break;

      default:
        error("Particle type not handled by VELOCIraptor.");
    }
  }

  /* Call VELOCIraptor. */
  if (!InvokeVelociraptor(nr_gparts, nr_hydro_parts, e->stf_output_count,
                          swift_parts, cell_node_ids, outputFileName))
    error("Exiting. Call to VELOCIraptor failed on rank: %d.", e->nodeID);

  /* Reset the pthread affinity mask after VELOCIraptor returns. */
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), engine_entry_affinity());

  /* Free cell node ids after VELOCIraptor has copied them. */
  free(cell_node_ids);
  free(swift_parts);

  /* Increase output counter (if not linked with snapshots) */
  if (!e->snapshot_invoke_stf) e->stf_output_count++;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT not configure to run with VELOCIraptor.");
#endif /* HAVE_VELOCIRAPTOR */
}
