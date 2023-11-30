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
#include <config.h>

/* Some standard headers. */
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/* This object's header. */
#include "velociraptor_interface.h"

/* Local includes. */
#include "black_holes_io.h"
#include "common_io.h"
#include "cooling.h"
#include "engine.h"
#include "gravity_io.h"
#include "hydro.h"
#include "hydro_io.h"
#include "stars_io.h"
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

  /*! Neutrino density parameter at z = 0 (cosmology.Omega_nu_0) */
  double Omega_nu;

  /*! Neutrino density parameter (cosmology.Omega_k) */
  double Omega_k;

  /*! Baryon density parameter (cosmology.Omega_b) */
  double Omega_b;

  /*! Radiation constant density parameter (cosmology.Omega_lambda) */
  double Omega_Lambda;

  /*! Dark matter density parameter (cosmology.Omega_cdm) */
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
 * @brief Structure for group information to return to swift
 */
struct groupinfo {

  /*! Index of a #gpart in the global array on this MPI rank */
  int index;

  /*! Group number of the #gpart. */
  long long groupID;
};

/**
 * @brief Structure for all information to return from VR invocation
 */
struct vr_return_data {

  /*! Total number of gparts in all groups on this MPI rank */
  int num_gparts_in_groups;

  /*! Assignment of particles to groups (must be freed by Swift, may be NULL) */
  struct groupinfo *group_info;

  /*! Number of most bound particles returned */
  int num_most_bound;

  /*! Swift gpart indexes of most bound particles (must be freed by Swift, may
   * be NULL) */
  int *most_bound_index;
};

int InitVelociraptor(char *config_name, struct unitinfo unit_info,
                     struct siminfo sim_info, const int numthreads);

struct vr_return_data InvokeVelociraptor(
    const int snapnum, char *output_name, struct cosmoinfo cosmo_info,
    struct siminfo sim_info, const size_t num_gravity_parts,
    const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, const int *cell_node_ids,
    const int numthreads, const int return_group_flags,
    const int return_most_bound);

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
  const ptrdiff_t index_offset = gparts - s->gparts;

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

  /* Convert particle properties into VELOCIraptor units.
   * VELOCIraptor wants:
   * - Un-dithered co-moving positions,
   * - Peculiar velocities,
   * - Co-moving potential,
   * - Physical internal energy (for the gas),
   * - Temperatures (for the gas).
   */
  for (int i = 0; i < nr_gparts; i++) {

#ifndef HAVE_VELOCIRAPTOR_WITH_NOMASS
    swift_parts[i].mass = gravity_get_mass(&gparts[i]);
#endif

    swift_parts[i].potential = gravity_get_comoving_potential(&gparts[i]);

    swift_parts[i].type = gparts[i].type;

    swift_parts[i].index = i + index_offset;
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

        convert_part_pos(e, p, xp, swift_parts[i].x);
        convert_part_vel(e, p, xp, swift_parts[i].v);
        swift_parts[i].id = parts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = hydro_get_drifted_physical_internal_energy(p, cosmo);
        swift_parts[i].T = cooling_get_temperature(phys_const, hydro_props, us,
                                                   cosmo, cool_func, p, xp);
      } break;

      case swift_type_stars: {
        const struct spart *sp = &sparts[-gparts[i].id_or_neg_offset];

        convert_spart_pos(e, sp, swift_parts[i].x);
        convert_spart_vel(e, sp, swift_parts[i].v);
        swift_parts[i].id = sparts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = 0.f;
        swift_parts[i].T = 0.f;
      } break;

      case swift_type_black_hole: {
        const struct bpart *bp = &bparts[-gparts[i].id_or_neg_offset];

        convert_bpart_pos(e, bp, swift_parts[i].x);
        convert_bpart_vel(e, bp, swift_parts[i].v);
        swift_parts[i].id = bparts[-gparts[i].id_or_neg_offset].id;
        swift_parts[i].u = 0.f;
        swift_parts[i].T = 0.f;
      } break;

      case swift_type_dark_matter:
      case swift_type_dark_matter_background:
      case swift_type_neutrino:

        convert_gpart_pos(e, &(gparts[i]), swift_parts[i].x);
        convert_gpart_vel(e, &(gparts[i]), swift_parts[i].v);
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

#ifdef HAVE_VELOCIRAPTOR_ORPHANS

/**
 * @brief Write an array to the output HDF5 file
 *
 * @param h_file HDF5 file handle of the file to write to
 * @param name Name of the dataset to write
 * @param buf The data to write out
 * @param dtype_id HDF5 data type to write
 * @param ndims Number of dimensions of the dataset
 * @param dims Total size of the data over all MPI ranks
 * @param start Offset into the dataset for data from this rank
 * @param count Number of particles to write on this rank
 *
 */
void write_orphan_particle_array(hid_t h_file, const char *name,
                                 const void *buf, const hid_t dtype_id,
                                 const int ndims, const hsize_t *dims,
                                 const hsize_t *start, const hsize_t *count,
                                 size_t nr_flagged_all,
                                 const struct unit_system *snapshot_units,
                                 const struct io_props props) {

  /* Make dataset transfer property list */
  hid_t xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
#if defined(HAVE_PARALLEL_HDF5) && defined(WITH_MPI)
  H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  /* Set up dataspaces */
  hid_t file_dspace_id = H5Screate_simple(ndims, dims, NULL);
  if (dims[0] > 0) {
    if (H5Sselect_hyperslab(file_dspace_id, H5S_SELECT_SET, start, NULL, count,
                            NULL) < 0) {
      error(
          "Failed to select hyperslab to write while writing orphan "
          "particles.");
    }
  } else {
    H5Sselect_none(file_dspace_id);
  }
  hid_t mem_dspace_id = H5Screate_simple(ndims, count, NULL);

  /* Create the dataset */
  hid_t dset_id = H5Dcreate(h_file, name, dtype_id, file_dspace_id, H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
  if (dset_id < 0)
    error("Failed to create dataset while writing orphan particles.");

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, props.units,
                              props.scale_factor_exponent);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, props.units);
  io_write_attribute_f(dset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(dset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(dset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(dset_id, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(dset_id, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(dset_id, "h-scale exponent", 0.f);
  io_write_attribute_f(dset_id, "a-scale exponent",
                       props.scale_factor_exponent);
  io_write_attribute_s(dset_id, "Expression for physical CGS units", buffer);

  /* Write data, if there is any */
  if (nr_flagged_all > 0) {
    if (H5Dwrite(dset_id, dtype_id, mem_dspace_id, file_dspace_id,
                 xfer_plist_id, buf) < 0) {
      error("Failed to write dataset while writing orphan particles.");
    }
  }

  /* Tidy up */
  H5Dclose(dset_id);
  H5Sclose(mem_dspace_id);
  H5Sclose(file_dspace_id);
  H5Pclose(xfer_plist_id);
}

/**
 * @brief Write all particles which have ever been most bound to a file
 *
 * @param e The #engine.
 */
void velociraptor_dump_orphan_particles(struct engine *e,
                                        char *outputFileName) {

  const struct space *s = e->s;
  const size_t nr_gparts = s->nr_gparts;
  const struct gpart *gparts = e->s->gparts;

  /* Handle on the other particle types */
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  /* Units */
  const struct unit_system *internal_units = e->internal_units;
  const struct unit_system *snapshot_units = e->snapshot_units;

  /* Count how many particles we need to write out */
  size_t nr_flagged = 0;
  for (size_t i = 0; i < nr_gparts; i += 1) {
    if (gparts[i].has_been_most_bound) nr_flagged += 1;
  }

  /* Allocate write buffers */
  double *pos;
  if (swift_memalign("VR.pos", (void **)&pos, SWIFT_STRUCT_ALIGNMENT,
                     3 * nr_flagged * sizeof(double)) != 0)
    error("Failed to allocate pos buffer for orphan particles.");
  float *vel;
  if (swift_memalign("VR.vel", (void **)&vel, SWIFT_STRUCT_ALIGNMENT,
                     3 * nr_flagged * sizeof(float)) != 0)
    error("Failed to allocate vel buffer for orphan particles.");
  long long *ids;
  if (swift_memalign("VR.ids", (void **)&ids, SWIFT_STRUCT_ALIGNMENT,
                     nr_flagged * sizeof(long long)) != 0)
    error("Failed to allocate ids buffer for orphan particles.");

  /* Populate write buffers */
  for (size_t i = 0, offset = 0; i < nr_gparts; i += 1) {
    if (gparts[i].has_been_most_bound) {
      switch (gparts[i].type) {
        case swift_type_gas: {
          const struct part *p = &parts[-gparts[i].id_or_neg_offset];
          const struct xpart *xp = &xparts[-gparts[i].id_or_neg_offset];
          convert_part_pos(e, p, xp, &pos[3 * offset]);
          convert_part_vel(e, p, xp, &vel[3 * offset]);
          ids[offset] = parts[-gparts[i].id_or_neg_offset].id;
        } break;
        case swift_type_stars: {
          const struct spart *sp = &sparts[-gparts[i].id_or_neg_offset];
          convert_spart_pos(e, sp, &pos[3 * offset]);
          convert_spart_vel(e, sp, &vel[3 * offset]);
          ids[offset] = sparts[-gparts[i].id_or_neg_offset].id;
        } break;
        case swift_type_black_hole: {
          const struct bpart *bp = &bparts[-gparts[i].id_or_neg_offset];
          convert_bpart_pos(e, bp, &pos[3 * offset]);
          convert_bpart_vel(e, bp, &vel[3 * offset]);
          ids[offset] = bparts[-gparts[i].id_or_neg_offset].id;
        } break;
        case swift_type_dark_matter:
        case swift_type_dark_matter_background:
        case swift_type_neutrino: {
          convert_gpart_pos(e, &gparts[i], &pos[3 * offset]);
          convert_gpart_vel(e, &gparts[i], &vel[3 * offset]);
          ids[offset] = gparts[i].id_or_neg_offset;
        } break;
        default:
          error("Particle type not handled by VELOCIraptor.");
      }
      offset += 1;
    }
  }

  /* Determine output file index for this rank:
   * this is the nodeID, unless we're doing collective I/O */
#if defined(HAVE_PARALLEL_HDF5) && defined(WITH_MPI)
  const int filenum = 0;
  const int write_metadata = e->nodeID == 0;
  const int num_files = 1;
#else
  const int filenum = e->nodeID;
  const int write_metadata = 1;
  const int num_files = e->nr_nodes;
#endif

  /* Determine output file name */
  char orphansFileName[FILENAME_BUFFER_SIZE];
  if (num_files > 1) {
    if (snprintf(orphansFileName, FILENAME_BUFFER_SIZE, "%s.orphans.%d.hdf5",
                 outputFileName, filenum) >= FILENAME_BUFFER_SIZE) {
      error(
          "FILENAME_BUFFER_SIZE is too small for orphan particles file name!");
    }
  } else {
    if (snprintf(orphansFileName, FILENAME_BUFFER_SIZE, "%s.orphans.hdf5",
                 outputFileName) >= FILENAME_BUFFER_SIZE) {
      error(
          "FILENAME_BUFFER_SIZE is too small for orphan particles file name!");
    }
  }

  /* Create output file and write metadata */
  hid_t h_file;
  if (write_metadata) {
    h_file =
        H5Fcreate(orphansFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Failed to open file for orphan particles.");
    io_write_meta_data(h_file, e, internal_units, snapshot_units);
    hid_t h_grp =
        H5Gcreate(h_file, "Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    io_write_attribute(h_grp, "NumFilesPerSnapshot", INT, &num_files, 1);
    H5Gclose(h_grp);
    H5Fclose(h_file);
  }

  /* Reopen the output file in MPI mode if necessary */
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
#if defined(HAVE_PARALLEL_HDF5) && defined(WITH_MPI)
  if (H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
    error("Unable to set MPI mode opening file for orphan particles.");
  }
#endif
  h_file = H5Fopen(orphansFileName, H5F_ACC_RDWR, fapl_id);
  if (h_file < 0) error("Failed to open file for orphan particles.");
  H5Pclose(fapl_id);

  /* Determine offsets and lengths to write in output file on this MPI rank */
  long long offset_ll = 0;
  long long count_ll = (long long)nr_flagged;
  long long ntot_ll = nr_flagged;
#if defined(HAVE_PARALLEL_HDF5) && defined(WITH_MPI)
  MPI_Exscan(&count_ll, &offset_ll, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ntot_ll, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  hsize_t start[2] = {(hsize_t)offset_ll, (hsize_t)0};
  hsize_t count[2] = {(hsize_t)count_ll, (hsize_t)3};
  hsize_t dims[2] = {(hsize_t)ntot_ll, (hsize_t)3};
  size_t nr_flagged_all = (size_t)ntot_ll;

  /* Get list of DM output fields - need this to get metadata for pos/vel/ids.
   * Note that this will be wrong if we have non-DM particles and they use
   * different position and velocity units from the DM particles. */
  int num_fields = 0;
  struct io_props list[100];
  darkmatter_write_particles(gparts, list, &num_fields);

  /* Write all particles as PartType1 */
  hid_t h_grp =
      H5Gcreate(h_file, "PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write the particle data */
  for (int i = 0; i < num_fields; i += 1) {
    if (strcmp(list[i].name, "Coordinates") == 0) {
      /* Convert length units if necessary */
      const double factor = units_conversion_factor(
          internal_units, snapshot_units, list[i].units);
      if (factor != 1.0) {
        for (size_t j = 0; j < 3 * nr_gparts; j += 1) pos[j] *= factor;
      }
      /* Write out the coordinates */
      write_orphan_particle_array(h_grp, list[i].name, pos, H5T_NATIVE_DOUBLE,
                                  2, dims, start, count, nr_flagged_all,
                                  snapshot_units, list[i]);
    } else if (strcmp(list[i].name, "Velocities") == 0) {
      /* Convert velocity units if necessary */
      const double factor = units_conversion_factor(
          internal_units, snapshot_units, list[i].units);
      if (factor != 1.0) {
        for (size_t j = 0; j < 3 * nr_gparts; j += 1) vel[j] *= factor;
      }
      /* Write out the velocities */
      write_orphan_particle_array(h_grp, list[i].name, vel, H5T_NATIVE_FLOAT, 2,
                                  dims, start, count, nr_flagged_all,
                                  snapshot_units, list[i]);
    } else if (strcmp(list[i].name, "ParticleIDs") == 0) {
      /* Write out the particle IDs */
      write_orphan_particle_array(h_grp, list[i].name, ids, H5T_NATIVE_LLONG, 1,
                                  dims, start, count, nr_flagged_all,
                                  snapshot_units, list[i]);
    }
  }

  /* Close output file */
  H5Gclose(h_grp);
  H5Fclose(h_file);

  /* Free write buffers */
  swift_free("VR.pos", pos);
  swift_free("VR.vel", vel);
  swift_free("VR.ids", ids);
}
#endif /* HAVE_VELOCIRAPTOR_ORPHANS */

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
  cosmo_info.Omega_m = e->cosmology->Omega_cdm + e->cosmology->Omega_b;
  cosmo_info.Omega_b = e->cosmology->Omega_b;
  cosmo_info.Omega_r = e->cosmology->Omega_r;
  cosmo_info.Omega_k = e->cosmology->Omega_k;
  cosmo_info.Omega_nu = e->cosmology->Omega_nu_0;
  cosmo_info.Omega_Lambda = e->cosmology->Omega_lambda;
  cosmo_info.Omega_cdm = e->cosmology->Omega_cdm;
  cosmo_info.w_de = e->cosmology->w;

  /* Report the cosmo info we use */
  if (e->verbose) {
    message("VELOCIraptor conf: Scale factor: %e", cosmo_info.atime);
    message("VELOCIraptor conf: Little h: %e", cosmo_info.littleh);
    message("VELOCIraptor conf: Omega_m: %e", cosmo_info.Omega_m);
    message("VELOCIraptor conf: Omega_b: %e", cosmo_info.Omega_b);
    message("VELOCIraptor conf: Omega_Lambda: %e", cosmo_info.Omega_Lambda);
    message("VELOCIraptor conf: Omega_cdm: %e", cosmo_info.Omega_cdm);
    message("VELOCIraptor conf: Omega_nu: %e", cosmo_info.Omega_nu);
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

    const double Omega_cdm = e->cosmology->Omega_cdm;
    const double Omega_b = e->cosmology->Omega_b;
    const double Omega_m = Omega_cdm + Omega_b;
    const double critical_density_0 = e->cosmology->critical_density_0;

    /* Linking length based on the mean DM inter-particle separation
     * in the zoom region and assuming the mean density of the Universe
     * is used in the zoom region. */
    double mean_matter_density;
    if (s->with_hydro)
      mean_matter_density = Omega_cdm * critical_density_0;
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
          box_wrap(cells_top[i].loc[0], 0.0, s->dim[0]);
      sim_info.cell_loc[i].loc[1] =
          box_wrap(cells_top[i].loc[1], 0.0, s->dim[1]);
      sim_info.cell_loc[i].loc[2] =
          box_wrap(cells_top[i].loc[2], 0.0, s->dim[2]);
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

#ifdef SWIFT_MEMUSE_REPORTS
  char report_filename[60];
  sprintf(report_filename, "memuse-VR-report-rank%d-step%d.txt", e->nodeID,
          e->step);
  memuse_log_dump(report_filename);
#endif

  /* Determine if we're writing out orphan particles */
#ifdef HAVE_VELOCIRAPTOR_ORPHANS
  const int return_most_bound = 1;
#else
  const int return_most_bound = 0;
#endif

  /* Call VELOCIraptor. */
  struct vr_return_data return_data = InvokeVelociraptor(
      e->stf_output_count, outputFileName, cosmo_info, sim_info, nr_gparts,
      nr_parts, nr_sparts, swift_parts, cell_node_ids, e->nr_threads,
      linked_with_snap, return_most_bound);

  /* Unpack returned data */
  int num_gparts_in_groups = return_data.num_gparts_in_groups;
  struct groupinfo *group_info = return_data.group_info;
  int num_most_bound = return_data.num_most_bound;
  int *most_bound_index = return_data.most_bound_index;

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
  if (return_most_bound && most_bound_index == NULL && num_most_bound > 0) {
    error("VELOCIraptor failed to return most bound particle indexes.");
  }
  if (!return_most_bound && most_bound_index != NULL) {
    error(
        "VELOCIraptor returned most bound particle indexes when it should not "
        "have.");
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

#ifdef HAVE_VELOCIRAPTOR_ORPHANS

  /* Flag most bound particles */
  if (most_bound_index) {
    for (int i = 0; i < num_most_bound; i++) {
      struct gpart *const gp = &(e->s->gparts[most_bound_index[i]]);
      gp->has_been_most_bound = 1;
    }
  }

  /* Output flagged particles (including those flagged in previous invocations)
   */
  if (e->verbose) message("Writing out orphan particles");
  velociraptor_dump_orphan_particles(e, outputFileName);

#endif

  /* Deallocate most bound particle indexes if necessary (may be allocated by
   * VELOCIraptor) */
  if (most_bound_index) free(most_bound_index);

  /* Reset the pthread affinity mask after VELOCIraptor returns. */
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), engine_entry_affinity());

  /* Increase output counter (if not linked with snapshot) */
  if (!linked_with_snap) e->stf_output_count++;

  /* Record we have ran stf this timestep */
  e->stf_this_timestep = 1;

#else
  error("SWIFT not configured to run with VELOCIraptor.");
#endif /* HAVE_VELOCIRAPTOR */
}

/* Dummy VELOCIraptor interface for testing compilation without linking the
 * actual VELOCIraptor library. Uses --enable-dummy-velociraptor configure
 * option. */
#ifdef HAVE_DUMMY_VELOCIRAPTOR

int InitVelociraptor(char *config_name, struct unitinfo unit_info,
                     struct siminfo sim_info, const int numthreads) {
  error("This is only a dummy. Call the real one!");
  return 0;
}

struct vr_return_data InvokeVelociraptor(
    const int snapnum, char *output_name, struct cosmoinfo cosmo_info,
    struct siminfo sim_info, const size_t num_gravity_parts,
    const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, const int *cell_node_ids,
    const int numthreads, const int return_group_flags,
    const int return_most_bound) {
  error("This is only a dummy. Call the real one!");
  struct vr_return_data data = {0};
  return data;
}

#endif /* HAVE_DUMMY_VELOCIRAPTOR */
