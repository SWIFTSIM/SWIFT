/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#if defined(HAVE_HDF5) && defined(WITH_MPI)

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

/* This object's header. */
#include "distributed_io.h"

/* Local includes. */
#include "black_holes_io.h"
#include "chemistry_io.h"
#include "common_io.h"
#include "cooling_io.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "fof_io.h"
#include "gravity_io.h"
#include "gravity_properties.h"
#include "hydro_io.h"
#include "hydro_properties.h"
#include "io_properties.h"
#include "memuse.h"
#include "output_list.h"
#include "output_options.h"
#include "part.h"
#include "part_type.h"
#include "star_formation_io.h"
#include "stars_io.h"
#include "tools.h"
#include "tracers_io.h"
#include "units.h"
#include "velociraptor_io.h"
#include "xmf.h"

/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param e The #engine we are writing from.
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param partTypeGroupName The name of the group containing the particles in
 * the HDF5 file.
 * @param props The #io_props of the field to read
 * @param N The number of particles to write.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 * the part array will be written once the structures have been stabilized.
 */
void write_distributed_array(const struct engine* e, hid_t grp,
                             const char* fileName,
                             const char* partTypeGroupName,
                             const struct io_props props, const size_t N,
                             const struct unit_system* internal_units,
                             const struct unit_system* snapshot_units) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t num_elements = N * props.dimension;

  /* message("Writing '%s' array...", props.name); */

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy the particle data to the temporary buffer */
  io_copy_temp_buffer(temp, e, props, N, internal_units, snapshot_units);

  /* Create data space */
  hid_t h_space;
  if (N > 0)
    h_space = H5Screate(H5S_SIMPLE);
  else
    h_space = H5Screate(H5S_NULL);

  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  int rank;
  hsize_t shape[2];
  hsize_t chunk_shape[2];

  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = props.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = 0;
  }

  /* Make sure the chunks are not larger than the dataset */
  if (chunk_shape[0] > N) chunk_shape[0] = N;

  /* Change shape of data space */
  hid_t h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
  if (h_err < 0)
    error("Error while changing data space shape for field '%s'.", props.name);

  /* Dataset properties */
  const hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Create filters and set compression level if we have something to write */
  if (N > 0) {

    /* Set chunk size */
    h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
    if (h_err < 0)
      error("Error while setting chunk size (%llu, %llu) for field '%s'.",
            chunk_shape[0], chunk_shape[1], props.name);

    /* Impose check-sum to verify data corruption */
    h_err = H5Pset_fletcher32(h_prop);
    if (h_err < 0)
      error("Error while setting checksum options for field '%s'.", props.name);

    /* Impose data compression */
    if (e->snapshot_compression > 0) {
      h_err = H5Pset_shuffle(h_prop);
      if (h_err < 0)
        error("Error while setting shuffling options for field '%s'.",
              props.name);

      h_err = H5Pset_deflate(h_prop, e->snapshot_compression);
      if (h_err < 0)
        error("Error while setting compression options for field '%s'.",
              props.name);
    }
  }

  /* Create dataset */
  const hid_t h_data = H5Dcreate(grp, props.name, io_hdf5_type(props.type),
                                 h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", props.name);

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(props.type), h_space, H5S_ALL,
                   H5P_DEFAULT, temp);
  if (h_err < 0) error("Error while writing data array '%s'.", props.name);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, props.units,
                              props.scale_factor_exponent);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, props.units);
  io_write_attribute_f(h_data, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(h_data, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(h_data, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(h_data, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(h_data, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(h_data, "h-scale exponent", 0.f);
  io_write_attribute_f(h_data, "a-scale exponent", props.scale_factor_exponent);
  io_write_attribute_s(h_data, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double factor =
      units_cgs_conversion_factor(snapshot_units, props.units);
  io_write_attribute_d(
      h_data,
      "Conversion factor to CGS (not including cosmological corrections)",
      factor);
  io_write_attribute_d(
      h_data,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor * pow(e->cosmology->a, props.scale_factor_exponent));

#ifdef SWIFT_DEBUG_CHECKS
  if (strlen(props.description) == 0)
    error("Invalid (empty) description of the field '%s'", props.name);
#endif

  /* Write the full description */
  io_write_attribute_s(h_data, "Description", props.description);

  /* Free and close everything */
  swift_free("writebuff", temp);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);
}

/**
 * @brief Writes a snapshot distributed into multiple files.
 *
 * @param e The engine containing all the system.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 * @param mpi_rank The rank number of the calling MPI rank.
 * @param mpi_size the number of MPI ranks.
 * @param comm The communicator used by the MPI ranks.
 * @param info The MPI information object.
 *
 * Creates a series of HDF5 output files (1 per MPI node) as a snapshot.
 * Writes the particles contained in the engine.
 * If such files already exist, it is erased and replaced by the new one.
 * The companion XMF file is also updated accordingly.
 */
void write_output_distributed(struct engine* e,
                              const struct unit_system* internal_units,
                              const struct unit_system* snapshot_units,
                              const int mpi_rank, const int mpi_size,
                              MPI_Comm comm, MPI_Info info) {

  hid_t h_file = 0, h_grp = 0;
  int numFiles = mpi_size;
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct spart* sparts = e->s->sparts;
  const struct bpart* bparts = e->s->bparts;
  struct output_options* output_options = e->output_options;
  struct output_list* output_list = e->output_list_snapshots;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int with_fof = e->policy & engine_policy_fof;
  const int with_DM_background = e->s->with_DM_background;
#ifdef HAVE_VELOCIRAPTOR
  const int with_stf = (e->policy & engine_policy_structure_finding) &&
                       (e->s->gpart_group_data != NULL);
#else
  const int with_stf = 0;
#endif

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Nblackholes = e->s->nr_bparts;

  size_t Ndm_background = 0;
  if (with_DM_background) {
    Ndm_background = io_count_dm_background_gparts(gparts, Ntot);
  }

  /* Number of particles that we will write in this file.
   * Recall that background particles are never inhibited and have no extras */
  const size_t Ntot_written =
      e->s->nr_gparts - e->s->nr_inhibited_gparts - e->s->nr_extra_gparts;
  const size_t Ngas_written =
      e->s->nr_parts - e->s->nr_inhibited_parts - e->s->nr_extra_parts;
  const size_t Nstars_written =
      e->s->nr_sparts - e->s->nr_inhibited_sparts - e->s->nr_extra_sparts;
  const size_t Nblackholes_written =
      e->s->nr_bparts - e->s->nr_inhibited_bparts - e->s->nr_extra_bparts;
  const size_t Nbaryons_written =
      Ngas_written + Nstars_written + Nblackholes_written;
  const size_t Ndm_written =
      Ntot_written > 0 ? Ntot_written - Nbaryons_written - Ndm_background : 0;

  int snap_count = -1;
  if (e->snapshot_int_time_label_on)
    snap_count = (int)round(e->time);
  else if (e->snapshot_invoke_stf)
    snap_count = e->stf_output_count;
  else
    snap_count = e->snapshot_output_count;

  int number_digits = -1;
  if (e->snapshot_int_time_label_on)
    number_digits = 6;
  else
    number_digits = 4;

  /* Directory and file name */
  char dirName[1024];
  char fileName[1024];

  if (strnlen(e->snapshot_subdir, PARSER_MAX_LINE_SIZE) > 0) {
    sprintf(dirName, "%s/%s_%0*d", e->snapshot_subdir, e->snapshot_base_name,
            number_digits, snap_count);

    sprintf(fileName, "%s/%s_%0*d/%s_%0*d.%d.hdf5", e->snapshot_subdir,
            e->snapshot_base_name, number_digits, snap_count,
            e->snapshot_base_name, number_digits, snap_count, mpi_rank);

  } else {
    sprintf(dirName, "%s_%0*d", e->snapshot_base_name, number_digits,
            snap_count);

    sprintf(fileName, "%s_%0*d/%s_%0*d.%d.hdf5", e->snapshot_base_name,
            number_digits, snap_count, e->snapshot_base_name, number_digits,
            snap_count, mpi_rank);
  }

  /* Create the directory */
  if (mpi_rank == 0) safe_checkdir(dirName, /*create=*/1);
  MPI_Barrier(comm);

  /* Compute offset in the file and total number of particles */
  const long long N[swift_type_count] = {Ngas_written,   Ndm_written,
                                         Ndm_background, 0,
                                         Nstars_written, Nblackholes_written};

  /* Gather the total number of particles to write */
  long long N_total[swift_type_count] = {0};
  MPI_Allreduce(N, N_total, swift_type_count, MPI_LONG_LONG_INT, MPI_SUM, comm);

  /* Open file */
  /* message("Opening file '%s'.", fileName); */
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fileName);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_TIME);
  const double factor_length =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  /* Determine if we are writing a reduced snapshot, and if so which
   * output selection type to use */
  char current_selection_name[FIELD_BUFFER_SIZE] =
      select_output_header_default_name;
  if (output_list) {
    /* Users could have specified a different Select Output scheme for each
     * snapshot. */
    output_list_get_current_select_output(output_list, current_selection_name);
  }

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
  const int dimension = (int)hydro_dimension;
  io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);
  io_write_attribute(h_grp, "Redshift", DOUBLE, &e->cosmology->z, 1);
  io_write_attribute(h_grp, "Scale-factor", DOUBLE, &e->cosmology->a, 1);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);

  /* Store the time at which the snapshot was written */
  time_t tm = time(NULL);
  struct tm* timeinfo = localtime(&tm);
  char snapshot_date[64];
  strftime(snapshot_date, 64, "%T %F %Z", timeinfo);
  io_write_attribute_s(h_grp, "Snapshot date", snapshot_date);

  /* GADGET-2 legacy values:  Number of particles of each type */
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};

  /* Total number of fields to write per ptype */
  int numFields[swift_type_count] = {0};

  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);

    numFields[ptype] = output_options_get_num_fields_to_write(
        output_options, current_selection_name, ptype);
  }

  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, N, swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total", UINT, numParticles,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total_HighWord", UINT,
                     numParticlesHighWord, swift_type_count);
  double MassTable[swift_type_count] = {0};
  io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                     swift_type_count);
  io_write_attribute_i(h_grp, "NumFilesPerSnapshot", numFiles);
  io_write_attribute_i(h_grp, "ThisFile", mpi_rank);
  io_write_attribute_s(h_grp, "OutputType", "FullVolume");
  io_write_attribute_s(h_grp, "SelectOutput", current_selection_name);

  /* Close header */
  H5Gclose(h_grp);

  /* Write all the meta-data */
  io_write_meta_data(h_file, e, internal_units, snapshot_units);

  /* Now write the top-level cell structure
   * We use a global offset of 0 here. This means that the cells will write
   * their offset with respect to the start of the file they belong to and
   * not a global offset */
  long long global_offsets[swift_type_count] = {0};
  h_grp = H5Gcreate(h_file, "/Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cells group");

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp, e->s->cdim, e->s->dim, e->s->pos_dithering,
                        e->s->cells_top, e->s->nr_cells, e->s->width, mpi_rank,
                        /*distributed=*/1, N_total, global_offsets, numFields,
                        internal_units, snapshot_units);
  H5Gclose(h_grp);

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if there are (a) no particles of this kind, or (b)
     * if we have disabled every field of this particle type. */
    if (numParticles[ptype] == 0 || numFields[ptype] == 0) continue;

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating particle group.\n");

    /* Add an alias name for convenience */
    char aliasName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(aliasName, PARTICLE_GROUP_BUFFER_SIZE, "/%sParticles",
             part_type_names[ptype]);
    hid_t h_err = H5Lcreate_soft(partTypeGroupName, h_grp, aliasName,
                                 H5P_DEFAULT, H5P_DEFAULT);
    if (h_err < 0) error("Error while creating alias for particle group.\n");

    /* Write the number of particles as an attribute */
    io_write_attribute_l(h_grp, "NumberOfParticles", N[ptype]);

    int num_fields = 0;
    struct io_props list[100];
    size_t Nparticles = 0;

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;
    struct spart* sparts_written = NULL;
    struct bpart* bparts_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas: {
        if (Ngas == Ngas_written) {

          /* No inhibted particles: easy case */
          Nparticles = Ngas;
          hydro_write_particles(parts, xparts, list, &num_fields);
          num_fields += chemistry_write_particles(parts, list + num_fields);
          if (with_cooling || with_temperature) {
            num_fields += cooling_write_particles(
                parts, xparts, list + num_fields, e->cooling_func);
          }
          if (with_fof) {
            num_fields += fof_write_parts(parts, xparts, list + num_fields);
          }
          if (with_stf) {
            num_fields +=
                velociraptor_write_parts(parts, xparts, list + num_fields);
          }
          num_fields += tracers_write_particles(
              parts, xparts, list + num_fields, with_cosmology);
          num_fields +=
              star_formation_write_particles(parts, xparts, list + num_fields);

        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Ngas_written;

          /* Allocate temporary arrays */
          if (swift_memalign("parts_written", (void**)&parts_written,
                             part_align,
                             Ngas_written * sizeof(struct part)) != 0)
            error("Error while allocating temporary memory for parts");
          if (swift_memalign("xparts_written", (void**)&xparts_written,
                             xpart_align,
                             Ngas_written * sizeof(struct xpart)) != 0)
            error("Error while allocating temporary memory for xparts");

          /* Collect the particles we want to write */
          io_collect_parts_to_write(parts, xparts, parts_written,
                                    xparts_written, Ngas, Ngas_written);

          /* Select the fields to write */
          hydro_write_particles(parts_written, xparts_written, list,
                                &num_fields);
          num_fields +=
              chemistry_write_particles(parts_written, list + num_fields);
          if (with_cooling || with_temperature) {
            num_fields +=
                cooling_write_particles(parts_written, xparts_written,
                                        list + num_fields, e->cooling_func);
          }
          if (with_fof) {
            num_fields += fof_write_parts(parts_written, xparts_written,
                                          list + num_fields);
          }
          if (with_stf) {
            num_fields += velociraptor_write_parts(
                parts_written, xparts_written, list + num_fields);
          }
          num_fields += tracers_write_particles(
              parts_written, xparts_written, list + num_fields, with_cosmology);
          num_fields += star_formation_write_particles(
              parts_written, xparts_written, list + num_fields);
        }
      } break;

      case swift_type_dark_matter: {
        if (Ntot == Ndm_written) {

          /* This is a DM-only run without background or inhibited particles */
          Nparticles = Ntot;
          darkmatter_write_particles(gparts, list, &num_fields);
          if (with_fof) {
            num_fields += fof_write_gparts(gparts, list + num_fields);
          }
          if (with_stf) {
            num_fields += velociraptor_write_gparts(e->s->gpart_group_data,
                                                    list + num_fields);
          }
        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Ndm_written;

          /* Allocate temporary array */
          if (swift_memalign("gparts_written", (void**)&gparts_written,
                             gpart_align,
                             Ndm_written * sizeof(struct gpart)) != 0)
            error("Error while allocating temporary memory for gparts");

          if (with_stf) {
            if (swift_memalign(
                    "gpart_group_written", (void**)&gpart_group_data_written,
                    gpart_align,
                    Ndm_written * sizeof(struct velociraptor_gpart_data)) != 0)
              error(
                  "Error while allocating temporary memory for gparts STF "
                  "data");
          }

          /* Collect the non-inhibited DM particles from gpart */
          io_collect_gparts_to_write(gparts, e->s->gpart_group_data,
                                     gparts_written, gpart_group_data_written,
                                     Ntot, Ndm_written, with_stf);

          /* Select the fields to write */
          darkmatter_write_particles(gparts_written, list, &num_fields);
          if (with_fof) {
            num_fields += fof_write_gparts(gparts_written, list + num_fields);
          }
          if (with_stf) {
            num_fields += velociraptor_write_gparts(gpart_group_data_written,
                                                    list + num_fields);
          }
        }
      } break;

      case swift_type_dark_matter_background: {

        /* Ok, we need to fish out the particles we want */
        Nparticles = Ndm_background;

        /* Allocate temporary array */
        if (swift_memalign("gparts_written", (void**)&gparts_written,
                           gpart_align,
                           Ndm_background * sizeof(struct gpart)) != 0)
          error("Error while allocating temporart memory for gparts");

        if (with_stf) {
          if (swift_memalign(
                  "gpart_group_written", (void**)&gpart_group_data_written,
                  gpart_align,
                  Ndm_background * sizeof(struct velociraptor_gpart_data)) != 0)
            error(
                "Error while allocating temporart memory for gparts STF "
                "data");
        }

        /* Collect the non-inhibited DM particles from gpart */
        io_collect_gparts_background_to_write(
            gparts, e->s->gpart_group_data, gparts_written,
            gpart_group_data_written, Ntot, Ndm_background, with_stf);

        /* Select the fields to write */
        darkmatter_write_particles(gparts_written, list, &num_fields);
        if (with_fof) {
          num_fields += fof_write_gparts(gparts_written, list + num_fields);
        }
        if (with_stf) {
          num_fields += velociraptor_write_gparts(gpart_group_data_written,
                                                  list + num_fields);
        }
      } break;

      case swift_type_stars: {
        if (Nstars == Nstars_written) {

          /* No inhibted particles: easy case */
          Nparticles = Nstars;
          stars_write_particles(sparts, list, &num_fields, with_cosmology);
          num_fields += chemistry_write_sparticles(sparts, list + num_fields);
          num_fields += tracers_write_sparticles(sparts, list + num_fields,
                                                 with_cosmology);
          if (with_fof) {
            num_fields += fof_write_sparts(sparts, list + num_fields);
          }
          if (with_stf) {
            num_fields += velociraptor_write_sparts(sparts, list + num_fields);
          }
        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nstars_written;

          /* Allocate temporary arrays */
          if (swift_memalign("sparts_written", (void**)&sparts_written,
                             spart_align,
                             Nstars_written * sizeof(struct spart)) != 0)
            error("Error while allocating temporary memory for sparts");

          /* Collect the particles we want to write */
          io_collect_sparts_to_write(sparts, sparts_written, Nstars,
                                     Nstars_written);

          /* Select the fields to write */
          stars_write_particles(sparts_written, list, &num_fields,
                                with_cosmology);
          num_fields +=
              chemistry_write_sparticles(sparts_written, list + num_fields);
          num_fields += tracers_write_sparticles(
              sparts_written, list + num_fields, with_cosmology);
          if (with_fof) {
            num_fields += fof_write_sparts(sparts_written, list + num_fields);
          }
          if (with_stf) {
            num_fields +=
                velociraptor_write_sparts(sparts_written, list + num_fields);
          }
        }
      } break;

      case swift_type_black_hole: {
        if (Nblackholes == Nblackholes_written) {

          /* No inhibted particles: easy case */
          Nparticles = Nblackholes;
          black_holes_write_particles(bparts, list, &num_fields,
                                      with_cosmology);
          num_fields += chemistry_write_bparticles(bparts, list + num_fields);
          if (with_fof) {
            num_fields += fof_write_bparts(bparts, list + num_fields);
          }
          if (with_stf) {
            num_fields += velociraptor_write_bparts(bparts, list + num_fields);
          }
        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nblackholes_written;

          /* Allocate temporary arrays */
          if (swift_memalign("bparts_written", (void**)&bparts_written,
                             bpart_align,
                             Nblackholes_written * sizeof(struct bpart)) != 0)
            error("Error while allocating temporary memory for bparts");

          /* Collect the particles we want to write */
          io_collect_bparts_to_write(bparts, bparts_written, Nblackholes,
                                     Nblackholes_written);

          /* Select the fields to write */
          black_holes_write_particles(bparts_written, list, &num_fields,
                                      with_cosmology);
          num_fields +=
              chemistry_write_bparticles(bparts_written, list + num_fields);
          if (with_fof) {
            num_fields += fof_write_bparts(bparts_written, list + num_fields);
          }
          if (with_stf) {
            num_fields +=
                velociraptor_write_bparts(bparts_written, list + num_fields);
          }
        }
      } break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Did the user specify a non-standard default for the entire particle
     * type? */
    const enum compression_levels compression_level_current_default =
        output_options_get_ptype_default(output_options->select_output,
                                         current_selection_name,
                                         (enum part_type)ptype);

    /* Write everything that is not cancelled */
    int num_fields_written = 0;
    for (int i = 0; i < num_fields; ++i) {

      /* Did the user cancel this field? */
      const int should_write = output_options_should_write_field(
          output_options, current_selection_name, list[i].name,
          (enum part_type)ptype, compression_level_current_default);

      if (should_write) {
        write_distributed_array(e, h_grp, fileName, partTypeGroupName, list[i],
                                Nparticles, internal_units, snapshot_units);
        num_fields_written++;
      }
    }

    /* Only write this now that we know exactly how many fields there are. */
    io_write_attribute_i(h_grp, "NumberOfFields", num_fields_written);

    /* Free temporary arrays */
    if (parts_written) swift_free("parts_written", parts_written);
    if (xparts_written) swift_free("xparts_written", xparts_written);
    if (gparts_written) swift_free("gparts_written", gparts_written);
    if (gpart_group_data_written)
      swift_free("gpart_group_written", gpart_group_data_written);
    if (sparts_written) swift_free("sparts_written", sparts_written);
    if (bparts_written) swift_free("bparts_written", bparts_written);

    /* Close particle group */
    H5Gclose(h_grp);
  }

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

  e->snapshot_output_count++;
  if (e->snapshot_invoke_stf) e->stf_output_count++;
}

#endif /* HAVE_HDF5 && WITH_MPI */
