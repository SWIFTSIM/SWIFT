/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#include "io_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
void hydro_read_particles(struct part* parts, struct io_props* list,
                          int* num_fields) {

  *num_fields = 8;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                parts, mass);
  list[3] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[4] = io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_ENERGY, parts, entropy);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, rho);
}

/**
 * @brief Writes the different particles to the HDF5 file
 *
 * @param h_grp The HDF5 group in which to write the arrays.
 * @param fileName The name of the file (unsued in MPI mode).
 * @param partTypeGroupName The name of the group containing the particles in
 * the HDF5 file.
 * @param xmfFile The XMF file to write to (unused in MPI mode).
 * @param N The number of particles on that MPI rank.
 * @param N_total The total number of particles (only used in MPI mode)
 * @param mpi_rank The MPI rank of this node (only used in MPI mode)
 * @param offset The offset of the particles for this MPI rank (only used in MPI
 * mode)
 * @param parts The particle array
 * @param internal_units The #UnitSystem used internally
 * @param snapshot_units The #UnitSystem used in the snapshots
 *
 */
void hydro_write_particles(struct part* parts, struct io_props* list,
                           int* num_fields) {

  *num_fields = 1;

  /* List what we want to read */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 parts, x);

  /* Write arrays */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Coordinates",
   * DOUBLE, */
  /*            N, 3, parts, N_total, mpi_rank, offset, x, internal_units, */
  /*            snapshot_units, UNIT_CONV_LENGTH); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Velocities",
   * FLOAT, */
  /*            N, 3, parts, N_total, mpi_rank, offset, v, internal_units, */
  /*            snapshot_units, UNIT_CONV_SPEED); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Masses", FLOAT, N,
   * 1, */
  /*            parts, N_total, mpi_rank, offset, mass, internal_units, */
  /*            snapshot_units, UNIT_CONV_MASS); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "SmoothingLength",
   */
  /*            FLOAT, N, 1, parts, N_total, mpi_rank, offset, h,
   * internal_units, */
  /*            snapshot_units, UNIT_CONV_LENGTH); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Entropy", FLOAT,
   * N, */
  /*            1, parts, N_total, mpi_rank, offset, entropy, internal_units, */
  /*            snapshot_units, UNIT_CONV_ENTROPY_PER_UNIT_MASS); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "ParticleIDs", */
  /*            ULONGLONG, N, 1, parts, N_total, mpi_rank, offset, id, */
  /*            internal_units, snapshot_units, UNIT_CONV_NO_UNITS); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Acceleration",
   * FLOAT, */
  /*            N, 3, parts, N_total, mpi_rank, offset, a_hydro, internal_units,
   */
  /*            snapshot_units, UNIT_CONV_ACCELERATION); */
  /* writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Density", FLOAT,
   * N, */
  /*            1, parts, N_total, mpi_rank, offset, rho, internal_units, */
  /*            snapshot_units, UNIT_CONV_DENSITY); */
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {

  /* Kernel function description */
  writeAttribute_s(h_grpsph, "Kernel", kernel_name);

  /* Viscosity and thermal conduction */
  writeAttribute_s(h_grpsph, "Thermal Conductivity Model",
                   "(No treatment) Legacy Gadget-2 as in Springel (2005)");
  writeAttribute_s(h_grpsph, "Viscosity Model",
                   "Legacy Gadget-2 as in Springel (2005)");
  writeAttribute_f(h_grpsph, "Viscosity alpha", const_viscosity_alpha);
  writeAttribute_f(h_grpsph, "Viscosity beta", 3.f);
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
int writeEntropyFlag() { return 1; }
