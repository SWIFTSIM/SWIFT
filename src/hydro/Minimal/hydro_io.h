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

/**
 * @brief Reads the different particles to the HDF5 file
 *
 * @param h_grp The HDF5 group in which to read the arrays.
 * @param N The number of particles on that MPI rank.
 * @param N_total The total number of particles (only used in MPI mode)
 * @param offset The offset of the particles for this MPI rank (only used in MPI
 *mode)
 * @param parts The particle array
 *
 */
__attribute__((always_inline)) INLINE static void hydro_read_particles(
    hid_t h_grp, int N, long long N_total, long long offset,
    struct part* parts) {

  /* Read arrays */
  readArray(h_grp, "Coordinates", DOUBLE, N, 3, parts, N_total, offset, x,
            COMPULSORY);
  readArray(h_grp, "Velocities", FLOAT, N, 3, parts, N_total, offset, v,
            COMPULSORY);
  readArray(h_grp, "Masses", FLOAT, N, 1, parts, N_total, offset, mass,
            COMPULSORY);
  readArray(h_grp, "SmoothingLength", FLOAT, N, 1, parts, N_total, offset, h,
            COMPULSORY);
  readArray(h_grp, "InternalEnergy", FLOAT, N, 1, parts, N_total, offset, u,
            COMPULSORY);
  readArray(h_grp, "ParticleIDs", ULONGLONG, N, 1, parts, N_total, offset, id,
            COMPULSORY);
  readArray(h_grp, "Acceleration", FLOAT, N, 3, parts, N_total, offset, a_hydro,
            OPTIONAL);
  readArray(h_grp, "Density", FLOAT, N, 1, parts, N_total, offset, rho,
            OPTIONAL);
}

/**
 * @brief Writes the different particles to the HDF5 file
 *
 * @param h_grp The HDF5 group in which to write the arrays.
 * @param fileName The name of the file (unsued in MPI mode).
 * @param partTypeGroupName The name of the group containing the particles in
 *the HDF5 file.
 * @param xmfFile The XMF file to write to (unused in MPI mode).
 * @param N The number of particles on that MPI rank.
 * @param N_total The total number of particles (only used in MPI mode)
 * @param mpi_rank The MPI rank of this node (only used in MPI mode)
 * @param offset The offset of the particles for this MPI rank (only used in MPI
 *mode)
 * @param parts The particle array
 * @param us The unit system to use
 *
 */
__attribute__((always_inline)) INLINE static void hydro_write_particles(
    hid_t h_grp, char* fileName, char* partTypeGroupName, FILE* xmfFile, int N,
    long long N_total, int mpi_rank, long long offset, struct part* parts,
    struct UnitSystem* us) {

  /* Write arrays */
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Coordinates", DOUBLE,
             N, 3, parts, N_total, mpi_rank, offset, x, us, UNIT_CONV_LENGTH);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Velocities", FLOAT,
             N, 3, parts, N_total, mpi_rank, offset, v, us, UNIT_CONV_SPEED);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Masses", FLOAT, N, 1,
             parts, N_total, mpi_rank, offset, mass, us, UNIT_CONV_MASS);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "SmoothingLength",
             FLOAT, N, 1, parts, N_total, mpi_rank, offset, h, us,
             UNIT_CONV_LENGTH);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "InternalEnergy",
             FLOAT, N, 1, parts, N_total, mpi_rank, offset, u, us,
             UNIT_CONV_ENERGY_PER_UNIT_MASS);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "ParticleIDs",
             ULONGLONG, N, 1, parts, N_total, mpi_rank, offset, id, us,
             UNIT_CONV_NO_UNITS);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Acceleration", FLOAT,
             N, 3, parts, N_total, mpi_rank, offset, a_hydro, us,
             UNIT_CONV_ACCELERATION);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Density", FLOAT, N,
             1, parts, N_total, mpi_rank, offset, rho, us, UNIT_CONV_DENSITY);
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {

  /* Kernel function description */
  writeAttribute_s(h_grpsph, "Kernel", kernel_name);
  writeAttribute_f(h_grpsph, "Hydro gamma", const_hydro_gamma);

  /* Viscosity and thermal conduction */
  /* Nothing in this minimal model... */
  writeAttribute_s(h_grpsph, "Thermal Conductivity Model", "No model");
  writeAttribute_s(h_grpsph, "Viscosity Model", "No model");

  /* Time integration properties */
  writeAttribute_f(h_grpsph, "Maximal Delta u change over dt",
                   const_max_u_change);
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
int writeEntropyFlag() { return 0; }
