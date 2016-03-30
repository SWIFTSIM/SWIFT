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
 * @brief Read dark matter particles from HDF5.
 *
 * @param h_grp The HDF5 group in which to read the arrays.
 * @param N The number of particles on that MPI rank.
 * @param N_total The total number of particles (only used in MPI mode)
 * @param offset The offset of the particles for this MPI rank (only used in MPI
 *mode)
 * @param gparts The particle array
 *
 */
__attribute__((always_inline)) INLINE static void darkmatter_read_particles(
    hid_t h_grp, int N, long long N_total, long long offset,
    struct gpart* gparts) {

  /* Read arrays */
  readArray(h_grp, "Coordinates", DOUBLE, N, 3, gparts, N_total, offset, x,
            COMPULSORY);
  readArray(h_grp, "Masses", FLOAT, N, 1, gparts, N_total, offset, mass,
            COMPULSORY);
  readArray(h_grp, "Velocities", FLOAT, N, 3, gparts, N_total, offset, v_full,
            COMPULSORY);
  readArray(h_grp, "ParticleIDs", ULONGLONG, N, 1, gparts, N_total, offset,
            id_or_neg_offset, COMPULSORY);
}

/**
 * @brief Writes the different particles to the HDF5 file
 *
 * @param h_grp The HDF5 group in which to write the arrays.
 * @param fileName The name of the file (unsued in MPI mode).
 * @param partTypeGroupName The name of the group containing the particles in
 *the HDF5 file.
 * @param xmfFile The XMF file to write to (unused in MPI mode).
 * @param Ndm The number of DM particles on that MPI rank.
 * @param Ndm_total The total number of g-particles (only used in MPI mode)
 * @param mpi_rank The MPI rank of this node (only used in MPI mode)
 * @param offset The offset of the particles for this MPI rank (only used in MPI
 *mode)
 * @param gparts The #gpart array
 * @param us The unit system to use
 *
 */
__attribute__((always_inline)) INLINE static void darkmatter_write_particles(
    hid_t h_grp, char* fileName, char* partTypeGroupName, FILE* xmfFile,
    int Ndm, long long Ndm_total, int mpi_rank, long long offset,
    struct gpart* gparts, struct UnitSystem* us) {

  /* Write arrays */
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Coordinates", DOUBLE,
             Ndm, 3, gparts, Ndm_total, mpi_rank, offset, x, us,
             UNIT_CONV_LENGTH);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Masses", FLOAT, Ndm,
             1, gparts, Ndm_total, mpi_rank, offset, mass, us, UNIT_CONV_MASS);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Velocities", FLOAT,
             Ndm, 3, gparts, Ndm_total, mpi_rank, offset, v_full, us,
             UNIT_CONV_SPEED);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "ParticleIDs",
             ULONGLONG, Ndm, 1, gparts, Ndm_total, mpi_rank, offset,
             id_or_neg_offset, us, UNIT_CONV_NO_UNITS);
}
