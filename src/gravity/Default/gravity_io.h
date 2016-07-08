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
 * mode)
 * @param gparts The particle array
 * @param internal_units The #UnitSystem used internally
 * @param ic_units The #UnitSystem used in the snapshots
 *
 */
__attribute__((always_inline)) INLINE static void darkmatter_read_particles(
    hid_t h_grp, int N, long long N_total, long long offset,
    struct gpart* gparts, const struct UnitSystem* internal_units,
    struct UnitSystem* ic_units) {

  const int num_fields = 4;
  struct io_props list[num_fields];

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, gparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, gparts, v_full);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                gparts, mass);
  list[3] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, gparts, id);

  /* Read arrays */
  /* readArray(h_grp, "Coordinates", DOUBLE, N, 3, gparts, N_total, offset, x,
   */
  /*           COMPULSORY, internal_units, ic_units, UNIT_CONV_LENGTH); */
  /* readArray(h_grp, "Masses", FLOAT, N, 1, gparts, N_total, offset, mass, */
  /*           COMPULSORY, internal_units, ic_units, UNIT_CONV_MASS); */
  /* readArray(h_grp, "Velocities", FLOAT, N, 3, gparts, N_total, offset,
   * v_full, */
  /*           COMPULSORY, internal_units, ic_units, UNIT_CONV_SPEED); */
  /* readArray(h_grp, "ParticleIDs", ULONGLONG, N, 1, gparts, N_total, offset,
   * id, */
  /*           COMPULSORY, internal_units, ic_units, UNIT_CONV_NO_UNITS); */

  /* And read everything */
  for (int i = 0; i < num_fields; ++i)
    readArray(h_grp, list[i], N, N_total, offset, internal_units, ic_units);
}

/**
 * @brief Writes the different particles to the HDF5 file
 *
 * @param h_grp The HDF5 group in which to write the arrays.
 * @param fileName The name of the file (unsued in MPI mode).
 * @param partTypeGroupName The name of the group containing the particles in
 * the HDF5 file.
 * @param xmfFile The XMF file to write to (unused in MPI mode).
 * @param Ndm The number of DM particles on that MPI rank.
 * @param Ndm_total The total number of g-particles (only used in MPI mode)
 * @param mpi_rank The MPI rank of this node (only used in MPI mode)
 * @param offset The offset of the particles for this MPI rank (only used in MPI
 * mode)
 * @param gparts The #gpart array
 * @param internal_units The #UnitSystem used internally
 * @param snapshot_units The #UnitSystem used in the snapshots
 *
 */
__attribute__((always_inline)) INLINE static void darkmatter_write_particles(
    hid_t h_grp, char* fileName, char* partTypeGroupName, FILE* xmfFile,
    int Ndm, long long Ndm_total, int mpi_rank, long long offset,
    struct gpart* gparts, const struct UnitSystem* internal_units,
    const struct UnitSystem* snapshot_units) {

  /* Write arrays */
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Coordinates", DOUBLE,
             Ndm, 3, gparts, Ndm_total, mpi_rank, offset, x, internal_units,
             snapshot_units, UNIT_CONV_LENGTH);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Masses", FLOAT, Ndm,
             1, gparts, Ndm_total, mpi_rank, offset, mass, internal_units,
             snapshot_units, UNIT_CONV_MASS);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "Velocities", FLOAT,
             Ndm, 3, gparts, Ndm_total, mpi_rank, offset, v_full,
             internal_units, snapshot_units, UNIT_CONV_SPEED);
  writeArray(h_grp, fileName, xmfFile, partTypeGroupName, "ParticleIDs",
             ULONGLONG, Ndm, 1, gparts, Ndm_total, mpi_rank, offset, id,
             internal_units, snapshot_units, UNIT_CONV_NO_UNITS);
}
