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


__attribute__((always_inline))
INLINE static void hydro_read_particles(hid_t h_grp, int N, struct part* parts) {

  /* Read arrays */
  readArray(h_grp, "Coordinates", DOUBLE, N, 3, parts, x, COMPULSORY);
  readArray(h_grp, "Velocities", FLOAT, N, 3, parts, v, COMPULSORY);
  readArray(h_grp, "Masses", FLOAT, N, 1, parts, mass, COMPULSORY);
  readArray(h_grp, "SmoothingLength", FLOAT, N, 1, parts, h, COMPULSORY);
  readArray(h_grp, "InternalEnergy", FLOAT, N, 1, parts, u, COMPULSORY); 
  readArray(h_grp, "ParticleIDs", ULONGLONG, N, 1, parts, id, COMPULSORY);
  readArray(h_grp, "Acceleration", FLOAT, N, 3, parts, a, OPTIONAL);
  readArray(h_grp, "Density", FLOAT, N, 1, parts, rho, OPTIONAL);  
}






__attribute__((always_inline))
INLINE static void hydro_write_particles(hid_t h_grp, char* fileName, FILE* xmfFile,
					 int N, struct part* parts,
					 struct UnitSystem* us) {

  /* Write arrays */
  writeArray(h_grp, fileName, xmfFile, "Coordinates", DOUBLE, N, 3, parts, x,
             us, UNIT_CONV_LENGTH);
  writeArray(h_grp, fileName, xmfFile, "Velocities", FLOAT, N, 3, parts, v, us,
             UNIT_CONV_SPEED);
  writeArray(h_grp, fileName, xmfFile, "Masses", FLOAT, N, 1, parts, mass, us,
             UNIT_CONV_MASS);
  writeArray(h_grp, fileName, xmfFile, "SmoothingLength", FLOAT, N, 1, parts, h,
             us, UNIT_CONV_LENGTH);
  writeArray(h_grp, fileName, xmfFile, "InternalEnergy", FLOAT, N, 1, parts,
             u, us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  writeArray(h_grp, fileName, xmfFile, "ParticleIDs", ULONGLONG, N, 1, parts,
             id, us, UNIT_CONV_NO_UNITS);
  writeArray(h_grp, fileName, xmfFile, "Acceleration", FLOAT, N, 3, parts, a,
             us, UNIT_CONV_ACCELERATION);
  writeArray(h_grp, fileName, xmfFile, "Density", FLOAT, N, 1, parts, rho, us,
             UNIT_CONV_DENSITY);

  
}
