/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 *               2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_SPHM1RT_H
#define SWIFT_RT_IO_SPHM1RT_H

#define RT_LABELS_SIZE 10

#include "rt.h"

/**
 * @file src/rt/SPHM1RT/rt_io.h
 * @brief Main header file for no radiative transfer scheme IO routines.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_particles(const struct part* parts,
                                    struct io_props* list) {

  /* List what we want to read */

  /* Note that in the input, we read radiation energy and flux
   * then we convert these quantities to radiation energy per mass and flux per
   * mass in rt_convert_quantities */

  char fieldname[30];
  int count = 0;
  for (int phg = 0; phg < RT_NGROUPS; phg++) {
    sprintf(fieldname, "PhotonEnergiesGroup%d", phg + 1);
    list[count++] =
        io_make_input_field(fieldname, FLOAT, 1, OPTIONAL, UNIT_CONV_ENERGY,
                            parts, rt_data.conserved[phg].urad);
    sprintf(fieldname, "PhotonFluxesGroup%d", phg + 1);
    list[count++] = io_make_input_field(fieldname, FLOAT, 3, OPTIONAL,
                                        UNIT_CONV_ENERGY_VELOCITY, parts,
                                        rt_data.conserved[phg].frad);
  }

  /* Read quantities for thermo-chemistry */
  list[count++] = io_make_input_field(
      "RtElementMassFractions", FLOAT, rt_chemistry_element_count, OPTIONAL,
      UNIT_CONV_NO_UNITS, parts, rt_data.tchem.metal_mass_fraction);

  list[count++] = io_make_input_field(
      "RtSpeciesAbundances", FLOAT, rt_species_count, OPTIONAL,
      UNIT_CONV_NO_UNITS, parts, rt_data.tchem.abundances);

  return count;
}

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_stars(const struct spart* sparts,
                                struct io_props* list) {
  return 0;
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 * we convert radiation energy per mass to radiation energy
 */
INLINE static void rt_convert_conserved_photon_energies(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[g] = part->rt_data.conserved[g].urad * part->mass;
  }
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 * we convert radiation flux per mass to radiation flux
 */
INLINE static void rt_convert_conserved_photon_fluxes(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  int i = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[i++] = part->rt_data.conserved[g].frad[0] * part->mass;
    ret[i++] = part->rt_data.conserved[g].frad[1] * part->mass;
    ret[i++] = part->rt_data.conserved[g].frad[2] * part->mass;
  }
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  /* Note that in the output, we write radiation energy and flux
   * then we convert these quantities from radiation energy per mass and flux
   * per mass
   * */
  int num_elements = 4;

  list[0] = io_make_output_field_convert_part(
      "PhotonEnergies", FLOAT, RT_NGROUPS, UNIT_CONV_ENERGY, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_energies,
      "Photon Energies (all groups)");

  list[1] = io_make_output_field_convert_part(
      "PhotonFluxes", FLOAT, 3 * RT_NGROUPS, UNIT_CONV_ENERGY_VELOCITY, 0,
      parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_fluxes,
      "Photon Fluxes (all groups; x, y, and z coordinates)");

  list[2] = io_make_output_field(
      "RtElementMassFractions", FLOAT, rt_chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, parts, rt_data.tchem.metal_mass_fraction,
      "Fractions of the particles' masses that are in the given element");

  list[3] = io_make_output_field(
      "RtSpeciesAbundances", FLOAT, rt_species_count, UNIT_CONV_NO_UNITS, 0.f,
      parts, rt_data.tchem.abundances,
      "Species Abundances in unit of hydrogen number density");

  return num_elements;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {
  return 0;
}

/**
 * @brief Write the RT model properties to the snapshot.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The engine
 * @param internal_units The internal unit system
 * @param snapshot_units Units used for the snapshot
 * @param rtp The #rt_props
 */
INLINE static void rt_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                    const struct engine* e,
                                    const struct unit_system* internal_units,
                                    const struct unit_system* snapshot_units,
                                    const struct rt_props* rtp) {

#if defined(HAVE_HDF5)

  /* Write scheme name */
  /* ----------------- */
  io_write_attribute_s(h_grp, "RT Scheme", RT_IMPLEMENTATION);

  /* Write photon group counts */
  /* ------------------------- */
  io_write_attribute_i(h_grp, "PhotonGroupNumber", RT_NGROUPS);

  /* Write photon group bin edges */
  /* ---------------------------- */

  /* Note: photon frequency bin edges are kept in cgs. Convert them here to
   * internal units so we're still compatible with swiftsimio. */
  /* TK comment: I think rtp->photon_groups is already in internal unit */
  // const float Hz_internal =
  //    units_cgs_conversion_factor(internal_units, UNIT_CONV_INV_TIME);
  // const float Hz_internal_inv = 1.f / Hz_internal;
  float photon_groups_internal[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++)
    photon_groups_internal[g] = rtp->photon_groups[g];
  // photon_groups_internal[g] = rtp->photon_groups[g] * Hz_internal_inv;

  hid_t type_float = H5Tcopy(io_hdf5_type(FLOAT));

  hsize_t dims[1] = {RT_NGROUPS};
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp, "PhotonGroupEdges", type_float, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type_float, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           photon_groups_internal);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, UNIT_CONV_INV_TIME,
                              /*scale_factor_exponent=*/0);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, UNIT_CONV_INV_TIME);
  io_write_attribute_f(dset, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(dset, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(dset, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(dset, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(dset, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(dset, "h-scale exponent", 0.f);
  io_write_attribute_f(dset, "a-scale exponent", 0.f);
  io_write_attribute_s(dset, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double factor =
      units_cgs_conversion_factor(snapshot_units, UNIT_CONV_INV_TIME);
  io_write_attribute_d(
      dset, "Conversion factor to CGS (not including cosmological corrections)",
      factor);
  io_write_attribute_d(
      dset,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor * pow(e->cosmology->a, 0.f));

  H5Dclose(dset);
  /* H5Tclose(type_float); [> close this later <] */

  /* If without RT, we have nothing more to do. */
  const int with_rt = e->policy & engine_policy_rt;
  if (!with_rt) return;

  /* Write photon group names */
  /* -------------------------*/

  /* Generate Energy Group names */
  char names_energy[RT_NGROUPS * RT_LABELS_SIZE];
  for (int g = 0; g < RT_NGROUPS; g++) {
    char newEname[RT_LABELS_SIZE];
    sprintf(newEname, "Group%d", g + 1);
    strcpy(names_energy + g * RT_LABELS_SIZE, newEname);
  }

  /* Now write them down */
  hid_t type_string_label = H5Tcopy(H5T_C_S1);
  H5Tset_size(type_string_label, RT_LABELS_SIZE);

  hsize_t dimsE[1] = {RT_NGROUPS};
  hid_t spaceE = H5Screate_simple(1, dimsE, NULL);
  hid_t dsetE = H5Dcreate(h_grp_columns, "PhotonEnergies", type_string_label,
                          spaceE, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dsetE, type_string_label, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           names_energy);
  H5Dclose(dsetE);

  /* Generate Fluxes Group Names */
  char names_fluxes[3 * RT_NGROUPS * RT_LABELS_SIZE];
  int i = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    char newFnameX[RT_LABELS_SIZE];
    sprintf(newFnameX, "Group%dX", g + 1);
    strcpy(names_fluxes + i * RT_LABELS_SIZE, newFnameX);
    i++;
    char newFnameY[RT_LABELS_SIZE];
    sprintf(newFnameY, "Group%dY", g + 1);
    strcpy(names_fluxes + i * RT_LABELS_SIZE, newFnameY);
    i++;
    char newFnameZ[RT_LABELS_SIZE];
    sprintf(newFnameZ, "Group%dZ", g + 1);
    strcpy(names_fluxes + i * RT_LABELS_SIZE, newFnameZ);
    i++;
  }

  /* Now write them down */
  hsize_t dimsF[1] = {3 * RT_NGROUPS};
  hid_t spaceF = H5Screate_simple(1, dimsF, NULL);
  hid_t dsetF = H5Dcreate(h_grp_columns, "PhotonFluxes", type_string_label,
                          spaceF, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dsetF, type_string_label, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           names_fluxes);
  H5Dclose(dsetF);

  /* H5Tclose(type_string_label); [> close this later <] */

  /* Write reduced speed of light */
  /* ---------------------------- */
  /* hid_t type2 = H5Tcopy(io_hdf5_type(FLOAT)); */

  hsize_t dims_cred[1] = {1};
  hid_t space_cred = H5Screate_simple(1, dims_cred, NULL);
  hid_t dset_cred =
      H5Dcreate(h_grp, "ReducedLightspeed", type_float, space_cred, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset_cred, type_float, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &rtp->cred_phys);

  /* Write unit conversion factors for this data set */
  char buffer_cred[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer_cred, snapshot_units, UNIT_CONV_VELOCITY,
                              /*scale_factor_exponent=*/0);
  float baseUnitsExp_cred[5];
  units_get_base_unit_exponents_array(baseUnitsExp_cred, UNIT_CONV_VELOCITY);
  io_write_attribute_f(dset_cred, "U_M exponent", baseUnitsExp_cred[UNIT_MASS]);
  io_write_attribute_f(dset_cred, "U_L exponent",
                       baseUnitsExp_cred[UNIT_LENGTH]);
  io_write_attribute_f(dset_cred, "U_t exponent", baseUnitsExp_cred[UNIT_TIME]);
  io_write_attribute_f(dset_cred, "U_I exponent",
                       baseUnitsExp_cred[UNIT_CURRENT]);
  io_write_attribute_f(dset_cred, "U_T exponent",
                       baseUnitsExp_cred[UNIT_TEMPERATURE]);
  io_write_attribute_f(dset_cred, "h-scale exponent", 0.f);
  io_write_attribute_f(dset_cred, "a-scale exponent", 0.f);
  io_write_attribute_s(dset_cred, "Expression for physical CGS units",
                       buffer_cred);

  /* Write the actual number this conversion factor corresponds to */
  /* TODO Mladen: check cosmology. reduced_speed_of_light is physical only for
   * now. */
  const double factor_cred =
      units_cgs_conversion_factor(snapshot_units, UNIT_CONV_VELOCITY);
  io_write_attribute_d(
      dset_cred,
      "Conversion factor to CGS (not including cosmological corrections)",
      factor_cred);
  io_write_attribute_d(
      dset_cred,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor_cred * pow(e->cosmology->a, 0.f));

  H5Dclose(dset_cred);

  /* Clean up after yourself */
  /* ----------------------- */

  /* Close up the types */
  H5Tclose(type_float);
  H5Tclose(type_string_label);

  /* Create an array of element names */
  const int rt_element_name_length = 32;
  char rt_element_names[rt_chemistry_element_count][rt_element_name_length];
  for (int elem = 0; elem < rt_chemistry_element_count; ++elem) {
    sprintf(rt_element_names[elem], "%s",
            rt_chemistry_get_element_name((enum rt_chemistry_element)elem));
  }

  /* Add to the named columns */
  hsize_t rt_dims[1] = {rt_chemistry_element_count};
  hid_t rt_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(rt_type, rt_element_name_length);
  hid_t rt_space = H5Screate_simple(1, rt_dims, NULL);
  hid_t rt_dset = H5Dcreate(h_grp_columns, "RtElementMassFractions", rt_type,
                            rt_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(rt_dset, rt_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           rt_element_names[0]);
  H5Dclose(rt_dset);

  H5Tclose(rt_type);
  H5Sclose(rt_space);

  /* Add the species names to the named columns */
  const int rt_species_name_length = 32;
  char rt_species_names[rt_species_count][rt_species_name_length];
  for (int spec = 0; spec < rt_species_count; ++spec) {
    sprintf(rt_species_names[spec], "%s",
            rt_get_species_name((enum rt_cooling_species)spec));
  }

  /* Add to the named columns */
  hsize_t rts_dims[1] = {rt_species_count};
  hid_t rts_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(rts_type, rt_species_name_length);
  hid_t rts_space = H5Screate_simple(1, rts_dims, NULL);
  hid_t rts_dset = H5Dcreate(h_grp_columns, "RtSpeciesAbundances", rts_type,
                             rts_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(rts_dset, rts_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           rt_species_names[0]);
  H5Dclose(rts_dset);

  H5Tclose(rts_type);
  H5Sclose(rts_space);

#endif /* HAVE_HDF5 */
}

#endif /* SWIFT_RT_IO_SPHM1RT_H */
