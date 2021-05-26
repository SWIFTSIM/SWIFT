/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_GEAR_H
#define SWIFT_RT_IO_GEAR_H

#define RT_LABELS_SIZE 10

/**
 * @file src/rt/GEAR/rt_io.h
 * @brief Main header file for GEAR M1 Closure radiative transfer
 * scheme IO routines.
 */



/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_particles(const struct part* parts, struct io_props* list) {
  return 0;
}

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_stars(const struct spart* sparts, struct io_props* list) {
  return 0;
}




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

  char fieldname[30];
  int count = 0;
  for (int phg = 0; phg < RT_NGROUPS; phg++) {
    sprintf(fieldname, "PhotonEnergiesGroup%d", phg + 1);
    list[count++] =
        io_make_input_field(fieldname, FLOAT, 1, OPTIONAL, UNIT_CONV_ENERGY,
                            parts, rt_data.conserved[phg].energy);
    sprintf(fieldname, "PhotonFluxesGroup%d", phg + 1);
    list[count++] = io_make_input_field(fieldname, FLOAT, 3, OPTIONAL,
                                        UNIT_CONV_RADIATION_FLUX, parts,
                                        rt_data.conserved[phg].flux);
  }

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
 */
INLINE static void rt_convert_conserved_photon_energies(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[g] = part->rt_data.conserved[g].energy;
  }
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 */
INLINE static void rt_convert_conserved_photon_fluxes(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  int i = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[i++] = part->rt_data.conserved[g].flux[0];
    ret[i++] = part->rt_data.conserved[g].flux[1];
    ret[i++] = part->rt_data.conserved[g].flux[2];
  }
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  int num_elements = 2;

  list[0] = io_make_output_field_convert_part(
      "PhotonEnergies", FLOAT, RT_NGROUPS, UNIT_CONV_ENERGY, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_energies,
      "Photon Energies (all groups)");
  list[1] = io_make_output_field_convert_part(
      "PhotonFluxes", FLOAT, 3 * RT_NGROUPS, UNIT_CONV_RADIATION_FLUX, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_fluxes,
      "Photon Fluxes (all groups; x, y, and z coordinates)");

#ifdef SWIFT_RT_DEBUG_CHECKS
  num_elements += 10;
  list[2] =
      io_make_output_field("RTDebugStarIact", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.debug_iact_stars_inject,
                           "number of interactions between this hydro particle"
                           " and any star particle during injection step");
  list[3] =
      io_make_output_field("RTDebugInjectionDone", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.debug_injection_done,
                           "How many times rt_injection_update_photon_density "
                           "has been called");
  list[4] = io_make_output_field(
      "RTDebugCallsIactGradient", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.debug_calls_iact_gradient,
      "number of calls to this particle during the gradient interaction loop");
  list[5] = io_make_output_field("RTDebugCallsIactTransport", INT, 1,
                                 UNIT_CONV_NO_UNITS, 0, parts,
                                 rt_data.debug_calls_iact_transport,
                                 "number of calls to this particle during the"
                                 "transport interaction loop");
  list[6] = io_make_output_field(
      "RTDebugCallsIactGradientInteractions", INT, 1, UNIT_CONV_NO_UNITS, 0,
      parts, rt_data.debug_calls_iact_gradient_interaction,
      "number of calls to this particle during the gradient interaction loop "
      "from the actual interaction function");
  list[7] = io_make_output_field("RTDebugCallsIactTransportInteractions", INT,
                                 1, UNIT_CONV_NO_UNITS, 0, parts,
                                 rt_data.debug_calls_iact_transport_interaction,
                                 "number of calls to this particle during the "
                                 "transport interaction loop from the actual "
                                 "interaction function");
  list[8] =
      io_make_output_field("RTDebugGradientsDone", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.debug_gradients_done,
                           "How many times finalise_gradients was called");
  list[9] =
      io_make_output_field("RTDebugTransportDone", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.debug_transport_done,
                           "How many times finalise_transport was called");
  list[10] = io_make_output_field(
      "RTDebugThermochemistryDone", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.debug_thermochem_done, "How many times rt_tchem was called");
  list[11] = io_make_output_field(
      "RTDebugRadAbsorbedTot", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.debug_radiation_absorbed_tot,
      "Radiation absorbed by this part during its lifetime");
#endif

  return num_elements;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {
  int num_elements = 0;

#ifdef SWIFT_RT_DEBUG_CHECKS
  num_elements += 4;
  list[0] = io_make_output_field("RTDebugHydroIact", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, sparts, rt_data.debug_iact_hydro_inject,
                                 "number of interactions between this hydro "
                                 "particle and any star particle");
  list[1] = io_make_output_field(
      "RTDebugEmissionRateSet", INT, 1, UNIT_CONV_NO_UNITS, 0, sparts,
      rt_data.debug_emission_rate_set, "Stellar photon emission rates set?");
  list[2] = io_make_output_field(
      "RTDebugRadEmittedTot", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0, sparts,
      rt_data.debug_radiation_emitted_tot,
      "Total radiation emitted during the lifetime of this star");
  list[3] = io_make_output_field("RTDebugInjectedPhotonEnergy", FLOAT,
                                 RT_NGROUPS, UNIT_CONV_ENERGY, 0, sparts,
                                 rt_data.debug_injected_energy_tot,
                                 "Total radiation actually injected into gas");
#endif

  return num_elements;
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
  if (rtp->hydro_controlled_injection) {
    io_write_attribute_s(h_grp, "RT Scheme",
                         RT_IMPLEMENTATION ", hydro controlled injection");
  } else {
    io_write_attribute_s(h_grp, "RT Scheme", RT_IMPLEMENTATION);
  }

  /* Write photon group counts */
  io_write_attribute_i(h_grp, "PhotonGroupNumber", RT_NGROUPS);

  /* Write photon group bin edges */
  hid_t type = H5Tcopy(io_hdf5_type(FLOAT));
  H5Tset_size(type, RT_NGROUPS);

  hsize_t dims[1] = {RT_NGROUPS};
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp, "PhotonGroupEdges", type, space, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rtp->photon_groups);

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
  H5Tclose(type);

  /* If without RT, we have nothing more to do. */
  const int with_rt = e->policy & engine_policy_rt;
  if (!with_rt) return;

  /* Write photon group names now */

  /* Generate Energy Group names */
  char names_energy[RT_NGROUPS * RT_LABELS_SIZE];
  for (int g = 0; g < RT_NGROUPS; g++) {
    char newEname[RT_LABELS_SIZE];
    sprintf(newEname, "Group%d", g + 1);
    strcpy(names_energy + g * RT_LABELS_SIZE, newEname);
  }

  /* Now write them down */
  type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, RT_LABELS_SIZE);

  hsize_t dimsE[1] = {RT_NGROUPS};
  hid_t spaceE = H5Screate_simple(1, dimsE, NULL);
  hid_t dsetE = H5Dcreate(h_grp_columns, "PhotonEnergies", type, spaceE,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dsetE, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names_energy);
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
  hid_t dsetF = H5Dcreate(h_grp_columns, "PhotonFluxes", type, spaceF,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dsetF, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names_fluxes);
  H5Dclose(dsetF);

  H5Tclose(type);

#endif /* HAVE_HDF5 */
}

#endif /* SWIFT_RT_IO_GEAR_H */
