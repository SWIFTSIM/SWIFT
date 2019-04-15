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
#ifndef SWIFT_EAGLE_STARS_YIELD_TABLES_H
#define SWIFT_EAGLE_STARS_YIELD_TABLES_H

#include "chemistry.h"

static const float log10_min_metallicity = -20;

/**
 * @brief returns index of element_name within array of element names
 * (element_array)
 *
 * @param element_name name of element
 * @param element_array array of element names
 * @param n_elements size of element_array
 */
inline static int get_element_index(const char *element_name,
                                    char **element_array, int n_elements) {

  /* Compare element name we are trying to index to names in element_array */
  for (int i = 0; i < n_elements; i++) {
    if (strcmp(element_array[i], element_name) == 0) return i;
  }

  /* If we don't find the index return flag  */
  return -1;
}

/**
 * @brief reads yield tables, flattens and stores them in stars_props data
 * struct
 *
 * @param stars the #stars_props data structure
 */
inline static void read_yield_tables(struct feedback_props *feedback_props) {

#ifdef HAVE_HDF5

  int i, j, k, flat_index;

  /* filenames to read HDF5 files */
  char fname[256], setname[100];

  hid_t file_id, dataset, datatype;
  herr_t status;

  /* Open SNIa tables for reading */
  sprintf(fname, "%s/SNIa.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->SNIa_element_names);
  if (status < 0) error("error reading SNIa element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read SNIa yields */
  dataset = H5Dopen(file_id, "Yield", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->yields_SNIa);
  if (status < 0) error("error reading SNIa yields");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read SNIa total metals released */
  dataset = H5Dopen(file_id, "Total_Metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &feedback_props->yield_SNIa_total_metals_IMF_resampled);
  if (status < 0) error("error reading SNIa total metal");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  status = H5Fclose(file_id);
  if (status < 0) error("error closing SNIa file");

  /* Open SNII tables for reading */
  sprintf(fname, "%s/SNII.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->SNII_element_names);
  if (status < 0) error("error reading SNII element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read array of masses */
  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->yield_SNII.mass);
  if (status < 0) error("error reading SNII masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read array of metallicities */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->yield_SNII.metallicity);
  if (status < 0) error("error reading SNII metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* declare temporary arrays to read data from HDF5 files */
  double temp_yield_SNII[feedback_props->SNII_n_elements]
                        [feedback_props->SNII_n_mass];
  double temp_ejecta_SNII[feedback_props->SNII_n_mass],
      tempmet1[feedback_props->SNII_n_mass];
  char *metallicity_yield_table_name_SNII[feedback_props->SNII_n_z];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   metallicity_yield_table_name_SNII);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read SNII yield tables */
  for (i = 0; i < feedback_props->SNII_n_z; i++) {
    /* read yields to temporary array */
    sprintf(setname, "/Yields/%s/Yield", metallicity_yield_table_name_SNII[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_yield_SNII);
    if (status < 0) error("error reading SNII yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read mass ejected table to temporary array */
    sprintf(setname, "/Yields/%s/Ejected_mass",
            metallicity_yield_table_name_SNII[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_ejecta_SNII);
    if (status < 0) error("error reading SNII ejected masses");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read total metals table to temporary array */
    sprintf(setname, "/Yields/%s/Total_Metals",
            metallicity_yield_table_name_SNII[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     tempmet1);
    if (status < 0) error("error reading SNII total metals");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* Flatten the temporary tables that were read, store in stars_props */
    for (k = 0; k < feedback_props->SNII_n_mass; k++) {
      flat_index = row_major_index_2d(i, k, feedback_props->SNII_n_z,
                                               feedback_props->SNII_n_mass);
      feedback_props->yield_SNII.ejecta[flat_index] = temp_ejecta_SNII[k];
      feedback_props->yield_SNII.total_metals[flat_index] = tempmet1[k];

      for (j = 0; j < feedback_props->SNII_n_elements; j++) {
        flat_index = row_major_index_3d(
            i, j, k, feedback_props->SNII_n_z, feedback_props->SNII_n_elements,
            feedback_props->SNII_n_mass);
        feedback_props->yield_SNII.yield[flat_index] = temp_yield_SNII[j][k];
      }
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* Read AGB tables */
  sprintf(fname, "%s/AGB.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->AGB_element_names);
  if (status < 0) error("error reading AGB element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read array of masses */
  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->yield_AGB.mass);
  if (status < 0) error("error reading AGB masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read array of metallicities */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->yield_AGB.metallicity);
  if (status < 0) error("error reading AGB metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* declare temporary arrays to read data from HDF5 files */
  double temp_yield_AGB[feedback_props->AGB_n_elements]
                       [feedback_props->AGB_n_mass];
  double temp_ejecta_AGB[feedback_props->AGB_n_mass],
      tempmet2[feedback_props->AGB_n_mass];
  char *metallicity_yield_table_name_AGB[feedback_props->AGB_n_z];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   metallicity_yield_table_name_AGB);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read AGB yield tables */
  for (i = 0; i < feedback_props->AGB_n_z; i++) {
    /* read yields to temporary array */
    sprintf(setname, "/Yields/%s/Yield", metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_yield_AGB);
    if (status < 0) error("error reading AGB yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read mass ejected table to temporary array */
    sprintf(setname, "/Yields/%s/Ejected_mass",
            metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_ejecta_AGB);
    if (status < 0) error("error reading AGB ejected masses");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read total metals table to temporary array */
    sprintf(setname, "/Yields/%s/Total_Metals",
            metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     tempmet2);
    if (status < 0) error("error reading AGB total metals");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* Flatten the temporary tables that were read, store in stars_props */
    for (k = 0; k < feedback_props->AGB_n_mass; k++) {
      flat_index = row_major_index_2d(i, k, feedback_props->AGB_n_z,
                                               feedback_props->AGB_n_mass);
      feedback_props->yield_AGB.ejecta[flat_index] = temp_ejecta_AGB[k];
      feedback_props->yield_AGB.total_metals[flat_index] = tempmet2[k];

      for (j = 0; j < feedback_props->AGB_n_elements; j++) {
        flat_index = row_major_index_3d(
            i, j, k, feedback_props->AGB_n_z, feedback_props->AGB_n_elements,
            feedback_props->AGB_n_mass);
        feedback_props->yield_AGB.yield[flat_index] = temp_yield_AGB[j][k];
      }
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* open lifetimes table */
  sprintf(fname, "%s/Lifetimes.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read lifetimes mass bins */
  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->lifetimes.mass);
  if (status < 0) error("error reading lifetime table masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read metallicity bins */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->lifetimes.metallicity);
  if (status < 0) error("error reading lifetimes metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* allocate temporary array to read lifetimes */
  double temp_lifetimes[feedback_props->lifetimes.n_z]
                       [feedback_props->lifetimes.n_mass];

  dataset = H5Dopen(file_id, "Lifetimes", H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          temp_lifetimes);
  H5Dclose(dataset);

  for (i = 0; i < feedback_props->lifetimes.n_z; i++) {
    for (j = 0; j < feedback_props->lifetimes.n_mass; j++) {
      feedback_props->lifetimes.dyingtime[i][j] = log10(temp_lifetimes[i][j]);
    }
  }

  H5Fclose(file_id);

#endif
}

/**
 * @brief allocates space for the yield tables
 *
 * @param stars the #stars_props data struct to store the tables in
 */
inline static void allocate_yield_tables(struct feedback_props *feedback_props) {

  /* Allocate array to store SNIa yield tables */
  if (swift_memalign("feedback-tables", (void **)&feedback_props->yields_SNIa,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNIa_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNIa yields array");
  }

  /* Allocate array to store SNIa yield table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNIa_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNIa_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNIa IMF resampled yields array");
  }

  /* Allocate array for AGB mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.mass,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->AGB_n_mass * sizeof(double)) != 0) {
    error("Failed to allocate AGB mass array");
  }

  /* Allocate array for AGB metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->AGB_n_z * sizeof(double)) != 0) {
    error("Failed to allocate AGB metallicity array");
  }

  /* Allocate array to store AGB yield tables */
  if (swift_memalign(
          "feedback-tables", (void **)&feedback_props->yield_AGB.yield,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->AGB_n_z * feedback_props->AGB_n_mass *
              feedback_props->AGB_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate AGB yield array");
  }

  /* Allocate array to store AGB yield table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.yield_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->AGB_n_z * feedback_props->n_imf_mass_bins *
                         chemistry_element_count * sizeof(double)) != 0) {
    error("Failed to allocate AGB IMF resampled array");
  }

  /* Allocate array to store AGB ejecta tables */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.ejecta,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->AGB_n_z * feedback_props->AGB_n_mass *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB ejecta array");
  }

  /* Allocate array to store AGB ejecta table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.ejecta_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->AGB_n_z * feedback_props->n_imf_mass_bins *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB ejecta IMF resampled array");
  }

  /* Allocate array to store table of total metals released by AGB */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.total_metals,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->AGB_n_z * feedback_props->AGB_n_mass *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB total metals array");
  }

  /* Allocate array to store table of total metals released by AGB resampled by
   * IMF mass bins */
  if (swift_memalign(
          "feedback-tables",
          (void **)&feedback_props->yield_AGB.total_metals_IMF_resampled,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->AGB_n_z * feedback_props->n_imf_mass_bins *
              sizeof(double)) != 0) {
    error("Failed to allocate AGB total metals IMF resampled array");
  }

  /* Allocate array for SNII mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.mass,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNII_n_mass * sizeof(double)) != 0) {
    error("Failed to allocate SNII mass array");
  }

  /* Allocate array for SNII metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNII_n_z * sizeof(double)) != 0) {
    error("Failed to allocate SNII metallicity array");
  }

  /* Allocate array to store SNII yield tables */
  if (swift_memalign(
          "feedback-tables", (void **)&feedback_props->yield_SNII.yield,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->SNII_n_z * feedback_props->SNII_n_mass *
              feedback_props->SNII_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNII yield array");
  }

  /* Allocate array to store SNII yield table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.yield_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNII_n_z *
                         feedback_props->n_imf_mass_bins *
                         chemistry_element_count * sizeof(double)) != 0) {
    error("Failed to allocate SNII IMF resampled array");
  }

  /* Allocate array to store SNII ejecta tables */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.ejecta,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNII_n_z * feedback_props->SNII_n_mass *
                         sizeof(double)) != 0) {
    error("Failed to allocate SNII ejecta array");
  }

  /* Allocate array to store SNII ejecta table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.ejecta_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNII_n_z *
                         feedback_props->n_imf_mass_bins * sizeof(double)) !=
      0) {
    error("Failed to allocate SNII ejecta IMF resampled array");
  }

  /* Allocate array to store table of total metals released by SNII */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.total_metals,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->SNII_n_z * feedback_props->SNII_n_mass *
                         sizeof(double)) != 0) {
    error("Failed to allocate SNII total metals array");
  }

  /* Allocate array to store table of total metals released by SNII resampled by
   * IMF mass bins */
  if (swift_memalign(
          "feedback-tables",
          (void **)&feedback_props->yield_SNII.total_metals_IMF_resampled,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->SNII_n_z * feedback_props->n_imf_mass_bins *
              sizeof(double)) != 0) {
    error("Failed to allocate SNII total metals IMF resampled array");
  }

  /* Allocate array for lifetimes mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->lifetimes.mass,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->lifetimes.n_mass * sizeof(double)) != 0) {
    error("Failed to allocate lifetime mass array");
  }

  /* Allocate array for lifetimes metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->lifetimes.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->lifetimes.n_z * sizeof(double)) != 0) {
    error("Failed to allocate lifetime metallicity array");
  }

  /* Allocate lifetimes array */
  feedback_props->lifetimes.dyingtime =
      (double **)malloc(feedback_props->lifetimes.n_z * sizeof(double *));
  for (int i = 0; i < feedback_props->lifetimes.n_z; i++) {
    feedback_props->lifetimes.dyingtime[i] =
        (double *)malloc(feedback_props->lifetimes.n_mass * sizeof(double));
  }

  /* Allocate SNII factor array  */
  if (swift_memalign("feedback-tables", (void **)&feedback_props->typeII_factor,
                     SWIFT_STRUCT_ALIGNMENT,
                     chemistry_element_count * sizeof(double)) != 0) {
    error("Failed to allocate SNII factor array");
  }

  /* Allocate arrays to store names of elements tracked for SNIa, SNII, AGB  */
  feedback_props->SNIa_element_names =
      (char **)malloc(feedback_props->SNIa_n_elements * sizeof(char *));
  for (int i = 0; i < feedback_props->SNIa_n_elements; i++) {
    feedback_props->SNIa_element_names[i] =
        (char *)malloc(feedback_props->element_name_length * sizeof(char));
  }
  feedback_props->SNII_element_names =
      (char **)malloc(feedback_props->SNII_n_elements * sizeof(char *));
  for (int i = 0; i < feedback_props->SNII_n_elements; i++) {
    feedback_props->SNII_element_names[i] =
        (char *)malloc(feedback_props->element_name_length * sizeof(char));
  }
  feedback_props->AGB_element_names =
      (char **)malloc(feedback_props->AGB_n_elements * sizeof(char *));
  for (int i = 0; i < feedback_props->AGB_n_elements; i++) {
    feedback_props->AGB_element_names[i] =
        (char *)malloc(feedback_props->element_name_length * sizeof(char));
  }

  /* Allocate array of IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_mass_bins,
                     SWIFT_STRUCT_ALIGNMENT,
                     feedback_props->n_imf_mass_bins * sizeof(double)) != 0) {
    error("Failed to allocate imf mass bins array");
  }
}

/**
 * @brief resamples yields based on IMF mass bins
 *
 * @param stars the #stars_props data structure
 */
inline static void compute_yields(struct feedback_props *feedback_props) {

  int flat_index_3d, flat_index_2d;

  /* convert SNII tables to log10  */
  for (int i = 0; i < feedback_props->SNII_n_mass; i++) {
    feedback_props->yield_SNII.mass[i] =
        log10(feedback_props->yield_SNII.mass[i]);
  }
  for (int i = 0; i < feedback_props->SNII_n_z; i++) {
    if (feedback_props->yield_SNII.metallicity[i] > 0) {
      feedback_props->yield_SNII.metallicity[i] =
          log10(feedback_props->yield_SNII.metallicity[i]);
    } else {
      feedback_props->yield_SNII.metallicity[i] = log10_min_metallicity;
    }
  }

  /* convert AGB tables to log10  */
  for (int i = 0; i < feedback_props->AGB_n_mass; i++) {
    feedback_props->yield_AGB.mass[i] =
        log10(feedback_props->yield_AGB.mass[i]);
  }
  for (int i = 0; i < feedback_props->AGB_n_z; i++) {
    if (feedback_props->yield_AGB.metallicity[i] > 0) {
      feedback_props->yield_AGB.metallicity[i] =
          log10(feedback_props->yield_AGB.metallicity[i]);
    } else {
      feedback_props->yield_AGB.metallicity[i] = log10_min_metallicity;
    }
  }

  /* Declare temporary tables to accumulate yields */
  double SNII_yield[feedback_props->SNII_n_mass];
  double AGB_yield[feedback_props->AGB_n_mass];
  float result;

  /* Resample yields for each element tracked in EAGLE */
  int element_index = 0;
  for (enum chemistry_element elem = chemistry_element_H;
       elem < chemistry_element_count; elem++) {
    /* SNIa  */
    element_index = get_element_index(chemistry_get_element_name(elem),
                                      feedback_props->SNIa_element_names,
                                      feedback_props->SNIa_n_elements);
    feedback_props->yield_SNIa_IMF_resampled[elem] =
        feedback_props->yields_SNIa[element_index];

    /* SNII  */
    element_index = get_element_index(chemistry_get_element_name(elem),
                                      feedback_props->SNII_element_names,
                                      feedback_props->SNII_n_elements);
    for (int i = 0; i < feedback_props->SNII_n_z; i++) {
      for (int j = 0; j < feedback_props->SNII_n_mass; j++) {
        flat_index_3d = row_major_index_3d(
            i, element_index, j, feedback_props->SNII_n_z,
            feedback_props->SNII_n_elements, feedback_props->SNII_n_mass);
        SNII_yield[j] = feedback_props->yield_SNII.yield[flat_index_3d] *
                        exp(M_LN10 * (-feedback_props->yield_SNII.mass[j]));
      }

      for (int k = 0; k < feedback_props->n_imf_mass_bins; k++) {
        if (feedback_props->yield_mass_bins[k] <
            feedback_props->yield_SNII.mass[0])
          result = SNII_yield[0];
        else if (feedback_props->yield_mass_bins[k] >
                 feedback_props->yield_SNII
                     .mass[feedback_props->SNII_n_mass - 1])
          result = SNII_yield[feedback_props->SNII_n_mass - 1];
        else {
          result = interpolate_1D_non_uniform(
              feedback_props->yield_SNII.mass, SNII_yield,
              feedback_props->SNII_n_mass, feedback_props->yield_mass_bins[k]);
        }

        flat_index_3d = row_major_index_3d(
            i, elem, k, feedback_props->SNII_n_z, chemistry_element_count,
            feedback_props->n_imf_mass_bins);
        feedback_props->yield_SNII.yield_IMF_resampled[flat_index_3d] =
            exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
      }
    }

    for (int i = 0; i < feedback_props->SNII_n_z; i++) {
      for (int k = 0; k < feedback_props->n_imf_mass_bins; k++) {
        flat_index_2d = row_major_index_2d(
            i, k, feedback_props->SNII_n_z, feedback_props->n_imf_mass_bins);
        flat_index_3d = row_major_index_3d(
            i, elem, k, feedback_props->SNII_n_z, chemistry_element_count,
            feedback_props->n_imf_mass_bins);
        if (strcmp(chemistry_get_element_name(elem), "Hydrogen") != 0 ||
            strcmp(chemistry_get_element_name(elem), "Helium") != 0) {
          feedback_props->yield_SNII
              .total_metals_IMF_resampled[flat_index_2d] +=
              (feedback_props->typeII_factor[elem] - 1) *
              feedback_props->yield_SNII.yield_IMF_resampled[flat_index_3d];
        }

        feedback_props->yield_SNII.yield_IMF_resampled[flat_index_3d] *=
            feedback_props->typeII_factor[elem];
      }
    }

    /* AGB  */
    element_index = get_element_index(chemistry_get_element_name(elem),
                                      feedback_props->AGB_element_names,
                                      feedback_props->AGB_n_elements);

    if (element_index < 0) {
      error("element not tracked for AGB");
    } else {
      for (int i = 0; i < feedback_props->AGB_n_z; i++) {
        for (int j = 0; j < feedback_props->AGB_n_mass; j++) {
          flat_index_3d = row_major_index_3d(
              i, element_index, j, feedback_props->AGB_n_z,
              feedback_props->AGB_n_elements, feedback_props->AGB_n_mass);
          AGB_yield[j] = feedback_props->yield_AGB.yield[flat_index_3d] *
                         exp(M_LN10 * (-feedback_props->yield_AGB.mass[j]));
        }

        for (int j = 0; j < feedback_props->n_imf_mass_bins; j++) {
          if (feedback_props->yield_mass_bins[j] <
              feedback_props->yield_AGB.mass[0])
            result = AGB_yield[0];
          else if (feedback_props->yield_mass_bins[j] >
                   feedback_props->yield_AGB
                       .mass[feedback_props->AGB_n_mass - 1])
            result = AGB_yield[feedback_props->AGB_n_mass - 1];
          else
            result = interpolate_1D_non_uniform(
                feedback_props->yield_AGB.mass, AGB_yield,
                feedback_props->AGB_n_mass, feedback_props->yield_mass_bins[j]);

          flat_index_3d = row_major_index_3d(
              i, elem, j, feedback_props->AGB_n_z, chemistry_element_count,
              feedback_props->n_imf_mass_bins);
          feedback_props->yield_AGB.yield_IMF_resampled[flat_index_3d] =
              exp(M_LN10 * feedback_props->yield_mass_bins[j]) * result;
        }
      }
    }
  }
}

/**
 * @brief resamples ejecta based on IMF mass bins
 *
 * @param stars the #stars_props data structure
 */
inline static void compute_ejecta(struct feedback_props *feedback_props) {

  /* Declare temporary tables to accumulate yields */
  double SNII_ejecta[feedback_props->SNII_n_mass];
  double AGB_ejecta[feedback_props->AGB_n_mass];
  float result;

  int flat_index;

  /* Resample SNII ejecta */
  for (int i = 0; i < feedback_props->SNII_n_z; i++) {
    for (int k = 0; k < feedback_props->SNII_n_mass; k++) {
      flat_index = row_major_index_2d(i, k, feedback_props->SNII_n_z,
                                               feedback_props->SNII_n_mass);
      SNII_ejecta[k] = feedback_props->yield_SNII.ejecta[flat_index] *
                       exp(M_LN10 * (-feedback_props->yield_SNII.mass[k]));
    }

    for (int k = 0; k < feedback_props->n_imf_mass_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_SNII.mass[0])
        result = SNII_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_SNII.mass[feedback_props->SNII_n_mass - 1])
        result = SNII_ejecta[feedback_props->SNII_n_mass - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_SNII.mass, SNII_ejecta,
            feedback_props->SNII_n_mass, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, feedback_props->SNII_n_z,
                                               feedback_props->n_imf_mass_bins);
      feedback_props->yield_SNII.ejecta_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }

  /* resample SNII total metals released */
  for (int i = 0; i < feedback_props->SNII_n_z; i++) {
    for (int k = 0; k < feedback_props->SNII_n_mass; k++) {
      flat_index = row_major_index_2d(i, k, feedback_props->SNII_n_z,
                                               feedback_props->SNII_n_mass);
      SNII_ejecta[k] = feedback_props->yield_SNII.total_metals[flat_index] *
                       exp(M_LN10 * (-feedback_props->yield_SNII.mass[k]));
    }

    for (int k = 0; k < feedback_props->n_imf_mass_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_SNII.mass[0])
        result = SNII_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_SNII.mass[feedback_props->SNII_n_mass - 1])
        result = SNII_ejecta[feedback_props->SNII_n_mass - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_SNII.mass, SNII_ejecta,
            feedback_props->SNII_n_mass, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, feedback_props->SNII_n_z,
                                               feedback_props->n_imf_mass_bins);
      feedback_props->yield_SNII.total_metals_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }

  /* AGB yields */
  for (int i = 0; i < feedback_props->AGB_n_z; i++) {
    for (int k = 0; k < feedback_props->AGB_n_mass; k++) {
      flat_index = row_major_index_2d(i, k, feedback_props->AGB_n_z,
                                               feedback_props->AGB_n_mass);
      AGB_ejecta[k] = feedback_props->yield_AGB.ejecta[flat_index] /
                      exp(M_LN10 * feedback_props->yield_AGB.mass[k]);
    }

    for (int k = 0; k < feedback_props->n_imf_mass_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_AGB.mass[0])
        result = AGB_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_AGB.mass[feedback_props->AGB_n_mass - 1])
        result = AGB_ejecta[feedback_props->AGB_n_mass - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_AGB.mass, AGB_ejecta,
            feedback_props->AGB_n_mass, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, feedback_props->AGB_n_z,
                                               feedback_props->n_imf_mass_bins);
      feedback_props->yield_AGB.ejecta_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }

  for (int i = 0; i < feedback_props->AGB_n_z; i++) {
    for (int k = 0; k < feedback_props->AGB_n_mass; k++) {
      flat_index = row_major_index_2d(i, k, feedback_props->AGB_n_z,
                                               feedback_props->AGB_n_mass);
      AGB_ejecta[k] = feedback_props->yield_AGB.total_metals[flat_index] *
                      exp(M_LN10 * (-feedback_props->yield_AGB.mass[k]));
    }

    for (int k = 0; k < feedback_props->n_imf_mass_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_AGB.mass[0])
        result = AGB_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_AGB.mass[feedback_props->AGB_n_mass - 1])
        result = AGB_ejecta[feedback_props->AGB_n_mass - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_AGB.mass, AGB_ejecta,
            feedback_props->AGB_n_mass, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, feedback_props->AGB_n_z,
                                               feedback_props->n_imf_mass_bins);
      feedback_props->yield_AGB.total_metals_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }
}

#endif
