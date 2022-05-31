/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_FEEDBACK_YIELD_TABLES_H
#define SWIFT_EAGLE_FEEDBACK_YIELD_TABLES_H

/* Local includes. */
#include "chemistry.h"
#include "inline.h"

static const float log10_min_metallicity = -20;

/*! Length of the name fields in the yields tables */
#define eagle_feedback_element_name_length 15

/*! Number of bins used to define the IMF */
#define eagle_feedback_N_imf_bins 200

/*! Number of elements considered for the SNIa yields */
#define eagle_feedback_SNIa_N_elements 42

/*! Number of elements considered for the SNII yields */
#define eagle_feedback_SNII_N_elements 11

/*! Number of mass bins considered for the SNII yields */
#define eagle_feedback_SNII_N_masses 11

/*! Number of metallicity bins considered for the SNII yields */
#define eagle_feedback_SNII_N_metals 5

/*! Number of elements considered for the AGB yields */
#define eagle_feedback_AGB_N_elements 11

/*! Number of mass bins considered for the AGB yields */
#define eagle_feedback_AGB_N_masses 23

/*! Number of metallicity bins considered for the AGB yields */
#define eagle_feedback_AGB_N_metals 3

/*! Number od mass bins along the mass axis of the lifetime table */
#define eagle_feedback_lifetime_N_masses 30

/*! Number od mass bins along the metal axis of the lifetime table */
#define eagle_feedback_lifetime_N_metals 6

/**
 * @brief returns index of element_name within array of element names
 * (element_array)
 *
 * @param element_name name of element
 * @param element_array array of element names
 * @param n_elements size of element_array
 */
INLINE static int get_element_index(const char *element_name,
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
 * @param feedback_props the #feedback_props data struct to read the table into.
 */
INLINE static void read_yield_tables(struct feedback_props *feedback_props) {

#ifdef HAVE_HDF5

  /* filenames to read HDF5 files */
  char fname[256], setname[100];
  char **temp;

  hid_t file_id, dataset, dataset2, datatype, dataspace;
  herr_t status;

  /* Open SNIa tables for reading */
  sprintf(fname, "%s/SNIa.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array into temporary array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  temp = (char **)malloc(eagle_feedback_SNIa_N_elements * sizeof(char *));
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
  if (status < 0) error("error reading SNIa element names");

  /* Copy the element names into their final destination */
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    memcpy(feedback_props->SNIa_element_names[i], temp[i], strlen(temp[i]));
  }

  /* Release the memory allocated by HDF5 for the strings */
  status = H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, temp);
  if (status < 0) error("error freeing string memory");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");
  status = H5Sclose(dataspace);
  if (status < 0) error("error closing dataspace");

  /* Free the temporary memory */
  free(temp);

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

  /**************************************************************************/

  /* Open SNII tables for reading */
  sprintf(fname, "%s/SNII.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array into temporary array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  temp = (char **)malloc(eagle_feedback_SNII_N_elements * sizeof(char *));
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
  if (status < 0) error("error reading SNII element names");

  /* Copy the element names into their final destination */
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    memcpy(feedback_props->SNII_element_names[i], temp[i], strlen(temp[i]));
  }

  /* Release the memory allocated by HDF5 for the strings */
  status = H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, temp);
  if (status < 0) error("error freeing string memory");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");
  status = H5Sclose(dataspace);
  if (status < 0) error("error closing dataspace");

  /* Free the temporary memory */
  free(temp);

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
  double temp_yield_SNII[eagle_feedback_SNII_N_elements]
                        [eagle_feedback_SNII_N_masses];
  double temp_ejecta_SNII[eagle_feedback_SNII_N_masses],
      tempmet1[eagle_feedback_SNII_N_masses];
  char *metallicity_yield_table_name_SNII[eagle_feedback_SNII_N_metals];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset2 = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset2);

  status = H5Dread(dataset2, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   metallicity_yield_table_name_SNII);
  if (status < 0) error("error reading yield table names");

  /* read SNII yield tables */
  for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {

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
    for (int k = 0; k < eagle_feedback_SNII_N_masses; k++) {
      const int flat_index = row_major_index_2d(
          i, k, eagle_feedback_SNII_N_metals, eagle_feedback_SNII_N_masses);

      feedback_props->yield_SNII.ejecta[flat_index] = temp_ejecta_SNII[k];
      feedback_props->yield_SNII.total_metals[flat_index] = tempmet1[k];

      for (int j = 0; j < eagle_feedback_SNII_N_elements; j++) {

        const int flat_index_Z = row_major_index_3d(
            i, j, k, eagle_feedback_SNII_N_metals,
            eagle_feedback_SNII_N_elements, eagle_feedback_SNII_N_masses);

        feedback_props->yield_SNII.yield[flat_index_Z] = temp_yield_SNII[j][k];
      }
    }
  }

  /* Release the memory allocated by HDF5 for the strings */
  status = H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT,
                           metallicity_yield_table_name_SNII);
  if (status < 0) error("error freeing string memory");
  status = H5Dclose(dataset2);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");
  status = H5Sclose(dataspace);
  if (status < 0) error("error closing dataspace");

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /**************************************************************************/

  /* Read AGB tables */
  sprintf(fname, "%s/AGB.hdf5", feedback_props->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  temp = (char **)malloc(eagle_feedback_AGB_N_elements * sizeof(char *));
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
  if (status < 0) error("error reading AGB element names");

  /* Copy the element names into their final destination */
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    memcpy(feedback_props->AGB_element_names[i], temp[i], strlen(temp[i]));
  }

  /* Release the memory allocated by HDF5 for the strings */
  status = H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, temp);
  if (status < 0) error("error freeing string memory");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");
  status = H5Sclose(dataspace);
  if (status < 0) error("error closing dataspace");

  /* Free the temporary memory */
  free(temp);

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
  double temp_yield_AGB[eagle_feedback_AGB_N_elements]
                       [eagle_feedback_AGB_N_masses];
  double temp_ejecta_AGB[eagle_feedback_AGB_N_masses],
      tempmet2[eagle_feedback_AGB_N_masses];
  char *metallicity_yield_table_name_AGB[eagle_feedback_AGB_N_metals];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset2 = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset2);

  status = H5Dread(dataset2, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   metallicity_yield_table_name_AGB);
  if (status < 0) error("error reading yield table names");

  /* read AGB yield tables */
  for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
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
    for (int k = 0; k < eagle_feedback_AGB_N_masses; k++) {

      const int flat_index = row_major_index_2d(
          i, k, eagle_feedback_AGB_N_metals, eagle_feedback_AGB_N_masses);

      feedback_props->yield_AGB.ejecta[flat_index] = temp_ejecta_AGB[k];
      feedback_props->yield_AGB.total_metals[flat_index] = tempmet2[k];

      for (int j = 0; j < eagle_feedback_AGB_N_elements; j++) {
        const int flat_index_Z = row_major_index_3d(
            i, j, k, eagle_feedback_AGB_N_metals, eagle_feedback_AGB_N_elements,
            eagle_feedback_AGB_N_masses);

        feedback_props->yield_AGB.yield[flat_index_Z] = temp_yield_AGB[j][k];
      }
    }
  }

  /* Release the memory allocated by HDF5 for the strings */
  status = H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT,
                           metallicity_yield_table_name_AGB);
  if (status < 0) error("error freeing string memory");
  status = H5Dclose(dataset2);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");
  status = H5Sclose(dataspace);
  if (status < 0) error("error closing dataspace");

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
  double temp_lifetimes[eagle_feedback_lifetime_N_metals]
                       [eagle_feedback_lifetime_N_masses];

  dataset = H5Dopen(file_id, "Lifetimes", H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          temp_lifetimes);
  H5Dclose(dataset);

  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    for (int j = 0; j < eagle_feedback_lifetime_N_masses; j++) {
      feedback_props->lifetimes.dyingtime[i][j] = log10(temp_lifetimes[i][j]);
    }
  }

  H5Fclose(file_id);

#endif
}

/**
 * @brief allocates space for the yield tables
 *
 * @param feedback_props the #feedback_props data struct to store the tables in
 */
INLINE static void allocate_yield_tables(
    struct feedback_props *feedback_props) {

  /* Allocate array to store SNIa yield tables */
  if (swift_memalign("feedback-tables", (void **)&feedback_props->yields_SNIa,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNIa_N_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNIa yields array");
  }

  /* Allocate array to store SNIa yield table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNIa_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNIa_N_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNIa IMF resampled yields array");
  }

  /* Allocate array for AGB mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.mass,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_masses * sizeof(double)) != 0) {
    error("Failed to allocate AGB mass array");
  }

  /* Allocate array for AGB metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * sizeof(double)) != 0) {
    error("Failed to allocate AGB metallicity array");
  }

  /* Allocate array to store AGB yield tables */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.yield,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * eagle_feedback_AGB_N_masses *
                         eagle_feedback_AGB_N_elements * sizeof(double)) != 0) {
    error("Failed to allocate AGB yield array");
  }

  /* Allocate array to store AGB yield table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.yield_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * eagle_feedback_N_imf_bins *
                         chemistry_element_count * sizeof(double)) != 0) {
    error("Failed to allocate AGB IMF resampled array");
  }

  /* Allocate array to store AGB ejecta tables */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.ejecta,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * eagle_feedback_AGB_N_masses *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB ejecta array");
  }

  /* Allocate array to store AGB ejecta table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.ejecta_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * eagle_feedback_N_imf_bins *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB ejecta IMF resampled array");
  }

  /* Allocate array to store table of total metals released by AGB */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_AGB.total_metals,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * eagle_feedback_AGB_N_masses *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB total metals array");
  }

  /* Allocate array to store table of total metals released by AGB resampled by
   * IMF mass bins */
  if (swift_memalign(
          "feedback-tables",
          (void **)&feedback_props->yield_AGB.total_metals_IMF_resampled,
          SWIFT_STRUCT_ALIGNMENT,
          eagle_feedback_AGB_N_metals * eagle_feedback_N_imf_bins *
              sizeof(double)) != 0) {
    error("Failed to allocate AGB total metals IMF resampled array");
  }

  /* Allocate array for SNII mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.mass,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNII_N_masses * sizeof(double)) != 0) {
    error("Failed to allocate SNII mass array");
  }

  /* Allocate array for SNII metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNII_N_metals * sizeof(double)) != 0) {
    error("Failed to allocate SNII metallicity array");
  }

  /* Allocate array to store SNII yield tables */
  if (swift_memalign(
          "feedback-tables", (void **)&feedback_props->yield_SNII.yield,
          SWIFT_STRUCT_ALIGNMENT,
          eagle_feedback_SNII_N_metals * eagle_feedback_SNII_N_masses *
              eagle_feedback_SNII_N_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNII yield array");
  }

  /* Allocate array to store SNII yield table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.yield_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNII_N_metals * eagle_feedback_N_imf_bins *
                         chemistry_element_count * sizeof(double)) != 0) {
    error("Failed to allocate SNII IMF resampled array");
  }

  /* Allocate array to store SNII ejecta tables */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.ejecta,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNII_N_metals *
                         eagle_feedback_SNII_N_masses * sizeof(double)) != 0) {
    error("Failed to allocate SNII ejecta array");
  }

  /* Allocate array to store SNII ejecta table resampled by IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.ejecta_IMF_resampled,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNII_N_metals * eagle_feedback_N_imf_bins *
                         sizeof(double)) != 0) {
    error("Failed to allocate SNII ejecta IMF resampled array");
  }

  /* Allocate array to store table of total metals released by SNII */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_SNII.total_metals,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_SNII_N_metals *
                         eagle_feedback_SNII_N_masses * sizeof(double)) != 0) {
    error("Failed to allocate SNII total metals array");
  }

  /* Allocate array to store table of total metals released by SNII resampled by
   * IMF mass bins */
  if (swift_memalign(
          "feedback-tables",
          (void **)&feedback_props->yield_SNII.total_metals_IMF_resampled,
          SWIFT_STRUCT_ALIGNMENT,
          eagle_feedback_SNII_N_metals * eagle_feedback_N_imf_bins *
              sizeof(double)) != 0) {
    error("Failed to allocate SNII total metals IMF resampled array");
  }

  /* Allocate array for lifetimes mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->lifetimes.mass,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_lifetime_N_masses * sizeof(double)) != 0) {
    error("Failed to allocate lifetime mass array");
  }

  /* Allocate array for lifetimes metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->lifetimes.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_lifetime_N_metals * sizeof(double)) != 0) {
    error("Failed to allocate lifetime metallicity array");
  }

  /* Allocate lifetimes array */
  feedback_props->lifetimes.dyingtime =
      (double **)malloc(eagle_feedback_lifetime_N_metals * sizeof(double *));
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    feedback_props->lifetimes.dyingtime[i] =
        (double *)malloc(eagle_feedback_lifetime_N_masses * sizeof(double));
  }

  /* Allocate arrays to store names of elements tracked for SNIa, SNII, AGB  */
  feedback_props->SNIa_element_names =
      (char **)malloc(eagle_feedback_SNIa_N_elements * sizeof(char *));
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    feedback_props->SNIa_element_names[i] =
        (char *)malloc(eagle_feedback_element_name_length * sizeof(char));
    memset(feedback_props->SNIa_element_names[i], 0,
           eagle_feedback_element_name_length);
  }
  feedback_props->SNII_element_names =
      (char **)malloc(eagle_feedback_SNII_N_elements * sizeof(char *));
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    feedback_props->SNII_element_names[i] =
        (char *)malloc(eagle_feedback_element_name_length * sizeof(char));
    memset(feedback_props->SNII_element_names[i], 0,
           eagle_feedback_element_name_length);
  }
  feedback_props->AGB_element_names =
      (char **)malloc(eagle_feedback_AGB_N_elements * sizeof(char *));
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    feedback_props->AGB_element_names[i] =
        (char *)malloc(eagle_feedback_element_name_length * sizeof(char));
    memset(feedback_props->AGB_element_names[i], 0,
           eagle_feedback_element_name_length);
  }

  /* Allocate array of IMF mass bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->yield_mass_bins,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_N_imf_bins * sizeof(double)) != 0) {
    error("Failed to allocate imf mass bins array");
  }
}

/**
 * @brief resamples yields based on IMF mass bins
 *
 * @param feedback_props the #feedback_props data struct.
 */
INLINE static void compute_yields(struct feedback_props *feedback_props) {

  int flat_index_3d, flat_index_2d;

  /* convert SNII tables to log10  */
  for (int i = 0; i < eagle_feedback_SNII_N_masses; i++) {
    feedback_props->yield_SNII.mass[i] =
        log10(feedback_props->yield_SNII.mass[i]);
  }
  for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
    if (feedback_props->yield_SNII.metallicity[i] > 0) {
      feedback_props->yield_SNII.metallicity[i] =
          log10(feedback_props->yield_SNII.metallicity[i]);
    } else {
      feedback_props->yield_SNII.metallicity[i] = log10_min_metallicity;
    }
  }

  /* convert AGB tables to log10  */
  for (int i = 0; i < eagle_feedback_AGB_N_masses; i++) {
    feedback_props->yield_AGB.mass[i] =
        log10(feedback_props->yield_AGB.mass[i]);
  }
  for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
    if (feedback_props->yield_AGB.metallicity[i] > 0) {
      feedback_props->yield_AGB.metallicity[i] =
          log10(feedback_props->yield_AGB.metallicity[i]);
    } else {
      feedback_props->yield_AGB.metallicity[i] = log10_min_metallicity;
    }
  }

  /* Declare temporary tables to accumulate yields */
  double SNII_yield[eagle_feedback_SNII_N_masses];
  double AGB_yield[eagle_feedback_AGB_N_masses];
  float result;

  /* Resample yields for each element tracked in EAGLE */
  int element_index = 0;
  for (int elem_nr = chemistry_element_H; elem_nr < chemistry_element_count;
       elem_nr++) {

    enum chemistry_element elem = (enum chemistry_element)elem_nr;

    /* SNIa  */
    element_index = get_element_index(chemistry_get_element_name(elem),
                                      feedback_props->SNIa_element_names,
                                      eagle_feedback_SNIa_N_elements);
    feedback_props->yield_SNIa_IMF_resampled[elem] =
        feedback_props->yields_SNIa[element_index];

    /* SNII  */
    element_index = get_element_index(chemistry_get_element_name(elem),
                                      feedback_props->SNII_element_names,
                                      eagle_feedback_SNII_N_elements);
    for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
      for (int j = 0; j < eagle_feedback_SNII_N_masses; j++) {
        flat_index_3d = row_major_index_3d(
            i, element_index, j, eagle_feedback_SNII_N_metals,
            eagle_feedback_SNII_N_elements, eagle_feedback_SNII_N_masses);
        SNII_yield[j] = feedback_props->yield_SNII.yield[flat_index_3d] *
                        exp(M_LN10 * (-feedback_props->yield_SNII.mass[j]));
      }

      for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
        if (feedback_props->yield_mass_bins[k] <
            feedback_props->yield_SNII.mass[0])
          result = SNII_yield[0];
        else if (feedback_props->yield_mass_bins[k] >
                 feedback_props->yield_SNII
                     .mass[eagle_feedback_SNII_N_masses - 1])
          result = SNII_yield[eagle_feedback_SNII_N_masses - 1];
        else {
          result = interpolate_1D_non_uniform(
              feedback_props->yield_SNII.mass, SNII_yield,
              eagle_feedback_SNII_N_masses, feedback_props->yield_mass_bins[k]);
        }

        flat_index_3d = row_major_index_3d(
            i, elem, k, eagle_feedback_SNII_N_metals, chemistry_element_count,
            eagle_feedback_N_imf_bins);
        feedback_props->yield_SNII.yield_IMF_resampled[flat_index_3d] =
            exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
      }
    }

    for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
      for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
        flat_index_2d = row_major_index_2d(i, k, eagle_feedback_SNII_N_metals,
                                           eagle_feedback_N_imf_bins);
        flat_index_3d = row_major_index_3d(
            i, elem, k, eagle_feedback_SNII_N_metals, chemistry_element_count,
            eagle_feedback_N_imf_bins);
        if (strcmp(chemistry_get_element_name(elem), "Hydrogen") != 0 ||
            strcmp(chemistry_get_element_name(elem), "Helium") != 0) {
          feedback_props->yield_SNII
              .total_metals_IMF_resampled[flat_index_2d] +=
              (feedback_props->SNII_yield_factor[elem] - 1) *
              feedback_props->yield_SNII.yield_IMF_resampled[flat_index_3d];
        }

        feedback_props->yield_SNII.yield_IMF_resampled[flat_index_3d] *=
            feedback_props->SNII_yield_factor[elem];
      }
    }

    /* AGB  */
    element_index = get_element_index(chemistry_get_element_name(elem),
                                      feedback_props->AGB_element_names,
                                      eagle_feedback_AGB_N_elements);

    if (element_index < 0) {
      error("element not tracked for AGB");
    } else {
      for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
        for (int j = 0; j < eagle_feedback_AGB_N_masses; j++) {
          flat_index_3d = row_major_index_3d(
              i, element_index, j, eagle_feedback_AGB_N_metals,
              eagle_feedback_AGB_N_elements, eagle_feedback_AGB_N_masses);
          AGB_yield[j] = feedback_props->yield_AGB.yield[flat_index_3d] *
                         exp(M_LN10 * (-feedback_props->yield_AGB.mass[j]));
        }

        for (int j = 0; j < eagle_feedback_N_imf_bins; j++) {
          if (feedback_props->yield_mass_bins[j] <
              feedback_props->yield_AGB.mass[0])
            result = AGB_yield[0];
          else if (feedback_props->yield_mass_bins[j] >
                   feedback_props->yield_AGB
                       .mass[eagle_feedback_AGB_N_masses - 1])
            result = AGB_yield[eagle_feedback_AGB_N_masses - 1];
          else
            result = interpolate_1D_non_uniform(
                feedback_props->yield_AGB.mass, AGB_yield,
                eagle_feedback_AGB_N_masses,
                feedback_props->yield_mass_bins[j]);

          flat_index_3d = row_major_index_3d(
              i, elem, j, eagle_feedback_AGB_N_metals, chemistry_element_count,
              eagle_feedback_N_imf_bins);
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
 * @param feedback_props the #feedback_props data struct.
 */
INLINE static void compute_ejecta(struct feedback_props *feedback_props) {

  /* Declare temporary tables to accumulate yields */
  double SNII_ejecta[eagle_feedback_SNII_N_masses];
  double AGB_ejecta[eagle_feedback_AGB_N_masses];
  float result;

  int flat_index;

  /* Resample SNII ejecta */
  for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
    for (int k = 0; k < eagle_feedback_SNII_N_masses; k++) {
      flat_index = row_major_index_2d(i, k, eagle_feedback_SNII_N_metals,
                                      eagle_feedback_SNII_N_masses);
      SNII_ejecta[k] = feedback_props->yield_SNII.ejecta[flat_index] *
                       exp(M_LN10 * (-feedback_props->yield_SNII.mass[k]));
    }

    for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_SNII.mass[0])
        result = SNII_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_SNII
                   .mass[eagle_feedback_SNII_N_masses - 1])
        result = SNII_ejecta[eagle_feedback_SNII_N_masses - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_SNII.mass, SNII_ejecta,
            eagle_feedback_SNII_N_masses, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, eagle_feedback_SNII_N_metals,
                                      eagle_feedback_N_imf_bins);
      feedback_props->yield_SNII.ejecta_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }

  /* resample SNII total metals released */
  for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
    for (int k = 0; k < eagle_feedback_SNII_N_masses; k++) {
      flat_index = row_major_index_2d(i, k, eagle_feedback_SNII_N_metals,
                                      eagle_feedback_SNII_N_masses);
      SNII_ejecta[k] = feedback_props->yield_SNII.total_metals[flat_index] *
                       exp(M_LN10 * (-feedback_props->yield_SNII.mass[k]));
    }

    for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_SNII.mass[0])
        result = SNII_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_SNII
                   .mass[eagle_feedback_SNII_N_masses - 1])
        result = SNII_ejecta[eagle_feedback_SNII_N_masses - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_SNII.mass, SNII_ejecta,
            eagle_feedback_SNII_N_masses, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, eagle_feedback_SNII_N_metals,
                                      eagle_feedback_N_imf_bins);
      feedback_props->yield_SNII.total_metals_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }

  /* AGB yields */
  for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
    for (int k = 0; k < eagle_feedback_AGB_N_masses; k++) {
      flat_index = row_major_index_2d(i, k, eagle_feedback_AGB_N_metals,
                                      eagle_feedback_AGB_N_masses);
      AGB_ejecta[k] = feedback_props->yield_AGB.ejecta[flat_index] /
                      exp(M_LN10 * feedback_props->yield_AGB.mass[k]);
    }

    for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_AGB.mass[0])
        result = AGB_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_AGB.mass[eagle_feedback_AGB_N_masses - 1])
        result = AGB_ejecta[eagle_feedback_AGB_N_masses - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_AGB.mass, AGB_ejecta,
            eagle_feedback_AGB_N_masses, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, eagle_feedback_AGB_N_metals,
                                      eagle_feedback_N_imf_bins);
      feedback_props->yield_AGB.ejecta_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }

  for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
    for (int k = 0; k < eagle_feedback_AGB_N_masses; k++) {
      flat_index = row_major_index_2d(i, k, eagle_feedback_AGB_N_metals,
                                      eagle_feedback_AGB_N_masses);
      AGB_ejecta[k] = feedback_props->yield_AGB.total_metals[flat_index] *
                      exp(M_LN10 * (-feedback_props->yield_AGB.mass[k]));
    }

    for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
      if (feedback_props->yield_mass_bins[k] <
          feedback_props->yield_AGB.mass[0])
        result = AGB_ejecta[0];
      else if (feedback_props->yield_mass_bins[k] >
               feedback_props->yield_AGB.mass[eagle_feedback_AGB_N_masses - 1])
        result = AGB_ejecta[eagle_feedback_AGB_N_masses - 1];
      else
        result = interpolate_1D_non_uniform(
            feedback_props->yield_AGB.mass, AGB_ejecta,
            eagle_feedback_AGB_N_masses, feedback_props->yield_mass_bins[k]);

      flat_index = row_major_index_2d(i, k, eagle_feedback_AGB_N_metals,
                                      eagle_feedback_N_imf_bins);
      feedback_props->yield_AGB.total_metals_IMF_resampled[flat_index] =
          exp(M_LN10 * feedback_props->yield_mass_bins[k]) * result;
    }
  }
}

#endif /* SWIFT_EAGLE_FEEDBACK_YIELD_TABLES_H */
