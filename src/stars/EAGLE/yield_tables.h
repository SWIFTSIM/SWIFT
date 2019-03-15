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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "chemistry.h"

static const float log_min_metallicity = -20;
static const int n_mass_bins =
    200;  // temporary, put in correct value and move elsewhere.

// Temporary, these functions need to be somewhere else
inline static char *mystrdup(const char *s) {
  char *p;

  p = (char *)malloc(strlen(s) + 1);
  strcpy(p, s);
  return p;
}

/**
 * @brief Returns the 1d index of element with 2d indices i,j
 * from a flattened 2d array in row major order
 *
 * @param i, j Indices of element of interest
 * @param nx, ny Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_2d(int i, int j,
                                                             int nx, int ny) {
  return i * ny + j;
}

/**
 * @brief Returns the 1d index of element with 3d indices i,j,k
 * from a flattened 3d array in row major order
 *
 * @param i, j, k Indices of element of interest
 * @param nx, ny, nz Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_3d(int i, int j,
                                                             int k, int nx,
                                                             int ny, int nz) {
  return i * ny * nz + j * nz + k;
}

/**
 * @brief returns index of element_name within array of element names (element_array)
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
 * @brief reads yield tables, flattens and stores them in stars_props data struct
 *
 * @param stars the #stars_props data structure 
 */
inline static void read_yield_tables(struct stars_props *restrict stars) {
#ifdef HAVE_HDF5

  // ALEXEI: sort out index names
  int i, j, k, index;

  /* filenames to read HDF5 files */
  char fname[256], setname[100];

  hid_t file_id, dataset, datatype;
  herr_t status;

  /* Open SNIa tables for reading */
  sprintf(fname, "%s/SNIa.hdf5", stars->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->SNIa_element_names);
  if (status < 0) error("error reading SNIa element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read SNIa yields */
  dataset = H5Dopen(file_id, "Yield", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->yields_SNIa);
  if (status < 0) error("error reading SNIa yields");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read SNIa total metals released */
  dataset = H5Dopen(file_id, "Total_Metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &stars->yield_SNIa_total_metals_SPH);
  if (status < 0) error("error reading SNIa total metal");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  status = H5Fclose(file_id);
  if (status < 0) error("error closing SNIa file");

  /* Open SNII tables for reading */
  sprintf(fname, "%s/SNII.hdf5", stars->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->SNII_element_names);
  if (status < 0) error("error reading SNII element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read array of masses */
  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->yield_SNII.mass);
  if (status < 0) error("error reading SNII masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read array of metallicities */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->yield_SNII.metallicity);
  if (status < 0) error("error reading SNII metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* declare temporary arrays to read data from HDF5 files */
  double temp_yield_SNII[stars->SNII_n_elements][stars->SNII_n_mass];
  double temp_ejecta_SNII[stars->SNII_n_mass], tempmet1[stars->SNII_n_mass];
  char *metallicity_yield_table_name_SNII[stars->SNII_n_z];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, metallicity_yield_table_name_SNII);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read SNII yield tables */
  for (i = 0; i < stars->SNII_n_z; i++) {
    /* read yields to temporary array */
    sprintf(setname, "/Yields/%s/Yield", metallicity_yield_table_name_SNII[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_yield_SNII);
    if (status < 0) error("error reading SNII yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read mass ejected table to temporary array */
    sprintf(setname, "/Yields/%s/Ejected_mass", metallicity_yield_table_name_SNII[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_ejecta_SNII);
    if (status < 0) error("error reading SNII ejected masses");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read total metals table to temporary array */
    sprintf(setname, "/Yields/%s/Total_Metals", metallicity_yield_table_name_SNII[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     tempmet1);
    if (status < 0) error("error reading SNII total metals");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* Flatten the temporary tables that were read, store in stars_props */
    for (k = 0; k < stars->SNII_n_mass; k++) {
      index = row_major_index_2d(i, k, stars->SNII_n_z, stars->SNII_n_mass);
      stars->yield_SNII.ejecta[index] = temp_ejecta_SNII[k];
      stars->yield_SNII.total_metals[index] = tempmet1[k];

      for (j = 0; j < stars->SNII_n_elements; j++) {
        index = row_major_index_3d(i, j, k, stars->SNII_n_z,
                                   stars->SNII_n_elements, stars->SNII_n_mass);
        stars->yield_SNII.yield[index] = temp_yield_SNII[j][k];
      }
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* Read AGB tables */
  sprintf(fname, "%s/AGB.hdf5", stars->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->AGB_element_names);
  if (status < 0) error("error reading AGB element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read array of masses */
  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->yield_AGB.mass);
  if (status < 0) error("error reading AGB masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read array of metallicities */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->yield_AGB.metallicity);
  if (status < 0) error("error reading AGB metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* declare temporary arrays to read data from HDF5 files */
  double temp_yield_AGB[stars->AGB_n_elements][stars->AGB_n_mass];
  double temp_ejecta_AGB[stars->AGB_n_mass], tempmet2[stars->AGB_n_mass];
  char *metallicity_yield_table_name_AGB[stars->AGB_n_z];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, metallicity_yield_table_name_AGB);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  /* read AGB yield tables */
  for (i = 0; i < stars->AGB_n_z; i++) {
    /* read yields to temporary array */
    sprintf(setname, "/Yields/%s/Yield", metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_yield_AGB);
    if (status < 0) error("error reading AGB yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read mass ejected table to temporary array */
    sprintf(setname, "/Yields/%s/Ejected_mass", metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_ejecta_AGB);
    if (status < 0) error("error reading AGB ejected masses");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* read total metals table to temporary array */
    sprintf(setname, "/Yields/%s/Total_Metals", metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     tempmet2);
    if (status < 0) error("error reading AGB total metals");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    /* Flatten the temporary tables that were read, store in stars_props */
    for (k = 0; k < stars->AGB_n_mass; k++) {
      index = row_major_index_2d(i, k, stars->AGB_n_z, stars->AGB_n_mass);
      stars->yield_AGB.ejecta[index] = temp_ejecta_AGB[k];
      stars->yield_AGB.total_metals[index] = tempmet2[k];

      for (j = 0; j < stars->AGB_n_elements; j++) {
        index = row_major_index_3d(i, j, k, stars->AGB_n_z,
                                   stars->AGB_n_elements, stars->AGB_n_mass);
        stars->yield_AGB.yield[index] = temp_yield_AGB[j][k];
      }
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* open lifetimes table */
  sprintf(fname, "%s/Lifetimes.hdf5", stars->yield_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read lifetimes mass bins */
  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->lifetimes.mass);
  if (status < 0) error("error reading lifetime table masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read metallicity bins */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   stars->lifetimes.metallicity);
  if (status < 0) error("error reading lifetimes metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* allocate temporary array to read lifetimes */
  double temp_lifetimes[stars->lifetimes.n_z][stars->lifetimes.n_mass];

  dataset = H5Dopen(file_id, "Lifetimes", H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_lifetimes);
  H5Dclose(dataset);

  for (i = 0; i < stars->lifetimes.n_z; i++) {
    for (j = 0; j < stars->lifetimes.n_mass; j++) {
      stars->lifetimes.dyingtime[i][j] = log10(temp_lifetimes[i][j]);
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
inline static void allocate_yield_tables(struct stars_props *restrict stars) {

  /* Allocate array to store SNIa yield tables */
  if (posix_memalign((void **)&stars->yields_SNIa, SWIFT_STRUCT_ALIGNMENT,
                     stars->SNIa_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNIa yields array");
  }

  /* Allocate array to store SNIa yield table resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_SNIa_SPH, SWIFT_STRUCT_ALIGNMENT,
                     stars->SNIa_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNIa SPH yields array");
  }

  /* Allocate array for AGB mass bins */
  if (posix_memalign((void **)&stars->yield_AGB.mass, SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_mass * sizeof(double)) != 0) {
    error("Failed to allocate AGB mass array");
  }
  
  /* Allocate array for AGB metallicity bins */
  if (posix_memalign((void **)&stars->yield_AGB.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_z * sizeof(double)) != 0) {
    error("Failed to allocate AGB metallicity array");
  }
  
  /* Allocate array to store AGB yield tables */
  if (posix_memalign((void **)&stars->yield_AGB.yield, SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_z * stars->AGB_n_mass *
                         stars->AGB_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate AGB yield array");
  }

  /* Allocate array to store AGB yield table resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_AGB.SPH, SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_z * n_mass_bins * chemistry_element_count *
                         sizeof(double)) != 0) {
    error("Failed to allocate AGB SPH array");
  }

  /* Allocate array to store AGB ejecta tables */
  if (posix_memalign((void **)&stars->yield_AGB.ejecta, SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_z * stars->AGB_n_mass * sizeof(double)) !=
      0) {
    error("Failed to allocate AGB ejecta array");
  }

  /* Allocate array to store AGB ejecta table resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_AGB.ejecta_SPH,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_z * n_mass_bins * sizeof(double)) != 0) {
    error("Failed to allocate AGB ejecta SPH array");
  }
  
  /* Allocate array to store table of total metals released by AGB */
  if (posix_memalign(
          (void **)&stars->yield_AGB.total_metals, SWIFT_STRUCT_ALIGNMENT,
          stars->AGB_n_z * stars->AGB_n_mass * sizeof(double)) != 0) {
    error("Failed to allocate AGB total metals array");
  }
  
  /* Allocate array to store table of total metals released by AGB resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_AGB.total_metals_SPH,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->AGB_n_z * n_mass_bins * sizeof(double)) != 0) {
    error("Failed to allocate AGB total metals SPH array");
  }

  /* Allocate array for SNII mass bins */
  if (posix_memalign((void **)&stars->yield_SNII.mass, SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_mass * sizeof(double)) != 0) {
    error("Failed to allocate SNII mass array");
  }

  /* Allocate array for SNII metallicity bins */
  if (posix_memalign((void **)&stars->yield_SNII.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_z * sizeof(double)) != 0) {
    error("Failed to allocate SNII metallicity array");
  }
  
  /* Allocate array to store SNII yield tables */
  if (posix_memalign((void **)&stars->yield_SNII.yield, SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_z * stars->SNII_n_mass *
                         stars->SNII_n_elements * sizeof(double)) != 0) {
    error("Failed to allocate SNII yield array");
  }

  /* Allocate array to store SNII yield table resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_SNII.SPH, SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_z * n_mass_bins * chemistry_element_count *
                         sizeof(double)) != 0) {
    error("Failed to allocate SNII SPH array");
  }

  /* Allocate array to store SNII ejecta tables */
  if (posix_memalign((void **)&stars->yield_SNII.ejecta, SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_z * stars->SNII_n_mass * sizeof(double)) !=
      0) {
    error("Failed to allocate SNII ejecta array");
  }

  /* Allocate array to store SNII ejecta table resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_SNII.ejecta_SPH,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_z * n_mass_bins * sizeof(double)) != 0) {
    error("Failed to allocate SNII ejecta SPH array");
  }

  /* Allocate array to store table of total metals released by SNII */
  if (posix_memalign(
          (void **)&stars->yield_SNII.total_metals, SWIFT_STRUCT_ALIGNMENT,
          stars->SNII_n_z * stars->SNII_n_mass * sizeof(double)) != 0) {
    error("Failed to allocate SNII total metals array");
  }
  
  /* Allocate array to store table of total metals released by SNII resampled by IMF mass bins */
  if (posix_memalign((void **)&stars->yield_SNII.total_metals_SPH,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->SNII_n_z * n_mass_bins * sizeof(double)) != 0) {
    error("Failed to allocate SNII total metals SPH array");
  }

  /* Allocate array for lifetimes mass bins */
  if (posix_memalign((void **)&stars->lifetimes.mass, SWIFT_STRUCT_ALIGNMENT,
                     stars->lifetimes.n_mass * sizeof(double)) != 0) {
    error("Failed to allocate lifetime mass array");
  }

  /* Allocate array for lifetimes metallicity bins */
  if (posix_memalign((void **)&stars->lifetimes.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     stars->lifetimes.n_z * sizeof(double)) != 0) {
    error("Failed to allocate lifetime metallicity array");
  }
  
  /* Allocate lifetimes array */
  stars->lifetimes.dyingtime =
      (double **)malloc(stars->lifetimes.n_z * sizeof(double *));
  for (int i = 0; i < stars->lifetimes.n_z; i++) {
    stars->lifetimes.dyingtime[i] =
        (double *)malloc(stars->lifetimes.n_mass * sizeof(double));
  }

  /* Allocate SNII factor array  */
  if (posix_memalign((void **)&stars->typeII_factor, SWIFT_STRUCT_ALIGNMENT,
                     chemistry_element_count * sizeof(double)) != 0) {
    error("Failed to allocate SNII factor array");
  }

  /* Allocate arrays to store names of elements tracked for SNIa, SNII, AGB  */
  stars->SNIa_element_names =
      (char **)malloc(stars->SNIa_n_elements * sizeof(char *));
  for (int i = 0; i < stars->SNIa_n_elements; i++) {
    stars->SNIa_element_names[i] =
        (char *)malloc(stars->element_name_length * sizeof(char));
  }
  stars->SNII_element_names =
      (char **)malloc(stars->SNII_n_elements * sizeof(char *));
  for (int i = 0; i < stars->SNII_n_elements; i++) {
    stars->SNII_element_names[i] =
        (char *)malloc(stars->element_name_length * sizeof(char));
  }
  stars->AGB_element_names =
      (char **)malloc(stars->AGB_n_elements * sizeof(char *));
  for (int i = 0; i < stars->AGB_n_elements; i++) {
    stars->AGB_element_names[i] =
        (char *)malloc(stars->element_name_length * sizeof(char));
  }

  /* Allocate array of IMF mass bins */
  if (posix_memalign((void **)&stars->yield_mass_bins, SWIFT_STRUCT_ALIGNMENT,
                     n_mass_bins * sizeof(double)) != 0) {
    error("Failed to allocate imf mass bins array");
  }
}

/**
 * @brief resamples yields based on IMF mass bins
 *
 * @param stars the #stars_props data structure 
 */
inline static void compute_yields(struct stars_props *restrict stars) {

  int index, index_2d;

  /* convert SNII tables to log10  */
  for (int i = 0; i < stars->SNII_n_mass; i++) {
    stars->yield_SNII.mass[i] = log10(stars->yield_SNII.mass[i]);
  }
  for (int i = 0; i < stars->SNII_n_z; i++) {
    if (stars->yield_SNII.metallicity[i] > 0) {
      stars->yield_SNII.metallicity[i] =
          log10(stars->yield_SNII.metallicity[i]);
    } else {
      stars->yield_SNII.metallicity[i] = log_min_metallicity;
    }
  }

  /* convert AGB tables to log10  */
  for (int i = 0; i < stars->AGB_n_mass; i++) {
    stars->yield_AGB.mass[i] = log10(stars->yield_AGB.mass[i]);
  }
  for (int i = 0; i < stars->AGB_n_z; i++) {
    if (stars->yield_AGB.metallicity[i] > 0) {
      stars->yield_AGB.metallicity[i] = log10(stars->yield_AGB.metallicity[i]);
    } else {
      stars->yield_AGB.metallicity[i] = log_min_metallicity;
    }
  }

  /* Declare temporary table to accumulate yields */
  double SNII_yield[stars->SNII_n_mass];
  double AGB_yield[stars->AGB_n_mass];
  float result;

  /* Resample yields for each element tracked in EAGLE */
  int element_index = 0;
  for (enum chemistry_element elem = chemistry_element_H;
       elem < chemistry_element_count; elem++) {
    /* SNIa  */
    element_index =
        get_element_index(chemistry_get_element_name(elem),
                          stars->SNIa_element_names, stars->SNIa_n_elements);
    stars->yield_SNIa_SPH[elem] = stars->yields_SNIa[element_index];

    /* SNII  */
    element_index =
        get_element_index(chemistry_get_element_name(elem),
                          stars->SNII_element_names, stars->SNII_n_elements);
    for (int i = 0; i < stars->SNII_n_z; i++) {
      for (int j = 0; j < stars->SNII_n_mass; j++) {
        index = row_major_index_3d(i, element_index, j, stars->SNII_n_z,
                                   stars->SNII_n_elements, stars->SNII_n_mass);
        SNII_yield[j] = stars->yield_SNII.yield[index] *
                        exp(M_LN10 * (-stars->yield_SNII.mass[j]));
      }

      for (int k = 0; k < n_mass_bins; k++) {
        if (stars->yield_mass_bins[k] < stars->yield_SNII.mass[0])
          result = SNII_yield[0];
        else if (stars->yield_mass_bins[k] >
                 stars->yield_SNII.mass[stars->SNII_n_mass - 1])
          result = SNII_yield[stars->SNII_n_mass - 1];
        else {
          result =
              interpolate_1D_non_uniform(stars->yield_SNII.mass, SNII_yield,
                             stars->SNII_n_mass, stars->yield_mass_bins[k]);
        }

        index = row_major_index_3d(i, elem, k, stars->SNII_n_z,
                                   chemistry_element_count, n_mass_bins);
        stars->yield_SNII.SPH[index] =
            exp(M_LN10 * stars->yield_mass_bins[k]) * result;
      }
    }

    for (int i = 0; i < stars->SNII_n_z; i++) {
      for (int k = 0; k < n_mass_bins; k++) {
        index_2d = row_major_index_2d(i, k, stars->SNII_n_z, n_mass_bins);
        index = row_major_index_3d(i, elem, k, stars->SNII_n_z,
                                   chemistry_element_count, n_mass_bins);
        if (strcmp(chemistry_get_element_name(elem), "Hydrogen") != 0 ||
            strcmp(chemistry_get_element_name(elem), "Helium") != 0) {
          stars->yield_SNII.total_metals_SPH[index_2d] +=
              (stars->typeII_factor[elem] - 1) *
              stars->yield_SNII.SPH[index];
        }

        stars->yield_SNII.SPH[index] *= stars->typeII_factor[elem];
      }
    }

    /* AGB  */
    element_index =
        get_element_index(chemistry_get_element_name(elem),
                          stars->AGB_element_names, stars->AGB_n_elements);

    if (element_index < 0) {
      error("element not tracked for AGB");
    } else {
      for (int i = 0; i < stars->AGB_n_z; i++) {
        for (int j = 0; j < stars->AGB_n_mass; j++) {
          index = row_major_index_3d(i, element_index, j, stars->AGB_n_z,
                                     stars->AGB_n_elements, stars->AGB_n_mass);
          AGB_yield[j] = stars->yield_AGB.yield[index] *
                         exp(M_LN10 * (-stars->yield_AGB.mass[j]));
        }

        for (int j = 0; j < n_mass_bins; j++) {
          if (stars->yield_mass_bins[j] < stars->yield_AGB.mass[0])
            result = AGB_yield[0];
          else if (stars->yield_mass_bins[j] >
                   stars->yield_AGB.mass[stars->AGB_n_mass - 1])
            result = AGB_yield[stars->AGB_n_mass - 1];
          else
            result =
                interpolate_1D_non_uniform(stars->yield_AGB.mass, AGB_yield,
                               stars->AGB_n_mass, stars->yield_mass_bins[j]);

          index = row_major_index_3d(i, elem, j, stars->AGB_n_z,
                                     chemistry_element_count, n_mass_bins);
          stars->yield_AGB.SPH[index] =
              exp(M_LN10 * stars->yield_mass_bins[j]) * result;
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
inline static void compute_ejecta(struct stars_props *restrict stars) {

  gsl_interp_accel *accel_ptr;

  gsl_spline *SNII_spline_ptr, *AGB_spline_ptr;

  // Do we really need SNII_ejecta and AGB_ejecta, they're not used
  // simultaneously, so can use only one?
  double SNII_ejecta[stars->SNII_n_mass];
  double AGB_ejecta[stars->AGB_n_mass];
  float result;

  accel_ptr = gsl_interp_accel_alloc();
  SNII_spline_ptr = gsl_spline_alloc(gsl_interp_linear, stars->SNII_n_mass);

  int index;

  /* Resample SNII ejecta */
  for (int i = 0; i < stars->SNII_n_z; i++) {
    for (int k = 0; k < stars->SNII_n_mass; k++) {
      index = row_major_index_2d(i, k, stars->SNII_n_z, stars->SNII_n_mass);
      SNII_ejecta[k] = stars->yield_SNII.ejecta[index] *
                      exp(M_LN10 * (-stars->yield_SNII.mass[k]));
    }

    gsl_spline_init(SNII_spline_ptr, stars->yield_SNII.mass, SNII_ejecta,
                    stars->SNII_n_mass);

    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_SNII.mass[0])
        result = SNII_ejecta[0];
      else if (stars->yield_mass_bins[k] >
               stars->yield_SNII.mass[stars->SNII_n_mass - 1])
        result = SNII_ejecta[stars->SNII_n_mass - 1];
      else
        result = gsl_spline_eval(SNII_spline_ptr, stars->yield_mass_bins[k],
                                 accel_ptr);

      index = row_major_index_2d(i, k, stars->SNII_n_z, n_mass_bins);
      stars->yield_SNII.ejecta_SPH[index] =
          exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    }
  }

  /* resample SNII total metals released */
  for (int i = 0; i < stars->SNII_n_z; i++) {
    for (int k = 0; k < stars->SNII_n_mass; k++) {
      index = row_major_index_2d(i, k, stars->SNII_n_z, stars->SNII_n_mass);
      SNII_ejecta[k] = stars->yield_SNII.total_metals[index] *
                      exp(M_LN10 * (-stars->yield_SNII.mass[k]));
    }

    gsl_spline_init(SNII_spline_ptr, stars->yield_SNII.mass, SNII_ejecta,
                    stars->SNII_n_mass);

    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_SNII.mass[0])
        result = SNII_ejecta[0];
      else if (stars->yield_mass_bins[k] >
               stars->yield_SNII.mass[stars->SNII_n_mass - 1])
        result = SNII_ejecta[stars->SNII_n_mass - 1];
      else
        result = gsl_spline_eval(SNII_spline_ptr, stars->yield_mass_bins[k],
                                 accel_ptr);

      index = row_major_index_2d(i, k, stars->SNII_n_z, n_mass_bins);
      stars->yield_SNII.total_metals_SPH[index] =
          exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    }
  }

  gsl_spline_free(SNII_spline_ptr);

  /* AGB yields */
  AGB_spline_ptr = gsl_spline_alloc(gsl_interp_linear, stars->AGB_n_mass);

  for (int i = 0; i < stars->AGB_n_z; i++) {
    for (int k = 0; k < stars->AGB_n_mass; k++) {
      index = row_major_index_2d(i, k, stars->AGB_n_z, stars->AGB_n_mass);
      AGB_ejecta[k] = stars->yield_AGB.ejecta[index] /
                     exp(M_LN10 * stars->yield_AGB.mass[k]);
    }

    gsl_spline_init(AGB_spline_ptr, stars->yield_AGB.mass, AGB_ejecta,
                    stars->AGB_n_mass);

    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_AGB.mass[0])
        result = AGB_ejecta[0];
      else if (stars->yield_mass_bins[k] >
               stars->yield_AGB.mass[stars->AGB_n_mass - 1])
        result = AGB_ejecta[stars->AGB_n_mass - 1];
      else
        result = gsl_spline_eval(AGB_spline_ptr, stars->yield_mass_bins[k],
                                 accel_ptr);

      index = row_major_index_2d(i, k, stars->AGB_n_z, n_mass_bins);
      stars->yield_AGB.ejecta_SPH[index] =
          exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    }
  }

  for (int i = 0; i < stars->AGB_n_z; i++) {
    for (int k = 0; k < stars->AGB_n_mass; k++) {
      index = row_major_index_2d(i, k, stars->AGB_n_z, stars->AGB_n_mass);
      AGB_ejecta[k] = stars->yield_AGB.total_metals[index] *
                     exp(M_LN10 * (-stars->yield_AGB.mass[k]));
    }

    gsl_spline_init(AGB_spline_ptr, stars->yield_AGB.mass, AGB_ejecta,
                    stars->AGB_n_mass);

    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_AGB.mass[0])
        result = AGB_ejecta[0];
      else if (stars->yield_mass_bins[k] >
               stars->yield_AGB.mass[stars->AGB_n_mass - 1])
        result = AGB_ejecta[stars->AGB_n_mass - 1];
      else
        result = gsl_spline_eval(AGB_spline_ptr, stars->yield_mass_bins[k],
                                 accel_ptr);

      index = row_major_index_2d(i, k, stars->AGB_n_z, n_mass_bins);
      stars->yield_AGB.total_metals_SPH[index] =
          exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    }
  }

  gsl_spline_free(AGB_spline_ptr);
  gsl_interp_accel_free(accel_ptr);
}

#endif
