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

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "chemistry.h"

static const float log_min_metallicity = -20;
static const int n_mass_bins = 200; // temporary, put in correct value and move elsewhere.

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

inline static int get_element_index(const char *element_name, char **element_array, int n_elements){
  
  /* Compare element name we are trying to index to every name in element array  */
  for (int i = 0; i < n_elements; i++) {
    if (strcmp(element_array[i], element_name) == 0) return i;
  }

  /* If we don't find the index return flag  */
  return -1;
}

inline static void read_yield_tables(struct stars_props *restrict stars){
#ifdef HAVE_HDF5

  int i, j, k, index;

  char fname[256], setname[100];
  
  hid_t file_id, dataset, datatype;
  herr_t status;

  /* Read SNIa tables */
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

  // What is this for? Copied from EAGLE
  //for (i = 0; i < stars->SNIa_n_elements; i++) {
  //  stars->SNIa_element_names[i] = mystrdup(stars->SNIa_element_names[i]);
  //}

  dataset = H5Dopen(file_id, "Yield", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->yields_SNIa);
  if (status < 0) error("error reading SNIa yields");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  dataset = H5Dopen(file_id, "Total_Metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &stars->yield_SNIa_total_metals_SPH);
  if (status < 0) error("error reading SNIa total metal");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  status = H5Fclose(file_id);
  if (status < 0) error("error closing SNIa file");

  /* Read SNII tables */

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

  // What is this for (again)? copied from EAGLE
  //for (i = 0; i < stars->SNII_n_elements; i++)
  //  yieldsSNII.ElementName[i] = mystrdup(yieldsSNII.ElementName[i]);

  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->yield_SNII.mass);
  if (status < 0) error("error reading SNII masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->yield_SNII.metallicity);
  if (status < 0) error("error reading SNII metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  double tempyield1[stars->SNII_n_elements][stars->SNII_n_mass];

  double tempej1[stars->SNII_n_mass], tempmet1[stars->SNII_n_mass];

  char *tempname1[stars->SNII_n_z];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname1);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  for (i = 0; i < stars->SNII_n_z; i++) {
    sprintf(setname, "/Yields/%s/Yield", tempname1[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            tempyield1);
    if (status < 0) error("error reading SNII yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    sprintf(setname, "/Yields/%s/Ejected_mass", tempname1[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempej1);
    if (status < 0) error("error reading SNII ejected masses");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    sprintf(setname, "/Yields/%s/Total_Metals", tempname1[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            tempmet1);
    if (status < 0) error("error reading SNII total metals");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    for (k = 0; k < stars->SNII_n_mass; k++) {
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      stars->yield_SNII.ejecta[index] = tempej1[k];
      stars->yield_SNII.total_metals[index] = tempmet1[k];

      for (j = 0; j < stars->SNII_n_elements; j++) {
        index = row_major_index_3d(i,j,k,stars->SNII_n_z,stars->SNII_n_elements,stars->SNII_n_mass);
        stars->yield_SNII.yield[index] = tempyield1[j][k];
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

  // What is this for (again)? copied from EAGLE
  //for (i = 0; i < stars->AGB_n_elements; i++)
  //  yieldsAGB.ElementName[i] = mystrdup(yieldsAGB.ElementName[i]);

  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->yield_AGB.mass);
  if (status < 0) error("error reading AGB masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->yield_AGB.metallicity);
  if (status < 0) error("error reading AGB metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  double tempyield2[stars->AGB_n_elements][stars->AGB_n_mass];

  double tempej2[stars->AGB_n_mass], tempmet2[stars->AGB_n_mass];

  char *tempname2[stars->AGB_n_z];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname2);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  for (i = 0; i < stars->AGB_n_z; i++) {
    sprintf(setname, "/Yields/%s/Yield", tempname2[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            tempyield2);
    if (status < 0) error("error reading AGB yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    sprintf(setname, "/Yields/%s/Ejected_mass", tempname2[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempej2);
    if (status < 0) error("error reading AGB ejected masses");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    sprintf(setname, "/Yields/%s/Total_Metals", tempname2[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            tempmet2);
    if (status < 0) error("error reading AGB total metals");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");

    for (k = 0; k < stars->AGB_n_mass; k++) {
      index = row_major_index_2d(i,k,stars->AGB_n_z,stars->AGB_n_mass);
      stars->yield_AGB.ejecta[index] = tempej2[k];
      stars->yield_AGB.total_metals[index] = tempmet2[k];

      for (j = 0; j < stars->AGB_n_elements; j++) {
        index = row_major_index_3d(i,j,k,stars->AGB_n_z,stars->AGB_n_elements,stars->AGB_n_mass);
        stars->yield_AGB.yield[index] = tempyield2[j][k];
      }
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* Read lifetimes table */

  sprintf(fname, "%s/Lifetimes.hdf5", stars->yield_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) error("unable to open file %s\n", fname);

  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->lifetimes.mass);
  if (status < 0) error("error reading lifetime table masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          stars->lifetimes.metallicity);
  if (status < 0) error("error reading AGB metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  double temptime[stars->lifetimes.n_z][stars->lifetimes.n_mass];

  dataset = H5Dopen(file_id, "Lifetimes", H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temptime);
  H5Dclose(dataset);

  for (i = 0; i < stars->lifetimes.n_z; i++) {
    for (j = 0; j < stars->lifetimes.n_mass; j++) {
      //index = row_major_index_2d(i,j,stars->lifetimes.n_z,stars->lifetimes.n_mass);
      //stars->lifetimes.dyingtime[index] = log10(temptime[i][j]);
      stars->lifetimes.dyingtime[i][j] = log10(temptime[i][j]);
    }
  }

  H5Fclose(file_id);


#endif
}

inline static void allocate_yield_tables(struct stars_props *restrict stars){

  /* Allocate SNIa arrays */
  if (posix_memalign((void **)&stars->yields_SNIa, SWIFT_STRUCT_ALIGNMENT, stars->SNIa_n_elements * sizeof(double)) !=0) {
    error("Failed to allocate SNIa yields array");
  }
  if (posix_memalign((void **)&stars->yield_SNIa_SPH, SWIFT_STRUCT_ALIGNMENT, stars->SNIa_n_elements * sizeof(double)) !=0) {
    error("Failed to allocate SNIa SPH yields array");
  }

  /* Allocate AGB arrays  */
  if (posix_memalign((void **)&stars->yield_AGB.mass, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_mass * sizeof(double)) !=0) {
    error("Failed to allocate AGB mass array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.metallicity, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * sizeof(double)) !=0) {
    error("Failed to allocate AGB metallicity array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.SPH, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * n_mass_bins * stars->AGB_n_elements * sizeof(double)) !=0) {
    error("Failed to allocate AGB SPH array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.yield, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * stars->AGB_n_mass * stars->AGB_n_elements * sizeof(double)) !=0) {
    error("Failed to allocate AGB yield array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.ejecta_SPH, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * n_mass_bins * sizeof(double)) !=0) {
    error("Failed to allocate AGB ejecta SPH array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.ejecta, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * stars->AGB_n_mass * sizeof(double)) !=0) {
    error("Failed to allocate AGB ejecta array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.total_metals_SPH, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * n_mass_bins * sizeof(double)) !=0) {
    error("Failed to allocate AGB total metals SPH array");
  }
  if (posix_memalign((void **)&stars->yield_AGB.total_metals, SWIFT_STRUCT_ALIGNMENT, stars->AGB_n_z * stars->AGB_n_mass * sizeof(double)) !=0) {
    error("Failed to allocate AGB total metals array");
  }

  /* Allocate SNII arrays  */
  if (posix_memalign((void **)&stars->yield_SNII.mass, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_mass * sizeof(double)) !=0) {
    error("Failed to allocate SNII mass array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.metallicity, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * sizeof(double)) !=0) {
    error("Failed to allocate SNII metallicity array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.SPH, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * n_mass_bins * stars->SNII_n_elements * sizeof(double)) !=0) {
    error("Failed to allocate SNII SPH array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.yield, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * stars->SNII_n_mass * stars->SNII_n_elements * sizeof(double)) !=0) {
    error("Failed to allocate SNII yield array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.ejecta_SPH, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * n_mass_bins * sizeof(double)) !=0) {
    error("Failed to allocate SNII ejecta SPH array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.ejecta, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * stars->SNII_n_mass * sizeof(double)) !=0) {
    error("Failed to allocate SNII ejecta array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.total_metals_SPH, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * n_mass_bins * sizeof(double)) !=0) {
    error("Failed to allocate SNII total metals SPH array");
  }
  if (posix_memalign((void **)&stars->yield_SNII.total_metals, SWIFT_STRUCT_ALIGNMENT, stars->SNII_n_z * stars->SNII_n_mass * sizeof(double)) !=0) {
    error("Failed to allocate SNII total metals array");
  }

  /* Allocate Lifetime arrays  */
  if (posix_memalign((void **)&stars->lifetimes.mass, SWIFT_STRUCT_ALIGNMENT, stars->lifetimes.n_mass * sizeof(double)) !=0) {
    error("Failed to allocate lifetime mass array");
  }
  if (posix_memalign((void **)&stars->lifetimes.metallicity, SWIFT_STRUCT_ALIGNMENT, stars->lifetimes.n_z * sizeof(double)) !=0) {
    error("Failed to allocate lifetime metallicity array");
  }
  //if (posix_memalign((void **)&stars->lifetimes.dyingtime, SWIFT_STRUCT_ALIGNMENT, stars->lifetimes.n_z * stars->lifetimes.n_mass * sizeof(double)) !=0) {
  //  error("Failed to allocate dyingtime array");
  //}
  stars->lifetimes.dyingtime = (double **)malloc(stars->lifetimes.n_z * sizeof(double *));
  for(int i = 0; i < stars->lifetimes.n_z; i++){
    stars->lifetimes.dyingtime[i] = (double *)malloc(stars->lifetimes.n_mass * sizeof(double));
  }

  /* Allocate SNII factor array  */
  if (posix_memalign((void **)&stars->typeII_factor, SWIFT_STRUCT_ALIGNMENT, chemistry_element_count * sizeof(double)) !=0) {
    error("Failed to allocate SNII factor array");
  }

  /* Allocate element name arrays  */
  stars->SNIa_element_names = (char **)malloc(stars->SNIa_n_elements * sizeof(char*));
  for(int i = 0; i < stars->SNIa_n_elements; i++){
    stars->SNIa_element_names[i] = (char *)malloc(stars->element_name_length * sizeof(char));
  }
  stars->SNII_element_names = (char **)malloc(stars->SNII_n_elements * sizeof(char*));
  for(int i = 0; i < stars->SNII_n_elements; i++){
    stars->SNII_element_names[i] = (char *)malloc(stars->element_name_length * sizeof(char));
  }
  stars->AGB_element_names = (char **)malloc(stars->AGB_n_elements * sizeof(char*));
  for(int i = 0; i < stars->AGB_n_elements; i++){
    stars->AGB_element_names[i] = (char *)malloc(stars->element_name_length * sizeof(char));
  }
  
  /* Allocate array of mass bins for yield and ejecta calculation */
  if (posix_memalign((void **)&stars->yield_mass_bins, SWIFT_STRUCT_ALIGNMENT, n_mass_bins * sizeof(double)) !=0) {
    error("Failed to allocate yield mass bins array");
  }
}

inline static double interpolate_1D(double *array_x, double *array_y, int size, double value){
  
  double result;

  if (value < array_x[0]) 
    error("interpolating value less than array min. value %.5e array min %.5e", value, array_x[0]);
  else if (value > array_x[size-1]) 
    error("interpolating value greater than array max. value %.5e array max %.5e", value, array_x[size-1]);
  else {
    int index = 0;
    while (array_x[index] <= value) index++;
    double offset = (array_x[index] - value)/(array_x[index] - array_x[index-1]);
    result = offset*array_y[index-1] + (1 - offset)*array_y[index];
  }

  return result;
}

inline static void compute_yields(struct stars_props *restrict stars){

  int index, index_2d;
  gsl_interp_accel *accel_ptr;

  gsl_spline *SNII_spline_ptr, *AGB_spline_ptr;

  // check whether imf is properly normalised? present in EAGLE

  /* convert SNII tables to log10  */
  for (int i = 0; i < stars->SNII_n_mass; i++) {
    stars->yield_SNII.mass[i] = log10(stars->yield_SNII.mass[i]);
  }
  for (int i = 0; i < stars->SNII_n_z; i++) {
    if (stars->yield_SNII.metallicity[i] > 0) {
      stars->yield_SNII.metallicity[i] = log10(stars->yield_SNII.metallicity[i]);
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

  accel_ptr = gsl_interp_accel_alloc();
  SNII_spline_ptr = gsl_spline_alloc(gsl_interp_linear, stars->SNII_n_mass);
  AGB_spline_ptr = gsl_spline_alloc(gsl_interp_linear, stars->AGB_n_mass);
  double SNII_yield[stars->SNII_n_mass];
  double AGB_yield[stars->AGB_n_mass];
  float result;


  /* Loop over elements tracked in EAGLE  */
  int element_index = 0;
  // ALEXEI: better name for eagle_elem?
  for (enum chemistry_element eagle_elem = chemistry_element_H; eagle_elem < chemistry_element_count; eagle_elem++) {
    /* SNIa  */
    element_index = get_element_index(chemistry_get_element_name(eagle_elem), stars->SNIa_element_names, stars->SNIa_n_elements);
    stars->yield_SNIa_SPH[eagle_elem] = stars->yields_SNIa[element_index];

    /* SNII  */
    element_index = get_element_index(chemistry_get_element_name(eagle_elem), stars->SNII_element_names, stars->SNII_n_elements);
    for (int i = 0; i < stars->SNII_n_z; i++) {
      for (int j = 0; j < stars->SNII_n_mass; j++) {
        index = row_major_index_3d(i, element_index, j, stars->SNII_n_z, stars->SNII_n_elements, stars->SNII_n_mass);
	SNII_yield[j] = stars->yield_SNII.yield[index] * exp(M_LN10 * (-stars->yield_SNII.mass[j]));
	//if (SNII_yield[j] < 0) error("SNII_yield negative %.5e z_i %d element_index %d mass_i %d table yield %.5e factor %.5e", SNII_yield[j], i, element_index, j, stars->yield_SNII.yield[index], exp(M_LN10 * (-stars->yield_SNII.mass[j])));
      }

      // ALEXEI: for some reason need to define another accel_ptr so that it doesn't crash due to thinking we're trying to interpolate out of bounds when we're not. Investigate more?
      //gsl_interp_accel *SNII_accel_ptr;
      //gsl_spline *SNII_spline_ptr2;
      //SNII_accel_ptr = gsl_interp_accel_alloc();
      //SNII_spline_ptr2 = gsl_spline_alloc(gsl_interp_linear, stars->SNII_n_mass);

      gsl_spline_init(SNII_spline_ptr, stars->yield_SNII.mass, SNII_yield, stars->SNII_n_mass);

      for (int k = 0; k < n_mass_bins; k++) {
        if (stars->yield_mass_bins[k] < stars->yield_SNII.mass[0])
          result = SNII_yield[0];
        else if (stars->yield_mass_bins[k] > stars->yield_SNII.mass[stars->SNII_n_mass - 1])
          result = SNII_yield[stars->SNII_n_mass - 1];
        else {
	  //message("k %d yield_mass %.5e min max mass %.5e %.5e n_mass %d", k, stars->yield_mass_bins[k], stars->yield_SNII.mass[0], stars->yield_SNII.mass[stars->SNII_n_mass - 1], stars->SNII_n_mass);
          //result =
          //    gsl_spline_eval(SNII_spline_ptr, stars->yield_mass_bins[k], accel_ptr);
	  result = interpolate_1D(stars->yield_SNII.mass,SNII_yield, stars->SNII_n_mass,stars->yield_mass_bins[k]);
	  //message("result %.5e", result);
	}

        index = row_major_index_3d(i,eagle_elem,k,stars->SNII_n_z,chemistry_element_count,n_mass_bins);
        stars->yield_SNII.SPH[index] = exp(M_LN10 * stars->yield_mass_bins[k]) * result;
      }

      //gsl_interp_accel_free(SNII_accel_ptr);
      //gsl_spline_free(SNII_spline_ptr2);
    }

    for (int i = 0; i < stars->SNII_n_z; i++) {
      for (int k = 0; k < n_mass_bins; k++) {
        index_2d = row_major_index_2d(i,k,stars->SNII_n_z,n_mass_bins);
        index = row_major_index_3d(i,eagle_elem,k,stars->SNII_n_z,chemistry_element_count,n_mass_bins);
        if (strcmp(chemistry_get_element_name(eagle_elem), "Hydrogen") != 0 ||
            strcmp(chemistry_get_element_name(eagle_elem), "Helium") != 0) {
          stars->yield_SNII.total_metals_SPH[index_2d] +=
              (stars->typeII_factor[eagle_elem] - 1) * stars->yield_SNII.SPH[index];
        }

        stars->yield_SNII.SPH[index] *= stars->typeII_factor[eagle_elem];
      }
    }

    /* AGB  */
    element_index = get_element_index(chemistry_get_element_name(eagle_elem), stars->AGB_element_names, stars->AGB_n_elements);
    if (element_index < 0) {
      for (int i = 0; i < stars->AGB_n_z; i++) {
        for (int j = 0; j < n_mass_bins; j++) {
	  index = row_major_index_3d(i,eagle_elem,j,stars->AGB_n_z,chemistry_element_count,n_mass_bins);
          stars->yield_AGB.SPH[index] = 0.0;
	}
      }
    } else {
      for (int i = 0; i < stars->AGB_n_z; i++) {
        for (int j = 0; j < stars->AGB_n_mass; j++) {
          index = row_major_index_3d(i, element_index, j, stars->AGB_n_z, stars->AGB_n_elements, stars->AGB_n_mass);
          AGB_yield[j] = stars->yield_AGB.yield[index] / exp(M_LN10 * stars->yield_AGB.mass[j]);
        }

        gsl_spline_init(AGB_spline_ptr, stars->yield_AGB.mass, AGB_yield, stars->AGB_n_mass);

        for (int j = 0; j < n_mass_bins; j++) {
          if (stars->yield_mass_bins[j] < stars->yield_AGB.mass[0])
            result = AGB_yield[0];
          else if (stars->yield_mass_bins[j] > stars->yield_AGB.mass[stars->AGB_n_mass - 1])
            result = AGB_yield[stars->AGB_n_mass - 1];
          else
            result =
                gsl_spline_eval(AGB_spline_ptr, stars->yield_mass_bins[j], accel_ptr);

          index = row_major_index_3d(i,eagle_elem,j,stars->AGB_n_z,chemistry_element_count,n_mass_bins);
          stars->yield_AGB.SPH[index] = exp(M_LN10 * stars->yield_mass_bins[j]) * result;
        }
      }
    }
  }
  gsl_spline_free(SNII_spline_ptr);
  gsl_spline_free(AGB_spline_ptr);
  gsl_interp_accel_free(accel_ptr);

}

inline static void compute_ejecta(struct stars_props *restrict stars) {
  
  gsl_interp_accel *accel_ptr;

  gsl_spline *SNII_spline_ptr, *AGB_spline_ptr;

  // Do we really need SNII_yield and AGB_yield, they're not used simultaneously, so can use only one?
  double SNII_yield[stars->SNII_n_mass];
  double AGB_yield[stars->AGB_n_mass];
  float result;
  
  accel_ptr = gsl_interp_accel_alloc();                                                                            
  SNII_spline_ptr = gsl_spline_alloc(gsl_interp_linear, stars->SNII_n_mass);                                             

  int index;
    
  for (int i = 0; i < stars->SNII_n_z; i++) {
    for (int k = 0; k < stars->SNII_n_mass; k++) {
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      SNII_yield[k] = stars->yield_SNII.ejecta[index] * exp(M_LN10 * (-stars->yield_SNII.mass[k]));                                             
    }

    gsl_spline_init(SNII_spline_ptr, stars->yield_SNII.mass, SNII_yield, stars->SNII_n_mass);                                        
    
    for (int k = 0; k < n_mass_bins; k++) {                                                                                
      if (stars->yield_mass_bins[k] < stars->yield_SNII.mass[0])                                                                     
        result = SNII_yield[0];
      else if (stars->yield_mass_bins[k] > stars->yield_SNII.mass[stars->SNII_n_mass - 1])                                            
        result = SNII_yield[stars->SNII_n_mass - 1];
      else
        result =                                                                                                      
            gsl_spline_eval(SNII_spline_ptr, stars->yield_mass_bins[k], accel_ptr);
            
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      stars->yield_SNII.ejecta_SPH[index] = exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    } 
  }       
    
  for (int i = 0; i < stars->SNII_n_z; i++) {                                                      
    for (int k = 0; k < stars->SNII_n_mass; k++) {
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      SNII_yield[k] = stars->yield_SNII.total_metals[index] * exp(M_LN10 * (-stars->yield_SNII.mass[k]));                                        
    }
    
    gsl_spline_init(SNII_spline_ptr, stars->yield_SNII.mass, SNII_yield, stars->SNII_n_mass);                                        
      
    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_SNII.mass[0])                                                                     
        result = SNII_yield[0];                                                                                            
      else if (stars->yield_mass_bins[k] > stars->yield_SNII.mass[stars->SNII_n_mass - 1])                                            
        result = SNII_yield[stars->SNII_n_mass - 1];                                                                        
      else 
        result =
            gsl_spline_eval(SNII_spline_ptr, stars->yield_mass_bins[k], accel_ptr);                                          
          
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      stars->yield_SNII.total_metals_SPH[index] = exp(M_LN10 * stars->yield_mass_bins[k]) * result;                                       
    }       
  }                                                                                                                   

  gsl_spline_free(SNII_spline_ptr);

  /* AGB yields */
  AGB_spline_ptr = gsl_spline_alloc(gsl_interp_linear, stars->AGB_n_mass);

  for (int i = 0; i < stars->AGB_n_z; i++) {
    for (int k = 0; k < stars->AGB_n_mass; k++) {
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      AGB_yield[k] = stars->yield_AGB.ejecta[index] / exp(M_LN10 * stars->yield_AGB.mass[k]);
    }

    gsl_spline_init(AGB_spline_ptr, stars->yield_AGB.mass, AGB_yield, stars->AGB_n_mass);

    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_AGB.mass[0])
        result = AGB_yield[0];
      else if (stars->yield_mass_bins[k] > stars->yield_AGB.mass[stars->AGB_n_mass - 1])
        result = AGB_yield[stars->AGB_n_mass - 1];
      else
        result =
            gsl_spline_eval(AGB_spline_ptr, stars->yield_mass_bins[k], accel_ptr);

      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      stars->yield_AGB.ejecta_SPH[index] = exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    }
  }

  for (int i = 0; i < stars->AGB_n_z; i++) {
    for (int k = 0; k < stars->AGB_n_mass; k++) {
      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      AGB_yield[k] = stars->yield_AGB.total_metals[index] * exp(M_LN10 * (-stars->yield_AGB.mass[k]));
    }

    gsl_spline_init(AGB_spline_ptr, stars->yield_AGB.mass, AGB_yield, stars->AGB_n_mass);

    for (int k = 0; k < n_mass_bins; k++) {
      if (stars->yield_mass_bins[k] < stars->yield_AGB.mass[0])
        result = AGB_yield[0];
      else if (stars->yield_mass_bins[k] > stars->yield_AGB.mass[stars->AGB_n_mass - 1])
        result = AGB_yield[stars->AGB_n_mass - 1];
      else
        result =
            gsl_spline_eval(AGB_spline_ptr, stars->yield_mass_bins[k], accel_ptr);

      index = row_major_index_2d(i,k,stars->SNII_n_z,stars->SNII_n_mass);
      stars->yield_AGB.total_metals_SPH[index] = exp(M_LN10 * stars->yield_mass_bins[k]) * result;
    }
  }

  gsl_spline_free(AGB_spline_ptr);
  gsl_interp_accel_free(accel_ptr);
}


#endif
