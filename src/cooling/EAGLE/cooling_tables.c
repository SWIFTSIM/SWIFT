/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 * @file src/cooling/EAGLE/cooling_tables.c
 * @brief Functions to read EAGLE tables
 */

/* Config parameters. */
#include "../config.h"

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "cooling_tables.h"
#include "error.h"
#include "interpolate.h"

/**
 * @brief Names of the elements in the order they are stored in the files
 */
static const char *eagle_tables_element_names[eagle_cooling_N_metal] = {
    "Carbon",  "Nitrogen", "Oxygen",  "Neon", "Magnesium",
    "Silicon", "Sulphur",  "Calcium", "Iron"};

/*! Number of elements in a z-slice of the H+He cooling rate tables */
static const size_t num_elements_cooling_rate =
    eagle_cooling_N_temperature * eagle_cooling_N_density;

/*! Number of elements in a z-slice of the metal cooling rate tables */
static const size_t num_elements_metal_heating = eagle_cooling_N_metal *
                                                 eagle_cooling_N_temperature *
                                                 eagle_cooling_N_density;

/*! Number of elements in a z-slice of the metal electron abundance tables */
static const size_t num_elements_electron_abundance =
    eagle_cooling_N_temperature * eagle_cooling_N_density;

/*! Number of elements in a z-slice of the temperature tables */
static const size_t num_elements_temperature = eagle_cooling_N_He_frac *
                                               eagle_cooling_N_temperature *
                                               eagle_cooling_N_density;

/*! Number of elements in a z-slice of the H+He cooling rate tables */
static const size_t num_elements_HpHe_heating = eagle_cooling_N_He_frac *
                                                eagle_cooling_N_temperature *
                                                eagle_cooling_N_density;

/*! Number of elements in a z-slice of the H+He electron abundance tables */
static const size_t num_elements_HpHe_electron_abundance =
    eagle_cooling_N_He_frac * eagle_cooling_N_temperature *
    eagle_cooling_N_density;

/**
 * @brief Reads in EAGLE table of redshift values
 *
 * @param cooling #cooling_function_data structure
 */
void get_cooling_redshifts(struct cooling_function_data *cooling) {

  /* Read the list of table redshifts */
  char redshift_filename[eagle_table_path_name_length + 16];
  sprintf(redshift_filename, "%s/redshifts.dat", cooling->cooling_table_path);

  FILE *infile = fopen(redshift_filename, "r");
  if (infile == NULL) {
    error("Cannot open the list of cooling table redshifts (%s)",
          redshift_filename);
  }

  int N_Redshifts = -1;

  /* Read the file */
  if (!feof(infile)) {

    char buffer[50];

    /* Read the number of redshifts (1st line in the file) */
    if (fgets(buffer, 50, infile) != NULL)
      N_Redshifts = atoi(buffer);
    else
      error("Impossible to read the number of redshifts");

    /* Be verbose about it */
    message("Found cooling tables at %d redhsifts", N_Redshifts);

    /* Check value */
    if (N_Redshifts != eagle_cooling_N_redshifts)
      error("Invalid redshift lenght array.");

    /* Allocate the list of redshifts */
    if (posix_memalign((void **)&cooling->Redshifts, SWIFT_STRUCT_ALIGNMENT,
                       eagle_cooling_N_redshifts * sizeof(float)) != 0)
      error("Failed to allocate redshift table");

    /* Read all the redshift values */
    int count = 0;
    while (!feof(infile)) {
      if (fgets(buffer, 50, infile) != NULL) {
        cooling->Redshifts[count] = atof(buffer);
        count++;
      }
    }

    /* Verify that the file was self-consistent */
    if (count != N_Redshifts) {
      error(
          "Redshift file (%s) does not contain the correct number of redshifts "
          "(%d vs. %d)",
          redshift_filename, count, N_Redshifts);
    }
  } else {
    error("Redshift file (%s) is empty!", redshift_filename);
  }

  /* We are done with this file */
  fclose(infile);

  /* EAGLE cooling assumes cooling->Redshifts table is in increasing order. Test
   * this. */
  for (int i = 0; i < N_Redshifts - 1; i++) {
    if (cooling->Redshifts[i + 1] < cooling->Redshifts[i]) {
      error("table should be in increasing order\n");
    }
  }
}

/**
 * @brief Reads in EAGLE cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, helium fraction
 * solar element abundances, and elements used to index the cooling tables.
 *
 * @param fname Filepath for cooling table from which to read header
 * @param cooling Cooling data structure
 */
void read_cooling_header(const char *fname,
                         struct cooling_function_data *cooling) {

#ifdef HAVE_HDF5

  int N_Temp, N_nH, N_He, N_SolarAbundances, N_Elements;

  /* read sizes of array dimensions */
  hid_t tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0) error("unable to open file %s\n", fname);

  /* read size of each table of values */
  hid_t dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
  herr_t status =
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_Temp);
  if (status < 0) error("error reading number of temperature bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Check value */
  if (N_Temp != eagle_cooling_N_temperature)
    error("Invalid temperature array length.");

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins", H5P_DEFAULT);
  status =
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_nH);
  if (status < 0) error("error reading number of density bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Check value */
  if (N_nH != eagle_cooling_N_density) error("Invalid density array length.");

  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
  status =
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_He);
  if (status < 0) error("error reading number of He fraction bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Check value */
  if (N_He != eagle_cooling_N_He_frac)
    error("Invalid Helium fraction array length.");

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &N_SolarAbundances);
  if (status < 0) error("error reading number of solar abundance bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Check value */
  if (N_SolarAbundances != eagle_cooling_N_abundances)
    error("Invalid solar abundances array length.");

  /* Check value */
  if (N_SolarAbundances != chemistry_element_count + 2)
    error("Number of abundances not compatible with the chemistry model.");

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &N_Elements);
  if (status < 0) error("error reading number of metal bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Check value */
  if (N_Elements != eagle_cooling_N_metal) error("Invalid metal array length.");

  /* allocate arrays of values for each of the above quantities */
  if (posix_memalign((void **)&cooling->Temp, SWIFT_STRUCT_ALIGNMENT,
                     N_Temp * sizeof(float)) != 0)
    error("Failed to allocate temperature table");
  if (posix_memalign((void **)&cooling->Therm, SWIFT_STRUCT_ALIGNMENT,
                     N_Temp * sizeof(float)) != 0)
    error("Failed to allocate internal energy table");
  if (posix_memalign((void **)&cooling->nH, SWIFT_STRUCT_ALIGNMENT,
                     N_nH * sizeof(float)) != 0)
    error("Failed to allocate nH table");
  if (posix_memalign((void **)&cooling->HeFrac, SWIFT_STRUCT_ALIGNMENT,
                     N_He * sizeof(float)) != 0)
    error("Failed to allocate HeFrac table");
  if (posix_memalign((void **)&cooling->SolarAbundances, SWIFT_STRUCT_ALIGNMENT,
                     N_SolarAbundances * sizeof(float)) != 0)
    error("Failed to allocate Solar abundances table");
  if (posix_memalign((void **)&cooling->SolarAbundances_inv,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_SolarAbundances * sizeof(float)) != 0)
    error("Failed to allocate Solar abundances inverses table");

  /* read in values for each of the arrays */
  dataset = H5Dopen(tempfile_id, "/Solar/Temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Temp);
  if (status < 0) error("error reading temperature bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Solar/Hydrogen_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->nH);
  if (status < 0) error("error reading H density bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Metal_free/Helium_mass_fraction_bins",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->HeFrac);
  if (status < 0) error("error reading He fraction bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Solar_mass_fractions",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->SolarAbundances);
  if (status < 0) error("error reading solar mass fraction bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Metal_free/Temperature/Energy_density_bins",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Therm);
  if (status < 0) error("error reading internal energy bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Convert to temperature, density and internal energy arrays to log10 */
  for (int i = 0; i < N_Temp; i++) {
    cooling->Temp[i] = log10(cooling->Temp[i]);
    cooling->Therm[i] = log10(cooling->Therm[i]);
  }
  for (int i = 0; i < N_nH; i++) {
    cooling->nH[i] = log10(cooling->nH[i]);
  }
  /* Compute inverse of solar mass fractions */
  for (int i = 0; i < N_SolarAbundances; ++i) {
    cooling->SolarAbundances_inv[i] = 1.f / cooling->SolarAbundances[i];
  }

#else
  error("Need HDF5 to read cooling tables");
#endif
}

/**
 * @brief Allocate space for cooling tables.
 *
 * @param cooling #cooling_function_data structure
 */
void allocate_cooling_tables(struct cooling_function_data *restrict cooling) {

  /* Allocate arrays to store cooling tables. Arrays contain two tables of
   * cooling rates with one table being for the redshift above current redshift
   * and one below. */

  if (posix_memalign((void **)&cooling->table.metal_heating,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_cooling_N_loaded_redshifts *
                         num_elements_metal_heating * sizeof(float)) != 0)
    error("Failed to allocate metal_heating array");

  if (posix_memalign((void **)&cooling->table.electron_abundance,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_cooling_N_loaded_redshifts *
                         num_elements_electron_abundance * sizeof(float)) != 0)
    error("Failed to allocate electron_abundance array");

  if (posix_memalign((void **)&cooling->table.temperature,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_cooling_N_loaded_redshifts *
                         num_elements_temperature * sizeof(float)) != 0)
    error("Failed to allocate temperature array");

  if (posix_memalign((void **)&cooling->table.H_plus_He_heating,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_cooling_N_loaded_redshifts *
                         num_elements_HpHe_heating * sizeof(float)) != 0)
    error("Failed to allocate H_plus_He_heating array");

  if (posix_memalign((void **)&cooling->table.H_plus_He_electron_abundance,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_cooling_N_loaded_redshifts *
                         num_elements_HpHe_electron_abundance *
                         sizeof(float)) != 0)
    error("Failed to allocate H_plus_He_electron_abundance array");
}

/**
 * @brief Get the redshift invariant table of cooling rates (before reionization
 * at redshift ~9) Reads in table of cooling rates and electron abundances due
 * to metals (depending on temperature, hydrogen number density), cooling rates
 * and electron abundances due to hydrogen and helium (depending on temperature,
 * hydrogen number density and helium fraction), and temperatures (depending on
 * internal energy, hydrogen number density and helium fraction; note: this is
 * distinct from table of temperatures read in ReadCoolingHeader, as that table
 * is used to index the cooling, electron abundance tables, whereas this one is
 * used to obtain temperature of particle)
 *
 * @param cooling #cooling_function_data structure
 * @param photodis Are we loading the photo-dissociation table?
 */
void get_redshift_invariant_table(
    struct cooling_function_data *restrict cooling, const int photodis) {
#ifdef HAVE_HDF5

  /* Temporary tables */
  float *net_cooling_rate = NULL;
  float *electron_abundance = NULL;
  float *temperature = NULL;
  float *he_net_cooling_rate = NULL;
  float *he_electron_abundance = NULL;

  /* Allocate arrays for reading in cooling tables.  */
  if (posix_memalign((void **)&net_cooling_rate, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_cooling_rate * sizeof(float)) != 0)
    error("Failed to allocate net_cooling_rate array");
  if (posix_memalign((void **)&electron_abundance, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_electron_abundance * sizeof(float)) != 0)
    error("Failed to allocate electron_abundance array");
  if (posix_memalign((void **)&temperature, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_temperature * sizeof(float)) != 0)
    error("Failed to allocate temperature array");
  if (posix_memalign((void **)&he_net_cooling_rate, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_HpHe_heating * sizeof(float)) != 0)
    error("Failed to allocate he_net_cooling_rate array");
  if (posix_memalign((void **)&he_electron_abundance, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_HpHe_electron_abundance * sizeof(float)) != 0)
    error("Failed to allocate he_electron_abundance array");

  /* Decide which high redshift table to read. Indices set in cooling_update */
  char filename[eagle_table_path_name_length + 21];
  if (photodis) {
    sprintf(filename, "%sz_photodis.hdf5", cooling->cooling_table_path);
    message("Reading cooling table 'z_photodis.hdf5'");
  } else {
    sprintf(filename, "%sz_8.989nocompton.hdf5", cooling->cooling_table_path);
    message("Reading cooling table 'z_8.989nocompton.hdf5'");
  }

  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", filename);

  char set_name[64];

  /* read in cooling rates due to metals */
  for (int specs = 0; specs < eagle_cooling_N_metal; specs++) {

    /* Read in the cooling rate for this metal */
    sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
    hid_t dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, net_cooling_rate);
    if (status < 0) error("error reading metal cooling rate table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    /* Transpose from order tables are stored in (temperature, nH)
     * to (metal species, nH, temperature) where fastest
     * varying index is on right. Tables contain cooling rates but we
     * want rate of change of internal energy, hence minus sign. */
    for (int j = 0; j < eagle_cooling_N_temperature; j++) {
      for (int k = 0; k < eagle_cooling_N_density; k++) {

        /* Index in the HDF5 table */
        const int hdf5_index = row_major_index_2d(
            j, k, eagle_cooling_N_temperature, eagle_cooling_N_density);

        /* Index in the internal table */
        const int internal_index = row_major_index_3d(
            specs, k, j, eagle_cooling_N_metal, eagle_cooling_N_density,
            eagle_cooling_N_temperature);

        /* Change the sign and transpose */
        cooling->table.metal_heating[internal_index] =
            -net_cooling_rate[hdf5_index];
      }
    }
  }

  /* read in cooling rates due to H + He */
  strcpy(set_name, "/Metal_free/Net_Cooling");
  hid_t dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, he_net_cooling_rate);
  if (status < 0) error("error reading metal free cooling rate table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* read in Temperatures */
  strcpy(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  if (status < 0) error("error reading temperature table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* read in H + He electron abundances */
  strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  if (status < 0) error("error reading electron density table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Transpose from order tables are stored in (helium fraction, temperature,
   * nH) to (nH, helium fraction, temperature) where fastest
   * varying index is on right. Tables contain cooling rates but we
   * want rate of change of internal energy, hence minus sign. */
  for (int i = 0; i < eagle_cooling_N_He_frac; i++) {
    for (int j = 0; j < eagle_cooling_N_temperature; j++) {
      for (int k = 0; k < eagle_cooling_N_density; k++) {

        /* Index in the HDF5 table */
        const int hdf5_index = row_major_index_3d(
            i, j, k, eagle_cooling_N_He_frac, eagle_cooling_N_temperature,
            eagle_cooling_N_density);

        /* Index in the internal table */
        const int internal_index = row_major_index_3d(
            k, i, j, eagle_cooling_N_density, eagle_cooling_N_He_frac,
            eagle_cooling_N_temperature);

        /* Change the sign and transpose */
        cooling->table.H_plus_He_heating[internal_index] =
            -he_net_cooling_rate[hdf5_index];

        /* Convert to log T and transpose */
        cooling->table.temperature[internal_index] =
            log10(temperature[hdf5_index]);

        /* Just transpose */
        cooling->table.H_plus_He_electron_abundance[internal_index] =
            he_electron_abundance[hdf5_index];
      }
    }
  }

  /* read in electron densities due to metals */
  strcpy(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  if (status < 0) error("error reading solar electron density table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Transpose from order tables are stored in (temperature, nH) to
   * (nH, temperature) where fastest varying index is on right. */
  for (int i = 0; i < eagle_cooling_N_temperature; i++) {
    for (int j = 0; j < eagle_cooling_N_density; j++) {

      /* Index in the HDF5 table */
      const int hdf5_index = row_major_index_2d(
          i, j, eagle_cooling_N_temperature, eagle_cooling_N_density);

      /* Index in the internal table */
      const int internal_index = row_major_index_2d(
          j, i, eagle_cooling_N_density, eagle_cooling_N_temperature);

      /* Just transpose */
      cooling->table.electron_abundance[internal_index] =
          electron_abundance[hdf5_index];
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

#ifdef SWIFT_DEBUG_CHECKS
  message("done reading in redshift invariant table");
#endif

#else
  error("Need HDF5 to read cooling tables");
#endif
}

/**
 * @brief Get redshift dependent table of cooling rates.
 * Reads in table of cooling rates and electron abundances due to
 * metals (depending on temperature, hydrogen number density), cooling rates and
 * electron abundances due to hydrogen and helium (depending on temperature,
 * hydrogen number density and helium fraction), and temperatures (depending on
 * internal energy, hydrogen number density and helium fraction; note: this is
 * distinct from table of temperatures read in ReadCoolingHeader, as that table
 * is used to index the cooling, electron abundance tables, whereas this one is
 * used to obtain temperature of particle)
 *
 * @param cooling #cooling_function_data structure
 * @param low_z_index Index of the lowest redshift table to load.
 * @param high_z_index Index of the highest redshift table to load.
 */
void get_cooling_table(struct cooling_function_data *restrict cooling,
                       const int low_z_index, const int high_z_index) {

#ifdef HAVE_HDF5

  /* Temporary tables */
  float *net_cooling_rate = NULL;
  float *electron_abundance = NULL;
  float *temperature = NULL;
  float *he_net_cooling_rate = NULL;
  float *he_electron_abundance = NULL;

  /* Allocate arrays for reading in cooling tables.  */
  if (posix_memalign((void **)&net_cooling_rate, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_cooling_rate * sizeof(float)) != 0)
    error("Failed to allocate net_cooling_rate array");
  if (posix_memalign((void **)&electron_abundance, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_electron_abundance * sizeof(float)) != 0)
    error("Failed to allocate electron_abundance array");
  if (posix_memalign((void **)&temperature, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_temperature * sizeof(float)) != 0)
    error("Failed to allocate temperature array");
  if (posix_memalign((void **)&he_net_cooling_rate, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_HpHe_heating * sizeof(float)) != 0)
    error("Failed to allocate he_net_cooling_rate array");
  if (posix_memalign((void **)&he_electron_abundance, SWIFT_STRUCT_ALIGNMENT,
                     num_elements_HpHe_electron_abundance * sizeof(float)) != 0)
    error("Failed to allocate he_electron_abundance array");

  /* Read in tables, transpose so that values for indices which vary most are
   * adjacent. Repeat for redshift above and redshift below current value.  */
  for (int z_index = low_z_index; z_index <= high_z_index; z_index++) {

    /* Index along redhsift dimension for the subset of tables we read */
    const int local_z_index = z_index - low_z_index;

#ifdef SWIFT_DEBUG_CHECKS
    if (local_z_index >= eagle_cooling_N_loaded_redshifts)
      error("Reading invalid number of tables along z axis.");
#endif

    /* Open table for this redshift index */
    char fname[eagle_table_path_name_length + 12];
    sprintf(fname, "%sz_%1.3f.hdf5", cooling->cooling_table_path,
            cooling->Redshifts[z_index]);
    message("Reading cooling table 'z_%1.3f.hdf5'",
            cooling->Redshifts[z_index]);

    hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) error("unable to open file %s", fname);

    char set_name[64];

    /* read in cooling rates due to metals */
    for (int specs = 0; specs < eagle_cooling_N_metal; specs++) {

      sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
      hid_t dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
      herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, net_cooling_rate);
      if (status < 0) error("error reading metal cooling rate table");
      status = H5Dclose(dataset);
      if (status < 0) error("error closing cooling dataset");

      /* Transpose from order tables are stored in (temperature, nH)
       * to (metal species, redshift, nH, temperature) where fastest
       * varying index is on right. Tables contain cooling rates but we
       * want rate of change of internal energy, hence minus sign. */
      for (int i = 0; i < eagle_cooling_N_density; i++) {
        for (int j = 0; j < eagle_cooling_N_temperature; j++) {

          /* Index in the HDF5 table */
          const int hdf5_index = row_major_index_2d(
              j, i, eagle_cooling_N_temperature, eagle_cooling_N_density);

          /* Index in the internal table */
          const int internal_index = row_major_index_4d(
              specs, local_z_index, i, j, eagle_cooling_N_metal,
              eagle_cooling_N_loaded_redshifts, eagle_cooling_N_density,
              eagle_cooling_N_temperature);

          /* Change the sign and transpose */
          cooling->table.metal_heating[internal_index] =
              -net_cooling_rate[hdf5_index];
        }
      }
    }

    /* read in cooling rates due to H + He */
    strcpy(set_name, "/Metal_free/Net_Cooling");
    hid_t dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, he_net_cooling_rate);
    if (status < 0) error("error reading metal free cooling rate table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    /* read in Temperature */
    strcpy(set_name, "/Metal_free/Temperature/Temperature");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temperature);
    if (status < 0) error("error reading temperature table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    /* Read in H + He electron abundance */
    strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_electron_abundance);
    if (status < 0) error("error reading electron density table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    /* Transpose from order tables are stored in (helium fraction, temperature,
     * nH) to (redshift, nH, helium fraction, temperature) where fastest
     * varying index is on right. */
    for (int i = 0; i < eagle_cooling_N_He_frac; i++) {
      for (int j = 0; j < eagle_cooling_N_temperature; j++) {
        for (int k = 0; k < eagle_cooling_N_density; k++) {

          /* Index in the HDF5 table */
          const int hdf5_index = row_major_index_3d(
              i, j, k, eagle_cooling_N_He_frac, eagle_cooling_N_temperature,
              eagle_cooling_N_density);

          /* Index in the internal table */
          const int internal_index = row_major_index_4d(
              local_z_index, k, i, j, eagle_cooling_N_loaded_redshifts,
              eagle_cooling_N_density, eagle_cooling_N_He_frac,
              eagle_cooling_N_temperature);

          /* Change the sign and transpose */
          cooling->table.H_plus_He_heating[internal_index] =
              -he_net_cooling_rate[hdf5_index];

          /* Convert to log T and transpose */
          cooling->table.temperature[internal_index] =
              log10(temperature[hdf5_index]);

          /* Just transpose */
          cooling->table.H_plus_He_electron_abundance[internal_index] =
              he_electron_abundance[hdf5_index];
        }
      }
    }

    /* read in electron densities due to metals */
    strcpy(set_name, "/Solar/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     electron_abundance);
    if (status < 0) error("error reading solar electron density table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    /* Transpose from order tables are stored in (temperature, nH) to
     * (redshift, nH, temperature) where fastest varying index is on right. */
    for (int i = 0; i < eagle_cooling_N_temperature; i++) {
      for (int j = 0; j < eagle_cooling_N_density; j++) {

        /* Index in the HDF5 table */
        const int hdf5_index = row_major_index_2d(
            i, j, eagle_cooling_N_temperature, eagle_cooling_N_density);

        /* Index in the internal table */
        const int internal_index = row_major_index_3d(
            local_z_index, j, i, eagle_cooling_N_loaded_redshifts,
            eagle_cooling_N_density, eagle_cooling_N_temperature);

        /* Just transpose */
        cooling->table.electron_abundance[internal_index] =
            electron_abundance[hdf5_index];
      }
    }

    status = H5Fclose(file_id);
    if (status < 0) error("error closing file");
  }

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

#ifdef SWIFT_DEBUG_CHECKS
  message("Done reading in general cooling table");
#endif

#else
  error("Need HDF5 to read cooling tables");
#endif
}
