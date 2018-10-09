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
 * @file src/cooling/EAGLE/eagle_cool_tables.c
 * @brief Functions to read EAGLE tables
 */

/* Config parameters. */
#include "../config.h"

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cooling_struct.h"
#include "eagle_cool_tables.h"
#include "error.h"
#include "interpolate.h"

/* Names of the elements in the order they are stored in the files */
static const char *eagle_tables_element_names[9] = {
    "Carbon",  "Nitrogen", "Oxygen",  "Neon", "Magnesium",
    "Silicon", "Sulphur",  "Calcium", "Iron"};

/*
 * @brief Constructs the data structure containting the relevant cooling tables
 * for the redshift index (set in cooling_update)
 *
 * @param cooling Cooling data structure
 */
struct cooling_tables eagle_readtable(
    struct cooling_function_data *restrict cooling) {

  struct cooling_tables table;
  if (cooling->z_index < 0) {
    // z_index is set to < 0 in cooling_update if need
    // to read any of the high redshift tables
    table = get_redshift_invariant_table(cooling);
  } else {
    table = get_cooling_table(cooling);
  }

  return table;
}

/**
 * @brief Checks the tables that are currently loaded in memory and read
 * new ones if necessary.
 *
 * @param cooling The #cooling_function_data we play with.
 * @param index_z The index along the redshift axis of the tables of the current
 * z.
 */
void eagle_check_cooling_tables(struct cooling_function_data *restrict cooling,
                                int index_z) {

  /* Do we already have the right table in memory? */
  if (cooling->low_z_index == index_z) return;

  /* Record the table indices */
  cooling->low_z_index = index_z;
  cooling->high_z_index = index_z + 1;

  /* Load the damn thing */
  cooling->table = eagle_readtable(cooling);
}

/*
 * @brief Reads in EAGLE table of redshift values
 *
 * @param cooling Cooling data structure
 */
void GetCoolingRedshifts(struct cooling_function_data *cooling) {
  FILE *infile;

  int i = 0;

  char buffer[500], redfilename[516];

  sprintf(redfilename, "%s/redshifts.dat", cooling->cooling_table_path);
  infile = fopen(redfilename, "r");
  if (infile == NULL) puts("GetCoolingRedshifts can't open a file");

  if (fscanf(infile, "%s", buffer) != EOF) {
    cooling->N_Redshifts = atoi(buffer);
    if (posix_memalign((void **)&cooling->Redshifts, SWIFT_STRUCT_ALIGNMENT, cooling->N_Redshifts*sizeof(float)) != 0)
      error("Failed to allocate redshift table");

    while (fscanf(infile, "%s", buffer) != EOF) {
      cooling->Redshifts[i] = atof(buffer);
      i += 1;
    }
  }
  fclose(infile);

  /* EAGLE cooling assumes cooling->Redshifts table is in increasing order. Test
   * this. */
  for (i = 0; i < cooling->N_Redshifts - 2; i++)
    if (cooling->Redshifts[i + 1] < cooling->Redshifts[i]) {
      error("table should be in increasing order\n");
    }
}

/*
 * @brief Reads in EAGLE cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, helium fraction
 * solar element abundances, and elements used to index the cooling tables.
 *
 * @param fname Filepath for cooling table from which to read header
 * @param cooling Cooling data structure
 */
void ReadCoolingHeader(char *fname, struct cooling_function_data *cooling) {
#ifdef HAVE_HDF5
  int i;

  hid_t tempfile_id, dataset;

  herr_t status;

  // read sizes of array dimensions
  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0) error("unable to open file %s\n", fname);

  // read size of each table of values
  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Temp);
  if (status < 0) error("error reading number of temperature bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_nH);
  if (status < 0) error("error reading number of density bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_He);
  if (status < 0) error("error reading number of He fraction bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_SolarAbundances);
  if (status < 0) error("error reading number of solar abundance bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Elements);
  if (status < 0) error("error reading number of metal bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  // allocate arrays of values for each of the above quantities
  if (posix_memalign((void **)&cooling->Temp, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * sizeof(float)) !=0) 
    error("Failed to allocate temperature table");
  if (posix_memalign((void **)&cooling->Therm, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * sizeof(float)) !=0)
    error("Failed to allocate internal energy table");
  if (posix_memalign((void **)&cooling->nH, SWIFT_STRUCT_ALIGNMENT, cooling->N_nH * sizeof(float)) !=0) 
    error("Failed to allocate nH table");
  if (posix_memalign((void **)&cooling->HeFrac, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * sizeof(float)) !=0)
    error("Failed to allocate HeFrac table");
  if (posix_memalign((void **)&cooling->SolarAbundances, SWIFT_STRUCT_ALIGNMENT, cooling->N_SolarAbundances * sizeof(float)) !=0)
    error("Failed to allocate Solar abundances table");

  // read in values for each of the arrays
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

  // Convert to temperature, density and internal energy arrays to log10
  for (i = 0; i < cooling->N_Temp; i++) {
    cooling->Temp[i] = log10(cooling->Temp[i]);
    cooling->Therm[i] = log10(cooling->Therm[i]);
  }

  for (i = 0; i < cooling->N_nH; i++) cooling->nH[i] = log10(cooling->nH[i]);
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
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
 * @param cooling Cooling data structure
 */
struct cooling_tables get_redshift_invariant_table(
    struct cooling_function_data *restrict cooling) {
#ifdef HAVE_HDF5

  struct cooling_tables cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[521], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate = NULL;
  float *electron_abundance = NULL;
  float *temperature = NULL;
  float *he_net_cooling_rate = NULL;
  float *he_electron_abundance = NULL;

  // Allocate arrays for reading in cooling tables. 
  if (posix_memalign((void **)&net_cooling_rate, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate net_cooling_rate array");
  if (posix_memalign((void **)&electron_abundance, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate electron_abundance array");
  if (posix_memalign((void **)&temperature, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate temperature array");
  if (posix_memalign((void **)&he_net_cooling_rate, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate he_net_cooling_rate array");
  if (posix_memalign((void **)&he_electron_abundance, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate he_electron_abundance array");

  // Allocate arrays to store cooling tables. 
  if (posix_memalign((void **)&cooling_table.metal_heating, SWIFT_STRUCT_ALIGNMENT, cooling->N_Elements * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate metal_heating array");
  if (posix_memalign((void **)&cooling_table.electron_abundance, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate electron_abundance array");
  if (posix_memalign((void **)&cooling_table.temperature, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate temperature array");
  if (posix_memalign((void **)&cooling_table.H_plus_He_heating, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate H_plus_He_heating array");
  if (posix_memalign((void **)&cooling_table.H_plus_He_electron_abundance, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate H_plus_He_electron_abundance array");

  // Decide which high redshift table to read. Indices set in cooling_update
  if (cooling->low_z_index == -1) {
    sprintf(fname, "%sz_8.989nocompton.hdf5", cooling->cooling_table_path);
  } else if (cooling->low_z_index == -2) {
    sprintf(fname, "%sz_photodis.hdf5", cooling->cooling_table_path);
  }

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  // read in cooling rates due to metals
  for (specs = 0; specs < cooling->N_Elements; specs++) {
    sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     net_cooling_rate);
    if (status < 0) error("error reading metal cooling rate table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    // Transpose from order tables are stored in (temperature, nH)
    // to (nH, temperature, metal species) where fastest 
    // varying index is on right. Tables contain cooling rates but we 
    // want rate of change of internal energy, hence minus sign.
    for (j = 0; j < cooling->N_Temp; j++) {
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_2d(j, k, cooling->N_Temp, cooling->N_nH);
        cooling_index = row_major_index_3d(
            k, j, specs, cooling->N_nH, cooling->N_Temp, cooling->N_Elements);
        cooling_table.metal_heating[cooling_index] =
            -net_cooling_rate[table_index];
      }
    }
  }

  // read in cooling rates due to hydrogen and helium, H + He electron
  // abundances, temperatures
  strcpy(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_net_cooling_rate);
  if (status < 0) error("error reading metal free cooling rate table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  strcpy(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  if (status < 0) error("error reading temperature table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  if (status < 0) error("error reading electron density table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  // Transpose from order tables are stored in (helium fraction, temperature,
  // nH) to (nH, helium fraction, temperature) where fastest 
  // varying index is on right. Tables contain cooling rates but we 
  // want rate of change of internal energy, hence minus sign.
  for (i = 0; i < cooling->N_He; i++) {
    for (j = 0; j < cooling->N_Temp; j++) {
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_3d(i, j, k, cooling->N_He,
                                         cooling->N_Temp, cooling->N_nH);
        cooling_index = row_major_index_3d(k, i, j, cooling->N_nH,
                                           cooling->N_He, cooling->N_Temp);
        cooling_table.H_plus_He_heating[cooling_index] =
            -he_net_cooling_rate[table_index];
        cooling_table.H_plus_He_electron_abundance[cooling_index] =
            he_electron_abundance[table_index];
        cooling_table.temperature[cooling_index] =
            log10(temperature[table_index]);
      }
    }
  }

  // read in electron densities due to metals
  strcpy(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  if (status < 0) error("error reading solar electron density table");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  // Transpose from order tables are stored in (temperature, nH) to 
  // (nH, temperature) where fastest varying index is on right.
  for (i = 0; i < cooling->N_Temp; i++) {
    for (j = 0; j < cooling->N_nH; j++) {
      table_index = row_major_index_2d(i, j, cooling->N_Temp, cooling->N_nH);
      cooling_index = row_major_index_2d(j, i, cooling->N_nH, cooling->N_Temp);
      cooling_table.electron_abundance[cooling_index] =
          electron_abundance[table_index];
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

  return cooling_table;
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
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
 * @param cooling Cooling data structure
 */

struct cooling_tables get_cooling_table(
    struct cooling_function_data *restrict cooling) {
#ifdef HAVE_HDF5

  struct cooling_tables cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[1024], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate = NULL;
  float *electron_abundance = NULL;
  float *temperature = NULL;
  float *he_net_cooling_rate = NULL;
  float *he_electron_abundance = NULL;

  // Allocate arrays for reading in cooling tables. 
  if (posix_memalign((void **)&net_cooling_rate, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate net_cooling_rate array");
  if (posix_memalign((void **)&electron_abundance, SWIFT_STRUCT_ALIGNMENT, cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate electron_abundance array");
  if (posix_memalign((void **)&temperature, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate temperature array");
  if (posix_memalign((void **)&he_net_cooling_rate, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate he_net_cooling_rate array");
  if (posix_memalign((void **)&he_electron_abundance, SWIFT_STRUCT_ALIGNMENT, cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate he_electron_abundance array");

  // Allocate arrays to store cooling tables. Arrays contain two tables of 
  // cooling rates with one table being for the redshift above current redshift and one below.
  if (posix_memalign((void **)&cooling_table.metal_heating, SWIFT_STRUCT_ALIGNMENT, 2 * cooling->N_Elements * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate metal_heating array");
  if (posix_memalign((void **)&cooling_table.electron_abundance, SWIFT_STRUCT_ALIGNMENT, 2 * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate electron_abundance array");
  if (posix_memalign((void **)&cooling_table.temperature, SWIFT_STRUCT_ALIGNMENT, 2 * cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate temperature array");
  if (posix_memalign((void **)&cooling_table.H_plus_He_heating, SWIFT_STRUCT_ALIGNMENT, 2 * cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate H_plus_He_heating array");
  if (posix_memalign((void **)&cooling_table.H_plus_He_electron_abundance, SWIFT_STRUCT_ALIGNMENT, 2 * cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float)) !=0)
    error("Failed to allocate H_plus_He_electron_abundance array");

  // Read in tables, transpose so that values for indices which vary most are
  // adjacent. Repeat for redshift above and redshift below current value. 
  for (int z_index = cooling->low_z_index; z_index <= cooling->high_z_index;
       z_index++) {
    sprintf(fname, "%sz_%1.3f.hdf5", cooling->cooling_table_path,
            cooling->Redshifts[z_index]);
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) error("unable to open file %s", fname);

    // read in cooling rates due to metals
    for (specs = 0; specs < cooling->N_Elements; specs++) {
      sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
      dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       net_cooling_rate);
      if (status < 0) error("error reading metal cooling rate table");
      status = H5Dclose(dataset);
      if (status < 0) error("error closing cooling dataset");

      // Transpose from order tables are stored in (temperature, nH)
      // to (redshift, nH, temperature, metal species) where fastest 
      // varying index is on right. Tables contain cooling rates but we 
      // want rate of change of internal energy, hence minus sign.
      for (i = 0; i < cooling->N_nH; i++) {
        for (j = 0; j < cooling->N_Temp; j++) {
          table_index =
              row_major_index_2d(j, i, cooling->N_Temp, cooling->N_nH);
          cooling_index = row_major_index_4d(
              z_index - cooling->low_z_index, i, j, specs, 2, cooling->N_nH,
              cooling->N_Temp, cooling->N_Elements);
          cooling_table.metal_heating[cooling_index] =
              -net_cooling_rate[table_index];
        }
      }
    }

    // read in cooling rates due to hydrogen and helium, H + He electron
    // abundances, temperatures
    strcpy(set_name, "/Metal_free/Net_Cooling");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_net_cooling_rate);
    if (status < 0) error("error reading metal free cooling rate table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    strcpy(set_name, "/Metal_free/Temperature/Temperature");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temperature);
    if (status < 0) error("error reading temperature table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_electron_abundance);
    if (status < 0) error("error reading electron density table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    // Transpose from order tables are stored in (helium fraction, temperature,
    // nH) to (redshift, nH, helium fraction, temperature) where fastest 
    // varying index is on right.
    for (i = 0; i < cooling->N_He; i++) {
      for (j = 0; j < cooling->N_Temp; j++) {
        for (k = 0; k < cooling->N_nH; k++) {
          table_index = row_major_index_3d(i, j, k, cooling->N_He,
                                           cooling->N_Temp, cooling->N_nH);
          cooling_index =
              row_major_index_4d(z_index - cooling->low_z_index, k, i, j, 2,
                                 cooling->N_nH, cooling->N_He, cooling->N_Temp);
          cooling_table.H_plus_He_heating[cooling_index] =
              -he_net_cooling_rate[table_index];
          cooling_table.H_plus_He_electron_abundance[cooling_index] =
              he_electron_abundance[table_index];
          cooling_table.temperature[cooling_index] =
              log10(temperature[table_index]);
        }
      }
    }

    // read in electron densities due to metals
    strcpy(set_name, "/Solar/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     electron_abundance);
    if (status < 0) error("error reading solar electron density table");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing cooling dataset");

    // Transpose from order tables are stored in (temperature, nH) to 
    // (redshift, nH, temperature) where fastest varying index is on right.
    for (i = 0; i < cooling->N_Temp; i++) {
      for (j = 0; j < cooling->N_nH; j++) {
        table_index = row_major_index_2d(i, j, cooling->N_Temp, cooling->N_nH);
        cooling_index = row_major_index_3d(z_index - cooling->low_z_index, j, i,
                                           2, cooling->N_nH, cooling->N_Temp);
        cooling_table.electron_abundance[cooling_index] =
            electron_abundance[table_index];
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
  message("done reading in general cooling table");
#endif

  return cooling_table;
#else
  error("Need HDF5 to read cooling tables");
#endif
}
