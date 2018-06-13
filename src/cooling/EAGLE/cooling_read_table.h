#ifndef SWIFT_COOLING_EAGLE_READ_TABLES_H
#define SWIFT_COOLING_EAGLE_READ_TABLES_H

/*
 * @brief Reads cooling header data into cooling data structure
 *
 * @param fname Filepath
 * @param cooling Cooling data structure
 */
INLINE void ReadCoolingHeader(char *fname,
                              struct cooling_function_data *cooling) {
  int i;

  hid_t tempfile_id, dataset, datatype;

  herr_t status;

  /* fill the constants */
  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (tempfile_id < 0) {
    error("[ReadCoolingHeader()]: unable to open file %s\n", fname);
  }

  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_nH);
  status = H5Dclose(dataset);

  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_He);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_SolarAbundances);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Elements);
  status = H5Dclose(dataset);

  /* allocate arrays for cooling table header */
  // allocate_header_arrays();

  cooling->Temp = malloc(cooling->N_Temp * sizeof(float));
  cooling->nH = malloc(cooling->N_Temp * sizeof(float));
  cooling->HeFrac = malloc(cooling->N_Temp * sizeof(float));
  cooling->SolarAbundances = malloc(cooling->N_Temp * sizeof(float));

  /* fill the arrays */
  dataset = H5Dopen(tempfile_id, "/Solar/Temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Solar/Hydrogen_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->nH);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Helium_mass_fraction_bins",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->HeFrac);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Solar_mass_fractions",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->SolarAbundances);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Temperature/Energy_density_bins",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Therm);
  status = H5Dclose(dataset);

  char element_names[cooling->N_Elements][element_name_length];
  hsize_t string_length = element_name_length;

  /* names of chemical elements stored in table */
  datatype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(datatype, string_length);
  dataset = H5Dopen(tempfile_id, "/Header/Metal_names", H5P_DEFAULT);
  status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_Elements; i++)
    cooling->ElementNames[i] = mystrdup(element_names[i]);

  // char solar_abund_names[cooling_N_SolarAbundances][EL_NAME_LENGTH];

  /* assumed solar abundances used in constructing the tables, and corresponding
   * names */
  // dataset = H5Dopen(tempfile_id, "/Header/Abundances/Abund_names",
  // H5P_DEFAULT);
  // status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
  //                 solar_abund_names);
  // status = H5Dclose(dataset);
  // H5Tclose(datatype);

  // for (i = 0; i < cooling_N_SolarAbundances; i++)
  //  cooling_SolarAbundanceNames[i] = mystrdup(solar_abund_names[i]);

  // status = H5Fclose(tempfile_id);

  /* Convert to temperature, density and internal energy arrays to log10 */
  for (i = 0; i < cooling->N_Temp; i++) {
    cooling->Temp[i] = log10(cooling->Temp[i]);
    cooling->Therm[i] = log10(cooling->Therm[i]);
  }

  for (i = 0; i < cooling->N_nH; i++) cooling->nH[i] = log10(cooling->nH[i]);

  printf("Done with cooling table header.\n");
  fflush(stdout);
}

#endif /* SWIFT_COOLING_EAGLE_READ_TABLES_H */
