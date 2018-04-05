#ifndef SWIFT_EAGLE_COOL_TABLES_H
#define SWIFT_EAGLE_COOL_TABLES_H

#include "cooling_struct.h"
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>
#include "error.h"

inline int whatisopen(hid_t fid) {
        ssize_t cnt;
        int howmany;
        int i;
        H5I_type_t ot;
        hid_t anobj;
        hid_t *objs;
        char name[1024];
        herr_t status;

        cnt = H5Fget_obj_count(fid, H5F_OBJ_ALL);

        if (cnt <= 0) return cnt;

        printf("%d object(s) open\n", (int) cnt);

        objs = malloc(cnt * sizeof(hid_t));

        howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, cnt, objs);

        printf("open objects:\n");

        for (i = 0; i < howmany; i++ ) {
             anobj = *objs++;
             ot = H5Iget_type(anobj);
             status = H5Iget_name(anobj, name, 1024);
             printf(" %d: type %d, name %s\n",i,ot,name);
        }

        return howmany;
}


inline char *mystrdup(const char *s) {
  char *p;

  p = (char *)malloc((strlen(s) + 1)*sizeof(char));
  strcpy(p, s);
  return p;
}

int row_major_index_2d(int, int, int, int);

int row_major_index_3d(int, int, int, int, int, int);

int row_major_index_4d(int, int, int, int, int, int, int, int);

inline void GetCoolingRedshifts(struct cooling_function_data *cooling) {
  FILE *infile;

  int i = 0;

  char buffer[500], redfilename[500];

  sprintf(redfilename, "%s/redshifts.dat", cooling->cooling_table_path);
  infile = fopen(redfilename, "r");
  if (infile == NULL) puts("GetCoolingRedshifts can't open a file");

  if (fscanf(infile, "%s", buffer) != EOF) {
    cooling->N_Redshifts = atoi(buffer);
    cooling->Redshifts =
        (float *)malloc(cooling->N_Redshifts * sizeof(float));

    while (fscanf(infile, "%s", buffer) != EOF) {
      cooling->Redshifts[i] = atof(buffer);
      i += 1;
    }
  }
  fclose(infile);

  //printf("eagle_cool_tables.h redshift max, min, N_Redshifts: %.5e, %.5e, %d\n", cooling->Redshifts[cooling->N_Redshifts-1], cooling->Redshifts[0], cooling->N_Redshifts);

}

inline void ReadCoolingHeader(char *fname, struct cooling_function_data *cooling) {
  int i;

  hid_t tempfile_id, dataset, datatype;

  herr_t status;

  /* fill the constants */
  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (tempfile_id < 0) {
    error("[ReadCoolingHeader()]: unable to open file %s\n", fname);
  }

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_nH);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_He);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_SolarAbundances);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Elements);
  status = H5Dclose(dataset);

  /* allocate arrays for cooling table header */
  //allocate_header_arrays();

  cooling->Temp = malloc(cooling->N_Temp*sizeof(float));
  cooling->nH = malloc(cooling->N_nH*sizeof(float));
  cooling->HeFrac = malloc(cooling->N_He*sizeof(float));
  cooling->SolarAbundances = malloc(cooling->N_SolarAbundances*sizeof(float));
  cooling->Therm = malloc(cooling->N_Temp*sizeof(float));
  cooling->ElementNames = malloc(cooling->N_Elements*eagle_element_name_length*sizeof(char));
  cooling->SolarAbundanceNames = malloc(cooling->N_SolarAbundances*eagle_element_name_length*sizeof(char));
  
  /* fill the arrays */
  dataset = H5Dopen(tempfile_id, "/Solar/Temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Solar/Hydrogen_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->nH);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Helium_mass_fraction_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->HeFrac);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Solar_mass_fractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->SolarAbundances);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Temperature/Energy_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Therm);
  status = H5Dclose(dataset);

  char element_names[cooling->N_Elements][eagle_element_name_length];
  hsize_t string_length = eagle_element_name_length;

  /* names of chemical elements stored in table */
  datatype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(datatype, string_length);
  dataset = H5Dopen(tempfile_id, "/Header/Metal_names", H5P_DEFAULT);
  status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  status = H5Dclose(dataset);
  H5Tclose(datatype);
  
  for (i = 0; i < cooling->N_Elements; i++)
    cooling->ElementNames[i] = mystrdup(element_names[i]);
  
  char solar_abund_names[cooling->N_SolarAbundances][eagle_element_name_length];

  /* assumed solar abundances used in constructing the tables, and corresponding
   * names */
  datatype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(datatype, string_length);
  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Abund_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   solar_abund_names);
  status = H5Dclose(dataset);
  H5Tclose(datatype);

  for (i = 0; i < cooling->N_SolarAbundances; i++)
    cooling->SolarAbundanceNames[i] = mystrdup(solar_abund_names[i]);

  status = H5Fclose(tempfile_id);

  /* Convert to temperature, density and internal energy arrays to log10 */
  for (i = 0; i < cooling->N_Temp; i++) {
    cooling->Temp[i] = log10(cooling->Temp[i]);
    cooling->Therm[i] = log10(cooling->Therm[i]);
  }

  for (i = 0; i < cooling->N_nH; i++) cooling->nH[i] = log10(cooling->nH[i]);
  
  //printf("eagle_cooling_tables.h temp max, min, N_Temp: %.5e, %.5e, %d\n",cooling->Temp[cooling->N_Temp-1], cooling->Temp[0],cooling->N_Temp);
  //printf("eagle_cooling_tables.h internal energy max, min, N_Temp: %.5e, %.5e, %d\n",cooling->Therm[cooling->N_Temp-1], cooling->Therm[0],cooling->N_Temp);
  //printf("eagle_cooling_tables.h H max, min, N_nH: %.5e, %.5e, %d\n",cooling->nH[cooling->N_nH-1], cooling->nH[0],cooling->N_nH);
  //printf("eagle_cooling_tables.h He max, min, N_He: %.5e, %.5e, %d\n",cooling->HeFrac[cooling->N_He-1], cooling->HeFrac[0],cooling->N_He);
  //printf("eagle_cooling_tables.h Solar abundances max, min, N_SolarAbundances: %.5e, %.5e, %d\n",cooling->SolarAbundances[cooling->N_SolarAbundances-1], cooling->SolarAbundances[0],cooling->N_SolarAbundances);
  //printf("eagle_cooling_tables.h N_Elements: %d\n",cooling->N_Elements);


  printf("Done with cooling table header.\n");
  fflush(stdout);

}


/*
 * ----------------------------------------------------------------------
 * Get the cooling table for photoionized cooling (before redshift ~9)
 * ----------------------------------------------------------------------
 */

inline struct cooling_tables_redshift_invariant get_no_compt_table(char *cooling_table_path, const struct cooling_function_data *restrict cooling) {

  struct cooling_tables_redshift_invariant cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  net_cooling_rate = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  electron_abundance = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  temperature = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  he_net_cooling_rate = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  he_electron_abundance = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  
  cooling_table.metal_heating = (float *)malloc(cooling->N_Elements*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.electron_abundance = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.temperature = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.H_plus_He_heating = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.H_plus_He_electron_abundance = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));

  sprintf(fname, "%sz_8.989nocompton.hdf5", cooling_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  //printf("GetNoCompTable Redshift 1 %ld %s\n", (long int)file_id, fname);
  //fflush(stdout);

  /* For normal elements */
  for (specs = 0; specs < cooling->N_Elements; specs++) {
    sprintf(set_name, "/%s/Net_Cooling", cooling->ElementNames[specs]);
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     net_cooling_rate);
    status = H5Dclose(dataset);
    //for(i = 0; i < cooling->N_Temp*cooling->N_nH; i++) printf("eagle_cool_tables.h i, net_cooling_rate %d, %.5e\n",i,net_cooling_rate[i]);

    for (j = 0; j < cooling->N_Temp; j++){
      for (k = 0; k < cooling->N_nH; k++){
        table_index = row_major_index_2d(j,k,cooling->N_Temp,cooling->N_nH);
        cooling_index = row_major_index_4d(0,specs,k,j,1,cooling->N_Elements,cooling->N_nH,cooling->N_Temp); //Redshift invariant table!!!
        cooling_table.metal_heating[cooling_index] = -net_cooling_rate[table_index];
	//printf("eagle_cool_tables.h j,k,j*cooling->N_nH+k,table_index,cooling_index,net cooling rate, cooling_table value %d, %d, %d, %d, %d %.5e, %.5e\n",j,k,j*cooling->N_nH+k,table_index,cooling_index,-net_cooling_rate[table_index],cooling_table.metal_heating[cooling_index]);
      }
    }
  }
  //for (i = 0; i < cooling->N_nH; i++){
  //  table_index = row_major_index_2d(0,i,cooling->N_Temp,cooling->N_nH);
  //  cooling_index = row_major_index_4d(0,0,i,0,1,cooling->N_Elements,cooling->N_nH,cooling->N_Temp);
  //  printf("eagle_cool_tables.h i, cooling_index, table_index, cooling_table.metal_heating, net heating = %d %d %d %.5e %.5e\n",i,cooling_index,table_index,cooling_table.metal_heating[cooling_index],net_cooling_rate[table_index]);
  //}
  //for (i = 0; i < cooling->N_Elements*cooling->N_Temp*cooling->N_nH; i++) printf("eagle cooling.h i, cooling_table %d, %.5e\n",i,cooling_table.metal_heating[i]);
  //error("stop in eagle_cool_tables.h");


  /* Helium */
  sprintf(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_net_cooling_rate);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_He; i++){
    for (j = 0; j < cooling->N_Temp; j++){
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_3d(i,j,k,cooling->N_He,cooling->N_Temp,cooling->N_nH);
        cooling_index = row_major_index_4d(0,i,k,j,1,cooling->N_He,cooling->N_nH,cooling->N_Temp); //Redshift invariant table!!!
        cooling_table.H_plus_He_heating[cooling_index] = -he_net_cooling_rate[table_index];
        cooling_table.H_plus_He_electron_abundance[cooling_index] =
            he_electron_abundance[table_index];
        cooling_table.temperature[cooling_index] = log10(temperature[table_index]);
      }
    }
  }

  sprintf(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_Temp; i++){
    for (j = 0; j < cooling->N_nH; j++){
      table_index = row_major_index_2d(i,j,cooling->N_Temp,cooling->N_nH);
      cooling_index = row_major_index_3d(0,j,i,1,cooling->N_nH,cooling->N_Temp); //Redshift invariant table!!!
      cooling_table.electron_abundance[cooling_index] = electron_abundance[table_index];
    }
  }

  status = H5Fclose(file_id);

  //cooling_table.metals_heating = cooling_MetalsNetHeating;
  //cooling_table.H_plus_He_heating = cooling_HplusHeNetHeating;
  //cooling_table.H_plus_He_electron_abundance = cooling_HplusHeElectronAbundance;
  //cooling_table.temperature = cooling_ThermalToTemp;
  //cooling_table.electron_abundance = cooling_SolarElectronAbundance;

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

  printf("eagle_cool_tables.h done reading in no compton table\n");

  return cooling_table;
}

/*
 * ----------------------------------------------------------------------
 * Get the cooling table for collisional cooling (before reionisation)
 * ----------------------------------------------------------------------
 */

inline struct cooling_tables_redshift_invariant get_collisional_table(char *cooling_table_path, const struct cooling_function_data* restrict cooling) {

  struct cooling_tables_redshift_invariant cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  net_cooling_rate = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  electron_abundance = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  temperature = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  he_net_cooling_rate = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  he_electron_abundance = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));

  cooling_table.metal_heating = (float *)malloc(cooling->N_Elements*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.electron_abundance = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.temperature = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.H_plus_He_heating = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.H_plus_He_electron_abundance = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));

  sprintf(fname, "%sz_photodis.hdf5", cooling_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* For normal elements */
  for (specs = 0; specs < cooling->N_Elements; specs++) {
    sprintf(set_name, "/%s/Net_Cooling", cooling->ElementNames[specs]);
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     net_cooling_rate);
    status = H5Dclose(dataset);

    for (j = 0; j < cooling->N_Temp; j++){
      for (k = 0; k < cooling->N_nH; k++){
        table_index = row_major_index_2d(j,k,cooling->N_Temp,cooling->N_nH);
        cooling_index = row_major_index_4d(0,specs,k,j,1,cooling->N_Elements,cooling->N_nH,cooling->N_Temp); //Redshift invariant table!!!
        cooling_table.metal_heating[cooling_index] = -net_cooling_rate[table_index];
      }
    }
  }

  /* Helium */
  sprintf(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_net_cooling_rate);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  status = H5Dclose(dataset);

  sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_He; i++){
    for (j = 0; j < cooling->N_Temp; j++){
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_3d(i,j,k,cooling->N_He,cooling->N_Temp,cooling->N_nH);
        cooling_index = row_major_index_4d(0,i,k,j,1,cooling->N_He,cooling->N_nH,cooling->N_Temp); //Redshift invariant table!!!
        cooling_table.H_plus_He_heating[cooling_index] = -he_net_cooling_rate[table_index];
        cooling_table.H_plus_He_electron_abundance[cooling_index] =
            he_electron_abundance[table_index];
        cooling_table.temperature[cooling_index] = log10(temperature[table_index]);
      }
    }
  }

  sprintf(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_Temp; i++){
    for (j = 0; j < cooling->N_nH; j++){
      table_index = row_major_index_2d(i,j,cooling->N_Temp,cooling->N_nH);
      cooling_index = row_major_index_3d(0,j,i,1,cooling->N_nH,cooling->N_Temp); //Redshift invariant table!!!
      cooling_table.electron_abundance[cooling_index] = electron_abundance[table_index];
    }
  }

  status = H5Fclose(file_id);

  //cooling_table.metal_heating = cooling_MetalsNetHeating;
  //cooling_table.H_plus_He_heating = cooling_HplusHeNetHeating;
  //cooling_table.H_plus_He_electron_abundance = cooling_HplusHeElectronAbundance;
  //cooling_table.temperature = cooling_ThermalToTemp;
  //cooling_table.electron_abundance = cooling_SolarElectronAbundance;

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

  printf("eagle_cool_tables.h done reading in collisional table\n");
  return cooling_table;
}

/*
 * ----------------------------------------------------------------------
 * Get the cooling tables that bound the given redshift
 * ----------------------------------------------------------------------
 */

inline struct cooling_tables get_cooling_table(char *cooling_table_path, const struct cooling_function_data* restrict cooling) {

  struct cooling_tables cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  net_cooling_rate = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  electron_abundance = (float *)malloc(cooling->N_Temp*cooling->N_nH*sizeof(float));
  temperature = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  he_net_cooling_rate = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  he_electron_abundance = (float *)malloc(cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  
  cooling_table.metal_heating = (float *)malloc(cooling->N_Redshifts*cooling->N_Elements*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.electron_abundance = (float *)malloc(cooling->N_Redshifts*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.temperature = (float *)malloc(cooling->N_Redshifts*cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.H_plus_He_heating = (float *)malloc(cooling->N_Redshifts*cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));
  cooling_table.H_plus_He_electron_abundance = (float *)malloc(cooling->N_Redshifts*cooling->N_He*cooling->N_Temp*cooling->N_nH*sizeof(float));

  /* For normal elements */
  for (int z_index = 0; z_index < cooling->N_Redshifts; z_index++){
    sprintf(fname, "%sz_%1.3f.hdf5", cooling_table_path,
            cooling->Redshifts[z_index]);
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0) {
      error("[GetCoolingTables()]: unable to open file %s\n", fname);
    }

    for (specs = 0; specs < cooling->N_Elements; specs++) {
      sprintf(set_name, "/%s/Net_Cooling", cooling->ElementNames[specs]);
      dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       net_cooling_rate);
      status = H5Dclose(dataset);

      for (i = 0; i < cooling->N_nH; i++){
        for (j = 0; j < cooling->N_Temp; j++){
          table_index = row_major_index_2d(j,i,cooling->N_Temp,cooling->N_nH);
          cooling_index = row_major_index_4d(z_index,specs,i,j,cooling->N_Redshifts,cooling->N_Elements,cooling->N_nH,cooling->N_Temp); 
          cooling_table.metal_heating[cooling_index] = -net_cooling_rate[table_index];
        }
      }
    }

    /* Helium */
    sprintf(set_name, "/Metal_free/Net_Cooling");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_net_cooling_rate);
    status = H5Dclose(dataset);

    sprintf(set_name, "/Metal_free/Temperature/Temperature");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temperature);
    status = H5Dclose(dataset);

    sprintf(set_name, "/Metal_free/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_electron_abundance);
    status = H5Dclose(dataset);

    for (i = 0; i < cooling->N_He; i++){
      for (j = 0; j < cooling->N_Temp; j++){
        for (k = 0; k < cooling->N_nH; k++) {
          table_index = row_major_index_3d(i,j,k,cooling->N_He,cooling->N_Temp,cooling->N_nH);
          cooling_index = row_major_index_4d(z_index,i,k,j,cooling->N_Redshifts,cooling->N_He,cooling->N_nH,cooling->N_Temp); 
          cooling_table.H_plus_He_heating[cooling_index] = -he_net_cooling_rate[table_index];
          cooling_table.H_plus_He_electron_abundance[cooling_index] =
              he_electron_abundance[table_index];
          cooling_table.temperature[cooling_index] = log10(temperature[table_index]);
        }
      }
    }

    sprintf(set_name, "/Solar/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     electron_abundance);
    status = H5Dclose(dataset);

    for (i = 0; i < cooling->N_Temp; i++){
      for (j = 0; j < cooling->N_nH; j++){
        table_index = row_major_index_2d(i,j,cooling->N_Temp,cooling->N_nH);
        cooling_index = row_major_index_3d(z_index,j,i,cooling->N_Redshifts,cooling->N_nH,cooling->N_Temp); 
        cooling_table.electron_abundance[cooling_index] = electron_abundance[table_index];
      }
    }

    status = H5Fclose(file_id);
  }

  //cooling_table.metal_heating = cooling_MetalsNetHeating;
  //cooling_table.H_plus_He_heating = cooling_HplusHeNetHeating;
  //cooling_table.H_plus_He_electron_abundance = cooling_HplusHeElectronAbundance;
  //cooling_table.temperature = cooling_ThermalToTemp;
  //cooling_table.electron_abundance = cooling_SolarElectronAbundance;
  
  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

  printf("eagle_cool_tables.h done reading in general cooling table\n");

  return cooling_table;
}

inline struct eagle_cooling_table eagle_readtable(char *cooling_table_path, const struct cooling_function_data* restrict cooling){

  struct eagle_cooling_table table;

  table.photoionisation_cooling = get_no_compt_table(cooling_table_path, cooling);
  table.collisional_cooling = get_collisional_table(cooling_table_path, cooling);
  table.element_cooling = get_cooling_table(cooling_table_path, cooling);

  return table;
}

inline int element_index(char *element_name, const struct cooling_function_data* restrict cooling) {
  int i;

  for (i = 0; i < cooling->N_Elements; i++)
    if (strcmp(cooling->ElementNames[i], element_name) == 0) return i;

  /* element not found */
  return -1;
}


inline int get_element_index(char *table[20], int size, char *element_name) {
  int i;

  for (i = 0; i < size; i++)
    if (strcmp(table[i], element_name) == 0) return i;

  /* element not found */
  return -1;
}

inline void MakeNamePointers(struct cooling_function_data* cooling) {
  int i, j, sili_index = 0;
  char ElementNames[cooling->N_Elements][eagle_element_name_length];

  /* This is ridiculous, way too many element name arrays. Needs to be changed */
  //ElementNames = malloc(cooling->N_Elements*eagle_element_name_length*sizeof(char));
  strcpy(ElementNames[0], "Hydrogen");
  strcpy(ElementNames[1], "Helium");
  strcpy(ElementNames[2], "Carbon");
  strcpy(ElementNames[3], "Nitrogen");
  strcpy(ElementNames[4], "Oxygen");
  strcpy(ElementNames[5], "Neon");
  strcpy(ElementNames[6], "Magnesium");
  strcpy(ElementNames[7], "Silicon");
  strcpy(ElementNames[8], "Iron");

  cooling->ElementNamePointers = malloc(cooling->N_Elements * sizeof(int));
  cooling->SolarAbundanceNamePointers = malloc(cooling->N_Elements * sizeof(int));

  for (i = 0; i < cooling->N_Elements; i++) {
    if (strcmp(ElementNames[i], "Silicon") == 0) sili_index = i;
  }

  for (i = 0; i < cooling->N_Elements; i++) {
    cooling->SolarAbundanceNamePointers[i] = -999;
    cooling->ElementNamePointers[i] = -999;

    for (j = 0; j < cooling->N_SolarAbundances; j++) {
      if (strcmp(cooling->ElementNames[i], cooling->SolarAbundanceNames[j]) == 0)
        cooling->SolarAbundanceNamePointers[i] = j;
    }

    if (strcmp(cooling->ElementNames[i], "Sulphur") == 0 ||
        strcmp(cooling->ElementNames[i], "Calcium") ==
            0) /* These elements are tracked! */
      cooling->ElementNamePointers[i] = -1 * sili_index;
    else {
      for (j = 0; j < cooling->N_Elements; j++) {
        if (strcmp(cooling->ElementNames[i], ElementNames[j]) == 0)
          cooling->ElementNamePointers[i] = j;
      }
    }
  }
}


//inline int set_cooling_SolarAbundances(const float *element_abundance,
//                                float *cooling_element_abundance, 
//				const struct cooling_function_data* restrict cooling) {
//  int i, index;
//
//  for(i = 0; i < cooling->N_Elements; i++){
//    index = get_element_index(cooling->SolarAbundanceNames,
//                          cooling->N_SolarAbundances, cooling->ElementNames[i]);
//    if (cooling->SolarAbundances[index] != 0.0) cooling_element_abundance[i] = element_abundance[i]/cooling->SolarAbundances[index];
//    else cooling_element_abundance[i] = 0.0;
//    printf ("eagle_cool_tables.h element, name, abundance, solar abundance, cooling abundance %d %s %.5e %.5e %.5e\n",index,cooling->ElementNames[i], element_abundance[i],cooling->SolarAbundances[index], cooling_element_abundance[i]);
//  }
//
//  return 0;
//}

//inline int set_cooling_SolarAbundances(const float *element_abundance,
//                                float *cooling_element_abundance,
//				const struct cooling_function_data* restrict cooling) {
//  int i, index;
//
//  int static Silicon_SPH_Index = -1;
//  int static Calcium_SPH_Index = -1;
//  int static Sulphur_SPH_Index = -1;
//
//  int static Silicon_CoolHeat_Index = -1;
//  int static Calcium_CoolHeat_Index = -1;
//  int static Sulphur_CoolHeat_Index = -1;
//
//  static int first_call = 0;
//
//  if (first_call == 0) {
//    /* determine (inverse of) solar abundance of these elements */
//    for (i = 0; i < cooling_N_Elements; i++) {
//      index =
//          get_element_index(cooling_SolarAbundanceNames,
//                            cooling_N_SolarAbundances, cooling_ElementNames[i]);
//
//      if (index < 0) endrun(-12345);
//
//      index = SolarAbundanceNamePointers[i];
//
//      cooling_ElementAbundance_SOLARM1[i] = 1. / cooling_SolarAbundances[index];
//
//      index = ElementNamePointers[i];
//
//      if (index < 0 && ThisTask == 0)
//        printf("[bg_cooling] element not found %s\n", cooling_ElementNames[i]);
//    }
//
//    /* Sulphur tracks Silicon: may choose not to follow Sulphur as SPH element
//     */
//    /* Same is true for Calcium */
//    /* We will assume the code tracks Silicon, and may need to scale Calcium and
//     * Sulphur accordingly */
//
//    Silicon_SPH_Index = element_index("Silicon");
//    Calcium_SPH_Index = element_index("Calcium");
//    Sulphur_SPH_Index = element_index("Sulphur");
//
//    Silicon_CoolHeat_Index =
//        get_element_index(cooling_ElementNames, cooling_N_Elements, "Silicon");
//    Calcium_CoolHeat_Index =
//        get_element_index(cooling_ElementNames, cooling_N_Elements, "Calcium");
//    Sulphur_CoolHeat_Index =
//        get_element_index(cooling_ElementNames, cooling_N_Elements, "Sulphur");
//
//    if (Silicon_CoolHeat_Index == -1 || Calcium_CoolHeat_Index == -1 ||
//        Sulphur_CoolHeat_Index == -1) {
//      if (ThisTask == 0)
//        printf("[bg_cooling] error: did not find Si or Ca or S??\n");
//      endrun(-1233);
//    }
//
//    first_call = 1;
//  }
//  for (i = 0; i < cooling_N_Elements; i++) {
//    if (i == Calcium_CoolHeat_Index && Calcium_SPH_Index == -1)
//      /* SPH does not track Calcium: use Si abundance */
//      if (Silicon_SPH_Index == -1)
//        cooling_element_abundance[i] = 0.0;
//      else
//        cooling_element_abundance[i] =
//            element_abundance[Silicon_SPH_Index] *
//            cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
//    else if (i == Sulphur_CoolHeat_Index && Sulphur_SPH_Index == -1)
//      /* SPH does not track Sulphur: use Si abundance */
//      if (Silicon_SPH_Index == -1)
//        cooling_element_abundance[i] = 0.0;
//      else
//        cooling_element_abundance[i] =
//            element_abundance[Silicon_SPH_Index] *
//            cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
//    else
//      cooling_element_abundance[i] = element_abundance[ElementNamePointers[i]] *
//                                     cooling_ElementAbundance_SOLARM1[i];
//  }
//
//  return 0;
//}


#endif
