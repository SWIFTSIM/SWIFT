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

inline void read_yield_tables(struct stars_props *restrict stars){
#ifdef HAVE_HDF5

  int i, j, k;

  char fname[256], setname[100];
  
  // move these definitions somewhere else...
  stars->yields_SNIa_n_elements = 42;
  stars->yields_SNII_n_mass = 11;
  stars->yields_SNII_n_elements = 11;
  stars->yields_SNII_n_z = 5;
  stars->yields_AGB_n_mass = 23;
  stars->yields_AGB_n_elements = 11;
  stars->yields_AGB_n_z = 3;
  stars->lifetimes.n_mass = 30;
  stars->lifetimes.n_z = 6;

  hid_t file_id, dataset, datatype;
  herr_t status;

  /* Read SNIa tables */
  sprintf(fname, "%s/SNIa.hdf5", stars->yield_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) error("unable to open file %s\n", fname);

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->SNIa_element_names));
  if (status < 0) error("error reading SNIa element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  // What is this for? Copied from EAGLE
  //for (i = 0; i < stars->yields_SNIa_n_elements; i++)
  //  yieldsSNIa.ElementName[i] = mystrdup(yieldsSNIa.ElementName[i]);

  dataset = H5Dopen(file_id, "Yield", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->yields_SNIa));
  if (status < 0) error("error reading SNIa yields");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Total_Metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->yield_SNIa_total_metals_SPH));
  if (status < 0) error("error reading SNIa total metal");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  status = H5Fclose(file_id);
  if (status < 0) error("error closing SNIa file");

  /* Read SNII tables */

  sprintf(fname, "%s/SNII.hdf5", stars->yield_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) error("unable to open file %s\n", fname);

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->SNII_element_names));
  if (status < 0) error("error reading SNII element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  // What is this for (again)? copied from EAGLE
  //for (i = 0; i < stars->yields_SNII_n_elements; i++)
  //  yieldsSNII.ElementName[i] = mystrdup(yieldsSNII.ElementName[i]);

  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->yield_SNII.mass));
  if (status < 0) error("error reading SNII masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->yield_SNII.metallicity));
  if (status < 0) error("error reading SNII metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  float tempyield1[stars->yields_SNII_n_elements][stars->yields_SNII_n_mass];

  float tempej1[stars->yields_SNII_n_mass], tempmet1[stars->yields_SNII_n_mass];

  char *tempname1[stars->yields_SNII_n_z];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname1);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  for (i = 0; i < stars->yields_SNII_n_z; i++) {
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

    for (k = 0; k < stars->yields_SNII_n_mass; k++) {
      stars->yield_SNII.ejecta[i][k] = tempej1[k];
      stars->yield_SNII.total_metals[i][k] = tempmet1[k];

      for (j = 0; j < stars->yields_SNII_n_elements; j++)
        stars->yield_SNII.yield[i][j][k] = tempyield1[j][k];
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* Read AGB tables */

  sprintf(fname, "%s/AGB.hdf5", stars->yield_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) error("unable to open file %s\n", fname);

  /* allocate element name array */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Species_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->AGB_element_names));
  if (status < 0) error("error reading AGB element names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  // What is this for (again)? copied from EAGLE
  //for (i = 0; i < stars->yields_AGB_n_elements; i++)
  //  yieldsAGB.ElementName[i] = mystrdup(yieldsAGB.ElementName[i]);

  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->yield_AGB.mass));
  if (status < 0) error("error reading AGB masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->yield_AGB.metallicity));
  if (status < 0) error("error reading AGB metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  float tempyield2[stars->yields_AGB_n_elements][stars->yields_AGB_n_mass];

  float tempej2[stars->yields_AGB_n_mass], tempmet2[stars->yields_AGB_n_mass];

  char *tempname2[stars->yields_AGB_n_z];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname2);
  if (status < 0) error("error reading yield table names");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");

  for (i = 0; i < stars->yields_AGB_n_z; i++) {
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

    for (k = 0; k < stars->yields_AGB_n_mass; k++) {
      stars->yield_AGB.ejecta[i][k] = tempej2[k];
      stars->yield_AGB.total_metals[i][k] = tempmet2[k];

      for (j = 0; j < stars->yields_AGB_n_elements; j++)
        stars->yield_AGB.yield[i][j][k] = tempyield2[j][k];
    }
  }

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

  /* Read lifetimes table */

  sprintf(fname, "%s/AGB.hdf5", stars->yield_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) error("unable to open file %s\n", fname);

  dataset = H5Dopen(file_id, "Masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->lifetimes.mass));
  if (status < 0) error("error reading lifetime table masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(stars->lifetimes.metallicity));
  if (status < 0) error("error reading AGB metallicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  float temptime[stars->lifetimes.n_z][stars->lifetimes.n_mass];

  dataset = H5Dopen(file_id, "Lifetimes", H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temptime);
  H5Dclose(dataset);

  for (i = 0; i < stars->lifetimes.n_z; i++) {
    for (j = 0; j < stars->lifetimes.n_mass; j++) {
      stars->lifetimes.dyingtime[i][j] = log10(temptime[i][j]);
    }
  }

  H5Fclose(file_id);


#endif
}

#endif
