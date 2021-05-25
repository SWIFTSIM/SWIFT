/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

/* Config parameters. */
#include "../config.h"
#include "parser.h"
 
/* Local headers */
#include "neutrino.h"

void renderer_init(struct swift_params *params, const struct unit_system *us,
                  const struct phys_const *phys_const,
                  const struct cosmology *c,
                  const struct neutrino_props *np, struct neutrino_renderer *rend) {

  /* Do we need to do anything? */
  if (!np->use_linear_response) return;

  /* Load the data */
  char filename[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "Neutrino:boltzmann_filename", filename);
  
  message("Loading perturbation file from %s", filename);
  
  /* Open the hdf5 file (file exists error handled by HDF5) */
  hid_t h_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Open the Header group */
  hid_t h_grp = H5Gopen(h_file, "Header", H5P_DEFAULT);

  /* Read the size of the perturbation */
  hid_t h_attr, h_err, h_tp;

  /* Number of wavenumbers */
  h_attr = H5Aopen(h_grp, "k_size", H5P_DEFAULT);
  h_err = H5Aread(h_attr, H5T_NATIVE_INT, &rend->N_wavenumbers);
  if (h_err) error("Error reading attribute");
  H5Aclose(h_attr);

  /* Number of time steps */
  h_attr = H5Aopen(h_grp, "tau_size", H5P_DEFAULT);
  h_err = H5Aread(h_attr, H5T_NATIVE_INT, &rend->N_redshifts);
  if (h_err) error("Error reading attribute");
  H5Aclose(h_attr);

  /* Number of transfer functions T(k,tau) */
  h_attr = H5Aopen(h_grp, "n_functions", H5P_DEFAULT);
  h_err = H5Aread(h_attr, H5T_NATIVE_INT, &rend->N_functions);
  if (h_err) error("Error reading attribute");
  H5Aclose(h_attr);
  
  const int N_k = rend->N_wavenumbers;
  const int N_z = rend->N_redshifts;
  const int N_f = rend->N_functions;
  
  message("Found %d %d %d", rend->N_wavenumbers, rend->N_redshifts, rend->N_functions);

  /* Allocate memory for the transfer function titles */
  char **titles = malloc(rend->N_functions * sizeof(char*));
    
  /* Read the titles of the transfer functions */
  h_attr = H5Aopen(h_grp, "FunctionTitles", H5P_DEFAULT);
  h_tp = H5Aget_type(h_attr);
  h_err = H5Aread(h_attr, h_tp, titles);
  if (h_err) error("Error reading transfer function titles");
  H5Aclose(h_attr);
  H5Tclose(h_tp);

  /* Titles of the cdm and neutrino density transfer functions */
  char d_cdm_name[PARSER_MAX_LINE_SIZE];
  char d_ncdm_name[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "Neutrino:title_d_cdm", filename);
  parser_get_param_string(params, "Neutrino:title_d_ncdm", filename);

  /* Check if we have the cdm and neutrino density transfer functions */
  int found_cdm = -1, found_ncdm = -1;
  for (int i=0; i<rend->N_functions; i++) {
      if (strcmp(titles[i], d_cdm_name) == 0)
      found_cdm = i;
      if (strcmp(titles[i], d_ncdm_name) == 0)
      found_ncdm = i;
  }
  
  if (found_cdm < 0 || found_ncdm < 0) {
    error("Did not find cdm or ncdm density transfer functions.");
  }
  
  /* Done with the titles */
  free (titles);
  
/* Read the units of the perturbation data file */
double UnitLengthCGS, UnitTimeCGS, UnitMassCGS;
double UnitLengthMetres, UnitTimeSeconds, UnitMassKilogram;

/* CGS Length units */
h_attr = H5Aopen(h_grp, "Unit length in cgs (U_L)", H5P_DEFAULT);
h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitLengthCGS);
H5Aclose(h_attr);

/* CGS Time units */
h_attr = H5Aopen(h_grp, "Unit time in cgs (U_t)", H5P_DEFAULT);
h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitTimeCGS);
H5Aclose(h_attr);

/* CGS Mass units */
h_attr = H5Aopen(h_grp, "Unit mass in cgs (U_M)", H5P_DEFAULT);
h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitMassCGS);
H5Aclose(h_attr);

/* Convert to SI units */
UnitLengthMetres = UnitLengthCGS/100;
UnitTimeSeconds = UnitTimeCGS;
UnitMassKilogram = UnitMassCGS/1000;

  /* Close the Header group */
  H5Gclose(h_grp);

  /* Open the data group */
  h_grp = H5Gopen(h_file, "Perturb", H5P_DEFAULT);

  /* Allocate memory */
  rend->redshifts = calloc(N_z, sizeof(double));
  rend->wavenumbers = calloc(N_k, sizeof(double));
  rend->ncdm_over_cdm = calloc(N_z * N_k, sizeof(double));
  
  /* Temporary array for all the transfer function data */
  double *transfer_functions = malloc(N_k * N_z * N_f * sizeof(double));
    
  /* Dataspace */
  hid_t h_data;
    
  /* Allocation successful? */
  if (rend->redshifts == NULL || rend->wavenumbers == NULL || rend->ncdm_over_cdm == NULL || transfer_functions == NULL) {
    error("ERROR: unable to allocate memory for perturbation data.");
  }
  
  /* Read the redshifts */
  h_data = H5Dopen2(h_grp, "Redshifts", H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rend->redshifts);
  if (h_err) error("Error reading redshifts");
  H5Dclose(h_data);
  
  /* Read the wavenumbers */
  h_data = H5Dopen2(h_grp, "Wavenumbers", H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rend->wavenumbers);
  if (h_err) error("Error reading wavenumbers");
  H5Dclose(h_data);

  /* Read the transfer functions */
  h_data = H5Dopen2(h_grp, "Transfer functions", H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, transfer_functions);
  if (h_err) error("Error reading transfer functions");
  H5Dclose(h_data);
  
  
    // 
    // if (h_err < 0 || pt->k[0] == 0 || pt->log_tau[0] == 0) {
    //     printf("ERROR: problem with reading the perturbation data.\n");
    //     return 1;
    // }
    // 
    // /* Close the data group */
    // H5Gclose(h_grp);
    // 
    // /* Close the file */
    // H5Fclose(h_file);
    // 
    // /* Perform unit conversions for the wavenumbers */
    // for (int i=0; i<pt->k_size; i++) {
    //     pt->k[i] *= us->UnitLengthMetres / UnitLengthMetres;
    // }
    // 
    // /* Perform unit conversions for the conformal times */
    // for (int i=0; i<pt->tau_size; i++) {
    //     pt->log_tau[i] += log(UnitTimeSeconds/us->UnitTimeSeconds);
    // }
    // 
    // const double unit_length_factor = UnitLengthMetres / us->UnitLengthMetres;
    // const double unit_time_factor = UnitTimeSeconds / us->UnitTimeSeconds;
    // 
    // if (fabs(1./unit_time_factor - 1) > 1e-5 ) {
    //   message(pars->rank, "Velocity factor = %e\n", 1./unit_time_factor);
    // }
    // 
    // /* Perform unit conversions for the Hubble rates */
    // for (int i=0; i<pt->tau_size; i++) {
    //     pt->Hubble_H[i] /= unit_time_factor;
    // }
    // 
    // if (fabs(unit_time_factor - 1) > 1e-5 ) {
    //   message(pars->rank, "Unit conversion factor for '%s' is %f\n", "Hubble_H", unit_time_factor);
    // }
    // 
    // /* Perform unit conversions for the transfer functions */
    // for (int i=0; i<pt->n_functions; i++) {
    //     /* Determine the unit conversion factor */
    //     const char *title = pt->titles[i];
    //     double unit_factor = unitConversionFactor(title, unit_length_factor, unit_time_factor);
    // 
    //     if (fabs(unit_factor - 1) > 1e-5 ) {
    //       message(pars->rank, "Unit conversion factor for '%s' is %f\n", title, unit_factor);
    //     }
    // 
    //     /* Convert from input units to internal units */
    //     for (int index_k=0; index_k<pt->k_size; index_k++) {
    //         for (int index_tau=0; index_tau<pt->tau_size; index_tau++) {
    //             int index = pt->tau_size * pt->k_size * i + pt->k_size * index_tau + index_k;
    //             pt->delta[index] *= unit_factor;
    //         }
    //     }
    // }
  
  
}
                  