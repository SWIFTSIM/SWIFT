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
                   const struct cosmology *c, const struct neutrino_props *np,
                   struct neutrino_renderer *rend) {

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

  message("Found %d %d %d", rend->N_wavenumbers, rend->N_redshifts,
          rend->N_functions);

  /* Allocate memory for the transfer function titles */
  char **titles = malloc(rend->N_functions * sizeof(char *));

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
  parser_get_param_string(params, "Neutrino:title_d_cdm", d_cdm_name);
  parser_get_param_string(params, "Neutrino:title_d_ncdm", d_ncdm_name);

  /* Check if we have the cdm and neutrino density transfer functions */
  int index_cdm = -1, index_ncdm = -1;
  for (int i = 0; i < rend->N_functions; i++) {
    if (strcmp(titles[i], d_cdm_name) == 0) index_cdm = i;
    if (strcmp(titles[i], d_ncdm_name) == 0) index_ncdm = i;
  }

  if (index_cdm < 0 || index_ncdm < 0) {
    error("Did not find cdm or ncdm density transfer functions.");
  }

  /* Done with the titles */
  free(titles);

  /* Read the units of the perturbation data file */
  double UnitLengthCGS;

  /* CGS Length units */
  h_attr = H5Aopen(h_grp, "Unit length in cgs (U_L)", H5P_DEFAULT);
  h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitLengthCGS);
  H5Aclose(h_attr);

  /* Close the Header group */
  H5Gclose(h_grp);

  /* Open the data group */
  h_grp = H5Gopen(h_file, "Perturb", H5P_DEFAULT);

  /* Allocate memory */
  rend->redshifts = calloc(N_z, sizeof(double));
  rend->log_wavenumbers = calloc(N_k, sizeof(double));
  rend->ncdm_over_cdm = calloc(N_z * N_k, sizeof(double));

  /* Temporary array for the wavenumbers */
  double *wavenumbers = calloc(N_k, sizeof(double));

  /* Temporary array for all the transfer function data */
  double *transfer_functions = malloc(N_k * N_z * N_f * sizeof(double));

  /* Dataspace */
  hid_t h_data;

  /* Allocation successful? */
  if (rend->redshifts == NULL || rend->log_wavenumbers == NULL ||
      rend->ncdm_over_cdm == NULL || wavenumbers == NULL || 
      transfer_functions == NULL) {
    error("ERROR: unable to allocate memory for perturbation data.");
  }

  /* Read the redshifts */
  h_data = H5Dopen2(h_grp, "Redshifts", H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  rend->redshifts);
  if (h_err) error("Error reading redshifts");
  H5Dclose(h_data);

  /* Read the wavenumbers */
  h_data = H5Dopen2(h_grp, "Wavenumbers", H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  wavenumbers);
  if (h_err) error("Error reading wavenumbers");
  H5Dclose(h_data);

  /* Read the transfer functions */
  h_data = H5Dopen2(h_grp, "Transfer functions", H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  transfer_functions);
  if (h_err) error("Error reading transfer functions");
  H5Dclose(h_data);

  /* Close the data group */
  H5Gclose(h_grp);

  /* Close the file */
  H5Fclose(h_file);

  /* Perform unit conversions for the wavenumbers and store the logarithm */
  for (int i = 0; i < N_k; i++) {
    wavenumbers[i] *= us->UnitLength_in_cgs / UnitLengthCGS;
    rend->log_wavenumbers[i] = log(wavenumbers[i]);
  }

  /* Redshifts and transfer function ratios are dimensionless */

  /* Compute the transfer function ratio */
  for (int i = 0; i < N_z; i++) {
    for (int j = 0; j < N_k; j++) {
      /* Indices of the cdm and ncdm densities at this (z, k) pair */
      int k_cdm = N_z * N_k * index_cdm + N_k * i + j;
      int k_ncdm = N_z * N_k * index_ncdm + N_k * i + j;

      /* Fetch the data */
      double d_cdm = transfer_functions[k_cdm];
      double d_ncdm = transfer_functions[k_ncdm];

      /* Store the ratio */
      rend->ncdm_over_cdm[N_k * i + j] = d_ncdm / d_cdm;
    }
  }

  /* Done with the transfer functions */
  free(transfer_functions);
}

void renderer_clean(struct neutrino_renderer *rend) {
  free(rend->redshifts);
  free(rend->log_wavenumbers);
  free(rend->ncdm_over_cdm);
}