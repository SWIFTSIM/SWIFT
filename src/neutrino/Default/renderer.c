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

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

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

/* Write binary boxes in HDF5 format */
int writeGRF_H5(const double *box, int N, double boxlen, const char *fname) {
  /* Create the hdf5 file */
  hid_t h_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Create the Header group */
  hid_t h_grp =
      H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Create dataspace for BoxSize attribute */
  const hsize_t arank = 1;
  const hsize_t adims[1] = {3};  // 3D space
  hid_t h_aspace = H5Screate_simple(arank, adims, NULL);

  /* Create the BoxSize attribute and write the data */
  hid_t h_attr =
      H5Acreate1(h_grp, "BoxSize", H5T_NATIVE_DOUBLE, h_aspace, H5P_DEFAULT);
  double boxsize[3] = {boxlen, boxlen, boxlen};
  H5Awrite(h_attr, H5T_NATIVE_DOUBLE, boxsize);

  /* Close the attribute, corresponding dataspace, and the Header group */
  H5Aclose(h_attr);
  H5Sclose(h_aspace);
  H5Gclose(h_grp);

  /* Create the Field group */
  h_grp = H5Gcreate(h_file, "/Field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Create dataspace for the field */
  const hsize_t frank = 3;
  const hsize_t fdims[3] = {N, N, N};  // 3D space
  hid_t h_fspace = H5Screate_simple(frank, fdims, NULL);

  /* Create the dataset for the field */
  hid_t h_data = H5Dcreate(h_grp, "Field", H5T_NATIVE_DOUBLE, h_fspace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write the data */
  H5Dwrite(h_data, H5T_NATIVE_DOUBLE, h_fspace, h_fspace, H5P_DEFAULT, box);

  /* Close the dataset, corresponding dataspace, and the Field group */
  H5Dclose(h_data);
  H5Sclose(h_fspace);
  H5Gclose(h_grp);

  /* Close the file */
  H5Fclose(h_file);

  return 0;
}

/* Find index along the time direction (0 <= u <= 1 between index and index+1)
 */
void interp_redshift(const struct neutrino_renderer *rend, const double z,
                     int *index, double *u) {

  /* Number of bins */
  int N_z = rend->N_redshifts;

  /* Quickly return if we are in the first or last bin */
  if (z <= rend->redshifts[N_z - 1]) {
    *index = N_z - 2;
    *u = 1.0;
    return;
  } else if (z >= rend->redshifts[0]) {
    *index = 0;
    *u = 0.0;
    return;
  }

  /* Find i such that z[i] <= z */
  int i;
  for (i = 1; i < N_z && rend->redshifts[i] >= z; i++)
    ;

  *index = i - 1;

  /* Find the bounding values */
  double left = rend->redshifts[*index];
  double right = rend->redshifts[*index + 1];

  /* Calculate the ratio (X - X_left) / (X_right - X_left) */
  *u = (z - left) / (right - left);
}

/* Find index along the k-direction (0 <= u <= 1 between index and index+1) */
void interp_wavenumber(const struct neutrino_renderer *rend, const float log_k,
                       int *index, double *u) {

  /* Number of bins */
  int N_k = rend->N_wavenumbers;

  /* Quickly return if we are in the first or last bin */
  if (log_k <= rend->log_wavenumbers[0]) {
    *index = 0;
    *u = 0.0;
    return;
  } else if (log_k >= rend->log_wavenumbers[N_k - 1]) {
    *index = N_k - 2;
    *u = 1.0;
    return;
  }

  /* Find i such that log_k[i] <= log_k*/
  int i;
  for (i = 1; i < N_k && rend->log_wavenumbers[i] < log_k; i++)
    ;

  *index = i - 1;

  /* Find the bounding values */
  double left = rend->log_wavenumbers[*index];
  double right = rend->log_wavenumbers[*index + 1];

  /* Calculate the ratio (X - X_left) / (X_right - X_left) */
  *u = (log_k - left) / (right - left);
}

void renderer_compute(const struct cosmology *c,
                      const struct neutrino_props *np,
                      struct neutrino_renderer *rend, struct pm_mesh *mesh,
                      int verbose) {
#ifdef HAVE_FFTW

  /* Grid size */
  const int N = mesh->N;
  const int N_half = N / 2;
  const double boxlen = mesh->dim[0];
  const double boxvol = boxlen * boxlen * boxlen;
  const double delta_k = 2 * M_PI / boxlen;  // U_L^-1

  /* Calculate the background neutrino density at the present time */
  const double Omega_nu = cosmology_get_neutrino_density(c, c->a);
  const double Omega_m = c->Omega_m;  // does not include neutrinos
  /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
  const double bg_density_ratio = (Omega_nu / c->a) / Omega_m;

  /* Bilinear interpolation indices in (z, log_k) space */
  int z_index = 0, k_index = 0;
  double u_z = 0.0, u_k = 0.0;

  /* Number of wavenumbers in the perturbation vector */
  const int N_k = rend->N_wavenumbers;

  interp_redshift(rend, c->z, &z_index, &u_z);

  message("index = %d, u = %f (%f)", z_index, u_z,
          rend->redshifts[z_index + 1]);

  /* Long-range potential smoothing length */
  const double r_s = mesh->r_s;

  /* Boxes in configuration and momentum space */
  double *restrict potential;
  fftw_complex *restrict fp;

  /* Allocate memory for the rendered field */
  potential = (double *)fftw_malloc(sizeof(double) * N * N * N);
  fp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N * (N / 2 + 1));

  if (potential == NULL || fp == NULL) {
    error("Error allocating memory for rendering.");
  }

  memuse_log_allocation("potential", potential, 1, sizeof(double) * N * N * N);
  memuse_log_allocation("f", fp, 1, sizeof(fftw_complex) * N * N * (N / 2 + 1));

  /* Prepare the FFTW plans */
  fftw_plan pr2c = fftw_plan_dft_r2c_3d(N, N, N, potential, fp, FFTW_ESTIMATE);
  fftw_plan pc2r = fftw_plan_dft_c2r_3d(N, N, N, fp, potential, FFTW_ESTIMATE);

  /* Start a timer */
  ticks tic = getticks();

  /* First, copy the matter potential from mesh_gravity into the array */
  memcpy(potential, mesh->potential, N * N * N * sizeof(double));

  /* Transform to momentum space */
  fftw_execute(pr2c);

  /* Normalization */
  for (int i = 0; i < N * N * (N / 2 + 1); i++) {
    fp[i][0] *= boxvol / (N * N * N);
    fp[i][1] *= boxvol / (N * N * N);
  }

  if (verbose)
    message("Forward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Apply the transfer function */
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z <= N / 2; z++) {
        /* Compute the wavenumber */
        float k_x = (x > N / 2) ? (x - N) * delta_k : x * delta_k;  // U_L^-1
        float k_y = (y > N / 2) ? (y - N) * delta_k : y * delta_k;  // U_L^-1
        float k_z = (z > N / 2) ? (z - N) * delta_k : z * delta_k;  // U_L^-1
        float k = sqrtf(k_x * k_x + k_y * k_y + k_z * k_z);

        int index = N * (N_half + 1) * x + (N_half + 1) * y + z;

        /* Ignore the DC mode */
        if (k > 0) {
          /* Find the k-space interpolation index */
          interp_wavenumber(rend, logf(k), &k_index, &u_k);

          /* Retrieve the bounding values */
          double T11 = rend->ncdm_over_cdm[N_k * z_index + k_index];
          double T21 = rend->ncdm_over_cdm[N_k * z_index + k_index + 1];
          double T12 = rend->ncdm_over_cdm[N_k * (z_index + 1) + k_index];
          double T22 = rend->ncdm_over_cdm[N_k * (z_index + 1) + k_index + 1];

          /* Bilinear interpolation of the tranfer function ratio */
          double ncdm_over_cdm = (1.0 - u_z) * ((1.0 - u_k) * T11 + u_k * T21) +
                                 u_z * ((1.0 - u_k) * T12 + u_k * T22);

          /* The long-range kernel */
          double K = 1.;
          fourier_kernel_long_grav_eval(k * k * r_s * r_s, &K);

          fp[index][0] *= ncdm_over_cdm * bg_density_ratio / K;
          fp[index][1] *= ncdm_over_cdm * bg_density_ratio / K;
        } else {
          fp[index][0] = 0;
          fp[index][1] = 0;
        }
      }
    }
  }

  if (verbose)
    message("Applying transfer function took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Transform back */
  fftw_execute(pc2r);

  /* Normalization */
  for (int i = 0; i < N * N * N; i++) {
    potential[i] /= boxvol;
  }

  if (verbose)
    message("Backward Fourier transform took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Time the next stage */
  tic = getticks();

  /* Export the potentials if necessary */
  writeGRF_H5(potential, mesh->N, mesh->dim[0], "scaled_potential.hdf5");
  writeGRF_H5(mesh->potential, mesh->N, mesh->dim[0], "m_potential.hdf5");

  /* Print some statistics */
  double rms_matter = 0.f;
  double rms_nu = 0.f;
  for (int i = 0; i < N * N * N; i++) {
    rms_matter += mesh->potential[i] * mesh->potential[i];
    rms_nu += potential[i] * potential[i];
  }

  message("Dumping render boxes took %.3f %s.",
          clocks_from_ticks(getticks() - tic), clocks_getunit());
  message("[Phi_m, Phi_nu] = [%e, %e].", rms_matter, rms_nu);

  /* Add the contribution to the gravity mesh */
  for (int i = 0; i < N * N * N; i++) {
    mesh->potential[i] += potential[i];
  }

  /* Free memory */
  fftw_free(potential);
  fftw_free(fp);
  fftw_destroy_plan(pr2c);
  fftw_destroy_plan(pc2r);

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

/**
 * @brief Write a renderer struct to the given FILE as a stream of bytes.
 *
 * @param renderer the struct
 * @param stream the file stream
 */
void renderer_struct_dump(const struct neutrino_renderer *rend, FILE *stream) {
  restart_write_blocks((void *)rend, sizeof(struct neutrino_renderer), 1,
                       stream, "rend", "neutrino renderer");

  /* Store the perturbation data */
  restart_write_blocks((double *)rend->redshifts, sizeof(double),
                       rend->N_redshifts, stream, "rend->redshifts",
                       "redshifts");
  restart_write_blocks((double *)rend->log_wavenumbers, sizeof(double),
                       rend->N_wavenumbers, stream, "rend->log_wavenumbers",
                       "log_wavenumbers");
  restart_write_blocks((double *)rend->ncdm_over_cdm, sizeof(double),
                       rend->N_wavenumbers * rend->N_redshifts, stream,
                       "rend->ncdm_over_cdm", "ncdm_over_cdm");
}

/**
 * @brief Restore a neutrino renderer struct from the given FILE as a stream of
 * bytes.
 *
 * @param rend the struct
 * @param stream the file stream
 */
void renderer_struct_restore(struct neutrino_renderer *rend, FILE *stream) {
  restart_read_blocks((void *)rend, sizeof(struct neutrino_renderer), 1, stream,
                      NULL, "neutrino renderer");

  /* Restore the perturbation data */
  rend->redshifts =
      (double *)swift_malloc("redshifts", rend->N_redshifts * sizeof(double));
  restart_read_blocks((double *)rend->redshifts, sizeof(double),
                      rend->N_redshifts, stream, NULL, "redshifts");
  rend->log_wavenumbers = (double *)swift_malloc(
      "log_wavenumbers", rend->N_wavenumbers * sizeof(double));
  restart_read_blocks((double *)rend->log_wavenumbers, sizeof(double),
                      rend->N_wavenumbers, stream, NULL, "log_wavenumbers");
  rend->ncdm_over_cdm = (double *)swift_malloc(
      "ncdm_over_cdm",
      rend->N_redshifts * rend->N_wavenumbers * sizeof(double));
  restart_read_blocks((double *)rend->ncdm_over_cdm, sizeof(double),
                      rend->N_redshifts * rend->N_wavenumbers, stream, NULL,
                      "ncdm_over_cdm");
}