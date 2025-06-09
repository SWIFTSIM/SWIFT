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
#include <config.h>

/* Local includes. */
#include "parser.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline2d.h>
#endif

/* Local headers */
#include "neutrino.h"

/*! Number of time steps in the neutrino density interpolation table */
const hsize_t timestep_length = 10000;
/*! Number of wavenumbers in the neutrino density interpolation table */
const hsize_t wavenumber_length = 10000;

/**
 * @brief Determine the length of an HDF5 dataset
 *
 * @param h_file The file object
 * @param title The title of the dataset
 * @param length (Output) The vector length
 */
void read_vector_length(hid_t h_file, char title[PARSER_MAX_LINE_SIZE],
                        hsize_t *length) {

  /* Open the dataset */
  hid_t h_data = H5Dopen2(h_file, title, H5P_DEFAULT);

  /* Open the dataspace and fetch the dimensions */
  hid_t h_space = H5Dget_space(h_data);

  /* Fetch the rank and size of the dataset */
  int ndims = H5Sget_simple_extent_ndims(h_space);
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(h_space, dims, NULL);

  /* Close the dataspace */
  H5Sclose(h_space);

  /* Verify that this is a vector and return the length */
  if (ndims != 1 || dims[0] <= 0) error("We expected a non-empty vector.");
  length[0] = dims[0];

  /* Close the dataset */
  H5Dclose(h_data);
}

/**
 * @brief Read transfer funtion from an HDF5 dataset with expected size N_z*N_k
 *
 * @param h_file The file object
 * @param title The title of the dataset
 * @param length (Output) Memory fot the transfer function data
 * @param N_z Expected number of timesteps
 * @param N_k Expected number of wavenumbers
 */
void read_transfer_function(hid_t h_file, char title[PARSER_MAX_LINE_SIZE],
                            double *dest, hsize_t N_z, hsize_t N_k) {

  /* Open the dataset */
  hid_t h_data = H5Dopen2(h_file, title, H5P_DEFAULT);

  /* Open the dataspace and fetch the dimensions */
  hid_t h_space = H5Dget_space(h_data);
  int ndims = H5Sget_simple_extent_ndims(h_space);
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(h_space, dims, NULL);

  /* Close the dataspace */
  H5Sclose(h_space);

  /* Verify the dimensions */
  if (ndims != 2) error("We expected a rank-2 tensor.");
  if (dims[0] != N_z)
    error("Number of redshifts does not match array dimensions.");
  if (dims[1] != N_k)
    error("Number of wavenumbers does not match array dimensions.");

  /* Read out the data */
  hid_t h_err =
      H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dest);
  if (h_err < 0) error("Error reading dataset '%s'.\n", title);

  /* Close the dataset */
  H5Dclose(h_data);
}

/**
 * @brief Initialize the #neutrino_response object by reading an HDF5 file
 * with transfer functions and pre-computing a transfer function ratio.
 *
 * @param numesh The #neutrino_mesh to be initialized
 * @param params The parsed parameter file.
 * @param us The system of units used internally.
 * @param dim Spatial dimensions of the domain.
 * @param c The #cosmology used for this run.
 * @param np The #neutrino_props used for this run.
 * @param gp The #gravity_props used for this run.
 * @param verbose Are we talkative ?
 */
void neutrino_response_init(struct neutrino_response *numesh,
                            struct swift_params *params,
                            const struct unit_system *us, const double dim[3],
                            const struct cosmology *c,
                            const struct neutrino_props *np,
                            const struct gravity_props *gp, int rank,
                            int verbose) {

  /* Do we need to do anything? */
  if (!np->use_linear_response) return;

  /* Check that we have a degenerate neutrino mass spectrum */
  if (c->N_nu > 1)
    error(
        "Non-degenerate neutrino mass spectra not supported with the linear "
        "response method.");

#ifdef HAVE_LIBGSL

  /* Parse file name parameter */
  char filename[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "Neutrino:transfer_functions_filename",
                          filename);

  /* Titles of the redshift, wavenumber, and transfer functions datasets */
  char z_name[PARSER_MAX_LINE_SIZE];
  char k_name[PARSER_MAX_LINE_SIZE];
  char d_cdm_name[PARSER_MAX_LINE_SIZE];
  char d_b_name[PARSER_MAX_LINE_SIZE];
  char d_ncdm_name[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "Neutrino:dataset_redshifts", z_name);
  parser_get_param_string(params, "Neutrino:dataset_wavenumbers", k_name);
  parser_get_param_string(params, "Neutrino:dataset_delta_cdm", d_cdm_name);
  parser_get_param_string(params, "Neutrino:dataset_delta_baryon", d_b_name);
  parser_get_param_string(params, "Neutrino:dataset_delta_nu", d_ncdm_name);

  /* Read additional parameters */
  numesh->fixed_bg_density =
      parser_get_param_int(params, "Neutrino:fixed_bg_density");

  if (rank == 0 && verbose)
    message("Reading transfer functions file '%s'", filename);

  /* Open the hdf5 file */
  hid_t h_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t h_data, h_err, h_status;

  /* Obtain the number of redshifts and wavenumbers in the file */
  hsize_t N_z;  // redshifts
  hsize_t N_k;  // wavenumbers
  read_vector_length(h_file, z_name, &N_z);
  read_vector_length(h_file, k_name, &N_k);

  /* The transfer function should be a rank-2 array of size (N_z, N_k) */
  const hsize_t tf_size = N_z * N_k;

  /* Allocate temporary memory for the transfer functions */
  double *delta_cdm;
  double *delta_b;
  double *delta_ncdm;
  double *ncdm_over_cb;
  if ((delta_cdm = (double *)swift_malloc("delta_cdm",
                                          sizeof(double) * tf_size)) == NULL)
    error("Failed to allocate memory for delta_cdm.");
  if ((delta_b = (double *)swift_malloc("delta_b", sizeof(double) * tf_size)) ==
      NULL)
    error("Failed to allocate memory for delta_b.");
  if ((delta_ncdm = (double *)swift_malloc("delta_ncdm",
                                           sizeof(double) * tf_size)) == NULL)
    error("Failed to allocate memory for delta_ncdm.");
  if ((ncdm_over_cb = (double *)swift_malloc("ncdm_over_cb",
                                             sizeof(double) * tf_size)) == NULL)
    error("Failed to allocate memory for ncdm_over_cb.");

  /* Read the necessary transfer functions */
  read_transfer_function(h_file, d_cdm_name, delta_cdm, N_z, N_k);
  read_transfer_function(h_file, d_b_name, delta_b, N_z, N_k);
  read_transfer_function(h_file, d_ncdm_name, delta_ncdm, N_z, N_k);

  /* Compute the background mass ratio of cdm to cb (= cdm + baryon) */
  const double f_cdm = c->Omega_cdm / (c->Omega_cdm + c->Omega_b);

  /* Compute the transfer function ratio */
  for (hsize_t i = 0; i < N_z; i++) {
    for (hsize_t j = 0; j < N_k; j++) {
      /* Fetch the data */
      double d_cdm = delta_cdm[N_k * i + j];
      double d_b = delta_b[N_k * i + j];
      double d_ncdm = delta_ncdm[N_k * i + j];

      /* Weighted baryon-cdm density */
      double d_cb = f_cdm * d_cdm + (1.0 - f_cdm) * d_b;

      /* Store the ratio */
      ncdm_over_cb[N_k * i + j] = d_ncdm / d_cb;
    }
  }

  /* Free temporary transfer function memory */
  swift_free("delta_cdm", delta_cdm);
  swift_free("delta_b", delta_b);
  swift_free("delta_ncdm", delta_ncdm);

  /* The length unit used by the transfer function file */
  double UnitLengthCGS;

  /* Check if a Units group exists */
  h_status = H5Eset_auto1(NULL, NULL);  // turn off error printing
  h_status = H5Gget_objinfo(h_file, "/Units", 0, NULL);

  /* If the group exists */
  if (h_status == 0) {
    hid_t h_attr, h_grp;

    /* Open the Units group */
    h_grp = H5Gopen(h_file, "Units", H5P_DEFAULT);

    /* Read the length unit in CGS units */
    h_attr = H5Aopen(h_grp, "Unit length in cgs (U_L)", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitLengthCGS);
    H5Aclose(h_attr);

    /* Close the Units group */
    H5Gclose(h_grp);
  } else {
    /* Assume the internal unit system */
    UnitLengthCGS = us->UnitLength_in_cgs;
  }

  /* Determine the range of wavenumbers needed for this simulation */
  const double boxlen = dim[0];
  const int mesh_size = gp->mesh_size;
  const double k_margin = 1.1;  // safety margin around the edges
  const double k_min = 2.0 * M_PI / boxlen / k_margin;
  const double k_max = sqrt(3.0) * k_min * mesh_size * k_margin * k_margin;
  const double k_unit = UnitLengthCGS / us->UnitLength_in_cgs;

  /* Determine the range of redshifts needed for this simulation */
  const double a_min = c->a_begin;
  const double a_max = c->a_end;
  const double z_max = 1.0 / a_min - 1.0;
  const double z_min = 1.0 / a_max - 1.0;

  /* Allocate temporary memory for the redshifts and wavenumbers */
  double *redshifts;
  double *wavenumbers;
  if ((redshifts = (double *)swift_malloc("redshifts", sizeof(double) * N_z)) ==
      NULL)
    error("Failed to allocate memory for redshifts.");
  if ((wavenumbers =
           (double *)swift_malloc("wavenumbers", sizeof(double) * N_k)) == NULL)
    error("Failed to allocate memory for wavenumbers.");

  /* Open the redshifts dataset, read the data, then close it */
  h_data = H5Dopen2(h_file, z_name, H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  redshifts);
  if (h_err < 0) error("Error reading dataset '%s'.\n", z_name);
  H5Dclose(h_data);

  /* Open the wavenumbers dataset, read the data, then close it */
  h_data = H5Dopen2(h_file, k_name, H5P_DEFAULT);
  h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  wavenumbers);
  if (h_err < 0) error("Error reading dataset '%s'.\n", k_name);
  H5Dclose(h_data);

  /* Ensure that the data is ascending in time and wavenumber */
  for (hsize_t i = 1; i < N_z; i++)
    if (redshifts[i] > redshifts[i - 1]) error("Redshifts not descending.");
  for (hsize_t i = 1; i < N_k; i++)
    if (wavenumbers[i] < wavenumbers[i - 1])
      error("Wavenumbers not ascending.");

  /* Ensure that we have data covering the required domain */
  if (redshifts[0] < z_max || redshifts[N_z - 1] > z_min)
    error("Redshifts do not cover the interval (%g, %g)", z_min, z_max);
  if (wavenumbers[0] * k_unit > k_min || wavenumbers[N_k - 1] * k_unit < k_max)
    error("Wavenumbers do not cover the interval (%g, %g)", k_min, k_max);

  if (rank == 0 && verbose) {
    message("We have transfer functions covering the domain:");
    message("(k_min, k_max) = (%g, %g) U_L^-1", k_min, k_max);
    message("(a_min, a_max) = (%g, %g)", a_min, a_max);
  }

  /* We will remap the data such that it just covers the required domain
   * with a constant log spacing, allowing for faster interpolation.
   * We only use the slower GSL interpolation for the remapping. */

  /* Determine the dimensions of the remapped data */
  const hsize_t remap_tf_size = timestep_length * wavenumber_length;

  /* Determine the constant log spacing and bounding values */
  const double log_a_min = log(a_min);
  const double log_a_max = log(a_max);
  const double log_k_min = log(k_min);
  const double log_k_max = log(k_max);
  const double delta_log_a = (log_a_max - log_a_min) / timestep_length;
  const double delta_log_k = (log_k_max - log_k_min) / wavenumber_length;

  /* Store the remapped bounding values and log spacing */
  numesh->log_a_min = log_a_min;
  numesh->log_a_max = log_a_max;
  numesh->log_k_min = log_k_min;
  numesh->log_k_max = log_k_max;
  numesh->delta_log_a = delta_log_a;
  numesh->delta_log_k = delta_log_k;

  /* Allocate temporary memory for the log of scale factors and wavenumbers */
  double *log_scale_factors;
  double *log_wavenumbers;
  if ((log_scale_factors = (double *)swift_malloc(
           "log_scale_factors", sizeof(double) * N_z)) == NULL)
    error("Failed to allocate memory for log_scale_factors.");
  if ((log_wavenumbers = (double *)swift_malloc("log_wavenumbers",
                                                sizeof(double) * N_k)) == NULL)
    error("Failed to allocate memory for log_wavenumbers.");

  /* Convert units and compute logarithms */
  for (hsize_t i = 0; i < N_z; i++) {
    log_scale_factors[i] = -log(1.0 + redshifts[i]);
  }
  for (hsize_t i = 0; i < N_k; i++) {
    log_wavenumbers[i] = log(wavenumbers[i] * k_unit);
  }

  /* Free the temporary buffers for the original redshifts and wavenumbers */
  swift_free("redshifts", redshifts);
  swift_free("wavenumbers", wavenumbers);

  /* Initialize GSL interpolation */
  /* NB: GSL uses column major indices, but we use row major. This is why we
   * feed the transpose into spline2d_init. */
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gsl_spline2d *spline = gsl_spline2d_alloc(T, N_k, N_z);
  gsl_interp_accel *a_acc = gsl_interp_accel_alloc();
  gsl_interp_accel *k_acc = gsl_interp_accel_alloc();
  gsl_spline2d_init(spline, log_wavenumbers, log_scale_factors, ncdm_over_cb,
                    N_k, N_z);

  /* Allocate memory for the transfer function ratio with constant spacing */
  if ((numesh->pt_density_ratio = (double *)swift_malloc(
           "numesh.pt_density_ratio", sizeof(double) * remap_tf_size)) == NULL)
    error("Failed to allocate memory for numesh.pt_density_ratio.");

  /* Remap the transfer function ratio */
  for (hsize_t i = 0; i < timestep_length; i++) {
    for (hsize_t j = 0; j < wavenumber_length; j++) {
      double log_a = log_a_min + i * delta_log_a;
      double log_k = log_k_min + j * delta_log_k;
      numesh->pt_density_ratio[wavenumber_length * i + j] =
          gsl_spline2d_eval(spline, log_k, log_a, k_acc, a_acc);
    }
  }

  /* Store the array size */
  numesh->tf_size = wavenumber_length * timestep_length;

  /* Clean up GSL interpolation */
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(a_acc);
  gsl_interp_accel_free(k_acc);

  /* Free the temporary buffers */
  swift_free("log_scale_factors", log_scale_factors);
  swift_free("log_wavenumbers", log_wavenumbers);
  swift_free("ncdm_over_cb", ncdm_over_cb);
#else
  error("No GSL library found. Cannot remap the transfer functions.");
#endif
}

void neutrino_response_clean(struct neutrino_response *numesh) {
  swift_free("numesh.pt_density_ratio", numesh->pt_density_ratio);
}

#ifdef HAVE_FFTW

/**
 * @brief Shared information about the neutrino response used by the threads.
 */
struct neutrino_response_tp_data {

  /* Mesh properties */
  int N;
  fftw_complex *frho;
  double boxlen;
  int slice_offset;
  int slice_width;

  /* Interpolation properties */
  double inv_delta_log_k;
  double log_k_min;
  int a_index;
  double u_a;

  /* Background and perturbed density ratios */
  double bg_density_ratio;
  const double *pt_density_ratio;
};

/**
 * @brief Mapper function for the application of the linear neutrino response.
 *
 * @param map_data The array of the density field Fourier transform.
 * @param num The number of elements to iterate on (along the x-axis).
 * @param extra The properties of the neutrino mesh.
 */
void neutrino_response_apply_neutrino_response_mapper(void *map_data,
                                                      const int num,
                                                      void *extra) {

  struct neutrino_response_tp_data *data =
      (struct neutrino_response_tp_data *)extra;

  /* Unpack the mesh properties */
  fftw_complex *const frho = data->frho;
  const int N = data->N;
  const int N_half = N / 2;
  const double delta_k = 2.0 * M_PI / data->boxlen;

  /* Unpack the interpolation properties */
  const hsize_t a_index = data->a_index;
  const double u_a = data->u_a;
  const double inv_delta_log_k = data->inv_delta_log_k;
  const double log_k_min = data->log_k_min;
  const int N_k = wavenumber_length;

  /* Unpack the density ratios (background & perturbation) */
  const double bg_density_ratio = data->bg_density_ratio;
  const double *pt_density_ratio = data->pt_density_ratio;

  /* Find what slice of the full mesh is stored on this MPI rank */
  const int slice_offset = data->slice_offset;

  /* Range of x coordinates in the full mesh handled by this call */
  const int x_start = ((fftw_complex *)map_data - frho) + slice_offset;
  const int x_end = x_start + num;

  /* Loop over the x range corresponding to this thread */
  for (int x = x_start; x < x_end; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z < N_half + 1; z++) {

        /* Compute the wavevector (U_L^-1) */
        const double k_x = (x > N_half) ? (x - N) * delta_k : x * delta_k;
        const double k_y = (y > N_half) ? (y - N) * delta_k : y * delta_k;
        const double k_z = (z > N_half) ? (z - N) * delta_k : z * delta_k;
        const double k2 = k_x * k_x + k_y * k_y + k_z * k_z;

        /* Skip the DC mode */
        if (k2 == 0.) continue;

        /* Interpolate along the k-axis */
        const double log_k = 0.5 * log(k2);
        const double log_k_steps = (log_k - log_k_min) * inv_delta_log_k;
        const hsize_t k_index = (hsize_t)log_k_steps;
        const double u_k = log_k_steps - k_index;

        /* Retrieve the bounding values */
        const double T11 = pt_density_ratio[N_k * a_index + k_index];
        const double T21 = pt_density_ratio[N_k * a_index + k_index + 1];
        const double T12 = pt_density_ratio[N_k * (a_index + 1) + k_index];
        const double T22 = pt_density_ratio[N_k * (a_index + 1) + k_index + 1];

        /* Bilinear interpolation of the tranfer function ratio */
        double pt_ratio_interp = (1.0 - u_a) * ((1.0 - u_k) * T11 + u_k * T21) +
                                 u_a * ((1.0 - u_k) * T12 + u_k * T22);
        double correction = 1.0 + pt_ratio_interp * bg_density_ratio;

#ifdef SWIFT_DEBUG_CHECKS
        if (u_k < 0 || u_a < 0 || u_k > 1 || u_a > 1 ||
            k_index > wavenumber_length || a_index > timestep_length)
          error("Interpolation out of bounds error: %g %g %g %g %llu %llu\n",
                u_k, u_a, sqrt(k2), pt_ratio_interp,
                (unsigned long long)k_index, (unsigned long long)a_index);
#endif

        /* Apply to the mesh */
        const int index =
            N * (N_half + 1) * (x - slice_offset) + (N_half + 1) * y + z;
        frho[index][0] *= correction;
        frho[index][1] *= correction;
      }
    }
  }
}

/**
 * @brief Apply the linear neutrino response to the Fourier transform of the
 * gravitational potential.
 *
 * @param s The current #space
 * @param mesh The #pm_mesh used to store the potential
 * @param tp The #threadpool object used for parallelisation
 * @param frho The NxNx(N/2) complex array of the Fourier transform of the
 * density field
 * @param slice_offset The x coordinate of the start of the slice on this MPI
 * rank
 * @param slice_width The width of the local slice on this MPI rank
 * @param verbose Are we talkative?
 */
void neutrino_response_compute(const struct space *s, struct pm_mesh *mesh,
                               struct threadpool *tp, fftw_complex *frho,
                               const int slice_offset, const int slice_width,
                               int verbose) {
#ifdef HAVE_FFTW

  const struct cosmology *c = s->e->cosmology;
  struct neutrino_response *numesh = s->e->neutrino_response;

  /* Grid size */
  const int N = mesh->N;
  const double boxlen = mesh->dim[0];

  /* Calculate the background neutrino density */
  const double a = numesh->fixed_bg_density ? 1.0 : c->a;
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_m = c->Omega_cdm + c->Omega_b;  // does not include nu's
  /* The comoving density is (Omega_nu * a^-4) * a^3  = Omega_nu / a */
  const double bg_density_ratio = (Omega_nu / a) / Omega_m;

  /* Transfer function bounds and spacing */
  const double inv_delta_log_a = 1.0 / numesh->delta_log_a;
  const double inv_delta_log_k = 1.0 / numesh->delta_log_k;
  const double log_a_min = numesh->log_a_min;
  const double log_k_min = numesh->log_k_min;

  /* Interpolate along the a-axis */
  const double log_a = log(c->a);
  const double log_a_steps = (log_a - log_a_min) * inv_delta_log_a;
  hsize_t a_index = (hsize_t)log_a_steps;
  double u_a = log_a_steps - a_index;

  /* Perform bounds checks for a */
  if (log_a_steps < 0) {
    a_index = 0;
    u_a = 0;
  } else if (a_index > timestep_length - 2) {
    a_index = timestep_length - 2;
    u_a = 1.0;
  }

  /* Some common factors */
  struct neutrino_response_tp_data data;
  data.frho = frho;
  data.N = N;
  data.boxlen = boxlen;
  data.slice_offset = slice_offset;
  data.slice_width = slice_width;
  data.inv_delta_log_k = inv_delta_log_k;
  data.log_k_min = log_k_min;
  data.a_index = a_index;
  data.u_a = u_a;
  data.bg_density_ratio = bg_density_ratio;
  data.pt_density_ratio = numesh->pt_density_ratio;

  /* Parallelize the neutrino linear response application using the threadpool
     to split the x-axis loop over the threads. The array is N x N x (N/2).
     We use the thread to each deal with a range [i_min, i_max[ x N x (N/2) */
  threadpool_map(tp, neutrino_response_apply_neutrino_response_mapper, frho,
                 slice_width, sizeof(fftw_complex), threadpool_auto_chunk_size,
                 &data);

  /* Correct singularity at (0,0,0) */
  if (slice_offset == 0 && slice_width > 0) {
    frho[0][0] = 0.;
    frho[0][1] = 0.;
  }

#else
  error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
}

#endif /* HAVE FFTW */

/**
 * @brief Write a neutrino response struct to the given FILE as a stream of
 * bytes.
 *
 * @param numesh the struct
 * @param stream the file stream
 */
void neutrino_response_struct_dump(const struct neutrino_response *numesh,
                                   FILE *stream) {
  restart_write_blocks((void *)numesh, sizeof(struct neutrino_response), 1,
                       stream, "numesh", "neutrino mesh");

  /* Store the perturbation data */
  if (numesh->tf_size > 0) {
    restart_write_blocks((double *)numesh->pt_density_ratio, sizeof(double),
                         numesh->tf_size, stream, "pt_density_ratio",
                         "pt_density_ratio");
  }
}

/**
 * @brief Restore a neutrino response struct from the given FILE as a stream
 * of bytes.
 *
 * @param numesh the struct
 * @param stream the file stream
 */
void neutrino_response_struct_restore(struct neutrino_response *numesh,
                                      FILE *stream) {
  restart_read_blocks((void *)numesh, sizeof(struct neutrino_response), 1,
                      stream, NULL, "neutrino mesh");

  /* Restore the perturbation data */
  if (numesh->tf_size > 0) {
    numesh->pt_density_ratio = (double *)swift_malloc(
        "numesh.pt_density_ratio", numesh->tf_size * sizeof(double));
    restart_read_blocks((double *)numesh->pt_density_ratio, sizeof(double),
                        numesh->tf_size, stream, NULL, "pt_density_ratio");
  }
}
