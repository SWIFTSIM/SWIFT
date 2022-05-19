/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Some standard headers. */
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/* This object's header. */
#include "gravity.h"

/* Local headers. */
#include "active.h"
#include "error.h"
#include "kernel_gravity.h"
#include "kernel_long_gravity.h"
#include "threadpool.h"
#include "version.h"

struct exact_force_data {
  const struct engine *e;
  const struct space *s;
  int counter_global;
  double const_G;
};

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

/* Size of the Ewald table */
#define Newald 64

/* Components of the Ewald correction */
static float fewald_x[Newald + 1][Newald + 1][Newald + 1];
static float fewald_y[Newald + 1][Newald + 1][Newald + 1];
static float fewald_z[Newald + 1][Newald + 1][Newald + 1];
static float potewald[Newald + 1][Newald + 1][Newald + 1];

/* Factor used to normalize the access to the Ewald table */
float ewald_fac;
#endif

/**
 * @brief Allocates the memory and computes one octant of the
 * Ewald correction table.
 *
 * We follow Hernquist, Bouchet & Suto, 1991, ApJS, Volume 75, p.231-240,
 * equations (2.14a) and (2.14b) with alpha = 2. We consider all terms with
 * |x - nL| < 4L and |h|^2 < 16.
 *
 * @param boxSize The side-length (L) of the volume.
 */
void gravity_exact_force_ewald_init(const double boxSize) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  const float boxSize_inv = 1.f / boxSize;
  const float boxSize_inv2 = 1.f / (boxSize * boxSize);

  int use_file = 0;

#ifdef HAVE_HDF5
  if (access("Ewald.hdf5", R_OK) != -1) use_file = 1;
#endif

  /* Can we use the stored HDF5 file? */
  if (use_file) {

    const ticks tic = getticks();
    message("Reading Ewald correction table from file...");

#ifdef HAVE_HDF5
    const hid_t h_file = H5Fopen("Ewald.hdf5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h_file < 0) error("Error opening the old 'Ewald.hdf5' file.");

    /* Check the table size */
    const hid_t h_grp = H5Gopen1(h_file, "Info");
    const hid_t h_attr = H5Aopen(h_grp, "Ewald_size", H5P_DEFAULT);
    int size;
    H5Aread(h_attr, H5T_NATIVE_INT, &size);
    if (size != Newald)
      error("File 'Ewald.hdf5' contains arrays of the wrong size");
    H5Aclose(h_attr);
    H5Gclose(h_grp);

    /* Now read the tables themselves */
    hid_t h_data;
    h_data = H5Dopen(h_file, "Ewald_x", H5P_DEFAULT);
    H5Dread(h_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(fewald_x[0][0][0]));
    H5Dclose(h_data);
    h_data = H5Dopen(h_file, "Ewald_y", H5P_DEFAULT);
    H5Dread(h_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(fewald_y[0][0][0]));
    H5Dclose(h_data);
    h_data = H5Dopen(h_file, "Ewald_z", H5P_DEFAULT);
    H5Dread(h_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(fewald_z[0][0][0]));
    H5Dclose(h_data);
    h_data = H5Dopen(h_file, "Ewald_pot", H5P_DEFAULT);
    H5Dread(h_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(potewald[0][0][0]));
    H5Dclose(h_data);

    /* Done */
    H5Fclose(h_file);
#endif

    /* Report time this took */
    message("Ewald correction table read in (took %.3f %s). ",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
  } else {

    /* Ok.. let's recompute everything */
    const ticks tic = getticks();
    message("Computing Ewald correction table...");

    /* Level of correction  (Hernquist et al. 1991)*/
    const float alpha = 2.f;

    /* some useful constants */
    const float alpha2 = alpha * alpha;
    const float factor_exp1 = 2.f * alpha / sqrt(M_PI);
    const float factor_exp2 = -M_PI * M_PI / alpha2;
    const float factor_sin = 2.f * M_PI;
    const float factor_cos = 2.f * M_PI;
    const float factor_pot = M_PI / alpha2;

    /* Zero everything */
    bzero(fewald_x, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));
    bzero(fewald_y, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));
    bzero(fewald_z, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));
    bzero(potewald, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));

    /* Hernquist, Bouchet & Suto, 1991, Eq. 2.10 and just below Eq. 2.15 */
    potewald[0][0][0] = 2.8372975f;

    /* Compute the values in one of the octants */
    for (int i = 0; i <= Newald; ++i) {
      for (int j = 0; j <= Newald; ++j) {
        for (int k = 0; k <= Newald; ++k) {

          if (i == 0 && j == 0 && k == 0) continue;

          /* Distance vector */
          const float r_x = 0.5f * ((float)i) / Newald;
          const float r_y = 0.5f * ((float)j) / Newald;
          const float r_z = 0.5f * ((float)k) / Newald;

          /* Norm of distance vector */
          const float r2 = r_x * r_x + r_y * r_y + r_z * r_z;
          const float r_inv = 1.f / sqrtf(r2);
          const float r_inv3 = r_inv * r_inv * r_inv;

          /* Normal gravity potential term */
          float f_x = r_x * r_inv3;
          float f_y = r_y * r_inv3;
          float f_z = r_z * r_inv3;
          float pot = r_inv + factor_pot;

          for (int n_i = -4; n_i <= 4; ++n_i) {
            for (int n_j = -4; n_j <= 4; ++n_j) {
              for (int n_k = -4; n_k <= 4; ++n_k) {

                const float d_x = r_x - n_i;
                const float d_y = r_y - n_j;
                const float d_z = r_z - n_k;

                /* Discretised distance */
                const float r_tilde2 = d_x * d_x + d_y * d_y + d_z * d_z;
                const float r_tilde_inv = 1.f / sqrtf(r_tilde2);
                const float r_tilde = r_tilde_inv * r_tilde2;
                const float r_tilde_inv3 =
                    r_tilde_inv * r_tilde_inv * r_tilde_inv;

                const float val_pot = erfcf(alpha * r_tilde);

                const float val_f =
                    val_pot + factor_exp1 * r_tilde * expf(-alpha2 * r_tilde2);

                /* First correction term */
                const float f = val_f * r_tilde_inv3;
                f_x -= f * d_x;
                f_y -= f * d_y;
                f_z -= f * d_z;
                pot -= val_pot * r_tilde_inv;
              }
            }
          }

          for (int h_i = -4; h_i <= 4; ++h_i) {
            for (int h_j = -4; h_j <= 4; ++h_j) {
              for (int h_k = -4; h_k <= 4; ++h_k) {

                const float h2 = h_i * h_i + h_j * h_j + h_k * h_k;

                if (h2 == 0.f) continue;

                const float h2_inv = 1.f / (h2 + FLT_MIN);
                const float h_dot_x = h_i * r_x + h_j * r_y + h_k * r_z;

                const float common = h2_inv * expf(h2 * factor_exp2);

                const float val_pot =
                    (float)M_1_PI * common * cosf(factor_cos * h_dot_x);

                const float val_f = 2.f * common * sinf(factor_sin * h_dot_x);

                /* Second correction term */
                f_x -= val_f * h_i;
                f_y -= val_f * h_j;
                f_z -= val_f * h_k;
                pot -= val_pot;
              }
            }
          }

          /* Save back to memory */
          fewald_x[i][j][k] = f_x;
          fewald_y[i][j][k] = f_y;
          fewald_z[i][j][k] = f_z;
          potewald[i][j][k] = pot;
        }
      }
    }

/* Dump the Ewald table to a file */
#ifdef HAVE_HDF5
    hid_t h_file =
        H5Fcreate("Ewald.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Error while opening file for Ewald dump.");

    /* Write the Ewald table size */
    const int size = Newald;
    const hid_t h_grp = H5Gcreate1(h_file, "Info", 0);
    const hid_t h_aspace = H5Screate(H5S_SCALAR);
    hid_t h_att =
        H5Acreate1(h_grp, "Ewald_size", H5T_NATIVE_INT, h_aspace, H5P_DEFAULT);
    H5Awrite(h_att, H5T_NATIVE_INT, &size);
    H5Aclose(h_att);
    H5Gclose(h_grp);
    H5Sclose(h_aspace);

    /* Create dataspace and write arrays */
    hsize_t dim[3] = {Newald + 1, Newald + 1, Newald + 1};
    hid_t h_space = H5Screate_simple(3, dim, NULL);
    hid_t h_data;
    h_data = H5Dcreate(h_file, "Ewald_x", H5T_NATIVE_FLOAT, h_space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
             &(fewald_x[0][0][0]));
    H5Dclose(h_data);
    h_data = H5Dcreate(h_file, "Ewald_y", H5T_NATIVE_FLOAT, h_space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
             &(fewald_y[0][0][0]));
    H5Dclose(h_data);
    h_data = H5Dcreate(h_file, "Ewald_z", H5T_NATIVE_FLOAT, h_space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
             &(fewald_z[0][0][0]));
    H5Dclose(h_data);
    h_data = H5Dcreate(h_file, "Ewald_pot", H5T_NATIVE_FLOAT, h_space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
             &(potewald[0][0][0]));
    H5Dclose(h_data);
    H5Sclose(h_space);
    H5Fclose(h_file);
#endif

    /* Report time this took */
    message("Ewald correction table computed (took %.3f %s). ",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

  /* Ewald factor to access the table */
  ewald_fac = (double)(2 * Newald) / boxSize;

  /* Apply the box-size correction */
  for (int i = 0; i <= Newald; ++i) {
    for (int j = 0; j <= Newald; ++j) {
      for (int k = 0; k <= Newald; ++k) {
        fewald_x[i][j][k] *= boxSize_inv2;
        fewald_y[i][j][k] *= boxSize_inv2;
        fewald_z[i][j][k] *= boxSize_inv2;
        potewald[i][j][k] *= boxSize_inv;
      }
    }
  }
#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Compute the Ewald correction for a given distance vector r.
 *
 * We interpolate the Ewald correction tables using a tri-linear interpolation
 * similar to a CIC.
 *
 * @param rx x-coordinate of distance vector.
 * @param ry y-coordinate of distance vector.
 * @param rz z-coordinate of distance vector.
 * @param corr_f (return) The Ewald correction for the force.
 * @param corr_p (return) The Ewald correction for the potential.
 */
void gravity_exact_force_ewald_evaluate(double rx, double ry, double rz,
                                        double corr_f[3], double *corr_p) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  const double s_x = (rx < 0.) ? 1. : -1.;
  const double s_y = (ry < 0.) ? 1. : -1.;
  const double s_z = (rz < 0.) ? 1. : -1.;
  rx = fabs(rx);
  ry = fabs(ry);
  rz = fabs(rz);

  int i = (int)(rx * ewald_fac);
  if (i >= Newald) i = Newald - 1;
  const double dx = rx * ewald_fac - i;
  const double tx = 1. - dx;

  int j = (int)(ry * ewald_fac);
  if (j >= Newald) j = Newald - 1;
  const double dy = ry * ewald_fac - j;
  const double ty = 1. - dy;

  int k = (int)(rz * ewald_fac);
  if (k >= Newald) k = Newald - 1;
  const double dz = rz * ewald_fac - k;
  const double tz = 1. - dz;

  /* Interpolation in X */
  corr_f[0] = 0.;
  corr_f[0] += fewald_x[i + 0][j + 0][k + 0] * tx * ty * tz;
  corr_f[0] += fewald_x[i + 0][j + 0][k + 1] * tx * ty * dz;
  corr_f[0] += fewald_x[i + 0][j + 1][k + 0] * tx * dy * tz;
  corr_f[0] += fewald_x[i + 0][j + 1][k + 1] * tx * dy * dz;
  corr_f[0] += fewald_x[i + 1][j + 0][k + 0] * dx * ty * tz;
  corr_f[0] += fewald_x[i + 1][j + 0][k + 1] * dx * ty * dz;
  corr_f[0] += fewald_x[i + 1][j + 1][k + 0] * dx * dy * tz;
  corr_f[0] += fewald_x[i + 1][j + 1][k + 1] * dx * dy * dz;
  corr_f[0] *= s_x;

  /* Interpolation in Y */
  corr_f[1] = 0.;
  corr_f[1] += fewald_y[i + 0][j + 0][k + 0] * tx * ty * tz;
  corr_f[1] += fewald_y[i + 0][j + 0][k + 1] * tx * ty * dz;
  corr_f[1] += fewald_y[i + 0][j + 1][k + 0] * tx * dy * tz;
  corr_f[1] += fewald_y[i + 0][j + 1][k + 1] * tx * dy * dz;
  corr_f[1] += fewald_y[i + 1][j + 0][k + 0] * dx * ty * tz;
  corr_f[1] += fewald_y[i + 1][j + 0][k + 1] * dx * ty * dz;
  corr_f[1] += fewald_y[i + 1][j + 1][k + 0] * dx * dy * tz;
  corr_f[1] += fewald_y[i + 1][j + 1][k + 1] * dx * dy * dz;
  corr_f[1] *= s_y;

  /* Interpolation in Z */
  corr_f[2] = 0.;
  corr_f[2] += fewald_z[i + 0][j + 0][k + 0] * tx * ty * tz;
  corr_f[2] += fewald_z[i + 0][j + 0][k + 1] * tx * ty * dz;
  corr_f[2] += fewald_z[i + 0][j + 1][k + 0] * tx * dy * tz;
  corr_f[2] += fewald_z[i + 0][j + 1][k + 1] * tx * dy * dz;
  corr_f[2] += fewald_z[i + 1][j + 0][k + 0] * dx * ty * tz;
  corr_f[2] += fewald_z[i + 1][j + 0][k + 1] * dx * ty * dz;
  corr_f[2] += fewald_z[i + 1][j + 1][k + 0] * dx * dy * tz;
  corr_f[2] += fewald_z[i + 1][j + 1][k + 1] * dx * dy * dz;
  corr_f[2] *= s_z;

  /* Interpolation for potential */
  *corr_p = 0.;
  *corr_p += potewald[i + 0][j + 0][k + 0] * tx * ty * tz;
  *corr_p += potewald[i + 0][j + 0][k + 1] * tx * ty * dz;
  *corr_p += potewald[i + 0][j + 1][k + 0] * tx * dy * tz;
  *corr_p += potewald[i + 0][j + 1][k + 1] * tx * dy * dz;
  *corr_p += potewald[i + 1][j + 0][k + 0] * dx * ty * tz;
  *corr_p += potewald[i + 1][j + 0][k + 1] * dx * ty * dz;
  *corr_p += potewald[i + 1][j + 1][k + 0] * dx * dy * tz;
  *corr_p += potewald[i + 1][j + 1][k + 1] * dx * dy * dz;

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Checks whether the file containing the exact accelerations for
 * the current choice of parameters already exists.
 *
 * @param e The #engine.
 */
int gravity_exact_force_file_exits(const struct engine *e) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  /* File name */
  char file_name[100];
  if (e->s->periodic)
    sprintf(file_name, "gravity_checks_exact_periodic_step%.4d.dat", e->step);
  else
    sprintf(file_name, "gravity_checks_exact_step%.4d.dat", e->step);

  /* Does the file exist ? */
  if (access(file_name, R_OK | W_OK) == 0) {

    /* Let's check whether the header matches the parameters of this run */
    FILE *file = fopen(file_name, "r");
    if (file == NULL) error("Problem reading gravity_check file");

    char line[100];
    char dummy1[10], dummy2[10];
    double newton_G;
    int N, periodic;
    /* Reads file header */
    if (fgets(line, 100, file) != line) error("Problem reading title");
    if (fgets(line, 100, file) != line) error("Problem reading G");
    sscanf(line, "%s %s %le", dummy1, dummy2, &newton_G);
    if (fgets(line, 100, file) != line) error("Problem reading N");
    sscanf(line, "%s %s %d", dummy1, dummy2, &N);
    if (fgets(line, 100, file) != line) error("Problem reading BC");
    sscanf(line, "%s %s %d", dummy1, dummy2, &periodic);
    fclose(file);

    /* Check whether it matches the current parameters */
    if (N == SWIFT_GRAVITY_FORCE_CHECKS && periodic == e->s->periodic &&
        (fabs(newton_G - e->physical_constants->const_newton_G) / newton_G <
         1e-5)) {
      return 1;
    }
  }
  return 0;
#else
  error("Gravity checking function called without the corresponding flag.");
  return 0;
#endif
}

/**
 * @brief Mapper function for the exact gravity calculation.
 */
void gravity_exact_force_compute_mapper(void *map_data, int nr_gparts,
                                        void *extra_data) {
#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  /* Unpack the data */
  struct gpart *restrict gparts = (struct gpart *)map_data;
  struct exact_force_data *data = (struct exact_force_data *)extra_data;
  const struct space *s = data->s;
  const struct part *parts = s->parts;
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double const_G = data->const_G;
  int counter = 0;

  for (int i = 0; i < nr_gparts; ++i) {

    struct gpart *gpi = &gparts[i];

    /* Get the particle ID */
    long long id = 0;
    if (gpi->type == swift_type_gas)
      id = parts[-gpi->id_or_neg_offset].id;
    else if (gpi->type == swift_type_stars)
      id = sparts[-gpi->id_or_neg_offset].id;
    else if (gpi->type == swift_type_black_hole)
      id = bparts[-gpi->id_or_neg_offset].id;
    else
      id = gpi->id_or_neg_offset;

    /* Is the particle active and part of the subset to be tested ? */
    if (id % SWIFT_GRAVITY_FORCE_CHECKS == 0 && gpart_is_active(gpi, e)) {

      /* Get some information about the particle */
      const double pix[3] = {gpi->x[0], gpi->x[1], gpi->x[2]};
      const double hi = gravity_get_softening(gpi, e->gravity_properties);
      const double hi_inv = 1. / hi;
      const double hi_inv3 = hi_inv * hi_inv * hi_inv;

      /* Be ready for the calculation */
      double a_grav[3] = {0., 0., 0.};
      double a_grav_short[3] = {0., 0., 0.};
      double a_grav_long[3] = {0., 0., 0.};
      double pot = 0.;

      /* Interact it with all other particles in the space.*/
      for (int j = 0; j < (int)s->nr_gparts; ++j) {

        const struct gpart *gpj = &s->gparts[j];

#ifdef SWIFT_DEBUG_CHECKS
        if (gpj->time_bin == time_bin_not_created) {
          error("Found an extra particle in the gravity check.");
        }
#endif

        /* No self interaction */
        if (gpi == gpj) continue;

        /* Compute the pairwise distance. */
        double dx = gpj->x[0] - pix[0];
        double dy = gpj->x[1] - pix[1];
        double dz = gpj->x[2] - pix[2];

        /* Now apply periodic BC */
        if (periodic) {
          dx = nearest(dx, dim[0]);
          dy = nearest(dy, dim[1]);
          dz = nearest(dz, dim[2]);
        }

        const double r2 = dx * dx + dy * dy + dz * dz;
        const double r_inv = 1. / sqrt(r2);
        const double r = r2 * r_inv;
        const double mj = gpj->mass;
        double f, phi;

        if (r >= hi) {

          /* Get Newtonian gravity */
          f = mj * r_inv * r_inv * r_inv;
          phi = -mj * r_inv;

        } else {

          const double ui = r * hi_inv;
          double Wf, Wp;

          kernel_grav_eval_force_double(ui, &Wf);
          kernel_grav_eval_pot_double(ui, &Wp);

          /* Get softened gravity */
          f = mj * hi_inv3 * Wf;
          phi = mj * hi_inv * Wp;
        }

        a_grav[0] += f * dx;
        a_grav[1] += f * dy;
        a_grav[2] += f * dz;
        pot += phi;

        /* Apply Ewald correction for periodic BC
         *
         * We also want to check what the tree and mesh do so we want to mimic
         * that:
         * - a_grav_short is the total acceleration multiplied by the
         * short-range correction.
         * - a_grav_long is the total acceleration (including Ewald correction)
         * minus the short-range acceleration.
         */
        if (periodic && r > 1e-5 * hi) {

          /* Compute trunctation for long and short range forces */
          const double r_s_inv = e->mesh->r_s_inv;
          const double u_lr = r * r_s_inv;
          double corr_f_lr;
          kernel_long_grav_force_eval_double(u_lr, &corr_f_lr);

          a_grav_short[0] += f * dx * corr_f_lr;
          a_grav_short[1] += f * dy * corr_f_lr;
          a_grav_short[2] += f * dz * corr_f_lr;

          a_grav_long[0] += f * dx * (1. - corr_f_lr);
          a_grav_long[1] += f * dy * (1. - corr_f_lr);
          a_grav_long[2] += f * dz * (1. - corr_f_lr);

          /* Ewald correction. */
          double corr_f[3], corr_pot;
          gravity_exact_force_ewald_evaluate(dx, dy, dz, corr_f, &corr_pot);

          a_grav[0] += mj * corr_f[0];
          a_grav[1] += mj * corr_f[1];
          a_grav[2] += mj * corr_f[2];
          pot += mj * corr_pot;

          a_grav_long[0] += mj * corr_f[0];
          a_grav_long[1] += mj * corr_f[1];
          a_grav_long[2] += mj * corr_f[2];
        }
      }

      /* Store the exact answer */
      for (int k = 0; k < 3; k++) {
        gpi->a_grav_exact[k] = a_grav[k] * const_G;
        gpi->a_grav_exact_short[k] = a_grav_short[k] * const_G;
        gpi->a_grav_exact_long[k] = a_grav_long[k] * const_G;
      }
      gpi->potential_exact = pot * const_G;

      counter++;
    }
  }
  atomic_add(&data->counter_global, counter);

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Run a brute-force gravity calculation for a subset of particles.
 *
 * All gpart with ID modulo SWIFT_GRAVITY_FORCE_CHECKS will get their forces
 * computed.
 *
 * @param s The #space to use.
 * @param e The #engine (to access the current time).
 */
void gravity_exact_force_compute(struct space *s, const struct engine *e) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  const ticks tic = getticks();

  /* Let's start by checking whether we already computed these forces */
  if (gravity_exact_force_file_exits(e)) {
    message("Exact accelerations already computed. Skipping calculation.");
    fflush(stdout);
    return;
  }

  /* No matching file present ? Do it then */
  struct exact_force_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;
  data.const_G = e->physical_constants->const_newton_G;

  threadpool_map(&s->e->threadpool, gravity_exact_force_compute_mapper,
                 s->gparts, s->nr_gparts, sizeof(struct gpart),
                 threadpool_auto_chunk_size, &data);

  message("Computed exact gravity for %d gparts (took %.3f %s). ",
          data.counter_global, clocks_from_ticks(getticks() - tic),
          clocks_getunit());

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Check the accuracy of the gravity calculation by comparing the
 * accelerations
 * to the brute-force computed ones.
 *
 * All gpart with ID modulo SWIFT_GRAVITY_FORCE_CHECKS will be checked.
 *
 * @param s The #space to use.
 * @param e The #engine (to access the current time).
 * @param rel_tol The maximal relative error. Will call error() if one particle
 * has a larger error.
 */
void gravity_exact_force_check(struct space *s, const struct engine *e,
                               float rel_tol) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  const struct part *parts = s->parts;
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  /* File name */
  char file_name_swift[100];
  sprintf(file_name_swift, "gravity_checks_swift_step%.4d_order%d.dat", e->step,
          SELF_GRAVITY_MULTIPOLE_ORDER);

  /* Creare files and write header */
  FILE *file_swift = fopen(file_name_swift, "w");
  if (file_swift == NULL) error("Could not create file '%s'.", file_name_swift);
  fprintf(file_swift, "# Gravity accuracy test - SWIFT FORCES\n");
  fprintf(file_swift, "# G= %16.8e\n", e->physical_constants->const_newton_G);
  fprintf(file_swift, "# N= %d\n", SWIFT_GRAVITY_FORCE_CHECKS);
  fprintf(file_swift, "# periodic= %d\n", s->periodic);
  fprintf(file_swift, "# theta= %16.8e\n", e->gravity_properties->theta_crit);
  fprintf(file_swift, "# Git Branch: %s\n", git_branch());
  fprintf(file_swift, "# Git Revision: %s\n", git_revision());
  fprintf(file_swift,
          "# %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s "
          "%16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "a_swift[0]", "a_swift[1]",
          "a_swift[2]", "potential", "a_PM[0]", "a_PM[1]", "a_PM[2]",
          "potentialPM", "a_p2p[0]", "a_p2p[1]", "a_p2p[2]", "a_m2p[0]",
          "a_m2p[1]", "a_m2p[2]", "a_m2l[0]", "a_m2l[1]", "a_m2l[2]", "n_p2p",
          "n_m2p", "n_m2l", "n_PM");

  /* Output particle SWIFT accelerations  */
  for (size_t i = 0; i < s->nr_gparts; ++i) {

    struct gpart *gpi = &s->gparts[i];

    /* Get the particle ID */
    long long id = 0;
    if (gpi->type == swift_type_gas)
      id = parts[-gpi->id_or_neg_offset].id;
    else if (gpi->type == swift_type_stars)
      id = sparts[-gpi->id_or_neg_offset].id;
    else if (gpi->type == swift_type_black_hole)
      id = bparts[-gpi->id_or_neg_offset].id;
    else
      id = gpi->id_or_neg_offset;

    /* Is the particle was active and part of the subset to be tested ? */
    if (id % SWIFT_GRAVITY_FORCE_CHECKS == 0 && gpart_is_starting(gpi, e)) {

      fprintf(file_swift,
              "%18lld %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e "
              "%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e "
              "%16.8e %16.8e %16.8e %16lld %16lld %16lld %16lld\n",
              id, gpi->x[0], gpi->x[1], gpi->x[2],
              gpi->a_grav[0] + gpi->a_grav_mesh[0],
              gpi->a_grav[1] + gpi->a_grav_mesh[1],
              gpi->a_grav[2] + gpi->a_grav_mesh[2],
              gravity_get_comoving_potential(gpi), gpi->a_grav_mesh[0],
              gpi->a_grav_mesh[1], gpi->a_grav_mesh[2],
              gravity_get_comoving_mesh_potential(gpi), gpi->a_grav_p2p[0],
              gpi->a_grav_p2p[1], gpi->a_grav_p2p[2], gpi->a_grav_m2p[0],
              gpi->a_grav_m2p[1], gpi->a_grav_m2p[2], gpi->a_grav_m2l[0],
              gpi->a_grav_m2l[1], gpi->a_grav_m2l[2], gpi->num_interacted_p2p,
              gpi->num_interacted_m2p, gpi->num_interacted_m2l,
              gpi->num_interacted_pm);
    }
  }

  message("Written SWIFT accelerations in file '%s'.", file_name_swift);

  /* Be nice */
  fclose(file_swift);

  if (!gravity_exact_force_file_exits(e)) {

    char file_name_exact[100];
    if (s->periodic)
      sprintf(file_name_exact, "gravity_checks_exact_periodic_step%.4d.dat",
              e->step);
    else
      sprintf(file_name_exact, "gravity_checks_exact_step%.4d.dat", e->step);

    FILE *file_exact = fopen(file_name_exact, "w");
    if (file_exact == NULL)
      error("Could not create file '%s'.", file_name_exact);
    fprintf(file_exact, "# Gravity accuracy test - EXACT FORCES\n");
    fprintf(file_exact, "# G= %16.8e\n", e->physical_constants->const_newton_G);
    fprintf(file_exact, "# N= %d\n", SWIFT_GRAVITY_FORCE_CHECKS);
    fprintf(file_exact, "# periodic= %d\n", s->periodic);
    fprintf(file_exact, "# Git Branch: %s\n", git_branch());
    fprintf(file_exact, "# Git Revision: %s\n", git_revision());
    fprintf(file_exact,
            "# %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s "
            "%16s %16s %16s\n",
            "id", "pos[0]", "pos[1]", "pos[2]", "a_exact[0]", "a_exact[1]",
            "a_exact[2]", "potential", "a_exact_short[0]", "a_exact_short[1]",
            "a_exact_short[2]", "a_exact_long[0]", "a_exact_long[1]",
            "a_exact_long[2]");

    /* Output particle exact accelerations  */
    for (size_t i = 0; i < s->nr_gparts; ++i) {

      struct gpart *gpi = &s->gparts[i];

      long long id = 0;
      if (gpi->type == swift_type_gas)
        id = parts[-gpi->id_or_neg_offset].id;
      else if (gpi->type == swift_type_stars)
        id = sparts[-gpi->id_or_neg_offset].id;
      else if (gpi->type == swift_type_black_hole)
        id = bparts[-gpi->id_or_neg_offset].id;
      else
        id = gpi->id_or_neg_offset;

      /* Is the particle was active and part of the subset to be tested ? */
      if (id % SWIFT_GRAVITY_FORCE_CHECKS == 0 && gpart_is_starting(gpi, e)) {
        fprintf(file_exact,
                "%18lld %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e "
                "%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
                id, gpi->x[0], gpi->x[1], gpi->x[2], gpi->a_grav_exact[0],
                gpi->a_grav_exact[1], gpi->a_grav_exact[2],
                gpi->potential_exact, gpi->a_grav_exact_short[0],
                gpi->a_grav_exact_short[1], gpi->a_grav_exact_short[2],
                gpi->a_grav_exact_long[0], gpi->a_grav_exact_long[1],
                gpi->a_grav_exact_long[2]);
      }
    }

    message("Written exact accelerations in file '%s'.", file_name_exact);

    /* Be nice */
    fclose(file_exact);
  }
#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}
