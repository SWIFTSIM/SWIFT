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
void gravity_exact_force_ewald_init(double boxSize) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  const ticks tic = getticks();
  message("Computing Ewald correction table...");

  /* Level of correction  (Hernquist et al. 1991)*/
  const float alpha = 2.f;

  /* some useful constants */
  const float alpha2 = alpha * alpha;
  const float factor_exp1 = 2.f * alpha / sqrt(M_PI);
  const float factor_exp2 = -M_PI * M_PI / alpha2;
  const float factor_sin = 2.f * M_PI;
  const float boxSize_inv2 = 1.f / (boxSize * boxSize);

  /* Ewald factor to access the table */
  ewald_fac = (double)(2 * Newald) / boxSize;

  /* Zero everything */
  bzero(fewald_x, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));
  bzero(fewald_y, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));
  bzero(fewald_z, (Newald + 1) * (Newald + 1) * (Newald + 1) * sizeof(float));

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

              const float val =
                  erfcf(alpha * r_tilde) +
                  factor_exp1 * r_tilde * expf(-alpha2 * r_tilde2);

              /* First correction term */
              const float f = val * r_tilde_inv3;
              f_x -= f * d_x;
              f_y -= f * d_y;
              f_z -= f * d_z;
            }
          }
        }

        for (int h_i = -4; h_i <= 4; ++h_i) {
          for (int h_j = -4; h_j <= 4; ++h_j) {
            for (int h_k = -4; h_k <= 4; ++h_k) {

              const float h2 = h_i * h_i + h_j * h_j + h_k * h_k;

              const float h2_inv = 1.f / (h2 + FLT_MIN);
              const float h_dot_x = h_i * r_x + h_j * r_y + h_k * r_z;

              const float val = 2.f * h2_inv * expf(h2 * factor_exp2) *
                                sinf(factor_sin * h_dot_x);

              /* Second correction term */
              f_x -= val * h_i;
              f_y -= val * h_j;
              f_z -= val * h_k;
            }
          }
        }

        /* Save back to memory */
        fewald_x[i][j][k] = f_x;
        fewald_y[i][j][k] = f_y;
        fewald_z[i][j][k] = f_z;
      }
    }
  }

/* Dump the Ewald table to a file */
#ifdef HAVE_HDF5
  hid_t h_file =
      H5Fcreate("Ewald.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file for Ewald dump.");

  /* Create dataspace */
  hsize_t dim[3] = {Newald + 1, Newald + 1, Newald + 1};
  hid_t h_space = H5Screate_simple(3, dim, NULL);
  hid_t h_data;
  h_data = H5Dcreate(h_file, "Ewald_x", H5T_NATIVE_FLOAT, h_space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
           &(fewald_x[0][0][0]));
  H5Dclose(h_data);
  h_data = H5Dcreate(h_file, "Ewald_y", H5T_NATIVE_FLOAT, h_space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
           &(fewald_y[0][0][0]));
  H5Dclose(h_data);
  h_data = H5Dcreate(h_file, "Ewald_z", H5T_NATIVE_FLOAT, h_space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(h_data, H5T_NATIVE_FLOAT, h_space, H5S_ALL, H5P_DEFAULT,
           &(fewald_z[0][0][0]));
  H5Dclose(h_data);
  H5Sclose(h_space);
  H5Fclose(h_file);
#endif

  /* Apply the box-size correction */
  for (int i = 0; i <= Newald; ++i) {
    for (int j = 0; j <= Newald; ++j) {
      for (int k = 0; k <= Newald; ++k) {
        fewald_x[i][j][k] *= boxSize_inv2;
        fewald_y[i][j][k] *= boxSize_inv2;
        fewald_z[i][j][k] *= boxSize_inv2;
      }
    }
  }

  /* Say goodbye */
  message("Ewald correction table computed (took %.3f %s). ",
          clocks_from_ticks(getticks() - tic), clocks_getunit());
#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
/**
 * @brief Compute the Ewald correction for a given distance vector r.
 *
 * We interpolate the Ewald correction tables using a tri-linear interpolation
 * similar to a CIC.
 *
 * @param rx x-coordinate of distance vector.
 * @param ry y-coordinate of distance vector.
 * @param rz z-coordinate of distance vector.
 * @param corr (return) The Ewald correction.
 */
__attribute__((always_inline)) INLINE static void
gravity_exact_force_ewald_evaluate(double rx, double ry, double rz,
                                   double corr[3]) {

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
  corr[0] = 0.;
  corr[0] += fewald_x[i + 0][j + 0][k + 0] * tx * ty * tz;
  corr[0] += fewald_x[i + 0][j + 0][k + 1] * tx * ty * dz;
  corr[0] += fewald_x[i + 0][j + 1][k + 0] * tx * dy * tz;
  corr[0] += fewald_x[i + 0][j + 1][k + 1] * tx * dy * dz;
  corr[0] += fewald_x[i + 1][j + 0][k + 0] * dx * ty * tz;
  corr[0] += fewald_x[i + 1][j + 0][k + 1] * dx * ty * dz;
  corr[0] += fewald_x[i + 1][j + 1][k + 0] * dx * dy * tz;
  corr[0] += fewald_x[i + 1][j + 1][k + 1] * dx * dy * dz;
  corr[0] *= s_x;

  /* Interpolation in Y */
  corr[1] = 0.;
  corr[1] += fewald_y[i + 0][j + 0][k + 0] * tx * ty * tz;
  corr[1] += fewald_y[i + 0][j + 0][k + 1] * tx * ty * dz;
  corr[1] += fewald_y[i + 0][j + 1][k + 0] * tx * dy * tz;
  corr[1] += fewald_y[i + 0][j + 1][k + 1] * tx * dy * dz;
  corr[1] += fewald_y[i + 1][j + 0][k + 0] * dx * ty * tz;
  corr[1] += fewald_y[i + 1][j + 0][k + 1] * dx * ty * dz;
  corr[1] += fewald_y[i + 1][j + 1][k + 0] * dx * dy * tz;
  corr[1] += fewald_y[i + 1][j + 1][k + 1] * dx * dy * dz;
  corr[1] *= s_y;

  /* Interpolation in Z */
  corr[2] = 0.;
  corr[2] += fewald_z[i + 0][j + 0][k + 0] * tx * ty * tz;
  corr[2] += fewald_z[i + 0][j + 0][k + 1] * tx * ty * dz;
  corr[2] += fewald_z[i + 0][j + 1][k + 0] * tx * dy * tz;
  corr[2] += fewald_z[i + 0][j + 1][k + 1] * tx * dy * dz;
  corr[2] += fewald_z[i + 1][j + 0][k + 0] * dx * ty * tz;
  corr[2] += fewald_z[i + 1][j + 0][k + 1] * dx * ty * dz;
  corr[2] += fewald_z[i + 1][j + 1][k + 0] * dx * dy * tz;
  corr[2] += fewald_z[i + 1][j + 1][k + 1] * dx * dy * dz;
  corr[2] *= s_z;
}
#endif

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
    sprintf(file_name, "gravity_checks_exact_periodic_step%d.dat", e->step);
  else
    sprintf(file_name, "gravity_checks_exact_step%d.dat", e->step);

  /* Does the file exist ? */
  if (access(file_name, R_OK | W_OK) == 0) {

    /* Let's check whether the header matches the parameters of this run */
    FILE *file = fopen(file_name, "r");
    if (!file) error("Problem reading gravity_check file");

    char line[100];
    char dummy1[10], dummy2[10];
    double epsilon, newton_G;
    int N, periodic;
    /* Reads file header */
    if (fgets(line, 100, file) != line) error("Problem reading title");
    if (fgets(line, 100, file) != line) error("Problem reading G");
    sscanf(line, "%s %s %le", dummy1, dummy2, &newton_G);
    if (fgets(line, 100, file) != line) error("Problem reading N");
    sscanf(line, "%s %s %d", dummy1, dummy2, &N);
    if (fgets(line, 100, file) != line) error("Problem reading epsilon");
    sscanf(line, "%s %s %le", dummy1, dummy2, &epsilon);
    if (fgets(line, 100, file) != line) error("Problem reading BC");
    sscanf(line, "%s %s %d", dummy1, dummy2, &periodic);
    fclose(file);

    /* Check whether it matches the current parameters */
    if (N == SWIFT_GRAVITY_FORCE_CHECKS && periodic == e->s->periodic &&
        (fabs(epsilon - e->gravity_properties->epsilon) / epsilon < 1e-5) &&
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
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double const_G = data->const_G;
  int counter = 0;

  for (int i = 0; i < nr_gparts; ++i) {

    struct gpart *gpi = &gparts[i];

    /* Is the particle active and part of the subset to be tested ? */
    if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
        gpart_is_active(gpi, e)) {

      /* Get some information about the particle */
      const double pix[3] = {gpi->x[0], gpi->x[1], gpi->x[2]};
      const double hi = gpi->epsilon;
      const double hi_inv = 1. / hi;
      const double hi_inv3 = hi_inv * hi_inv * hi_inv;

      /* Be ready for the calculation */
      double a_grav[3] = {0., 0., 0.};

      /* Interact it with all other particles in the space.*/
      for (int j = 0; j < (int)s->nr_gparts; ++j) {

        struct gpart *gpj = &s->gparts[j];

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
        double f;

        if (r >= hi) {

          /* Get Newtonian gravity */
          f = mj * r_inv * r_inv * r_inv;

        } else {

          const double ui = r * hi_inv;
          double W;

          kernel_grav_eval_double(ui, &W);

          /* Get softened gravity */
          f = mj * hi_inv3 * W;
        }

        a_grav[0] += f * dx;
        a_grav[1] += f * dy;
        a_grav[2] += f * dz;

        /* Apply Ewald correction for periodic BC */
        if (periodic && r > 1e-5 * hi) {

          double corr[3];
          gravity_exact_force_ewald_evaluate(dx, dy, dz, corr);

          a_grav[0] += mj * corr[0];
          a_grav[1] += mj * corr[1];
          a_grav[2] += mj * corr[2];
        }
      }

      /* Store the exact answer */
      gpi->a_grav_exact[0] = a_grav[0] * const_G;
      gpi->a_grav_exact[1] = a_grav[1] * const_G;
      gpi->a_grav_exact[2] = a_grav[2] * const_G;

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
    return;
  }

  /* No matching file present ? Do it then */
  struct exact_force_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;
  data.const_G = e->physical_constants->const_newton_G;

  threadpool_map(&s->e->threadpool, gravity_exact_force_compute_mapper,
                 s->gparts, s->nr_gparts, sizeof(struct gpart), 0, &data);

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

  /* File name */
  char file_name_swift[100];
  sprintf(file_name_swift, "gravity_checks_swift_step%d_order%d.dat", e->step,
          SELF_GRAVITY_MULTIPOLE_ORDER);

  /* Creare files and write header */
  FILE *file_swift = fopen(file_name_swift, "w");
  fprintf(file_swift, "# Gravity accuracy test - SWIFT FORCES\n");
  fprintf(file_swift, "# G= %16.8e\n", e->physical_constants->const_newton_G);
  fprintf(file_swift, "# N= %d\n", SWIFT_GRAVITY_FORCE_CHECKS);
  fprintf(file_swift, "# epsilon= %16.8e\n", e->gravity_properties->epsilon);
  fprintf(file_swift, "# periodic= %d\n", s->periodic);
  fprintf(file_swift, "# theta= %16.8e\n", e->gravity_properties->theta_crit);
  fprintf(file_swift, "# Git Branch: %s\n", git_branch());
  fprintf(file_swift, "# Git Revision: %s\n", git_revision());
  fprintf(file_swift, "# %16s %16s %16s %16s %16s %16s %16s\n", "id", "pos[0]",
          "pos[1]", "pos[2]", "a_swift[0]", "a_swift[1]", "a_swift[2]");

  /* Output particle SWIFT accelerations  */
  for (size_t i = 0; i < s->nr_gparts; ++i) {

    struct gpart *gpi = &s->gparts[i];

    /* Is the particle was active and part of the subset to be tested ? */
    if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
        gpart_is_starting(gpi, e)) {

      fprintf(file_swift, "%18lld %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n",
              gpi->id_or_neg_offset, gpi->x[0], gpi->x[1], gpi->x[2],
              gpi->a_grav[0], gpi->a_grav[1], gpi->a_grav[2]);
    }
  }

  message("Written SWIFT accelerations in file '%s'.", file_name_swift);

  /* Be nice */
  fclose(file_swift);

  if (!gravity_exact_force_file_exits(e)) {

    char file_name_exact[100];
    if (s->periodic)
      sprintf(file_name_exact, "gravity_checks_exact_periodic_step%d.dat",
              e->step);
    else
      sprintf(file_name_exact, "gravity_checks_exact_step%d.dat", e->step);

    FILE *file_exact = fopen(file_name_exact, "w");
    fprintf(file_exact, "# Gravity accuracy test - EXACT FORCES\n");
    fprintf(file_exact, "# G= %16.8e\n", e->physical_constants->const_newton_G);
    fprintf(file_exact, "# N= %d\n", SWIFT_GRAVITY_FORCE_CHECKS);
    fprintf(file_exact, "# epsilon=%16.8e\n", e->gravity_properties->epsilon);
    fprintf(file_exact, "# periodic= %d\n", s->periodic);
    fprintf(file_exact, "# Git Branch: %s\n", git_branch());
    fprintf(file_exact, "# Git Revision: %s\n", git_revision());
    fprintf(file_exact, "# %16s %16s %16s %16s %16s %16s %16s\n", "id",
            "pos[0]", "pos[1]", "pos[2]", "a_exact[0]", "a_exact[1]",
            "a_exact[2]");

    /* Output particle exact accelerations  */
    for (size_t i = 0; i < s->nr_gparts; ++i) {

      struct gpart *gpi = &s->gparts[i];

      /* Is the particle was active and part of the subset to be tested ? */
      if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
          gpart_is_starting(gpi, e)) {

        fprintf(
            file_exact, "%18lld %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n",
            gpi->id_or_neg_offset, gpi->x[0], gpi->x[1], gpi->x[2],
            gpi->a_grav_exact[0], gpi->a_grav_exact[1], gpi->a_grav_exact[2]);
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
