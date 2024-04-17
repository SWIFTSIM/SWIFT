/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* local headers */
#include "initial_mass_function.h"

#include "hdf5_functions.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Get the IMF exponent in between mass_min and mass_max.
 */
float initial_mass_function_get_exponent(
    const struct initial_mass_function *imf, float mass_min, float mass_max) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mass_max > imf->mass_max)
    error("Cannot have mass larger than the largest one in the IMF");
  if (mass_min < imf->mass_min)
    error("Cannot have mass smaller than the smallest one in the IMF");
  if (mass_max < mass_min) error("Cannot have mass_min larger than mass_max");
#endif

  for (int i = 0; i < imf->n_parts; i++) {

    /* Check if in the correct part of the IMF */
    if (mass_min < imf->mass_limits[i + 1]) {

      /* Check if in only one segment */
      if (mass_max > imf->mass_limits[i + 1]) {
        error("Cannot get a single exponent for the interval [%g, %g]",
              mass_min, mass_max);
      }

      return imf->exp[i];
    }
  }

  error("Masses outside IMF ranges");

  return -1;
}

/** @brief Print the initial mass function */
void initial_mass_function_print(const struct initial_mass_function *imf) {

  message("Number of parts: %i", imf->n_parts);
  message("Number of stars per mass units: %g", imf->N_tot);
  message("Mass interval: [%g, %g]", imf->mass_min, imf->mass_max);
  for (int i = 0; i < imf->n_parts; i++) {
    message("[%7.3f, %7.3f]: %5.2g * m^{%g}", imf->mass_limits[i],
            imf->mass_limits[i + 1], imf->coef[i], imf->exp[i]);
  }

  message("Mass fractions");
  for (int i = 0; i < imf->n_parts + 1; i++)
    message("m=%7.3f: x=%5.3f", imf->mass_limits[i], imf->mass_fraction[i]);
}

/** @brief Sample the initial mass function */
float initial_mass_function_sample(const struct initial_mass_function *imf,
                                   float f) {

  for (int i = 0; i < imf->n_parts; i++)
    if (f < imf->mass_fraction[i + 1]) {
      float pmin = pow(imf->mass_limits[i], imf->exp[i]);
      float exponent = 1. / imf->exp[i];
      float base_part_1 = imf->N_tot * imf->exp[i] / imf->coef[i];
      base_part_1 *= (f - imf->mass_fraction[i]);
      float base = base_part_1 + pmin;

      /* The mathematical expression is:
           (N_tot * exp_{imf, i} / coeff_{imf, i} * (f - mass_fraction_{imf, i})
             + pmin)**(1/exp_{imf, i})
      */
      return pow(base, exponent);
    }

  return -1;
}

/**
 * @brief Integrate the #interpolation_1d data with the initial mass function.
 *
 * The x are supposed to be linear in log.
 *
 * @param imf The #initial_mass_function.
 * @param data The data to integrate.
 * @param count The number of element in data.
 * @param log_mass_min The value of the first element.
 * @param step_size The distance between two points.
 */
void initial_mass_function_integrate(const struct initial_mass_function *imf,
                                     float *data, size_t count,
                                     float log_mass_min, float step_size) {

  /* Index in the data */
  size_t j = 1;
  const float mass_min = exp10(log_mass_min);
  const float mass_max = exp10(log_mass_min + (count - 1) * step_size);

  float m = mass_min;

  float *tmp = (float *)malloc(sizeof(float) * count);

  /* Set lower limit */
  tmp[0] = 0;
  for (int i = 0; i < imf->n_parts; i++) {

    /* Check if already in the correct part */
    if (mass_min > imf->mass_limits[i + 1]) {
      continue;
    }

    /* Check if already above the maximal mass */
    if (mass_max < imf->mass_limits[i]) {
      break;
    }

    /* Integrate the data */
    while ((m < imf->mass_limits[i + 1] || i == imf->n_parts - 1) &&
           j < count) {

      /* Compute the masses */
      const float log_m1 = log_mass_min + (j - 1) * step_size;
      const float m1 = exp10(log_m1);
      const float log_m2 = log_mass_min + j * step_size;
      float m2 = exp10(log_m2);

      /* Ensure that we stay within the limits */
      if (m2 > imf->mass_max) {
        m2 = imf->mass_max;
      }

      const float dm = m2 - m1;
      const float imf_1 = imf->coef[i] * pow(m1, imf->exp[i]);

      /* Get the imf of the upper limit  */
      float imf_2;
      if (m2 > imf->mass_limits[i + 1]) {
        imf_2 = imf->coef[i + 1] * pow(m2, imf->exp[i + 1]);
      } else {
        imf_2 = imf->coef[i] * pow(m2, imf->exp[i]);
      }

      /* Compute the integral */
      tmp[j] = tmp[j - 1] + 0.5 * (imf_1 * data[j - 1] + imf_2 * data[j]) * dm;

      /* Update j and m */
      j += 1;
      m = m2;
    }
  }

  /* The rest is extrapolated with 0 */
  for (size_t k = j; k < count; k++) {
    tmp[k] = tmp[k - 1];
  }

  /* Copy temporary array */
  memcpy(data, tmp, count * sizeof(float));

  /* clean everything */
  free(tmp);
}

/**
 * @brief Get the IMF coefficient in between mass_min and mass_max.
 *
 * @param imf The #initial_mass_function.
 * @param mass_min The minimal mass of the requested interval.
 * @param mass_max The maximal mass of the requested interval.
 *
 * @return The imf's coefficient of the interval.
 */
float initial_mass_function_get_coefficient(
    const struct initial_mass_function *imf, float mass_min, float mass_max) {

  for (int i = 0; i < imf->n_parts; i++) {

    /* Check if in the correct part of the IMF */
    if (mass_min < imf->mass_limits[i + 1]) {

      /* Check if in only one segment */
      if (mass_max > imf->mass_limits[i + 1]) {
        error("Cannot get a single coefficient for the interval [%g, %g]",
              mass_min, mass_max);
      }

      return imf->coef[i];
    }
  }

  error("Masses outside IMF ranges");

  return -1;
}

/**
 * @brief Compute the integral of the fraction number of the initial mass
 * function.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lower mass to evaluate.
 * @param m2 The upper mass to evaluate.
 *
 * @return The number fraction.
 */
float initial_mass_function_get_integral_xi(
    const struct initial_mass_function *imf, float m1, float m2) {

  /* Ensure the masses to be withing the limits */
  m1 = min(m1, imf->mass_max);
  m1 = max(m1, imf->mass_min);

  m2 = min(m2, imf->mass_max);
  m2 = max(m2, imf->mass_min);

  int k = -1;
  /* Find the correct part */
  for (int i = 0; i < imf->n_parts; i++) {
    if (m1 <= imf->mass_limits[i + 1]) {
      k = i;
      break;
    }
  }

  /* Check if found a part */
  if (k == -1) {
    error("Failed to find the correct function part: %g %g", m1, m2);
  }

  /* Check if m2 is inside the part */
  if (m2 < imf->mass_limits[k] || m2 > imf->mass_limits[k + 1]) {
    error("This function is not able to integrate in two different parts %g %g",
          m1, m2);
  }

  /* Compute the integral */
  const float int_xi1 = pow(m1, imf->exp[k]);
  const float int_xi2 = pow(m2, imf->exp[k]);

  return imf->coef[k] * (int_xi2 - int_xi1) / imf->exp[k];
};

/**
 * @brief Compute the mass fraction of the initial mass function.
 *
 * @param imf The #initial_mass_function.
 * @param m The mass to evaluate.
 *
 * @return The mass fraction.
 */
float initial_mass_function_get_imf(const struct initial_mass_function *imf,
                                    float m) {

  /* Check the mass to be within the limits */

  if (m > imf->mass_max || m < imf->mass_min)
    error("Mass below or above limits expecting %g < %g < %g.", imf->mass_min,
          m, imf->mass_max);

  for (int i = 0; i < imf->n_parts; i++) {
    if (m <= imf->mass_limits[i + 1]) {
      return imf->coef[i] * pow(m, imf->exp[i]);
    }
  }

  error("Failed to find correct function part: %g larger than mass max %g.", m,
        imf->mass_max);
  return 0.;
};

/**
 * @brief Compute the the mass fraction (of stars) between m1 and m2 per mass
 * unit.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lower mass to evaluate.
 * @param m2 The upper mass to evaluate.
 *
 * @return The integral of the mass fraction.
 */
float initial_mass_function_get_imf_mass_fraction(
    const struct initial_mass_function *imf, float m1, float m2) {

  /* Check that m2 is > m1 */
  if (m1 > m2)
    error("Mass m1 (=%g) larger or equal to m2 (=%g). This is not allowed", m1,
          m2);

  /* Check the masses to be within the limits */

  if (m1 > imf->mass_max || m1 < imf->mass_min)
    error("Mass m1 below or above limits expecting %g < %g < %g.",
          imf->mass_min, m1, imf->mass_max);

  if (m2 > imf->mass_max || m2 < imf->mass_min)
    error("Mass m2 below or above limits expecting %g < %g < %g.",
          imf->mass_min, m2, imf->mass_max);

  const int n = imf->n_parts;
  float integral = 0;

  /* loop over all segments */
  for (int i = 0; i < n; i++) {
    float mmin = max(imf->mass_limits[i], m1);
    float mmax = min(imf->mass_limits[i + 1], m2);

    if (mmin < mmax) {
      float p = imf->exp[i] + 1;
      integral += (imf->coef[i] / p) * (pow(mmax, p) - pow(mmin, p));
    } else /* nothing in this segment go to the next one */
      continue;

    if (m2 == mmax) /* nothing after this segment, stop */
      break;
  }

  return integral;
};

/**
 * @brief Compute the number fraction (of stars) between m1 and m2 per mass
 * unit.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lower mass to evaluate.
 * @param m2 The upper mass to evaluate.
 *
 * @return The integral of the mass fraction.
 */
float initial_mass_function_get_imf_number_fraction(
    const struct initial_mass_function *imf, float m1, float m2) {

  /* Check that m2 is > m1 */
  if (m1 > m2)
    error("Mass m1 (=%g) larger or equal to m2 (=%g). This is not allowed", m1,
          m2);

  /* Check the masses to be within the limits */

  if (m1 > imf->mass_max || m1 < imf->mass_min)
    error("Mass m1 below or above limits expecting %g < %g < %g.",
          imf->mass_min, m1, imf->mass_max);

  if (m2 > imf->mass_max || m2 < imf->mass_min)
    error("Mass m2 below or above limits expecting %g < %g < %g.",
          imf->mass_min, m2, imf->mass_max);

  const int n = imf->n_parts;
  float integral = 0;

  /* loop over all segments */
  for (int i = 0; i < n; i++) {
    float mmin = max(imf->mass_limits[i], m1);
    float mmax = min(imf->mass_limits[i + 1], m2);

    if (mmin < mmax) {
      float p = imf->exp[i];
      integral += (imf->coef[i] / p) * (pow(mmax, p) - pow(mmin, p));
    } else /* nothing in this segment go to the next one */
      continue;

    if (m2 == mmax) /* nothing after this segment, stop */
      break;
  }

  return integral;
};

/**
 * @brief Compute the coefficients of the initial mass function
 * as well as the mass fraction at the interface between IMF segments.
 *
 * @param imf The #initial_mass_function.
 */
void initial_mass_function_compute_coefficients(
    struct initial_mass_function *imf) {

  /* Allocate memory */
  if ((imf->coef = (float *)malloc(sizeof(float) * imf->n_parts)) == NULL)
    error("Failed to allocate the IMF coefficients.");

  /* Suppose that the first coefficients is 1 (will be corrected later) */
  imf->coef[0] = 1.;

  /* Use the criterion of continuity for the IMF */
  for (int i = 1; i < imf->n_parts; i++) {
    float exp = imf->exp[i - 1] - imf->exp[i];
    imf->coef[i] = imf->coef[i - 1] * pow(imf->mass_limits[i], exp);
  }

  /* Use the criterion on the integral = 1 */
  float integral = 0;
  for (int i = 0; i < imf->n_parts; i++) {
    const float exp = imf->exp[i] + 1.;
    const float m_i = pow(imf->mass_limits[i], exp);
    const float m_i1 = pow(imf->mass_limits[i + 1], exp);
    integral += imf->coef[i] * (m_i1 - m_i) / exp;
  }

  /* Normalize the coefficients (fix initial supposition) */
  for (int i = 0; i < imf->n_parts; i++) {
    imf->coef[i] /= integral;
  }

  /* Compute the total number of stars per mass unit */
  imf->N_tot = initial_mass_function_get_imf_number_fraction(imf, imf->mass_min,
                                                             imf->mass_max);

  /* Allocate the memory for the mass fraction */
  if ((imf->mass_fraction =
           (float *)malloc(sizeof(float) * (imf->n_parts + 1))) == NULL)
    error("Failed to allocate the IMF mass_fraction.");

  for (int i = 0; i < imf->n_parts + 1; i++) {
    imf->mass_fraction[i] = initial_mass_function_get_imf_number_fraction(
                                imf, imf->mass_min, imf->mass_limits[i]) /
                            imf->N_tot;
  }
}

/**
 * @brief Reads the initial mass function parameters from the tables.
 *
 * @param imf The #initial_mass_function.
 * @param params The #swift_params.
 * @param filename The filename of the chemistry table.
 */
void initial_mass_function_read_from_table(struct initial_mass_function *imf,
                                           struct swift_params *params,
                                           const char *filename) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(filename, "Data/IMF", &file_id, &group_id);

  /* Read number of parts */
  io_read_attribute(group_id, "n", INT, &imf->n_parts);

  /* The tables have a different definition of n */
  imf->n_parts += 1;

  /* Allocate the memory for the exponents */
  if ((imf->exp = (float *)malloc(sizeof(float) * imf->n_parts)) == NULL)
    error("Failed to allocate the IMF exponents.");

  /* Read the exponents */
  io_read_array_attribute(group_id, "as", FLOAT, imf->exp, imf->n_parts);

  /* Allocate the memory for the temporary mass limits */
  if ((imf->mass_limits =
           (float *)malloc(sizeof(float) * (imf->n_parts + 1))) == NULL)
    error("Failed to allocate the IMF masses.");

  /* Read the mass limits */
  io_read_array_attribute(group_id, "ms", FLOAT, imf->mass_limits,
                          imf->n_parts - 1);

  /* Copy the data (need to shift for mass_min) */
  for (int i = imf->n_parts - 1; i > 0; i--) {
    imf->mass_limits[i] = imf->mass_limits[i - 1];
  }

  /* Read the minimal mass limit */
  io_read_attribute(group_id, "Mmin", FLOAT, &imf->mass_limits[0]);

  /* Read the maximal mass limit */
  io_read_attribute(group_id, "Mmax", FLOAT, &imf->mass_limits[imf->n_parts]);

  /* Close everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Initialize the initial mass function.
 *
 * @param imf The #initial_mass_function.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The #swift_params.
 * @param filename The filename of the chemistry table.
 */
void initial_mass_function_init(struct initial_mass_function *imf,
                                const struct phys_const *phys_const,
                                const struct unit_system *us,
                                struct swift_params *params,
                                const char *filename) {

  /* Read the parameters from the yields table */
  initial_mass_function_read_from_table(imf, params, filename);

  /* Write the masses in the correct attributes */
  imf->mass_min = imf->mass_limits[0];
  imf->mass_max = imf->mass_limits[imf->n_parts];

  /* Compute the coefficients */
  initial_mass_function_compute_coefficients(imf);

  /* Print info */
  initial_mass_function_print(imf);
}

/**
 * @brief Write a initial_mass_function struct to the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param imf the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void initial_mass_function_dump(const struct initial_mass_function *imf,
                                FILE *stream, const struct stellar_model *sm) {

  /* Dump the mass limits. */
  if (imf->mass_limits != NULL) {
    restart_write_blocks((void *)imf->mass_limits, sizeof(float),
                         imf->n_parts + 1, stream, "imf_mass_limits",
                         "imf_mass_limits");
  }

  /*! Dump the exponents. */
  if (imf->exp != NULL) {
    restart_write_blocks((void *)imf->exp, sizeof(float), imf->n_parts, stream,
                         "imf_exponents", "imf_exponents");
  }

  /*! Dump the coefficients. */
  if (imf->coef != NULL) {
    restart_write_blocks((void *)imf->coef, sizeof(float), imf->n_parts, stream,
                         "imf_coef", "imf_coef");
  }
}

/**
 * @brief Restore a initial_mass_function struct from the given FILE as a stream
 * of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param imf the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void initial_mass_function_restore(struct initial_mass_function *imf,
                                   FILE *stream,
                                   const struct stellar_model *sm) {

  /* Restore the mass limits */
  if (imf->mass_limits != NULL) {
    imf->mass_limits = (float *)malloc(sizeof(float) * imf->n_parts + 1);
    restart_read_blocks((void *)imf->mass_limits, sizeof(float),
                        imf->n_parts + 1, stream, NULL, "imf_mass_limits");
  }

  /* Restore the exponents */
  if (imf->exp != NULL) {
    imf->exp = (float *)malloc(sizeof(float) * imf->n_parts);
    restart_read_blocks((void *)imf->exp, sizeof(float), imf->n_parts, stream,
                        NULL, "imf_exponents");
  }

  /* Restore the coefficients */
  if (imf->coef != NULL) {
    imf->coef = (float *)malloc(sizeof(float) * imf->n_parts);
    restart_read_blocks((void *)imf->coef, sizeof(float), imf->n_parts, stream,
                        NULL, "imf_coef");
  }
}

/**
 * @brief Clean the allocated memory.
 *
 * @param imf the #initial_mass_function.
 */
void initial_mass_function_clean(struct initial_mass_function *imf) {

  /* Free the pointers */
  free(imf->mass_limits);
  imf->mass_limits = NULL;

  free(imf->exp);
  imf->exp = NULL;

  free(imf->coef);
  imf->coef = NULL;

  free(imf->mass_fraction);
  imf->mass_fraction = NULL;
}

/** @brief Sample a power law distribution (IMF)
 *
 * @param min_mass : the minimal IMF mass.
 * @param max_mass : the maximal IMF mass.
 * @param exp : the power law slope.
 * @param x : a random number in the range [0, 1].
 */
INLINE double initial_mass_function_sample_power_law(double min_mass,
                                                     double max_mass,
                                                     double exp, double x) {

  double pmin = pow(min_mass, exp);
  double pmax = pow(max_mass, exp);
  return pow(x * (pmax - pmin) + pmin, 1. / exp);
}

/** @brief
 *
 * Note: This function does not verify if it computes the masses for the first
 * stars or not. You need to verify this before this function and pass the
 * correct values to 'minimal_discrete_mass' and 'stellar_particle_mass'.
 *
 * @param imf The #initial_mass_function.
 * @param minimal_discrete_mass
 * @param stellar_particle_mass
 * @param (return) M_continuous Mass of the continous part of the IMF.
 * @param (return) M_discrete Mass of the discrete part of the IMF.
 * @param (return) M_tot Total mass of the IMF.
 */
void initial_mass_function_compute_Mc_Md_Mtot(const struct initial_mass_function* imf,
					      const double minimal_discrete_mass,
					      const double stellar_particle_mass,
					      double* M_continuous, double* M_discrete,
					      double* M_tot) {

  /* Get the IMF mass limits (all in M_sun) */
  const float mass_min = imf->mass_min;

  double f_continuous = 0.0;

  /* f_continuous is the imf mass fraction of the continuous part (of the IMF). */
  f_continuous = initial_mass_function_get_imf_mass_fraction(imf, mass_min,
						   minimal_discrete_mass);

  /* Determine Mc and Md the masses of the continuous and discrete parts of the
     IMF, as well as Mtot the total mass of the IMF. */
  if (f_continuous > 0) {
    *M_tot = stellar_particle_mass / f_continuous;
    *M_discrete =* M_tot - stellar_particle_mass;
    *M_continuous = stellar_particle_mass;
  } else {
    *M_tot = stellar_particle_mass;
    *M_discrete = *M_tot;
    *M_continuous = 0;
  }
}
