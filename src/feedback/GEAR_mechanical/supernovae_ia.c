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

/* Include header */
#include "supernovae_ia.h"

/* Local headers */
#include "hdf5_functions.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Print the supernovae Ia model.
 *
 * @param snia The #supernovae_ia.
 */
void supernovae_ia_print(const struct supernovae_ia *snia) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  message("Mass of the white dwarf = %g", snia->mass_white_dwarf);
  message("Mass range of the progenitor = [%g, %g]", snia->mass_min_progenitor,
          snia->mass_max_progenitor);

  for (int i = 0; i < GEAR_NUMBER_TYPE_OF_COMPANION; i++) {
    message("Mass range of the companion %i = [%g, %g]", i,
            snia->companion[i].mass_min, snia->companion[i].mass_max);
  }
}

/**
 * @brief Check if the given mass is able to produce a SNIa.
 *
 * @param snia The #supernovae_ia model.
 * @param m_low The lower mass.
 * @param m_high The higher mass.
 *
 * @return If the mass is in the range of SNIa.
 */
int supernovae_ia_can_explode(const struct supernovae_ia *snia, float m_low,
                              float m_high) {

  if (m_low > snia->mass_max_progenitor) return 0;

  for (int i = 0; i < GEAR_NUMBER_TYPE_OF_COMPANION; i++) {
    if (m_low < snia->companion[i].mass_max &&
        m_high > snia->companion[i].mass_min) {
      return 1;
    }
  }

  return 0;
}

/**
 * @brief Get the yields of a supernovae Ia.
 *
 * @param snia The #supernovae_ia model.
 */
const float *supernovae_ia_get_yields(const struct supernovae_ia *snia) {
  return snia->yields;
}

/**
 * @brief Get the processed mass ejected of a supernovae Ia.
 *
 * @param snia The #supernovae_ia model.
 */
float supernovae_ia_get_ejected_mass_processed(
    const struct supernovae_ia *snia) {
  return snia->mass_white_dwarf;
}

/**
 * @brief Compute the companion integral (second integral in equation 3.46 in
 * Poirier 2004)
 *
 * @param snia The #supernovae_ia model.
 * @param m1 The lower mass limit.
 * @param m2 The upper mass limit.
 * @param companion_type The type of companion (e.g. index of snia->companion).
 *
 * @return The fraction of companion.
 */
float supernovae_ia_get_companion_fraction(const struct supernovae_ia *snia,
                                           float m1, float m2,
                                           int companion_type) {
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > m2) error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  const float tmp =
      pow(m2, snia->companion_exponent) - pow(m1, snia->companion_exponent);
  return snia->companion[companion_type].coef * tmp;
}

/**
 * @brief Compute the number of supernovae Ia per unit of mass (equation 3.46 in
 * Poirier 2004).
 *
 * @param snia The #supernovae_ia model.
 * @param m1 The lower mass limit.
 * @param m2 The upper mass limit.
 *
 * @return The number of supernovae Ia per unit of mass.
 */
float supernovae_ia_get_number_per_unit_mass(const struct supernovae_ia *snia,
                                             float m1, float m2) {

#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > m2) error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  /* Do we have white dwarf? */
  if (m1 > snia->mass_max_progenitor) {
    return 0.;
  }

  float number_companion = 0.;
  for (int i = 0; i < GEAR_NUMBER_TYPE_OF_COMPANION; i++) {
    /* Check if we are in the possible interval */
    if (m1 > snia->companion[i].mass_max || m2 < snia->companion[i].mass_min)
      continue;

    /* Get mass limits */
    const float mass_min = max(m1, snia->companion[i].mass_min);
    const float mass_max = min(m2, snia->companion[i].mass_max);

    /* Compute number of companions */
    number_companion +=
        supernovae_ia_get_companion_fraction(snia, mass_min, mass_max, i);
  }

  /* Use only the white dwarf already created */
  const float mass_min = max(m1, snia->mass_min_progenitor);

  /* Compute number of white dwarf */
  float number_white_dwarf =
      pow(snia->mass_max_progenitor, snia->progenitor_exponent);
  number_white_dwarf -= pow(mass_min, snia->progenitor_exponent);
  number_white_dwarf *= snia->progenitor_coef_exp;

  return number_companion * number_white_dwarf;
};

/**
 * @brief Read the SNIa yields from the table.
 *
 * @param snia The #supernovae_ia model.
 * @param params The #swift_params.
 * @param sm The #stellar_model.
 * @param filename The filename of the chemistry table.
 */
void supernovae_ia_read_yields(struct supernovae_ia *snia,
                               struct swift_params *params,
                               const struct stellar_model *sm,
                               const char *filename) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(filename, "Data/SNIa/Metals", &file_id, &group_id);

  /* Get the number of elements. */
  const hid_t attr = H5Aopen(group_id, "elts", H5P_DEFAULT);
  if (attr < 0) error("Error while opening attribute elts in Data/SNIa/Metals");

  size_t nval = io_get_number_element_in_attribute(attr);
  H5Aclose(attr);

  /* Read the yields */
  float *yields = (float *)malloc(sizeof(float) * nval);
  io_read_array_attribute(group_id, "data", FLOAT, yields, nval);

  /* Read the labels */
  char *labels = malloc(nval * GEAR_LABELS_SIZE);
  io_read_string_array_attribute(group_id, "elts", labels, nval,
                                 GEAR_LABELS_SIZE);

  /* Save the yields */
  /* Loop over the elements in sm */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    int found = 0;
    /* Loop over SNIa yields labels */
    for (size_t j = 0; j < nval; j++) {
      const char *s1 = labels + j * GEAR_LABELS_SIZE;
      const char *s2 = stellar_evolution_get_element_name(sm, i);
      if (strcmp(s1, s2) == 0) {
        found = 1;
        snia->yields[i] = yields[j];
        break;
      }
    }

    /* Check if found an element */
    if (!found) {
      error("Cannot find element '%s' in SNIa yields",
            stellar_evolution_get_element_name(sm, i));
    }
  }

  /* Cleanup everything */
  free(yields);
  free(labels);
  h5_close_group(file_id, group_id);
};

/**
 * @brief Initialize the companion structure in the #supernovae_ia.
 */
void supernovae_ia_init_companion(struct supernovae_ia *snia) {

  /* Get the exponent for the integral */
  const float exp = snia->companion_exponent;

  for (int i = 0; i < GEAR_NUMBER_TYPE_OF_COMPANION; i++) {

    /* Compute the integral */
    float integral = pow(snia->companion[i].mass_max, exp + 1.);
    integral -= pow(snia->companion[i].mass_min, exp + 1.);
    integral /= exp + 1.;

    /* Update the coefficient for a normalization to 1 of the IMF */
    snia->companion[i].coef /= exp * integral;
  }
}

/**
 * @brief Reads the supernovae Ia parameters from the tables.
 *
 * @param snia The #supernovae_ia model.
 * @param params The simulation parameters.
 * @param filename The filename of the chemistry table.
 */
void supernovae_ia_read_from_tables(struct supernovae_ia *snia,
                                    struct swift_params *params,
                                    const char *filename) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(filename, "Data/SNIa", &file_id, &group_id);

  /* Read the exponent of the IMF for companion */
  io_read_attribute(group_id, "a", FLOAT, &snia->companion_exponent);

  /* Read the minimal mass for a white dwarf progenitor */
  io_read_attribute(group_id, "Mpl", FLOAT, &snia->mass_min_progenitor);

  /* Read the maximal mass for a white dwarf */
  io_read_attribute(group_id, "Mpu", FLOAT, &snia->mass_max_progenitor);

  /* Read the maximal mass of a red giant companion */
  io_read_attribute(group_id, "Mdu1", FLOAT, &snia->companion[0].mass_max);

  /* Read the minimal mass of a red giant companion */
  io_read_attribute(group_id, "Mdl1", FLOAT, &snia->companion[0].mass_min);

  /* Read the coefficient of the main sequence companion */
  io_read_attribute(group_id, "bb1", FLOAT, &snia->companion[0].coef);

  /* Read the maximal mass of a main sequence companion */
  io_read_attribute(group_id, "Mdu2", FLOAT, &snia->companion[1].mass_max);

  /* Read the minimal mass of a main sequence companion */
  io_read_attribute(group_id, "Mdl2", FLOAT, &snia->companion[1].mass_min);

  /* Read the coefficient of the main sequence companion */
  io_read_attribute(group_id, "bb2", FLOAT, &snia->companion[1].coef);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);

  /* Read the white dwarf mass */

  /* Open IMF group */
  h5_open_group(filename, "Data", &file_id, &group_id);

  /* Read the white dwarf mass */
  io_read_attribute(group_id, "MeanWDMass", FLOAT, &snia->mass_white_dwarf);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Initialize the #supernovae_ia structure.
 *
 * @param snia The #supernovae_ia model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 */
void supernovae_ia_init(struct supernovae_ia *snia,
                        const struct phys_const *phys_const,
                        const struct unit_system *us,
                        struct swift_params *params,
                        const struct stellar_model *sm) {

  /* Read the parameters from the tables */
  supernovae_ia_read_from_tables(snia, params, sm->yields_table);

  /* Read the yields */
  supernovae_ia_read_yields(snia, params, sm, sm->yields_table);

  /* Get the IMF parameters */
  snia->progenitor_exponent = initial_mass_function_get_exponent(
      &sm->imf, snia->mass_min_progenitor, snia->mass_max_progenitor);
  snia->progenitor_coef_exp = initial_mass_function_get_coefficient(
      &sm->imf, snia->mass_min_progenitor, snia->mass_max_progenitor);
  snia->progenitor_coef_exp /= snia->progenitor_exponent;

  /* Compute the normalization coefficients of the companion IMF */
  supernovae_ia_init_companion(snia);

  /* Supernovae energy */
  double e_feedback =
      parser_get_param_double(params, "GEARFeedback:supernovae_energy_erg");
  e_feedback /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  snia->energy_per_supernovae = e_feedback;
}

/**
 * @brief Write a supernovae_ia struct to the given FILE as a stream of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param snia the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void supernovae_ia_dump(const struct supernovae_ia *snia, FILE *stream,
                        const struct stellar_model *sm) {

  /* Nothing to do here */
}

/**
 * @brief Restore a supernovae_ia struct from the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param snia the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void supernovae_ia_restore(struct supernovae_ia *snia, FILE *stream,
                           const struct stellar_model *sm) {

  /* Nothing to do here */
}

/**
 * @brief Clean the allocated memory.
 *
 * @param snia the #supernovae_ia.
 */
void supernovae_ia_clean(struct supernovae_ia *snia) {}
