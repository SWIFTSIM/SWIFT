/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
/**
 * @file src/feedback/GEAR/radiation_table_io.c
 * @brief HDF5 reading and interpolation-table building for GEAR radiation
 * feedback.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "engine.h"
#include "error.h"
#include "hdf5_functions.h"
#include "interpolation.h"
#include "minmax.h"
#include "radiation.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"
#include "units.h"

#include <float.h>
#include <string.h>

/**
 * @brief Read a scalar HDF5 string attribute (fixed- or variable-length)
 * into a NUL-terminated buffer.
 *
 * io_read_attribute() cannot read the variable-length UTF-8 strings that
 * pychem writes.
 *
 * @param group_id Open HDF5 group id.
 * @param name Attribute name.
 * @param out Output buffer, NUL-terminated on return.
 * @param out_size Size of out, including the terminating NUL.
 * @param truncate If 0, an attribute too long for out_size is fatal; if
 * non-zero, it is silently truncated to fit.
 */
static void radiation_read_string_attribute_impl(hid_t group_id,
                                                 const char *name, char *out,
                                                 size_t out_size,
                                                 const int truncate) {

  const hid_t h_attr = H5Aopen(group_id, name, H5P_DEFAULT);
  if (h_attr < 0) error("Error while opening attribute '%s'", name);

  const hid_t h_type = H5Aget_type(h_attr);
  if (h_type < 0) error("Error while getting the type of attribute '%s'", name);

  if (H5Tis_variable_str(h_type) > 0) {
    /* Read into the attribute's native type: a fresh H5T_C_S1 type differs in
       character set (ASCII vs UTF-8) and has no registered conversion. */
    char *tmp = NULL;
    if (H5Aread(h_attr, h_type, &tmp) < 0)
      error("Error while reading string attribute '%s'", name);

    const size_t len = strlen(tmp);
    if (len >= out_size && !truncate) {
      error(
          "String attribute '%s' (%zu bytes) does not fit in the %zu-byte "
          "buffer.",
          name, len, out_size);
    }

    strncpy(out, tmp, out_size - 1);
    out[out_size - 1] = '\0';

    const hid_t h_space = H5Aget_space(h_attr);
    H5Dvlen_reclaim(h_type, h_space, H5P_DEFAULT, &tmp);
    H5Sclose(h_space);
  } else {
    const size_t fixed_size = H5Tget_size(h_type);
    if (fixed_size >= out_size && !truncate)
      error(
          "String attribute '%s' (%zu bytes) does not fit in the %zu-byte "
          "buffer.",
          name, fixed_size, out_size);

    char *tmp = (char *)calloc(fixed_size + 1, sizeof(char));
    if (tmp == NULL) error("Failed to allocate string attribute buffer.");
    if (H5Aread(h_attr, h_type, tmp) < 0)
      error("Error while reading string attribute '%s'", name);
    const size_t kept = min(fixed_size, out_size - 1);
    memcpy(out, tmp, kept);
    out[kept] = '\0';
    free(tmp);
  }

  H5Tclose(h_type);
  H5Aclose(h_attr);
}

/**
 * @brief Read a string attribute, failing if it does not fit.
 *
 * The caller parses the result, so truncation is not acceptable.
 *
 * @param group_id Open HDF5 group id.
 * @param name The attribute's name.
 * @param out (output) The buffer to fill.
 * @param out_size The size of @p out, in bytes.
 */
static void radiation_read_string_attribute(hid_t group_id, const char *name,
                                            char *out, size_t out_size) {
  radiation_read_string_attribute_impl(group_id, name, out, out_size,
                                       /*truncate=*/0);
}

/**
 * @brief Read a string attribute, truncating it if it does not fit.
 *
 * For values that are only reported, never parsed.
 *
 * @param group_id Open HDF5 group id.
 * @param name The attribute's name.
 * @param out (output) The buffer to fill.
 * @param out_size The size of @p out, in bytes.
 */
static void radiation_read_string_attribute_truncating(hid_t group_id,
                                                       const char *name,
                                                       char *out,
                                                       size_t out_size) {
  radiation_read_string_attribute_impl(group_id, name, out, out_size,
                                       /*truncate=*/1);
}

const char *const radiation_table_source_keys[RADIATION_TABLE_SOURCE_COUNT] = {
    "qh_source", "lwpe_source", "stellar_evolution_source", "source"};

/**
 * @brief Read the table's own provenance into the #radiation model.
 *
 * Reported by #radiation_print. Every attribute is optional: a table with
 * none leaves every #radiation.table_source entry empty.
 *
 * @param rad (output) The #radiation model to fill.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The #radiation_grid_metadata already read from @p group_id.
 */
static void radiation_read_table_identity(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid) {

  for (int i = 0; i < RADIATION_TABLE_SOURCE_COUNT; i++) {
    rad->table_source[i][0] = '\0';

    if (H5Aexists(group_id, radiation_table_source_keys[i]) <= 0) continue;

    radiation_read_string_attribute_truncating(
        group_id, radiation_table_source_keys[i], rad->table_source[i],
        RADIATION_TABLE_SOURCE_SIZE);
  }

  rad->table_n_mass = grid->n_mass;
  rad->table_mass_min = (float)exp10((double)grid->log_mass_min);
  rad->table_mass_max =
      (float)exp10((double)grid->log_mass_min +
                   (grid->n_mass - 1) * (double)grid->mass_step);

  if (grid->is_2d) {
    rad->table_n_metallicity = grid->n_metallicity;
    rad->table_metallicity_min = grid->metallicity[0];
    rad->table_metallicity_max = grid->metallicity[grid->n_metallicity - 1];
  } else {
    rad->table_n_metallicity = 0;
    rad->table_metallicity_min = 0.f;
    rad->table_metallicity_max = 0.f;
  }
}

/**
 * @brief Check that a Data/Radiation dataset's "units" attribute matches
 * the unit the reader assumes.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to check.
 * @param expected_units The expected unit string (e.g. "erg/s"), compared
 * verbatim.
 */
static void radiation_check_dataset_units(hid_t group_id,
                                          const char *dataset_name,
                                          const char *expected_units) {

  const hid_t h_dataset = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
  if (h_dataset < 0)
    error("Error while opening dataset '%s' to check its units.", dataset_name);

  char actual_units[64];
  radiation_read_string_attribute(h_dataset, "units", actual_units,
                                  sizeof(actual_units));

  H5Dclose(h_dataset);

  if (strcmp(actual_units, expected_units) != 0) {
    error(
        "Data/Radiation/%s declares units='%s', but SWIFT's reader assumes "
        "'%s'. This table's unit convention no longer matches what this "
        "code converts. Aborting rather than silently misinterpreting "
        "the physics.",
        dataset_name, actual_units, expected_units);
  }
}

/**
 * @brief Turn one field's edge_policy_<field>_below/above attribute pair
 * into an #interpolate_boundary_condition on the mass axis.
 *
 * "linear" is not supported, and "constant below / zero above" has no
 * #interpolate_boundary_condition. Both stop the run.
 *
 * @param below The field's edge_policy_<field>_below value.
 * @param above The field's edge_policy_<field>_above value.
 * @param field_name Name of the field these came from, for the error
 * message only.
 * @return The matching #interpolate_boundary_condition.
 */
static enum interpolate_boundary_condition radiation_parse_edge_policy(
    const char *below, const char *above, const char *field_name) {

  int zero_below = 0;
  if (strcmp(below, "zero") == 0) {
    zero_below = 1;
  } else if (strcmp(below, "constant") != 0) {
    error(
        "Data/Radiation field '%s': edge_policy_%s_below='%s' has no SWIFT "
        "#interpolate_boundary_condition equivalent (only \"zero\" and "
        "\"constant\" are supported).",
        field_name, field_name, below);
  }

  int zero_above = 0;
  if (strcmp(above, "zero") == 0) {
    zero_above = 1;
  } else if (strcmp(above, "constant") != 0) {
    error(
        "Data/Radiation field '%s': edge_policy_%s_above='%s' has no SWIFT "
        "#interpolate_boundary_condition equivalent (only \"zero\" and "
        "\"constant\" are supported).",
        field_name, field_name, above);
  }

  if (zero_below && zero_above) return boundary_condition_zero;
  if (zero_below && !zero_above) return boundary_condition_zero_const;
  if (!zero_below && !zero_above) return boundary_condition_const;

  error(
      "Data/Radiation field '%s': edge policy 'constant below / zero "
      "above' has no SWIFT #interpolate_boundary_condition equivalent.",
      field_name);
  return boundary_condition_error;
}

/**
 * @brief Read the grid metadata of the Data/Radiation group.
 *
 * Reads "dimensionality" ("M" or "M,Z"), the "m0"/"dm"/"nm" mass grid and,
 * for "M,Z", "nz", "Metallicity" and the mass-axis edge_policy_* attributes,
 * see radiation_parse_edge_policy(). The edge_policy_* of L_PE, L_LW and Teff
 * are read only if their dataset exists.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid (output) The #radiation_grid_metadata to fill in.
 */
void radiation_read_grid_metadata(hid_t group_id,
                                  struct radiation_grid_metadata *grid) {

  radiation_read_string_attribute(group_id, "dimensionality",
                                  grid->dimensionality,
                                  sizeof(grid->dimensionality));

  io_read_attribute(group_id, "m0", FLOAT, &grid->log_mass_min);
  io_read_attribute(group_id, "dm", FLOAT, &grid->mass_step);
  io_read_attribute(group_id, "nm", INT, &grid->n_mass);

  if (!(grid->mass_step > 0.f))
    error(
        "Data/Radiation's 'dm' attribute is %.4g; it must be a strictly "
        "positive mass-grid step.",
        (double)grid->mass_step);

  if (grid->n_mass < 2)
    error(
        "Data/Radiation's 'nm' attribute is %d; at least 2 mass points are "
        "needed to interpolate.",
        grid->n_mass);

  grid->n_metallicity = 0;
  grid->metallicity = NULL;

  if (strcmp(grid->dimensionality, "M,Z") == 0) {
    grid->is_2d = 1;

    io_read_attribute(group_id, "nz", INT, &grid->n_metallicity);

    if (grid->n_metallicity < 2)
      error(
          "Data/Radiation's 'nz' attribute is %d; at least 2 metallicity "
          "points are needed to interpolate the log10(Z) axis.",
          grid->n_metallicity);

    grid->metallicity = (float *)malloc(sizeof(float) * grid->n_metallicity);
    if (grid->metallicity == NULL)
      error("Failed to allocate the RAD metallicity grid.");

    io_read_array_dataset(group_id, "Metallicity", FLOAT, grid->metallicity,
                          grid->n_metallicity);

    for (int i = 0; i < grid->n_metallicity; i++) {
      if (grid->metallicity[i] <= 0.f)
        error(
            "Data/Radiation's 'Metallicity' dataset entry %d is %.4g "
            "(<= 0): the metallicity axis is interpolated in log10(Z), "
            "which requires every entry to be strictly positive.",
            i, (double)grid->metallicity[i]);
      if (i > 0 && !(grid->metallicity[i] > grid->metallicity[i - 1]))
        error(
            "Data/Radiation's 'Metallicity' dataset is not strictly "
            "increasing at entry %d (%.4g after %.4g): the metallicity axis "
            "is bracketed node by node, which requires a sorted grid.",
            i, (double)grid->metallicity[i], (double)grid->metallicity[i - 1]);
    }

    char luminosity_below[16], luminosity_above[16];
    char q_h_below[16], q_h_above[16];
    char mean_excess_energy_below[16], mean_excess_energy_above[16];
    radiation_read_string_attribute(group_id, "edge_policy_luminosity_below",
                                    luminosity_below, sizeof(luminosity_below));
    radiation_read_string_attribute(group_id, "edge_policy_luminosity_above",
                                    luminosity_above, sizeof(luminosity_above));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_below",
                                    q_h_below, sizeof(q_h_below));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_above",
                                    q_h_above, sizeof(q_h_above));
    radiation_read_string_attribute(
        group_id, "edge_policy_mean_excess_energy_below",
        mean_excess_energy_below, sizeof(mean_excess_energy_below));
    radiation_read_string_attribute(
        group_id, "edge_policy_mean_excess_energy_above",
        mean_excess_energy_above, sizeof(mean_excess_energy_above));

    grid->edge_policy_luminosity = radiation_parse_edge_policy(
        luminosity_below, luminosity_above, "luminosity");
    grid->edge_policy_q_h =
        radiation_parse_edge_policy(q_h_below, q_h_above, "q_h");
    grid->edge_policy_dot_e_excess = radiation_parse_edge_policy(
        mean_excess_energy_below, mean_excess_energy_above,
        "mean_excess_energy");

    /* Teff is optional: older tables carry neither dataset nor attributes. */
    grid->edge_policy_teff = boundary_condition_error;
    if (H5Lexists(group_id, "Teff", H5P_DEFAULT) > 0) {
      char teff_below[16], teff_above[16];
      radiation_read_string_attribute(group_id, "edge_policy_teff_below",
                                      teff_below, sizeof(teff_below));
      radiation_read_string_attribute(group_id, "edge_policy_teff_above",
                                      teff_above, sizeof(teff_above));
      grid->edge_policy_teff =
          radiation_parse_edge_policy(teff_below, teff_above, "teff");
    }

    /* L_PE/L_LW are optional, see #with_ISRF: guard on the dataset so an old
       table still loads. */
    grid->edge_policy_l_pe = boundary_condition_error;
    if (H5Lexists(group_id, "L_PE", H5P_DEFAULT) > 0) {
      char l_pe_below[16], l_pe_above[16];
      radiation_read_string_attribute(group_id, "edge_policy_l_pe_below",
                                      l_pe_below, sizeof(l_pe_below));
      radiation_read_string_attribute(group_id, "edge_policy_l_pe_above",
                                      l_pe_above, sizeof(l_pe_above));
      grid->edge_policy_l_pe =
          radiation_parse_edge_policy(l_pe_below, l_pe_above, "l_pe");
    }

    grid->edge_policy_l_lw = boundary_condition_error;
    if (H5Lexists(group_id, "L_LW", H5P_DEFAULT) > 0) {
      char l_lw_below[16], l_lw_above[16];
      radiation_read_string_attribute(group_id, "edge_policy_l_lw_below",
                                      l_lw_below, sizeof(l_lw_below));
      radiation_read_string_attribute(group_id, "edge_policy_l_lw_above",
                                      l_lw_above, sizeof(l_lw_above));
      grid->edge_policy_l_lw =
          radiation_parse_edge_policy(l_lw_below, l_lw_above, "l_lw");
    }
  } else if (strcmp(grid->dimensionality, "M") == 0) {
    grid->is_2d = 0;
    grid->edge_policy_luminosity = boundary_condition_error;
    grid->edge_policy_q_h = boundary_condition_error;
    grid->edge_policy_dot_e_excess = boundary_condition_error;
    grid->edge_policy_teff = boundary_condition_error;
    grid->edge_policy_l_pe = boundary_condition_error;
    grid->edge_policy_l_lw = boundary_condition_error;
  } else {
    error(
        "Data/Radiation has an unrecognised 'dimensionality' attribute "
        "'%s' (expected 'M' or 'M,Z').",
        grid->dimensionality);
  }
}

/**
 * @brief Check a table's IMF attributes ("imf_a_s", "imf_m_s",
 * "mass_min_msun", "mass_max_msun") against the run's IMF.
 *
 * A table without "Integrated_Q_H" is skipped. "imf_m_s" holds only the
 * interior breakpoints, compared with #initial_mass_function.mass_limits[1 ..
 * n_parts - 1]; mass_min_msun and mass_max_msun stand for the outer two.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param sm The #stellar_model (its imf must already be initialised).
 */
static void radiation_check_imf_consistency(hid_t group_id,
                                            const struct stellar_model *sm) {

  const htri_t exists = H5Lexists(group_id, "Integrated_Q_H", H5P_DEFAULT);
  if (exists <= 0) return;

  const struct initial_mass_function *imf = &sm->imf;

  double mass_min_msun, mass_max_msun;
  io_read_attribute(group_id, "mass_min_msun", DOUBLE, &mass_min_msun);
  io_read_attribute(group_id, "mass_max_msun", DOUBLE, &mass_max_msun);

  /* Check the segment count first, for a clearer error. */
  const hid_t attr_a_s = H5Aopen(group_id, "imf_a_s", H5P_DEFAULT);
  if (attr_a_s < 0) error("Error while opening attribute 'imf_a_s'");
  const hsize_t n_a_s = io_get_number_element_in_attribute(attr_a_s);
  H5Aclose(attr_a_s);
  if (n_a_s != (hsize_t)imf->n_parts)
    error(
        "Data/Radiation's 'imf_a_s' attribute has %llu segment(s), but "
        "sm->imf (read from this same file's Data/IMF group) has %d: "
        "Data/Radiation's Integrated_* datasets were precomputed against a "
        "different IMF than this file's own Data/IMF. Regenerate this "
        "table with pychem so Data/Radiation and Data/IMF agree, or point "
        "GEARFeedback:yields_table (or yields_table_first_stars) at a "
        "table whose Data/IMF already matches its own Integrated_* "
        "datasets.",
        (unsigned long long)n_a_s, imf->n_parts);

  double *imf_a_s = (double *)malloc(sizeof(double) * imf->n_parts);
  if (imf_a_s == NULL) error("Failed to allocate the 'imf_a_s' buffer.");
  io_read_array_attribute(group_id, "imf_a_s", DOUBLE, imf_a_s,
                          (hsize_t)imf->n_parts);

  const int n_interior = imf->n_parts - 1;
  double *imf_m_s = NULL;
  if (n_interior > 0) {
    const hid_t attr_m_s = H5Aopen(group_id, "imf_m_s", H5P_DEFAULT);
    if (attr_m_s < 0) error("Error while opening attribute 'imf_m_s'");
    const hsize_t n_m_s = io_get_number_element_in_attribute(attr_m_s);
    H5Aclose(attr_m_s);
    if (n_m_s != (hsize_t)n_interior)
      error(
          "Data/Radiation's 'imf_m_s' attribute has %llu segment(s), but "
          "sm->imf (read from this same file's Data/IMF group) has %d "
          "interior breakpoint(s): Data/Radiation's Integrated_* datasets "
          "were precomputed against a different IMF than this file's own "
          "Data/IMF. Regenerate this table with pychem so Data/Radiation "
          "and Data/IMF agree, or point GEARFeedback:yields_table (or "
          "yields_table_first_stars) at a table whose Data/IMF already "
          "matches its own Integrated_* datasets.",
          (unsigned long long)n_m_s, n_interior);

    imf_m_s = (double *)malloc(sizeof(double) * n_interior);
    if (imf_m_s == NULL) error("Failed to allocate the 'imf_m_s' buffer.");
    io_read_array_attribute(group_id, "imf_m_s", DOUBLE, imf_m_s,
                            (hsize_t)n_interior);
  }

  for (int k = 0; k <= imf->n_parts; k++) {
    const double expected = (k == 0)              ? mass_min_msun
                            : (k == imf->n_parts) ? mass_max_msun
                                                  : imf_m_s[k - 1];
    const float a = imf->mass_limits[k];
    const float b = (float)expected;
    if (fabsf(a - b) > 1e-4f * fmaxf(1.f, fabsf(a)))
      error(
          "Data/Radiation's IMF attrs disagree with sm->imf.mass_limits[%d] "
          "(table=%.8g Msun, SWIFT=%.8g Msun, read from this same file's "
          "Data/IMF group): Data/Radiation's Integrated_* datasets were "
          "precomputed against a different IMF than this file's own "
          "Data/IMF. Regenerate this table with pychem so Data/Radiation "
          "and Data/IMF agree, or point GEARFeedback:yields_table (or "
          "yields_table_first_stars) at a table whose Data/IMF already "
          "matches its own Integrated_* datasets.",
          k, b, a);
  }

  for (int k = 0; k < imf->n_parts; k++) {
    const float a = imf->exp[k];
    const float b = (float)imf_a_s[k];
    if (fabsf(a - b) > 1e-4f * fmaxf(1.f, fabsf(a)))
      error(
          "Data/Radiation's 'imf_a_s' attribute disagrees with "
          "sm->imf.exp[%d] (table=%.8g, SWIFT=%.8g, read from this same "
          "file's Data/IMF group): Data/Radiation's Integrated_* datasets "
          "were precomputed against a different IMF than this file's own "
          "Data/IMF. Regenerate this table with pychem so Data/Radiation "
          "and Data/IMF agree, or point GEARFeedback:yields_table (or "
          "yields_table_first_stars) at a table whose Data/IMF already "
          "matches its own Integrated_* datasets.",
          k, b, a);
  }

  free(imf_a_s);
  free(imf_m_s);
}

/**
 * @brief Read one CGS-valued dataset of Data/Radiation, convert it to
 * internal units (and #RADIATION_DOT_N_ION_TABLE_SCALING) as float, and
 * optionally compute its log10.
 *
 * Aborts on float overflow. With debug checks, warns when a nonzero CGS value
 * collapses to exactly zero.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to read.
 * @param count Number of elements ("nm", or "nm" * "nz" for a 2D table).
 * @param conversion_factor CGS values are divided by this to reach internal
 * units.
 * @param extra_scaling Additional divisor (#RADIATION_DOT_N_ION_TABLE_SCALING
 * for Q_H/DotEExcess, 1 for Luminosity).
 * @param expected_units Asserted against the dataset's "units" attribute.
 * @param log_data_internal (output, optional) If not NULL, array of length
 * @p count filled with log10 of the internal-unit value, floored at
 * #RADIATION_LOG_FLOOR_CGS on the CGS side.
 * @return Newly malloc'd float array of length @p count, in linear (not
 * logged) internal units. Caller must free().
 */
static float *radiation_read_cgs_array(hid_t group_id, const char *dataset_name,
                                       hsize_t count, double conversion_factor,
                                       double extra_scaling,
                                       const char *expected_units,
                                       float *log_data_internal) {

  radiation_check_dataset_units(group_id, dataset_name, expected_units);

  double *data_cgs = (double *)malloc(sizeof(double) * count);
  if (data_cgs == NULL)
    error("Failed to allocate the RAD yields for %s.", dataset_name);

  io_read_array_dataset(group_id, dataset_name, DOUBLE, data_cgs, count);

  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for %s.", dataset_name);

  /* log10(internal value) = log10(cgs value) - log_conversion */
  const double log_conversion = log10(conversion_factor) + log10(extra_scaling);

  for (hsize_t j = 0; j < count; j++) {
    const double value_internal =
        data_cgs[j] / conversion_factor / extra_scaling;

    if (fabs(value_internal) > (double)FLT_MAX) {
      error(
          "Radiation table '%s' entry %llu (%e cgs) converts to %e in "
          "internal units, exceeding FLT_MAX. This is a units/scaling bug; "
          "aborting rather than silently corrupting the physics.",
          dataset_name, (unsigned long long)j, data_cgs[j], value_internal);
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (data_cgs[j] != 0. && (float)value_internal == 0.0f) {
      message(
          "WARNING: radiation table '%s' entry %llu (%e cgs, nonzero) "
          "collapsed to exactly 0 in internal units after conversion. "
          "Check RADIATION_DOT_N_ION_TABLE_SCALING and the unit system.",
          dataset_name, (unsigned long long)j, data_cgs[j]);
    }
#endif

    data[j] = (float)value_internal;

    if (log_data_internal != NULL) {
      const double floored_cgs = max(data_cgs[j], RADIATION_LOG_FLOOR_CGS);
      const double log_value_internal = log10(floored_cgs) - log_conversion;
      log_data_internal[j] = (float)log_value_internal;

#ifdef SWIFT_DEBUG_CHECKS
      /* Check the log10/exp10 round trip against the linear value, away from
         the floor. */
      if (data_cgs[j] > RADIATION_LOG_FLOOR_CGS * 1e10) {
        const double round_trip = exp10(log_value_internal);
        const double rel_diff =
            fabs(round_trip - value_internal) / fabs(value_internal);
        if (rel_diff > 1e-4) {
          error(
              "Radiation table '%s' entry %llu: log10/exp10 round-trip "
              "mismatch (internal=%e, round-trip=%e, rel_diff=%e). This "
              "indicates a bug in the log-log conversion, not the data.",
              dataset_name, (unsigned long long)j, value_internal, round_trip,
              rel_diff);
        }
      }
#endif
    }
  }

  free(data_cgs);
  return data;
}

/**
 * @brief Build the raw and IMF-integrated interpolation tables of one
 * Data/Radiation quantity, for a 1D or 2D table.
 *
 * The raw table holds log10(value) so that the linear interpolators act as
 * log-log ones: every raw getter must exponentiate the result. The
 * integrated table stays in linear space. It is read from the table's
 * "Integrated_<dataset_name>" dataset (number-weighted, cumulative from
 * Mmin) and shares the output-grid bounds of @p raw_2d, so a two-point
 * subtraction query never sees a Z-axis mismatch. The 2D metallicity axis
 * uses the "Metallicity" nodes and is linear in log10(Z).
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to read.
 * @param grid The group's own grid metadata (see
 * radiation_read_grid_metadata()).
 * @param sm The #stellar_model (for the output mass-grid bounds).
 * @param interpolation_size_mass Number of points in the mass
 * interpolation output grid.
 * @param conversion_factor See radiation_read_cgs_array().
 * @param extra_scaling See radiation_read_cgs_array().
 * @param expected_units See radiation_read_cgs_array().
 * @param raw_1d (output) Raw 1D interpolation table (1D tables only),
 * holding log10(value in internal units), pychem-floored.
 * @param integrated_1d (output, optional) IMF-integrated 1D table (linear
 * values). NULL for a dataset with no IMF-integrated concept.
 * @param raw_2d (output) Raw 2D interpolation table (2D tables only),
 * holding log10(value in internal units), pychem-floored.
 * @param integrated_2d (output, optional) IMF-integrated 2D table (linear
 * values). NULL for a dataset with no IMF-integrated concept
 * (MainSequenceLifetime).
 * @param boundary_condition_mass Mass-axis boundary condition of @p raw_2d
 * (ignored for a 1D table, which clamps). The integrated tables and the
 * metallicity axis always clamp.
 */
static void radiation_build_tables(
    hid_t group_id, const char *dataset_name,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    int interpolation_size_mass, double conversion_factor, double extra_scaling,
    const char *expected_units, struct interpolation_1d *raw_1d,
    struct interpolation_1d *integrated_1d, struct interpolation_2d *raw_2d,
    struct interpolation_2d *integrated_2d,
    enum interpolate_boundary_condition boundary_condition_mass) {

  const float log_mass_min_out = log10f(sm->imf.mass_min);
  const float log_mass_max_out = log10f(sm->imf.mass_max);

  if (grid->is_2d) {

    const hsize_t count = (hsize_t)grid->n_mass * (hsize_t)grid->n_metallicity;
    float *log_data = (float *)malloc(sizeof(float) * count);
    if (log_data == NULL)
      error("Failed to allocate the RAD 2D log-value yields for %s.",
            dataset_name);
    float *data = radiation_read_cgs_array(group_id, dataset_name, count,
                                           conversion_factor, extra_scaling,
                                           expected_units, log_data);

    float *log_z_nodes = (float *)malloc(sizeof(float) * grid->n_metallicity);
    if (log_z_nodes == NULL)
      error("Failed to allocate the RAD 2D log10(Z) axis for %s.",
            dataset_name);
    for (int i = 0; i < grid->n_metallicity; i++)
      log_z_nodes[i] = log10f(grid->metallicity[i]);

    /* interpolate_2d_init() takes a double source array. */
    double *log_data_double = (double *)malloc(sizeof(double) * count);
    if (log_data_double == NULL)
      error("Failed to allocate the RAD 2D log-value yields for %s.",
            dataset_name);
    for (hsize_t i = 0; i < count; i++)
      log_data_double[i] = (double)log_data[i];

    interpolate_2d_init(raw_2d, log_z_nodes, grid->n_metallicity, log_z_nodes,
                        grid->n_metallicity, log_mass_min_out, log_mass_max_out,
                        interpolation_size_mass, grid->log_mass_min,
                        grid->mass_step, grid->n_mass, log_data_double,
                        boundary_condition_const, boundary_condition_mass);

    free(log_data_double);
    free(log_data);
    free(data);

    if (integrated_2d == NULL) {
      free(log_z_nodes);
      return;
    }

    /* integrated_2d is read from the precomputed integral, see above. */
    char integrated_dataset_name[64];
    int written =
        snprintf(integrated_dataset_name, sizeof(integrated_dataset_name),
                 "Integrated_%s", dataset_name);
    if (written < 0 || (size_t)written >= sizeof(integrated_dataset_name))
      error("Dataset name 'Integrated_%s' does not fit in the buffer.",
            dataset_name);

    const htri_t integrated_exists =
        H5Lexists(group_id, integrated_dataset_name, H5P_DEFAULT);
    if (integrated_exists <= 0) {
      error(
          "This Data/Radiation group has no '%s' dataset. This table was "
          "generated before pychem added the precomputed IMF-integrated "
          "datasets and needs regenerating: run pychem's "
          "pychem_generate_hdf5_parameters on this table's own chimieparam "
          "file, then point GEARFeedback:yields_table (or "
          "yields_table_first_stars, for the PopIII model) at the "
          "regenerated file.",
          integrated_dataset_name);
    }

    char integrated_expected_units[32];
    written =
        snprintf(integrated_expected_units, sizeof(integrated_expected_units),
                 "%s/Msun", expected_units);
    if (written < 0 || (size_t)written >= sizeof(integrated_expected_units))
      error("Units string '%s/Msun' does not fit in the buffer.",
            expected_units);

    /* NULL: the integral stays in linear space. */
    float *integrated_data = radiation_read_cgs_array(
        group_id, integrated_dataset_name, count, conversion_factor,
        extra_scaling, integrated_expected_units, NULL);

    double *integrated_data_double = (double *)malloc(sizeof(double) * count);
    if (integrated_data_double == NULL)
      error("Failed to allocate the RAD 2D integrated yields for %s.",
            dataset_name);
    for (hsize_t i = 0; i < count; i++)
      integrated_data_double[i] = (double)integrated_data[i];

    /* Both axes clamp. */
    interpolate_2d_init(integrated_2d, log_z_nodes, grid->n_metallicity,
                        log_z_nodes, grid->n_metallicity, log_mass_min_out,
                        log_mass_max_out, interpolation_size_mass,
                        grid->log_mass_min, grid->mass_step, grid->n_mass,
                        integrated_data_double, boundary_condition_const,
                        boundary_condition_const);

    free(integrated_data_double);
    free(integrated_data);
    free(log_z_nodes);
    return;
  }

  float *log_data = (float *)malloc(sizeof(float) * grid->n_mass);
  if (log_data == NULL)
    error("Failed to allocate the RAD log-value yields for %s.", dataset_name);
  float *data = radiation_read_cgs_array(
      group_id, dataset_name, (hsize_t)grid->n_mass, conversion_factor,
      extra_scaling, expected_units, log_data);

  interpolate_1d_init(raw_1d, log_mass_min_out, log_mass_max_out,
                      interpolation_size_mass, grid->log_mass_min,
                      grid->mass_step, grid->n_mass, log_data,
                      boundary_condition_const);
  free(data);
  free(log_data);

  if (integrated_1d == NULL) return;

  /* integrated_1d is read from the precomputed integral, see above. */
  char integrated_dataset_name[64];
  int written =
      snprintf(integrated_dataset_name, sizeof(integrated_dataset_name),
               "Integrated_%s", dataset_name);
  if (written < 0 || (size_t)written >= sizeof(integrated_dataset_name))
    error("Dataset name 'Integrated_%s' does not fit in the buffer.",
          dataset_name);

  const htri_t integrated_exists =
      H5Lexists(group_id, integrated_dataset_name, H5P_DEFAULT);
  if (integrated_exists <= 0) {
    error(
        "This Data/Radiation group has no '%s' dataset. This table was "
        "generated before pychem added the precomputed IMF-integrated "
        "datasets and needs regenerating: run pychem's "
        "pychem_generate_hdf5_parameters on this table's own chimieparam "
        "file, then point GEARFeedback:yields_table (or "
        "yields_table_first_stars, for the PopIII model) at the "
        "regenerated file.",
        integrated_dataset_name);
  }

  /* "per Msun" is a fixed physical mass, so the conversion of the raw
     sibling applies unchanged. */
  char integrated_expected_units[32];
  written =
      snprintf(integrated_expected_units, sizeof(integrated_expected_units),
               "%s/Msun", expected_units);
  if (written < 0 || (size_t)written >= sizeof(integrated_expected_units))
    error("Units string '%s/Msun' does not fit in the buffer.", expected_units);

  /* NULL: the integral stays in linear space. */
  float *integrated_data = radiation_read_cgs_array(
      group_id, integrated_dataset_name, (hsize_t)grid->n_mass,
      conversion_factor, extra_scaling, integrated_expected_units, NULL);

  interpolate_1d_init(integrated_1d, log_mass_min_out, log_mass_max_out,
                      interpolation_size_mass, grid->log_mass_min,
                      grid->mass_step, grid->n_mass, integrated_data,
                      boundary_condition_const);

  free(integrated_data);
}

/**
 * @brief Read an array of luminosities data from the table.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_luminosities_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(
      group_id, "Luminosity", grid, sm, rad->interpolation_size,
      units_cgs_conversion_factor(us, UNIT_CONV_POWER), 1., "erg/s",
      &rad->raw.luminosities, &rad->integrated.luminosities,
      &rad->raw.luminosities_2d, &rad->integrated.luminosities_2d,
      grid->edge_policy_luminosity);
}

/**
 * @brief Read an array of ionizing emission rates data from the table.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_ionization_rate_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(
      group_id, "Q_H", grid, sm, rad->interpolation_size,
      units_cgs_conversion_factor(us, UNIT_CONV_PHOTONS_PER_TIME),
      RADIATION_DOT_N_ION_TABLE_SCALING, "1/s", &rad->raw.dot_N_ion,
      &rad->integrated.dot_N_ion, &rad->raw.dot_N_ion_2d,
      &rad->integrated.dot_N_ion_2d, grid->edge_policy_q_h);
}

/**
 * @brief Read an array of excess-photon-energy emission rate data from the
 * table: DotEExcess(m) = Q_H(m) * MeanExcessPhotonEnergyHI(m).
 *
 * Converted with the rate-only #UNIT_CONV_PHOTONS_PER_TIME factor, as for
 * Q_H, although the table stores erg/s. The scaling then cancels in
 * #radiation_get_mean_excess_photon_energy_HI_from_integral, which returns
 * cgs erg. The units check still asserts "erg/s".
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_mean_excess_photon_energy_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(
      group_id, "DotEExcess", grid, sm, rad->interpolation_size,
      units_cgs_conversion_factor(us, UNIT_CONV_PHOTONS_PER_TIME),
      RADIATION_DOT_N_ION_TABLE_SCALING, "erg/s", &rad->raw.dot_E_excess,
      &rad->integrated.dot_E_excess, &rad->raw.dot_E_excess_2d,
      &rad->integrated.dot_E_excess_2d, grid->edge_policy_dot_e_excess);
}

/**
 * @brief Read the Teff (photospheric effective temperature) array from the
 * table.
 *
 * Only called when the group carries a "Teff" dataset (#radiation.has_teff).
 * Raw-only: there is no "Integrated_Teff" dataset.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_teff_array(struct radiation *rad, hid_t group_id,
                               const struct radiation_grid_metadata *grid,
                               const struct stellar_model *sm,
                               const struct unit_system *us) {

  radiation_build_tables(group_id, "Teff", grid, sm, rad->interpolation_size,
                         units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE),
                         1., "K", &rad->raw.teff, NULL, &rad->raw.teff_2d, NULL,
                         grid->edge_policy_teff);
}

/**
 * @brief Check the provenance of a band-edge spectral-rate dataset.
 *
 * Its "lower_edge_energy" must match the band edge SWIFT assumes, and
 * "native_grid_spacing_dlnE" must be present (a value of 0.0 is legitimate).
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the raw (non-"Integrated_") dataset to check.
 * @param expected_lower_edge_ev The band lower edge this reader is about to
 * assume, in eV.
 */
static void radiation_check_band_edge_provenance(
    hid_t group_id, const char *dataset_name, double expected_lower_edge_ev) {

  const hid_t h_dataset = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
  if (h_dataset < 0)
    error("Error while opening dataset '%s' to check its provenance.",
          dataset_name);

  if (H5Aexists(h_dataset, "native_grid_spacing_dlnE") <= 0) {
    H5Dclose(h_dataset);
    error(
        "Data/Radiation/%s has no 'native_grid_spacing_dlnE' attribute: this "
        "table was generated before pychem recorded the edge value's own "
        "grid provenance and needs regenerating with pychem's "
        "pychem_generate_hdf5_parameters.",
        dataset_name);
  }

  char units[16];
  radiation_read_string_attribute(h_dataset, "lower_edge_energy_units", units,
                                  sizeof(units));
  if (strcmp(units, "eV") != 0) {
    H5Dclose(h_dataset);
    error(
        "Data/Radiation/%s declares lower_edge_energy_units='%s', but SWIFT "
        "assumes 'eV'.",
        dataset_name, units);
  }

  double lower_edge_energy_ev = 0.;
  io_read_attribute(h_dataset, "lower_edge_energy", DOUBLE,
                    &lower_edge_energy_ev);
  H5Dclose(h_dataset);

  if (fabs(lower_edge_energy_ev - expected_lower_edge_ev) >
      1e-6 * expected_lower_edge_ev) {
    error(
        "Data/Radiation/%s declares lower_edge_energy=%.6g eV, but SWIFT's "
        "own band split assumes %.6g eV: this table was generated for a "
        "different PE/LW band split than the code's, or is corrupted.",
        dataset_name, lower_edge_energy_ev, expected_lower_edge_ev);
  }
}

/**
 * @brief Read the L_PE (non-ionizing PE band emission rate) array from the
 * table.
 *
 * Only called when #radiation.with_ISRF is set.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_luminosity_pe_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(group_id, "L_PE", grid, sm, rad->interpolation_size,
                         units_cgs_conversion_factor(us, UNIT_CONV_POWER), 1.,
                         "erg/s", &rad->raw.l_pe, &rad->integrated.l_pe,
                         &rad->raw.l_pe_2d, &rad->integrated.l_pe_2d,
                         grid->edge_policy_l_pe);
}

/**
 * @brief Read the L_LW (Lyman-Werner band emission rate) array from the
 * table. Same as #radiation_read_luminosity_pe_array.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_luminosity_lw_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(group_id, "L_LW", grid, sm, rad->interpolation_size,
                         units_cgs_conversion_factor(us, UNIT_CONV_POWER), 1.,
                         "erg/s", &rad->raw.l_lw, &rad->integrated.l_lw,
                         &rad->raw.l_lw_2d, &rad->integrated.l_lw_2d,
                         grid->edge_policy_l_lw);
}

/**
 * @brief Read the SpectralPhotonRateAtPEEdge (PE band lower-edge spectral
 * photon rate dQ/dE) array from the table.
 *
 * Only called when #radiation.with_ISRF is set.
 *
 * The built table is E_lo(PE)^2 * dQ/dE in internal power units, comparable
 * to #radiation.raw.l_pe. The conversion factor includes
 * #RADIATION_PE_BAND_LOWER_EDGE_CGS^2. The expected units stay the dataset's
 * own "1/s/erg".
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_luminosity_edge_pe_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_check_band_edge_provenance(group_id, "SpectralPhotonRateAtPEEdge",
                                       RADIATION_PE_BAND_LOWER_EDGE_EV);

  const double conversion_factor =
      units_cgs_conversion_factor(us, UNIT_CONV_POWER) /
      (RADIATION_PE_BAND_LOWER_EDGE_CGS * RADIATION_PE_BAND_LOWER_EDGE_CGS);

  radiation_build_tables(group_id, "SpectralPhotonRateAtPEEdge", grid, sm,
                         rad->interpolation_size, conversion_factor, 1.,
                         "1/s/erg", &rad->raw.l_edge_pe,
                         &rad->integrated.l_edge_pe, &rad->raw.l_edge_pe_2d,
                         &rad->integrated.l_edge_pe_2d, grid->edge_policy_l_pe);
}

/**
 * @brief Read the SpectralPhotonRateAtLWEdge (LW band lower-edge spectral
 * photon rate dQ/dE) array from the table. Same as
 * #radiation_read_luminosity_edge_pe_array.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_luminosity_edge_lw_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_check_band_edge_provenance(group_id, "SpectralPhotonRateAtLWEdge",
                                       RADIATION_LW_BAND_LOWER_EDGE_EV);

  const double conversion_factor =
      units_cgs_conversion_factor(us, UNIT_CONV_POWER) /
      (RADIATION_LW_BAND_LOWER_EDGE_CGS * RADIATION_LW_BAND_LOWER_EDGE_CGS);

  radiation_build_tables(group_id, "SpectralPhotonRateAtLWEdge", grid, sm,
                         rad->interpolation_size, conversion_factor, 1.,
                         "1/s/erg", &rad->raw.l_edge_lw,
                         &rad->integrated.l_edge_lw, &rad->raw.l_edge_lw_2d,
                         &rad->integrated.l_edge_lw_2d, grid->edge_policy_l_lw);
}

/**
 * @brief Read the MeanPhotonEnergyLW / Integrated_MeanPhotonEnergyLW
 * (photon-number-weighted mean Lyman-Werner photon energy, L_LW/Q_LW over
 * 11.2-13.6 eV) arrays from the table.
 *
 * Both datasets are kept in cgs erg (no unit conversion or scaling), like
 * #radiation_get_mean_excess_photon_energy_HI_from_integral and
 * radiation_set_lw_photon_energy_cgs().
 *
 * Unlike every other Integrated_* dataset, Integrated_MeanPhotonEnergyLW
 * is not per Msun and not a cumulative integral: it is the ratio
 * Integrated_L_LW/Integrated_Q_LW. Both datasets are therefore read as
 * independent raw tables, in log10 space.
 *
 * Both axes clamp: the fallback below the LW mass floor is the 12.4 eV band
 * midpoint, so L_LW's "zero below" policy would give a zero photon energy.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 */
void radiation_read_mean_photon_energy_lw_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid,
    const struct stellar_model *sm) {

  radiation_build_tables(
      group_id, "MeanPhotonEnergyLW", grid, sm, rad->interpolation_size, 1., 1.,
      "erg", &rad->raw.mean_photon_energy_lw, NULL,
      &rad->raw.mean_photon_energy_lw_2d, NULL, boundary_condition_const);

  radiation_build_tables(group_id, "Integrated_MeanPhotonEnergyLW", grid, sm,
                         rad->interpolation_size, 1., 1., "erg",
                         &rad->integrated.mean_photon_energy_lw, NULL,
                         &rad->integrated.mean_photon_energy_lw_2d, NULL,
                         boundary_condition_const);
}

/**
 * @brief Read the main-sequence lifetime table (2D "M,Z" tables only).
 *
 * MainSequenceLifetime has no 1D analogue, so @p grid must be 2D.
 *
 * Read in Myr, not internal units, see
 * #radiation.raw.main_sequence_lifetime_2d.
 *
 * Both axes clamp: a star below (above) the tabulated mass range is very
 * long-lived (short-lived).
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata; must be a 2D ("M,Z") table.
 * @param sm The #stellar_model.
 * @param us The unit system (unused).
 */
void radiation_read_main_sequence_lifetime_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  if (!grid->is_2d) {
    error(
        "radiation_read_main_sequence_lifetime_array() called on a 1D "
        "('M') table: MainSequenceLifetime has no 1D analogue.");
  }

  const htri_t exists =
      H5Lexists(group_id, "MainSequenceLifetime", H5P_DEFAULT);
  if (exists <= 0) {
    error(
        "This 2D ('M,Z') Data/Radiation group has no 'MainSequenceLifetime' "
        "dataset. This table was generated before pychem added it and "
        "needs regenerating: run pychem's pychem_generate_hdf5_parameters "
        "on this table's own chimieparam file, then point "
        "GEARFeedback:yields_table (or yields_table_first_stars, for the "
        "PopIII model) at the regenerated file.");
  }

  radiation_build_tables(group_id, "MainSequenceLifetime", grid, sm,
                         rad->interpolation_size, 1., 1., "Myr", NULL, NULL,
                         &rad->raw.main_sequence_lifetime_2d, NULL,
                         boundary_condition_const);
}

/**
 * @brief Read the main-sequence-lifetime-inverse table (2D "M,Z" tables
 * only): "Age", "MainSequenceLifetimeInverse" and
 * "MainSequenceLifetimeInverseExcluded".
 *
 * Does not use #radiation_build_tables(): the output axis is age, not mass.
 * The Z axis keeps the native log10(Z) nodes and the age axis is the native
 * "Age" grid (identity resample, so no blending with the placeholder value
 * of Excluded cells). Each row's Excluded cells must be contiguous to the
 * end of the row. They are reduced to #radiation.longest_ms_lifetime_myr per
 * native Z row: the largest non-Excluded age, or FLT_MAX if none is Excluded.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata; must be a 2D ("M,Z") table.
 * @param sm The #stellar_model (unused).
 * @param us The unit system (unused).
 */
void radiation_read_main_sequence_lifetime_inverse_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  if (!grid->is_2d) {
    error(
        "radiation_read_main_sequence_lifetime_inverse_array() called on a "
        "1D ('M') table: MainSequenceLifetimeInverse has no 1D analogue.");
  }

  if (grid->n_metallicity > RADIATION_MAX_METALLICITY_ROWS) {
    error(
        "Data/Radiation's 'nz' attribute is %d, exceeding this build's "
        "RADIATION_MAX_METALLICITY_ROWS=%d. Increase that compile-time cap "
        "(stellar_evolution_struct.h) and rebuild to load this table.",
        grid->n_metallicity, RADIATION_MAX_METALLICITY_ROWS);
  }

  const htri_t age_exists = H5Lexists(group_id, "Age", H5P_DEFAULT);
  const htri_t msli_exists =
      H5Lexists(group_id, "MainSequenceLifetimeInverse", H5P_DEFAULT);
  const htri_t excl_exists =
      H5Lexists(group_id, "MainSequenceLifetimeInverseExcluded", H5P_DEFAULT);
  if (age_exists <= 0 || msli_exists <= 0 || excl_exists <= 0) {
    error(
        "This 2D ('M,Z') Data/Radiation group is missing one of 'Age'/"
        "'MainSequenceLifetimeInverse'/'MainSequenceLifetimeInverseExcluded'. "
        "This table was generated before pychem added the population "
        "main-sequence-lifetime cap and needs regenerating: run pychem's "
        "pychem_generate_hdf5_parameters on this table's own chimieparam "
        "file, then point GEARFeedback:yields_table (or "
        "yields_table_first_stars, for the PopIII model) at the "
        "regenerated file.");
  }

  double a0, da, age_max_myr;
  int na;
  io_read_attribute(group_id, "a0", DOUBLE, &a0);
  io_read_attribute(group_id, "da", DOUBLE, &da);
  io_read_attribute(group_id, "na", INT, &na);
  io_read_attribute(group_id, "age_max_myr", DOUBLE, &age_max_myr);

  if (!(da > 0.))
    error(
        "Data/Radiation's 'da' attribute is %.4g; it must be a strictly "
        "positive log10(age) grid step.",
        da);
  if (na < 2)
    error(
        "Data/Radiation's 'na' attribute is %d; at least 2 age points are "
        "needed to interpolate.",
        na);

  rad->age_max_myr = (float)age_max_myr;

  /* Native log10(Z) axis, as in radiation_build_tables(). */
  float *log_z_nodes = (float *)malloc(sizeof(float) * grid->n_metallicity);
  if (log_z_nodes == NULL)
    error(
        "Failed to allocate the RAD MainSequenceLifetimeInverse log10(Z) "
        "axis.");
  for (int i = 0; i < grid->n_metallicity; i++)
    log_z_nodes[i] = log10f(grid->metallicity[i]);

  /* Identity resample of the native age grid. */
  const float log_age_min = (float)a0;
  const float log_age_max = (float)(a0 + (na - 1) * da);

  const hsize_t count = (hsize_t)grid->n_metallicity * (hsize_t)na;

  float *log_data = (float *)malloc(sizeof(float) * count);
  if (log_data == NULL)
    error(
        "Failed to allocate the RAD MainSequenceLifetimeInverse log-value "
        "buffer.");
  float *data = radiation_read_cgs_array(
      group_id, "MainSequenceLifetimeInverse", count, 1., 1., "Msun", log_data);

  double *log_data_double = (double *)malloc(sizeof(double) * count);
  if (log_data_double == NULL)
    error(
        "Failed to allocate the RAD MainSequenceLifetimeInverse log-value "
        "double buffer.");
  for (hsize_t i = 0; i < count; i++) log_data_double[i] = (double)log_data[i];

  interpolate_2d_init(&rad->raw.main_sequence_lifetime_inverse_2d, log_z_nodes,
                      grid->n_metallicity, log_z_nodes, grid->n_metallicity,
                      log_age_min, log_age_max, na, log_age_min, (float)da, na,
                      log_data_double, boundary_condition_const,
                      boundary_condition_const);

  free(log_z_nodes);
  free(log_data_double);
  free(log_data);
  free(data);

  /* Reduce the Excluded mask to one scalar per native Z row. */
  hbool_t *excluded = (hbool_t *)malloc(sizeof(hbool_t) * count);
  if (excluded == NULL)
    error(
        "Failed to allocate the RAD MainSequenceLifetimeInverseExcluded "
        "buffer.");
  io_read_array_dataset(group_id, "MainSequenceLifetimeInverseExcluded", BOOL,
                        excluded, count);

  double *age = (double *)malloc(sizeof(double) * na);
  if (age == NULL) error("Failed to allocate the RAD Age buffer.");
  io_read_array_dataset(group_id, "Age", DOUBLE, age, na);

  for (int z = 0; z < grid->n_metallicity; z++) {
    int longest_index = -1;
    int seen_excluded = 0;
    for (int a = 0; a < na; a++) {
      const int is_excluded = excluded[(hsize_t)z * na + a] != 0;
      if (!is_excluded) {
        if (seen_excluded)
          error(
              "Data/Radiation's 'MainSequenceLifetimeInverseExcluded' row "
              "%d is not contiguous-True-to-the-end (a non-Excluded age "
              "follows an Excluded one at age index %d): the MS-lifetime "
              "cap's min()-gate safety proof requires this shape. "
              "Regenerate the table with pychem, or report this as a "
              "genuine pychem bug.",
              z, a);
        longest_index = a;
      } else {
        seen_excluded = 1;
      }
    }
    if (!seen_excluded) {
      rad->longest_ms_lifetime_myr[z] = FLT_MAX;
    } else if (longest_index >= 0) {
      rad->longest_ms_lifetime_myr[z] = (float)age[longest_index];
    } else {
      /* Every age of this row is Excluded: a threshold of 0 rejects every
         query for this row. */
      rad->longest_ms_lifetime_myr[z] = 0.f;
    }
  }

  free(age);
  free(excluded);
}

/**
 * @brief Open the "Data/Radiation" group of a yields table. The error names
 * the file and the regeneration step if the group is missing.
 *
 * @param filename The yields table filename to open (sm->yields_table).
 * @param file_id (output) The opened HDF5 file id.
 * @param group_id (output) The opened "Data/Radiation" group id.
 */
static void radiation_open_data_group(const char *filename, hid_t *file_id,
                                      hid_t *group_id) {

  *file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (*file_id < 0) error("unable to open file %s.\n", filename);

  const htri_t exists = H5Lexists(*file_id, "Data/Radiation", H5P_DEFAULT);
  if (exists <= 0) {
    error(
        "'%s' has no 'Data/Radiation' group. This yields table was "
        "generated before the radiation-table migration and needs "
        "regenerating: run pychem's pychem_generate_hdf5_parameters on "
        "this table's own chimieparam file, then point "
        "GEARFeedback:yields_table (or yields_table_first_stars, for the "
        "PopIII model) at the regenerated file.",
        filename);
  }

  *group_id = H5Gopen(*file_id, "Data/Radiation", H5P_DEFAULT);
  if (*group_id < 0)
    error("unable to open group 'Data/Radiation' in %s.\n", filename);
}

/**
 * @brief Read the RAD yields from the table.
 *
 * The tables are in internal units, except for a 2D table:
 * raw.main_sequence_lifetime_2d stays in Myr and
 * raw.main_sequence_lifetime_inverse_2d in Msun.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void radiation_read_data(struct radiation *rad, struct swift_params *params,
                         const struct stellar_model *sm,
                         const struct unit_system *us,
                         const struct phys_const *phys_const,
                         const int restart) {

  if (!restart) {
    rad->interpolation_size = parser_get_opt_param_int(
        params, "GEARRadiation:interpolation_size_mass", 500);
    if (rad->interpolation_size < 2) {
      error(
          "GEARRadiation:interpolation_size_mass must be >= 2; got "
          "%d.",
          rad->interpolation_size);
    }
  }

  /* radiation_zero_pointers() clears these, which are needed below: save and
     restore them. */
  const int interpolation_size_before = rad->interpolation_size;
  const int n_HII_pixels_before = rad->n_HII_pixels;
  /* with_ISRF is saved for the same reason. */
  const char with_ISRF_before = rad->with_ISRF;

  /* Zero every table so that radiation_clean() is safe whichever of the 1D
     and 2D variants is built. */
  radiation_zero_pointers(rad);

  rad->interpolation_size = interpolation_size_before;
  rad->n_HII_pixels = n_HII_pixels_before;
  rad->with_ISRF = with_ISRF_before;

  hid_t file_id, group_id;
  radiation_open_data_group(sm->yields_table, &file_id, &group_id);

  struct radiation_grid_metadata grid;
  radiation_read_grid_metadata(group_id, &grid);
  rad->is_2d = grid.is_2d;

  radiation_read_table_identity(rad, group_id, &grid);

  /* The ISRF needs L_PE, L_LW and their integrated datasets. */
  if (rad->with_ISRF) {
    const int has_l_pe = H5Lexists(group_id, "L_PE", H5P_DEFAULT) > 0;
    const int has_l_lw = H5Lexists(group_id, "L_LW", H5P_DEFAULT) > 0;
    const int has_integrated_l_pe =
        H5Lexists(group_id, "Integrated_L_PE", H5P_DEFAULT) > 0;
    const int has_integrated_l_lw =
        H5Lexists(group_id, "Integrated_L_LW", H5P_DEFAULT) > 0;
    if (!(has_l_pe && has_l_lw && has_integrated_l_pe && has_integrated_l_lw)) {
      error(
          "'%s': GEARFeedback:with_interstellar_radiation_field is on but this "
          "Data/Radiation group is missing%s%s%s%s. Regenerate the table "
          "with pychem's pychem_generate_hdf5_parameters on its own "
          "chimieparam file.",
          sm->yields_table, has_l_pe ? "" : " 'L_PE'",
          has_l_lw ? "" : " 'L_LW'",
          has_integrated_l_pe ? "" : " 'Integrated_L_PE'",
          has_integrated_l_lw ? "" : " 'Integrated_L_LW'");
    }
  }

  /* The band-edge datasets are required, like the datasets above. */
  if (rad->with_ISRF) {
    const int has_l_edge_pe =
        H5Lexists(group_id, "SpectralPhotonRateAtPEEdge", H5P_DEFAULT) > 0;
    const int has_l_edge_lw =
        H5Lexists(group_id, "SpectralPhotonRateAtLWEdge", H5P_DEFAULT) > 0;
    const int has_integrated_l_edge_pe =
        H5Lexists(group_id, "Integrated_SpectralPhotonRateAtPEEdge",
                  H5P_DEFAULT) > 0;
    const int has_integrated_l_edge_lw =
        H5Lexists(group_id, "Integrated_SpectralPhotonRateAtLWEdge",
                  H5P_DEFAULT) > 0;
    if (!(has_l_edge_pe && has_l_edge_lw && has_integrated_l_edge_pe &&
          has_integrated_l_edge_lw)) {
      error(
          "'%s': GEARFeedback:with_interstellar_radiation_field is on but "
          "this Data/Radiation group is missing%s%s%s%s. Regenerate the "
          "table with pychem's pychem_generate_hdf5_parameters on its own "
          "chimieparam file.",
          sm->yields_table,
          has_l_edge_pe ? "" : " 'SpectralPhotonRateAtPEEdge'",
          has_l_edge_lw ? "" : " 'SpectralPhotonRateAtLWEdge'",
          has_integrated_l_edge_pe ? ""
                                   : " 'Integrated_SpectralPhotonRateAtPEEdge'",
          has_integrated_l_edge_lw
              ? ""
              : " 'Integrated_SpectralPhotonRateAtLWEdge'");
    }
  }

  /* Required for every table. */
  {
    const int has_mean_photon_energy_lw =
        H5Lexists(group_id, "MeanPhotonEnergyLW", H5P_DEFAULT) > 0;
    const int has_integrated_mean_photon_energy_lw =
        H5Lexists(group_id, "Integrated_MeanPhotonEnergyLW", H5P_DEFAULT) > 0;
    if (!(has_mean_photon_energy_lw && has_integrated_mean_photon_energy_lw)) {
      error(
          "'%s': Data/Radiation group is missing%s%s. Regenerate the table "
          "with pychem's pychem_generate_hdf5_parameters on its own "
          "chimieparam file.",
          sm->yields_table,
          has_mean_photon_energy_lw ? "" : " 'MeanPhotonEnergyLW'",
          has_integrated_mean_photon_energy_lw
              ? ""
              : " 'Integrated_MeanPhotonEnergyLW'");
    }
  }

  /* A no-op on a table without precomputed IMF-integrated datasets. */
  radiation_check_imf_consistency(group_id, sm);

  /* Table-coverage check: a star outside the table's mass grid would
     silently get the nearest edge's value. The tolerance of half a grid cell
     absorbs float rounding of the two independently computed bounds. */
  const float log_mass_min_imf = log10f(sm->imf.mass_min);
  const float log_mass_max_imf = log10f(sm->imf.mass_max);
  const double log_mass_max_table =
      (double)grid.log_mass_min + (grid.n_mass - 1) * (double)grid.mass_step;
  const double tol = 0.5 * (double)grid.mass_step;
  if ((double)log_mass_min_imf < (double)grid.log_mass_min - tol ||
      (double)log_mass_max_imf > log_mass_max_table + tol) {
    error(
        "'%s': Data/Radiation's native mass grid [%.4g, %.4g] Msun does not "
        "cover the IMF's mass range [%.4g, %.4g] Msun. A star outside the "
        "table's own range would silently receive the nearest edge's "
        "photon budget (boundary_condition_const) instead of its own "
        "value. Regenerate the table over (at least) the IMF's mass range "
        "with pychem, or adjust the IMF's own mass_min/mass_max to fit "
        "inside the table.",
        sm->yields_table, (double)exp10(grid.log_mass_min),
        (double)exp10(log_mass_max_table), (double)sm->imf.mass_min,
        (double)sm->imf.mass_max);
  }

  if (rad->is_2d && engine_rank == 0) {
    message(
        "'%s' is a mass x metallicity (2D) radiation table; both "
        "individual-star and population (SSP/continuous_IMF) feedback are "
        "supported.",
        sm->yields_table);
  }

  /* Read the luminosities */
  radiation_read_luminosities_array(rad, group_id, &grid, sm, us);

  /* Read the ionization emission rates */
  radiation_read_ionization_rate_array(rad, group_id, &grid, sm, us);

  /* Read the excess-photon-energy emission rates */
  radiation_read_mean_excess_photon_energy_array(rad, group_id, &grid, sm, us);

  /* Optional: without Teff every star reports 0. */
  rad->has_teff = (char)(H5Lexists(group_id, "Teff", H5P_DEFAULT) > 0);
  if (rad->has_teff) radiation_read_teff_array(rad, group_id, &grid, sm, us);

  /* Validated above to exist. */
  if (rad->with_ISRF) {
    radiation_read_luminosity_pe_array(rad, group_id, &grid, sm, us);
    radiation_read_luminosity_lw_array(rad, group_id, &grid, sm, us);
    radiation_read_luminosity_edge_pe_array(rad, group_id, &grid, sm, us);
    radiation_read_luminosity_edge_lw_array(rad, group_id, &grid, sm, us);
  }

  /* Validated above to exist. */
  radiation_read_mean_photon_energy_lw_array(rad, group_id, &grid, sm);

  /* These datasets exist only for 2D tables. */
  if (grid.is_2d) {
    radiation_read_main_sequence_lifetime_array(rad, group_id, &grid, sm, us);
    radiation_read_main_sequence_lifetime_inverse_array(rad, group_id, &grid,
                                                        sm, us);
  }

  free(grid.metallicity);
  h5_close_group(file_id, group_id);

  /* Mark the tables valid. */
  rad->is_active = 1;
}
