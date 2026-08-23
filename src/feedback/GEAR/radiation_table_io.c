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
 * feedback: the Data/Radiation group's attributes and grid metadata, the
 * generic CGS-to-internal-units array reader, table construction (raw and
 * IMF-integrated), and the four public radiation_read_*_array() entry
 * points plus radiation_read_data() itself.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "engine.h"
#include "error.h"
#include "hdf5_functions.h"
#include "inline.h"
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
 * SWIFT's generic io_read_attribute() (common_io.h) only supports the
 * numeric/char/bool #IO_DATA_TYPE variants; there is no read counterpart to
 * common_io.c's write-only io_writeStringAttribute() for a python/h5py
 * variable-length UTF-8 string, which is what pychem writes for
 * Data/Radiation's "dimensionality" attribute. This is the minimal reader
 * needed to dispatch on it.
 *
 * @param group_id Open HDF5 group id.
 * @param name Attribute name.
 * @param out Output buffer, NUL-terminated on return.
 * @param out_size Size of out, including the terminating NUL.
 */
static void radiation_read_string_attribute(hid_t group_id, const char *name,
                                            char *out, size_t out_size) {

  const hid_t h_attr = H5Aopen(group_id, name, H5P_DEFAULT);
  if (h_attr < 0) error("Error while opening attribute '%s'", name);

  const hid_t h_type = H5Aget_type(h_attr);
  if (h_type < 0) error("Error while getting the type of attribute '%s'", name);

  if (H5Tis_variable_str(h_type) > 0) {
    /* Read into the attribute's own native type (variable-length,
       H5T_CSET_UTF8, as pychem/h5py writes it) rather than a freshly
       crafted H5T_C_S1/H5T_VARIABLE memory type: the two differ in
       character set (ASCII vs UTF-8), and this HDF5 build has no
       registered ASCII<->UTF-8 conversion path, so H5Aread() into the
       mismatched type fails ("no appropriate function for conversion
       path") -- caught by actually running this against a real pychem
       table, not just compiling it. */
    char *tmp = NULL;
    if (H5Aread(h_attr, h_type, &tmp) < 0)
      error("Error while reading string attribute '%s'", name);

    const size_t len = strlen(tmp);
    if (len >= out_size) {
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
    if (fixed_size >= out_size)
      error(
          "String attribute '%s' (%zu bytes) does not fit in the %zu-byte "
          "buffer.",
          name, fixed_size, out_size);

    char *tmp = (char *)calloc(fixed_size + 1, sizeof(char));
    if (tmp == NULL) error("Failed to allocate string attribute buffer.");
    if (H5Aread(h_attr, h_type, tmp) < 0)
      error("Error while reading string attribute '%s'", name);
    memcpy(out, tmp, fixed_size);
    out[fixed_size] = '\0';
    free(tmp);
  }

  H5Tclose(h_type);
  H5Aclose(h_attr);
}

/**
 * @brief Assert that a Data/Radiation dataset's own "units" attribute
 * matches what this reader is about to assume.
 *
 * pychem writes an explicit "units" attribute on every dataset in
 * Data/Radiation specifically so a unit mismatch is never silently
 * implicit (PyChemInitTable/libradiation.py's own write_h5_table()
 * docstring: "Grackle's unit conventions are a known source of silent
 * errors, so units are never left implicit here"). SWIFT's reader
 * (radiation_read_cgs_array()) converts every dataset with a fixed,
 * hardcoded physical-dimension assumption (CGS erg/s for Luminosity, CGS
 * 1/s for Q_H, CGS erg/s for DotEExcess); this check makes that assumption
 * self-verifying against the table itself instead of trusting it blindly,
 * the same fail-loud-on-mismatch policy this file already applies to
 * float overflow (radiation_read_cgs_array()) and table-coverage (the
 * inline check in radiation_read_data()).
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset whose "units" attribute to check.
 * @param expected_units The unit string this reader is about to assume
 * (e.g. "erg/s"), compared verbatim (case-sensitive) against the table's
 * own attribute.
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
        "code converts -- aborting rather than silently misinterpreting "
        "the physics.",
        dataset_name, actual_units, expected_units);
  }
}

/**
 * @brief Turn one field's edge_policy_<field>_below/above attribute pair
 * into the #interpolate_boundary_condition SWIFT's 2D interpolator can
 * apply on the mass axis.
 *
 * pychem's schema allows "constant", "zero" or "linear" independently per
 * side (PyChemInitTable/libradiation.py's ExtrapolationPolicy); SWIFT has
 * no "linear" extrapolation mode, and #interpolate_boundary_condition can
 * only express "zero both sides", "zero below/constant above" or "constant
 * both sides" -- not "constant below/zero above". Every real shipped
 * table so far only uses the three representable combinations (see
 * radiation_parsec_popIII.hdf5's edge_policy_* attributes); this errors
 * loudly on anything else rather than silently misapplying a policy.
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
 * @brief Read the Data/Radiation group's own grid metadata: the
 * "dimensionality" attribute ("M" or "M,Z") and the group-level
 * "m0"/"dm"/"nm" mass-grid attributes shared by every dataset in the
 * group, plus, for a 2D ("M,Z") table, "nz", the "Metallicity" dataset,
 * and the mass-axis edge_policy_* attributes (dispatched through "source"
 * for Q_H; see radiation_parse_edge_policy()).
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid (output) The #radiation_grid_metadata to fill in.
 */
static void radiation_read_grid_metadata(hid_t group_id,
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
    }

    char source[32];
    radiation_read_string_attribute(group_id, "source", source, sizeof(source));

    char luminosity_below[16], luminosity_above[16];
    char q_h_blackbody_below[16], q_h_blackbody_above[16];
    char q_h_parsec_below[16], q_h_parsec_above[16];
    radiation_read_string_attribute(group_id, "edge_policy_luminosity_below",
                                    luminosity_below, sizeof(luminosity_below));
    radiation_read_string_attribute(group_id, "edge_policy_luminosity_above",
                                    luminosity_above, sizeof(luminosity_above));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_blackbody_below",
                                    q_h_blackbody_below,
                                    sizeof(q_h_blackbody_below));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_blackbody_above",
                                    q_h_blackbody_above,
                                    sizeof(q_h_blackbody_above));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_parsec_below",
                                    q_h_parsec_below, sizeof(q_h_parsec_below));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_parsec_above",
                                    q_h_parsec_above, sizeof(q_h_parsec_above));

    grid->edge_policy_luminosity = radiation_parse_edge_policy(
        luminosity_below, luminosity_above, "luminosity");
    const enum interpolate_boundary_condition edge_policy_q_h_blackbody =
        radiation_parse_edge_policy(q_h_blackbody_below, q_h_blackbody_above,
                                    "q_h_blackbody");
    const enum interpolate_boundary_condition edge_policy_q_h_parsec =
        radiation_parse_edge_policy(q_h_parsec_below, q_h_parsec_above,
                                    "q_h_parsec");
    grid->edge_policy_dot_e_excess = edge_policy_q_h_blackbody;

    if (strcmp(source, "parsec_blackbody") == 0) {
      grid->edge_policy_q_h = edge_policy_q_h_blackbody;
    } else if (strcmp(source, "parsec_qtable") == 0) {
      grid->edge_policy_q_h = edge_policy_q_h_parsec;
    } else {
      error(
          "Data/Radiation's 'source' attribute is '%s' (expected "
          "'parsec_blackbody' or 'parsec_qtable'): cannot tell which "
          "Q_H edge policy applies to the table's primary 'Q_H' dataset.",
          source);
    }
  } else if (strcmp(grid->dimensionality, "M") == 0) {
    grid->is_2d = 0;
    grid->edge_policy_luminosity = boundary_condition_error;
    grid->edge_policy_q_h = boundary_condition_error;
    grid->edge_policy_dot_e_excess = boundary_condition_error;
  } else {
    error(
        "Data/Radiation has an unrecognised 'dimensionality' attribute "
        "'%s' (expected 'M' or 'M,Z').",
        grid->dimensionality);
  }
}

/**
 * @brief Cross-check a table's "imf_a_s"/"imf_m_s"/"mass_min_msun"/
 * "mass_max_msun" attributes against @p sm's own #initial_mass_function,
 * for a table carrying pychem's precomputed IMF-integrated datasets.
 *
 * pychem's Integrated_* cumulative datasets (see radiation_build_tables())
 * were computed against a specific IMF; if that IMF differs from the one
 * this run configured (#sm->imf), the precomputed integral is silently
 * wrong for this run -- exactly the class of bug this migration exists to
 * close. Only tables with the precomputed datasets write these attrs, so a
 * table without "Integrated_Q_H" is skipped here (not an error: an
 * old-format table is instead caught when a reader tries to open
 * "Integrated_Q_H" directly, in radiation_build_tables()). "imf_m_s" holds
 * only the IMF's *interior* breakpoints, so it is compared against
 * #initial_mass_function.mass_limits[1 .. n_parts - 1], with
 * mass_min_msun/mass_max_msun standing in for its outer two entries.
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

  /* A segment-count mismatch is not a rounding question: check it here,
     with the same actionable framing as the value mismatches below,
     before io_read_array_attribute()'s own generic "different number of
     elements than expected" error would fire instead. */
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
 * @brief Read one CGS-valued dataset from an open Data/Radiation group,
 * convert it to internal units and (for Q_H/DotEExcess)
 * #RADIATION_DOT_N_ION_TABLE_SCALING, narrow it to float, and (optionally)
 * also compute its log10, pychem-style, for the caller's raw (log-log)
 * interpolation table.
 *
 * Shared by radiation_read_luminosities_array(),
 * radiation_read_ionization_rate_array() and
 * radiation_read_mean_excess_photon_energy_array() to avoid tripling the
 * read/convert/guard boilerplate.
 *
 * Guards against float overflow the way stellar_evolution.c:676-695 does:
 * error() aborts (MPI_Abort/swift_abort, src/error.h) rather than capping
 * -- a units/scaling bug should stop the run, not silently corrupt the
 * physics. Also flags (debug-checks only) an implausible collapse to
 * exactly zero for a CGS input that was not itself zero: pychem bakes a
 * literal 0 into Q_H/DotEExcess below its own ionization threshold, so an
 * exact-zero result is only suspicious when the source value was nonzero.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to read.
 * @param count Number of elements to read (the group's "nm", or "nm" *
 * "nz" for a 2D table).
 * @param conversion_factor units_cgs_conversion_factor() for this
 * dataset's physical dimension; CGS values are divided by this to reach
 * internal units (SWIFT's convention).
 * @param extra_scaling Additional SWIFT-side-only divisor applied after
 * unit conversion (#RADIATION_DOT_N_ION_TABLE_SCALING for Q_H/DotEExcess,
 * 1 for Luminosity).
 * @param expected_units The dataset's own "units" attribute is asserted to
 * equal this string before any conversion happens (see
 * radiation_check_dataset_units()).
 * @param log_data_internal (output, optional) If not NULL, a caller-owned
 * float array of length @p count filled with log10 of the same
 * internal-unit value #RADIATION_LOG_FLOOR_CGS-floored on the CGS side
 * before conversion (see that macro's own doxygen) -- exactly pychem's
 * log-log convention, used to build a raw table's #interpolation_1d /
 * #interpolation_2d in log-value space instead of raw-value space. Left
 * untouched if NULL (the integrated-table caller has no use for it: see
 * radiation_build_tables()'s own doxygen for why the IMF-integrated table
 * stays in linear space).
 * @return Newly malloc'd float array of length count, in internal
 * (optionally rescaled) units, NOT logged -- this is the value a raw
 * dataset read needs (@p log_data_internal built alongside it) or, with
 * @p log_data_internal NULL, the linear-space value an "Integrated_*"
 * cumulative-table read needs. Caller must free().
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

  /* log10(internal value) = log10(cgs value) - log10(conversion_factor *
     extra_scaling); computed once here rather than per-entry below. */
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
          "collapsed to exactly 0 in internal units after conversion -- "
          "check RADIATION_DOT_N_ION_TABLE_SCALING and the unit system.",
          dataset_name, (unsigned long long)j, data_cgs[j]);
    }
#endif

    data[j] = (float)value_internal;

    if (log_data_internal != NULL) {
      const double floored_cgs = max(data_cgs[j], RADIATION_LOG_FLOOR_CGS);
      const double log_value_internal = log10(floored_cgs) - log_conversion;
      log_data_internal[j] = (float)log_value_internal;

#ifdef SWIFT_DEBUG_CHECKS
      /* Self-check the log10/exp10 round-trip against the independently
         computed linear value above, away from the floor (where the two
         are expected to diverge by construction): a sign error or a
         dropped extra_scaling term in log_conversion would silently make
         every raw getter wrong by many orders of magnitude while still
         running to completion, so this is checked at load time on every
         run rather than trusted from inspection alone. */
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
 * @brief Build the raw (and, for a 1D table, IMF-integrated) interpolation
 * table for one Data/Radiation quantity, dispatching on the table's
 * dimensionality.
 *
 * The raw table (1D or 2D) is built in log10(value) space, floored and
 * exponentiated pychem-style (see #RADIATION_LOG_FLOOR_CGS and
 * radiation_read_cgs_array()): #interpolate_1d_init()/#interpolate_2d_init()
 * are otherwise-unmodified generic linear interpolators, so feeding them
 * already-logged data makes their existing linear interpolation a log-log
 * interpolation for free. Every raw getter (radiation_get_*_from_raw())
 * must exponentiate the result back; see their own doxygen.
 *
 * The IMF-integrated table (1D only) is deliberately left in linear
 * (un-logged) value space, unlike the raw table above: it is not the same
 * kind of quantity pychem's log-log scheme governs. It is read directly
 * from the table's own "Integrated_<dataset_name>" dataset -- pychem's
 * precomputed, number-weighted (n(m), not #initial_mass_function_integrate()'s
 * mass-weighted m*n(m)), cumulative-from-Mmin integral -- rather than
 * integrated on the SWIFT side, so no logged/un-logged ordering constraint
 * applies to it the way it does to the raw table above.
 *
 * The 2D ("M,Z") branch only builds the raw table: no caller integrates a
 * 2D table over the IMF yet (see #radiation's own doxygen), so
 * @p integrated_1d is left untouched -- radiation_read_data() zeroes it
 * beforehand so radiation_clean() stays safe either way. It also
 * approximates the metallicity axis as log-uniformly spaced from the
 * "Metallicity" dataset's first/last values: pychem does not guarantee
 * this (it is the curated set of PARSEC metallicities actually collapsed
 * into the table, not a synthetic grid), but interpolate_2d_init()
 * requires a uniform grid. Acceptable only because this 2D path has zero
 * test coverage until the PARSEC validation track (plan Phase 6)
 * exercises it.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to read.
 * @param grid The group's own grid metadata (see
 * radiation_read_grid_metadata()).
 * @param sm The #stellar_model (for the output mass-grid bounds).
 * @param interpolation_size_mass Number of points in the mass
 * interpolation output grid.
 * @param interpolation_size_metallicity Number of points in the
 * metallicity interpolation output grid (2D tables only).
 * @param conversion_factor See radiation_read_cgs_array().
 * @param extra_scaling See radiation_read_cgs_array().
 * @param expected_units See radiation_read_cgs_array().
 * @param raw_1d (output) Raw 1D interpolation table (1D tables only),
 * holding log10(value in internal units), pychem-floored.
 * @param integrated_1d (output) IMF-integrated 1D interpolation table (1D
 * tables only), holding the linear (un-logged) cumulative value -- see
 * this function's own doxygen for why.
 * @param raw_2d (output) Raw 2D interpolation table (2D tables only),
 * holding log10(value in internal units), pychem-floored.
 * @param boundary_condition_mass Mass-axis #interpolate_boundary_condition
 * for @p dataset_name (2D tables only; ignored for a 1D table, which always
 * clamps -- see radiation_read_grid_metadata()'s edge_policy_* fields for
 * where this comes from). The metallicity axis always clamps
 * (boundary_condition_const), matching pychem's "clamp to nearest grid Z,
 * never extrapolated" convention.
 */
static void radiation_build_tables(
    hid_t group_id, const char *dataset_name,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    int interpolation_size_mass, int interpolation_size_metallicity,
    double conversion_factor, double extra_scaling, const char *expected_units,
    struct interpolation_1d *raw_1d, struct interpolation_1d *integrated_1d,
    struct interpolation_2d *raw_2d,
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

    const float log_z_min = log10f(grid->metallicity[0]);
    const float log_z_max = log10f(grid->metallicity[grid->n_metallicity - 1]);
    const float log_z_step =
        grid->n_metallicity > 1
            ? (log_z_max - log_z_min) / (grid->n_metallicity - 1)
            : 0.f;

    /* interpolate_2d_init() takes a double source array (its internal
       storage is float; see interpolation.h); re-widen the already
       guarded/narrowed/logged float data rather than duplicating the guard
       for a double codepath. */
    double *log_data_double = (double *)malloc(sizeof(double) * count);
    if (log_data_double == NULL)
      error("Failed to allocate the RAD 2D log-value yields for %s.",
            dataset_name);
    for (hsize_t i = 0; i < count; i++)
      log_data_double[i] = (double)log_data[i];

    interpolate_2d_init(raw_2d, log_z_min, log_z_max,
                        interpolation_size_metallicity, log_mass_min_out,
                        log_mass_max_out, interpolation_size_mass, log_z_min,
                        grid->log_mass_min, log_z_step, grid->mass_step,
                        grid->n_metallicity, grid->n_mass, log_data_double,
                        boundary_condition_const, boundary_condition_mass);

    free(log_data_double);
    free(log_data);
    free(data);
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

  /* integrated_1d is built from pychem's own precomputed, number-weighted,
     cumulative-from-Mmin "Integrated_<dataset_name>" dataset, not from
     integrating the raw values above -- see this function's own doxygen. */
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

  /* "per Msun of stars formed" is a fixed physical mass unit, not a SWIFT
     internal-unit-system quantity, so the SAME conversion_factor/
     extra_scaling this dataset's raw sibling already uses apply unchanged
     to the "erg/s" or "1/s" part of the stored value. */
  char integrated_expected_units[32];
  written =
      snprintf(integrated_expected_units, sizeof(integrated_expected_units),
               "%s/Msun", expected_units);
  if (written < 0 || (size_t)written >= sizeof(integrated_expected_units))
    error("Units string '%s/Msun' does not fit in the buffer.", expected_units);

  /* log_data_internal = NULL: the cumulative integral stays in linear
     (un-logged) value space, unlike the raw table -- see this function's
     own doxygen for why. */
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
      rad->interpolation_size_metallicity,
      units_cgs_conversion_factor(us, UNIT_CONV_POWER), 1., "erg/s",
      &rad->raw.luminosities, &rad->integrated.luminosities,
      &rad->raw.luminosities_2d, grid->edge_policy_luminosity);
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
      rad->interpolation_size_metallicity,
      units_cgs_conversion_factor(us, UNIT_CONV_PHOTONS_PER_TIME),
      RADIATION_DOT_N_ION_TABLE_SCALING, "1/s", &rad->raw.dot_N_ion,
      &rad->integrated.dot_N_ion, &rad->raw.dot_N_ion_2d,
      grid->edge_policy_q_h);
}

/**
 * @brief Read an array of excess-photon-energy emission rate data from the
 * table: DotEExcess(m) = Q_H(m) * MeanExcessPhotonEnergyHI(m).
 *
 * Converted with the SAME (rate-only) #UNIT_CONV_PHOTONS_PER_TIME factor
 * used for Q_H, not #UNIT_CONV_POWER -- deliberately, even though the file
 * stores DotEExcess in erg/s (a power): dividing only the rate part by
 * #RADIATION_DOT_N_ION_TABLE_SCALING and the unit conversion, while
 * leaving the "erg" part of the product in cgs, reproduces the mixed-unit
 * convention #radiation_get_mean_excess_photon_energy_HI_from_integral (an
 * existing, unmodified function/call site) already relies on: both raw
 * tables share the same rate-only scaling, so it cancels exactly in that
 * ratio, and the result comes out in cgs erg -- matching
 * feedback_struct.h's documented cgs-erg convention for
 * mean_excess_photon_energy_HI, without needing a unit_system argument on
 * the getters. Using #UNIT_CONV_POWER here instead would leave that ratio
 * in internal energy units, silently changing the existing getter's
 * output. The units check below still asserts "erg/s" (not the rate-only
 * factor's implied "1/s"): it verifies the table's stored unit, which the
 * deliberate mismatch above depends on staying exactly "erg/s" for the
 * cancellation to hold.
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
      rad->interpolation_size_metallicity,
      units_cgs_conversion_factor(us, UNIT_CONV_PHOTONS_PER_TIME),
      RADIATION_DOT_N_ION_TABLE_SCALING, "erg/s", &rad->raw.dot_E_excess,
      &rad->integrated.dot_E_excess, &rad->raw.dot_E_excess_2d,
      grid->edge_policy_dot_e_excess);
}

/**
 * @brief Read the main-sequence lifetime table (2D "M,Z" tables only).
 *
 * MainSequenceLifetime has no 1D/"M"-table analogue -- pychem only writes
 * it for a PARSEC table -- so the caller must only call this on a 2D grid
 * (checked below); see radiation_read_data()'s own #radiation_grid_metadata
 * .is_2d gate around this call.
 *
 * Read in Myr, not internal units: #conversion_factor and #extra_scaling
 * are both the identity (1.0), so the table's own "Myr" values are stored
 * as-is rather than run through units_cgs_conversion_factor(us,
 * UNIT_CONV_TIME), which would treat the stored value as if it were
 * already in CGS seconds -- see #radiation's raw.main_sequence_lifetime_2d
 * doxygen for why staying in Myr is the deliberate choice here.
 *
 * The mass-axis boundary condition is hardcoded to boundary_condition_const
 * (clamp to the nearest tabulated lifetime at both ends) rather than
 * dispatched from a per-field edge_policy_* attribute: pychem's schema has
 * no edge_policy_main_sequence_lifetime_* pair (only Luminosity/Teff/
 * Radius/Q_H_Blackbody/Q_H_PARSEC have one). Clamping is the physically
 * sane default in both directions: a star below the table's lowest
 * tabulated mass is very long-lived, one above the highest is very
 * short-lived, and zero-extrapolating either would give a nonsensical
 * instantaneous cap.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata; must be a 2D ("M,Z") table.
 * @param sm The #stellar_model.
 * @param us The unit system; unused (this table is Myr-native, see above),
 * kept only for signature symmetry with the other three _array() readers.
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

  radiation_build_tables(
      group_id, "MainSequenceLifetime", grid, sm, rad->interpolation_size,
      rad->interpolation_size_metallicity, 1., 1., "Myr", NULL, NULL,
      &rad->raw.main_sequence_lifetime_2d, boundary_condition_const);
}

/**
 * @brief Open the "Data/Radiation" group of a yields table, with an
 * actionable error naming the file and the regeneration step if the group
 * is missing.
 *
 * Every yields table generated before this migration (the ones currently
 * committed/fetched for the 8 radiation examples included) lacks this
 * group, and #h5_open_group's own generic "unable to open group" message
 * gives no hint that regenerating the table, not fixing a typo, is the
 * actual fix -- mirrors the actionable-message convention radiation_init()
 * (radiation.c) already uses for GEARFeedback:HII_angular_nside.
 *
 * @param filename The yields table filename (sm->yields_table or
 * sm->yields_table for the first-stars model -- same check either way).
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
 * The tables are in internal units at the end of this function, with one
 * exception: for a 2D ("M,Z") table, raw.main_sequence_lifetime_2d stays
 * in Myr, deliberately not run through the unit system -- see its own
 * doxygen on #radiation's raw sub-struct for why.
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
    rad->interpolation_size_metallicity = parser_get_opt_param_int(
        params, "GEARRadiation:interpolation_size_metallicity", 110);
    if (rad->interpolation_size < 2) {
      error(
          "GEARRadiation:interpolation_size_mass must be >= 2; got "
          "%d.",
          rad->interpolation_size);
    }
    if (rad->interpolation_size_metallicity < 2) {
      error(
          "GEARRadiation:interpolation_size_metallicity must be >= "
          "2; got %d.",
          rad->interpolation_size_metallicity);
    }
  }

  /* radiation_zero_pointers() below also clears interpolation_size(_
     metallicity) and n_HII_pixels, so callers with radiation disabled
     don't inherit uninitialized garbage in them -- but on this call path
     those three either were just set above (!restart) or hold the values
     radiation_restore() flat-restored moments ago (restart), and
     radiation_build_tables() below needs them either way. Round-trip them
     around the call. */
  const int interpolation_size_before = rad->interpolation_size;
  const int interpolation_size_metallicity_before =
      rad->interpolation_size_metallicity;
  const int n_HII_pixels_before = rad->n_HII_pixels;

  /* Zero every table up front: radiation_build_tables() only populates the
     _1d or _2d variant matching this table's dimensionality, and only the
     2D path's raw tables at that (see its own doxygen) -- the rest must be
     safe no-ops for radiation_clean()'s interpolate_1d_free()/
     interpolate_2d_free() calls regardless of which branch ran. */
  radiation_zero_pointers(rad);

  rad->interpolation_size = interpolation_size_before;
  rad->interpolation_size_metallicity = interpolation_size_metallicity_before;
  rad->n_HII_pixels = n_HII_pixels_before;

  hid_t file_id, group_id;
  radiation_open_data_group(sm->yields_table, &file_id, &group_id);

  struct radiation_grid_metadata grid;
  radiation_read_grid_metadata(group_id, &grid);
  rad->is_2d = grid.is_2d;

  /* A no-op on a table without pychem's precomputed IMF-integrated
     datasets; see radiation_check_imf_consistency()'s own doxygen. Runs
     for both 1D and 2D tables. */
  radiation_check_imf_consistency(group_id, sm);

  /* Table-coverage check: the raw table's boundary condition is
     boundary_condition_const (radiation_build_tables()), so a star outside
     the table's own native mass grid silently gets the nearest edge's
     value instead of its own -- e.g. pychem's real PopIII table has a
     native floor of 13 Msun, well above a typical IMF's own mass_min.
     Fail loudly instead, matching this file's existing convention (the
     FLT_MAX guard above) and radiation_init()'s HII_angular_nside checks
     (radiation.c).

     Compared with a half-grid-cell tolerance, not exact equality: both
     sides are independently accumulated (log_mass_min_imf/log_mass_max_imf
     from log10f() of the IMF's own bounds, log_mass_max_table from the
     grid's own log_mass_min/mass_step/n_mass), so a star whose mass range
     was generated to sit exactly on the table's own edge can differ from
     it by a few ULP of float rounding. A tolerance-free comparison rejects
     the overwhelming majority of otherwise-valid (mass_min, mass_max, nm)
     combinations for no physical reason; a genuine shortfall (e.g. the
     documented 13 Msun PopIII floor case) is orders of magnitude outside
     this tolerance and still aborts. */
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

  /* Fail at load time, not at the first star's feedback computation: every
     radiation_get_*_from_raw()/_from_integral() getter aborts on a 2D
     table (see radiation_check_dimensionality()), since no caller passes a
     metallicity yet. */
  if (rad->is_2d && engine_rank == 0) {
    message(
        "WARNING: '%s' is a mass x metallicity (2D) radiation table; no "
        "individual-star or population feedback getter supports one yet -- "
        "any star reaching feedback will abort the run.",
        sm->yields_table);
  }

  /* Read the luminosities */
  radiation_read_luminosities_array(rad, group_id, &grid, sm, us);

  /* Read the ionization emission rates */
  radiation_read_ionization_rate_array(rad, group_id, &grid, sm, us);

  /* Read the excess-photon-energy emission rates */
  radiation_read_mean_excess_photon_energy_array(rad, group_id, &grid, sm, us);

  /* MainSequenceLifetime has no 1D ("M") table analogue: only read it for
     a 2D table, where the HDF5 dataset actually exists. */
  if (grid.is_2d) {
    radiation_read_main_sequence_lifetime_array(rad, group_id, &grid, sm, us);
  }

  free(grid.metallicity);
  h5_close_group(file_id, group_id);

  /* The tables above are now valid: mark this #radiation active so callers
     use them instead of skipping to their zeroed defaults. */
  rad->is_active = 1;
};
