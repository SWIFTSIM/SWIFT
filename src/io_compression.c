/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020  Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

/* This object's header. */
#include "io_compression.h"

/* Local includes. */
#include "error.h"

/* Some standard headers. */
#include <stdlib.h>
#include <string.h>

/**
 * @brief Names of the compression levels, used in the select_output.yml
 *        parameter file.
 **/
const char* lossy_compression_schemes_names[compression_level_count] = {
    "off",        "on",          "DScale1",   "DScale2",    "DScale3",
    "DScale4",    "DScale5",     "DScale6",   "DMantissa9", "DMantissa13",
    "FMantissa9", "FMantissa13", "HalfFloat", "BFloat16",   "Nbit32",
    "Nbit36",     "Nbit40",      "Nbit44",    "Nbit48",     "Nbit56"};

/**
 * @brief Returns the lossy compression scheme given its name
 *
 * Calls error if the name is not found in the list of filters.
 *
 * @param name The name of the filter
 * @return The #lossy_compression_schemes
 */
enum lossy_compression_schemes compression_scheme_from_name(const char* name) {

  for (int i = 0; i < compression_level_count; ++i) {
    if (strcmp(name, lossy_compression_schemes_names[i]) == 0)
      return (enum lossy_compression_schemes)i;
  }

  error("Invalid lossy compression scheme name: '%s'", name);
  return (enum lossy_compression_schemes)0;
}

#ifdef HAVE_HDF5

/**
 * @brief Sets the properties and type of an HDF5 dataspace to apply a given
 * lossy compression scheme.
 *
 * The filter name is only updated if a filter was applied. If no lossy
 * compression is chosen the string is not altered.
 *
 * @param h_prop The properties of the dataspace.
 * @param h_type The type of the dataspace.
 * @param comp The lossy compression scheme to apply to this dataspace.
 * @param field_name The name of the field to write in the dataspace.
 * @param filter_name (return) The name of the filter (if one is applied).
 */
void set_hdf5_lossy_compression(hid_t* h_prop, hid_t* h_type,
                                const enum lossy_compression_schemes comp,
                                const char* field_name, char filter_name[32]) {

  if (comp == compression_do_not_write) {
    error(
        "Applying a lossy compression filter to a field labelled as 'do not "
        "write'");
  }

  else if (comp == compression_write_d_scale_6) {

    /* Scale filter with a scaling by 10^6 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 6);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_5) {

    /* Scale filter with a scaling by 10^5 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 5);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_4) {

    /* Scale filter with a scaling by 10^4 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 4);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_3) {

    /* Scale filter with a scaling by 10^3 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 3);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_2) {

    /* Scale filter with a scaling by 10^2 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 2);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_1) {

    /* Scale filter with a scaling by 10^1 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 1);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_mantissa_9) {

    /* Double numbers with 9-bits mantissa and 11-bits exponent
     *
     * This has a relative accuracy of log10(2^(9+1)) = 3.01 decimal digits
     * and the same range as a regular float.
     *
     * This leads to a compression ratio of 3.05 */

    /* Note a regular IEEE-754 double has:
     * - size = 8
     * - m_size = 52
     * - e_size = 11
     * i.e. 52 + 11 + 1 (the sign bit) == 64 bits (== 8 bytes) */

    const int size = 8;
    const int m_size = 9;
    const int e_size = 11;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_FLOAT);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  else if (comp == compression_write_d_mantissa_13) {

    /* Double numbers with 13-bits mantissa and 11-bits exponent
     *
     * This has a relative accuracy of log10(2^(13+1)) = 4.21 decimal digits
     * and the same range as a regular float.
     *
     * This leads to a compression ratio of 2.56 */

    /* Note a regular IEEE-754 double has:
     * - size = 8
     * - m_size = 52
     * - e_size = 11
     * i.e. 52 + 11 + 1 (the sign bit) == 64 bits (== 8 bytes) */

    const int size = 8;
    const int m_size = 13;
    const int e_size = 11;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_FLOAT);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  else if (comp == compression_write_f_mantissa_9) {

    /* Float numbers with 9-bits mantissa and 8-bits exponent
     *
     * This has a relative accuracy of log10(2^(9+1)) = 3.01 decimal digits
     * and the same range as a regular float.
     *
     * This leads to a compression ratio of 1.778 */

    /* Note a regular IEEE-754 float has:
     * - size = 4
     * - m_size = 23
     * - e_size = 8
     * i.e. 23 + 8 + 1 (the sign bit) == 32 bits (== 4 bytes) */

    const int size = 4;
    const int m_size = 9;
    const int e_size = 8;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_FLOAT);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  else if (comp == compression_write_f_mantissa_13) {

    /* Float numbers with 13-bits mantissa and 8-bits exponent
     *
     * This has a relative accuracy of log10(2^(13+1)) = 4.21 decimal digits
     * and the same range as a regular float.
     *
     * This leads to a compression ratio of 1.455 */

    /* Note a regular IEEE-754 float has:
     * - size = 4
     * - m_size = 23
     * - e_size = 8
     * i.e. 23 + 8 + 1 (the sign bit) == 32 bits (== 4 bytes) */

    const int size = 4;
    const int m_size = 13;
    const int e_size = 8;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_FLOAT);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  else if (comp == compression_write_half_float) {

    /* Float numbers with 10-bits mantissa and 5-bits exponent
     *
     * This has a relative accuracy of log10(2^(10+1)) = 3.31 decimal digits
     * and can represent numbers in the range [6.1*10^-5 - 6.5*10^4].
     *
     * This leads to a compression ratio of 2 */

    /* Note a regular IEEE-754 float has:
     * - size = 4
     * - m_size = 23
     * - e_size = 8
     * i.e. 23 + 8 + 1 (the sign bit) == 32 bits (== 4 bytes) */

    const int size = 4;
    const int m_size = 10;
    const int e_size = 5;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_FLOAT);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  else if (comp == compression_write_bfloat_16) {

    /* Float numbers with 7-bits mantissa and 8-bits exponent
     *
     * This has a relative accuracy of log10(2^(7+1)) = 2.4 decimal digits
     * and the same range as a regular float.
     *
     * This leads to a compression ratio of 2 */

    /* Note a regular IEEE-754 float has:
     * - size = 4
     * - m_size = 23
     * - e_size = 8
     * i.e. 23 + 8 + 1 (the sign bit) == 32 bits (== 4 bytes) */

    const int size = 4;
    const int m_size = 7;
    const int e_size = 8;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_FLOAT);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  else if (comp == compression_write_Nbit_36 ||
           comp == compression_write_Nbit_32 ||
           comp == compression_write_Nbit_40 ||
           comp == compression_write_Nbit_44 ||
           comp == compression_write_Nbit_48 ||
           comp == compression_write_Nbit_56) {

    int n_bits = 0;
    if (comp == compression_write_Nbit_32)
      n_bits = 32;
    else if (comp == compression_write_Nbit_36)
      n_bits = 36;
    else if (comp == compression_write_Nbit_40)
      n_bits = 40;
    else if (comp == compression_write_Nbit_44)
      n_bits = 44;
    else if (comp == compression_write_Nbit_48)
      n_bits = 48;
    else if (comp == compression_write_Nbit_56)
      n_bits = 56;
    else
      error("Invalid Nbit size");

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_NATIVE_LLONG);
    hid_t h_err = H5Tset_precision(*h_type, n_bits);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_offset(*h_type, 0);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  /* Other case: Do nothing! */

  /* Finish by returning the filter name */
  if (comp != compression_write_lossless)
    snprintf(filter_name, 32, "%s", lossy_compression_schemes_names[comp]);
}

#endif /* HAVE_HDF5 */
