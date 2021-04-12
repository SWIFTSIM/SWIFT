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
#ifndef SWIFT_IO_COMPRESSION_H
#define SWIFT_IO_COMPRESSION_H

/* Config parameters. */
#include "../config.h"

/**
 * @brief Compression levels for snapshot fields
 */
enum lossy_compression_schemes {
  compression_do_not_write = 0,    /*!< Do not write that field */
  compression_write_lossless,      /*!< Do not apply any lossy compression */
  compression_write_d_scale_1,     /*!< D-scale filter of magnitude 10^1 */
  compression_write_d_scale_2,     /*!< D-scale filter of magnitude 10^2 */
  compression_write_d_scale_3,     /*!< D-scale filter of magnitude 10^3 */
  compression_write_d_scale_6,     /*!< D-scale filter of magnitude 10^6 */
  compression_write_f_mantissa_9,  /*!< Conversion to 9-bits mantissa float */
  compression_write_f_mantissa_13, /*!< Conversion to 13-bits mantissa float */
  compression_write_half_float,    /*!< Conversion to IEEE754 half-float */
  compression_write_bfloat_16,     /*!< Conversion to Bfloat16 */
  compression_write_Nbit_36, /*!< Conversion to 36-bit int (from long long) */
  compression_write_Nbit_40, /*!< Conversion to 40-bit int (from long long) */
  compression_write_Nbit_44, /*!< Conversion to 44-bit int (from long long) */
  compression_write_Nbit_48, /*!< Conversion to 48-bit int (from long long) */
  compression_write_Nbit_56, /*!< Conversion to 56-bit int (from long long) */
  /* Counter, always leave last */
  compression_level_count,
};

/**
 * @brief Names of the compression levels, used in the select_output.yml
 *        parameter file.
 **/
extern const char* lossy_compression_schemes_names[];

enum lossy_compression_schemes compression_scheme_from_name(const char* name);

#ifdef HAVE_HDF5

#include <hdf5.h>

void set_hdf5_lossy_compression(hid_t* h_prop, hid_t* h_type,
                                const enum lossy_compression_schemes comp,
                                const char* field_name);

#endif /* HAVE_HDF5 */

#endif /* SWIFT_IO_COMPRESSION_H */
