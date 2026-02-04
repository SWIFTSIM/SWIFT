/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026   Jacob Kegerreis (j.kegerreis@imperial.ac.uk).
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
#ifndef SWIFT_MIXED_HHE_HEAVY_EQUATION_OF_STATE_H
#define SWIFT_MIXED_HHE_HEAVY_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/mixed_HHe_heavy.h
 *
 * Contains the mixed HHe--heavy-elements EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 * 
 * Many of these should be replaced by generalised interpolation functions...
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "units.h"
#include "utilities.h"

// mixedHHeHeavy parameters
struct single_mixedHHeHeavy_params {
  float *table_mix;
  float *table_log_rho;
  float *table_log_T;
  float *table_log_u_mix_rho_T;
  float *table_P_mix_rho_T;
  float *table_c_mix_rho_T;
  float *table_log_s_mix_rho_T;
  int version_date, num_mix, num_rho, num_T, num_rho_T;
  enum eos_planetary_material_id mat_id;
};
struct mixedHHeHeavy_params {
  struct *single_mixedHHeHeavy_params single_params;
  int num_mat;
};

// Parameter values for each material
INLINE static void set_mixedHHeHeavy(struct mixedHHeHeavy_params *mat,
                                     enum eos_planetary_material_id mat_id,
                                     int num_mat) {
  mat->mat_id = mat_id;
  mat->version_date = 20260203;
  mat->num_mat = num_mat;
  mat->single_params = (single_mixedHHeHeavy_params *)malloc(
                          mat->num_mat * sizeof(single_mixedHHeHeavy_params));
}
INLINE static void set_mixedHHeHeavy_rock(struct single_mixedHHeHeavy_params *mat,
                                          enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 20260203;
}
INLINE static void set_mixedHHeHeavy_water(struct single_mixedHHeHeavy_params *mat,
                                           enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 20260203;
}
INLINE static void set_mixedHHeHeavy_iron(struct single_mixedHHeHeavy_params *mat,
                                          enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 20260203;
}

/*
    Read the table data from the HDF5 file.
*/
INLINE static void load_table_mixedHHeHeavy(struct mixedHHeHeavy_params *mat,
                                     char *table_file) {
                                      
#ifndef HAVE_HDF5
  error("Need HDF5 to read mixed EoS tables");
#endif

  // Load table file
  hid_t file_id = H5Fopen(table_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("Failed to open mixedHHeHeavy EoS file %s\n", table_file);

  // Load properties
  hid_t grp_header = H5Gopen(table_file, "Header", H5P_DEFAULT);
  io_read_attribute(grp_header, "version_date", INT, &mat->version_date);
  io_read_attribute(grp_header, "num_mix", INT, &mat->num_mix);
  io_read_attribute(grp_header, "num_rho", INT, &mat->num_rho);
  io_read_attribute(grp_header, "num_T", INT, &mat->num_T);
  mat->num_rho_T = mat->num_rho * mat->num_T;
  const int num_mix_rho_T = mat->num_mix * mat->num_rho_T;

  // Allocate table memory
  mat->table_mix = (float *)malloc(mat->num_mix * sizeof(float));
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_T = (float *)malloc(mat->num_T * sizeof(float));
  mat->table_log_u_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));
  mat->table_P_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));
  mat->table_c_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));
  mat->table_log_s_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));

  // Load table data
  hid_t grp_table = H5Gopen(table_file, "/table", H5P_DEFAULT);
  io_read_array_dataset(grp_table, "table_mix", FLOAT, 
                        &mat->table_mix, mat->num_mix);
  io_read_array_dataset(grp_table, "table_log_rho", FLOAT, 
                        &mat->table_log_rho, mat->num_rho);
  io_read_array_dataset(grp_table, "table_log_T", FLOAT, 
                        &mat->table_log_T, mat->num_T);
  io_read_array_dataset(grp_table, "table_log_u_mix_rho_T", FLOAT, 
                        &mat->table_log_u_mix_rho_T, num_mix_rho_T);
  io_read_array_dataset(grp_table, "table_P_mix_rho_T", FLOAT, 
                        &mat->table_P_mix_rho_T, num_mix_rho_T);
  io_read_array_dataset(grp_table, "table_c_mix_rho_T", FLOAT, 
                        &mat->table_c_mix_rho_T, num_mix_rho_T);
  io_read_array_dataset(grp_table, "table_log_s_mix_rho_T", FLOAT, 
                        &mat->table_log_s_mix_rho_T, num_mix_rho_T);
}

// Misc. modifications
INLINE static void prepare_table_SESAME(struct SESAME_params *mat) {

  // Convert to log
  // Density and temperature
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    mat->table_log_rho[i_rho] = logf(mat->table_log_rho[i_rho]);
  }  
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    mat->table_log_T[i_T] = logf(mat->table_log_T[i_T]);
  }
  // Sp. int. energy and entropy
  for (int i_mix = 0; i_mix < mat->num_mix; i_mix++) {
    for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
      for (int i_T = 0; i_T < mat->num_T; i_T++) {
        mat->table_log_u_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] =
            logf(mat->table_log_u_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T]);

        mat->table_log_s_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] =
            logf(mat->table_log_s_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T]);
      }
    }
  }
}

// Convert to internal units
INLINE static void convert_units_mixedHHeHeavy(struct mixedHHeHeavy_params *mat,
                                        const struct unit_system *us) {

  // Convert input table values from all-SI to internal units
  struct unit_system si;
  units_init_si(&si);

  // Densities (log)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    mat->table_log_rho[i_rho] +=
        logf(units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY) /
             units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  }

  // Temperatures (log)
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    mat->table_log_T[i_T] +=
        logf(units_cgs_conversion_factor(&si, UNIT_CONV_TEMPERATURE) /
             units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
  }

  // Sp. int. energies (log), pressures, sound speeds, and sp. entropies
  for (int i_mix = 0; i_mix < mat->num_mix; i_mix++) {
    for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
      for (int i_T = 0; i_T < mat->num_T; i_T++) {
        mat->table_log_u_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] += 
            logf(units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
                 units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
        mat->table_P_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] *=
            units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
            units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
        mat->table_c_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] *=
            units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
            units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
        mat->table_log_s_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] +=
            logf(units_cgs_conversion_factor(
                    &si, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS) /
                units_cgs_conversion_factor(
                    us, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS));
      }
    }
  }
}

// gas_internal_energy_from_entropy
INLINE static float mixedHHeHeavy_internal_energy_from_entropy(
    const float density, const float entropy, const struct mixedHHeHeavy_params *mat) {

  // Return zero if entropy is zero
  if (entropy <= 0.f) {
    return 0.f;
  }

  const float log_rho = logf(density);
  const float log_s = logf(entropy);

  // 2D interpolation (bilinear with log(rho), log(s) to find u(rho, s))

  // Density index
  int idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }

  // Sp. entropy at this and the next density (in relevant slice of s array)
  int idx_s_1 = find_value_in_monot_incr_array(
      log_s, mat->table_log_s_rho_T + idx_rho * mat->num_T, mat->num_T);
  int idx_s_2 = find_value_in_monot_incr_array(
      log_s, mat->table_log_s_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_s_1 <= -1) {
    idx_s_1 = 0;
  } else if (idx_s_1 >= mat->num_T) {
    idx_s_1 = mat->num_T - 2;
  }
  if (idx_s_2 <= -1) {
    idx_s_2 = 0;
  } else if (idx_s_2 >= mat->num_T) {
    idx_s_2 = mat->num_T - 2;
  }

  // Check for duplicates in mixedHHeHeavy tables before interpolation
  float intp_rho, intp_s_1, intp_s_2;
  if (mat->table_log_rho[idx_rho + 1] != mat->table_log_rho[idx_rho]) {
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
               (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.f;
  }
  if (mat->table_log_s_rho_T[idx_rho * mat->num_T + (idx_s_1 + 1)] !=
      mat->table_log_s_rho_T[idx_rho * mat->num_T + idx_s_1]) {
    intp_s_1 =
        (log_s - mat->table_log_s_rho_T[idx_rho * mat->num_T + idx_s_1]) /
        (mat->table_log_s_rho_T[idx_rho * mat->num_T + (idx_s_1 + 1)] -
         mat->table_log_s_rho_T[idx_rho * mat->num_T + idx_s_1]);
  } else {
    intp_s_1 = 1.f;
  }
  if (mat->table_log_s_rho_T[(idx_rho + 1) * mat->num_T + (idx_s_2 + 1)] !=
      mat->table_log_s_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2]) {
    intp_s_2 =
        (log_s - mat->table_log_s_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2]) /
        (mat->table_log_s_rho_T[(idx_rho + 1) * mat->num_T + (idx_s_2 + 1)] -
         mat->table_log_s_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2]);
  } else {
    intp_s_2 = 1.f;
  }

  // Table values
  const float log_u_1 = mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_s_1];
  const float log_u_2 =
      mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_s_1 + 1];
  const float log_u_3 =
      mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2];
  const float log_u_4 =
      mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2 + 1];

  // If below the minimum s at this rho then just use the lowest table values
  if ((idx_rho > 0.f) && ((intp_s_1 < 0.f) || (intp_s_2 < 0.f) ||
                          (log_u_1 > log_u_2) || (log_u_3 > log_u_4))) {
    intp_s_1 = 0;
    intp_s_2 = 0;
  }

  // Interpolate with the log values
  const float u =
      (1.f - intp_rho) * ((1.f - intp_s_1) * log_u_1 + intp_s_1 * log_u_2) +
      intp_rho * ((1.f - intp_s_2) * log_u_3 + intp_s_2 * log_u_4);

  // Convert back from log
  return expf(u);
}

// gas_pressure_from_entropy
INLINE static float mixedHHeHeavy_pressure_from_entropy(
    const float density, const float entropy, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float mixedHHeHeavy_entropy_from_pressure(
    const float density, const float pressure,
    const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float mixedHHeHeavy_soundspeed_from_entropy(
    const float density, const float entropy, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float mixedHHeHeavy_entropy_from_internal_energy(
    const float density, const float u, const struct mixedHHeHeavy_params *mat) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float single_mixedHHeHeavy_pressure_from_internal_energy(
    const float density, const float u, const float mix, 
    const struct single_mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 3D interpolation (linear with mix, log(rho), log(u)) to find P(mix, rho, u)

  // Mix
  int idx_mix =
      find_value_in_monot_incr_array(mix, mat->table_mix, mat->num_mix);
  float intp_mix = (mix - mat->table_mix[idx_mix]) / 
                   (mat->table_mix[idx_mix + 1] - mat->table_mix[idx_mix]);

  // Density index
  int idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }

  // 2D interpolate within first single-mix subtable

  // Sp. int. energy at this and the next density (in relevant slice of each u array)
  int idx_u_1 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_mix_rho_T + idx_mix * mat->num_rho_T
            + idx_rho * mat->num_T, mat->num_T);
  int idx_u_2 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_mix_rho_T + idx_mix * mat->num_rho_T 
            + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_u_1 <= -1) {
    idx_u_1 = 0;
  } else if (idx_u_1 >= mat->num_T) {
    idx_u_1 = mat->num_T - 2;
  }
  if (idx_u_2 <= -1) {
    idx_u_2 = 0;
  } else if (idx_u_2 >= mat->num_T) {
    idx_u_2 = mat->num_T - 2;
  }

  // Check for duplicates in the table before interpolation
  float intp_rho, intp_u_1, intp_u_2;
  int idx_u_1_mix_rho_T = idx_mix * mat->num_rho_T + idx_rho * mat->num_T + idx_u_1;
  int idx_u_2_mix_rho_T = idx_mix * mat->num_rho_T + (idx_rho + 1) * mat->num_T + idx_u_2;
  intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
              (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  intp_u_1 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]);
  intp_u_2 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]);

  // Table values
  float P_1 = mat->table_P_mix_rho_T[idx_u_1_mix_rho_T];
  float P_2 = mat->table_P_mix_rho_T[idx_u_1_mix_rho_T + 1];
  float P_3 = mat->table_P_mix_rho_T[idx_u_2_mix_rho_T];
  float P_4 = mat->table_P_mix_rho_T[idx_u_2_mix_rho_T + 1];

  // If below the minimum u at this rho then just use the lowest table values
  if ((idx_rho > 0.f) &&
      ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (P_1 > P_2) || (P_3 > P_4))) {
    intp_u_1 = 0;
    intp_u_2 = 0;
  }

  // Interpolate with the log values
  P_1 = logf(P_1);
  P_2 = logf(P_2);
  P_3 = logf(P_3);
  P_4 = logf(P_4);
  float P_1_2 = (1.f - intp_u_1) * P_1 + intp_u_1 * P_2;
  float P_3_4 = (1.f - intp_u_2) * P_3 + intp_u_2 * P_4;

  float P = (1.f - intp_rho) * P_1_2 + intp_rho * P_3_4;

  // 2D interpolate within second single-mix table, then combine
  // (todo: tidy into a generalised function!)
  if ((mix > 0.f) && (mix < 1.f)) {
    // Sp. int. energy at this and the next density (in relevant slice of each u array)
    int idx_u_1 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_mix_rho_T + (idx_mix + 1) * mat->num_rho_T
              + idx_rho * mat->num_T, mat->num_T);
    int idx_u_2 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_mix_rho_T + (idx_mix + 1) * mat->num_rho_T 
              + (idx_rho + 1) * mat->num_T, mat->num_T);

    // If outside the table then extrapolate from the edge and edge-but-one values
    if (idx_u_1 <= -1) {
      idx_u_1 = 0;
    } else if (idx_u_1 >= mat->num_T) {
      idx_u_1 = mat->num_T - 2;
    }
    if (idx_u_2 <= -1) {
      idx_u_2 = 0;
    } else if (idx_u_2 >= mat->num_T) {
      idx_u_2 = mat->num_T - 2;
    }

    // Check for duplicates in the table before interpolation
    idx_u_1_mix_rho_T = (idx_mix + 1) * mat->num_rho_T + idx_rho * mat->num_T + idx_u_1;
    idx_u_2_mix_rho_T = (idx_mix + 1) * mat->num_rho_T + (idx_rho + 1) * mat->num_T + idx_u_2;
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
                (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
    intp_u_1 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]);
    intp_u_2 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]);

    // Table values
    P_1 = mat->table_P_mix_rho_T[idx_u_1_mix_rho_T];
    P_2 = mat->table_P_mix_rho_T[idx_u_1_mix_rho_T + 1];
    P_3 = mat->table_P_mix_rho_T[idx_u_2_mix_rho_T];
    P_4 = mat->table_P_mix_rho_T[idx_u_2_mix_rho_T + 1];

    // If below the minimum u at this rho then just use the lowest table values
    if ((idx_rho > 0.f) &&
        ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (P_1 > P_2) || (P_3 > P_4))) {
      intp_u_1 = 0;
      intp_u_2 = 0;
    }

    // Interpolate with the log values
    P_1 = logf(P_1);
    P_2 = logf(P_2);
    P_3 = logf(P_3);
    P_4 = logf(P_4);
    P_1_2 = (1.f - intp_u_1) * P_1 + intp_u_1 * P_2;
    P_3_4 = (1.f - intp_u_2) * P_3 + intp_u_2 * P_4;

    const float P_mix2 = (1.f - intp_rho) * P_1_2 + intp_rho * P_3_4;

    // Interpolate between mix values
    P = (1 - intp_mix) * P + intp_mix * P_mix2;
  }

  // Convert back from log
  return expf(P);
}

INLINE static float mixedHHeHeavy_pressure_from_internal_energy(
    const float density, const float u, const float* mixes, 
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  // No heavy-element fraction
  float mix_tot = 0.f;
  for (i_mat = 0; i_mat < num_mat; i_mat++) {
    mix_tot += mixes[i_mat];
  }
  if (mix_tot == 0.f) {
    return single_mixedHHeHeavy_pressure_from_internal_energy(
      density, u, 0.f, mat->single_params[0])
  }

  // Accumulate contribution from each non-zero mix
  float P = 0.f;
  for (i_mat = 0; i_mat < num_mat; i_mat++) {
    float mix = mixes[i_mat];
    if (mix > 0.f) {
      // Evaluate for this single heavy mix
      float P_mat = single_mixedHHeHeavy_pressure_from_internal_energy(
        density, u, mix, mat->single_params[i_mat])

      P += P_mat * mix / mix_tot;
    }
  }

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float mixedHHeHeavy_internal_energy_from_pressure(
    const float density, const float P, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float single_mixedHHeHeavy_soundspeed_from_internal_energy(
    const float density, const float u, const float mix, 
    const struct single_mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 3D interpolation (linear with mix, log(rho), log(u)) to find c(mix, rho, u)

  // Mix
  int idx_mix =
      find_value_in_monot_incr_array(mix, mat->table_mix, mat->num_mix);
  float intp_mix = (mix - mat->table_mix[idx_mix]) / 
                    (mat->table_mix[idx_mix + 1] - mat->table_mix[idx_mix]);

  // Density index
  int idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }

  // 2D interpolate within first single-mix subtable

  // Sp. int. energy at this and the next density (in relevant slice of each u array)
  int idx_u_1 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_mix_rho_T + idx_mix * mat->num_rho_T
            + idx_rho * mat->num_T, mat->num_T);
  int idx_u_2 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_mix_rho_T + idx_mix * mat->num_rho_T 
            + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_u_1 <= -1) {
    idx_u_1 = 0;
  } else if (idx_u_1 >= mat->num_T) {
    idx_u_1 = mat->num_T - 2;
  }
  if (idx_u_2 <= -1) {
    idx_u_2 = 0;
  } else if (idx_u_2 >= mat->num_T) {
    idx_u_2 = mat->num_T - 2;
  }

  // Check for duplicates in the table before interpolation
  float intp_rho, intp_u_1, intp_u_2;
  int idx_u_1_mix_rho_T = idx_mix * mat->num_rho_T + idx_rho * mat->num_T + idx_u_1;
  int idx_u_2_mix_rho_T = idx_mix * mat->num_rho_T + (idx_rho + 1) * mat->num_T + idx_u_2;
  intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
              (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  intp_u_1 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]);
  intp_u_2 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]);

  // Table values
  float c_1 = mat->table_c_mix_rho_T[idx_u_1_mix_rho_T];
  float c_2 = mat->table_c_mix_rho_T[idx_u_1_mix_rho_T + 1];
  float c_3 = mat->table_c_mix_rho_T[idx_u_2_mix_rho_T];
  float c_4 = mat->table_c_mix_rho_T[idx_u_2_mix_rho_T + 1];

  // If below the minimum u at this rho then just use the lowest table values
  if ((idx_rho > 0.f) &&
      ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (c_1 > c_2) || (c_3 > c_4))) {
    intp_u_1 = 0;
    intp_u_2 = 0;
  }

  // Interpolate with the log values
  c_1 = logf(c_1);
  c_2 = logf(c_2);
  c_3 = logf(c_3);
  c_4 = logf(c_4);
  float c_1_2 = (1.f - intp_u_1) * c_1 + intp_u_1 * c_2;
  float c_3_4 = (1.f - intp_u_2) * c_3 + intp_u_2 * c_4;

  float c = (1.f - intp_rho) * c_1_2 + intp_rho * c_3_4;

  // 2D interpolate within second single-mix table, then combine
  // (todo: tidy into a generalised function!)
  if ((mix > 0.f) && (mix < 1.f)) {
    // Sp. int. energy at this and the next density (in relevant slice of each u array)
    int idx_u_1 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_mix_rho_T + (idx_mix + 1) * mat->num_rho_T
              + idx_rho * mat->num_T, mat->num_T);
    int idx_u_2 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_mix_rho_T + (idx_mix + 1) * mat->num_rho_T 
              + (idx_rho + 1) * mat->num_T, mat->num_T);

    // If outside the table then extrapolate from the edge and edge-but-one values
    if (idx_u_1 <= -1) {
      idx_u_1 = 0;
    } else if (idx_u_1 >= mat->num_T) {
      idx_u_1 = mat->num_T - 2;
    }
    if (idx_u_2 <= -1) {
      idx_u_2 = 0;
    } else if (idx_u_2 >= mat->num_T) {
      idx_u_2 = mat->num_T - 2;
    }

    // Check for duplicates in the table before interpolation
    idx_u_1_mix_rho_T = (idx_mix + 1) * mat->num_rho_T + idx_rho * mat->num_T + idx_u_1;
    idx_u_2_mix_rho_T = (idx_mix + 1) * mat->num_rho_T + (idx_rho + 1) * mat->num_T + idx_u_2;
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
                (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
    intp_u_1 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T1[idx_u_1_mix_rho_T]);
    intp_u_2 = (log_u - mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T1[idx_u_2_mix_rho_T]);

    // Table values
    c_1 = mat->table_c_mix_rho_T[idx_u_1_mix_rho_T];
    c_2 = mat->table_c_mix_rho_T[idx_u_1_mix_rho_T + 1];
    c_3 = mat->table_c_mix_rho_T[idx_u_2_mix_rho_T];
    c_4 = mat->table_c_mix_rho_T[idx_u_2_mix_rho_T + 1];

    // If below the minimum u at this rho then just use the lowest table values
    if ((idx_rho > 0.f) &&
        ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (c_1 > c_2) || (c_3 > c_4))) {
      intp_u_1 = 0;
      intp_u_2 = 0;
    }

    // Interpolate with the log values
    c_1 = logf(c_1);
    c_2 = logf(c_2);
    c_3 = logf(c_3);
    c_4 = logf(c_4);
    c_1_2 = (1.f - intp_u_1) * c_1 + intp_u_1 * c_2;
    c_3_4 = (1.f - intp_u_2) * c_3 + intp_u_2 * c_4;

    const float c_mix2 = (1.f - intp_rho) * c_1_2 + intp_rho * c_3_4;

    // Interpolate between mix values
    c = (1 - intp_mix) * c + intp_mix * c_mix2;
  }

  // Convert back from log
  return expf(c);
}
  
INLINE static float mixedHHeHeavy_soundspeed_from_internal_energy(
    const float density, const float u, const float* mixes, 
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  // No heavy-element fraction
  float mix_tot = 0.f;
  for (i_mat = 0; i_mat < num_mat; i_mat++) {
    mix_tot += mixes[i_mat];
  }
  if (mix_tot == 0.f) {
    return single_mixedHHeHeavy_pressure_from_internal_energy(
      density, u, 0.f, mat->single_params[0])
  }

  // Accumulate contribution from each non-zero mix
  float c = 0.f;
  for (i_mat = 0; i_mat < num_mat; i_mat++) {
    float mix = mixes[i_mat];
    if (mix > 0.f) {
      // Evaluate for this single heavy mix
      float c_mat = single_mixedHHeHeavy_pressure_from_internal_energy(
        density, u, mix, mat->single_params[i_mat])

      c += c_mat * mix / mix_tot;
    }
  }

  return c;
}

// gas_soundspeed_from_pressure
INLINE static float mixedHHeHeavy_soundspeed_from_pressure(
    const float density, const float P, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_temperature_from_internal_energy
INLINE static float mixedHHeHeavy_temperature_from_internal_energy(
    const float density, const float u, const struct mixedHHeHeavy_params *mat) {

  // Return zero if zero internal energy
  if (u <= 0.f) {
    return 0.f;
  }

  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u) to find T(rho, u)

  // Density index
  int idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  int idx_u_1 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_rho_T + idx_rho * mat->num_T, mat->num_T);
  int idx_u_2 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_u_1 <= -1) {
    idx_u_1 = 0;
  } else if (idx_u_1 >= mat->num_T) {
    idx_u_1 = mat->num_T - 2;
  }
  if (idx_u_2 <= -1) {
    idx_u_2 = 0;
  } else if (idx_u_2 >= mat->num_T) {
    idx_u_2 = mat->num_T - 2;
  }

  // Check for duplicates in mixedHHeHeavy tables before interpolation
  float intp_u_1, intp_u_2;
  const int idx_u_1_mix_rho_T = idx_rho * mat->num_T + idx_u_1;
  if (mat->table_log_u_rho_T[idx_u_1_mix_rho_T + 1] !=
      mat->table_log_u_rho_T[idx_u_1_mix_rho_T]) {
    intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_1_mix_rho_T]) /
               (mat->table_log_u_rho_T[idx_u_1_mix_rho_T + 1] -
                mat->table_log_u_rho_T[idx_u_1_mix_rho_T]);
  } else {
    intp_u_1 = 1.f;
  }
  if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
      mat->table_log_u_rho_T[idx_u_2_mix_rho_T]) {
    intp_u_2 =
        (log_u - mat->table_log_u_rho_T[idx_u_2_mix_rho_T]) /
        (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
         mat->table_log_u_rho_T[idx_u_2_mix_rho_T]);
  } else {
    intp_u_2 = 1.f;
  }

  // Compute line points
  float log_rho_1 = mat->table_log_rho[idx_rho];
  float log_rho_2 = mat->table_log_rho[idx_rho + 1];
  float log_T_1 = mat->table_log_T[idx_u_1];
  log_T_1 +=
      intp_u_1 * (mat->table_log_T[idx_u_1 + 1] - mat->table_log_T[idx_u_1]);
  float log_T_2 = mat->table_log_T[idx_u_2];
  log_T_2 +=
      intp_u_2 * (mat->table_log_T[idx_u_2 + 1] - mat->table_log_T[idx_u_2]);

  // Intersect line passing through (log_rho_1, log_T_1), (log_rho_2, log_T_2)
  // with line density = log_rho

  // Check for log_T_1 == log_T_2
  float log_T;
  if (log_T_1 == log_T_2) {
    log_T = log_T_1;
  } else {
    // log_rho = slope*log_T + intercept
    const float slope = (log_rho_1 - log_rho_2) / (log_T_1 - log_T_2);
    const float intercept = log_rho_1 - slope * log_T_1;
    log_T = (log_rho - intercept) / slope;
  }

  // Convert back from log
  return expf(log_T);
}

// gas_density_from_pressure_and_temperature
INLINE static float mixedHHeHeavy_density_from_pressure_and_temperature(
    float P, float T, const struct mixedHHeHeavy_params *mat) {

  // Return zero if pressure is non-positive
  if (P <= 0.f) {
    return 0.f;
  }

  const float log_T = logf(T);

  // 2D interpolation (bilinear with log(T), P to find rho(T, P))

  // Temperature index
  int idx_T =
      find_value_in_monot_incr_array(log_T, mat->table_log_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_T <= -1) {
    idx_T = 0;
  } else if (idx_T >= mat->num_T) {
    idx_T = mat->num_T - 2;
  }

  // Pressure at this and the next temperature (in relevant vertical slice of P
  // array)
  int idx_P_1 = vertical_find_value_in_monot_incr_array(
      P, mat->table_P_rho_T, mat->num_rho, mat->num_T, idx_T);
  int idx_P_2 = vertical_find_value_in_monot_incr_array(
      P, mat->table_P_rho_T, mat->num_rho, mat->num_T, idx_T + 1);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_P_1 <= -1) {
    idx_P_1 = 0;
  } else if (idx_P_1 >= mat->num_rho) {
    idx_P_1 = mat->num_rho - 2;
  }
  if (idx_P_2 <= -1) {
    idx_P_2 = 0;
  } else if (idx_P_2 >= mat->num_rho) {
    idx_P_2 = mat->num_rho - 2;
  }

  // Check for duplicates in mixedHHeHeavy tables before interpolation
  float intp_P_1, intp_P_2;
  if (mat->table_P_rho_T[(idx_P_1 + 1) * mat->num_T + idx_T] !=
      mat->table_P_rho_T[idx_P_1 * mat->num_T + idx_T]) {
    intp_P_1 = (P - mat->table_P_rho_T[idx_P_1 * mat->num_T + idx_T]) /
               (mat->table_P_rho_T[(idx_P_1 + 1) * mat->num_T + idx_T] -
                mat->table_P_rho_T[idx_P_1 * mat->num_T + idx_T]);
  } else {
    intp_P_1 = 1.f;
  }
  if (mat->table_P_rho_T[(idx_P_2 + 1) * mat->num_T + (idx_T + 1)] !=
      mat->table_P_rho_T[idx_P_2 * mat->num_T + (idx_T + 1)]) {
    intp_P_2 = (P - mat->table_P_rho_T[idx_P_2 * mat->num_T + (idx_T + 1)]) /
               (mat->table_P_rho_T[(idx_P_2 + 1) * mat->num_T + (idx_T + 1)] -
                mat->table_P_rho_T[idx_P_2 * mat->num_T + (idx_T + 1)]);
  } else {
    intp_P_2 = 1.f;
  }

  // Compute line points
  const float log_T_1 = mat->table_log_T[idx_T];
  const float log_T_2 = mat->table_log_T[idx_T + 1];
  const float log_rho_1 = (1.f - intp_P_1) * mat->table_log_rho[idx_P_1] +
                          intp_P_1 * mat->table_log_rho[idx_P_1 + 1];
  const float log_rho_2 = (1.f - intp_P_2) * mat->table_log_rho[idx_P_2] +
                          intp_P_2 * mat->table_log_rho[idx_P_2 + 1];

  // Intersect line passing through (log_rho_1, log_T_1), (log_rho_2, log_T_2)
  // with line temperature = log_T

  // Check for log_rho_1 == log_rho_2
  float log_rho;
  if (log_rho_1 == log_rho_2) {
    log_rho = log_rho_1;
  } else {
    // log_T = slope*log_rho + intercept
    const float slope = (log_T_1 - log_T_2) / (log_rho_1 - log_rho_2);
    const float intercept = log_T_1 - slope * log_rho_1;
    log_rho = (log_T - intercept) / slope;
  }

  // Convert back from log
  return expf(log_rho);
}

// gas_density_from_pressure_and_internal_energy
INLINE static float mixedHHeHeavy_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    const struct mixedHHeHeavy_params *mat) {

  // Return the unchanged density if u or P is non-positive
  if (u <= 0.f || P <= 0.f) {
    return rho_sph;
  }

  // Convert inputs to log
  const float log_u = logf(u);
  const float log_P = logf(P);
  const float log_rho_ref = logf(rho_ref);

  // Find rounded down index of reference density. This is where we start our
  // search
  int idx_rho_ref = find_value_in_monot_incr_array(
      log_rho_ref, mat->table_log_rho, mat->num_rho);

  // If no roots are found in the current search range, we increase search range
  // by search_factor_log_rho above and below the reference density each
  // iteration.
  const float search_factor_log_rho = logf(10.f);

  // Initialise the minimum and maximum densities we're searching to at the
  // reference density. These will change before the first iteration.
  float log_rho_min = log_rho_ref;
  float log_rho_max = log_rho_ref;

  // Initialise search indices around rho_ref
  int idx_rho_below_min, idx_rho_above_max;

  // If we find a root, it will get stored as closest_root
  float closest_root = -1.f;
  float root_below;

  // Initialise pressures
  float P_above_lower, P_above_upper;
  float P_below_lower, P_below_upper;
  P_above_upper = 0.f;
  P_below_lower = 0.f;

  // Increase search range by search_factor_log_rho
  log_rho_max += search_factor_log_rho;
  idx_rho_above_max = find_value_in_monot_incr_array(
      log_rho_max, mat->table_log_rho, mat->num_rho);
  log_rho_min -= search_factor_log_rho;
  idx_rho_below_min = find_value_in_monot_incr_array(
      log_rho_min, mat->table_log_rho, mat->num_rho);

  // Ensure not outside the table
  if (idx_rho_above_max >= mat->num_rho) {
    idx_rho_above_max = mat->num_rho - 1;
  }
  if (idx_rho_below_min <= -1) {
    idx_rho_below_min = 0;
  }

  // When searching above/below, we are looking for where the pressure P(rho, u)
  // of the table densities changes from being less than to more than, or vice
  // versa, the desired pressure. If this is the case, there is a root between
  // these table values of rho.

  // First look for roots above rho_ref
  int idx_u_1_mix_rho_T;
  float intp_rho;
  for (int idx_rho = idx_rho_ref; idx_rho <= idx_rho_above_max; idx_rho++) {

    // This is similar to P_u_rho, but we're not interested in intp_rho,
    // and instead calculate the pressure for both intp_rho=0 and intp_rho=1

    // Sp. int. energy at this and the next density (in relevant slice of u
    // array)
    int idx_u_1 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_rho_T + idx_rho * mat->num_T, mat->num_T);
    int idx_u_2 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

    // If outside the table then extrapolate from the edge and edge-but-one
    // values
    if (idx_u_1 <= -1) {
      idx_u_1 = 0;
    } else if (idx_u_1 >= mat->num_T) {
      idx_u_1 = mat->num_T - 2;
    }
    if (idx_u_2 <= -1) {
      idx_u_2 = 0;
    } else if (idx_u_2 >= mat->num_T) {
      idx_u_2 = mat->num_T - 2;
    }

    float intp_u_1, intp_u_2;
    idx_u_1_mix_rho_T = idx_rho * mat->num_T + idx_u_1;
    if (mat->table_log_u_rho_T[idx_u_1_mix_rho_T + 1] !=
        mat->table_log_u_rho_T[idx_u_1_mix_rho_T]) {
      intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_1_mix_rho_T]) /
                 (mat->table_log_u_rho_T[idx_u_1_mix_rho_T + 1] -
                  mat->table_log_u_rho_T[idx_u_1_mix_rho_T]);
    } else {
      intp_u_1 = 1.f;
    }
    if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
        mat->table_log_u_rho_T[idx_u_2_mix_rho_T]) {
      intp_u_2 =
          (log_u -
           mat->table_log_u_rho_T[idx_u_2_mix_rho_T]) /
          (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
           mat->table_log_u_rho_T[idx_u_2_mix_rho_T]);
    } else {
      intp_u_2 = 1.f;
    }

    // Table values
    float P_1 = mat->table_P_rho_T[idx_u_1_mix_rho_T];
    float P_2 = mat->table_P_rho_T[idx_u_1_mix_rho_T + 1];
    float P_3 = mat->table_P_rho_T[idx_u_2_mix_rho_T];
    float P_4 = mat->table_P_rho_T[idx_u_2_mix_rho_T + 1];

    // If below the minimum u at this rho then just use the lowest table values
    if ((idx_rho > 0.f) &&
        ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (P_1 > P_2) || (P_3 > P_4))) {
      intp_u_1 = 0;
      intp_u_2 = 0;
    }

    // Interpolate with the log values
    P_1 = logf(P_1);
    P_2 = logf(P_2);
    P_3 = logf(P_3);
    P_4 = logf(P_4);
    float P_1_2 = (1.f - intp_u_1) * P_1 + intp_u_1 * P_2;
    float P_3_4 = (1.f - intp_u_2) * P_3 + intp_u_2 * P_4;

    // Pressure for intp_rho = 0
    P_above_lower = expf(P_1_2);

    // Because of linear interpolation, pressures are not exactly continuous
    // as we go from one side of a grid point to another. See if there is
    // a root between the last P_above_upper and the new P_above_lower,
    // which are approx the same.
    if (idx_rho != idx_rho_ref) {
      if ((P_above_lower - P) * (P_above_upper - P) <= 0) {
        closest_root = expf(mat->table_log_rho[idx_rho]);
        break;
      }
    }

    // Pressure for intp_rho = 1
    P_above_upper = expf(P_3_4);

    // Does the pressure of the adjacent table densities switch from being
    // above to below the desired pressure, or vice versa? If so, there is a
    // root.
    if ((P_above_lower - P) * (P_above_upper - P) <= 0.f) {

      // If there is a root, interpolate between the table values:
      intp_rho = (log_P - P_1_2) / (P_3_4 - P_1_2);

      closest_root = expf(mat->table_log_rho[idx_rho] +
                          intp_rho * (mat->table_log_rho[idx_rho + 1] -
                                      mat->table_log_rho[idx_rho]));

      // If the root is between the same table values as the reference value,
      // then this is the closest root, so we can return it without further
      // searching
      if (idx_rho == idx_rho_ref) {
        return closest_root;
      }

      // Found a root, so no need to search higher densities
      break;
    }
  }

  // If we found a root above, change search range below so that we're only
  // looking for closer (in log) roots than the one we found
  if (closest_root > 0.f) {
    log_rho_min = log_rho_ref - (logf(closest_root) - log_rho_ref);
    idx_rho_below_min = find_value_in_monot_incr_array(
        log_rho_min, mat->table_log_rho, mat->num_rho);
  }

  // Now look for roots below rho_ref
  for (int idx_rho = idx_rho_ref; idx_rho >= idx_rho_below_min; idx_rho--) {

    // Sp. int. energy at this and the next density (in relevant slice of u
    // array)
    int idx_u_1 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_rho_T + idx_rho * mat->num_T, mat->num_T);
    int idx_u_2 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

    // If outside the table then extrapolate from the edge and edge-but-one
    // values
    if (idx_u_1 <= -1) {
      idx_u_1 = 0;
    } else if (idx_u_1 >= mat->num_T) {
      idx_u_1 = mat->num_T - 2;
    }
    if (idx_u_2 <= -1) {
      idx_u_2 = 0;
    } else if (idx_u_2 >= mat->num_T) {
      idx_u_2 = mat->num_T - 2;
    }

    idx_u_1_mix_rho_T = idx_rho * mat->num_T + idx_u_1;
    float intp_u_1, intp_u_2;
    if (mat->table_log_u_rho_T[idx_u_1_mix_rho_T + 1] !=
        mat->table_log_u_rho_T[idx_u_1_mix_rho_T]) {
      intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_1_mix_rho_T]) /
                 (mat->table_log_u_rho_T[idx_u_1_mix_rho_T + 1] -
                  mat->table_log_u_rho_T[idx_u_1_mix_rho_T]);
    } else {
      intp_u_1 = 1.f;
    }
    if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
        mat->table_log_u_rho_T[idx_u_2_mix_rho_T]) {
      intp_u_2 =
          (log_u -
           mat->table_log_u_rho_T[idx_u_2_mix_rho_T]) /
          (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
           mat->table_log_u_rho_T[idx_u_2_mix_rho_T]);
    } else {
      intp_u_2 = 1.f;
    }

    // Table values
    float P_1 = mat->table_P_rho_T[idx_u_1_mix_rho_T];
    float P_2 = mat->table_P_rho_T[idx_u_1_mix_rho_T + 1];
    float P_3 = mat->table_P_rho_T[idx_u_2_mix_rho_T];
    float P_4 = mat->table_P_rho_T[idx_u_2_mix_rho_T + 1];

    // If below the minimum u at this rho then just use the lowest table values
    if ((idx_rho > 0.f) &&
        ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (P_1 > P_2) || (P_3 > P_4))) {
      intp_u_1 = 0;
      intp_u_2 = 0;
    }

    // Interpolate with the log values
    P_1 = logf(P_1);
    P_2 = logf(P_2);
    P_3 = logf(P_3);
    P_4 = logf(P_4);
    const float P_1_2 = (1.f - intp_u_1) * P_1 + intp_u_1 * P_2;
    const float P_3_4 = (1.f - intp_u_2) * P_3 + intp_u_2 * P_4;

    // Pressure for intp_rho = 1
    P_below_upper = expf(P_3_4);
    // Because of linear interpolation, pressures are not exactly continuous
    // as we go from one side of a grid point to another. See if there is
    // a root between the last P_below_lower and the new P_below_upper,
    // which are approx the same.
    if (idx_rho != idx_rho_ref) {
      if ((P_below_lower - P) * (P_below_upper - P) <= 0) {
        closest_root = expf(mat->table_log_rho[idx_rho + 1]);
        break;
      }
    }
    // Pressure for intp_rho = 0
    P_below_lower = expf(P_1_2);

    // Does the pressure of the adjacent table densities switch from being
    // above to below the desired pressure, or vice versa? If so, there is a
    // root.
    if ((P_below_lower - P) * (P_below_upper - P) <= 0.f) {

      // If there is a root, interpolate between the table values:
      intp_rho = (log_P - P_1_2) / (P_3_4 - P_1_2);

      root_below = expf(mat->table_log_rho[idx_rho] +
                        intp_rho * (mat->table_log_rho[idx_rho + 1] -
                                    mat->table_log_rho[idx_rho]));

      // If we found a root above, which one is closer to the reference rho?
      if (closest_root > 0.f) {
        if (fabsf(logf(root_below) - logf(rho_ref)) <
            fabsf(logf(closest_root) - logf(rho_ref))) {
          closest_root = root_below;
        }
      } else {
        closest_root = root_below;
      }
      // Found a root, so no need to search higher densities
      break;
    }
  }

  // Return the root if we found one
  if (closest_root > 0.f) {
    return closest_root;
  }

  // If we don't find a root before we reach max_counter, return rho_ref. Maybe
  // we should give an error here?
  return rho_sph;
}
#endif /* SWIFT_MIXED_HHE_HEAVY_EQUATION_OF_STATE_H */
