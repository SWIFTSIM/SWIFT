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
  float *table_log_c_mix_rho_T;
  float *table_log_s_mix_rho_T;
  int version_date, num_mix, num_rho, num_T, num_rho_T;
  enum eos_planetary_material_id mat_id;
};
struct mixedHHeHeavy_params {
  struct single_mixedHHeHeavy_params *single_params;
  int num_mat;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material
INLINE static void set_mixedHHeHeavy(struct mixedHHeHeavy_params *mat,
                                     enum eos_planetary_material_id mat_id,
                                     int num_mat) {
  mat->mat_id = mat_id;
  mat->num_mat = num_mat;
  mat->single_params = malloc(mat->num_mat * sizeof(struct single_mixedHHeHeavy_params));
}
INLINE static void set_mixedHHe_rock(struct single_mixedHHeHeavy_params *mat,
                                          enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 20250128;
}
INLINE static void set_mixedHHe_water(struct single_mixedHHeHeavy_params *mat,
                                           enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 20250128;
}
INLINE static void set_mixedHHe_iron(struct single_mixedHHeHeavy_params *mat,
                                          enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 20260203;
}

/*
    Read the table data from the HDF5 file.
*/
INLINE static void load_table_mixedHHeHeavy(struct single_mixedHHeHeavy_params *mat,
                                     char *table_file) {
                                      
#ifndef HAVE_HDF5
  error("Need HDF5 to read mixed EoS tables");
#endif

  // Load table file
  hid_t h_file = H5Fopen(table_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) error("Failed to open mixedHHeHeavy EoS file %s\n", table_file);

  // Load header data
  hid_t grp_header = H5Gopen(h_file, "/Header", H5P_DEFAULT);
  int version_date;
  io_read_attribute(grp_header, "VersionDate", INT, &version_date);
  if ((version_date != mat->version_date) && (mat->version_date != 0))
    error(
        "EoS file %s version_date %d does not match expected %d (YYYYMMDD)."
        "\nPlease download the file using "
        "examples/Planetary/EoSTables/get_eos_tables.sh",
        table_file, version_date, mat->version_date);
  io_read_attribute(grp_header, "MixMassFractionCount", INT, &mat->num_mix);
  io_read_attribute(grp_header, "DensityCount", INT, &mat->num_rho);
  io_read_attribute(grp_header, "TemperatureCount", INT, &mat->num_T);
  mat->num_rho_T = mat->num_rho * mat->num_T;
  const int num_mix_rho_T = mat->num_mix * mat->num_rho_T;

  // Allocate table memory
  mat->table_mix = (float *)malloc(mat->num_mix * sizeof(float));
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_T = (float *)malloc(mat->num_T * sizeof(float));
  mat->table_log_u_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));
  mat->table_P_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));
  mat->table_log_c_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));
  mat->table_log_s_mix_rho_T = (float *)malloc(num_mix_rho_T * sizeof(float));

  // Load table data (not log yet)
  hid_t grp_table = H5Gopen(h_file, "/Table", H5P_DEFAULT);
  io_read_array_dataset(grp_table, "MixMassFraction", FLOAT, mat->table_mix, mat->num_mix);
  io_read_array_dataset(grp_table, "Density", FLOAT, mat->table_log_rho, mat->num_rho);
  io_read_array_dataset(grp_table, "Temperature", FLOAT, mat->table_log_T, mat->num_T);
  io_read_array_dataset(grp_table, "SpecificInternalEnergy", FLOAT, mat->table_log_u_mix_rho_T, num_mix_rho_T);
  io_read_array_dataset(grp_table, "Pressure", FLOAT, mat->table_P_mix_rho_T, num_mix_rho_T);
  io_read_array_dataset(grp_table, "SoundSpeed", FLOAT, mat->table_log_c_mix_rho_T, num_mix_rho_T);
  io_read_array_dataset(grp_table, "SpecificEntropy", FLOAT, mat->table_log_s_mix_rho_T, num_mix_rho_T);
}

// Misc. modifications
INLINE static void prepare_table_mixedHHeHeavy(struct single_mixedHHeHeavy_params *mat) {

  // Convert to log
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    mat->table_log_rho[i_rho] = logf(mat->table_log_rho[i_rho]);
  }  
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    mat->table_log_T[i_T] = logf(mat->table_log_T[i_T]);
  }
  for (int i_mix = 0; i_mix < mat->num_mix; i_mix++) {
    for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
      for (int i_T = 0; i_T < mat->num_T; i_T++) {
        mat->table_log_u_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] =
            logf(mat->table_log_u_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T]);

            mat->table_log_c_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] =
                logf(mat->table_log_c_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T]);

            mat->table_log_s_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] =
                logf(mat->table_log_s_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T]);
      }
    }
  }
}

// Convert to internal units
INLINE static void convert_units_mixedHHeHeavy(struct single_mixedHHeHeavy_params *mat,
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
        mat->table_log_c_mix_rho_T[i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T] *=
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
    const float density, const float entropy, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_pressure_from_entropy
INLINE static float mixedHHeHeavy_pressure_from_entropy(
    const float density, const float entropy, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float mixedHHeHeavy_entropy_from_pressure(
    const float density, const float pressure,
    const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float mixedHHeHeavy_soundspeed_from_entropy(
    const float density, const float entropy, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float mixedHHeHeavy_entropy_from_internal_energy(
    const float density, const float u, const float *mixes, const struct mixedHHeHeavy_params *mat) {

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
  intp_u_1 = (log_u - mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]);
  intp_u_2 = (log_u - mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]);

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
    idx_u_1 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_mix_rho_T + (idx_mix + 1) * mat->num_rho_T
              + idx_rho * mat->num_T, mat->num_T);
    idx_u_2 = find_value_in_monot_incr_array(
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
    intp_u_1 = (log_u - mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]);
    intp_u_2 = (log_u - mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]);

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
    const float density, const float u, const float *mixes, 
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  // No heavy-element fraction
  float mix_tot = 0.f;
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    mix_tot += mixes[i_mat];
  }
  if (mix_tot == 0.f) {
    return single_mixedHHeHeavy_pressure_from_internal_energy(
      density, u, 0.f, &mat->single_params[0]);
  }

  // Accumulate contribution from each non-zero mix
  float P = 0.f;
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    float mix = mixes[i_mat];
    if (mix > 0.f) {
      // Evaluate for this single heavy mix
      float P_mat = single_mixedHHeHeavy_pressure_from_internal_energy(
        density, u, mix, &mat->single_params[i_mat]);

      P += P_mat * mix / mix_tot;
    }
  }

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float mixedHHeHeavy_internal_energy_from_pressure(
    const float density, const float P, const float *mixes, const struct mixedHHeHeavy_params *mat) {

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
  intp_u_1 = (log_u - mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]);
  intp_u_2 = (log_u - mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]) /
              (mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T + 1] -
              mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]);

  // Table values
  float log_c_1 = mat->table_log_c_mix_rho_T[idx_u_1_mix_rho_T];
  float log_c_2 = mat->table_log_c_mix_rho_T[idx_u_1_mix_rho_T + 1];
  float log_c_3 = mat->table_log_c_mix_rho_T[idx_u_2_mix_rho_T];
  float log_c_4 = mat->table_log_c_mix_rho_T[idx_u_2_mix_rho_T + 1];

  // If below the minimum u at this rho then just use the lowest table values
  if ((idx_rho > 0.f) &&
      ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (log_c_1 > log_c_2) || (log_c_3 > log_c_4))) {
    intp_u_1 = 0;
    intp_u_2 = 0;
  }

  // Interpolate with the log values
  float log_c_1_2 = (1.f - intp_u_1) * log_c_1 + intp_u_1 * log_c_2;
  float log_c_3_4 = (1.f - intp_u_2) * log_c_3 + intp_u_2 * log_c_4;

  float log_c = (1.f - intp_rho) * log_c_1_2 + intp_rho * log_c_3_4;

  // 2D interpolate within second single-mix table, then combine
  // (todo: tidy into a generalised function!)
  if ((mix > 0.f) && (mix < 1.f)) {
    // Sp. int. energy at this and the next density (in relevant slice of each u array)
    idx_u_1 = find_value_in_monot_incr_array(
        log_u, mat->table_log_u_mix_rho_T + (idx_mix + 1) * mat->num_rho_T
              + idx_rho * mat->num_T, mat->num_T);
    idx_u_2 = find_value_in_monot_incr_array(
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
    intp_u_1 = (log_u - mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T[idx_u_1_mix_rho_T]);
    intp_u_2 = (log_u - mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]) /
                (mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T + 1] -
                mat->table_log_u_mix_rho_T[idx_u_2_mix_rho_T]);

    // Table values
    log_c_1 = mat->table_log_c_mix_rho_T[idx_u_1_mix_rho_T];
    log_c_2 = mat->table_log_c_mix_rho_T[idx_u_1_mix_rho_T + 1];
    log_c_3 = mat->table_log_c_mix_rho_T[idx_u_2_mix_rho_T];
    log_c_4 = mat->table_log_c_mix_rho_T[idx_u_2_mix_rho_T + 1];

    // If below the minimum u at this rho then just use the lowest table values
    if ((idx_rho > 0.f) &&
        ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (log_c_1 > log_c_2) || (log_c_3 > log_c_4))) {
      intp_u_1 = 0;
      intp_u_2 = 0;
    }

    // Interpolate with the log values
    log_c_1_2 = (1.f - intp_u_1) * log_c_1 + intp_u_1 * log_c_2;
    log_c_3_4 = (1.f - intp_u_2) * log_c_3 + intp_u_2 * log_c_4;

    const float log_c_mix2 = (1.f - intp_rho) * log_c_1_2 + intp_rho * log_c_3_4;

    // Interpolate between mix values
    log_c = (1 - intp_mix) * log_c + intp_mix * log_c_mix2;
  }

  // Convert back from log
  return expf(log_c);
}
  
INLINE static float mixedHHeHeavy_soundspeed_from_internal_energy(
    const float density, const float u, const float *mixes, 
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  // No heavy-element fraction
  float mix_tot = 0.f;
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    mix_tot += mixes[i_mat];
  }
  if (mix_tot == 0.f) {
    return single_mixedHHeHeavy_soundspeed_from_internal_energy(
      density, u, 0.f, &mat->single_params[0]);
  }

  // Accumulate contribution from each non-zero mix
  float c = 0.f;
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    float mix = mixes[i_mat];
    if (mix > 0.f) {
      // Evaluate for this single heavy mix
      float c_mat = single_mixedHHeHeavy_soundspeed_from_internal_energy(
        density, u, mix, &mat->single_params[i_mat]);

      c += c_mat * mix / mix_tot;
    }
  }

  return c;
}

// gas_soundspeed_from_pressure
INLINE static float mixedHHeHeavy_soundspeed_from_pressure(
    const float density, const float P, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_temperature_from_internal_energy
INLINE static float mixedHHeHeavy_temperature_from_internal_energy(
    const float density, const float u, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_density_from_pressure_and_temperature
INLINE static float mixedHHeHeavy_density_from_pressure_and_temperature(
    float P, float T, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_density_from_pressure_and_internal_energy
INLINE static float mixedHHeHeavy_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

#endif /* SWIFT_MIXED_HHE_HEAVY_EQUATION_OF_STATE_H */
