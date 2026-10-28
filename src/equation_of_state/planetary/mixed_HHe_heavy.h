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
  float *table_mix;
  float *table_log_rho;
  float *table_log_T;
  int num_mix, num_rho, num_T, num_rho_T;
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
        int idx = i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T;
        mat->table_log_u_mix_rho_T[idx] = logf(mat->table_log_u_mix_rho_T[idx]);
        mat->table_log_c_mix_rho_T[idx] = logf(mat->table_log_c_mix_rho_T[idx]);
        mat->table_log_s_mix_rho_T[idx] = logf(mat->table_log_s_mix_rho_T[idx]);
      }
    }
  }
}
INLINE static void prepare_mixedHHeHeavy(struct mixedHHeHeavy_params *mat) {
  // Same array dimensions for all mix components
  mat->num_mix = mat->single_params[0].num_mix;
  mat->num_rho = mat->single_params[0].num_rho;
  mat->num_T = mat->single_params[0].num_T;
  mat->num_rho_T = mat->single_params[0].num_rho_T;

  // Check
  for (int i = 1; i < mat->num_mat; i++) {
    if (mat->num_mix != mat->single_params[i].num_mix)
      error("Mixed EoS mix-table dimensions don't match");
    if (mat->num_rho != mat->single_params[i].num_rho)
      error("Mixed EoS density-table dimensions don't match");
    if (mat->num_T != mat->single_params[i].num_T)
      error("Mixed EoS temperature-table dimensions don't match");
  }

  // Same 1D axis arrays for all mix components
  mat->table_mix = (float *)malloc(mat->num_mix * sizeof(float));
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_T = (float *)malloc(mat->num_T * sizeof(float));
  for (int i = 0; i < mat->num_mix; i++) {
    mat->table_mix[i] = mat->single_params[0].table_mix[i];
  }
  for (int i = 0; i < mat->num_rho; i++) {
    mat->table_log_rho[i] = mat->single_params[0].table_log_rho[i];
  }
  for (int i = 0; i < mat->num_T; i++) {
    mat->table_log_T[i] = mat->single_params[0].table_log_T[i];
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
        int idx = i_mix * mat->num_rho_T + i_rho * mat->num_T + i_T;
        mat->table_log_u_mix_rho_T[idx] += 
            logf(units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
                 units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
        mat->table_P_mix_rho_T[idx] *=
            units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
            units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
        mat->table_log_c_mix_rho_T[idx] +=
            logf(units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
                units_cgs_conversion_factor(us, UNIT_CONV_SPEED));
        mat->table_log_s_mix_rho_T[idx] +=
            logf(units_cgs_conversion_factor(
                    &si, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS) /
                units_cgs_conversion_factor(
                    us, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS));
      }
    }
  }
}

/*
    Find a value Z from the density and one other variable Y:

    1. Combine the mix components' tables into a single table at this mix.
    2. Interpolate each temperature column linearly in log(rho) to get Y(T) at
       this particle's density.
    3. Invert that single Y(T) array to find T.
    4. Interpolate Z linearly in log(rho), log(T) at that same (rho, T).
*/

// The table parameter choices
enum mixedHHeHeavy_param_choice {
  param_choice_P,
  param_choice_u,
  param_choice_c,
  param_choice_s,
  param_choice_T,
};

// The mix-table position and each component's weight for one particle's mixes
struct mixedHHeHeavy_mix_weights {
  int idx_mix;
  float intp_mix;
  float f_mat[eos_planetary_mixed_num_mat];
};

// Set idx_mix, intp_mix, and the weight of each mix component's contribution
INLINE static void mixedHHeHeavy_set_mix_weights(
    struct mixedHHeHeavy_mix_weights *w, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  // Mix fraction with respect to HHe rather than between each heavy component
  float mix_tot = 0.f;
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    mix_tot += mixes[i_mat];
  }

  // Relative weight of each heavy component (equal if no heavy elements)
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    if (mix_tot > 0.f) {
      w->f_mat[i_mat] = mixes[i_mat] / mix_tot;
    } else {
      w->f_mat[i_mat] = 1.f / mat->num_mat;
    }
  }

  // Mix, clamped to the table
  if ((mix_tot >= 0.9999f) || (mix_tot >= mat->table_mix[mat->num_mix - 1])) {
    w->idx_mix = mat->num_mix - 2;
    w->intp_mix = 1.f;
  } else if (mix_tot <= mat->table_mix[0]) {
    w->idx_mix = 0;
    w->intp_mix = 0.f;
  } else {
    w->idx_mix =
        find_value_in_monot_incr_array(mix_tot, mat->table_mix, mat->num_mix);
    w->intp_mix = (mix_tot - mat->table_mix[w->idx_mix]) /
                  (mat->table_mix[w->idx_mix + 1] - mat->table_mix[w->idx_mix]);
  }
}

// log(Z) at a table (rho, T) entry, interpolated over the mix components
INLINE static float mixedHHeHeavy_log_Z_rho_T(
    const struct mixedHHeHeavy_params *mat,
    const struct mixedHHeHeavy_mix_weights *w, const enum mixedHHeHeavy_param_choice Z_choice,
    const int i_rho, const int i_T) {

  const int idx_1 = w->idx_mix * mat->num_rho_T + i_rho * mat->num_T + i_T;
  const int idx_2 = idx_1 + mat->num_rho_T;

  float Z = 0.f;
  for (int i_mat = 0; i_mat < mat->num_mat; i_mat++) {
    if (w->f_mat[i_mat] == 0.f) continue;

    const struct single_mixedHHeHeavy_params *mat_single = &mat->single_params[i_mat];
    float log_Z_1, log_Z_2;
    switch (Z_choice) {
      case param_choice_u:
        log_Z_1 = mat_single->table_log_u_mix_rho_T[idx_1];
        log_Z_2 = mat_single->table_log_u_mix_rho_T[idx_2];
        break;
      case param_choice_P:
        // Avoid log(-ve) -- though that shouldn't occur in these tables
        log_Z_1 = logf(fmaxf(mat_single->table_P_mix_rho_T[idx_1], FLT_MIN));
        log_Z_2 = logf(fmaxf(mat_single->table_P_mix_rho_T[idx_2], FLT_MIN));
        break;
      case param_choice_c:
        log_Z_1 = mat_single->table_log_c_mix_rho_T[idx_1];
        log_Z_2 = mat_single->table_log_c_mix_rho_T[idx_2];
        break;
      case param_choice_s:
        log_Z_1 = mat_single->table_log_s_mix_rho_T[idx_1];
        log_Z_2 = mat_single->table_log_s_mix_rho_T[idx_2];
        break;
      default:
        error("Invalid mixedHHeHeavy table variable %d", Z_choice);
        return 0.f;
    }

    Z += w->f_mat[i_mat] * expf((1.f - w->intp_mix) * log_Z_1 + w->intp_mix * log_Z_2);
  }

  return logf(Z);
}

// log(Z) from the density and another variable Y (Z can also be T)
INLINE static float mixedHHeHeavy_log_Z_rho_Y(
    const float density, const float Y, const float *mixes,
    const struct mixedHHeHeavy_params *mat, const enum mixedHHeHeavy_param_choice Y_choice,
    const enum mixedHHeHeavy_param_choice Z_choice) {

  struct mixedHHeHeavy_mix_weights w;
  mixedHHeHeavy_set_mix_weights(&w, mixes, mat);

  const float log_rho = logf(density);
  const float log_Y = logf(Y);

  // Density
  int idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);
  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho - 1) {
    idx_rho = mat->num_rho - 2;
  }
  const float intp_rho =
      (log_rho - mat->table_log_rho[idx_rho]) /
      (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);

  // Y(T) at this density
  float log_Y_T[mat->num_T];
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    log_Y_T[i_T] =
        (1.f - intp_rho) *
            mixedHHeHeavy_log_Z_rho_T(mat, &w, Y_choice, idx_rho, i_T) +
        intp_rho *
            mixedHHeHeavy_log_Z_rho_T(mat, &w, Y_choice, idx_rho + 1, i_T);
  }

  // Temperature
  int idx_T = find_value_in_monot_incr_array(log_Y, log_Y_T, mat->num_T);
  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_T <= -1) {
    idx_T = 0;
  } else if (idx_T >= mat->num_T - 1) {
    idx_T = mat->num_T - 2;
  }
  const float delta_log_Y = log_Y_T[idx_T + 1] - log_Y_T[idx_T];
  float intp_T = 0.f;
  if (delta_log_Y != 0.f) intp_T = (log_Y - log_Y_T[idx_T]) / delta_log_Y;

  // If below the minimum Y at this rho then just use the lowest table values
  if ((idx_rho > 0) && (intp_T < 0.f)) intp_T = 0.f;

  if (Z_choice == param_choice_T) {
    return (1.f - intp_T) * mat->table_log_T[idx_T] +
           intp_T * mat->table_log_T[idx_T + 1];
  }

  // Interpolate Z at the same (rho, T)
  const float log_Z_1 =
      mixedHHeHeavy_log_Z_rho_T(mat, &w, Z_choice, idx_rho, idx_T);
  const float log_Z_2 =
      mixedHHeHeavy_log_Z_rho_T(mat, &w, Z_choice, idx_rho, idx_T + 1);
  const float log_Z_3 =
      mixedHHeHeavy_log_Z_rho_T(mat, &w, Z_choice, idx_rho + 1, idx_T);
  const float log_Z_4 =
      mixedHHeHeavy_log_Z_rho_T(mat, &w, Z_choice, idx_rho + 1, idx_T + 1);

  return (1.f - intp_rho) * ((1.f - intp_T) * log_Z_1 + intp_T * log_Z_2) +
         intp_rho * ((1.f - intp_T) * log_Z_3 + intp_T * log_Z_4);
}

// gas_internal_energy_from_entropy
INLINE static float mixedHHeHeavy_internal_energy_from_entropy(
    const float density, const float entropy, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (entropy <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, entropy, mixes, mat,
                                        param_choice_s,
                                        param_choice_u));
}

// gas_pressure_from_entropy
INLINE static float mixedHHeHeavy_pressure_from_entropy(
    const float density, const float entropy, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (entropy <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, entropy, mixes, mat,
                                        param_choice_s,
                                        param_choice_P));
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
    const float density, const float entropy, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (entropy <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, entropy, mixes, mat,
                                        param_choice_s,
                                        param_choice_c));
}

// gas_entropy_from_internal_energy
INLINE static float mixedHHeHeavy_entropy_from_internal_energy(
    const float density, const float u, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, u, mixes, mat,
                                        param_choice_u,
                                        param_choice_s));
}

// gas_pressure_from_internal_energy
INLINE static float mixedHHeHeavy_pressure_from_internal_energy(
    const float density, const float u, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, u, mixes, mat,
                                        param_choice_u,
                                        param_choice_P));
}

// gas_internal_energy_from_pressure
INLINE static float mixedHHeHeavy_internal_energy_from_pressure(
    const float density, const float P, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float mixedHHeHeavy_soundspeed_from_internal_energy(
    const float density, const float u, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, u, mixes, mat,
                                        param_choice_u,
                                        param_choice_c));
}

// gas_soundspeed_from_pressure
INLINE static float mixedHHeHeavy_soundspeed_from_pressure(
    const float density, const float P, const float *mixes, const struct mixedHHeHeavy_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_temperature_from_internal_energy
INLINE static float mixedHHeHeavy_temperature_from_internal_energy(
    const float density, const float u, const float *mixes,
    const struct mixedHHeHeavy_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  return expf(mixedHHeHeavy_log_Z_rho_Y(density, u, mixes, mat,
                                        param_choice_u,
                                        param_choice_T));
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
