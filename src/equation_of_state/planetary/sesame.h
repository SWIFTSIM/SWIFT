/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_SESAME_EQUATION_OF_STATE_H
#define SWIFT_SESAME_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/sesame.h
 *
 * Contains the SESAME and ANEOS-in-SESAME-style EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 *
 */

/* Some standard headers. */
#include <math.h>
#include <float.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "units.h"
#include "utilities.h"

// SESAME parameters
struct SESAME_params {
  float *table_log_rho;
  float *table_log_u_rho_T;
  float *table_P_rho_T;
  float *table_c_rho_T;
  float *table_s_rho_T;
  int num_rho, num_T;
  float P_tiny, c_tiny;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material (cgs units)
INLINE static void set_SESAME_iron(struct SESAME_params *mat,
                                   enum eos_planetary_material_id mat_id) {
  // SESAME 2140
  mat->mat_id = mat_id;
}
INLINE static void set_SESAME_basalt(struct SESAME_params *mat,
                                     enum eos_planetary_material_id mat_id) {
  // SESAME 7530
  mat->mat_id = mat_id;
}
INLINE static void set_SESAME_water(struct SESAME_params *mat,
                                    enum eos_planetary_material_id mat_id) {
  // SESAME 7154
  mat->mat_id = mat_id;
}
INLINE static void set_SS08_water(struct SESAME_params *mat,
                                  enum eos_planetary_material_id mat_id) {
  // Senft & Stewart (2008)
  mat->mat_id = mat_id;
}
INLINE static void set_ANEOS_forsterite(struct SESAME_params *mat,
                                        enum eos_planetary_material_id mat_id) {
  // Stewart et al. (2019)
  mat->mat_id = mat_id;
}
INLINE static void set_ANEOS_iron(struct SESAME_params *mat,
                                  enum eos_planetary_material_id mat_id) {
  // Stewart (2020)
  mat->mat_id = mat_id;
}
INLINE static void set_ANEOS_Fe85Si15(struct SESAME_params *mat,
                                      enum eos_planetary_material_id mat_id) {
  // Stewart (2020)
  mat->mat_id = mat_id;
}

// Read the tables from file
INLINE static void load_table_SESAME(struct SESAME_params *mat,
                                     char *table_file) {

  // Load table contents from file
  FILE *f = fopen(table_file, "r");
  if (f == NULL) error("Failed to open the SESAME EoS file '%s'", table_file);

  // Ignore header lines
  char buffer[100];
  for (int i = 0; i < 5; i++) {
    if (fgets(buffer, 100, f) == NULL)
      error("Failed to read the SESAME EoS file header %s", table_file);
  }
  float ignore;

  // Table properties
  int c = fscanf(f, "%d %d", &mat->num_rho, &mat->num_T);
  if (c != 2) error("Failed to read the SESAME EoS table %s", table_file);

  // Ignore the first elements of rho = 0, T = 0
  mat->num_rho--;
  mat->num_T--;

  // Allocate table memory
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_u_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_P_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_c_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_s_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));

  // Densities (not log yet)
  for (int i_rho = -1; i_rho < mat->num_rho; i_rho++) {
    // Ignore the first elements of rho = 0, T = 0
    if (i_rho == -1) {
      c = fscanf(f, "%f", &ignore);
      if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
    } else {
      c = fscanf(f, "%f", &mat->table_log_rho[i_rho]);
      if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
    }
  }

  // Temperatures (ignored)
  for (int i_T = -1; i_T < mat->num_T; i_T++) {
    c = fscanf(f, "%f", &ignore);
    if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
  }

  // Sp. int. energies (not log yet), pressures, sound speeds, and entropies
  for (int i_T = -1; i_T < mat->num_T; i_T++) {
    for (int i_rho = -1; i_rho < mat->num_rho; i_rho++) {
      // Ignore the first elements of rho = 0, T = 0
      if ((i_T == -1) || (i_rho == -1)) {
        c = fscanf(f, "%f %f %f %f", &ignore, &ignore, &ignore, &ignore);
        if (c != 4) error("Failed to read the SESAME EoS table %s", table_file);
      } else {
        c = fscanf(f, "%f %f %f %f",
                   &mat->table_log_u_rho_T[i_rho * mat->num_T + i_T],
                   &mat->table_P_rho_T[i_rho * mat->num_T + i_T],
                   &mat->table_c_rho_T[i_rho * mat->num_T + i_T],
                   &mat->table_s_rho_T[i_rho * mat->num_T + i_T]);
        if (c != 4) error("Failed to read the SESAME EoS table %s", table_file);
      }
    }
  }

  fclose(f);
}

// Misc. modifications
INLINE static void prepare_table_SESAME(struct SESAME_params *mat) {

  // Convert densities to log(density)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    mat->table_log_rho[i_rho] = logf(mat->table_log_rho[i_rho]);
  }

  // Convert sp. int. energies to log(sp. int. energy)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = 0; i_T < mat->num_T; i_T++) {
      // If not positive then set very small for the log
      if (mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] <= 0) {
        mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] = 1.f;
      }

      mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] =
          logf(mat->table_log_u_rho_T[i_rho * mat->num_T + i_T]);
    }
  }

  // Initialise tiny pressure and soundspeed
  mat->P_tiny = FLT_MAX;
  mat->c_tiny = FLT_MAX;

  // Enforce that the 1D arrays of u (at each rho) are monotonic
  // This is necessary because, for some high-density u slices at very low T,
  // u decreases (very slightly) with T, which makes the interpolation fail
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = mat->num_T - 1; i_T > 0; i_T--) {

      // If the one-lower-T u is greater than this u
      if (mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] <
          mat->table_log_u_rho_T[i_rho * mat->num_T + i_T - 1]) {

        // Replace it and all elements below it with that value
        for (int j_u = 0; j_u < i_T; j_u++) {
          mat->table_log_u_rho_T[i_rho * mat->num_T + j_u] =
              mat->table_log_u_rho_T[i_rho * mat->num_T + i_T];
        }
        break;
      }

      // Smallest positive pressure and sound speed
      if ((mat->table_P_rho_T[i_rho * mat->num_T + i_T] < mat->P_tiny) &&
          (mat->table_P_rho_T[i_rho * mat->num_T + i_T] > 0)) {
        mat->P_tiny = mat->table_P_rho_T[i_rho * mat->num_T + i_T];
      }
      if ((mat->table_c_rho_T[i_rho * mat->num_T + i_T] < mat->c_tiny) &&
          (mat->table_c_rho_T[i_rho * mat->num_T + i_T] > 0)) {
        mat->c_tiny = mat->table_c_rho_T[i_rho * mat->num_T + i_T];
      }
    }
  }

  // Tiny pressure to allow interpolation near non-positive values
  mat->P_tiny *= 1e-3f;
  mat->c_tiny *= 1e-3f;
}

// Convert to internal units
INLINE static void convert_units_SESAME(struct SESAME_params *mat,
                                        const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  // All table values in SI
  // Densities (log)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    mat->table_log_rho[i_rho] +=
        logf(units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY) /
             units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  }

  // Sp. Int. Energies (log), pressures, and sound speeds
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = 0; i_T < mat->num_T; i_T++) {
      mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] += logf(
          units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
      mat->table_P_rho_T[i_rho * mat->num_T + i_T] *=
          units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
          units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
      mat->table_c_rho_T[i_rho * mat->num_T + i_T] *=
          units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
          units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
      mat->table_s_rho_T[i_rho * mat->num_T + i_T] *=
          units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
          units_cgs_conversion_factor(us, UNIT_CONV_ENTROPY);
    }
  }

  // Tiny pressure and sound speed
  mat->P_tiny *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
                 units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->c_tiny *= units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
                 units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
}

// gas_internal_energy_from_entropy
INLINE static float SESAME_internal_energy_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_pressure_from_entropy
INLINE static float SESAME_pressure_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float SESAME_entropy_from_pressure(
    float density, float pressure, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float SESAME_soundspeed_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float SESAME_entropy_from_internal_energy(
    float density, float u, const struct SESAME_params *mat) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float SESAME_pressure_from_internal_energy(
    float density, float u, const struct SESAME_params *mat) {

  float P, P_1, P_2, P_3, P_4;

  if (u <= 0.f) {
    return 0.f;
  }

  int idx_rho, idx_u_1, idx_u_2;
  float intp_rho, intp_u_1, intp_u_2;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u)) to find P(rho, u)
  // Density index
  idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  idx_u_1 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_rho_T + idx_rho * mat->num_T, mat->num_T);
  idx_u_2 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }
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

  // Check for duplicates in SESAME tables before interpolation
  if (mat->table_log_rho[idx_rho + 1] != mat->table_log_rho[idx_rho]) {
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
               (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] !=
      mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) /
        (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] -
         mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
      mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u - mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) /
        (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
         mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.;
  }

  // Table values
  P_1 = mat->table_P_rho_T[idx_rho * mat->num_T + idx_u_1];
  P_2 = mat->table_P_rho_T[idx_rho * mat->num_T + idx_u_1 + 1];
  P_3 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
  P_4 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

  // If more than two table values are non-positive then return zero
  int num_non_pos = 0;
  if (P_1 <= 0.f) num_non_pos++;
  if (P_2 <= 0.f) num_non_pos++;
  if (P_3 <= 0.f) num_non_pos++;
  if (P_4 <= 0.f) num_non_pos++;
  if (num_non_pos > 0) {
    // If just one or two are non-positive then replace them with a tiny value
    // Unless already trying to extrapolate in which case return zero
    if ((num_non_pos > 2) || (mat->P_tiny == 0.f) || (intp_rho < 0.f) || 
        (intp_u_1 < 0.f) || (intp_u_2 < 0.f)) {
      return 0.f;
    }
    if (P_1 <= 0.f) P_1 = mat->P_tiny;
    if (P_2 <= 0.f) P_2 = mat->P_tiny;
    if (P_3 <= 0.f) P_3 = mat->P_tiny;
    if (P_4 <= 0.f) P_4 = mat->P_tiny;
  }

  // Interpolate with the log values
  P_1 = logf(P_1);
  P_2 = logf(P_2);
  P_3 = logf(P_3);
  P_4 = logf(P_4);

  P = (1.f - intp_rho) * ((1.f - intp_u_1) * P_1 + intp_u_1 * P_2) +
      intp_rho * ((1.f - intp_u_2) * P_3 + intp_u_2 * P_4);

  // Convert back from log
  P = expf(P);

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float SESAME_internal_energy_from_pressure(
    float density, float P, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float SESAME_soundspeed_from_internal_energy(
    float density, float u, const struct SESAME_params *mat) {

  float c, c_1, c_2, c_3, c_4;

  if (u <= 0.f) {
    return 0.f;
  }

  int idx_rho, idx_u_1, idx_u_2;
  float intp_rho, intp_u_1, intp_u_2;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u)) to find c(rho, u)
  // Density index
  idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  idx_u_1 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_rho_T + idx_rho * mat->num_T, mat->num_T);
  idx_u_2 = find_value_in_monot_incr_array(
      log_u, mat->table_log_u_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }
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

  // Check for duplicates in SESAME tables before interpolation
  if (mat->table_log_rho[idx_rho + 1] != mat->table_log_rho[idx_rho]) {
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
               (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] !=
      mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) /
        (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] -
         mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
      mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u - mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) /
        (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
         mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.;
  }

  // Table values
  c_1 = mat->table_c_rho_T[idx_rho * mat->num_T + idx_u_1];
  c_2 = mat->table_c_rho_T[idx_rho * mat->num_T + idx_u_1 + 1];
  c_3 = mat->table_c_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
  c_4 = mat->table_c_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

  // If more than two table values are non-positive then return zero
  int num_non_pos = 0;
  if (c_1 <= 0.f) num_non_pos++;
  if (c_2 <= 0.f) num_non_pos++;
  if (c_3 <= 0.f) num_non_pos++;
  if (c_4 <= 0.f) num_non_pos++;
  if (num_non_pos > 2) {
    return mat->c_tiny;
  }
  // If just one or two are non-positive then replace them with a tiny value
  else if (num_non_pos > 0) {
    // Unless already trying to extrapolate in which case return zero
    if ((intp_rho < 0.f) || (intp_u_1 < 0.f) || (intp_u_2 < 0.f)) {
      return mat->c_tiny;
    }
    if (c_1 <= 0.f) c_1 = mat->c_tiny;
    if (c_2 <= 0.f) c_2 = mat->c_tiny;
    if (c_3 <= 0.f) c_3 = mat->c_tiny;
    if (c_4 <= 0.f) c_4 = mat->c_tiny;
  }

  // Interpolate with the log values
  c_1 = logf(c_1);
  c_2 = logf(c_2);
  c_3 = logf(c_3);
  c_4 = logf(c_4);

  c = (1.f - intp_rho) * ((1.f - intp_u_1) * c_1 + intp_u_1 * c_2) +
      intp_rho * ((1.f - intp_u_2) * c_3 + intp_u_2 * c_4);

  // Convert back from log
  c = expf(c);

  return c;
}

// gas_soundspeed_from_pressure
INLINE static float SESAME_soundspeed_from_pressure(
    float density, float P, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

#endif /* SWIFT_SESAME_EQUATION_OF_STATE_H */
