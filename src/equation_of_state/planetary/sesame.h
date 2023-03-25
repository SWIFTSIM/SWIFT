/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

// SESAME parameters
struct SESAME_params {
  float *table_log_rho;
  float *table_log_u_rho_T;
  float *table_P_rho_T;
  float *table_c_rho_T;
  float *table_log_s_rho_T;
  int version_date, num_rho, num_T;
  float u_tiny, P_tiny, c_tiny, s_tiny;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material
INLINE static void set_SESAME_iron(struct SESAME_params *mat,
                                   enum eos_planetary_material_id mat_id) {
  // SESAME 2140
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_SESAME_basalt(struct SESAME_params *mat,
                                     enum eos_planetary_material_id mat_id) {
  // SESAME 7530
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_SESAME_water(struct SESAME_params *mat,
                                    enum eos_planetary_material_id mat_id) {
  // SESAME 7154
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_SS08_water(struct SESAME_params *mat,
                                  enum eos_planetary_material_id mat_id) {
  // Senft & Stewart (2008)
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_ANEOS_forsterite(struct SESAME_params *mat,
                                        enum eos_planetary_material_id mat_id) {
  // Stewart et al. (2019)
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_ANEOS_iron(struct SESAME_params *mat,
                                  enum eos_planetary_material_id mat_id) {
  // Stewart (2020)
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_ANEOS_Fe85Si15(struct SESAME_params *mat,
                                      enum eos_planetary_material_id mat_id) {
  // Stewart (2020)
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}

// Generic user-provided custom materials
INLINE static void set_custom(struct SESAME_params *mat,
                              enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->version_date = 0;
}

/*
    Skip a line while reading a file.
*/
INLINE static int skip_line(FILE *f) {
  int c;

  // Read each character until reaching the end of the line or file
  do {
    c = fgetc(f);
  } while ((c != '\n') && (c != EOF));

  return c;
}

/*
    Skip n lines while reading a file.
*/
INLINE static int skip_lines(FILE *f, int n) {
  int c;

  for (int i = 0; i < n; i++) c = skip_line(f);

  return c;
}

/*
    Read the data from the table file.

    File contents (SESAME-like format plus header info etc)
    -------------
    # header (12 lines)
    version_date                                                (YYYYMMDD)
    num_rho  num_T
    rho[0]   rho[1]  ...  rho[num_rho]                          (kg/m^3)
    T[0]     T[1]    ...  T[num_T]                              (K)
    u[0, 0]                 P[0, 0]     c[0, 0]     s[0, 0]     (J/kg, Pa, m/s,
   J/K/kg) u[1, 0]                 ...         ...         ...
    ...                     ...         ...         ...
    u[num_rho-1, 0]         ...         ...         ...
    u[0, 1]                 ...         ...         ...
    ...                     ...         ...         ...
    u[num_rho-1, num_T-1]   ...         ...         s[num_rho-1, num_T-1]
*/
INLINE static void load_table_SESAME(struct SESAME_params *mat,
                                     char *table_file) {

  // Load table contents from file
  FILE *f = fopen(table_file, "r");
  if (f == NULL) error("Failed to open the SESAME EoS file '%s'", table_file);

  // Skip header lines
  skip_lines(f, 12);

  // Table properties
  int version_date;
  int c = fscanf(f, "%d", &version_date);
  if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
  if ((version_date != mat->version_date) && (mat->version_date != 0))
    error(
        "EoS file %s version_date %d does not match expected %d (YYYYMMDD)."
        "\nPlease download the file using "
        "examples/Planetary/EoSTables/get_eos_tables.sh",
        table_file, version_date, mat->version_date);
  c = fscanf(f, "%d %d", &mat->num_rho, &mat->num_T);
  if (c != 2) error("Failed to read the SESAME EoS table %s", table_file);

  // Ignore the first elements of rho = 0, T = 0
  mat->num_rho--;
  mat->num_T--;
  float ignore;

  // Allocate table memory
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_u_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_P_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_c_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_log_s_rho_T =
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

  // Sp. int. energies (not log yet), pressures, sound speeds, and sp.
  // entropies (not log yet)
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
                   &mat->table_log_s_rho_T[i_rho * mat->num_T + i_T]);
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

  // Initialise tiny values
  mat->u_tiny = FLT_MAX;
  mat->P_tiny = FLT_MAX;
  mat->c_tiny = FLT_MAX;
  mat->s_tiny = FLT_MAX;

  // Enforce that the 1D arrays of u (at each rho) are monotonic
  // This is necessary because, for some high-density u slices at very low T,
  // u decreases (very slightly) with T, which makes the interpolation fail
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = mat->num_T - 1; i_T > 0; i_T--) {

      // If the one-lower-T u is greater than this u
      if (mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] <
          mat->table_log_u_rho_T[i_rho * mat->num_T + i_T - 1]) {

        // Replace with this lower value
        mat->table_log_u_rho_T[i_rho * mat->num_T + i_T - 1] =
            mat->table_log_u_rho_T[i_rho * mat->num_T + i_T];
      }

      // Smallest positive values
      if ((mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] < mat->u_tiny) &&
          (mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] > 0)) {
        mat->u_tiny = mat->table_log_u_rho_T[i_rho * mat->num_T + i_T];
      }
      if ((mat->table_P_rho_T[i_rho * mat->num_T + i_T] < mat->P_tiny) &&
          (mat->table_P_rho_T[i_rho * mat->num_T + i_T] > 0)) {
        mat->P_tiny = mat->table_P_rho_T[i_rho * mat->num_T + i_T];
      }
      if ((mat->table_c_rho_T[i_rho * mat->num_T + i_T] < mat->c_tiny) &&
          (mat->table_c_rho_T[i_rho * mat->num_T + i_T] > 0)) {
        mat->c_tiny = mat->table_c_rho_T[i_rho * mat->num_T + i_T];
      }
      if ((mat->table_log_s_rho_T[i_rho * mat->num_T + i_T] < mat->s_tiny) &&
          (mat->table_log_s_rho_T[i_rho * mat->num_T + i_T] > 0)) {
        mat->s_tiny = mat->table_log_s_rho_T[i_rho * mat->num_T + i_T];
      }
    }
  }

  // Tiny values to allow interpolation near non-positive values
  mat->u_tiny *= 1e-3f;
  mat->P_tiny *= 1e-3f;
  mat->c_tiny *= 1e-3f;
  mat->s_tiny *= 1e-3f;

  // Convert sp. int. energies to log(sp. int. energy), same for sp. entropies
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = 0; i_T < mat->num_T; i_T++) {
      // If not positive then set very small for the log
      if (mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] <= 0) {
        mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] = mat->u_tiny;
      }

      mat->table_log_u_rho_T[i_rho * mat->num_T + i_T] =
          logf(mat->table_log_u_rho_T[i_rho * mat->num_T + i_T]);

      if (mat->table_log_s_rho_T[i_rho * mat->num_T + i_T] <= 0) {
        mat->table_log_s_rho_T[i_rho * mat->num_T + i_T] = mat->s_tiny;
      }

      mat->table_log_s_rho_T[i_rho * mat->num_T + i_T] =
          logf(mat->table_log_s_rho_T[i_rho * mat->num_T + i_T]);
    }
  }
}

// Convert to internal units
INLINE static void convert_units_SESAME(struct SESAME_params *mat,
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

  // Sp. int. energies (log), pressures, sound speeds, and sp. entropies
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
      mat->table_log_s_rho_T[i_rho * mat->num_T + i_T] +=
          logf(units_cgs_conversion_factor(
                   &si, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS) /
               units_cgs_conversion_factor(
                   us, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS));
    }
  }

  // Tiny values
  mat->u_tiny *=
      units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->P_tiny *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
                 units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->c_tiny *= units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
                 units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
  mat->s_tiny *=
      units_cgs_conversion_factor(&si,
                                  UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS);
}

// gas_internal_energy_from_entropy
INLINE static float SESAME_internal_energy_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  float u, log_u_1, log_u_2, log_u_3, log_u_4;

  if (entropy <= 0.f) {
    return 0.f;
  }

  int idx_rho, idx_s_1, idx_s_2;
  float intp_rho, intp_s_1, intp_s_2;
  const float log_rho = logf(density);
  const float log_s = logf(entropy);

  // 2D interpolation (bilinear with log(rho), log(s)) to find u(rho, s))
  // Density index
  idx_rho =
      find_value_in_monot_incr_array(log_rho, mat->table_log_rho, mat->num_rho);

  // Sp. entropy at this and the next density (in relevant slice of s array)
  idx_s_1 = find_value_in_monot_incr_array(
      log_s, mat->table_log_s_rho_T + idx_rho * mat->num_T, mat->num_T);
  idx_s_2 = find_value_in_monot_incr_array(
      log_s, mat->table_log_s_rho_T + (idx_rho + 1) * mat->num_T, mat->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho) {
    idx_rho = mat->num_rho - 2;
  }
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

  // Check for duplicates in SESAME tables before interpolation
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
  log_u_1 = mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_s_1];
  log_u_2 = mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_s_1 + 1];
  log_u_3 = mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2];
  log_u_4 = mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_s_2 + 1];

  // If below the minimum s at this rho then just use the lowest table values
  if ((idx_rho > 0.f) && ((intp_s_1 < 0.f) || (intp_s_2 < 0.f) ||
                          (log_u_1 > log_u_2) || (log_u_3 > log_u_4))) {
    intp_s_1 = 0;
    intp_s_2 = 0;
  }

  // Interpolate with the log values
  u = (1.f - intp_rho) * ((1.f - intp_s_1) * log_u_1 + intp_s_1 * log_u_2) +
      intp_rho * ((1.f - intp_s_2) * log_u_3 + intp_s_2 * log_u_4);

  // Convert back from log
  u = expf(u);

  return u;
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

  // 2D interpolation (bilinear with log(rho), log(u)) to find P(rho, u))
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
    intp_rho = 1.f;
  }
  if (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] !=
      mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) /
        (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] -
         mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.f;
  }
  if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
      mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u - mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) /
        (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
         mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.f;
  }

  // Table values
  P_1 = mat->table_P_rho_T[idx_rho * mat->num_T + idx_u_1];
  P_2 = mat->table_P_rho_T[idx_rho * mat->num_T + idx_u_1 + 1];
  P_3 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
  P_4 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

  // If below the minimum u at this rho then just use the lowest table values
  if ((idx_rho > 0.f) &&
      ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (P_1 > P_2) || (P_3 > P_4))) {
    intp_u_1 = 0;
    intp_u_2 = 0;
  }

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

  // 2D interpolation (bilinear with log(rho), log(u)) to find c(rho, u))
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
    intp_rho = 1.f;
  }
  if (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] !=
      mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]) /
        (mat->table_log_u_rho_T[idx_rho * mat->num_T + (idx_u_1 + 1)] -
         mat->table_log_u_rho_T[idx_rho * mat->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.f;
  }
  if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
      mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u - mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) /
        (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
         mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.f;
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

  // If below the minimum u at this rho then just use the lowest table values
  if ((idx_rho > 0.f) &&
      ((intp_u_1 < 0.f) || (intp_u_2 < 0.f) || (c_1 > c_2) || (c_3 > c_4))) {
    intp_u_1 = 0;
    intp_u_2 = 0;
  }

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
