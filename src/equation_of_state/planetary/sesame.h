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
  float *table_log_T;
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
INLINE static void set_AQUA(struct SESAME_params *mat,
                            enum eos_planetary_material_id mat_id) {
  // Haldemann et al. (2020)
  mat->mat_id = mat_id;
  mat->version_date = 20220714;
}
INLINE static void set_CMS19_H(struct SESAME_params *mat,
                               enum eos_planetary_material_id mat_id) {
  // Chabrier et al. (2019)
  mat->mat_id = mat_id;
  mat->version_date = 20220905;
}
INLINE static void set_CMS19_He(struct SESAME_params *mat,
                                enum eos_planetary_material_id mat_id) {
  // Chabrier et al. (2019)
  mat->mat_id = mat_id;
  mat->version_date = 20220905;
}
INLINE static void set_CD21_HHe(struct SESAME_params *mat,
                                enum eos_planetary_material_id mat_id) {
  // Chabrier & Debras (2021)
  mat->mat_id = mat_id;
  mat->version_date = 20220905;
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

  // Allocate table memory
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_T = (float *)malloc(mat->num_T * sizeof(float));
  mat->table_log_u_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_P_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_c_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));
  mat->table_log_s_rho_T =
      (float *)malloc(mat->num_rho * mat->num_T * sizeof(float));

  // Densities (not log yet)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    c = fscanf(f, "%f", &mat->table_log_rho[i_rho]);
    if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
  }

  // Temperatures (not log yet)
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    c = fscanf(f, "%f", &mat->table_log_T[i_T]);
    if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
  }

  // Ensure first density and temperature aren't zero before taking logs
  if (mat->table_log_rho[0] == 0) {
    mat->table_log_rho[0] = mat->table_log_rho[1] * 1e-5;
  }
  if (mat->table_log_T[0] == 0) {
    mat->table_log_T[0] = mat->table_log_T[1] * 1e-5;
  }

  // Sp. int. energies (not log yet), pressures, sound speeds, and sp.
  // entropies (not log yet)
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
      c = fscanf(f, "%f %f %f %f",
                 &mat->table_log_u_rho_T[i_rho * mat->num_T + i_T],
                 &mat->table_P_rho_T[i_rho * mat->num_T + i_T],
                 &mat->table_c_rho_T[i_rho * mat->num_T + i_T],
                 &mat->table_log_s_rho_T[i_rho * mat->num_T + i_T]);
      if (c != 4) error("Failed to read the SESAME EoS table %s", table_file);
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
  // Convert temperatures to log(temperature)
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    mat->table_log_T[i_T] = logf(mat->table_log_T[i_T]);
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

      // Ensure P > 0
      if (mat->table_P_rho_T[i_rho * mat->num_T + i_T] <= 0) {
        mat->table_P_rho_T[i_rho * mat->num_T + i_T] = mat->P_tiny;
      }
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

  // Temperatures (log)
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    mat->table_log_T[i_T] +=
        logf(units_cgs_conversion_factor(&si, UNIT_CONV_TEMPERATURE) /
             units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
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
    const float density, const float entropy, const struct SESAME_params *mat) {

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

  // Check for duplicates in SESAME tables before interpolation
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
INLINE static float SESAME_pressure_from_entropy(
    const float density, const float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float SESAME_entropy_from_pressure(
    const float density, const float pressure,
    const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float SESAME_soundspeed_from_entropy(
    const float density, const float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float SESAME_entropy_from_internal_energy(
    const float density, const float u, const struct SESAME_params *mat) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float SESAME_pressure_from_internal_energy(
    const float density, const float u, const struct SESAME_params *mat) {

  if (u <= 0.f) {
    return 0.f;
  }

  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u) to find P(rho, u))

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

  // Check for duplicates in SESAME tables before interpolation
  float intp_rho, intp_u_1, intp_u_2;
  const int idx_u_rho_T = idx_rho * mat->num_T + idx_u_1;
  if (mat->table_log_rho[idx_rho + 1] != mat->table_log_rho[idx_rho]) {
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
               (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.f;
  }
  if (mat->table_log_u_rho_T[idx_u_rho_T + 1] !=
      mat->table_log_u_rho_T[idx_u_rho_T]) {
    intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_rho_T]) /
               (mat->table_log_u_rho_T[idx_u_rho_T + 1] -
                mat->table_log_u_rho_T[idx_u_rho_T]);
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
  float P_1 = mat->table_P_rho_T[idx_u_rho_T];
  float P_2 = mat->table_P_rho_T[idx_u_rho_T + 1];
  float P_3 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
  float P_4 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

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
  const float P_1_2 = (1.f - intp_u_1) * P_1 + intp_u_1 * P_2;
  const float P_3_4 = (1.f - intp_u_2) * P_3 + intp_u_2 * P_4;

  const float P = (1.f - intp_rho) * P_1_2 + intp_rho * P_3_4;

  // Convert back from log
  return expf(P);
}

// gas_internal_energy_from_pressure
INLINE static float SESAME_internal_energy_from_pressure(
    const float density, const float P, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float SESAME_soundspeed_from_internal_energy(
    const float density, const float u, const struct SESAME_params *mat) {

  // Return zero if internal energy is non-positive
  if (u <= 0.f) {
    return 0.f;
  }

  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u) to find c(rho, u))

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

  // Check for duplicates in SESAME tables before interpolation
  float intp_rho, intp_u_1, intp_u_2;
  if (mat->table_log_rho[idx_rho + 1] != mat->table_log_rho[idx_rho]) {
    intp_rho = (log_rho - mat->table_log_rho[idx_rho]) /
               (mat->table_log_rho[idx_rho + 1] - mat->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.f;
  }
  const int idx_u_rho_T = idx_rho * mat->num_T + idx_u_1;
  if (mat->table_log_u_rho_T[idx_u_rho_T + 1] !=
      mat->table_log_u_rho_T[idx_u_rho_T]) {
    intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_rho_T]) /
               (mat->table_log_u_rho_T[idx_u_rho_T + 1] -
                mat->table_log_u_rho_T[idx_u_rho_T]);
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
  float c_1 = mat->table_c_rho_T[idx_u_rho_T];
  float c_2 = mat->table_c_rho_T[idx_u_rho_T + 1];
  float c_3 = mat->table_c_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
  float c_4 = mat->table_c_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

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

  const float c = (1.f - intp_rho) * ((1.f - intp_u_1) * c_1 + intp_u_1 * c_2) +
                  intp_rho * ((1.f - intp_u_2) * c_3 + intp_u_2 * c_4);

  // Convert back from log
  return expf(c);
}

// gas_soundspeed_from_pressure
INLINE static float SESAME_soundspeed_from_pressure(
    const float density, const float P, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_temperature_from_internal_energy
INLINE static float SESAME_temperature_from_internal_energy(
    const float density, const float u, const struct SESAME_params *mat) {

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

  // Check for duplicates in SESAME tables before interpolation
  float intp_u_1, intp_u_2;
  const int idx_u_rho_T = idx_rho * mat->num_T + idx_u_1;
  if (mat->table_log_u_rho_T[idx_u_rho_T + 1] !=
      mat->table_log_u_rho_T[idx_u_rho_T]) {
    intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_rho_T]) /
               (mat->table_log_u_rho_T[idx_u_rho_T + 1] -
                mat->table_log_u_rho_T[idx_u_rho_T]);
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
INLINE static float SESAME_density_from_pressure_and_temperature(
    float P, float T, const struct SESAME_params *mat) {

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

  // Check for duplicates in SESAME tables before interpolation
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
INLINE static float SESAME_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    const struct SESAME_params *mat) {

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
  int idx_u_rho_T;
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
    idx_u_rho_T = idx_rho * mat->num_T + idx_u_1;
    if (mat->table_log_u_rho_T[idx_u_rho_T + 1] !=
        mat->table_log_u_rho_T[idx_u_rho_T]) {
      intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_rho_T]) /
                 (mat->table_log_u_rho_T[idx_u_rho_T + 1] -
                  mat->table_log_u_rho_T[idx_u_rho_T]);
    } else {
      intp_u_1 = 1.f;
    }
    if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
        mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) {
      intp_u_2 =
          (log_u -
           mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) /
          (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
           mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]);
    } else {
      intp_u_2 = 1.f;
    }

    // Table values
    float P_1 = mat->table_P_rho_T[idx_u_rho_T];
    float P_2 = mat->table_P_rho_T[idx_u_rho_T + 1];
    float P_3 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
    float P_4 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

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
      if ((num_non_pos > 2) || (mat->P_tiny == 0.f) || (intp_u_1 < 0.f) ||
          (intp_u_2 < 0.f)) {
        break;  // return rho_sph;
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

    idx_u_rho_T = idx_rho * mat->num_T + idx_u_1;
    float intp_u_1, intp_u_2;
    if (mat->table_log_u_rho_T[idx_u_rho_T + 1] !=
        mat->table_log_u_rho_T[idx_u_rho_T]) {
      intp_u_1 = (log_u - mat->table_log_u_rho_T[idx_u_rho_T]) /
                 (mat->table_log_u_rho_T[idx_u_rho_T + 1] -
                  mat->table_log_u_rho_T[idx_u_rho_T]);
    } else {
      intp_u_1 = 1.f;
    }
    if (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] !=
        mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) {
      intp_u_2 =
          (log_u -
           mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]) /
          (mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + (idx_u_2 + 1)] -
           mat->table_log_u_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2]);
    } else {
      intp_u_2 = 1.f;
    }

    // Table values
    float P_1 = mat->table_P_rho_T[idx_u_rho_T];
    float P_2 = mat->table_P_rho_T[idx_u_rho_T + 1];
    float P_3 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2];
    float P_4 = mat->table_P_rho_T[(idx_rho + 1) * mat->num_T + idx_u_2 + 1];

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
      if ((num_non_pos > 2) || (mat->P_tiny == 0.f) || (intp_u_1 < 0.f) ||
          (intp_u_2 < 0.f)) {
        break;
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
#endif /* SWIFT_SESAME_EQUATION_OF_STATE_H */
