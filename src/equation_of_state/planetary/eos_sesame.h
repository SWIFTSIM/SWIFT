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
#include "eos_setup.h"
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
INLINE static void set_SESAME_iron(struct SESAME_params *SESAME,
                                   enum eos_planetary_material_id mat_id) {
  // SESAME 2140
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_SESAME_basalt(struct SESAME_params *SESAME,
                                     enum eos_planetary_material_id mat_id) {
  // SESAME 7530
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_SESAME_water(struct SESAME_params *SESAME,
                                    enum eos_planetary_material_id mat_id) {
  // SESAME 7154
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_SS08_water(struct SESAME_params *SESAME,
                                  enum eos_planetary_material_id mat_id) {
  // Senft & Stewart (2008)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_AQUA(struct SESAME_params *SESAME,
                            enum eos_planetary_material_id mat_id) {
  // Haldemann et al. (2020)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_CMS19_H(struct SESAME_params *SESAME,
                               enum eos_planetary_material_id mat_id) {
  // Chabrier et al. (2019)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220905;
}
INLINE static void set_CMS19_He(struct SESAME_params *SESAME,
                                enum eos_planetary_material_id mat_id) {
  // Chabrier et al. (2019)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220905;
}
INLINE static void set_CD21_HHe(struct SESAME_params *SESAME,
                                enum eos_planetary_material_id mat_id) {
  // Chabrier & Debras (2021)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220905;
}
INLINE static void set_ANEOS_forsterite(struct SESAME_params *SESAME,
                                        enum eos_planetary_material_id mat_id) {
  // Stewart et al. (2019)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_ANEOS_iron(struct SESAME_params *SESAME,
                                  enum eos_planetary_material_id mat_id) {
  // Stewart (2020)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}
INLINE static void set_ANEOS_Fe85Si15(struct SESAME_params *SESAME,
                                      enum eos_planetary_material_id mat_id) {
  // Stewart (2020)
  SESAME->mat_id = mat_id;
  SESAME->version_date = 20220714;
}

// Generic user-provided custom materials
INLINE static void set_custom(struct SESAME_params *SESAME,
                              enum eos_planetary_material_id mat_id) {
  SESAME->mat_id = mat_id;
  SESAME->version_date = 0;
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
INLINE static void load_table_SESAME(struct SESAME_params *SESAME,
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
  if ((version_date != SESAME->version_date) && (SESAME->version_date != 0))
    error(
        "EoS file %s version_date %d does not match expected %d (YYYYMMDD)."
        "\nPlease download the file using "
        "examples/Planetary/EoSTables/get_eos_tables.sh",
        table_file, version_date, SESAME->version_date);
  c = fscanf(f, "%d %d", &SESAME->num_rho, &SESAME->num_T);
  if (c != 2) error("Failed to read the SESAME EoS table %s", table_file);

  // Ignore the first elements of rho = 0, T = 0
  SESAME->num_rho--;
  SESAME->num_T--;
  float ignore;

  // Allocate table memory
  SESAME->table_log_rho = (float *)malloc(SESAME->num_rho * sizeof(float));
  SESAME->table_log_T = (float *)malloc(SESAME->num_T * sizeof(float));
  SESAME->table_log_u_rho_T =
      (float *)malloc(SESAME->num_rho * SESAME->num_T * sizeof(float));
  SESAME->table_P_rho_T =
      (float *)malloc(SESAME->num_rho * SESAME->num_T * sizeof(float));
  SESAME->table_c_rho_T =
      (float *)malloc(SESAME->num_rho * SESAME->num_T * sizeof(float));
  SESAME->table_log_s_rho_T =
      (float *)malloc(SESAME->num_rho * SESAME->num_T * sizeof(float));

  // Densities (not log yet)
  for (int i_rho = -1; i_rho < SESAME->num_rho; i_rho++) {
    // Ignore the first elements of rho = 0, T = 0
    if (i_rho == -1) {
      c = fscanf(f, "%f", &ignore);
      if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
    } else {
      c = fscanf(f, "%f", &SESAME->table_log_rho[i_rho]);
      if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
    }
  }

  // Temperatures (not log yet)
  for (int i_T = -1; i_T < SESAME->num_T; i_T++) {
    // Ignore the first elements of rho = 0, T = 0
    if (i_T == -1) {
      c = fscanf(f, "%f", &ignore);
      if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
    } else {
      c = fscanf(f, "%f", &SESAME->table_log_T[i_T]);
      if (c != 1) error("Failed to read the SESAME EoS table %s", table_file);
    }
  }

  // Sp. int. energies (not log yet), pressures, sound speeds, and sp.
  // entropies (not log yet)
  for (int i_T = -1; i_T < SESAME->num_T; i_T++) {
    for (int i_rho = -1; i_rho < SESAME->num_rho; i_rho++) {
      // Ignore the first elements of rho = 0, T = 0
      if ((i_T == -1) || (i_rho == -1)) {
        c = fscanf(f, "%f %f %f %f", &ignore, &ignore, &ignore, &ignore);
        if (c != 4) error("Failed to read the SESAME EoS table %s", table_file);
      } else {
        c = fscanf(f, "%f %f %f %f",
                   &SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T],
                   &SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T],
                   &SESAME->table_c_rho_T[i_rho * SESAME->num_T + i_T],
                   &SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T]);
        if (c != 4) error("Failed to read the SESAME EoS table %s", table_file);
      }
    }
  }

  fclose(f);
}

// Misc. modifications
INLINE static void prepare_table_SESAME(struct SESAME_params *SESAME) {

  // Convert densities to log(density)
  for (int i_rho = 0; i_rho < SESAME->num_rho; i_rho++) {
    SESAME->table_log_rho[i_rho] = logf(SESAME->table_log_rho[i_rho]);
  }
  // Convert temperatures to log(temperature)
  for (int i_T = 0; i_T < SESAME->num_T; i_T++) {
    SESAME->table_log_T[i_T] = logf(SESAME->table_log_T[i_T]);
  }

  // Initialise tiny values
  SESAME->u_tiny = FLT_MAX;
  SESAME->P_tiny = FLT_MAX;
  SESAME->c_tiny = FLT_MAX;
  SESAME->s_tiny = FLT_MAX;

  // Enforce that the 1D arrays of u (at each rho) are monotonic
  // This is necessary because, for some high-density u slices at very low T,
  // u decreases (very slightly) with T, which makes the interpolation fail
  for (int i_rho = 0; i_rho < SESAME->num_rho; i_rho++) {
    for (int i_T = SESAME->num_T - 1; i_T > 0; i_T--) {

      // If the one-lower-T u is greater than this u
      if (SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] <
          SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T - 1]) {

        // Replace with this lower value
        SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T - 1] =
            SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T];
      }

      // Smallest positive values
      if ((SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] <
           SESAME->u_tiny) &&
          (SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] > 0)) {
        SESAME->u_tiny = SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T];
      }
      if ((SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T] <
           SESAME->P_tiny) &&
          (SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T] > 0)) {
        SESAME->P_tiny = SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T];
      }
      if ((SESAME->table_c_rho_T[i_rho * SESAME->num_T + i_T] <
           SESAME->c_tiny) &&
          (SESAME->table_c_rho_T[i_rho * SESAME->num_T + i_T] > 0)) {
        SESAME->c_tiny = SESAME->table_c_rho_T[i_rho * SESAME->num_T + i_T];
      }
      if ((SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T] <
           SESAME->s_tiny) &&
          (SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T] > 0)) {
        SESAME->s_tiny = SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T];
      }
    }
  }

  // Tiny values to allow interpolation near non-positive values
  SESAME->u_tiny *= 1e-3f;
  SESAME->P_tiny *= 1e-3f;
  SESAME->c_tiny *= 1e-3f;
  SESAME->s_tiny *= 1e-3f;

  // Convert sp. int. energies to log(sp. int. energy), same for sp. entropies
  for (int i_rho = 0; i_rho < SESAME->num_rho; i_rho++) {
    for (int i_T = 0; i_T < SESAME->num_T; i_T++) {
      // If not positive then set very small for the log
      if (SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] <= 0) {
        SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] = SESAME->u_tiny;
      }

      SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] =
          logf(SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T]);

      if (SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T] <= 0) {
        SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T] = SESAME->s_tiny;
      }

      SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T] =
          logf(SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T]);

      // Ensure P > 0
      if (SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T] <= 0) {
        SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T] = SESAME->P_tiny;
      }
    }
  }
}

// Convert to internal units
INLINE static void convert_units_SESAME(struct SESAME_params *SESAME,
                                        const struct unit_system *us) {

  // Convert input table values from all-SI to internal units
  struct unit_system si;
  units_init_si(&si);

  // Densities (log)
  for (int i_rho = 0; i_rho < SESAME->num_rho; i_rho++) {
    SESAME->table_log_rho[i_rho] +=
        logf(units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY) /
             units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  }

  // Temperatures (log)
  for (int i_T = 0; i_T < SESAME->num_T; i_T++) {
    SESAME->table_log_T[i_T] +=
        logf(units_cgs_conversion_factor(&si, UNIT_CONV_TEMPERATURE) /
             units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
  }

  // Sp. int. energies (log), pressures, sound speeds, and sp. entropies
  for (int i_rho = 0; i_rho < SESAME->num_rho; i_rho++) {
    for (int i_T = 0; i_T < SESAME->num_T; i_T++) {
      SESAME->table_log_u_rho_T[i_rho * SESAME->num_T + i_T] += logf(
          units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
      SESAME->table_P_rho_T[i_rho * SESAME->num_T + i_T] *=
          units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
          units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
      SESAME->table_c_rho_T[i_rho * SESAME->num_T + i_T] *=
          units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
          units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
      SESAME->table_log_s_rho_T[i_rho * SESAME->num_T + i_T] +=
          logf(units_cgs_conversion_factor(
                   &si, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS) /
               units_cgs_conversion_factor(
                   us, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS));
    }
  }

  // Tiny values
  SESAME->u_tiny *=
      units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  SESAME->P_tiny *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
                    units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  SESAME->c_tiny *= units_cgs_conversion_factor(&si, UNIT_CONV_SPEED) /
                    units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
  SESAME->s_tiny *=
      units_cgs_conversion_factor(&si,
                                  UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS);
}

// gas_internal_energy_from_entropy
INLINE static float SESAME_internal_energy_from_entropy(
    float density, float entropy, const struct SESAME_params *SESAME) {

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
  idx_rho = find_value_in_monot_incr_array(log_rho, SESAME->table_log_rho,
                                           SESAME->num_rho);

  // Sp. entropy at this and the next density (in relevant slice of s array)
  idx_s_1 = find_value_in_monot_incr_array(
      log_s, SESAME->table_log_s_rho_T + idx_rho * SESAME->num_T,
      SESAME->num_T);
  idx_s_2 = find_value_in_monot_incr_array(
      log_s, SESAME->table_log_s_rho_T + (idx_rho + 1) * SESAME->num_T,
      SESAME->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= SESAME->num_rho) {
    idx_rho = SESAME->num_rho - 2;
  }
  if (idx_s_1 <= -1) {
    idx_s_1 = 0;
  } else if (idx_s_1 >= SESAME->num_T) {
    idx_s_1 = SESAME->num_T - 2;
  }
  if (idx_s_2 <= -1) {
    idx_s_2 = 0;
  } else if (idx_s_2 >= SESAME->num_T) {
    idx_s_2 = SESAME->num_T - 2;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (SESAME->table_log_rho[idx_rho + 1] != SESAME->table_log_rho[idx_rho]) {
    intp_rho =
        (log_rho - SESAME->table_log_rho[idx_rho]) /
        (SESAME->table_log_rho[idx_rho + 1] - SESAME->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.f;
  }
  if (SESAME->table_log_s_rho_T[idx_rho * SESAME->num_T + (idx_s_1 + 1)] !=
      SESAME->table_log_s_rho_T[idx_rho * SESAME->num_T + idx_s_1]) {
    intp_s_1 =
        (log_s - SESAME->table_log_s_rho_T[idx_rho * SESAME->num_T + idx_s_1]) /
        (SESAME->table_log_s_rho_T[idx_rho * SESAME->num_T + (idx_s_1 + 1)] -
         SESAME->table_log_s_rho_T[idx_rho * SESAME->num_T + idx_s_1]);
  } else {
    intp_s_1 = 1.f;
  }
  if (SESAME
          ->table_log_s_rho_T[(idx_rho + 1) * SESAME->num_T + (idx_s_2 + 1)] !=
      SESAME->table_log_s_rho_T[(idx_rho + 1) * SESAME->num_T + idx_s_2]) {
    intp_s_2 =
        (log_s -
         SESAME->table_log_s_rho_T[(idx_rho + 1) * SESAME->num_T + idx_s_2]) /
        (SESAME->table_log_s_rho_T[(idx_rho + 1) * SESAME->num_T +
                                   (idx_s_2 + 1)] -
         SESAME->table_log_s_rho_T[(idx_rho + 1) * SESAME->num_T + idx_s_2]);
  } else {
    intp_s_2 = 1.f;
  }

  // Table values
  log_u_1 = SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_s_1];
  log_u_2 = SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_s_1 + 1];
  log_u_3 = SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_s_2];
  log_u_4 =
      SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_s_2 + 1];

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
    float density, float entropy, const struct SESAME_params *SESAME) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float SESAME_entropy_from_pressure(
    float density, float pressure, const struct SESAME_params *SESAME) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float SESAME_soundspeed_from_entropy(
    float density, float entropy, const struct SESAME_params *SESAME) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float SESAME_entropy_from_internal_energy(
    float density, float u, const struct SESAME_params *SESAME) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float SESAME_pressure_from_internal_energy(
    float density, float u, const struct SESAME_params *SESAME) {

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
  idx_rho = find_value_in_monot_incr_array(log_rho, SESAME->table_log_rho,
                                           SESAME->num_rho);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  idx_u_1 = find_value_in_monot_incr_array(
      log_u, SESAME->table_log_u_rho_T + idx_rho * SESAME->num_T,
      SESAME->num_T);
  idx_u_2 = find_value_in_monot_incr_array(
      log_u, SESAME->table_log_u_rho_T + (idx_rho + 1) * SESAME->num_T,
      SESAME->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= SESAME->num_rho) {
    idx_rho = SESAME->num_rho - 2;
  }
  if (idx_u_1 <= -1) {
    idx_u_1 = 0;
  } else if (idx_u_1 >= SESAME->num_T) {
    idx_u_1 = SESAME->num_T - 2;
  }
  if (idx_u_2 <= -1) {
    idx_u_2 = 0;
  } else if (idx_u_2 >= SESAME->num_T) {
    idx_u_2 = SESAME->num_T - 2;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (SESAME->table_log_rho[idx_rho + 1] != SESAME->table_log_rho[idx_rho]) {
    intp_rho =
        (log_rho - SESAME->table_log_rho[idx_rho]) /
        (SESAME->table_log_rho[idx_rho + 1] - SESAME->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.f;
  }
  if (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] !=
      SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) /
        (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] -
         SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.f;
  }
  if (SESAME
          ->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + (idx_u_2 + 1)] !=
      SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u -
         SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) /
        (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                   (idx_u_2 + 1)] -
         SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.f;
  }

  // Table values
  P_1 = SESAME->table_P_rho_T[idx_rho * SESAME->num_T + idx_u_1];
  P_2 = SESAME->table_P_rho_T[idx_rho * SESAME->num_T + idx_u_1 + 1];
  P_3 = SESAME->table_P_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2];
  P_4 = SESAME->table_P_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2 + 1];

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
    if ((num_non_pos > 2) || (SESAME->P_tiny == 0.f) || (intp_rho < 0.f) ||
        (intp_u_1 < 0.f) || (intp_u_2 < 0.f)) {
      return 0.f;
    }
    if (P_1 <= 0.f) P_1 = SESAME->P_tiny;
    if (P_2 <= 0.f) P_2 = SESAME->P_tiny;
    if (P_3 <= 0.f) P_3 = SESAME->P_tiny;
    if (P_4 <= 0.f) P_4 = SESAME->P_tiny;
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
    float density, float P, const struct SESAME_params *SESAME) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float SESAME_soundspeed_from_internal_energy(
    float density, float u, const struct SESAME_params *SESAME) {

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
  idx_rho = find_value_in_monot_incr_array(log_rho, SESAME->table_log_rho,
                                           SESAME->num_rho);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  idx_u_1 = find_value_in_monot_incr_array(
      log_u, SESAME->table_log_u_rho_T + idx_rho * SESAME->num_T,
      SESAME->num_T);
  idx_u_2 = find_value_in_monot_incr_array(
      log_u, SESAME->table_log_u_rho_T + (idx_rho + 1) * SESAME->num_T,
      SESAME->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= SESAME->num_rho) {
    idx_rho = SESAME->num_rho - 2;
  }
  if (idx_u_1 <= -1) {
    idx_u_1 = 0;
  } else if (idx_u_1 >= SESAME->num_T) {
    idx_u_1 = SESAME->num_T - 2;
  }
  if (idx_u_2 <= -1) {
    idx_u_2 = 0;
  } else if (idx_u_2 >= SESAME->num_T) {
    idx_u_2 = SESAME->num_T - 2;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (SESAME->table_log_rho[idx_rho + 1] != SESAME->table_log_rho[idx_rho]) {
    intp_rho =
        (log_rho - SESAME->table_log_rho[idx_rho]) /
        (SESAME->table_log_rho[idx_rho + 1] - SESAME->table_log_rho[idx_rho]);
  } else {
    intp_rho = 1.f;
  }
  if (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] !=
      SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) /
        (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] -
         SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.f;
  }
  if (SESAME
          ->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + (idx_u_2 + 1)] !=
      SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u -
         SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) /
        (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                   (idx_u_2 + 1)] -
         SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.f;
  }

  // Table values
  c_1 = SESAME->table_c_rho_T[idx_rho * SESAME->num_T + idx_u_1];
  c_2 = SESAME->table_c_rho_T[idx_rho * SESAME->num_T + idx_u_1 + 1];
  c_3 = SESAME->table_c_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2];
  c_4 = SESAME->table_c_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2 + 1];

  // If more than two table values are non-positive then return zero
  int num_non_pos = 0;
  if (c_1 <= 0.f) num_non_pos++;
  if (c_2 <= 0.f) num_non_pos++;
  if (c_3 <= 0.f) num_non_pos++;
  if (c_4 <= 0.f) num_non_pos++;
  if (num_non_pos > 2) {
    return SESAME->c_tiny;
  }
  // If just one or two are non-positive then replace them with a tiny value
  else if (num_non_pos > 0) {
    // Unless already trying to extrapolate in which case return zero
    if ((intp_rho < 0.f) || (intp_u_1 < 0.f) || (intp_u_2 < 0.f)) {
      return SESAME->c_tiny;
    }
    if (c_1 <= 0.f) c_1 = SESAME->c_tiny;
    if (c_2 <= 0.f) c_2 = SESAME->c_tiny;
    if (c_3 <= 0.f) c_3 = SESAME->c_tiny;
    if (c_4 <= 0.f) c_4 = SESAME->c_tiny;
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
    float density, float P, const struct SESAME_params *SESAME) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_temperature_from_internal_energy
INLINE static float SESAME_temperature_from_internal_energy(
    float density, float u, const struct SESAME_params *SESAME) {

  float T, log_T, log_T_1, log_T_2, log_rho_1, log_rho_2;

  if (u <= 0.f) {
    return 0.f;
  }

  int idx_rho, idx_u_1, idx_u_2;
  float intp_u_1, intp_u_2;
  float slope, intercept;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u)) to find T(rho, u)
  // Density index
  idx_rho = find_value_in_monot_incr_array(log_rho, SESAME->table_log_rho,
                                           SESAME->num_rho);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  idx_u_1 = find_value_in_monot_incr_array(
      log_u, SESAME->table_log_u_rho_T + idx_rho * SESAME->num_T,
      SESAME->num_T);
  idx_u_2 = find_value_in_monot_incr_array(
      log_u, SESAME->table_log_u_rho_T + (idx_rho + 1) * SESAME->num_T,
      SESAME->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= SESAME->num_rho) {
    idx_rho = SESAME->num_rho - 2;
  }
  if (idx_u_1 <= -1) {
    idx_u_1 = 0;
  } else if (idx_u_1 >= SESAME->num_T) {
    idx_u_1 = SESAME->num_T - 2;
  }
  if (idx_u_2 <= -1) {
    idx_u_2 = 0;
  } else if (idx_u_2 >= SESAME->num_T) {
    idx_u_2 = SESAME->num_T - 2;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] !=
      SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) {
    intp_u_1 =
        (log_u - SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) /
        (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] -
         SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]);
  } else {
    intp_u_1 = 1.f;
  }
  if (SESAME
          ->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + (idx_u_2 + 1)] !=
      SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) {
    intp_u_2 =
        (log_u -
         SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) /
        (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                   (idx_u_2 + 1)] -
         SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]);
  } else {
    intp_u_2 = 1.f;
  }

  // Compute line points
  log_rho_1 = SESAME->table_log_rho[idx_rho];
  log_rho_2 = SESAME->table_log_rho[idx_rho + 1];
  log_T_1 = SESAME->table_log_T[idx_u_1];
  log_T_1 += intp_u_1 *
             (SESAME->table_log_T[idx_u_1 + 1] - SESAME->table_log_T[idx_u_1]);
  log_T_2 = SESAME->table_log_T[idx_u_2];
  log_T_2 += intp_u_2 *
             (SESAME->table_log_T[idx_u_2 + 1] - SESAME->table_log_T[idx_u_2]);

  // Intersect line passing through (log_rho_1, log_T_1), (log_rho_2, log_T_2)
  // with line density = log_rho

  // Check for log_T_1 == log_T_2
  if (log_T_1 == log_T_2) {
    log_T = log_T_1;
  } else {
    // log_rho = slope*log_T + intercept
    slope = (log_rho_1 - log_rho_2) / (log_T_1 - log_T_2);
    intercept = log_rho_1 - slope * log_T_1;
    log_T = (log_rho - intercept) / slope;
  }

  // Convert back from log
  T = expf(log_T);

  return T;
}

// gas_density_from_pressure_and_temperature
INLINE static float SESAME_density_from_pressure_and_temperature(
    float P, float T, const struct SESAME_params *SESAME) {

  float rho, log_rho, log_T_1, log_T_2, log_rho_1, log_rho_2;

  if (P <= 0.f) {
    return 0.f;
  }

  int idx_T, idx_P_1, idx_P_2;
  float intp_P_1, intp_P_2;
  float slope, intercept;
  const float log_T = logf(T);

  // 2D interpolation (bilinear with log(T), P to find rho(T, P)
  // Temperature index
  idx_T =
      find_value_in_monot_incr_array(log_T, SESAME->table_log_T, SESAME->num_T);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_T <= -1) {
    idx_T = 0;
  } else if (idx_T >= SESAME->num_T) {
    idx_T = SESAME->num_T - 2;
  }

  // Pressure at this and the next temperature (in relevant vertical slice of P
  // array)
  idx_P_1 = vertical_find_value_in_monot_incr_array(
      P, SESAME->table_P_rho_T, SESAME->num_rho, SESAME->num_T, idx_T);
  idx_P_2 = vertical_find_value_in_monot_incr_array(
      P, SESAME->table_P_rho_T, SESAME->num_rho, SESAME->num_T, idx_T + 1);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_P_1 <= -1) {
    idx_P_1 = 0;
  } else if (idx_P_1 >= SESAME->num_rho) {
    idx_P_1 = SESAME->num_rho - 2;
  }
  if (idx_P_2 <= -1) {
    idx_P_2 = 0;
  } else if (idx_P_2 >= SESAME->num_rho) {
    idx_P_2 = SESAME->num_rho - 2;
  }

  // Check for duplicates in SESAME tables before interpolation
  if (SESAME->table_P_rho_T[(idx_P_1 + 1) * SESAME->num_T + idx_T] !=
      SESAME->table_P_rho_T[idx_P_1 * SESAME->num_T + idx_T]) {
    intp_P_1 = (P - SESAME->table_P_rho_T[idx_P_1 * SESAME->num_T + idx_T]) /
               (SESAME->table_P_rho_T[(idx_P_1 + 1) * SESAME->num_T + idx_T] -
                SESAME->table_P_rho_T[idx_P_1 * SESAME->num_T + idx_T]);
  } else {
    intp_P_1 = 1.f;
  }
  if (SESAME->table_P_rho_T[(idx_P_2 + 1) * SESAME->num_T + (idx_T + 1)] !=
      SESAME->table_P_rho_T[idx_P_2 * SESAME->num_T + (idx_T + 1)]) {
    intp_P_2 =
        (P - SESAME->table_P_rho_T[idx_P_2 * SESAME->num_T + (idx_T + 1)]) /
        (SESAME->table_P_rho_T[(idx_P_2 + 1) * SESAME->num_T + (idx_T + 1)] -
         SESAME->table_P_rho_T[idx_P_2 * SESAME->num_T + (idx_T + 1)]);
  } else {
    intp_P_2 = 1.f;
  }

  // Compute line points
  log_T_1 = SESAME->table_log_T[idx_T];
  log_T_2 = SESAME->table_log_T[idx_T + 1];
  log_rho_1 = SESAME->table_log_rho[idx_P_1];
  log_rho_1 += intp_P_1 * (SESAME->table_log_rho[idx_P_1 + 1] -
                           SESAME->table_log_rho[idx_P_1]);
  log_rho_2 = SESAME->table_log_rho[idx_P_2];
  log_rho_2 += intp_P_2 * (SESAME->table_log_rho[idx_P_2 + 1] -
                           SESAME->table_log_rho[idx_P_2]);

  // Intersect line passing through (log_rho_1, log_T_1), (log_rho_2, log_T_2)
  // with line temperature = log_T

  // Check for log_rho_1 == log_rho_2
  if (log_rho_1 == log_rho_2) {
    log_rho = log_rho_1;
  } else {
    // log_T = slope*log_rho + intercept
    slope = (log_T_1 - log_T_2) / (log_rho_1 - log_rho_2);
    intercept = log_T_1 - slope * log_rho_1;
    log_rho = (log_T - intercept) / slope;
  }

  // Convert back from log
  rho = expf(log_rho);

  return rho;
}

// gas_density_from_pressure_and_internal_energy
INLINE static float SESAME_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    const struct SESAME_params *SESAME) {

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
      log_rho_ref, SESAME->table_log_rho, SESAME->num_rho);

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
  float closest_root = 0.f;
  float root_below;

  // Initialise pressures
  float P_above_lower, P_above_upper;
  float P_below_lower, P_below_upper;
  P_above_upper = 0.f;
  P_below_lower = 0.f;

  // Increase search range by search_factor_log_rho
  log_rho_max += search_factor_log_rho;
  idx_rho_above_max = find_value_in_monot_incr_array(
      log_rho_max, SESAME->table_log_rho, SESAME->num_rho);
  log_rho_min -= search_factor_log_rho;
  idx_rho_below_min = find_value_in_monot_incr_array(
      log_rho_min, SESAME->table_log_rho, SESAME->num_rho);

  float P_1, P_2, P_3, P_4;
  int idx_rho, idx_u_1, idx_u_2;
  float intp_rho, intp_u_1, intp_u_2;

  // When searching above/below, we are looking for where the pressure P(rho, u)
  // of the table densities changes from being less than to more than, or vice
  // versa, the desired pressure. If this is the case, there is a root between
  // these table values of rho.

  // First look for roots above rho_ref
  for (idx_rho = idx_rho_ref; idx_rho <= idx_rho_above_max; idx_rho++) {

    // This is similar to P_u_rho, but we're not interest in intp_rho,
    // but instead calculate the pressure for both intp_rho=0 and intp_rho=1

    // Sp. int. energy at this and the next density (in relevant slice of u
    // array)
    idx_u_1 = find_value_in_monot_incr_array(
        log_u, SESAME->table_log_u_rho_T + idx_rho * SESAME->num_T,
        SESAME->num_T);
    idx_u_2 = find_value_in_monot_incr_array(
        log_u, SESAME->table_log_u_rho_T + (idx_rho + 1) * SESAME->num_T,
        SESAME->num_T);

    // If outside the table then extrapolate from the edge and edge-but-one
    // values
    if (idx_rho <= -1) {
      idx_rho = 0;
    } else if (idx_rho >= SESAME->num_rho) {
      idx_rho = SESAME->num_rho - 2;
    }
    if (idx_u_1 <= -1) {
      idx_u_1 = 0;
    } else if (idx_u_1 >= SESAME->num_T) {
      idx_u_1 = SESAME->num_T - 2;
    }
    if (idx_u_2 <= -1) {
      idx_u_2 = 0;
    } else if (idx_u_2 >= SESAME->num_T) {
      idx_u_2 = SESAME->num_T - 2;
    }

    if (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] !=
        SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) {
      intp_u_1 =
          (log_u -
           SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) /
          (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] -
           SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]);
    } else {
      intp_u_1 = 1.f;
    }
    if (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                  (idx_u_2 + 1)] !=
        SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) {
      intp_u_2 =
          (log_u -
           SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) /
          (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                     (idx_u_2 + 1)] -
           SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]);
    } else {
      intp_u_2 = 1.f;
    }

    // Table values
    P_1 = SESAME->table_P_rho_T[idx_rho * SESAME->num_T + idx_u_1];
    P_2 = SESAME->table_P_rho_T[idx_rho * SESAME->num_T + idx_u_1 + 1];
    P_3 = SESAME->table_P_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2];
    P_4 = SESAME->table_P_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2 + 1];

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
      if ((num_non_pos > 2) || (SESAME->P_tiny == 0.f) || (intp_u_1 < 0.f) ||
          (intp_u_2 < 0.f)) {
        break;  // return rho_sph;
      }
      if (P_1 <= 0.f) P_1 = SESAME->P_tiny;
      if (P_2 <= 0.f) P_2 = SESAME->P_tiny;
      if (P_3 <= 0.f) P_3 = SESAME->P_tiny;
      if (P_4 <= 0.f) P_4 = SESAME->P_tiny;
    }

    // Interpolate with the log values
    P_1 = logf(P_1);
    P_2 = logf(P_2);
    P_3 = logf(P_3);
    P_4 = logf(P_4);

    // Pressure for intp_rho = 0
    P_above_lower = expf(((1.f - intp_u_1) * P_1 + intp_u_1 * P_2));

    // Because of linear interpolation, pressures are not exactly continuous
    // as we go from one side of a grid point to another. See if there is
    // a root between the last P_above_upper and the new P_above_lower,
    // which are approx the same.
    if (idx_rho != idx_rho_ref) {
      if ((P_above_lower - P) * (P_above_upper - P) <= 0) {
        closest_root = expf(SESAME->table_log_rho[idx_rho]);
        break;
      }
    }

    // Pressure for intp_rho = 1
    P_above_upper = expf(((1.f - intp_u_2) * P_3 + intp_u_2 * P_4));

    // Does the pressure of the adjacent table densities switch from being
    // above to below the desired pressure, or vice versa? If so, there is a
    // root.
    if ((P_above_lower - P) * (P_above_upper - P) <= 0.f) {

      // If there is a root, interpolate between the table values:
      intp_rho = (log_P - ((1 - intp_u_1) * P_1 + intp_u_1 * P_2)) /
                 (((1 - intp_u_2) * P_3 + intp_u_2 * P_4) -
                  ((1 - intp_u_1) * P_1 + intp_u_1 * P_2));

      closest_root = expf(SESAME->table_log_rho[idx_rho] +
                          intp_rho * (SESAME->table_log_rho[idx_rho + 1] -
                                      SESAME->table_log_rho[idx_rho]));

      // If the root is between the same table values as the reference value,
      // then this is the closest root, so we can return it without further
      // searching
      if (idx_rho == idx_rho_ref) {
        return closest_root;
      }

      break;
    }
  }

  // if we found a root above, change search range below so that we're only
  // looking for closer (in log) roots than the one we found
  if (closest_root) {
    log_rho_min = log_rho_ref - (logf(closest_root) - log_rho_ref);
    idx_rho_below_min = find_value_in_monot_incr_array(
        log_rho_min, SESAME->table_log_rho, SESAME->num_rho);
  }

  // Now look for roots below rho_ref
  for (idx_rho = idx_rho_ref; idx_rho >= idx_rho_below_min; idx_rho--) {

    // Sp. int. energy at this and the next density (in relevant slice of u
    // array)
    idx_u_1 = find_value_in_monot_incr_array(
        log_u, SESAME->table_log_u_rho_T + idx_rho * SESAME->num_T,
        SESAME->num_T);
    idx_u_2 = find_value_in_monot_incr_array(
        log_u, SESAME->table_log_u_rho_T + (idx_rho + 1) * SESAME->num_T,
        SESAME->num_T);

    // If outside the table then extrapolate from the edge and edge-but-one
    // values
    if (idx_rho <= -1) {
      idx_rho = 0;
    } else if (idx_rho >= SESAME->num_rho) {
      idx_rho = SESAME->num_rho - 2;
    }
    if (idx_u_1 <= -1) {
      idx_u_1 = 0;
    } else if (idx_u_1 >= SESAME->num_T) {
      idx_u_1 = SESAME->num_T - 2;
    }
    if (idx_u_2 <= -1) {
      idx_u_2 = 0;
    } else if (idx_u_2 >= SESAME->num_T) {
      idx_u_2 = SESAME->num_T - 2;
    }

    if (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] !=
        SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) {
      intp_u_1 =
          (log_u -
           SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]) /
          (SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + (idx_u_1 + 1)] -
           SESAME->table_log_u_rho_T[idx_rho * SESAME->num_T + idx_u_1]);
    } else {
      intp_u_1 = 1.f;
    }
    if (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                  (idx_u_2 + 1)] !=
        SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) {
      intp_u_2 =
          (log_u -
           SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]) /
          (SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T +
                                     (idx_u_2 + 1)] -
           SESAME->table_log_u_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2]);
    } else {
      intp_u_2 = 1.f;
    }

    // Table values
    P_1 = SESAME->table_P_rho_T[idx_rho * SESAME->num_T + idx_u_1];
    P_2 = SESAME->table_P_rho_T[idx_rho * SESAME->num_T + idx_u_1 + 1];
    P_3 = SESAME->table_P_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2];
    P_4 = SESAME->table_P_rho_T[(idx_rho + 1) * SESAME->num_T + idx_u_2 + 1];

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
      if ((num_non_pos > 2) || (SESAME->P_tiny == 0.f) || (intp_u_1 < 0.f) ||
          (intp_u_2 < 0.f)) {
        break;  // return rho_sph;
      }
      if (P_1 <= 0.f) P_1 = SESAME->P_tiny;
      if (P_2 <= 0.f) P_2 = SESAME->P_tiny;
      if (P_3 <= 0.f) P_3 = SESAME->P_tiny;
      if (P_4 <= 0.f) P_4 = SESAME->P_tiny;
    }

    // Interpolate with the log values
    P_1 = logf(P_1);
    P_2 = logf(P_2);
    P_3 = logf(P_3);
    P_4 = logf(P_4);

    // Pressure for intp_rho = 1
    P_below_upper = expf(((1.f - intp_u_2) * P_3 + intp_u_2 * P_4));
    // Because of linear interpolation, pressures are not exactly continuous
    // as we go from one side of a grid point to another. See if there is
    // a root between the last P_below_lower and the new P_below_upper,
    // which are approx the same.
    if (idx_rho != idx_rho_ref) {
      if ((P_below_lower - P) * (P_below_upper - P) <= 0) {
        closest_root = expf(SESAME->table_log_rho[idx_rho + 1]);
        break;
      }
    }
    // Pressure for intp_rho = 0
    P_below_lower = expf(((1.f - intp_u_1) * P_1 + intp_u_1 * P_2));

    // Does the pressure of the adjacent table densities switch from being
    // above to below the desired pressure, or vice versa? If so, there is a
    // root.
    if ((P_below_lower - P) * (P_below_upper - P) <= 0.f) {

      // If there is a root, interpolate between the table values:
      intp_rho = (log_P - ((1 - intp_u_1) * P_1 + intp_u_1 * P_2)) /
                 (((1 - intp_u_2) * P_3 + intp_u_2 * P_4) -
                  ((1 - intp_u_1) * P_1 + intp_u_1 * P_2));

      root_below = expf(SESAME->table_log_rho[idx_rho] +
                        intp_rho * (SESAME->table_log_rho[idx_rho + 1] -
                                    SESAME->table_log_rho[idx_rho]));

      // If we found a root above, which one is closer to the reference rho?
      if (closest_root) {
        if (fabs(logf(root_below) - logf(rho_ref)) <
            fabs(logf(closest_root) - logf(rho_ref))) {
          closest_root = root_below;
        }
      } else {
        closest_root = root_below;
      }
      break;
    }
  }

  // Return the root if we found one
  if (closest_root) {
    return closest_root;
  }

  // If we don't find a root before we reach max_counter, return rho_ref. Maybe
  // we should give an error here?
  return rho_sph;
}

// material_phase_state_from_internal_energy
INLINE static float SESAME_phase_state_from_internal_energy(
    float density, float u, const struct mat_params *SESAME,
    const struct SESAME_params *SESAME_eos) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

#endif /* SWIFT_SESAME_EQUATION_OF_STATE_H */