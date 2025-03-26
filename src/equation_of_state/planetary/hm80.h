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
#ifndef SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H
#define SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/hm80.h
 *
 * Contains the Hubbard & MacFarlane (1980) Uranus/Neptune EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 *
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "units.h"

// Hubbard & MacFarlane (1980) parameters
struct HM80_params {
  float *table_log_P_rho_u;
  int version_date, num_rho, num_u;
  float log_rho_min, log_rho_max, log_rho_step, inv_log_rho_step, log_u_min,
      log_u_max, log_u_step, inv_log_u_step, bulk_mod, P_min_for_c_min;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material (SI units)
INLINE static void set_HM80_HHe(struct HM80_params *mat,
                                enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->bulk_mod = 0.f;
  mat->P_min_for_c_min = 1e3f;
  mat->version_date = 20230710;
}
INLINE static void set_HM80_ice(struct HM80_params *mat,
                                enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->bulk_mod = 2.0e9f;
  mat->P_min_for_c_min = 0.f;
  mat->version_date = 20230710;
}
INLINE static void set_HM80_rock(struct HM80_params *mat,
                                 enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->bulk_mod = 3.49e10f;
  mat->P_min_for_c_min = 0.f;
  mat->version_date = 20230710;
}

// Read the table from file
INLINE static void load_table_HM80(struct HM80_params *mat, char *table_file) {

  /* File contents:
  header (11 lines)
  version_date
  log_rho_min  log_rho_max  num_rho  log_u_min  log_u_max  num_u  (SI)
  P_0_0   P_0_1   ...     P_0_num_u           # Array of pressures (Pa)
  P_1_0   ...     ...     P_1_num_u
  ...     ...     ...     ...
  P_num_rho_0     ...     P_num_rho_num_u
  T_0_0   T_0_1   ...     T_0_num_u           # Array of temperatures (K)
  T_1_0   ...     ...     T_1_num_u
  ...     ...     ...     ...
  T_num_rho_0     ...     T_num_rho_num_u
  */

  // Load table contents from file
  FILE *f = fopen(table_file, "r");
  if (f == NULL) error("Failed to open the HM80 EoS file '%s'", table_file);

  // Ignore header lines
  char buffer[100];
  for (int i = 0; i < 11; i++) {
    if (fgets(buffer, 100, f) == NULL)
      error("Failed to read the HM80 EoS file header %s", table_file);
  }

  // Table properties
  int version_date;
  int c = fscanf(f, "%d", &version_date);
  if (c != 1) error("Failed to read the HM80 EoS table %s", table_file);
  if (version_date != mat->version_date)
    error(
        "EoS file %s version_date %d does not match expected %d"
        "\nPlease download the file using "
        "examples/Planetary/EoSTables/get_eos_tables.sh",
        table_file, version_date, mat->version_date);
  c = fscanf(f, "%f %f %d %f %f %d", &mat->log_rho_min, &mat->log_rho_max,
             &mat->num_rho, &mat->log_u_min, &mat->log_u_max, &mat->num_u);
  if (c != 6) error("Failed to read the HM80 EoS table %s", table_file);
  mat->log_rho_step =
      (mat->log_rho_max - mat->log_rho_min) / (mat->num_rho - 1);
  mat->log_u_step = (mat->log_u_max - mat->log_u_min) / (mat->num_u - 1);
  mat->inv_log_rho_step = 1.f / mat->log_rho_step;
  mat->inv_log_u_step = 1.f / mat->log_u_step;

  // Allocate table memory
  mat->table_log_P_rho_u =
      (float *)malloc(mat->num_rho * mat->num_u * sizeof(float));

  // Pressures (not log yet)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_u = 0; i_u < mat->num_u; i_u++) {
      c = fscanf(f, "%f", &mat->table_log_P_rho_u[i_rho * mat->num_u + i_u]);
      if (c != 1) error("Failed to read the HM80 EoS table %s", table_file);
    }
  }
  fclose(f);
}

// Misc. modifications
INLINE static void prepare_table_HM80(struct HM80_params *mat) {

  // Convert pressures to log(pressure)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_u = 0; i_u < mat->num_u; i_u++) {
      mat->table_log_P_rho_u[i_rho * mat->num_u + i_u] =
          logf(mat->table_log_P_rho_u[i_rho * mat->num_u + i_u]);
    }
  }
}

// Convert to internal units
INLINE static void convert_units_HM80(struct HM80_params *mat,
                                      const struct unit_system *us) {
  struct unit_system si;
  units_init_si(&si);

  // All table values in SI
  mat->log_rho_min += logf(units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY) /
                           units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  mat->log_rho_max += logf(units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY) /
                           units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));

  mat->log_u_min +=
      logf(units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
           units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  mat->log_u_max +=
      logf(units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
           units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));

  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_u = 0; i_u < mat->num_u; i_u++) {
      mat->table_log_P_rho_u[i_rho * mat->num_u + i_u] +=
          logf(units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
               units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE));
    }
  }

  mat->bulk_mod *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
                   units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->P_min_for_c_min *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE) /
                          units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
}

// gas_internal_energy_from_entropy
INLINE static float HM80_internal_energy_from_entropy(
    float density, float entropy, const struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_pressure_from_entropy
INLINE static float HM80_pressure_from_entropy(float density, float entropy,
                                               const struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float HM80_entropy_from_pressure(float density, float pressure,
                                               const struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float HM80_soundspeed_from_entropy(
    float density, float entropy, const struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float HM80_entropy_from_internal_energy(
    float density, float u, const struct HM80_params *mat) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float HM80_pressure_from_internal_energy(
    float density, float u, const struct HM80_params *mat) {

  float log_P, log_P_1, log_P_2, log_P_3, log_P_4;

  if (u <= 0.f) {
    return 0.f;
  }

  int idx_rho, idx_u;
  float intp_rho, intp_u;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u)) to find P(rho, u)
  idx_rho = floor((log_rho - mat->log_rho_min) * mat->inv_log_rho_step);
  idx_u = floor((log_u - mat->log_u_min) * mat->inv_log_u_step);

  // If outside the table then extrapolate from the edge and edge-but-one values
  if (idx_rho <= -1) {
    idx_rho = 0;
  } else if (idx_rho >= mat->num_rho - 1) {
    idx_rho = mat->num_rho - 2;
  }
  if (idx_u <= -1) {
    idx_u = 0;
  } else if (idx_u >= mat->num_u - 1) {
    idx_u = mat->num_u - 2;
  }

  intp_rho = (log_rho - mat->log_rho_min - idx_rho * mat->log_rho_step) *
             mat->inv_log_rho_step;
  intp_u =
      (log_u - mat->log_u_min - idx_u * mat->log_u_step) * mat->inv_log_u_step;

  // Table values
  log_P_1 = mat->table_log_P_rho_u[idx_rho * mat->num_u + idx_u];
  log_P_2 = mat->table_log_P_rho_u[idx_rho * mat->num_u + idx_u + 1];
  log_P_3 = mat->table_log_P_rho_u[(idx_rho + 1) * mat->num_u + idx_u];
  log_P_4 = mat->table_log_P_rho_u[(idx_rho + 1) * mat->num_u + idx_u + 1];

  log_P = (1.f - intp_rho) * ((1.f - intp_u) * log_P_1 + intp_u * log_P_2) +
          intp_rho * ((1.f - intp_u) * log_P_3 + intp_u * log_P_4);

  return expf(log_P);
}

// gas_internal_energy_from_pressure
INLINE static float HM80_internal_energy_from_pressure(
    float density, float P, const struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float HM80_soundspeed_from_internal_energy(
    float density, float u, const struct HM80_params *mat) {

  float c, P;

  // Bulk modulus
  if (mat->bulk_mod != 0) {
    c = sqrtf(mat->bulk_mod / density);
  }
  // Ideal gas
  else {
    P = HM80_pressure_from_internal_energy(density, u, mat);
    c = sqrtf(hydro_gamma * P / density);

    if (c <= 0) {
      c = sqrtf(hydro_gamma * mat->P_min_for_c_min / density);
    }
  }

  return c;
}

// gas_soundspeed_from_pressure
INLINE static float HM80_soundspeed_from_pressure(
    float density, float P, const struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

#endif /* SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H */
