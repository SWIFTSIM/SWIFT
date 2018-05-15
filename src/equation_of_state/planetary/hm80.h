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

/* ------------------------------------------------------------------------- */

// Hubbard & MacFarlane (1980) parameters
struct HM80_params {
  int mat_id;
  int num_rho, num_u;
  float log_rho_min, log_rho_max, log_rho_step, inv_log_rho_step, log_u_min,
      log_u_max, log_u_step, inv_log_u_step, bulk_mod;
  float **table_P_rho_u;
};
// Table file names
/// to be read in from the parameter file instead once finished testing...
#define HM80_HHe_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_HHe.txt"
#define HM80_ice_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_ice.txt"
#define HM80_rock_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_roc.txt"

// Parameter values for each material (cgs units)
INLINE static void set_HM80_HHe(struct HM80_params *mat, int mat_id) {
  mat->mat_id = mat_id;
  mat->num_rho = 100;
  mat->num_u = 100;
  mat->log_rho_min = -9.2103404;
  mat->log_rho_max = 1.6094379;
  mat->log_rho_step = 0.1092907;
  mat->log_u_min = 9.2103404;
  mat->log_u_max = 22.3327037;
  mat->log_u_step = 0.1325491;
  mat->bulk_mod = 0;

  mat->inv_log_rho_step = 1.f / mat->log_rho_step;
  mat->inv_log_u_step = 1.f / mat->log_u_step;
}
INLINE static void set_HM80_ice(struct HM80_params *mat, int mat_id) {
  mat->mat_id = mat_id;
  mat->num_rho = 200;
  mat->num_u = 200;
  mat->log_rho_min = -6.9077553;
  mat->log_rho_max = 2.7080502;
  mat->log_rho_step = 0.0483206;
  mat->log_u_min = 6.9077553;
  mat->log_u_max = 22.3327037;
  mat->log_u_step = 0.0775123;
  mat->bulk_mod = 2.0e10;

  mat->inv_log_rho_step = 1.f / mat->log_rho_step;
  mat->inv_log_u_step = 1.f / mat->log_u_step;
}
INLINE static void set_HM80_rock(struct HM80_params *mat, int mat_id) {
  mat->mat_id = mat_id;
  mat->num_rho = 100;
  mat->num_u = 100;
  mat->log_rho_min = -6.9077553;
  mat->log_rho_max = 2.9957323;
  mat->log_rho_step = 0.1000352;
  mat->log_u_min = 9.2103404;
  mat->log_u_max = 20.7232658;
  mat->log_u_step = 0.1162922;
  mat->bulk_mod = 3.49e11;

  mat->inv_log_rho_step = 1.f / mat->log_rho_step;
  mat->inv_log_u_step = 1.f / mat->log_u_step;
}

// Read the table from file
INLINE static void load_HM80_table(struct HM80_params *mat, char *table_file) {
  // Allocate table memory
  mat->table_P_rho_u = (float **)malloc(mat->num_rho * sizeof(float *));
  for (int i = 0; i < mat->num_rho; i++) {
    mat->table_P_rho_u[i] = (float *)malloc(mat->num_u * sizeof(float));
  }

  // Load table contents from file
  FILE *f = fopen(table_file, "r");
  int c;
  for (int i = 0; i < mat->num_rho; i++) {
    for (int j = 0; j < mat->num_u; j++) {
      c = fscanf(f, "%f", &mat->table_P_rho_u[i][j]);
      if (c != 1) {
        error("Failed to read EOS table");
      }
    }
  }
  fclose(f);
}

// Convert from cgs to internal units
INLINE static void convert_units_HM80(struct HM80_params *mat,
                                      const struct unit_system *us) {
  const float Mbar_to_Ba = 1e12;    // Convert Megabar to Barye
  const float J_kg_to_erg_g = 1e4;  // Convert J/kg to erg/g

  // Table densities in cgs
  mat->log_rho_min -= logf(units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  mat->log_rho_max -= logf(units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));

  // Table energies in SI
  mat->log_u_min +=
      logf(J_kg_to_erg_g /
           units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  mat->log_u_max +=
      logf(J_kg_to_erg_g /
           units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));

  // Table Pressures in Mbar
  for (int i = 0; i < mat->num_rho; i++) {
    for (int j = 0; j < mat->num_u; j++) {
      mat->table_P_rho_u[i][j] *=
          Mbar_to_Ba / units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    }
  }

  mat->bulk_mod /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
}

// gas_internal_energy_from_entropy
INLINE static float HM80_internal_energy_from_entropy(float density,
                                                      float entropy,
                                                      struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_pressure_from_entropy
INLINE static float HM80_pressure_from_entropy(float density, float entropy,
                                               struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_pressure
INLINE static float HM80_entropy_from_pressure(float density, float pressure,
                                               struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_entropy
INLINE static float HM80_soundspeed_from_entropy(float density, float entropy,
                                                 struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_internal_energy
INLINE static float HM80_entropy_from_internal_energy(float density, float u,
                                                      struct HM80_params *mat) {

  return 0;
}

// gas_pressure_from_internal_energy
INLINE static float HM80_pressure_from_internal_energy(
    float density, float u, struct HM80_params *mat) {

  float P;

  if (u <= 0) {
    return 0;
  }

  int rho_idx, u_idx;
  float intp_rho, intp_u;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (linear in log(rho), log(u)) to find P(rho, u)
  rho_idx = floorf((log_rho - mat->log_rho_min) * mat->inv_log_rho_step);
  u_idx = floorf((log_u - mat->log_u_min) * mat->inv_log_u_step);

  intp_rho = (log_rho - mat->log_rho_min - rho_idx * mat->log_rho_step) *
             mat->inv_log_rho_step;
  intp_u =
      (log_u - mat->log_u_min - u_idx * mat->log_u_step) * mat->inv_log_u_step;

  // Return zero pressure if below the table minimum/a
  // Extrapolate the pressure for low densities
  if (rho_idx < 0) {  // Too-low rho
    P = expf(logf((1 - intp_u) * mat->table_P_rho_u[0][u_idx] +
                  intp_u * mat->table_P_rho_u[0][u_idx + 1]) +
             log_rho - mat->log_rho_min);
    if (u_idx < 0) {  // and too-low u
      P = 0;
    }
  } else if (u_idx < 0) {  // Too-low u
    P = 0;
  }
  // Return an edge value if above the table maximum/a
  else if (rho_idx >= mat->num_rho - 1) {  // Too-high rho
    if (u_idx >= mat->num_u - 1) {         // and too-high u
      P = mat->table_P_rho_u[mat->num_rho - 1][mat->num_u - 1];
    } else {
      P = mat->table_P_rho_u[mat->num_rho - 1][u_idx];
    }
  } else if (u_idx >= mat->num_u - 1) {  // Too-high u
    P = mat->table_P_rho_u[rho_idx][mat->num_u - 1];
  }
  // Normal interpolation within the table
  else {
    P = (1 - intp_rho) * ((1 - intp_u) * mat->table_P_rho_u[rho_idx][u_idx] +
                          intp_u * mat->table_P_rho_u[rho_idx][u_idx + 1]) +
        intp_rho * ((1 - intp_u) * mat->table_P_rho_u[rho_idx + 1][u_idx] +
                    intp_u * mat->table_P_rho_u[rho_idx + 1][u_idx + 1]);
  }

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float HM80_internal_energy_from_pressure(
    float density, float P, struct HM80_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_internal_energy
INLINE static float HM80_soundspeed_from_internal_energy(
    float density, float u, struct HM80_params *mat) {

  float c, P;

  // Bulk modulus
  if (mat->bulk_mod != 0) {
    c = sqrtf(mat->bulk_mod / density);
  }
  // Ideal gas
  else {
    P = HM80_pressure_from_internal_energy(density, u, mat);
    c = sqrtf(hydro_gamma * P / density);
  }

  return c;
}

// gas_soundspeed_from_pressure
INLINE static float HM80_soundspeed_from_pressure(float density, float P,
                                                  struct HM80_params *mat) {

  float c;

  // Bulk modulus
  if (mat->bulk_mod != 0) {
    c = sqrtf(mat->bulk_mod / density);
  }
  // Ideal gas
  else {
    c = sqrtf(hydro_gamma * P / density);
  }

  return c;
}

#endif /* SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H */
