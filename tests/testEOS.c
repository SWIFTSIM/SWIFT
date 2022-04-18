/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

/* Local headers. */
#include "equation_of_state.h"
#include "swift.h"

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif

/**
 * @brief Write a list of densities, energies, and resulting pressures to file
 *  from an equation of state.
 *
 *                      WORK IN PROGRESS
 *
 * So far only does the Hubbard & MacFarlane (1980) equations of state.
 *
 * Usage:
 *      $  ./testEOS  (mat_id)  (do_output)
 *
 * Sys args (optional):
 *      mat_id | int | Material ID, see equation_of_state.h for the options.
 *          Default: 201 (= id_HM80_ice).
 *
 *      do_output | int | Set 1 to write the output file of rho, u, P values,
 *          set 0 for no output. Default: 0.
 *
 * Output text file contains:
 *  header
 *  num_rho num_u   mat_id                      # Header values
 *  rho_0   rho_1   rho_2   ...   rho_num_rho   # Array of densities, rho
 *  u_0     u_1     u_2     ...   u_num_u       # Array of energies, u
 *  P_0_0   P_0_1   ...     P_0_num_u           # Array of pressures, P(rho, u)
 *  P_1_0   ...     ...     P_1_num_u
 *  ...     ...     ...     ...
 *  P_num_rho_0     ...     P_num_rho_num_u
 *  c_0_0   c_0_1   ...     c_0_num_u           # Array of sound speeds, c(rho,
 * u)
 *  c_1_0   ...     ...     c_1_num_u
 *  ...     ...     ...     ...
 *  c_num_rho_0     ...     c_num_rho_num_u
 *
 * Note that the values tested extend beyond the range that most EOS are
 * designed for (e.g. outside table limits), to help test the EOS in case of
 * unexpected particle behaviour.
 *
 */

#ifdef EOS_PLANETARY
int main(int argc, char *argv[]) {
  float rho, u, log_rho, log_u, P, c;
  struct unit_system us;
  struct swift_params *params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  if (params == NULL) error("Error allocating memory for the parameter file.");
  const struct phys_const *phys_const = 0;  // Unused placeholder
  const float J_kg_to_erg_g = 1e4;          // Convert J/kg to erg/g
  char filename[64];
  // Output table params
  const int num_rho = 100, num_u = 100;
  float log_rho_min = logf(1e-4f), log_rho_max = logf(1e3f),  // Densities (cgs)
      log_u_min = logf(1e4f),
        log_u_max = logf(1e13f),  // Sp. int. energies (SI)
      log_rho_step = (log_rho_max - log_rho_min) / (num_rho - 1.f),
        log_u_step = (log_u_max - log_u_min) / (num_u - 1.f);
  float A1_rho[num_rho], A1_u[num_u];
  // Sys args
  int mat_id_in, do_output;
  // Default sys args
  const int mat_id_def = eos_planetary_id_HM80_ice;
  const int do_output_def = 0;

  // Check the number of system arguments and use defaults if not provided
  switch (argc) {
    case 1:
      // Default both
      mat_id_in = mat_id_def;
      do_output = do_output_def;
      break;

    case 2:
      // Read mat_id, default do_output
      mat_id_in = atoi(argv[1]);
      do_output = do_output_def;
      break;

    case 3:
      // Read both
      mat_id_in = atoi(argv[1]);
      do_output = atoi(argv[2]);
      break;

    default:
      error("Invalid number of system arguments!\n");
      mat_id_in = mat_id_def;  // Ignored, just here to keep the compiler happy
      do_output = do_output_def;
  };

  enum eos_planetary_material_id mat_id =
      (enum eos_planetary_material_id)mat_id_in;

  /* Greeting message */
  printf("This is %s\n", package_description());

  // Check material ID
  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  // Select the material base type
  switch (type) {
    // Tillotson
    case eos_planetary_type_Til:
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          printf("  Tillotson iron \n");
          break;

        case eos_planetary_id_Til_granite:
          printf("  Tillotson granite \n");
          break;

        case eos_planetary_id_Til_water:
          printf("  Tillotson water \n");
          break;

        default:
          error("Unknown material ID! mat_id = %d \n", mat_id);
      };
      break;

    // Hubbard & MacFarlane (1980)
    case eos_planetary_type_HM80:
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          printf("  Hubbard & MacFarlane (1980) hydrogen-helium atmosphere \n");
          break;

        case eos_planetary_id_HM80_ice:
          printf("  Hubbard & MacFarlane (1980) ice mix \n");
          break;

        case eos_planetary_id_HM80_rock:
          printf("  Hubbard & MacFarlane (1980) rock mix \n");
          break;

        default:
          error("Unknown material ID! mat_id = %d \n", mat_id);
      };
      break;

    // SESAME
    case eos_planetary_type_SESAME:
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          printf("  SESAME basalt 7530 \n");
          break;

        case eos_planetary_id_SESAME_basalt:
          printf("  SESAME basalt 7530 \n");
          break;

        case eos_planetary_id_SESAME_water:
          printf("  SESAME water 7154 \n");
          break;

        case eos_planetary_id_SS08_water:
          printf("  Senft & Stewart (2008) SESAME-like water \n");
          break;

        default:
          error("Unknown material ID! mat_id = %d \n", mat_id);
      };
      break;

    default:
      error("Unknown material type! mat_id = %d \n", mat_id);
  }

  // Convert to internal units
  // Earth masses and radii
  //  units_init(&us, 5.9724e27, 6.3710e8, 1.f, 1.f, 1.f);
  // SI
  units_init(&us, 1000.f, 100.f, 1.f, 1.f, 1.f);
  log_rho_min -= logf(units_cgs_conversion_factor(&us, UNIT_CONV_DENSITY));
  log_rho_max -= logf(units_cgs_conversion_factor(&us, UNIT_CONV_DENSITY));
  log_u_min += logf(J_kg_to_erg_g / units_cgs_conversion_factor(
                                        &us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  log_u_max += logf(J_kg_to_erg_g / units_cgs_conversion_factor(
                                        &us, UNIT_CONV_ENERGY_PER_UNIT_MASS));

  // Set the input parameters
  // Which EOS to initialise
  parser_set_param(params, "EoS:planetary_use_Til:1");
  parser_set_param(params, "EoS:planetary_use_HM80:1");
  parser_set_param(params, "EoS:planetary_use_SESAME:1");
  // Table file names
  parser_set_param(params,
                   "EoS:planetary_HM80_HHe_table_file:"
                   "../examples/planetary_HM80_HHe.txt");
  parser_set_param(params,
                   "EoS:planetary_HM80_ice_table_file:"
                   "../examples/planetary_HM80_ice.txt");
  parser_set_param(params,
                   "EoS:planetary_HM80_rock_table_file:"
                   "../examples/planetary_HM80_rock.txt");
  parser_set_param(params,
                   "EoS:planetary_SESAME_iron_table_file:"
                   "../examples/planetary_SESAME_iron_2140.txt");
  parser_set_param(params,
                   "EoS:planetary_SESAME_basalt_table_file:"
                   "../examples/planetary_SESAME_basalt_7530.txt");
  parser_set_param(params,
                   "EoS:planetary_SESAME_water_table_file:"
                   "../examples/planetary_SESAME_water_7154.txt");
  parser_set_param(params,
                   "EoS:planetary_SS08_water_table_file:"
                   "../examples/planetary_SS08_water.txt");

  // Initialise the EOS materials
  eos_init(&eos, phys_const, &us, params);

  // Manual debug testing
  if (1) {
    printf("\n ### MANUAL DEBUG TESTING ### \n");

    rho = 5960;
    u = 1.7e8;
    P = gas_pressure_from_internal_energy(rho, u, eos_planetary_id_HM80_ice);
    printf("u = %.2e,    rho = %.2e,    P = %.2e \n", u, rho, P);

    return 0;
  }

  // Output file
  sprintf(filename, "testEOS_rho_u_P_c_%d.txt", mat_id);
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("Could not open output file!\n");
    exit(EXIT_FAILURE);
  }

  if (do_output == 1) {
    fprintf(f, "Density  Sp.Int.Energy  mat_id \n");
    fprintf(f, "%d      %d            %d \n", num_rho, num_u, mat_id);
  }

  // Densities
  log_rho = log_rho_min;
  for (int i = 0; i < num_rho; i++) {
    A1_rho[i] = exp(log_rho);
    log_rho += log_rho_step;

    if (do_output == 1)
      fprintf(f, "%.6e ",
              A1_rho[i] * units_cgs_conversion_factor(&us, UNIT_CONV_DENSITY));
  }
  if (do_output == 1) fprintf(f, "\n");

  // Sp. int. energies
  log_u = log_u_min;
  for (int i = 0; i < num_u; i++) {
    A1_u[i] = exp(log_u);
    log_u += log_u_step;

    if (do_output == 1)
      fprintf(f, "%.6e ",
              A1_u[i] * units_cgs_conversion_factor(
                            &us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  }
  if (do_output == 1) fprintf(f, "\n");

  // Pressures
  for (int i = 0; i < num_rho; i++) {
    rho = A1_rho[i];

    for (int j = 0; j < num_u; j++) {
      P = gas_pressure_from_internal_energy(rho, A1_u[j], mat_id);

      if (do_output == 1)
        fprintf(f, "%.6e ",
                P * units_cgs_conversion_factor(&us, UNIT_CONV_PRESSURE));
    }

    if (do_output == 1) fprintf(f, "\n");
  }

  // Sound speeds
  for (int i = 0; i < num_rho; i++) {
    rho = A1_rho[i];

    for (int j = 0; j < num_u; j++) {
      c = gas_soundspeed_from_internal_energy(rho, A1_u[j], mat_id);

      if (do_output == 1)
        fprintf(f, "%.6e ",
                c * units_cgs_conversion_factor(&us, UNIT_CONV_SPEED));
    }

    if (do_output == 1) fprintf(f, "\n");
  }
  fclose(f);

  return 0;
}
#else
int main(int argc, char *argv[]) { return 0; }
#endif
