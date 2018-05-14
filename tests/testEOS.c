/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 * So far only has the Hubbard & MacFarlane (1980) equations of state.
 *
 * Sys args:
 *      mat_id | int | Material ID, see equation_of_state.h for the options.
 *          Default: 21 (= HM80_ice).
 *
 *      do_output | int | Set 1 to write the output file of rho, u, P values,
 *          set 0 for no output. Default: 0.
 *
 * Output text file contains:
 *  header
 *  num_rho  num_u  mat_id                      # Header info
 *  rho_0   rho_1   rho_2   ...   rho_num_rho   # Array of densities, rho
 *  u_0     u_1     u_2     ...   u_num_u       # Array of energies, u
 *  P_0_0   P_0_1   ...     P_0_num_u           # Array of pressures, P(rho, u)
 *  P_1_0   ...     ...     P_1_num_u
 *  ...     ...     ...     ...
 *  P_num_rho_0     ...     P_num_rho_num_u
 *
 */

int main(int argc, char *argv[]) {
  int mat_id, do_output;
  struct HM80_params mat;
  float rho, log_rho, log_u, P;
  int num_rho, num_u;
  struct unit_system us;
  const struct phys_const *phys_const = 0;  // Unused placeholder
  const struct swift_params *params = 0;    // Unused placeholder

  // Check the number of system arguments and set defaults if not provided
  switch (argc) {
    case 1:
      // Default both
      mat_id = HM80_ice;
      do_output = 0;
      break;

    case 2:
      // Read mat_id, default do_output
      mat_id = atoi(argv[1]);
      do_output = 0;
      break;

    case 3:
      // Read both
      mat_id = atoi(argv[1]);
      do_output = atoi(argv[2]);
      break;

    default:
      error("Invalid number of system arguments!\n");
      mat_id = HM80_ice;  // Ignored, just here to keep the compiler happy
      do_output = 0;
  };

  /* Greeting message */
  printf("This is %s\n", package_description());

  // Select the material parameters
  switch (mat_id) {
    case HM80_HHe:
      printf("HM80_HHe \n");
      set_HM80_HHe(&mat);
      load_HM80_table(&mat, HM80_HHe_table_file);
      break;

    case HM80_ice:
      printf("HM80_ice \n");
      set_HM80_ice(&mat);
      load_HM80_table(&mat, HM80_ice_table_file);
      break;

    case HM80_rock:
      printf("HM80_rock \n");
      set_HM80_rock(&mat);
      load_HM80_table(&mat, HM80_rock_table_file);
      break;

    default:
      error("Unknown material ID! mat_id = %d", mat_id);
      set_HM80_rock(&mat);  // Ignored, just here to keep the compiler happy
      load_HM80_table(&mat, HM80_rock_table_file);
  };

  // Convert to internal units
  units_init(&us, 5.9724e27, 6.3710e8, 1, 1, 1);
  convert_units_HM80(&mat, &us);

  eos_init(&eos, phys_const, &us, params);

  // Output file
  FILE *f = fopen("testEOS_rho_u_P.txt", "w");
  if (f == NULL) {
    printf("Could not open output file!\n");
    exit(EXIT_FAILURE);
  }

  num_rho = (mat.log_rho_max - mat.log_rho_min) / mat.log_rho_step;
  num_u = (mat.log_u_max - mat.log_u_min) / mat.log_u_step;
  if (do_output == 1) {
    fprintf(f, "Density  Sp.Int.Energy  mat_id \n");
    fprintf(f, "%d      %d            %d \n", num_rho, num_u, mat_id);
  }

  // Arrays of densities and energies
  float A1_rho[num_rho], A1_u[num_u];

  log_rho = mat.log_rho_min;
  for (int i = 0; i < num_rho; i++) {
    A1_rho[i] = exp(log_rho);
    log_rho += mat.log_rho_step;

    if (do_output == 1)
      fprintf(f, "%.6e ",
              A1_rho[i] * units_cgs_conversion_factor(&us, UNIT_CONV_DENSITY));
  }

  if (do_output == 1) fprintf(f, "\n");
  log_u = mat.log_u_min;
  for (int i = 0; i < num_u; i++) {
    A1_u[i] = exp(log_u);
    log_u += mat.log_u_step;

    if (do_output == 1)
      fprintf(f, "%.6e ", A1_u[i] * units_cgs_conversion_factor(
                                        &us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
  }

  // Pressures P(rho, u)
  if (do_output == 1) fprintf(f, "\n");
  for (int i = 0; i < num_rho; i++) {
    rho = A1_rho[i];

    for (int j = 0; j < num_u; j++) {
      P = gas_pressure_from_internal_energy(rho, A1_u[j], mat.mat_id);

      if (do_output == 1)
        fprintf(f, "%.6e ",
                P * units_cgs_conversion_factor(&us, UNIT_CONV_PRESSURE));
    }

    if (do_output == 1) fprintf(f, "\n");
  }
  fclose(f);

  return 0;
}
