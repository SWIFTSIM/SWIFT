/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_MULTI_TAIT_EQUATION_OF_STATE_H
#define SWIFT_MULTI_TAIT_EQUATION_OF_STATE_H

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "common_io.h"
#include "inline.h"
#include "error.h"
#include "physical_constants.h"

extern struct eos_parameters eos;

/**
 * @brief The parameters of the equation of state for the gas.
 *
 * This equation of state takes a single argument, the soundspeed.
 */
struct eos_parameters {

  float soundspeed;
  float soundspeed_squared;
  float gamma;
};

__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_pressure(float density, float P) {

  return eos.soundspeed;
}


__attribute__((always_inline)) INLINE static float pressure_from_density( double density, double reference_density ){ 
   const double B = eos.soundspeed_squared * reference_density / eos.gamma;
   const double rho_over_rho0 = density / reference_density;
   const double rho_over_rho0_2 = rho_over_rho0 * rho_over_rho0;
   const double rho_over_rho0_4 = rho_over_rho0_2 * rho_over_rho0_2;
   const double rho_over_rho0_6 = rho_over_rho0_4 * rho_over_rho0_2;
   const double rho_over_rho0_7 = rho_over_rho0_6 * rho_over_rho0;

   return B * ( rho_over_rho0_7 - 1.0);
}

INLINE static void eos_init(struct eos_parameters *e,
                            const struct phys_const *phys_const,
                            const struct unit_system *us,
                            struct swift_params *params){
//  e->B = parser_get_opt_param_float(params, "EoS:B", 0.01);
  e->soundspeed = parser_get_param_float(params, "EoS:soundspeed");
  e->gamma = 7.0;
//  e->soundspeed = sqrt((double)e->gamma * (double)e->B / (double)e->density_reference);
  e->soundspeed_squared = e->soundspeed * e->soundspeed;
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
INLINE static void eos_print(const struct eos_parameters *e) {

  message("Equation of state: Multifluid Tait.");

  message("Soundspeed: %f.\n", eos.soundspeed);
}

#if defined(HAVE_HDF5)
/**
 * @brief Write equation of state information to the snapshot
 *
 * @param h_grpsph The HDF5 group in which to write
 * @param e The #eos_parameters
 */
INLINE static void eos_print_snapshot(hid_t h_grpsph,
                                      const struct eos_parameters *e) {

  io_write_attribute_s(h_grpsph, "Equation of state", "Weakly Compressible");
}
#endif

#endif /* SWIFT_IDEAL_GAS_EQUATION_OF_STATE_H */
