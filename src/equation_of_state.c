/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* This object's header. */
#include "equation_of_state.h"

/* local headers */
#include "common_io.h"

/* Equation of state for the physics model
 * (temporary ugly solution as a global variable) */
struct eos_parameters eos;

void eos_init(struct eos_parameters *e, const struct swift_params *params) {

#if defined(EOS_IDEAL_GAS)
/* nothing to do here */
#elif defined(EOS_ISOTHERMAL_GAS)
  e->isothermal_internal_energy =
      parser_get_param_float(params, "EoS:isothermal_internal_energy");
#endif
}

void eos_print(const struct eos_parameters *e) {

#if defined(EOS_IDEAL_GAS)
  message("Equation of state: Ideal gas.");
#elif defined(EOS_ISOTHERMAL_GAS)
  message(
      "Equation of state: Isothermal with internal energy "
      "per unit mass set to %f.",
      e->isothermal_internal_energy);
#endif

  message("Adiabatic index gamma: %f.", hydro_gamma);
}

#if defined(HAVE_HDF5)
void eos_print_snapshot(hid_t h_grpsph, const struct eos_parameters *e) {

  io_write_attribute_f(h_grpsph, "Adiabatic index", hydro_gamma);

#if defined(EOS_IDEAL_GAS)
  io_write_attribute_s(h_grpsph, "Equation of state", "Ideal gas");
#elif defined(EOS_ISOTHERMAL_GAS)
  io_write_attribute_s(h_grpsph, "Equation of state", "Isothermal gas");
  io_write_attribute_f(h_grpsph, "Thermal energy per unit mass",
                       e->isothermal_internal_energy);
#endif
}
#endif
