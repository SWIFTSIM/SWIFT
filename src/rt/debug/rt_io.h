/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_DEBUG_H
#define SWIFT_RT_IO_DEBUG_H

#include "io_properties.h"

/**
 * @file src/rt/debug/rt_io.h
 * @brief Main header file for the debug radiative transfer scheme IO routines.
 */

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  list[0] =
      io_make_output_field("RTStarIact", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
                           rt_data.iact_stars_inject,
                           "number of interactions between this hydro particle"
                           " and any star particle during injection step");
  list[1] = io_make_output_field("RTTotalCalls", INT, 1, UNIT_CONV_NO_UNITS, 0,
                                 parts, rt_data.calls_tot,
                                 "total number of calls to this "
                                 "particle during the run");
  list[2] = io_make_output_field("RTCallsThisStep", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, parts, rt_data.calls_per_step,
                                 "number of calls "
                                 "to this particle during one time step");
  list[3] =
      io_make_output_field("RTPhotonsUpdated", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.photon_number_updated,
                           "=1 if photon number has been updated in this step");
  list[4] =
      io_make_output_field("RTCallsIactGradient", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_iact_gradient,
                           "number of calls to this particle during the"
                           "gradient interaction loop");
  list[5] =
      io_make_output_field("RTCallsIactTransport", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_iact_transport,
                           "number of calls to this particle during the"
                           "transport interaction loop");
  list[6] = io_make_output_field(
      "RTGradientsDone", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.gradients_done, "How many times finalise_gradients was called");
  list[7] = io_make_output_field(
      "RTTransportDone", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.transport_done, "How many times finalise_transport was called");
  list[8] = io_make_output_field(
      "RTThermochemistryDone", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.thermochem_done, "How many times rt_tchem was called");

  return 9;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {

  list[0] = io_make_output_field("RTHydroIact", INT, 1, UNIT_CONV_NO_UNITS, 0,
                                 sparts, rt_data.iact_hydro_inject,
                                 "number of interactions between this hydro "
                                 "particle and any star particle");
  list[1] = io_make_output_field("RTTotalCalls", INT, 1, UNIT_CONV_NO_UNITS, 0,
                                 sparts, rt_data.calls_tot,
                                 "total number of calls "
                                 "to this particle during the run");
  list[2] = io_make_output_field("RTCallsThisStep", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, sparts, rt_data.calls_per_step,
                                 "number of calls to "
                                 "this particle during one time step");
  list[3] =
      io_make_output_field("RTEmissionRateSet", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           sparts, rt_data.emission_rate_set,
                           "Stellar photon "
                           "emission rates set?");

  return 4;
}
#endif /* SWIFT_RT_IO_DEBUG_H */
