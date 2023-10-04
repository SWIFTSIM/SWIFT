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
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_particles(const struct part* parts,
                                    struct io_props* list) {
  return 0;
}

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_stars(const struct spart* sparts,
                                struct io_props* list) {
  return 0;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  list[0] =
      io_make_output_field("RTDebugInjectionDone", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.debug_injection_done,
                           "How many times rt_injection_update_photon_density "
                           "has been called");
  list[1] = io_make_output_field(
      "RTDebugCallsIactGradientInteractions", INT, 1, UNIT_CONV_NO_UNITS, 0,
      parts, rt_data.debug_calls_iact_gradient_interaction,
      "number of calls to this particle during the gradient interaction loop "
      "from the actual interaction function");
  list[2] = io_make_output_field("RTDebugCallsIactTransportInteractions", INT,
                                 1, UNIT_CONV_NO_UNITS, 0, parts,
                                 rt_data.debug_calls_iact_transport_interaction,
                                 "number of calls to this particle during the "
                                 "transport interaction loop from the actual "
                                 "interaction function");
  list[3] =
      io_make_output_field("RTDebugGradientsDone", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.debug_gradients_done,
                           "How many times finalise_gradients was called");
  list[4] =
      io_make_output_field("RTDebugTransportDone", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.debug_transport_done,
                           "How many times finalise_transport was called");
  list[5] = io_make_output_field(
      "RTDebugThermochemistryDone", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.debug_thermochem_done, "How many times rt_tchem was called");
  list[6] = io_make_output_field(
      "RTDebugRadAbsorbedTot", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0, parts,
      rt_data.debug_radiation_absorbed_tot,
      "Radiation absorbed by this part during its lifetime");
  list[7] = io_make_output_field("RTDebugSubcycles", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, parts, rt_data.debug_nsubcycles,
                                 "How many times this part was subcycled");

  return 8;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {

  list[0] = io_make_output_field("RTDebugHydroIact", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, sparts, rt_data.debug_iact_hydro_inject,
                                 "number of interactions between this star "
                                 "particle and any particle during injection");
  list[1] = io_make_output_field(
      "RTDebugEmissionRateSet", INT, 1, UNIT_CONV_NO_UNITS, 0, sparts,
      rt_data.debug_emission_rate_set, "Stellar photon emission rates set?");
  list[2] = io_make_output_field(
      "RTDebugRadEmittedTot", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0, sparts,
      rt_data.debug_radiation_emitted_tot,
      "Total radiation emitted during the lifetime of this star");
  return 3;
}

/**
 * @brief Write the RT model properties to the snapshot.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The engine
 * @param internal_units The internal unit system
 * @param snapshot_units Units used for the snapshot
 * @param rtp The #rt_props
 */
INLINE static void rt_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                    const struct engine* e,
                                    const struct unit_system* internal_units,
                                    const struct unit_system* snapshot_units,
                                    const struct rt_props* rtp) {
#if defined(HAVE_HDF5)

  io_write_attribute_s(h_grp, "RT Scheme", RT_IMPLEMENTATION);

#endif
}

#endif /* SWIFT_RT_IO_DEBUG_H */
