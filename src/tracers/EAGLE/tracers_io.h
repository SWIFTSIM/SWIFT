/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2021 Edo Altamura (edoardo.altamura@manchester.ac.uk)
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
#ifndef SWIFT_TRACERS_EAGLE_IO_H
#define SWIFT_TRACERS_EAGLE_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "io_properties.h"
#include "tracers.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of tracers to the file.
 *
 * @param h_grp The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void tracers_write_flavour(
    hid_t h_grp) {

  io_write_attribute_s(h_grp, "Tracers", "EAGLE");
}
#endif

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int tracers_write_particles(
    const struct part* parts, const struct xpart* xparts, struct io_props* list,
    const int with_cosmology) {

  list[0] = io_make_output_field(
      "MaximalTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, xparts,
      tracers_data.maximum_temperature,
      "Maximal temperatures ever reached by the particles");

  if (with_cosmology) {
    list[1] = io_make_output_field(
        "MaximalTemperatureScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        xparts, tracers_data.maximum_temperature_scale_factor,
        "Scale-factors at which the maximal temperature was reached");

  } else {

    list[1] = io_make_output_field(
        "MaximalTemperatureTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, xparts,
        tracers_data.maximum_temperature_time,
        "Times at which the maximal temperature was reached");
  }

  list[2] =
      io_make_output_field("HeatedBySNIIFeedback", CHAR, 1, UNIT_CONV_NO_UNITS,
                           0.f, xparts, tracers_data.hit_by_SNII_feedback,
                           "Flags the particles that have been directly hit by "
                           "a SNII feedback event at some point in the past. "
                           "If > 0, contains the number of individual events.");

  list[3] =
      io_make_output_field("HeatedByAGNFeedback", CHAR, 1, UNIT_CONV_NO_UNITS,
                           0.f, xparts, tracers_data.hit_by_AGN_feedback,
                           "Flags the particles that have been directly hit by "
                           "an AGN feedback event at some point in the past. "
                           "If > 0, contains the number of individual events.");

  list[4] = io_make_output_field("EnergiesReceivedFromAGNFeedback", FLOAT, 1,
                                 UNIT_CONV_ENERGY, 0.f, xparts,
                                 tracers_data.AGN_feedback_energy,
                                 "Total amount of thermal energy from AGN "
                                 "feedback events received by the particles.");

  list[5] = io_make_output_field(
      "DensitiesBeforeLastAGNEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, xparts,
      tracers_data.density_before_last_AGN_feedback_event,
      "Physical density (not subgrid) of the gas fetched before the last AGN "
      "feedback event that hit the particles. -1 if the particles have never "
      "been heated.");

  list[6] = io_make_output_field(
      "EntropiesBeforeLastAGNEvent", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS,
      0.f, xparts, tracers_data.entropy_before_last_AGN_feedback_event,
      "Physical entropy (not subgrid) per unit mass of the gas fetched before "
      "the last AGN feedback event that hit the particles. -1 if the particles "
      "have never been heated.");

  list[7] = io_make_output_field(
      "DensitiesAtLastAGNEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, xparts,
      tracers_data.density_at_last_AGN_feedback_event,
      "Physical density (not subgrid) of the gas at the last AGN feedback "
      "event that hit the particles. -1 if the particles have never been "
      "heated.");

  list[8] = io_make_output_field(
      "EntropiesAtLastAGNEvent", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS, 0.f,
      xparts, tracers_data.entropy_at_last_AGN_feedback_event,
      "Physical entropy (not subgrid) per unit mass of the gas at the last AGN "
      "feedback event that hit the particles. -1 if the particles have never "
      "been heated.");

  if (with_cosmology) {

    list[9] = io_make_output_field(
        "LastAGNFeedbackScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        xparts, tracers_data.last_AGN_injection_scale_factor,
        "Scale-factors at which the particles were last hit by AGN feedback. "
        "-1 if a particle has never been hit by feedback");

  } else {

    list[9] = io_make_output_field(
        "LastAGNFeedbackTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, xparts,
        tracers_data.last_AGN_injection_time,
        "Times at which the particles were last hit by AGN feedback. -1 if a "
        "particle has never been hit by feedback");
  }

  return 10;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart* sparts, struct io_props* list,
    const int with_cosmology) {

  list[0] = io_make_output_field(
      "MaximalTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, sparts,
      tracers_data.maximum_temperature,
      "Maximal temperatures ever reached by the particles before they got "
      "converted to stars");

  if (with_cosmology) {
    list[1] = io_make_output_field(
        "MaximalTemperatureScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        sparts, tracers_data.maximum_temperature_scale_factor,
        "Scale-factors at which the maximal temperature was reached");

  } else {

    list[1] = io_make_output_field(
        "MaximalTemperatureTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.maximum_temperature_time,
        "Times at which the maximal temperature was reached");
  }

  list[2] =
      io_make_output_field("HeatedBySNIIFeedback", CHAR, 1, UNIT_CONV_NO_UNITS,
                           0.f, sparts, tracers_data.hit_by_SNII_feedback,
                           "Flags the particles that have been directly hit by "
                           "a SNII feedback event at some point in the past "
                           "when the particle was still a gas particle.");

  list[3] =
      io_make_output_field("HeatedByAGNFeedback", CHAR, 1, UNIT_CONV_NO_UNITS,
                           0.f, sparts, tracers_data.hit_by_AGN_feedback,
                           "Flags the particles that have been directly hit by "
                           "an AGN feedback event at some point in the past "
                           "when the particle was still a gas particle.");

  list[4] = io_make_output_field(
      "EnergiesReceivedFromAGNFeedback", FLOAT, 1, UNIT_CONV_ENERGY, 0.f,
      sparts, tracers_data.AGN_feedback_energy,
      "Total amount of thermal energy from AGN feedback events received by the "
      "particles when the particle was still a gas particle.");

  list[5] = io_make_output_field(
      "DensitiesBeforeLastAGNEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.density_before_last_AGN_feedback_event,
      "Physical density (not subgrid) of the gas fetched before the last AGN "
      "feedback "
      "event that hit the particles when they were still gas particles. -1 if "
      "the particles have never been heated.");

  list[6] = io_make_output_field(
      "EntropiesBeforeLastAGNEvent", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS,
      0.f, sparts, tracers_data.entropy_before_last_AGN_feedback_event,
      "Physical entropy (not subgrid) per unit mass of the gas fetched before "
      "the last AGN "
      "feedback event that hit the particles when they were still gas "
      "particles."
      " -1 if the particles have never been heated.");

  list[7] = io_make_output_field(
      "DensitiesAtLastAGNEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.density_at_last_AGN_feedback_event,
      "Physical density (not subgrid) of the gas at the last AGN feedback "
      "event that hit the particles when they were still gas particles. -1 if "
      "the particles have never been heated.");

  list[8] = io_make_output_field(
      "EntropiesAtLastAGNEvent", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS, 0.f,
      sparts, tracers_data.entropy_at_last_AGN_feedback_event,
      "Physical entropy (not subgrid) per unit mass of the gas at the last AGN "
      "feedback event that hit the particles when they were still gas "
      "particles."
      " -1 if the particles have never been heated.");

  if (with_cosmology) {

    list[9] = io_make_output_field(
        "LastAGNFeedbackScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        sparts, tracers_data.last_AGN_injection_scale_factor,
        "Scale-factors at which the particles were last hit by AGN feedback "
        "when they were still gas particles. -1 if a particle has never been "
        "hit by feedback");

  } else {

    list[9] =
        io_make_output_field("LastAGNFeedbackTimes", FLOAT, 1, UNIT_CONV_TIME,
                             0.f, sparts, tracers_data.last_AGN_injection_time,
                             "Times at which the particles were last hit by "
                             "AGN feedback when they were still gas particles. "
                             "-1 if a particle has never been hit by feedback");
  }

  return 10;
}

#endif /* SWIFT_TRACERS_EAGLE_IO_H */
