/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_TRACERS_NONE_IO_H
#define SWIFT_TRACERS_NONE_IO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "io_properties.h"
#include "tracers.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of tracers to the file.
 *
 * @param h_grp The HDF5 group in which to write
 * @param tracers The #tracers_function_data
 */
__attribute__((always_inline)) INLINE static void tracers_write_flavour(
    hid_t h_grp) {

  io_write_attribute_s(h_grp, "Tracers", "GEAR");
}
#endif

INLINE static void convert_sink_averaged_SFR(const struct engine *e,
                                             const struct sink *sink,
                                             float *ret) {

  for (int i = 0; i < num_snapshot_triggers_sink; ++i) {
    if (e->snapshot_recording_triggers_started_sink[i]) {
      ret[i] = sink->tracers_data.averaged_SFR[i] /
               e->snapshot_recording_triggers_sink[i];
    } else {
      ret[i] = 0.f;
    }
  }
}

INLINE static void convert_sink_averaged_accretion_rate(const struct engine *e,
                                                        const struct sink *sink,
                                                        float *ret) {

  for (int i = 0; i < num_snapshot_triggers_sink; ++i) {
    if (e->snapshot_recording_triggers_started_sink[i]) {
      ret[i] = sink->tracers_data.averaged_accretion_rate[i] /
               e->snapshot_recording_triggers_sink[i];
    } else {
      ret[i] = 0.f;
    }
  }
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int tracers_write_particles(
    const struct part *parts, const struct xpart *xparts, struct io_props *list,
    const int with_cosmology) {

  int num = 2;

  list[0] = io_make_output_field(
      "IsIonizedFlags", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      tracers_data.HII_region.is_ionized,
      "Were the particles flagged as ionized by HII ionzation subgrid model?");

  list[1] = io_make_output_field(
      "HIIStarIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      tracers_data.HII_region.star_id,
      "Star particle IDs that ionized these gas particles due to HII ionzation "
      "subgrid model?");

  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart *sparts, struct io_props *list,
    const int with_cosmology) {

  int num = 1;

  list[0] = io_make_output_field(
      "FinalHIIRegionRadii", FLOAT, 1, UNIT_CONV_LENGTH, 0.f, sparts,
      tracers_data.final_HII_radius,
      "Co-moving HII region radius of the star particles before they die or "
      "were not eligible to form HII regions anymore.");
  
  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_bparticles(
    const struct bpart *bparts, struct io_props *list,
    const int with_cosmology) {

  return 0;
}

__attribute__((always_inline)) INLINE static int tracers_write_sinkparticles(
    const struct sink *sinks, struct io_props *list, const int with_cosmology) {

  list[0] = io_make_output_field_convert_sink(
      "AveragedAccretionRates", FLOAT, num_snapshot_triggers_sink,
      UNIT_CONV_MASS_PER_UNIT_TIME, 0.f, sinks,
      convert_sink_averaged_accretion_rate,
      "Accretion rates of the sinks averaged over the period set by the "
      "first N snapshot triggers");

  list[1] = io_make_output_field_convert_sink(
      "AveragedStarFormationRates", FLOAT, num_snapshot_triggers_sink,
      UNIT_CONV_SFR, 0.f, sinks, convert_sink_averaged_SFR,
      "Star formation rates of the particles averaged over the period set by "
      "the first N snapshot triggers");

  return 2;
}

#endif /* SWIFT_TRACERS_NONE_IO_H */
