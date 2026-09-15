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

INLINE static void convert_bpart_averaged_accretion_rate(const struct engine *e,
                                                         const struct bpart *bp,
                                                         float *ret) {

  for (int i = 0; i < num_snapshot_triggers_bpart; ++i) {
    if (e->snapshot_recording_triggers_started_bpart[i]) {
      ret[i] = bp->tracers_data.averaged_accretion_rate[i] /
               e->snapshot_recording_triggers_bpart[i];
    } else {
      ret[i] = 0.f;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (ret[i] < 0.f)
      error(
          "Negative averaged accretion rate for black hole id=%lld "
          "trigger=%d value=%e",
          bp->id, i, ret[i]);
#endif
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

  int num = 6;

  list[0] = io_make_physical_output_field(
      "CumulativeMomentumFromSupernovae", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f,
      xparts, tracers_data.feedback_cumulative.momentum_supernovae,
      /*can convert to comoving=*/0,
      "Cumulative |delta_p| per event from supernovae over this particle's "
      "lifetime (scalar sum, not vector: isotropic kicks would else "
      "cancel).");

  list[1] = io_make_physical_output_field(
      "CumulativeMomentumFromWinds", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f, xparts,
      tracers_data.feedback_cumulative.momentum_winds,
      /*can convert to comoving=*/0,
      "Same convention as CumulativeMomentumFromSupernovae, for stellar "
      "winds.");

  list[2] = io_make_physical_output_field(
      "CumulativeEnergyFromSupernovae", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, xparts,
      tracers_data.feedback_cumulative.energy_supernovae,
      /*can convert to comoving=*/0,
      "Cumulative specific internal energy received from supernovae over "
      "this particle's lifetime.");

  list[3] = io_make_physical_output_field(
      "CumulativeEnergyFromWinds", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, xparts, tracers_data.feedback_cumulative.energy_winds,
      /*can convert to comoving=*/0,
      "Cumulative specific internal energy received from stellar winds. "
      "A conservation residual, not strictly positive: can go negative "
      "when the gas was already moving towards the star before the kick.");

  list[4] = io_make_physical_output_field(
      "MaxKickVelocityFromSupernovae", FLOAT, 1, UNIT_CONV_SPEED, 0.f, xparts,
      tracers_data.feedback_cumulative.max_kick_velocity_supernovae,
      /*can convert to comoving=*/0,
      "Largest single-event kick velocity this particle received from "
      "supernovae (outflow diagnostic).");

  list[5] = io_make_physical_output_field(
      "MaxKickVelocityFromWinds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, xparts,
      tracers_data.feedback_cumulative.max_kick_velocity_winds,
      /*can convert to comoving=*/0,
      "Same convention as MaxKickVelocityFromSupernovae, for stellar "
      "winds.");

  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart *sparts, struct io_props *list,
    const int with_cosmology) {

  int num = 6;

  list[0] = io_make_output_field(
      "NumberOfSNIIEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      tracers_data.snii_events.n_events,
      "Number of SNII events this star produced over its lifetime so far "
      "(fractional for a continuously-sampled population particle; always "
      "0 or 1 for a discrete star).");

  list[1] = io_make_physical_output_field(
      "DensityAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.snii_events.density_at_last_event,
      /*can convert to comoving=*/0,
      "Gas density at the star's location at its most recent SNII event. "
      "0 if it has never had one.");

  if (with_cosmology) {
    list[2] = io_make_physical_output_field(
        "ScaleFactorAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.snii_events.last_event_scale_factor,
        /*can convert to comoving=*/0,
        "Scale-factor at this star's most recent SNII event. 0 if it has "
        "never had one.");
  } else {
    list[2] = io_make_output_field(
        "TimeAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.snii_events.last_event_time,
        "Simulation time at this star's most recent SNII event. 0 if it "
        "has never had one.");
  }

  list[3] = io_make_output_field(
      "NumberOfSNIaEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      tracers_data.snia_events.n_events,
      "Number of SNIa events this star produced over its lifetime so far. "
      "Always 0 for a discrete (single_star) particle: SNIa is a "
      "population-level channel in this model.");

  list[4] = io_make_physical_output_field(
      "DensityAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.snia_events.density_at_last_event,
      /*can convert to comoving=*/0,
      "Same as DensityAtLastSNIIEvent, for the SNIa channel.");

  if (with_cosmology) {
    list[5] = io_make_physical_output_field(
        "ScaleFactorAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.snia_events.last_event_scale_factor,
        /*can convert to comoving=*/0,
        "Same as ScaleFactorAtLastSNIIEvent, for the SNIa channel.");
  } else {
    list[5] = io_make_output_field(
        "TimeAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.snia_events.last_event_time,
        "Same as TimeAtLastSNIIEvent, for the SNIa channel.");
  }

  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_bparticles(
    const struct bpart *bparts, struct io_props *list,
    const int with_cosmology) {

  list[0] = io_make_output_field_convert_bpart(
      "AveragedAccretionRates", FLOAT, num_snapshot_triggers_bpart,
      UNIT_CONV_MASS_PER_UNIT_TIME, 0.f, bparts,
      convert_bpart_averaged_accretion_rate,
      "Accretion rates of the black holes averaged over the period set by "
      "the first N snapshot triggers");

  return 1;
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
