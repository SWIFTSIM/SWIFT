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
#include "feedback.h"
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
 * @brief Snapshot converter for #IsIonizedFlags, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_is_ionized(const struct engine *e,
                                           const struct part *p,
                                           const struct xpart *xp, char *ret) {
  ret[0] = feedback_is_part_tagged_as_ionized(p, xp);
}

/**
 * @brief Snapshot converter for #HIIStarIDs, see #tracers_write_particles.
 */
INLINE static void convert_part_HII_star_id(const struct engine *e,
                                            const struct part *p,
                                            const struct xpart *xp,
                                            long long *ret) {
  ret[0] = feedback_get_part_ionized_star_id(p, xp);
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

  int num = 10;

  /* The tag core (is_ionized/star_id) lives on struct part's feedback_data,
     not tracers_xpart_data, so read it through the feedback-model dispatch
     wrapper: this must compile under any --with-feedback choice paired
     with --with-tracers=GEAR, not only --with-feedback=GEAR. */
  list[0] = io_make_output_field_convert_part(
      "IsIonizedFlags", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_is_ionized,
      "Were the particles flagged as ionized by HII ionzation subgrid model?");

  list[1] = io_make_output_field_convert_part(
      "HIIStarIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_HII_star_id,
      "Star particle IDs that ionized these gas particles due to HII ionzation "
      "subgrid model?");

  list[2] = io_make_output_field(
      "CumulativeMomentumFromSN", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f, xparts,
      tracers_data.feedback_cumulative.momentum_SN,
      "Cumulative momentum magnitude received from SN over this particle's "
      "lifetime (scalar sum of |delta_p| per event, not a vector sum: "
      "isotropic kicks would otherwise cancel out).");

  list[3] = io_make_output_field(
      "CumulativeMomentumFromWinds", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f, xparts,
      tracers_data.feedback_cumulative.momentum_winds,
      "Cumulative momentum magnitude received from stellar winds over this "
      "particle's lifetime. Same scalar-sum convention as "
      "CumulativeMomentumFromSN.");

  list[4] = io_make_output_field(
      "CumulativeMomentumFromRadiationPressure", FLOAT, 1, UNIT_CONV_MOMENTUM,
      0.f, xparts, tracers_data.feedback_cumulative.momentum_radiation,
      "Cumulative momentum magnitude received from radiation pressure over "
      "this particle's lifetime. Same scalar-sum convention as "
      "CumulativeMomentumFromSN.");

  list[5] = io_make_output_field(
      "CumulativeEnergyFromSN", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      xparts, tracers_data.feedback_cumulative.energy_SN,
      "Cumulative specific internal energy received from SN over this "
      "particle's lifetime.");

  list[6] = io_make_output_field(
      "CumulativeEnergyFromWinds", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, xparts, tracers_data.feedback_cumulative.energy_winds,
      "Cumulative specific internal energy received from stellar winds over "
      "this particle's lifetime. Unlike CumulativeEnergyFromSN this is a "
      "conservation residual (budgeted energy minus the actual kinetic "
      "energy change from the kick), not a strictly positive injected "
      "quantity: known to go negative when the gas was already moving "
      "towards the star before the kick.");

  list[7] = io_make_output_field(
      "MaxKickVelocityFromSN", FLOAT, 1, UNIT_CONV_SPEED, 0.f, xparts,
      tracers_data.feedback_cumulative.max_kick_velocity_SN,
      "Largest single-event kick velocity this particle received from SN "
      "(outflow diagnostic: peak coupling speed near the source, before "
      "deceleration).");

  list[8] = io_make_output_field(
      "MaxKickVelocityFromWinds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, xparts,
      tracers_data.feedback_cumulative.max_kick_velocity_winds,
      "Largest single-event kick velocity this particle received from "
      "stellar winds. Same convention as MaxKickVelocityFromSN.");

  list[9] = io_make_output_field(
      "MaxKickVelocityFromRadiationPressure", FLOAT, 1, UNIT_CONV_SPEED, 0.f,
      xparts, tracers_data.feedback_cumulative.max_kick_velocity_radiation,
      "Largest single-event kick velocity this particle received from "
      "radiation pressure. Same convention as MaxKickVelocityFromSN.");

  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart *sparts, struct io_props *list,
    const int with_cosmology) {

  int num = 8;

  list[0] = io_make_output_field(
      "FinalHIIRegionRadii", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sparts,
      tracers_data.final_HII_radius,
      "Co-moving HII region radius of the star particles before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionRadii it is retired from.");

  list[1] = io_make_output_field(
      "FinalHIIRegionMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      tracers_data.final_HII_mass,
      "Ionized gas mass of the star particles' HII region before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionMasses it is retired from.");

  list[2] = io_make_output_field(
      "NumberOfSNIIEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      tracers_data.snii_events.n_events,
      "Number of SNII events this star produced over its lifetime so far "
      "(fractional for a continuously-sampled population particle; always "
      "0 or 1 for a discrete star).");

  list[3] = io_make_output_field(
      "DensityAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.snii_events.density_at_last_event,
      "Physical gas density at the star's own location (kernel-averaged, "
      "same value the enrichment loop already computes) at its most recent "
      "SNII event. 0 if it has never had one. Compare NumberOfSNIIEvents "
      "between two snapshots to see if any events were missed between "
      "them.");

  if (with_cosmology) {
    list[4] = io_make_output_field(
        "ScaleFactorAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.snii_events.last_event_scale_factor,
        "Scale-factor at this star's most recent SNII event. 0 if it has "
        "never had one.");
  } else {
    list[4] = io_make_output_field(
        "TimeAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.snii_events.last_event_time,
        "Simulation time at this star's most recent SNII event. 0 if it "
        "has never had one.");
  }

  list[5] = io_make_output_field(
      "NumberOfSNIaEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      tracers_data.snia_events.n_events,
      "Number of SNIa events this star produced over its lifetime so far. "
      "Always 0 for a discrete (single_star) particle: SNIa is a "
      "population-level channel in this model.");

  list[6] = io_make_output_field(
      "DensityAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.snia_events.density_at_last_event,
      "Same as DensityAtLastSNIIEvent, for the SNIa channel.");

  if (with_cosmology) {
    list[7] = io_make_output_field(
        "ScaleFactorAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.snia_events.last_event_scale_factor,
        "Same as ScaleFactorAtLastSNIIEvent, for the SNIa channel.");
  } else {
    list[7] = io_make_output_field(
        "TimeAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.snia_events.last_event_time,
        "Same as TimeAtLastSNIIEvent, for the SNIa channel.");
  }

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
