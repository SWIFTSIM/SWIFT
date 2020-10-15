/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_BLACK_HOLES_IO_H
#define SWIFT_EAGLE_BLACK_HOLES_IO_H

#include "adiabatic_index.h"
#include "black_holes_part.h"
#include "black_holes_properties.h"
#include "io_properties.h"

/**
 * @brief Specifies which b-particle fields to read from a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void black_holes_read_particles(struct bpart* bparts,
                                              struct io_props* list,
                                              int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 7;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, bparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, bparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                bparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, bparts, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, bparts, h);
  list[5] = io_make_input_field("EnergyReservoir", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY, bparts, energy_reservoir);
  list[6] = io_make_input_field("SubgridMasses", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_MASS, bparts, subgrid_mass);
}

INLINE static void convert_bpart_pos(const struct engine* e,
                                     const struct bpart* bp, double* ret) {

  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(bp->x[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(bp->x[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(bp->x[2], 0.0, s->dim[2]);
  } else {
    ret[0] = bp->x[0];
    ret[1] = bp->x[1];
    ret[2] = bp->x[2];
  }
}

INLINE static void convert_bpart_vel(const struct engine* e,
                                     const struct bpart* bp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, bp->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, bp->time_bin);

  /* Get time-step since the last kick */
  float dt_kick_grav;
  if (with_cosmology) {
    dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_grav -=
        cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
  } else {
    dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
  }

  /* Extrapolate the velocites to the current time */
  const struct gpart* gp = bp->gpart;
  ret[0] = gp->v_full[0] + gp->a_grav[0] * dt_kick_grav;
  ret[1] = gp->v_full[1] + gp->a_grav[1] * dt_kick_grav;
  ret[2] = gp->v_full[2] + gp->a_grav[2] * dt_kick_grav;

  /* Conversion from internal to physical units */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

INLINE static void convert_bpart_gas_vel(const struct engine* e,
                                         const struct bpart* bp, float* ret) {

  const struct cosmology* cosmo = e->cosmology;

  /* Convert relative velocities to physical units */
  ret[0] = bp->velocity_gas[0] * cosmo->a_inv;
  ret[1] = bp->velocity_gas[1] * cosmo->a_inv;
  ret[2] = bp->velocity_gas[2] * cosmo->a_inv;
}

INLINE static void convert_bpart_gas_circular_vel(const struct engine* e,
                                                  const struct bpart* bp,
                                                  float* ret) {

  const struct cosmology* cosmo = e->cosmology;

  /* Conversion from internal to physical units */
  ret[0] = bp->circular_velocity_gas[0] * cosmo->a_inv;
  ret[1] = bp->circular_velocity_gas[1] * cosmo->a_inv;
  ret[2] = bp->circular_velocity_gas[2] * cosmo->a_inv;
}

INLINE static void convert_bpart_gas_temperatures(const struct engine* e,
                                                  const struct bpart* bp,
                                                  float* ret) {

  const struct black_holes_props* props = e->black_holes_properties;
  const struct cosmology* cosmo = e->cosmology;

  /* Conversion from specific internal energy to temperature */
  ret[0] = bp->internal_energy_gas * cosmo->a_factor_internal_energy /
           props->temp_to_u_factor;
}

/**
 * @brief Specifies which b-particle fields to write to a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * @param with_cosmology Are we running a cosmological simulation?
 */
INLINE static void black_holes_write_particles(const struct bpart* bparts,
                                               struct io_props* list,
                                               int* num_fields,
                                               int with_cosmology) {

  /* Say how much we want to write */
  *num_fields = 42;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_bpart(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, bparts,
      convert_bpart_pos, "Co-moving position of the particles");

  list[1] = io_make_output_field_convert_bpart(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts, convert_bpart_vel,
      "Peculiar velocities of the particles. This is a * dx/dt where x is the "
      "co-moving position of the particles.");

  list[2] =
      io_make_output_field("DynamicalMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, mass, "Dynamical masses of the particles");

  list[3] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           bparts, id, "Unique ID of the particles");

  list[4] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, bparts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[5] = io_make_output_field("SubgridMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 bparts, subgrid_mass,
                                 "Subgrid masses of the particles");

  if (with_cosmology) {
    list[6] = io_make_output_field(
        "FormationScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
        formation_scale_factor, "Scale-factors at which the BHs were formed");
  } else {
    list[6] = io_make_output_field("FormationTimes", FLOAT, 1, UNIT_CONV_TIME,
                                   0.f, bparts, formation_time,
                                   "Times at which the BHs were formed");
  }

  list[7] = io_make_output_field(
      "GasDensities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f, bparts, rho_gas,
      "Co-moving densities of the gas around the particles");

  list[8] = io_make_output_field(
      "GasSoundSpeeds", FLOAT, 1, UNIT_CONV_SPEED,
      -1.5f * hydro_gamma_minus_one, bparts, sound_speed_gas,
      "Co-moving sound-speeds of the gas around the particles");

  list[9] = io_make_output_field(
      "EnergyReservoirs", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, bparts,
      energy_reservoir,
      "Physcial energy contained in the feedback reservoir of the particles");

  list[10] = io_make_output_field(
      "AccretionRates", FLOAT, 1, UNIT_CONV_MASS_PER_UNIT_TIME, 0.f, bparts,
      accretion_rate,
      "Physical instantaneous accretion rates of the particles");

  list[11] = io_make_output_field(
      "TotalAccretedMasses", FLOAT, 1, UNIT_CONV_MASS_PER_UNIT_TIME, 0.f,
      bparts, total_accreted_mass,
      "Total mass accreted onto the particles since its birth");

  list[12] = io_make_output_field(
      "CumulativeNumberOfSeeds", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      cumulative_number_seeds,
      "Total number of BH seeds that have merged into this black hole");

  list[13] =
      io_make_output_field("NumberOfMergers", INT, 1, UNIT_CONV_NO_UNITS, 0.f,
                           bparts, number_of_mergers,
                           "Number of mergers the black holes went through. "
                           "This does not include the number of mergers "
                           "accumulated by any merged black hole.");

  if (with_cosmology) {
    list[14] = io_make_output_field(
        "LastHighEddingtonFractionScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS,
        0.f, bparts, last_high_Eddington_fraction_scale_factor,
        "Scale-factors at which the black holes last reached a large Eddington "
        "ratio. -1 if never reached.");
  } else {
    list[14] = io_make_output_field(
        "LastHighEddingtonFractionTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_high_Eddington_fraction_time,
        "Times at which the black holes last reached a large Eddington ratio. "
        "-1 if never reached.");
  }

  if (with_cosmology) {
    list[15] = io_make_output_field(
        "LastMinorMergerScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_minor_merger_scale_factor,
        "Scale-factors at which the black holes last had a minor merger.");
  } else {
    list[15] = io_make_output_field(
        "LastMinorMergerScaleTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_minor_merger_time,
        "Times at which the black holes last had a minor merger.");
  }

  if (with_cosmology) {
    list[16] = io_make_output_field(
        "LastMajorMergerScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_major_merger_scale_factor,
        "Scale-factors at which the black holes last had a major merger.");
  } else {
    list[16] = io_make_output_field(
        "LastMajorMergerScaleTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_major_merger_time,
        "Times at which the black holes last had a major merger.");
  }

  list[17] = io_make_output_field(
      "SwallowedAngularMomenta", FLOAT, 3, UNIT_CONV_ANGULAR_MOMENTUM, 0.f,
      bparts, swallowed_angular_momentum,
      "Physical angular momenta that the black holes have accumulated by "
      "swallowing gas particles.");

  list[18] = io_make_output_field_convert_bpart(
      "GasRelativeVelocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts,
      convert_bpart_gas_vel,
      "Peculiar relative velocities of the gas particles around the black "
      "holes. This is a * dx/dt where x is the co-moving position of the "
      "particles.");

  list[19] = io_make_output_field_convert_bpart(
      "GasCircularVelocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts,
      convert_bpart_gas_circular_vel,
      "Circular velocities of the gas around the black hole at the "
      "smoothing radius. This is j / h_BH, where j is the smoothed, peculiar "
      "specific angular momentum of gas around the black holes, and h_BH is "
      "the smoothing length of each black hole.");

  list[20] =
      io_make_output_field("TimeBins", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
                           time_bin, "Time-bins of the particles");

  list[21] = io_make_output_field(
      "NumberOfSwallows", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_gas_swallows,
      "Number of gas particles the black holes have swallowed. "
      "This includes the particles swallowed by any of the black holes that "
      "merged into this one.");

  list[22] = io_make_output_field(
      "NumberOfDirectSwallows", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_direct_gas_swallows,
      "Number of gas particles the black holes have swallowed. "
      "This does not include any particles swallowed by any of the black holes "
      "that merged into this one.");

  list[23] = io_make_output_field(
      "NumberOfRepositions", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_repositions,
      "Number of repositioning events the black holes went through. This does "
      "not include the number of reposition events accumulated by any merged "
      "black holes.");

  list[24] = io_make_output_field(
      "NumberOfRepositionAttempts", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_reposition_attempts,
      "Number of time steps in which the black holes had an eligible particle "
      "to reposition to. They may or may not have ended up moving there, "
      "depending on their subgrid mass and on whether these particles were at "
      "a lower or higher potential than the black holes themselves. It does "
      "not include attempted repositioning events accumulated by any merged "
      "black holes.");

  list[25] = io_make_output_field(
      "NumberOfTimeSteps", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_time_steps,
      "Total number of time steps at which the black holes were active.");

  list[26] = io_make_output_field(
      "ViscosityFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts, f_visc,
      "Multiplicative factors by which the Bondi-Hoyle-Lyttleton accretion "
      "rates have been suppressed by the Rosas-Guevara et al. (2015) "
      "accretion disc model.");

  list[27] = io_make_output_field(
      "SubgridDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, bparts,
      rho_subgrid_gas,
      "Physical subgrid densities used in the subgrid-Bondi model.");

  list[28] = io_make_output_field(
      "SubgridSoundSpeeds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, bparts,
      sound_speed_subgrid_gas,
      "Physical subgrid sound-speeds used in the subgrid-Bondi model.");

  list[29] = io_make_output_field(
      "BirthGasDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, bparts,
      formation_gas_density,
      "Physical densities of the converted part at the time of birth. "
      "We store the physical density at the birth redshift, no conversion is "
      "needed.");

  list[30] = io_make_output_field(
      "AccretedAngularMomenta", FLOAT, 3, UNIT_CONV_ANGULAR_MOMENTUM, 0.f,
      bparts, accreted_angular_momentum,
      "Physical angular momenta that the black holes have accumulated through "
      "subgrid accretion.");

  list[31] = io_make_output_field(
      "NumberOfGasNeighbours", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      num_ngbs,
      "Integer number of gas neighbour particles within the black hole "
      "kernels.");

  list[32] = io_make_output_field(
      "FeedbackDeltaT", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, bparts,
      AGN_delta_T,
      "Temperature by which gas particles have been heated by the black hole "
      "particles in the most recent feedback event.");

  list[33] = io_make_output_field(
      "NumberOfHeatingEvents", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      AGN_number_of_energy_injections,
      "Integer number of (thermal) energy injections the black hole has had "
      "so far");

  list[34] = io_make_output_field(
      "NumberOfAGNEvents", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      AGN_number_of_AGN_events,
      "Integer number of AGN events the black hole has had so far"
      " (the number of times the BH did AGN feedback)");

  if (with_cosmology) {
    list[35] = io_make_output_field(
        "LastAGNFeedbackScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_AGN_event_scale_factor,
        "Scale-factors at which the black holes last had an AGN event.");
  } else {
    list[35] = io_make_output_field(
        "LastAGNFeedbackTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_AGN_event_time,
        "Times at which the black holes last had an AGN event.");
  }

  list[36] = io_make_output_field(
      "AccretionLimitedTimeSteps", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
      dt_heat, "Accretion-limited time-steps of black holes.");

  list[37] = io_make_output_field(
      "AGNTotalInjectedEnergies", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, bparts,
      AGN_cumulative_energy,
      "Total (cumulative) physical energies injected into gas particles "
      "in AGN feedback.");

  list[38] = io_make_output_field(
      "AccretionBoostFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      accretion_boost_factor,
      "Multiplicative factors by which the Bondi-Hoyle-Lyttleton accretion "
      "rates have been increased by the density-dependent Booth & Schaye "
      "(2009) accretion model.");

  list[39] = io_make_output_field_convert_bpart(
      "GasTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, bparts,
      convert_bpart_gas_temperatures,
      "Temperature of the gas surrounding the black holes.");

  list[40] = io_make_output_field(
      "EnergyReservoirThresholds", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      num_ngbs_to_heat,
      "Minimum energy reservoir required for the black holes to do feedback, "
      "expressed in units of the (constant) target heating temperature "
      "increase.");

  list[41] = io_make_output_field(
      "EddingtonFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      eddington_fraction,
      "Accretion rates of black holes in units of their Eddington rates. "
      "This is based on the unlimited accretion rates, so these fractions "
      "can be above the limiting fEdd.");

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES

  list += *num_fields;
  *num_fields += 4;

  list[0] = io_make_output_field("Num_ngb_density", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_density);
  list[1] = io_make_output_field("Num_ngb_force", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_force);
  list[2] = io_make_output_field("Ids_ngb_density", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_density);
  list[3] = io_make_output_field("Ids_ngb_force", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_force);
#endif
}

#endif /* SWIFT_EAGLE_BLACK_HOLES_IO_H */
