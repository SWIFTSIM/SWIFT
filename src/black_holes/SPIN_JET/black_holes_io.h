/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 * Copyright (c) 2022 Filip Husko (filip.husko@durham.ac.uk)
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
#ifndef SWIFT_SPIN_JET_BLACK_HOLES_IO_H
#define SWIFT_SPIN_JET_BLACK_HOLES_IO_H

#include "adiabatic_index.h"
#include "black_holes_part.h"
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
  *num_fields = 9;

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
  list[7] = io_make_input_field("Spins", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, bparts, spin);
  list[8] = io_make_input_field("AngularMomentumDirections", FLOAT, 3,
                                COMPULSORY, UNIT_CONV_NO_UNITS, bparts,
                                angular_momentum_direction);
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

INLINE static void convert_bpart_potential(const struct engine* e,
                                           const struct bpart* bp, float* ret) {

  if (bp->gpart != NULL)
    ret[0] = gravity_get_comoving_potential(bp->gpart);
  else
    ret[0] = 0.f;
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

INLINE static void convert_bpart_gas_velocity_dispersion(const struct engine* e,
                                                         const struct bpart* bp,
                                                         float* ret) {

  const struct cosmology* cosmo = e->cosmology;

  /* Conversion from internal to physical units */
  ret[0] = bp->velocity_dispersion_gas * cosmo->a_inv;
}

INLINE static void convert_bpart_gas_velocity_curl(const struct engine* e,
                                                   const struct bpart* bp,
                                                   float* ret) {

  const struct cosmology* cosmo = e->cosmology;

  /* Conversion from internal to physical units */
  ret[0] = bp->curl_v_gas[0] * cosmo->a2_inv;
  ret[1] = bp->curl_v_gas[1] * cosmo->a2_inv;
  ret[2] = bp->curl_v_gas[2] * cosmo->a2_inv;
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
                                               const int with_cosmology) {

  /* Say how much we want to write */
  *num_fields = 56;

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
      "GasDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, bparts, rho_gas,
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
      "TotalAccretedMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      total_accreted_mass,
      "Total mass accreted onto the main progenitor of the black holes since "
      "birth. This does not include any mass accreted onto any merged black "
      "holes.");

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

  list[27] = io_make_output_field_convert_bpart(
      "GasVelocityDispersions", FLOAT, 1, UNIT_CONV_SPEED, 0.f, bparts,
      convert_bpart_gas_velocity_dispersion,
      "Velocity dispersion (3D) of the gas particles around the black "
      "holes. This is a * sqrt(<|dx/dt|^2> - <|dx/dt|>^2) where x is the "
      "co-moving position of the particles relative to the black holes.");

  list[28] = io_make_output_field_convert_bpart(
      "GasCurlVelocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts,
      convert_bpart_gas_velocity_curl,
      "Velocity curl (3D) of the gas particles around the black holes.");

  list[29] = io_make_output_field(
      "AccretedAngularMomenta", FLOAT, 3, UNIT_CONV_ANGULAR_MOMENTUM, 0.f,
      bparts, accreted_angular_momentum,
      "Physical angular momenta that the black holes have accumulated through "
      "subgrid accretion.");

  list[30] = io_make_output_field(
      "NumberOfGasNeighbours", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      num_ngbs,
      "Integer number of gas neighbour particles within the black hole "
      "kernels.");

  list[31] = io_make_output_field(
      "NumberOfHeatingEvents", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      AGN_number_of_energy_injections,
      "Integer number of (thermal) energy injections the black hole has had "
      "so far. This counts each heated gas particle separately, and so can "
      "increase by more than one during a single time step.");

  list[32] = io_make_output_field(
      "NumberOfAGNEvents", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      AGN_number_of_AGN_events,
      "Integer number of AGN events the black hole has had so far"
      " (the number of time steps the BH did AGN feedback)");

  if (with_cosmology) {
    list[33] = io_make_output_field(
        "LastAGNFeedbackScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_AGN_event_scale_factor,
        "Scale-factors at which the black holes last had an AGN event.");
  } else {
    list[33] = io_make_output_field(
        "LastAGNFeedbackTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_AGN_event_time,
        "Times at which the black holes last had an AGN event.");
  }

  list[34] = io_make_output_field(
      "AccretionLimitedTimeSteps", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
      dt_heat,
      "Accretion-limited time steps of black holes. The actual time step of "
      "the particles may differ due to the minimum allowed value.");

  list[35] = io_make_output_field(
      "AGNTotalInjectedEnergies", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, bparts,
      AGN_cumulative_energy,
      "Total (cumulative) physical energies injected into gas particles "
      "in AGN feedback.");

  list[36] = io_make_output_field_convert_bpart(
      "Potentials", FLOAT, 1, UNIT_CONV_POTENTIAL, -1.f, bparts,
      convert_bpart_potential, "Gravitational potentials of the particles");

  list[37] = io_make_output_field(
      "Spins", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts, spin,
      "Dimensionless spins of the black holes. "
      "Negative values indicate retrograde accretion.");

  list[38] = io_make_output_field(
      "AngularMomentumDirections", FLOAT, 3, UNIT_CONV_NO_UNITS, 0.f, bparts,
      angular_momentum_direction,
      "Direction of the black hole spin vector, normalised to unity.");

  list[39] = io_make_output_field(
      "AccretionDiscAspectRatios", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      aspect_ratio,
      "The aspect ratio, h/r, of the subgrid accretion disc "
      "around the black hole.");

  list[40] = io_make_output_field(
      "JetEfficiencies", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      jet_efficiency, "Jet power divided by accretion rate.");

  list[41] = io_make_output_field(
      "RadiativeEfficiencies", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      radiative_efficiency, "AGN luminosity divided by accretion rate.");

  list[42] = io_make_output_field("CosAccretionDiskAngle", FLOAT, 1,
                                  UNIT_CONV_NO_UNITS, 0.f, bparts,
                                  accretion_disk_angle,
                                  "Cosine of the angle between the spin vector "
                                  "and the accreting gas angular momentum.");

  list[43] = io_make_output_field(
      "AccretionModes", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts, accretion_mode,
      "Accretion flow regime. 0 - Thick disk, 1 - Thin disk, 2 - Slim disk");

  list[44] = io_make_output_field(
      "JetReservoir", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, bparts, jet_reservoir,
      "Total jet energy waiting to be released (once it "
      "grows large enough to kick a single particle).");

  list[45] = io_make_output_field(
      "InjectedJetEnergies", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, bparts,
      total_jet_energy, "Total jet energy injected into AGN surroundings.");

  list[46] = io_make_output_field(
      "JetTimeSteps", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts, dt_jet,
      "Jet-launching-limited time-steps of black holes.");

  list[47] = io_make_output_field(
      "NumberOfJetParticlesLaunched", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      AGN_number_of_jet_injections,
      "Integer number of (kinetic) energy injections the black hole has had "
      "so far");

  list[48] = io_make_output_field(
      "NumberOfAGNJetEvents", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      AGN_number_of_AGN_jet_events,
      "Integer number of AGN jet launching events the black hole has had"
      " (the number of times the BH did AGN jet feedback)");

  if (with_cosmology) {
    list[49] = io_make_output_field(
        "LastAGNJetScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
        last_AGN_jet_event_scale_factor,
        "Scale-factors at which the black holes last had an AGN jet event.");
  } else {
    list[49] = io_make_output_field(
        "LastAGNJetTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_AGN_jet_event_time,
        "Times at which the black holes last had an AGN jet event.");
  }

  list[50] = io_make_output_field(
      "EddingtonFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      eddington_fraction,
      "Accretion rates of black holes in units of their Eddington rates. "
      "This is based on the unlimited accretion rates, so these fractions "
      "can be above the limiting fEdd.");

  list[51] = io_make_output_field(
      "FOFGroupMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts, group_mass,
      "Parent halo masses of the black holes, as determined from the FOF  "
      "algorithm.");

  list[52] =
      io_make_output_field("JetVelocities", FLOAT, 1, UNIT_CONV_VELOCITY, 0.f,
                           bparts, v_jet, "The current jet velocities.");

  list[53] = io_make_output_field(
      "TotalAccretedMassesByMode", FLOAT, 3, UNIT_CONV_MASS, 0.f, bparts,
      accreted_mass_by_mode,
      "The total accreted mass in each accretion mode. The components to the "
      "mass accreted in the thick, thin and slim disc modes, respectively.");

  list[54] = io_make_output_field(
      "AGNTotalInjectedEnergiesByMode", FLOAT, 3, UNIT_CONV_ENERGY, 0.f, bparts,
      thermal_energy_by_mode,
      "The total energy injected in the thermal AGN feedback mode, split by "
      "accretion mode. The components correspond to the thermal energy dumped "
      "in the thick, thin and slim disc modes, respectively.");

  list[55] = io_make_output_field(
      "InjectedJetEnergiesByMode", FLOAT, 3, UNIT_CONV_ENERGY, 0.f, bparts,
      jet_energy_by_mode,
      "The total energy injected in the kinetic jet AGN feedback mode, split "
      "by accretion mode. The components correspond to the thermal energy "
      "dumped in the thick, thin and slim disc modes, respectively.");

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

#endif /* SWIFT_SPIN_JET_BLACK_HOLES_IO_H */
