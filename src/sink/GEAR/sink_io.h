/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_GEAR_SINK_IO_H
#define SWIFT_GEAR_SINK_IO_H

#include "io_properties.h"
#include "sink_part.h"

/**
 * @brief Specifies which sink-particle fields to read from a dataset
 *
 * @param sinks The sink-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void sink_read_particles(struct sink* sinks,
                                       struct io_props* list, int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 6;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, sinks, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, sinks, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sinks, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, sinks, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, sinks, h);
  list[5] = io_make_input_field("BirthTime", FLOAT, 1, OPTIONAL, UNIT_CONV_MASS,
                                sinks, birth_time);
}

INLINE static void convert_sink_pos(const struct engine* e,
                                    const struct sink* sp, double* ret) {

  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(sp->x[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(sp->x[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(sp->x[2], 0.0, s->dim[2]);
  } else {
    ret[0] = sp->x[0];
    ret[1] = sp->x[1];
    ret[2] = sp->x[2];
  }
}

INLINE static void convert_sink_vel(const struct engine* e,
                                    const struct sink* sp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, sp->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, sp->time_bin);

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
  const struct gpart* gp = sp->gpart;
  ret[0] = gp->v_full[0] + gp->a_grav[0] * dt_kick_grav;
  ret[1] = gp->v_full[1] + gp->a_grav[1] * dt_kick_grav;
  ret[2] = gp->v_full[2] + gp->a_grav[2] * dt_kick_grav;

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

INLINE static void convert_sink_target_mass(const struct engine* e,
                                            const struct sink* sink,
                                            float* ret) {
  /* Recall that the target_mass_Msun is in M_sun in the code. We nee to convert
     it to internal units for consistency in the output. */
  ret[0] = sink->target_mass_Msun * e->physical_constants->const_solar_mass;
}

INLINE static void convert_sink_swallowed_angular_momentum(
    const struct engine* e, const struct sink* sink, float* ret) {
  ret[0] = sink->swallowed_angular_momentum[0];
  ret[1] = sink->swallowed_angular_momentum[1];
  ret[2] = sink->swallowed_angular_momentum[2];
}

/**
 * @brief Specifies which sink-particle fields to write to a dataset
 *
 * @param sinks The sink-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * @param with_cosmology Are we running a cosmological simulation?
 */
INLINE static void sink_write_particles(const struct sink* sinks,
                                        struct io_props* list, int* num_fields,
                                        int with_cosmology) {

  /* Say how much we want to write */
  *num_fields = 11;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_sink(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, sinks, convert_sink_pos,
      "Co-moving position of the particles");

  list[1] = io_make_output_field_convert_sink(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, sinks, convert_sink_vel,
      "Peculiar velocities of the particles. This is a * dx/dt where x is the "
      "co-moving position of the particles.");

  list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f, sinks,
                                 mass, "Masses of the particles");

  list[3] = io_make_physical_output_field(
      "ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, sinks, id,
      /*can convert to comoving=*/0, "Unique ID of the particles");

  list[4] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sinks, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[5] = io_make_physical_output_field(
      "NumberOfSinkSwallows", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sinks,
      number_of_sink_swallows, /*can convert to comoving=*/0,
      "Total number of sink merger events");

  list[6] = io_make_physical_output_field(
      "NumberOfGasSwallows", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sinks,
      number_of_gas_swallows, /*can convert to comoving=*/0,
      "Total number of gas merger events");

  list[7] = io_make_output_field_convert_sink(
      "TargetMass", FLOAT, 1, UNIT_CONV_MASS, 0.f, sinks,
      convert_sink_target_mass, "Sink target mass to spawn star particles");

  list[8] = io_make_physical_output_field(
      "Nstars", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sinks, n_stars,
      /*can convert to comoving=*/0,
      "Number of stars spawned by the sink particles");

  /* Note: Since the swallowed momentum is computed with the physical velocity,
     i.e. including the Hubble flow term, it is not convertible to comoving
     frame. */
  list[9] = io_make_physical_output_field_convert_sink(
      "SwallowedAngularMomentum", FLOAT, 3, UNIT_CONV_ANGULAR_MOMENTUM, 0.f,
      sinks,
      /*can convert to comoving=*/0, convert_sink_swallowed_angular_momentum,
      "Physical swallowed angular momentum of the particles");

  if (with_cosmology) {
    list[10] = io_make_physical_output_field(
        "BirthScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sinks,
        birth_scale_factor, /*can convert to comoving=*/0,
        "Scale-factors at which the sinks were born");
  } else {
    list[10] =
        io_make_output_field("BirthTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, sinks,
                             birth_time, "Times at which the sinks were born");
  }

#ifdef DEBUG_INTERACTIONS_SINKS

  list += *num_fields;
  *num_fields += 4;

  list[0] =
      io_make_output_field("Num_ngb_formation", INT, 1, UNIT_CONV_NO_UNITS, 0.f,
                           sinks, num_ngb_formation, "Number of neighbors");
  list[1] =
      io_make_output_field("Ids_ngb_formation", LONGLONG,
                           MAX_NUM_OF_NEIGHBOURS_SINKS, UNIT_CONV_NO_UNITS, 0.f,
                           sinks, ids_ngbs_formation, "IDs of the neighbors");

  list[2] =
      io_make_output_field("Num_ngb_merger", INT, 1, UNIT_CONV_NO_UNITS, 0.f,
                           sinks, num_ngb_merger, "Number of merger");
  list[3] = io_make_output_field(
      "Ids_ngb_merger", LONGLONG, MAX_NUM_OF_NEIGHBOURS_SINKS,
      UNIT_CONV_NO_UNITS, 0.f, sinks, ids_ngbs_merger, "IDs of the neighbors");

#endif
}

#endif /* SWIFT_GEAR_SINK_IO_H */
