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
  *num_fields = 6;

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
}

INLINE static void convert_bpart_pos(const struct engine* e,
                                     const struct bpart* bp, double* ret) {

  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(bp->x[0] - s->pos_dithering[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(bp->x[1] - s->pos_dithering[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(bp->x[2] - s->pos_dithering[2], 0.0, s->dim[2]);
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

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
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
  *num_fields = 12;

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

  list[3] = io_make_output_field("ParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 0.f, bparts, id, "Unique ID of the particles");

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
      "GasSoundSpeeds", FLOAT, 1, UNIT_CONV_SPEED, 1.5f * hydro_gamma_minus_one,
      bparts, sound_speed_gas,
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
      "CumulativeNumberSeeds", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      cumulative_number_seeds,
      "Total number of BH seeds that have merged into this black hole");

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
