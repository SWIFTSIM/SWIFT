/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STARS_IO_H
#define SWIFT_EAGLE_STARS_IO_H

#include "io_properties.h"
#include "stars_part.h"

/**
 * @brief Specifies which s-particle fields to read from a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void stars_read_particles(struct spart *sparts,
                                        struct io_props *list,
                                        int *num_fields) {

  /* Say how much we want to read */
  *num_fields = 9;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, sparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, sparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, sparts, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, sparts, h);
  list[5] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass_init);
  list[6] =
      io_make_input_field_default("StellarFormationTime", FLOAT, 1, OPTIONAL,
                                  UNIT_CONV_NO_UNITS, sparts, birth_time, -1.);
  list[7] = io_make_input_field("BirthDensities", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, sparts, birth_density);
  list[8] =
      io_make_input_field("BirthTemperatures", FLOAT, 1, OPTIONAL,
                          UNIT_CONV_TEMPERATURE, sparts, birth_temperature);
}

INLINE static void convert_spart_pos(const struct engine *e,
                                     const struct spart *sp, double *ret) {

  const struct space *s = e->s;
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

INLINE static void convert_spart_vel(const struct engine *e,
                                     const struct spart *sp, float *ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology *cosmo = e->cosmology;
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
  const struct gpart *gp = sp->gpart;
  ret[0] = gp->v_full[0] + gp->a_grav[0] * dt_kick_grav;
  ret[1] = gp->v_full[1] + gp->a_grav[1] * dt_kick_grav;
  ret[2] = gp->v_full[2] + gp->a_grav[2] * dt_kick_grav;

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

/**
 * @brief Specifies which s-particle fields to write to a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * @param with_cosmology Are we running a cosmological simulation?
 */
INLINE static void stars_write_particles(const struct spart *sparts,
                                         struct io_props *list, int *num_fields,
                                         const int with_cosmology) {

  /* Say how much we want to write */
  *num_fields = 12;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_spart(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, sparts,
      convert_spart_pos, "Co-moving position of the particles");

  list[1] = io_make_output_field_convert_spart(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, sparts, convert_spart_vel,
      "Peculiar velocities of the particles. This is a * dx/dt where x is the "
      "co-moving position of the particles.");

  list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 sparts, mass,
                                 "Masses of the particles at the current point "
                                 "in time (i.e. after stellar losses");

  list[3] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           sparts, id, "Unique ID of the particles");

  list[4] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sparts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[5] = io_make_output_field("InitialMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 sparts, mass_init,
                                 "Masses of the star particles at birth time");

  if (with_cosmology) {
    list[6] = io_make_output_field(
        "BirthScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        birth_scale_factor, "Scale-factors at which the stars were born");
  } else {
    list[6] = io_make_output_field("BirthTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f,
                                   sparts, birth_time,
                                   "Times at which the stars were born");
  }

  list[7] = io_make_output_field(
      "FeedbackEnergyFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts, f_E,
      "Fractions of the canonical feedback energy that was used for the stars' "
      "SNII feedback events");

  list[8] = io_make_output_field(
      "NumberOfFeedbackEvents", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      number_of_SNII_events,
      "Number of SNII energy injection events the stars went through.");

  list[9] = io_make_output_field(
      "BirthDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts, birth_density,
      "Physical densities at the time of birth of the gas particles that "
      "turned into stars (note that "
      "we store the physical density at the birth redshift, no conversion is "
      "needed)");

  list[10] =
      io_make_output_field("BirthTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE,
                           0.f, sparts, birth_temperature,
                           "Temperatures at the time of birth of the gas "
                           "particles that turned into stars");

  list[11] = io_make_output_field(
      "FeedbackNumberOfHeatingEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      sparts, number_of_heating_events,
      "Expected number of particles that were heated by each star particle.");
}

/**
 * @brief Initialize the global properties of the stars scheme.
 *
 * By default, takes the values provided by the hydro.
 *
 * @param sp The #stars_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param p The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void stars_props_init(struct stars_props *sp,
                                    const struct phys_const *phys_const,
                                    const struct unit_system *us,
                                    struct swift_params *params,
                                    const struct hydro_props *p,
                                    const struct cosmology *cosmo) {

  /* Kernel properties */
  sp->eta_neighbours = parser_get_opt_param_float(
      params, "Stars:resolution_eta", p->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  sp->h_tolerance =
      parser_get_opt_param_float(params, "Stars:h_tolerance", p->h_tolerance);

  /* Get derived properties */
  sp->target_neighbours = pow_dimension(sp->eta_neighbours) * kernel_norm;
  const float delta_eta = sp->eta_neighbours * (1.f + sp->h_tolerance);
  sp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(sp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  sp->max_smoothing_iterations = parser_get_opt_param_int(
      params, "Stars:max_ghost_iterations", p->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "Stars:max_volume_change", -1);
  if (max_volume_change == -1)
    sp->log_max_h_change = p->log_max_h_change;
  else
    sp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Do we want to overwrite the stars' birth properties? */
  sp->overwrite_birth_time =
      parser_get_opt_param_int(params, "Stars:overwrite_birth_time", 0);
  sp->overwrite_birth_density =
      parser_get_opt_param_int(params, "Stars:overwrite_birth_density", 0);
  sp->overwrite_birth_temperature =
      parser_get_opt_param_int(params, "Stars:overwrite_birth_temperature", 0);

  /* Read birth time to set all stars in ICs */
  if (sp->overwrite_birth_time) {
    sp->spart_first_init_birth_time =
        parser_get_param_float(params, "Stars:birth_time");
  }

  /* Read birth density to set all stars in ICs */
  if (sp->overwrite_birth_density) {
    sp->spart_first_init_birth_density =
        parser_get_param_float(params, "Stars:birth_density");
  }

  /* Read birth temperature to set all stars in ICs */
  if (sp->overwrite_birth_temperature) {
    sp->spart_first_init_birth_temperature =
        parser_get_param_float(params, "Stars:birth_temperature");
  }

  /* Maximal time-step lengths */
  const double Myr = 1e6 * 365.25 * 24. * 60. * 60.;
  const double conv_fac = units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  const double max_time_step_young_Myr = parser_get_opt_param_float(
      params, "Stars:max_timestep_young_Myr", FLT_MAX);
  const double max_time_step_old_Myr =
      parser_get_opt_param_float(params, "Stars:max_timestep_old_Myr", FLT_MAX);
  const double age_threshold_Myr = parser_get_opt_param_float(
      params, "Stars:timestep_age_threshold_Myr", FLT_MAX);
  const double age_threshold_unlimited_Myr = parser_get_opt_param_float(
      params, "Stars:timestep_age_threshold_unlimited_Myr", 0.);

  /* Check for consistency */
  if (age_threshold_unlimited_Myr != 0. && age_threshold_Myr != FLT_MAX) {
    if (age_threshold_unlimited_Myr < age_threshold_Myr)
      error(
          "The age threshold for unlimited stellar time-step sizes (%e Myr) is "
          "smaller than the transition threshold from young to old ages (%e "
          "Myr)",
          age_threshold_unlimited_Myr, age_threshold_Myr);
  }

  /* Convert to internal units */
  sp->max_time_step_young = max_time_step_young_Myr * Myr / conv_fac;
  sp->max_time_step_old = max_time_step_old_Myr * Myr / conv_fac;
  sp->age_threshold = age_threshold_Myr * Myr / conv_fac;
  sp->age_threshold_unlimited = age_threshold_unlimited_Myr * Myr / conv_fac;
}

/**
 * @brief Print the global properties of the stars scheme.
 *
 * @param sp The #stars_props.
 */
INLINE static void stars_props_print(const struct stars_props *sp) {

  message("Stars kernel: %s with eta=%f (%.2f neighbours).", kernel_name,
          sp->eta_neighbours, sp->target_neighbours);

  message("Stars relative tolerance in h: %.5f (+/- %.4f neighbours).",
          sp->h_tolerance, sp->delta_neighbours);

  message(
      "Stars integration: Max change of volume: %.2f "
      "(max|dlog(h)/dt|=%f).",
      pow_dimension(expf(sp->log_max_h_change)), sp->log_max_h_change);

  message("Maximal iterations in ghost task set to %d",
          sp->max_smoothing_iterations);

  if (sp->overwrite_birth_time)
    message("Stars' birth time read from the ICs will be overwritten to %f",
            sp->spart_first_init_birth_time);

  message("Stars' age threshold for unlimited dt: %e [U_t]",
          sp->age_threshold_unlimited);
  message("Stars' young/old age threshold: %e [U_t]", sp->age_threshold);
  message("Max time-step size of young stars: %e [U_t]",
          sp->max_time_step_young);
  message("Max time-step size of old stars: %e [U_t]", sp->max_time_step_old);
}

#if defined(HAVE_HDF5)
INLINE static void stars_props_print_snapshot(hid_t h_grpstars,
                                              const struct stars_props *sp) {

  io_write_attribute_s(h_grpstars, "Kernel function", kernel_name);
  io_write_attribute_f(h_grpstars, "Kernel target N_ngb",
                       sp->target_neighbours);
  io_write_attribute_f(h_grpstars, "Kernel delta N_ngb", sp->delta_neighbours);
  io_write_attribute_f(h_grpstars, "Kernel eta", sp->eta_neighbours);
  io_write_attribute_f(h_grpstars, "Smoothing length tolerance",
                       sp->h_tolerance);
  io_write_attribute_f(h_grpstars, "Volume log(max(delta h))",
                       sp->log_max_h_change);
  io_write_attribute_f(h_grpstars, "Volume max change time-step",
                       pow_dimension(expf(sp->log_max_h_change)));
  io_write_attribute_i(h_grpstars, "Max ghost iterations",
                       sp->max_smoothing_iterations);
}
#endif

/**
 * @brief Write a #stars_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
INLINE static void stars_props_struct_dump(const struct stars_props *p,
                                           FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct stars_props), 1, stream,
                       "starsprops", "stars props");
}

/**
 * @brief Restore a stars_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
INLINE static void stars_props_struct_restore(const struct stars_props *p,
                                              FILE *stream) {
  restart_read_blocks((void *)p, sizeof(struct stars_props), 1, stream, NULL,
                      "stars props");
}

#endif /* SWIFT_EAGLE_STAR_IO_H */
