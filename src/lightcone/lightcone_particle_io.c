/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <hdf5.h>

/* This object's header. */
#include "lightcone/lightcone_particle_io.h"

/* Local headers */
#include "black_holes.h"
#include "chemistry.h"
#include "chemistry_struct.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "lightcone/lightcone.h"
#include "neutrino.h"
#include "particle_buffer.h"
#include "stars.h"

void lightcone_io_field_list_init(struct lightcone_io_field_list *list) {

  list->first = NULL;
  list->last = NULL;
  list->num_fields = 0;
}

void lightcone_io_field_list_clean(struct lightcone_io_field_list *list) {

  struct lightcone_io_field *current;
  struct lightcone_io_field *next;

  current = list->first;
  while (current) {
    next = current->next;
    free(current);
    current = next;
  }

  list->first = NULL;
  list->last = NULL;
  list->num_fields = 0;
}

void lightcone_io_field_list_append(struct lightcone_io_field_list *list,
                                    char *name, enum IO_DATA_TYPE type,
                                    int dimension, size_t offset,
                                    enum unit_conversion_factor units,
                                    float scale_factor_exponent,
                                    char *compression) {

  /* Make the new lightcone_io_field struct */
  struct lightcone_io_field *r =
      (struct lightcone_io_field *)malloc(sizeof(struct lightcone_io_field));
  bzero(r, sizeof(struct lightcone_io_field));
  strcpy(r->name, name);
  r->type = type;
  r->dimension = dimension;
  r->offset = offset;
  r->units = units;
  r->scale_factor_exponent = scale_factor_exponent;
  r->compression = compression_scheme_from_name(compression);
  r->next = NULL;

  /* Append to the linked list */
  if (list->last) {
    list->last->next = r;
  } else {
    list->first = r;
  }
  list->last = r;
  list->num_fields += 1;
}

/**
 * @brief Make a linked list of output fields for gas particles
 */
void lightcone_io_append_gas_output_fields(
    struct lightcone_io_field_list *list) {

#define OFFSET(x) offsetof(struct lightcone_gas_data, x)
  lightcone_io_field_list_append(list, "ParticleIDs", LONGLONG, 1, OFFSET(id),
                                 UNIT_CONV_NO_UNITS, 0.0, "Nbit40");
  lightcone_io_field_list_append(list, "Coordinates", DOUBLE, 3, OFFSET(x),
                                 UNIT_CONV_LENGTH, 1.0, "DScale5");
  lightcone_io_field_list_append(list, "Velocities", FLOAT, 3, OFFSET(vel),
                                 UNIT_CONV_SPEED, 0.0, "DScale1");
  lightcone_io_field_list_append(list, "Masses", FLOAT, 1, OFFSET(mass),
                                 UNIT_CONV_MASS, 0.0, "on");
  lightcone_io_field_list_append(list, "ExpansionFactors", FLOAT, 1, OFFSET(a),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "SmoothingLengths", FLOAT, 1, OFFSET(h),
                                 UNIT_CONV_LENGTH, 0.0, "on");
  lightcone_io_field_list_append(list, "Densities", FLOAT, 1, OFFSET(rho),
                                 UNIT_CONV_DENSITY, -3.0, "FMantissa9");
  lightcone_io_field_list_append(list, "Temperatures", FLOAT, 1,
                                 OFFSET(temperature), UNIT_CONV_TEMPERATURE,
                                 0.0, "FMantissa9");
#ifdef CHEMISTRY_EAGLE
  lightcone_io_field_list_append(list, "SmoothedElementMassFractions", FLOAT,
                                 chemistry_element_count,
                                 OFFSET(smoothed_metal_mass_fraction),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
  lightcone_io_field_list_append(list, "SmoothedMetalMassFractions", FLOAT, 1,
                                 OFFSET(smoothed_metal_mass_fraction_total),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
  lightcone_io_field_list_append(list, "MetalMassFractions", FLOAT, 1,
                                 OFFSET(metal_mass_fraction_total),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
#endif
#ifdef COOLING_PS2020
  lightcone_io_field_list_append(list, "ElectronNumberDensities", DOUBLE, 1,
                                 OFFSET(electron_density),
                                 UNIT_CONV_NUMBER_DENSITY, 0.0, "DMantissa9");
  lightcone_io_field_list_append(list, "ComptonYParameters", DOUBLE, 1,
                                 OFFSET(ycompton), UNIT_CONV_AREA, 0.0,
                                 "DMantissa9");
#endif
#ifdef WITH_FOF
  lightcone_io_field_list_append(list, "FOFGroupIDs", LONGLONG, 1,
                                 OFFSET(group_id), UNIT_CONV_NO_UNITS, 0.0,
                                 "on");
#endif
#if defined(TRACERS_EAGLE) || defined(TRACERS_FLAMINGO)
  lightcone_io_field_list_append(list, "LastAGNFeedbackScaleFactors", FLOAT, 1,
                                 OFFSET(last_AGN_injection_scale_factor),
                                 UNIT_CONV_NO_UNITS, 0.0, "BFloat16");
#endif
#ifdef STAR_FORMATION_EAGLE
  lightcone_io_field_list_append(list, "StarFormationRates", FLOAT, 1,
                                 OFFSET(sfr), UNIT_CONV_SFR, 0.0, "on");
#endif
#undef OFFSET
}

/**
 * @brief Make a linked list of output fields for DM particles
 */
void lightcone_io_append_dark_matter_output_fields(
    struct lightcone_io_field_list *list) {

#define OFFSET(x) offsetof(struct lightcone_dark_matter_data, x)
  lightcone_io_field_list_append(list, "ParticleIDs", LONGLONG, 1, OFFSET(id),
                                 UNIT_CONV_NO_UNITS, 0.0, "Nbit40");
  lightcone_io_field_list_append(list, "Coordinates", DOUBLE, 3, OFFSET(x),
                                 UNIT_CONV_LENGTH, 1.0, "DScale5");
  lightcone_io_field_list_append(list, "Velocities", FLOAT, 3, OFFSET(vel),
                                 UNIT_CONV_SPEED, 0.0, "DScale1");
  lightcone_io_field_list_append(list, "Masses", FLOAT, 1, OFFSET(mass),
                                 UNIT_CONV_MASS, 0.0, "on");
  lightcone_io_field_list_append(list, "ExpansionFactors", FLOAT, 1, OFFSET(a),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
#undef OFFSET
}

/**
 * @brief Make a linked list of output fields for DM background particles
 */
void lightcone_io_append_dark_matter_background_output_fields(
    struct lightcone_io_field_list *list) {

#define OFFSET(x)                             \
  offsetof(struct lightcone_dark_matter_data, \
           x) /* Uses same struct as dark matter */
  lightcone_io_field_list_append(list, "ParticleIDs", LONGLONG, 1, OFFSET(id),
                                 UNIT_CONV_NO_UNITS, 0.0, "Nbit40");
  lightcone_io_field_list_append(list, "Coordinates", DOUBLE, 3, OFFSET(x),
                                 UNIT_CONV_LENGTH, 1.0, "DScale5");
  lightcone_io_field_list_append(list, "Velocities", FLOAT, 3, OFFSET(vel),
                                 UNIT_CONV_SPEED, 0.0, "DScale1");
  lightcone_io_field_list_append(list, "Masses", FLOAT, 1, OFFSET(mass),
                                 UNIT_CONV_MASS, 0.0, "on");
  lightcone_io_field_list_append(list, "ExpansionFactors", FLOAT, 1, OFFSET(a),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
#undef OFFSET
}

/**
 * @brief Make a linked list of output fields for star particles
 */
void lightcone_io_append_stars_output_fields(
    struct lightcone_io_field_list *list) {

#define OFFSET(x) offsetof(struct lightcone_stars_data, x)
  lightcone_io_field_list_append(list, "ParticleIDs", LONGLONG, 1, OFFSET(id),
                                 UNIT_CONV_NO_UNITS, 0.0, "Nbit40");
  lightcone_io_field_list_append(list, "Coordinates", DOUBLE, 3, OFFSET(x),
                                 UNIT_CONV_LENGTH, 1.0, "DScale5");
  lightcone_io_field_list_append(list, "Velocities", FLOAT, 3, OFFSET(vel),
                                 UNIT_CONV_SPEED, 0.0, "DScale1");
  lightcone_io_field_list_append(list, "Masses", FLOAT, 1, OFFSET(mass),
                                 UNIT_CONV_MASS, 0.0, "on");
  lightcone_io_field_list_append(list, "ExpansionFactors", FLOAT, 1, OFFSET(a),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
#ifdef WITH_FOF
  lightcone_io_field_list_append(list, "FOFGroupIDs", LONGLONG, 1,
                                 OFFSET(group_id), UNIT_CONV_NO_UNITS, 0.0,
                                 "on");
#endif
#ifdef STARS_EAGLE
  lightcone_io_field_list_append(list, "InitialMasses", FLOAT, 1,
                                 OFFSET(mass_init), UNIT_CONV_MASS, 0.0,
                                 "FMantissa9");
  lightcone_io_field_list_append(list, "BirthScaleFactors", FLOAT, 1,
                                 OFFSET(birth_scale_factor), UNIT_CONV_NO_UNITS,
                                 0.0, "FMantissa9");
  lightcone_io_field_list_append(list, "BirthDensities", FLOAT, 1,
                                 OFFSET(birth_density), UNIT_CONV_DENSITY, 0.0,
                                 "BFloat16");
  lightcone_io_field_list_append(list, "Luminosities", FLOAT,
                                 luminosity_bands_count, OFFSET(luminosities),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
#endif
#ifdef CHEMISTRY_EAGLE
  lightcone_io_field_list_append(list, "SmoothedElementMassFractions", FLOAT,
                                 chemistry_element_count,
                                 OFFSET(smoothed_metal_mass_fraction),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
  lightcone_io_field_list_append(list, "SmoothedMetalMassFractions", FLOAT, 1,
                                 OFFSET(smoothed_metal_mass_fraction_total),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
  lightcone_io_field_list_append(list, "MetalMassFractions", FLOAT, 1,
                                 OFFSET(metal_mass_fraction_total),
                                 UNIT_CONV_NO_UNITS, 0.0, "FMantissa9");
#endif
#if defined(TRACERS_EAGLE) || defined(TRACERS_FLAMINGO)
  lightcone_io_field_list_append(list, "LastAGNFeedbackScaleFactors", FLOAT, 1,
                                 OFFSET(last_AGN_injection_scale_factor),
                                 UNIT_CONV_NO_UNITS, 0.0, "BFloat16");
#endif
#undef OFFSET
}

/**
 * @brief Make a linked list of output fields for black hole particles
 */
void lightcone_io_append_black_hole_output_fields(
    struct lightcone_io_field_list *list) {

#define OFFSET(x) offsetof(struct lightcone_black_hole_data, x)
  lightcone_io_field_list_append(list, "ParticleIDs", LONGLONG, 1, OFFSET(id),
                                 UNIT_CONV_NO_UNITS, 0.0, "Nbit40");
  lightcone_io_field_list_append(list, "Coordinates", DOUBLE, 3, OFFSET(x),
                                 UNIT_CONV_LENGTH, 1.0, "DScale5");
  lightcone_io_field_list_append(list, "Velocities", FLOAT, 3, OFFSET(vel),
                                 UNIT_CONV_SPEED, 0.0, "DScale1");
  lightcone_io_field_list_append(list, "DynamicalMasses", FLOAT, 1,
                                 OFFSET(mass), UNIT_CONV_MASS, 0.0, "on");
  lightcone_io_field_list_append(list, "ExpansionFactors", FLOAT, 1, OFFSET(a),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
#ifdef BLACK_HOLES_EAGLE
  lightcone_io_field_list_append(list, "SubgridMasses", FLOAT, 1,
                                 OFFSET(subgrid_mass), UNIT_CONV_MASS, 0.0,
                                 "on");
  lightcone_io_field_list_append(list, "FormationScaleFactors", FLOAT, 1,
                                 OFFSET(formation_scale_factor),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "AccretionRates", FLOAT, 1,
                                 OFFSET(accretion_rate),
                                 UNIT_CONV_MASS_PER_UNIT_TIME, 0.0, "on");
  lightcone_io_field_list_append(list, "TotalAccretedMasses", FLOAT, 1,
                                 OFFSET(total_accreted_mass), UNIT_CONV_MASS,
                                 0.0, "on");
  lightcone_io_field_list_append(list, "LastMinorMergerScaleFactors", FLOAT, 1,
                                 OFFSET(last_minor_merger_scale_factor),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "LastMajorMergerScaleFactors", FLOAT, 1,
                                 OFFSET(last_major_merger_scale_factor),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "NumberOfMergers", INT, 1,
                                 OFFSET(number_of_mergers), UNIT_CONV_NO_UNITS,
                                 0.0, "on");
  lightcone_io_field_list_append(list, "LastAGNFeedbackScaleFactors", FLOAT, 1,
                                 OFFSET(last_AGN_event_scale_factor),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "NumberOfAGNEvents", INT, 1,
                                 OFFSET(AGN_number_of_AGN_events),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "NumberOfHeatingEvents", INT, 1,
                                 OFFSET(AGN_number_of_energy_injections),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(
      list, "LastHighEddingtonFractionScaleFactors", FLOAT, 1,
      OFFSET(last_high_Eddington_fraction_scale_factor), UNIT_CONV_NO_UNITS,
      0.0, "on");
  lightcone_io_field_list_append(list, "CumulativeNumberOfSeeds", INT, 1,
                                 OFFSET(cumulative_number_seeds),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
#ifdef WITH_FOF
  lightcone_io_field_list_append(list, "FOFGroupIDs", LONGLONG, 1,
                                 OFFSET(group_id), UNIT_CONV_NO_UNITS, 0.0,
                                 "on");
#endif
#endif
#undef OFFSET
}

/**
 * @brief Make a linked list of output fields for neutrino particles
 */
void lightcone_io_append_neutrino_output_fields(
    struct lightcone_io_field_list *list) {

#define OFFSET(x) offsetof(struct lightcone_neutrino_data, x)
  lightcone_io_field_list_append(list, "ParticleIDs", LONGLONG, 1, OFFSET(id),
                                 UNIT_CONV_NO_UNITS, 0.0, "Nbit40");
  lightcone_io_field_list_append(list, "Coordinates", DOUBLE, 3, OFFSET(x),
                                 UNIT_CONV_LENGTH, 1.0, "DScale5");
  lightcone_io_field_list_append(list, "Velocities", FLOAT, 3, OFFSET(vel),
                                 UNIT_CONV_SPEED, 0.0, "DScale1");
  lightcone_io_field_list_append(list, "Masses", FLOAT, 1, OFFSET(mass),
                                 UNIT_CONV_MASS, 0.0, "on");
  lightcone_io_field_list_append(list, "Weights", FLOAT, 1, OFFSET(weight),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
  lightcone_io_field_list_append(list, "ExpansionFactors", FLOAT, 1, OFFSET(a),
                                 UNIT_CONV_NO_UNITS, 0.0, "on");
#undef OFFSET
}

/*
  Functions to store particle properties in the lightcone_*_data structs.

  These should determine whether the particle should be included in the
  lightcone and, if so, copy the needed quantities into the struct and
  return 1. If the particle should be discarded the function should
  return 0.

 */

/**
 * @brief Store gas properties to write to the lightcone
 *
 * If the particle should be included in the lightcone output this function
 * copies its information to the lightcone_gas_data struct and returns 1.
 * If the particle should not be output the function returns 0.
 *
 * @param e the #engine structure
 * @param gp the #gpart which crossed the lightcone
 * @param p the #part associated with this #gpart
 * @param xp the #xpart associated with this #gpart
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 * @param the #lightcone_gas_data struct to update
 */
int lightcone_store_gas(const struct engine *e, struct lightcone_props *props,
                        const struct gpart *gp, const struct part *p,
                        const struct xpart *xp, const double a_cross,
                        const double x_cross[3],
                        struct lightcone_gas_data *data) {

  /*! Check if we're filtering gas particles */
  if (props->gas_filtering_enabled) {
    if (a_cross < props->max_a_for_gas_filtering) {

      /* Check hydrogen number density of this particle */
#ifdef CHEMISTRY_EAGLE
      const double density = p->rho;
      const double proton_mass = e->physical_constants->const_proton_mass;
      const double hydrogen_fraction =
          p->chemistry_data.metal_mass_fraction[chemistry_element_H];
      const double nh = density * hydrogen_fraction / proton_mass;
      if (nh < props->min_nh_for_filtered_gas * pow(a_cross, -4.0)) return 0;
#else
      error(
          "Lightcone gas particle filtering is only implemented for EAGLE "
          "chemistry");
#endif
      /* Check temperature of this particle */
      const double T = cooling_get_temperature(
          e->physical_constants, e->hydro_properties, e->internal_units,
          e->cosmology, e->cooling_func, p, xp);
      if (T < props->min_temp_for_filtered_gas) return 0;
    }
  }

  data->id = p->id;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->vel[0] =
      xp->v_full[0] / a_cross;  // TODO: extrapolate velocities to a_cross?
  data->vel[1] = xp->v_full[1] / a_cross;
  data->vel[2] = xp->v_full[2] / a_cross;
  data->mass = hydro_get_mass(p);
  data->a = a_cross;
  data->h = p->h;
  data->rho = hydro_get_comoving_density(p);
  data->temperature = cooling_get_temperature(
      e->physical_constants, e->hydro_properties, e->internal_units,
      e->cosmology, e->cooling_func, p, xp);
#ifdef WITH_FOF
  data->group_id = (long long)gp->fof_data.group_id;
#endif

#ifdef CHEMISTRY_EAGLE
  for (int i = 0; i < chemistry_element_count; i += 1)
    data->smoothed_metal_mass_fraction[i] =
        p->chemistry_data.smoothed_metal_mass_fraction[i];
  data->metal_mass_fraction_total = p->chemistry_data.metal_mass_fraction_total;
  data->smoothed_metal_mass_fraction_total =
      p->chemistry_data.smoothed_metal_mass_fraction_total;
#endif

#ifdef COOLING_PS2020
  data->electron_density = cooling_get_electron_density(
      e->physical_constants, e->hydro_properties, e->internal_units,
      e->cosmology, e->cooling_func, p, xp);
  data->ycompton = cooling_get_ycompton(e->physical_constants,
                                        e->hydro_properties, e->internal_units,
                                        e->cosmology, e->cooling_func, p, xp);
#endif

#if defined(TRACERS_EAGLE) || defined(TRACERS_FLAMINGO)
  data->last_AGN_injection_scale_factor =
      xp->tracers_data.last_AGN_injection_scale_factor;
#endif

#ifdef STAR_FORMATION_EAGLE
  data->sfr = xp->sf_data.SFR;
#endif

  return 1;
}

/**
 * @brief Store dark matter properties to write to the lightcone
 *
 * If the particle should be included in the lightcone output this function
 * copies its information to the lightcone_dark_matter_data struct and returns
 * 1. If the particle should not be output the function returns 0.
 *
 * @param e the #engine structure
 * @param gp the #gpart which crossed the lightcone
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 * @param the #lightcone_dark_matter_data struct to update
 */
int lightcone_store_dark_matter(const struct engine *e,
                                struct lightcone_props *props,
                                const struct gpart *gp, const double a_cross,
                                const double x_cross[3],
                                struct lightcone_dark_matter_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->vel[0] =
      gp->v_full[0] / a_cross;  // TODO: extrapolate velocities to a_cross?
  data->vel[1] = gp->v_full[1] / a_cross;
  data->vel[2] = gp->v_full[2] / a_cross;
  data->mass = gp->mass;
  data->a = a_cross;

  return 1;
}

/**
 * @brief Store star properties to write to the lightcone
 *
 * If the particle should be included in the lightcone output this function
 * copies its information to the lightcone_star_data struct and returns
 * 1. If the particle should not be output the function returns 0.
 *
 * @param e the #engine structure
 * @param gp the #gpart which crossed the lightcone
 * @param sp the #spart associated with the #gpart
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 * @param the #lightcone_stars_data struct to update
 */
int lightcone_store_stars(const struct engine *e, struct lightcone_props *props,
                          const struct gpart *gp, const struct spart *sp,
                          const double a_cross, const double x_cross[3],
                          struct lightcone_stars_data *data) {
  data->id = sp->id;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->vel[0] =
      sp->v[0] / a_cross;  // TODO: extrapolate velocities to a_cross?
  data->vel[1] = sp->v[1] / a_cross;
  data->vel[2] = sp->v[2] / a_cross;
  data->mass = sp->mass;
  data->a = a_cross;

#ifdef WITH_FOF
  data->group_id = (long long)gp->fof_data.group_id;
#endif

#ifdef STARS_EAGLE
  data->mass_init = sp->mass_init;
  data->birth_scale_factor = sp->birth_scale_factor;
  data->birth_density = sp->birth_density;
  stars_get_luminosities(sp, e->policy & engine_policy_cosmology, e->cosmology,
                         e->time, e->physical_constants, e->stars_properties,
                         data->luminosities);
#endif

#ifdef CHEMISTRY_EAGLE
  for (int i = 0; i < chemistry_element_count; i += 1)
    data->smoothed_metal_mass_fraction[i] =
        sp->chemistry_data.smoothed_metal_mass_fraction[i];
  data->metal_mass_fraction_total =
      sp->chemistry_data.metal_mass_fraction_total;
  data->smoothed_metal_mass_fraction_total =
      sp->chemistry_data.smoothed_metal_mass_fraction_total;
#endif

#if defined(TRACERS_EAGLE) || defined(TRACERS_FLAMINGO)
  data->last_AGN_injection_scale_factor =
      sp->tracers_data.last_AGN_injection_scale_factor;
#endif

  return 1;
}

/**
 * @brief Store black hole properties to write to the lightcone
 *
 * If the particle should be included in the lightcone output this function
 * copies its information to the lightcone_black_hole_data struct and returns
 * 1. If the particle should not be output the function returns 0.
 *
 * @param e the #engine structure
 * @param gp the #gpart which crossed the lightcone
 * @param bp the #bpart associated with the #gpart
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 * @param the #lightcone_black_hole_data struct to update
 */
int lightcone_store_black_hole(const struct engine *e,
                               struct lightcone_props *props,
                               const struct gpart *gp, const struct bpart *bp,
                               const double a_cross, const double x_cross[3],
                               struct lightcone_black_hole_data *data) {
  data->id = bp->id;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->vel[0] =
      bp->v[0] / a_cross;  // TODO: extrapolate velocities to a_cross?
  data->vel[1] = bp->v[1] / a_cross;
  data->vel[2] = bp->v[2] / a_cross;
  data->mass = bp->mass;
  data->a = a_cross;
#ifdef BLACK_HOLES_EAGLE
  data->subgrid_mass = bp->subgrid_mass;
  data->formation_scale_factor = bp->formation_scale_factor;
  data->accretion_rate = bp->accretion_rate;
  data->total_accreted_mass = bp->total_accreted_mass;
  data->last_minor_merger_scale_factor = bp->last_minor_merger_scale_factor;
  data->last_major_merger_scale_factor = bp->last_major_merger_scale_factor;
  data->number_of_mergers = bp->number_of_mergers;
  data->last_AGN_event_scale_factor = bp->last_AGN_event_scale_factor;
  data->AGN_number_of_AGN_events = bp->AGN_number_of_AGN_events;
  data->AGN_number_of_energy_injections = bp->AGN_number_of_energy_injections;
  data->last_high_Eddington_fraction_scale_factor =
      bp->last_high_Eddington_fraction_scale_factor;
  data->cumulative_number_seeds = bp->cumulative_number_seeds;
#ifdef WITH_FOF
  data->group_id = (long long)gp->fof_data.group_id;
#endif
#endif
  return 1;
}

/**
 * @brief Store neutrino properties to write to the lightcone
 *
 * If the particle should be included in the lightcone output this function
 * copies its information to the lightcone_neutrino_data struct and returns
 * 1. If the particle should not be output the function returns 0.
 *
 * @param e the #engine structure
 * @param gp the #gpart which crossed the lightcone
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 * @param the #lightcone_neutrino_data struct to update
 */
int lightcone_store_neutrino(const struct engine *e,
                             struct lightcone_props *props,
                             const struct gpart *gp, const double a_cross,
                             const double x_cross[3],
                             struct lightcone_neutrino_data *data) {

  /* Compute neutrino weight */
  struct neutrino_model nu_model;
  bzero(&nu_model, sizeof(struct neutrino_model));
  if (e->neutrino_properties->use_delta_f_mesh_only)
    gather_neutrino_consts(e->s, &nu_model);
  double weight = 1.0;
  gpart_neutrino_weight_mesh_only(gp, &nu_model, &weight);

  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->vel[0] =
      gp->v_full[0] / a_cross;  // TODO: extrapolate velocities to a_cross?
  data->vel[1] = gp->v_full[1] / a_cross;
  data->vel[2] = gp->v_full[2] / a_cross;
  data->mass = gp->mass;
  data->weight = weight;
  data->a = a_cross;

  return 1;
}

/**
 * @brief Write data to a HDF5 dataset, appending along first axis if it already
 * exists
 */
void append_dataset(const struct unit_system *snapshot_units,
                    enum unit_conversion_factor units,
                    float scale_factor_exponent, hid_t loc_id, const char *name,
                    hid_t mem_type_id, hsize_t chunk_size,
                    int lossy_compression,
                    enum lossy_compression_schemes compression_scheme,
                    int gzip_level, const int rank, const hsize_t dims[2],
                    const hsize_t num_written, const void *data) {

  if (rank > 2) error("HDF5 dataset has too may dimensions.");
  if (rank < 1) error("HDF5 dataset must be at least one dimensional");

  /* If we have zero elements to append, there's nothing to do */
  if (dims[0] == 0) return;

  /* Determine size of the dataset after we append our data */
  hsize_t full_dims[2];
  for (int i = 0; i < rank; i += 1) full_dims[i] = dims[i];
  full_dims[0] += num_written;

  /* Determine maximum size in each dimension */
  hsize_t max_dims[2];
  for (int i = 1; i < rank; i += 1) max_dims[i] = full_dims[i];
  max_dims[0] = H5S_UNLIMITED;

  /* Determine chunk size in each dimension */
  hsize_t chunk_dims[2];
  for (int i = 1; i < rank; i += 1) chunk_dims[i] = full_dims[i];
  chunk_dims[0] = (hsize_t)chunk_size;

  /* Find offset to region to write in each dimension */
  hsize_t offset[2];
  for (int i = 1; i < rank; i += 1) offset[i] = 0;
  offset[0] = num_written;

  hid_t dataset_id;
  hid_t file_space_id;
  if (num_written == 0) {

    /* We need to create a new dataset */
    file_space_id = H5Screate_simple(rank, full_dims, max_dims);
    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);

    /* Type of the dataset to create - this is initially the same as the type
       in memory but may be modified by lossy compression. */
    hid_t file_type_id = H5Tcopy(mem_type_id);

    /* Set chunk size and lossy compression scheme, if any  */
    H5Pset_chunk(prop_id, rank, chunk_dims);
    char filter_name[32];
    if (lossy_compression && (compression_scheme != compression_write_lossless))
      set_hdf5_lossy_compression(&prop_id, &file_type_id, compression_scheme,
                                 name, filter_name);

    /* Set lossless compression, if any */
    if (gzip_level > 0) {
      H5Pset_shuffle(prop_id);
      H5Pset_deflate(prop_id, gzip_level);
    }

    /* Create the dataset */
    dataset_id = H5Dcreate(loc_id, name, file_type_id, file_space_id,
                           H5P_DEFAULT, prop_id, H5P_DEFAULT);
    if (dataset_id < 0) error("Failed to create new dataset: %s", name);
    H5Pclose(prop_id);
    H5Tclose(file_type_id);

    /* Write unit conversion factors for this data set */
    char buffer[FIELD_BUFFER_SIZE] = {0};
    units_cgs_conversion_string(buffer, snapshot_units, units,
                                scale_factor_exponent);
    float baseUnitsExp[5];
    units_get_base_unit_exponents_array(baseUnitsExp, units);
    io_write_attribute_f(dataset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
    io_write_attribute_f(dataset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
    io_write_attribute_f(dataset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
    io_write_attribute_f(dataset_id, "U_I exponent",
                         baseUnitsExp[UNIT_CURRENT]);
    io_write_attribute_f(dataset_id, "U_T exponent",
                         baseUnitsExp[UNIT_TEMPERATURE]);
    io_write_attribute_f(dataset_id, "h-scale exponent", 0.f);
    io_write_attribute_f(dataset_id, "a-scale exponent", scale_factor_exponent);
    io_write_attribute_s(dataset_id, "Expression for physical CGS units",
                         buffer);

    /* Write the actual number this conversion factor corresponds to */
    const double factor = units_cgs_conversion_factor(snapshot_units, units);
    io_write_attribute_d(
        dataset_id,
        "Conversion factor to CGS (not including cosmological corrections)",
        factor);

    /* Note that we can't write the conversion factor including cosmological
       corrections as an attribute because it will be different for each
       particle. */

  } else {

    /* We're appending to an existing dataset */
    dataset_id = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (dataset_id < 0) error("Failed to open existing dataset: %s", name);
    if (H5Dset_extent(dataset_id, full_dims) < 0)
      error("Unable to extend dataset: %s", name);
    file_space_id = H5Dget_space(dataset_id);
  }

  /* Create memory dataspace */
  hid_t mem_space_id = H5Screate_simple(rank, dims, NULL);

  /* Select region to write in the file */
  if (H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, dims,
                          NULL) < 0)
    error("Failed to select region in dataset: %s", name);

  /* Write the data */
  if (H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id,
               H5P_DEFAULT, data) < 0)
    error("Failed to write dataset: %s", name);

  /* Clean up*/
  H5Sclose(file_space_id);
  H5Sclose(mem_space_id);
  H5Dclose(dataset_id);
}

hid_t init_write(struct lightcone_props *props, hid_t file_id, int ptype,
                 size_t *num_written, size_t *num_to_write) {

  /* Number of particles already written to the file */
  *num_written = props->num_particles_written_to_file[ptype];

  /* Number of buffered particles */
  *num_to_write = particle_buffer_num_elements(&props->buffer[ptype]);

  /* Create or open the HDF5 group for this particle type */
  const char *name = part_type_names[ptype];
  hid_t group_id;
  if (*num_written > 0) {
    group_id = H5Gopen(file_id, name, H5P_DEFAULT);
    if (group_id < 0) error("Failed to open existing group: %s", name);
  } else {
    group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) error("Failed to create new group: %s", name);
  }
  return group_id;
}

/**
 * @brief Append buffered particles to the output file.
 */
void lightcone_write_particles(struct lightcone_props *props,
                               const struct unit_system *internal_units,
                               const struct unit_system *snapshot_units,
                               int ptype, hid_t file_id) {

  if (props->particle_fields[ptype].num_fields > 0) {

    /* Open group and get number and offset of particles to write */
    size_t num_written, num_to_write;
    hid_t group_id =
        init_write(props, file_id, ptype, &num_written, &num_to_write);

    /* Get size of the data struct for this type */
    const size_t data_struct_size = lightcone_io_struct_size(ptype);

    /* Loop over output fields */
    struct lightcone_io_field *f = props->particle_fields[ptype].first;
    while (f) {

      /* Find output field info */
      hid_t dtype_id = io_hdf5_type(f->type);     /* HDF5 data type */
      size_t type_size = io_sizeof_type(f->type); /* Bytes per value */
      const size_t field_size =
          f->dimension * type_size; /* Bytes per particle */
      const enum lossy_compression_schemes compression_scheme =
          f->compression; /* Compression scheme */

      /* Find unit conversion factor for this quantity */
      const double conversion_factor =
          units_conversion_factor(internal_units, snapshot_units, f->units);

      /* Allocate output buffer */
      char *outbuf = (char *)malloc(num_to_write * field_size);
      if (!outbuf) error("Unable to allocate lightcone output buffer");
      char *outptr = outbuf;

      /* Loop over blocks of buffered particles and copy to output array */
      size_t num_elements;
      struct particle_buffer_block *block = NULL;
      char *block_data;
      do {
        particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements,
                                (void **)&block_data);
        for (size_t i = 0; i < num_elements; i += 1) {
          char *src = block_data + i * data_struct_size + f->offset;
          char *dest = outptr;
          memcpy(dest, src, field_size);
          outptr += field_size;
        }
      } while (block);

      /* Convert units if necessary */
      if (conversion_factor != 1.0) {
        const size_t nr_values = num_to_write * f->dimension;
        switch (f->type) {
          case INT: {
            int *values = (int *)outbuf;
            for (size_t i = 0; i < nr_values; i += 1)
              values[i] *= conversion_factor;
          } break;
          case LONGLONG: {
            long long *values = (long long *)outbuf;
            for (size_t i = 0; i < nr_values; i += 1)
              values[i] *= conversion_factor;
          } break;
          case FLOAT: {
            float *values = (float *)outbuf;
            for (size_t i = 0; i < nr_values; i += 1)
              values[i] *= conversion_factor;
          } break;
          case DOUBLE: {
            double *values = (double *)outbuf;
            for (size_t i = 0; i < nr_values; i += 1)
              values[i] *= conversion_factor;
          } break;
          default:
            error("Unhandled data type");
        }
      }

      /* Write the data */
      const hsize_t chunk_size = props->hdf5_chunk_size;
      hsize_t dims[] = {(hsize_t)num_to_write, (hsize_t)f->dimension};
      int rank = 1;
      if (f->dimension > 1) rank = 2;
      append_dataset(snapshot_units, f->units, f->scale_factor_exponent,
                     group_id, f->name, dtype_id, chunk_size,
                     props->particles_lossy_compression, compression_scheme,
                     props->particles_gzip_level, rank, dims, num_written,
                     outbuf);

      /* Free the output buffer */
      free(outbuf);

      /* Advance to next output field */
      f = f->next;
    }

    /* If all fields are done, we can close the particle type group */
    H5Gclose(group_id);
  }
}
