/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_GRACKLE_H
#define SWIFT_COOLING_GRACKLE_H

/**
 * @file src/cooling/none/cooling.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/* Some standard headers. */
#include <float.h>
#include <grackle.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Cooling Model", "Grackle");
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct cooling_function_data* cooling) {

  xp->cooling_data.radiated_energy = 0.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle'.");
  message("Using Grackle           = %i", cooling->chemistry.use_grackle);
  message("Chemical network        = %i",
          cooling->chemistry.primordial_chemistry);
  message("Radiative cooling       = %i",
          cooling->chemistry.with_radiative_cooling);
  message("Metal cooling           = %i", cooling->chemistry.metal_cooling);

  message("CloudyTable             = %s", cooling->cloudy_table);
  message("UVbackground            = %d", cooling->uv_background);
  message("Redshift                = %g", cooling->redshift);
  message("Density Self Shielding  = %g", cooling->density_self_shielding);
  message("Units:");
  message("\tComoving     = %i", cooling->units.comoving_coordinates);
  message("\tLength       = %g", cooling->units.length_units);
  message("\tDensity      = %g", cooling->units.density_units);
  message("\tTime         = %g", cooling->units.time_units);
  message("\tScale Factor = %g", cooling->units.a_units);
}

/**
 * @brief Compute the cooling rate
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param dt The time-step of this particle.
 *
 * @return du / dt
 */
__attribute__((always_inline)) INLINE static double cooling_rate(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, float dt) {

  if (cooling->chemistry.primordial_chemistry > 1) error("Not implemented");

  /* set current time */
  code_units units = cooling->units;
  if (cooling->redshift == -1)
    error("TODO time dependant redshift");
  else
    units.a_value = 1. / (1. + cooling->redshift);

  /* initialize data */
  grackle_field_data data;

  /* set values */
  /* grid */
  int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  int grid_start[GRACKLE_RANK] = {0, 0, 0};
  int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};

  data.grid_rank = GRACKLE_RANK;
  data.grid_dimension = grid_dimension;
  data.grid_start = grid_start;
  data.grid_end = grid_end;

  /* general particle data */
  gr_float density = hydro_get_density(p);
  const double energy_before = hydro_get_internal_energy(p);
  gr_float energy = energy_before;
  gr_float vx = 0;
  gr_float vy = 0;
  gr_float vz = 0;

  data.density = &density;
  data.internal_energy = &energy;
  data.x_velocity = &vx;
  data.y_velocity = &vy;
  data.z_velocity = &vz;

  /* /\* primordial chemistry >= 1 *\/ */
  /* gr_float HI_density = density; */
  /* gr_float HII_density = 0.; */
  /* gr_float HeI_density = 0.; */
  /* gr_float HeII_density = 0.; */
  /* gr_float HeIII_density = 0.; */
  /* gr_float e_density = 0.; */

  /* data.HI_density = &HI_density; */
  /* data.HII_density = &HII_density; */
  /* data.HeI_density = &HeI_density; */
  /* data.HeII_density = &HeII_density; */
  /* data.HeIII_density = &HeIII_density; */
  /* data.e_density = &e_density; */

  /* /\* primordial chemistry >= 2 *\/ */
  /* gr_float HM_density = 0.; */
  /* gr_float H2I_density = 0.; */
  /* gr_float H2II_density = 0.; */

  /* data.HM_density = &HM_density; */
  /* data.H2I_density = &H2I_density; */
  /* data.H2II_density = &H2II_density; */

  /* /\* primordial chemistry >= 3 *\/ */
  /* gr_float DI_density = 0.; */
  /* gr_float DII_density = 0.; */
  /* gr_float HDI_density = 0.; */

  /* data.DI_density = &DI_density; */
  /* data.DII_density = &DII_density; */
  /* data.HDI_density = &HDI_density; */

  /* metal cooling = 1 */
  gr_float metal_density = density * grackle_data->SolarMetalFractionByMass;

  data.metal_density = &metal_density;

  /* /\* volumetric heating rate *\/ */
  /* gr_float volumetric_heating_rate = 0.; */

  /* data.volumetric_heating_rate = &volumetric_heating_rate; */

  /* /\* specific heating rate *\/ */
  /* gr_float specific_heating_rate = 0.; */

  /* data.specific_heating_rate = &specific_heating_rate; */

  /* solve chemistry with table */
  if (solve_chemistry(&units, &data, dt) == 0) {
    error("Error in solve_chemistry.");
  }

  return (energy - energy_before) / dt;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {

  if (dt == 0.) return;

  /* Current du_dt */
  const float hydro_du_dt = hydro_get_internal_energy_dt(p);

  /* compute cooling rate */
  const float du_dt = cooling_rate(phys_const, us, cooling, p, dt);

  /* record energy lost */
  xp->cooling_data.radiated_energy += -du_dt * dt * hydro_get_mass(p);

  /* Update the internal energy */
  hydro_set_internal_energy_dt(p, hydro_du_dt + du_dt);
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const struct part* restrict p) {

  return FLT_MAX;
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_backend(
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {

  /* read parameters */
  parser_get_param_string(parameter_file, "GrackleCooling:GrackleCloudyTable",
                          cooling->cloudy_table);
  cooling->uv_background =
      parser_get_param_int(parameter_file, "GrackleCooling:UVbackground");

  cooling->redshift =
      parser_get_param_double(parameter_file, "GrackleCooling:GrackleRedshift");

  cooling->density_self_shielding = parser_get_param_double(
      parameter_file, "GrackleCooling:GrackleHSShieldingDensityThreshold");

#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

  /* Set up the units system.
     These are conversions from code units to cgs. */

  /* first cosmo */
  cooling->units.a_units = 1.0;  // units for the expansion factor (1/1+zi)
  cooling->units.a_value = 1.0;

  /* We assume here all physical quantities to
     be in proper coordinate (not comobile)  */
  cooling->units.comoving_coordinates = 0;

  /* then units */
  cooling->units.density_units =
      us->UnitMass_in_cgs / pow(us->UnitLength_in_cgs, 3);
  cooling->units.length_units = us->UnitLength_in_cgs;
  cooling->units.time_units = us->UnitTime_in_cgs;
  cooling->units.velocity_units = cooling->units.a_units *
                                  cooling->units.length_units /
                                  cooling->units.time_units;

  chemistry_data* chemistry = &cooling->chemistry;

  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters(chemistry) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry.
  chemistry->use_grackle = 1;
  chemistry->with_radiative_cooling = 1;
  /* molecular network with H, He, D
   From Cloudy table */
  chemistry->primordial_chemistry = 0;
  chemistry->metal_cooling = 1;  // metal cooling on
  chemistry->UVbackground = cooling->uv_background;
  chemistry->grackle_data_file = cooling->cloudy_table;
  chemistry->use_radiative_transfer = 0;
  chemistry->use_volumetric_heating_rate = 0;
  chemistry->use_specific_heating_rate = 0;

  /* Initialize the chemistry object. */
  if (initialize_chemistry_data(&cooling->units) == 0) {
    error("Error in initialize_chemistry_data.");
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");
  float threshold = cooling->density_self_shielding;

  threshold /= phys_const->const_proton_mass;
  threshold /= pow(us->UnitLength_in_cgs, 3);

  message("***************************************");
  message("initializing grackle cooling function");
  message("");
  cooling_print_backend(cooling);
  message("Density Self Shielding = %g atom/cm3", threshold);

  message("");
  message("***************************************");
#endif
}

#endif /* SWIFT_COOLING_GRACKLE_H */
