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
#include <math.h>
#include <grackle.h>

/* Local includes. */
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1


/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_init_part(
    const struct part* restrict p, struct xpart* restrict xp) {

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

  /* set current time */
  float scale_factor;
  if (cooling->redshift == -1)
    error("TODO time dependant redshift");
  else
    scale_factor = 1. / (1. + cooling->redshift);

  /* Get current internal energy (dt=0) */
  const double energy_before = hydro_get_internal_energy(p);
  
  /* Get current density */
  const float rho = hydro_get_density(p);
  
  /* 0.02041 (= 1 Zsun in Grackle v2.0, but = 1.5761 Zsun in
     Grackle v2.1) */
  const double Z = 0.02041;

  /* create grackle struct */
  /* velocities */
  gr_float x_velocity[GRACKLE_NPART] = {0.0};
  gr_float y_velocity[GRACKLE_NPART] = {0.0};
  gr_float z_velocity[GRACKLE_NPART] = {0.0};

  /* particle data */
  gr_float density[GRACKLE_NPART] = {rho};
  gr_float metal_density[GRACKLE_NPART] = {Z * density[0]};
  gr_float energy[GRACKLE_NPART] = {energy_before};

  /* dimensions */
  int grid_dimension[3] = {GRACKLE_NPART, 0, 0};
  int grid_start[3] = {0, 0, 0};
  int grid_end[3] = {0, 0, 0};

  /* solve chemistry with table */
  if (solve_chemistry_table(&cooling->units, scale_factor, dt, grid_rank, grid_dimension,
                            grid_start, grid_end, density, energy, x_velocity,
                            y_velocity, z_velocity, metal_density) == 0) {
    error("Error in solve_chemistry.");
  }

  return (energy[0] - energy_before) / dt;
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
  xp->cooling_data.radiated_energy += - du_dt * dt * hydro_get_mass(p);

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
static INLINE void cooling_init_backend(
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

  /* We assume here all physical quantities to
     be in proper coordinate (not comobile)  */
  cooling->comoving_coordinates = 0;

  /* then units */
  cooling->units.density_units = us->UnitMass_in_cgs / pow(us->UnitLength_in_cgs, 3);
  cooling->units.length_units = us->UnitLength_in_cgs;
  cooling->units.time_units = us->UnitTime_in_cgs;
  cooling->units.velocity_units =
    cooling->units.a_units * cooling->units.length_units / cooling->units.time_units;
 
  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters() == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry.
  grackle_data.use_grackle = 1;
  grackle_data.with_radiative_cooling = 1;
  /* molecular network with H, He, D
   From Cloudy table */
  grackle_data.primordial_chemistry = 0;
  grackle_data.metal_cooling = 1;  // metal cooling on
  grackle_data.UVbackground = cooling->uv_background;
  grackle_data.grackle_data_file = cooling->cloudy_table;

  /* Initialize the chemistry object.
     a_value is not the true initial a
     This should get set before any computation */
  double a_value = 1.;

  if (initialize_chemistry_data(&cooling->units, a_value) == 0) {
    error("Error in initialize_chemistry_data.");
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");
  float threshold = cooling->GrackleHSShieldingDensityThreshold;

  threshold /= phys_const->const_proton_mass;
  threshold /= pow(us->UnitLength_in_cgs, 3);

  message("***************************************");
  message("initializing grackle cooling function");
  message("");
  cooling_print_backend(cooling);
  message("Density Self Shielding = %g atom/cm3", threshold);


  //grackle_print_data();
  message("");
  message("***************************************");
#endif
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle'.");
  message("CloudyTable             = %s",
          cooling->cloudy_table);
  message("UVbackground            = %d", cooling->uv_background);
  message("Redshift                = %g", cooling->redshift);
  message("Density Self Shielding  = %g",
          cooling->density_self_shielding);
  message("Units:");
  message("\tComoving     = %g", cooling->units.comoving_coordinates)
  message("\tLength  = %g", cooling->units.length_units);
  message("\tDensity = %g", cooling->units.density_units);
  message("\tTime         = %g", cooling->units.time_units);
  message("\tScale Factor = %g", cooling->units.a_units);  
}


/**
 * @brief print data in grackle struct
 *
 * Should only be used for debugging
 */
void grackle_print_data() {
  message("Grackle Data:");
  message("\t Data file: %s", grackle_data.grackle_data_file);
  message("\t With grackle: %i", grackle_data.use_grackle);
  message("\t With radiative cooling: %i", grackle_data.with_radiative_cooling);
  message("\t With UV background: %i", grackle_data.UVbackground);
  message("\t With primordial chemistry: %i",
          grackle_data.primordial_chemistry);
  message("\t Number temperature bins: %i",
          grackle_data.NumberOfTemperatureBins);
  message("\t T = (%g, ..., %g)", grackle_data.TemperatureStart,
          grackle_data.TemperatureEnd);

  message("Primordial Cloudy");
  cloudy_print_data(grackle_data.cloudy_primordial, 1);
  if (grackle_data.metal_cooling) {
    message("Metal Cooling");
    cloudy_print_data(grackle_data.cloudy_metal, 0);
  }

  message("\t Gamma: %g", grackle_data.Gamma);

  /* UVB */
  if (grackle_data.UVbackground && grackle_data.primordial_chemistry != 0) {
    struct UVBtable uvb = grackle_data.UVbackground_table;
    long long N = uvb.Nz;
    message("\t UV Background");
    message("\t\t Redshift from %g to %g with %lli steps", uvb.zmin, uvb.zmax,
            N);
    message("\t\t z = (%g, ..., %g)", uvb.z[0], uvb.z[N - 1]);
  }
}

/**
 * @brief print data in cloudy struct
 *
 * Should only be used for debugging
 */
void cloudy_print_data(const cloudy_data c, const int print_mmw) {
  long long N = c.data_size;
  message("\t Data size: %lli", N);
  message("\t Grid rank: %lli", c.grid_rank);

  char msg[200] = "\t Dimension: (";
  for (long long i = 0; i < c.grid_rank; i++) {
    char tmp[200] = "%lli%s";
    if (i == c.grid_rank - 1)
      sprintf(tmp, tmp, c.grid_dimension[i], ")");
    else
      sprintf(tmp, tmp, c.grid_dimension[i], ", ");

    strcat(msg, tmp);
  }
  message("%s", msg);

  if (c.heating_data)
    message("\t Heating: (%g, ..., %g)", c.heating_data[0],
            c.heating_data[N - 1]);
  if (c.cooling_data)
    message("\t Cooling: (%g, ..., %g)", c.cooling_data[0],
            c.cooling_data[N - 1]);
  if (c.mmw_data && print_mmw)
    message("\t Mean molecular weigth: (%g, ..., %g)", c.mmw_data[0],
            c.mmw_data[N - 1]);
}

#endif /* SWIFT_COOLING_GRACKLE_H */
