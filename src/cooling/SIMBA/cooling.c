/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
/**
 * @file src/cooling/SIMBA/cooling.c
 * @brief Cooling using the GRACKLE 3.x library.
 */

#include "../config.h"

/* Include header */
#include "cooling.h"

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <math.h>

/* The grackle library itself */
#include <grackle.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling_io.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
void cooling_update(const struct cosmology* cosmo,
                    struct cooling_function_data* cooling, struct space* s) {
  /* set current time */
  if (cooling->redshift == -1)
    cooling->units.a_value = cosmo->a;
  else
    cooling->units.a_value = 1. / (1. + cooling->redshift);
}

/**
 * @brief Print the chemical network
 *
 * @param xp The #xpart to print
 */
void cooling_print_fractions(const struct xpart* restrict xp) {

  const struct cooling_xpart_data tmp = xp->cooling_data;
#if COOLING_GRACKLE_MODE > 0
  message("HI %g, HII %g, HeI %g, HeII %g, HeIII %g, e %g", tmp.HI_frac,
          tmp.HII_frac, tmp.HeI_frac, tmp.HeII_frac, tmp.HeIII_frac,
          tmp.e_frac);
#endif

#if COOLING_GRACKLE_MODE > 1
  message("HM %g, H2I %g, H2II %g", tmp.HM_frac, tmp.H2I_frac, tmp.H2II_frac);
#endif

#if COOLING_GRACKLE_MODE > 2
  message("DI %g, DII %g, HDI %g", tmp.DI_frac, tmp.DII_frac, tmp.HDI_frac);
#endif
  message("Metal: %g", tmp.metal_frac);
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_first_init_part(const struct phys_const* restrict phys_const,
                             const struct unit_system* restrict us,
                             const struct hydro_props* hydro_props,
                             const struct cosmology* restrict cosmo,
                             const struct cooling_function_data* cooling,
                             const struct part* restrict p,
                             struct xpart* restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.time_last_event = -cooling->thermal_time;

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;

  /* primordial chemistry >= 1: Start with everything neutral (as in dark ages)
   */
  xp->cooling_data.HI_frac = grackle_data->HydrogenFractionByMass;
  xp->cooling_data.HII_frac = zero;
  xp->cooling_data.HeI_frac = 1. - grackle_data->HydrogenFractionByMass;
  xp->cooling_data.HeII_frac = zero;
  xp->cooling_data.HeIII_frac = zero;
  xp->cooling_data.e_frac = xp->cooling_data.HII_frac +
                            0.25 * xp->cooling_data.HeII_frac +
                            0.5 * xp->cooling_data.HeIII_frac;
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  xp->cooling_data.HM_frac = zero;
  xp->cooling_data.H2I_frac = zero;
  xp->cooling_data.H2II_frac = zero;
#endif  // MODE >= 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  xp->cooling_data.DI_frac = grackle_data->DeuteriumToHydrogenRatio *
                             grackle_data->HydrogenFractionByMass;
  xp->cooling_data.DII_frac = zero;
  xp->cooling_data.HDI_frac = zero;
#endif  // MODE >= 3
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
float cooling_get_radiated_energy(const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print_backend(const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle'.");
  message("Using Grackle = %i", cooling->chemistry.use_grackle);
  message("Chemical network = %i", cooling->chemistry.primordial_chemistry);
  message("CloudyTable = %s", cooling->cloudy_table);
  message("Redshift = %g", cooling->redshift);
  message("UV background = %d", cooling->with_uv_background);
  message("Metal cooling = %i", cooling->chemistry.metal_cooling);
  message("Self Shielding = %i", cooling->self_shielding_method);
  if (cooling->self_shielding_method == -1) {
    message("Self Shelding density = %g", cooling->self_shielding_threshold);
  }
  message("Thermal time = %g", cooling->thermal_time);
  message("Specific Heating Rates = %i",
          cooling->provide_specific_heating_rates);
  message("Volumetric Heating Rates = %i",
          cooling->provide_volumetric_heating_rates);
  message("Units:");
  message("\tComoving = %i", cooling->units.comoving_coordinates);
  message("\tLength = %g", cooling->units.length_units);
  message("\tDensity = %g", cooling->units.density_units);
  message("\tTime = %g", cooling->units.time_units);
  message("\tScale Factor = %g (units: %g)", cooling->units.a_value,
          cooling->units.a_units);
}

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 0
void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho,
                              gr_float species_densities[12]) {
  /* HI */
  species_densities[0] = xp->cooling_data.HI_frac * rho;
  data->HI_density = &species_densities[0];
  /* HII */
  species_densities[1] = xp->cooling_data.HII_frac * rho;
  data->HII_density = &species_densities[1];

  /* HeI */
  species_densities[2] = xp->cooling_data.HeI_frac * rho;
  data->HeI_density = &species_densities[2];

  /* HeII */
  species_densities[3] = xp->cooling_data.HeII_frac * rho;
  data->HeII_density = &species_densities[3];

  /* HeIII */
  species_densities[4] = xp->cooling_data.HeIII_frac * rho;
  data->HeIII_density = &species_densities[4];

  /* HeII */
  species_densities[5] = xp->cooling_data.e_frac * rho;
  data->e_density = &species_densities[5];
}
#else
void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho,
                              gr_float species_densities[12]) {
  data->HI_density = NULL;
  data->HII_density = NULL;
  data->HeI_density = NULL;
  data->HeII_density = NULL;
  data->HeIII_density = NULL;
  data->e_density = NULL;
}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 1
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho,
                              gr_float species_densities[12]) {
  /* HM */
  species_densities[6] = xp->cooling_data.HM_frac * rho;
  data->HM_density = &species_densities[6];

  /* H2I */
  species_densities[7] = xp->cooling_data.H2I_frac * rho;
  data->H2I_density = &species_densities[7];

  /* H2II */
  species_densities[8] = xp->cooling_data.H2II_frac * rho;
  data->H2II_density = &species_densities[8];
}
#else
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho,
                              gr_float species_densities[12]) {
  data->HM_density = NULL;
  data->H2I_density = NULL;
  data->H2II_density = NULL;
}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 2
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho,
                              gr_float species_densities[12]) {
  /* DI */
  species_densities[9] = xp->cooling_data.DI_frac * rho;
  data->DI_density = &species_densities[9];

  /* DII */
  species_densities[10] = xp->cooling_data.DII_frac * rho;
  data->DII_density = &species_densities[10];

  /* HDI */
  species_densities[11] = xp->cooling_data.HDI_frac * rho;
  data->HDI_density = &species_densities[11];
}
#else
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho,
                              gr_float species_densities[12]) {
  data->DI_density = NULL;
  data->DII_density = NULL;
  data->HDI_density = NULL;
}
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
#if COOLING_GRACKLE_MODE > 0
void cooling_copy_from_grackle1(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {

  /* HI */
  xp->cooling_data.HI_frac = *data->HI_density / rho;

  /* HII */
  xp->cooling_data.HII_frac = *data->HII_density / rho;

  /* HeI */
  xp->cooling_data.HeI_frac = *data->HeI_density / rho;

  /* HeII */
  xp->cooling_data.HeII_frac = *data->HeII_density / rho;

  /* HeIII */
  xp->cooling_data.HeIII_frac = *data->HeIII_density / rho;

  /* e */
  xp->cooling_data.e_frac = *data->e_density / rho;
}
#else
void cooling_copy_from_grackle1(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {}
#endif

/**
 * @brief copy the grackle data to a #xpart.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
#if COOLING_GRACKLE_MODE > 1
void cooling_copy_from_grackle2(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {
  /* HM */
  xp->cooling_data.HM_frac = *data->HM_density / rho;
  /* H2I */
  xp->cooling_data.H2I_frac = *data->H2I_density / rho;
  /* H2II */
  xp->cooling_data.H2II_frac = *data->H2II_density / rho;
}
#else
void cooling_copy_from_grackle2(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {}
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
#if COOLING_GRACKLE_MODE > 2
void cooling_copy_from_grackle3(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {

  /* DI */
  xp->cooling_data.DI_frac = *data->DI_density / rho;

  /* DII */
  xp->cooling_data.DII_frac = *data->DII_density / rho;

  /* HDI */
  xp->cooling_data.HDI_frac = *data->HDI_density / rho;
}
#else
void cooling_copy_from_grackle3(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
void cooling_copy_to_grackle(grackle_field_data* data,
                             const struct cosmology* restrict cosmo,
                             const struct part* p, struct xpart* xp,
                             gr_float species_densities[12]) {

  /* set values */
  /* grid */
  int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  int grid_start[GRACKLE_RANK] = {0, 0, 0};
  int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};

  data->grid_rank = GRACKLE_RANK;
  data->grid_dx = 0.f;
  data->grid_dimension = grid_dimension;
  data->grid_start = grid_start;
  data->grid_end = grid_end;

  /* get general particle data in physical cgs units */
  gr_float rho = hydro_get_physical_density(p, cosmo);
  gr_float energy = hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* load particle data into grackle structure */
  data->density = &rho;
  data->internal_energy = &energy;
  data->specific_heating_rate = &hydro_du_dt;

  /* velocity (maybe not needed?) */
  double v_tmp[3] = {p->v[0], p->v[1], p->v[2]};
  data->x_velocity = &v_tmp[0];
  data->y_velocity = &v_tmp[1];
  data->z_velocity = &v_tmp[2];

  cooling_copy_to_grackle1(data, p, xp, rho, species_densities);
  cooling_copy_to_grackle2(data, p, xp, rho, species_densities);
  cooling_copy_to_grackle3(data, p, xp, rho, species_densities);

  data->volumetric_heating_rate = NULL;
  data->RT_heating_rate = NULL;
  data->RT_HI_ionization_rate = NULL;
  data->RT_HeI_ionization_rate = NULL;
  data->RT_HeII_ionization_rate = NULL;
  data->RT_H2_dissociation_rate = NULL;

  gr_float* metal_density = (gr_float*)malloc(sizeof(gr_float));
  *metal_density = chemistry_get_total_metal_mass_fraction_for_cooling(p) * rho;
  data->metal_density = metal_density;
}

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
void cooling_copy_from_grackle(grackle_field_data* data, const struct part* p,
                               struct xpart* xp, gr_float rho) {

  cooling_copy_from_grackle1(data, p, xp, rho);
  cooling_copy_from_grackle2(data, p, xp, rho);
  cooling_copy_from_grackle3(data, p, xp, rho);

  free(data->metal_density);
}

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param mode 0=energy, 1=cooling time, 2=temperature, 3=pressure, 4=gamma
 *
 * @return desired quantity based on mode
 */
gr_float cooling_grackle_driver(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp, double dt,
    int mode) {

  /* set current units for conversion to physical quantities */
  code_units units = cooling->units;

  /* initialize data to send to grackle */
  grackle_field_data data;
  gr_float species_densities[12];

  /* copy species_densities from particle to grackle data */
  cooling_copy_to_grackle(&data, cosmo, p, xp, species_densities);

  /* Run Grackle in desired mode */
  gr_float return_value = 0.f;
  switch (mode) {
    case 0:
      /* solve chemistry, advance thermal energy by dt */
      if (solve_chemistry(&units, &data, dt) == 0) {
        error("Error in Grackle solve_chemistry.");
      }
      /* copy from grackle data to particle */
      gr_float rho = hydro_get_physical_density(p, cosmo);
      cooling_copy_from_grackle(&data, p, xp, rho);
      return_value = data.internal_energy[0];
      break;
    case 1:
      /* compute cooling time */
      if (calculate_cooling_time(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_cooling_time.");
      }
      break;
    case 2:
      /* compute temperature */
      if (calculate_temperature(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_temperature.");
      }
      break;
    case 3:
      /* compute pressure */
      if (calculate_pressure(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_pressure.");
      }
      break;
    case 4:
      /* compute gamma */
      if (calculate_gamma(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_gamma.");
      }
      break;
  }

  return return_value;
}

/**
 * @brief Compute the cooling time
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 *
 * @return cooling time
 */
gr_float cooling_time(const struct phys_const* restrict phys_const,
                      const struct unit_system* restrict us,
                      const struct hydro_props* hydro_properties,
                      const struct cosmology* restrict cosmo,
                      const struct cooling_function_data* restrict cooling,
                      const struct part* restrict p,
                      struct xpart* restrict xp) {

  gr_float cooling_time = cooling_grackle_driver(
      phys_const, us, cosmo, hydro_properties, cooling, p, xp, 0., 1);
  return cooling_time;
}

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_properties,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  struct xpart xp_temp = *xp;  // gets rid of const in declaration
  float temperature = cooling_grackle_driver(
      phys_const, us, cosmo, hydro_properties, cooling, p, &xp_temp, 0., 2);
  return temperature;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
void cooling_cool_part(const struct phys_const* restrict phys_const,
                       const struct unit_system* restrict us,
                       const struct cosmology* restrict cosmo,
                       const struct hydro_props* hydro_props,
                       const struct entropy_floor_properties* floor_props,
                       const struct cooling_function_data* restrict cooling,
                       struct part* restrict p, struct xpart* restrict xp,
                       const double dt, const double dt_therm,
                       const double time) {

  /* Nothing to do here? */
  // message("GRACKLE: z=%g  dt=%g  dt_therm=%g", cosmo->z, dt, dt_therm);
  if (dt == 0.) return;

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Calculate energy after dt */
  gr_float u_new;

  /* Is the cooling turned off */
  if (time - xp->cooling_data.time_last_event < cooling->thermal_time) {
    u_new = u_old;
  } else {
    u_new = cooling_grackle_driver(phys_const, us, cosmo, hydro_props, cooling,
                                   p, xp, dt_therm, 0);
  }

  /* We now need to check that we are not going to go below any of the limits */
  u_new = max(u_new, hydro_props->minimal_internal_energy);

  /* Update the internal energy time derivative */
  float cool_du_dt = (u_new - u_old) / dt_therm;
  hydro_set_physical_internal_energy_dt(p, cosmo, cool_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cool_du_dt * dt_therm;
}

/**
 * @brief Compute the y-Compton contribution of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
double Cooling_get_ycompton(const struct phys_const* phys_const,
                            const struct hydro_props* hydro_props,
                            const struct unit_system* us,
                            const struct cosmology* cosmo,
                            const struct cooling_function_data* cooling,
                            const struct part* p, const struct xpart* xp) {

  return 0.;
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step,
 * since Grackle sub-cycles the cooling as needed.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The #cosmology.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 */
float cooling_timestep(const struct cooling_function_data* restrict cooling,
                       const struct phys_const* restrict phys_const,
                       const struct cosmology* restrict cosmo,
                       const struct unit_system* restrict us,
                       const struct hydro_props* hydro_props,
                       const struct part* restrict p,
                       const struct xpart* restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part* p, struct xpart* xp, double n) {

  xp->cooling_data.radiated_energy /= n;
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param phys_const The #phys_const.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_units(const struct unit_system* us,
                        const struct phys_const* phys_const,
                        struct cooling_function_data* cooling) {

  /* These are conversions from code units to cgs. */

  /* first cosmo */
  cooling->units.a_units = 1.0;  // units for the expansion factor
  cooling->units.a_value = 1.0;

  /* We assume here all physical quantities to
     be in proper coordinate (not comoving)  */
  cooling->units.comoving_coordinates = 0;

  /* then units */
  cooling->units.density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->units.length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  cooling->units.time_units = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->units.velocity_units =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  /* Self shielding */
  if (cooling->self_shielding_method == -1) {
    cooling->self_shielding_threshold *=
        phys_const->const_proton_mass *
        pow(units_cgs_conversion_factor(us, UNIT_CONV_LENGTH), 3.);
  }
}

/**
 * @brief Initialises Grackle.
 *
 * @param cooling The cooling properties to initialize
 */
void cooling_init_grackle(struct cooling_function_data* cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

  chemistry_data* chemistry = &cooling->chemistry;

  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters(chemistry) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry & cooling

  // Flag to activate the grackle machinery:
  chemistry->use_grackle = 1;  // grackle on (duh)
  // Flag to include radiative cooling and actually update the thermal energy
  // during the chemistry solver. If off, the chemistry species will still be
  // updated. The most common reason to set this to off is to iterate the
  // chemistry network to an equilibrium state. Default: 1.
  chemistry->with_radiative_cooling = 1;  // cooling on
  // Flag to control which primordial chemistry network is used (set by Config
  // file)
  chemistry->primordial_chemistry = COOLING_GRACKLE_MODE;
  // Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat
  // transfer follow Omukai (2000). This assumes that the dust to gas ratio
  // scales with the metallicity. Default: 0.
  chemistry->h2_on_dust = 0;  // dust cooling/chemistry on
  // Flag to enable metal cooling using the Cloudy tables. If enabled, the
  // cooling table to be used must be specified with the grackle_data_file
  // parameter. Default: 0.
  chemistry->metal_cooling = 1;  // metal cooling on
  // Flag to enable an effective CMB temperature floor. This is implemented by
  // subtracting the value of the cooling rate at TCMB from the total cooling
  // rate. Default: 1.
  chemistry->cmb_temperature_floor = 1;
  // Flag to enable a UV background. If enabled, the cooling table to be used
  // must be specified with the grackle_data_file parameter. Default: 0.
  chemistry->UVbackground = 1;  // UV background on
  // Path to the data file containing the metal cooling and UV background
  // tables:
  chemistry->grackle_data_file = cooling->cloudy_table;  // data file
  // The ratio of specific heats for an ideal gas. A direct calculation for the
  // molecular component is used if primordial_chemistry > 1. Default: 5/3.
  chemistry->Gamma = hydro_gamma;  // our eos set in Config.sh
  // Flag to control which three-body H2 formation rate is used.
  chemistry->three_body_rate = 0;
  // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel
  // (2004). Default: 0.
  chemistry->cie_cooling = 0;
  // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004).
  // Default: 0
  chemistry->h2_optical_depth_approximation = 0;
  // Flag to enable a spatially uniform heating term approximating
  // photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
  chemistry->photoelectric_heating =
      0;  // photo-electric on [but not adjusted to local background, beware!]
  chemistry->photoelectric_heating_rate = 8.5e-26;
  // Flag to enable Compton heating from an X-ray background following Madau &
  // Efstathiou (1999). Default: 0.
  chemistry->Compton_xray_heating = 0;
  // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
  // in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
  chemistry->LWbackground_intensity = 0;
  // Flag to enable suppression of Lyman-Werner flux due to Lyman-series
  // absorption
  //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000).
  //    Default: 0.
  chemistry->LWbackground_sawtooth_suppression = 0;
  // volumetric heating rates is being provided in the volumetric_heating_rate
  // field of grackle_field_data
  chemistry->use_volumetric_heating_rate = 0;
  // specific heating rates is being provided in the specific_heating_rate field
  // of grackle_field_data
  chemistry->use_specific_heating_rate = 1;
  // arrays of ionization and heating rates from radiative transfer solutions
  // are being provided
  chemistry->use_radiative_transfer = 0;
  // must be enabled to couple the passed radiative transfer fields to the
  // chemistry solver
  chemistry->radiative_transfer_coupled_rate_solver = 0;
  // enable intermediate stepping in applying radiative transfer fields to
  // chemistry solver.
  chemistry->radiative_transfer_intermediate_step = 0;
  // only use hydrogen ionization and heating rates from the radiative transfer
  // solutions.
  chemistry->radiative_transfer_hydrogen_only = 0;
  // Use Rahmati+13 self-shielding; 0=none, 1=HI only, 2=HI+HeI, 3=HI+HeI but
  // set HeII rates to 0
  chemistry->self_shielding_method = 0;
  chemistry->self_shielding_method = cooling->self_shielding_method;

  /* Initialize the chemistry object. */
  if (initialize_chemistry_data(&cooling->units) == 0) {
    error("Error in initialize_chemistry_data.");
  }
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_backend(struct swift_params* parameter_file,
                          const struct unit_system* us,
                          const struct phys_const* phys_const,
                          const struct hydro_props* hydro_props,
                          struct cooling_function_data* cooling) {

  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");

  /* read parameters */
  cooling_read_parameters(parameter_file, cooling, phys_const);

  /* Set up the units system. */
  cooling_init_units(us, phys_const, cooling);

  /* Set up grackle */
  cooling_init_grackle(cooling);
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data* cooling) {
  //_free_chemistry_data(&cooling->chemistry, &grackle_rates);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * Nothing to do beyond writing the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data* cooling,
                         FILE* stream) {
  restart_write_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
                       stream, "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Nothing to do beyond reading the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data* cooling, FILE* stream,
                            const struct cosmology* cosmo) {
  restart_read_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  /* Set up grackle */
  cooling_init_grackle(cooling);
}
