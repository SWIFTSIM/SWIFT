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
 * @file src/cooling/grackle/cooling.h
 * @brief Cooling using the GRACKLE 3.0 library.
 */

/* Some standard headers. */
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

/* prototype */
static gr_float cooling_time(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp);
static gr_float cooling_new_energy(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp, double dt);

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
INLINE static void cooling_update(const struct cosmology* cosmo,
                                  struct cooling_function_data* cooling,
                                  struct space* s) {
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
__attribute__((always_inline)) INLINE static void cooling_print_fractions(
    const struct xpart* restrict xp) {

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
 * @brief Check if the difference of a given field is lower than limit
 *
 * @param xp First #xpart
 * @param old Second #xpart
 * @param field The field to check
 * @param limit Difference limit
 *
 * @return 0 if diff > limit
 */
#define cooling_check_field(xp, old, field, limit)                \
  ({                                                              \
    float tmp = xp->cooling_data.field - old->cooling_data.field; \
    tmp = fabs(tmp) / xp->cooling_data.field;                     \
    if (tmp > limit) return 0;                                    \
  })

/**
 * @brief Check if difference between two particles is lower than a given value
 *
 * @param xp One of the #xpart
 * @param old The other #xpart
 * @param limit The difference limit
 */
__attribute__((always_inline)) INLINE static int cooling_converged(
    const struct xpart* restrict xp, const struct xpart* restrict old,
    const float limit) {

#if COOLING_GRACKLE_MODE > 0
  cooling_check_field(xp, old, HI_frac, limit);
  cooling_check_field(xp, old, HII_frac, limit);
  cooling_check_field(xp, old, HeI_frac, limit);
  cooling_check_field(xp, old, HeII_frac, limit);
  cooling_check_field(xp, old, HeIII_frac, limit);
  cooling_check_field(xp, old, e_frac, limit);
#endif
#if COOLING_GRACKLE_MODE > 1
  cooling_check_field(xp, old, HM_frac, limit);
  cooling_check_field(xp, old, H2I_frac, limit);
  cooling_check_field(xp, old, H2II_frac, limit);
#endif

#if COOLING_GRACKLE_MODE > 2
  cooling_check_field(xp, old, DI_frac, limit);
  cooling_check_field(xp, old, DII_frac, limit);
  cooling_check_field(xp, old, HDI_frac, limit);
#endif

  return 1;
}

/**
 * @brief Compute equilibrium fraction
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_compute_equilibrium(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp) {

  /* get temporary data */
  struct part p_tmp = *p;
  struct cooling_function_data cooling_tmp = *cooling;
  cooling_tmp.chemistry.with_radiative_cooling = 0;
  /* need density for computation, therefore quick estimate */
  p_tmp.rho = 0.2387 * p_tmp.mass / pow(p_tmp.h, 3);

  /* compute time step */
  const double alpha = 0.01;
  double dt =
      fabs(cooling_time(phys_const, us, cosmo, &cooling_tmp, &p_tmp, xp));
  cooling_new_energy(phys_const, us, cosmo, &cooling_tmp, &p_tmp, xp, dt);
  dt = alpha *
       fabs(cooling_time(phys_const, us, cosmo, &cooling_tmp, &p_tmp, xp));

  /* init simple variables */
  int step = 0;
  const int max_step = cooling_tmp.max_step;
  const float conv_limit = cooling_tmp.convergence_limit;
  struct xpart old;

  do {
    /* update variables */
    step += 1;
    old = *xp;

    /* update chemistry */
    cooling_new_energy(phys_const, us, cosmo, &cooling_tmp, &p_tmp, xp, dt);
  } while (step < max_step && !cooling_converged(xp, &old, conv_limit));

  if (step == max_step)
    error(
        "A particle element fraction failed to converge."
        "You can change 'GrackleCooling:MaxSteps' or "
        "'GrackleCooling:ConvergenceLimit' to avoid this problem");
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
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* cooling, const struct part* restrict p,
    struct xpart* restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;

  /* primordial chemistry >= 1 */
  xp->cooling_data.HI_frac = zero;
  xp->cooling_data.HII_frac = grackle_data->HydrogenFractionByMass;
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

#if COOLING_GRACKLE_MODE > 0
  cooling_compute_equilibrium(phys_const, us, cosmo, cooling, p, xp);
#endif
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
  message("Using Grackle    = %i", cooling->chemistry.use_grackle);
  message("Chemical network = %i", cooling->chemistry.primordial_chemistry);
  message("CloudyTable      = %s", cooling->cloudy_table);
  message("Redshift         = %g", cooling->redshift);
  message("UV background    = %d", cooling->with_uv_background);
  message("Metal cooling    = %i", cooling->chemistry.metal_cooling);
  message("Self Shielding   = %i", cooling->self_shielding_method);
  message("Specific Heating Rates   = %i",
          cooling->provide_specific_heating_rates);
  message("Volumetric Heating Rates = %i",
          cooling->provide_volumetric_heating_rates);
  message("Units:");
  message("\tComoving     = %i", cooling->units.comoving_coordinates);
  message("\tLength       = %g", cooling->units.length_units);
  message("\tDensity      = %g", cooling->units.density_units);
  message("\tTime         = %g", cooling->units.time_units);
  message("\tScale Factor = %g (units: %g)", cooling->units.a_value,
          cooling->units.a_units);
}

/**
 * @brief copy a single field from the grackle data to a #xpart
 *
 * @param data The #grackle_field_data
 * @param xp The #xpart
 * @param rho Particle density
 * @param field The field to copy
 */
#define cooling_copy_field_from_grackle(data, xp, rho, field) \
  xp->cooling_data.field##_frac = *data.field##_density / rho;

/**
 * @brief copy a single field from a #xpart to the grackle data
 *
 * @param data The #grackle_field_data
 * @param xp The #xpart
 * @param rho Particle density
 * @param field The field to copy
 */
#define cooling_copy_field_to_grackle(data, xp, rho, field)       \
  gr_float grackle_##field = xp->cooling_data.field##_frac * rho; \
  data.field##_density = &grackle_##field;

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 0
#define cooling_copy_to_grackle1(data, p, xp, rho)     \
  cooling_copy_field_to_grackle(data, xp, rho, HI);    \
  cooling_copy_field_to_grackle(data, xp, rho, HII);   \
  cooling_copy_field_to_grackle(data, xp, rho, HeI);   \
  cooling_copy_field_to_grackle(data, xp, rho, HeII);  \
  cooling_copy_field_to_grackle(data, xp, rho, HeIII); \
  cooling_copy_field_to_grackle(data, xp, rho, e);
#else
#define cooling_copy_to_grackle1(data, p, xp, rho) \
  data.HI_density = NULL;                          \
  data.HII_density = NULL;                         \
  data.HeI_density = NULL;                         \
  data.HeII_density = NULL;                        \
  data.HeIII_density = NULL;                       \
  data.e_density = NULL;
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 1
#define cooling_copy_to_grackle2(data, p, xp, rho)   \
  cooling_copy_field_to_grackle(data, xp, rho, HM);  \
  cooling_copy_field_to_grackle(data, xp, rho, H2I); \
  cooling_copy_field_to_grackle(data, xp, rho, H2II);
#else
#define cooling_copy_to_grackle2(data, p, xp, rho) \
  data.HM_density = NULL;                          \
  data.H2I_density = NULL;                         \
  data.H2II_density = NULL;
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 2
#define cooling_copy_to_grackle3(data, p, xp, rho)   \
  cooling_copy_field_to_grackle(data, xp, rho, DI);  \
  cooling_copy_field_to_grackle(data, xp, rho, DII); \
  cooling_copy_field_to_grackle(data, xp, rho, HDI);
#else
#define cooling_copy_to_grackle3(data, p, xp, rho) \
  data.DI_density = NULL;                          \
  data.DII_density = NULL;                         \
  data.HDI_density = NULL;
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 0
#define cooling_copy_from_grackle1(data, p, xp, rho)     \
  cooling_copy_field_from_grackle(data, xp, rho, HI);    \
  cooling_copy_field_from_grackle(data, xp, rho, HII);   \
  cooling_copy_field_from_grackle(data, xp, rho, HeI);   \
  cooling_copy_field_from_grackle(data, xp, rho, HeII);  \
  cooling_copy_field_from_grackle(data, xp, rho, HeIII); \
  cooling_copy_field_from_grackle(data, xp, rho, e);
#else
#define cooling_copy_from_grackle1(data, p, xp, rho)
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 1
#define cooling_copy_from_grackle2(data, p, xp, rho)   \
  cooling_copy_field_from_grackle(data, xp, rho, HM);  \
  cooling_copy_field_from_grackle(data, xp, rho, H2I); \
  cooling_copy_field_from_grackle(data, xp, rho, H2II);
#else
#define cooling_copy_from_grackle2(data, p, xp, rho)
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 2
#define cooling_copy_from_grackle3(data, p, xp, rho)   \
  cooling_copy_field_from_grackle(data, xp, rho, DI);  \
  cooling_copy_field_from_grackle(data, xp, rho, DII); \
  cooling_copy_field_from_grackle(data, xp, rho, HDI);
#else
#define cooling_copy_from_grackle3(data, p, xp, rho)
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#define cooling_copy_to_grackle(data, p, xp, rho)                      \
  cooling_copy_to_grackle1(data, p, xp, rho);                          \
  cooling_copy_to_grackle2(data, p, xp, rho);                          \
  cooling_copy_to_grackle3(data, p, xp, rho);                          \
  data.volumetric_heating_rate = NULL;                                 \
  data.specific_heating_rate = NULL;                                   \
  data.RT_heating_rate = NULL;                                         \
  data.RT_HI_ionization_rate = NULL;                                   \
  data.RT_HeI_ionization_rate = NULL;                                  \
  data.RT_HeII_ionization_rate = NULL;                                 \
  data.RT_H2_dissociation_rate = NULL;                                 \
  gr_float metal_density = chemistry_metal_mass_fraction(p, xp) * rho; \
  data.metal_density = &metal_density;

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#define cooling_copy_from_grackle(data, p, xp, rho) \
  cooling_copy_from_grackle1(data, p, xp, rho);     \
  cooling_copy_from_grackle2(data, p, xp, rho);     \
  cooling_copy_from_grackle3(data, p, xp, rho);

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 *
 * @return du / dt
 */
__attribute__((always_inline)) INLINE static gr_float cooling_new_energy(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp, double dt) {

  /* set current time */
  code_units units = cooling->units;

  /* initialize data */
  grackle_field_data data;

  /* set values */
  /* grid */
  int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  int grid_start[GRACKLE_RANK] = {0, 0, 0};
  int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};

  data.grid_dx = 0.;
  data.grid_rank = GRACKLE_RANK;
  data.grid_dimension = grid_dimension;
  data.grid_start = grid_start;
  data.grid_end = grid_end;

  /* general particle data */
  gr_float density = hydro_get_physical_density(p, cosmo);
  const float energy_before = hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float energy = energy_before;

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  /* copy to grackle structure */
  cooling_copy_to_grackle(data, p, xp, density);

  /* solve chemistry */
  chemistry_data chemistry_grackle = cooling->chemistry;
  if (local_solve_chemistry(&chemistry_grackle, &grackle_rates, &units, &data,
                            dt) == 0) {
    error("Error in solve_chemistry.");
  }

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(data, p, xp, density);

  return energy;
}

/**
 * @brief Compute the cooling time
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 *
 * @return cooling time
 */
__attribute__((always_inline)) INLINE static gr_float cooling_time(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  /* set current time */
  code_units units = cooling->units;

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
  const gr_float energy_before =
      hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float density = hydro_get_physical_density(p, cosmo);
  gr_float energy = energy_before;

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  /* copy data from particle to grackle data */
  cooling_copy_to_grackle(data, p, xp, density);

  /* Compute cooling time */
  gr_float cooling_time;
  chemistry_data chemistry_grackle = cooling->chemistry;
  chemistry_data_storage chemistry_rates = grackle_rates;
  if (local_calculate_cooling_time(&chemistry_grackle, &chemistry_rates, &units,
                                   &data, &cooling_time) == 0) {
    error("Error in calculate_cooling_time.");
  }

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(data, p, xp, density);

  /* compute rate */
  return cooling_time;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param dt The time-step of this particle.
 * @param hydro_properties the hydro_props struct, used for
 * getting the minimal internal energy allowed in by SWIFT.
 * Read from yml file into engine struct.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, double dt,
    double dt_therm) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Current du_dt in physical coordinates (internal units) */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Calculate energy after dt */
  gr_float u_new =
      cooling_new_energy(phys_const, us, cosmo, cooling, p, xp, dt);

  float delta_u = u_new - u_old + hydro_du_dt * dt_therm;

  /* We now need to check that we are not going to go below any of the limits */

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_u. */
  if (u_old + 1.5 * delta_u < hydro_props->minimal_internal_energy) {
    delta_u = (hydro_props->minimal_internal_energy - u_old) / 1.5;
  }

  /* Second, check whether the energy used in the prediction could get negative.
   * We need to check for the 1/2 dt kick followed by a full time-step drift
   * that could potentially be for a time-step twice as big. We hence check
   * for 2.5 delta_u but this time against 0 energy not the minimum.
   * To avoid numerical rounding bringing us below 0., we add a tiny tolerance.
   */
  const float rounding_tolerance = 1.0e-4;

  if (u_old + 2.5 * delta_u < 0.) {
    delta_u = -u_old / (2.5 + rounding_tolerance);
  }

  /* Turn this into a rate of change (including cosmology term) */
  const float cooling_du_dt = delta_u / dt_therm;

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cooling_du_dt * dt;
}

static INLINE float cooling_get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  error("This function needs implementing!!");
  return 0.;
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props, const struct part* restrict p,
    const struct xpart* restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_units(
    const struct unit_system* us, struct cooling_function_data* cooling) {

  /* These are conversions from code units to cgs. */

  /* first cosmo */
  cooling->units.a_units = 1.0;  // units for the expansion factor
  cooling->units.a_value = 1.0;

  /* We assume here all physical quantities to
     be in proper coordinate (not comobile)  */
  cooling->units.comoving_coordinates = 0;

  /* then units */
  cooling->units.density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->units.length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  cooling->units.time_units = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->units.velocity_units =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);
}

/**
 * @brief Initialises Grackle.
 *
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_grackle(
    struct cooling_function_data* cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

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
  chemistry->primordial_chemistry = COOLING_GRACKLE_MODE;
  chemistry->metal_cooling = cooling->with_metal_cooling;
  chemistry->UVbackground = cooling->with_uv_background;
  chemistry->grackle_data_file = cooling->cloudy_table;

  /* radiative transfer */
  chemistry->use_radiative_transfer = cooling->provide_specific_heating_rates ||
                                      cooling->provide_volumetric_heating_rates;
  chemistry->use_volumetric_heating_rate =
      cooling->provide_volumetric_heating_rates;
  chemistry->use_specific_heating_rate =
      cooling->provide_specific_heating_rates;

  if (cooling->provide_specific_heating_rates &&
      cooling->provide_volumetric_heating_rates)
    message(
        "WARNING: You should specified either the specific or the volumetric "
        "heating rates, not both");

  /* self shielding */
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
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_backend(
    struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {

  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");

  /* read parameters */
  cooling_read_parameters(parameter_file, cooling);

  /* Set up the units system. */
  cooling_init_units(us, cooling);

  /* Set up grackle */
  cooling_init_grackle(cooling);
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
static INLINE void cooling_clean(struct cooling_function_data* cooling) {

  // MATTHIEU: To do: free stuff here
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * Nothing to do beyond writing the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
static INLINE void cooling_struct_dump(
    const struct cooling_function_data* cooling, FILE* stream) {
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
static INLINE void cooling_struct_restore(struct cooling_function_data* cooling,
                                          FILE* stream,
                                          const struct cosmology* cosmo) {
  restart_read_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");
}

#endif /* SWIFT_COOLING_GRACKLE_H */
