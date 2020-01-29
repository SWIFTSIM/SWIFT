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
/**
 * @file src/cooling/grackle/cooling.c
 * @brief Cooling using the GRACKLE 3.0 library.
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
int cooling_converged(const struct xpart* restrict xp,
                      const struct xpart* restrict old, const float limit) {

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
void cooling_compute_equilibrium(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_properties,
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
  double dt = fabs(cooling_time(phys_const, us, hydro_properties, cosmo,
                                &cooling_tmp, &p_tmp, xp));
  cooling_new_energy(phys_const, us, cosmo, hydro_properties, &cooling_tmp,
                     &p_tmp, xp, dt);
  dt = alpha * fabs(cooling_time(phys_const, us, hydro_properties, cosmo,
                                 &cooling_tmp, &p_tmp, xp));

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
    cooling_new_energy(phys_const, us, cosmo, hydro_properties, &cooling_tmp,
                       &p_tmp, xp, dt);
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
  cooling_compute_equilibrium(phys_const, us, hydro_props, cosmo, cooling, p,
                              xp);
#endif
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
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 0
void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho) {
  /* HI */
  xp->cooling_data.HI_frac *= rho;
  data->HI_density = &xp->cooling_data.HI_frac;
  /* HII */
  xp->cooling_data.HII_frac *= rho;
  data->HII_density = &xp->cooling_data.HII_frac;

  /* HeI */
  xp->cooling_data.HeI_frac *= rho;
  data->HeI_density = &xp->cooling_data.HeI_frac;

  /* HeII */
  xp->cooling_data.HeII_frac *= rho;
  data->HeII_density = &xp->cooling_data.HeII_frac;

  /* HeIII */
  xp->cooling_data.HeIII_frac *= rho;
  data->HeIII_density = &xp->cooling_data.HeIII_frac;

  /* HeII */
  xp->cooling_data.e_frac *= rho;
  data->e_density = &xp->cooling_data.e_frac;
}
#else
void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho) {
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
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 1
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho) {
  /* HM */
  xp->cooling_data.HM_frac *= rho;
  data->HM_density = &xp->cooling_data.HM_frac;

  /* H2I */
  xp->cooling_data.H2I_frac *= rho;
  data->H2I_density = &xp->cooling_data.H2I_frac;

  /* H2II */
  xp->cooling_data.H2II_frac *= rho;
  data->H2II_density = &xp->cooling_data.H2II_frac;
}
#else
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho) {
  data->HM_density = NULL;
  data->H2I_density = NULL;
  data->H2II_density = NULL;
}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE > 2
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho) {
  /* DI */
  xp->cooling_data.DI_frac *= rho;
  data->DI_density = &xp->cooling_data.DI_frac;

  /* DII */
  xp->cooling_data.DII_frac *= rho;
  data->DII_density = &xp->cooling_data.DII_frac;

  /* HDI */
  xp->cooling_data.HDI_frac *= rho;
  data->HDI_density = &xp->cooling_data.HDI_frac;
}
#else
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho) {
  data->DI_density = NULL;
  data->DII_density = NULL;
  data->HDI_density = NULL;
}
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
 * @brief copy the grackle data to a #xpart
 *
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
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
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
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
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
void cooling_copy_to_grackle(grackle_field_data* data, const struct part* p,
                             struct xpart* xp, gr_float rho) {

  cooling_copy_to_grackle1(data, p, xp, rho);
  cooling_copy_to_grackle2(data, p, xp, rho);
  cooling_copy_to_grackle3(data, p, xp, rho);

  data->volumetric_heating_rate = NULL;
  data->specific_heating_rate = NULL;
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
 * @param data The #grackle_field_data
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
void cooling_copy_from_grackle(grackle_field_data* data, const struct part* p,
                               struct xpart* xp, gr_float rho) {
  cooling_copy_from_grackle1(data, p, xp, rho);
  cooling_copy_from_grackle2(data, p, xp, rho);
  cooling_copy_from_grackle3(data, p, xp, rho);

  free(data->metal_density);
}

/**
 * @brief Apply the self shielding (if needed) by turning on/off the UV
 * background.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param chemistry The #chemistry_data from grackle.
 * @param p Pointer to the particle data.
 *
 */
void cooling_apply_self_shielding(
    const struct cooling_function_data* restrict cooling,
    chemistry_data* restrict chemistry, const struct part* restrict p,
    const struct cosmology* cosmo) {

  /* Are we using self shielding or UV background? */
  if (!cooling->with_uv_background || cooling->self_shielding_method >= 0) {
    return;
  }

  /* Are we in a self shielding regime? */
  const float rho = hydro_get_physical_density(p, cosmo);
  if (rho > cooling->self_shielding_threshold) {
    chemistry->UVbackground = 0;
  } else {
    chemistry->UVbackground = 1;
  }
}

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
gr_float cooling_new_energy(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp, double dt) {

  /* set current time */
  code_units units = cooling->units;
  chemistry_data chemistry_grackle = cooling->chemistry;

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
  gr_float energy = hydro_get_physical_internal_energy(p, xp, cosmo) +
                    dt * hydro_get_physical_internal_energy_dt(p, cosmo);

  /* We now need to check that we are not going to go below any of the limits */
  const double u_minimal = hydro_props->minimal_internal_energy;
  energy = max(energy, u_minimal);

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  /* copy to grackle structure */
  cooling_copy_to_grackle(&data, p, xp, density);

  /* Apply the self shielding if requested */
  cooling_apply_self_shielding(cooling, &chemistry_grackle, p, cosmo);

  /* solve chemistry */
  if (local_solve_chemistry(&chemistry_grackle, &grackle_rates, &units, &data,
                            dt) == 0) {
    error("Error in solve_chemistry.");
  }

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(&data, p, xp, density);

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
gr_float cooling_time(const struct phys_const* restrict phys_const,
                      const struct unit_system* restrict us,
                      const struct hydro_props* hydro_props,
                      const struct cosmology* restrict cosmo,
                      const struct cooling_function_data* restrict cooling,
                      const struct part* restrict p,
                      struct xpart* restrict xp) {

  error("TODO: use energy after adiabatic cooling");

  /* set current time */
  code_units units = cooling->units;

  /* initialize data */
  grackle_field_data data;
  chemistry_data chemistry_grackle = cooling->chemistry;

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
  gr_float density = hydro_get_physical_density(p, cosmo);
  gr_float energy = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  /* copy data from particle to grackle data */
  cooling_copy_to_grackle(&data, p, xp, density);

  /* Apply the self shielding if requested */
  cooling_apply_self_shielding(cooling, &chemistry_grackle, p, cosmo);

  /* Compute cooling time */
  gr_float cooling_time;
  chemistry_data_storage chemistry_rates = grackle_rates;
  if (local_calculate_cooling_time(&chemistry_grackle, &chemistry_rates, &units,
                                   &data, &cooling_time) == 0) {
    error("Error in calculate_cooling_time.");
  }

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(&data, p, xp, density);

  /* compute rate */
  return cooling_time;
}

/**
 * @brief Apply the cooling to a particle.
 *
 * Depending on the task order, you may wish to either
 * cool down the particle immediately or do it during the drift.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the xparticle data.
 * @param cosmo The current cosmological model.
 * @param cooling_du_dt Time derivative of the cooling.
 * @param u_new Internal energy after the cooling.
 */
void cooling_apply(struct part* restrict p, struct xpart* restrict xp,
                   const struct cosmology* restrict cosmo, float cooling_du_dt,
                   gr_float u_new) {

#ifdef TASK_ORDER_GEAR
  /* Cannot use du / dt as it will be erased before being used */
  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, u_new);
#else
  hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);
#endif
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param time The current time.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 */
void cooling_cool_part(const struct phys_const* restrict phys_const,
                       const struct unit_system* restrict us,
                       const struct cosmology* restrict cosmo,
                       const struct hydro_props* hydro_props,
                       const struct entropy_floor_properties* floor_props,
                       const struct cooling_function_data* restrict cooling,
                       struct part* restrict p, struct xpart* restrict xp,
                       const double time, const double dt,
                       const double dt_therm) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Is the cooling turn off */
  if (time - xp->cooling_data.time_last_event < cooling->thermal_time) {
    return;
  }

  /* Get the change in internal energy due to hydro forces */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Calculate energy after dt */
  gr_float u_new = cooling_new_energy(phys_const, us, cosmo, hydro_props,
                                      cooling, p, xp, dt);

  /* We now need to check that we are not going to go below any of the limits */
  const double u_minimal = hydro_props->minimal_internal_energy;
  u_new = max(u_new, u_minimal);

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u = u_new - u_old;

  /* Turn this into a rate of change (including cosmology term) */
  const float cooling_du_dt = delta_u / dt_therm;

  /* Update the internal energy time derivative */
  cooling_apply(p, xp, cosmo, cooling_du_dt, u_new);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -=
      hydro_get_mass(p) * (cooling_du_dt - hydro_du_dt) * dt;
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
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {
  // TODO use the grackle library

  /* Physical constants */
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  /* Gas properties */
  const double T_transition = hydro_props->hydrogen_ionization_temperature;
  const double mu_neutral = hydro_props->mu_neutral;
  const double mu_ionised = hydro_props->mu_ionised;

  /* Particle temperature */
  const double u = hydro_get_drifted_physical_internal_energy(p, cosmo);

  /* Temperature over mean molecular weight */
  const double T_over_mu = hydro_gamma_minus_one * u * m_H / k_B;

  /* Are we above or below the HII -> HI transition? */
  if (T_over_mu > (T_transition + 1.) / mu_ionised)
    return T_over_mu * mu_ionised;
  else if (T_over_mu < (T_transition - 1.) / mu_neutral)
    return T_over_mu * mu_neutral;
  else
    return T_transition;
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

  error("Loic: to be implemented");
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
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

  // Set parameter values for chemistry.
  chemistry->use_grackle = 1;
  chemistry->with_radiative_cooling = 1;

  /* molecular network with H, He, D
   From Cloudy table */
  chemistry->primordial_chemistry = cooling->primordial_chemistry;
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
  _free_chemistry_data(&cooling->chemistry, &grackle_rates);
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
