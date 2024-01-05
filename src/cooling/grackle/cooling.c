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
 * @file src/cooling/grackle/cooling.c
 * @brief Cooling using the GRACKLE 3.1.1 library.
 */

#include <config.h>

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

gr_float cooling_new_energy(const struct phys_const* phys_const,
                            const struct unit_system* us,
                            const struct cosmology* cosmo,
                            const struct hydro_props* hydro_properties,
                            const struct cooling_function_data* cooling,
                            const struct part* p, struct xpart* xp, double dt,
                            double dt_therm);

gr_float cooling_time(const struct phys_const* phys_const,
                      const struct unit_system* us,
                      const struct hydro_props* hydro_properties,
                      const struct cosmology* cosmo,
                      const struct cooling_function_data* cooling,
                      const struct part* p, struct xpart* xp);

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param phys_const The #phys_const.
 * @param cosmo The current cosmological model.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 * @param time The current system time
 */
void cooling_update(const struct phys_const* phys_const,
                    const struct cosmology* cosmo,
                    const struct pressure_floor_props* pressure_floor,
                    struct cooling_function_data* cooling, struct space* s,
                    const double time) {
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
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param hydro_properties The #hydro_props.
 * @param cosmo The #cosmology
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_compute_equilibrium(const struct phys_const* phys_const,
                                 const struct unit_system* us,
                                 const struct hydro_props* hydro_properties,
                                 const struct cosmology* cosmo,
                                 const struct cooling_function_data* cooling,
                                 const struct part* p, struct xpart* xp) {

  return;
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
void cooling_first_init_part(const struct phys_const* phys_const,
                             const struct unit_system* us,
                             const struct hydro_props* hydro_props,
                             const struct cosmology* cosmo,
                             const struct cooling_function_data* cooling,
                             const struct part* p, struct xpart* xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.time_last_event = -cooling->thermal_time;

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;

  /* NOTE: at this stage, we assume neutral gas
   * a better determination will be done in cooling_post_init_part */

  /* primordial chemistry >= 1 */
  /* assume neutral gas */
  xp->cooling_data.HI_frac = grackle_data->HydrogenFractionByMass;
  xp->cooling_data.HII_frac = zero;
  xp->cooling_data.HeI_frac = 1. - grackle_data->HydrogenFractionByMass;
  xp->cooling_data.HeII_frac = zero;
  xp->cooling_data.HeIII_frac = zero;
  xp->cooling_data.e_frac = zero;
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
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state. The function requires the density to be defined and thus must
 * be called after its computation.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_post_init_part(const struct phys_const* phys_const,
                            const struct unit_system* us,
                            const struct hydro_props* hydro_props,
                            const struct cosmology* cosmo,
                            const struct cooling_function_data* cooling,
                            const struct part* p, struct xpart* xp) {

  // const float rho = hydro_get_physical_density(p, cosmo);
  // const float energy = hydro_get_physical_internal_energy(p, xp, cosmo);
  // message("rho = %g energy = %g",rho,energy);

#if COOLING_GRACKLE_MODE > 0
  /* The function below currently does nothing. Will have to be updated. */
  cooling_compute_equilibrium(phys_const, us, hydro_props, cosmo, cooling, p,
                              xp);
#endif
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
float cooling_get_radiated_energy(const struct xpart* xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print_backend(const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle'.");
  message("Using Grackle = %i", cooling->chemistry_data.use_grackle);
  message("Chemical network = %i",
          cooling->chemistry_data.primordial_chemistry);
  message("CloudyTable = %s", cooling->cloudy_table);
  message("Redshift = %g", cooling->redshift);
  message("UV background = %d", cooling->with_uv_background);
  message("Metal cooling = %i", cooling->chemistry_data.metal_cooling);
  message("Self Shielding = %i", cooling->self_shielding_method);
  if (cooling->self_shielding_method == -1) {
    message("Self Shelding density = %g", cooling->self_shielding_threshold);
  }

  message("Thermal time = %g", cooling->thermal_time);
  message("Specific Heating Rates = %g", cooling->specific_heating_rates);
  message("Volumetric Heating Rates = %g", cooling->volumetric_heating_rates);

  message("grackle_chemistry_data.RT_heating_rate = %g",
          cooling->RT_heating_rate);
  message("grackle_chemistry_data.RT_HI_ionization_rate = %g",
          cooling->RT_HI_ionization_rate);
  message("grackle_chemistry_data.RT_HeI_ionization_rate = %g",
          cooling->RT_HeI_ionization_rate);
  message("grackle_chemistry_data.RT_HeII_ionization_rate = %g",
          cooling->RT_HeII_ionization_rate);
  message("grackle_chemistry_data.RT_H2_dissociation_rate = %g",
          cooling->RT_H2_dissociation_rate);

  message("Units:");
  message("\tComoving = %i", cooling->units.comoving_coordinates);
  message("\tLength = %g", cooling->units.length_units);
  message("\tDensity = %g", cooling->units.density_units);
  message("\tTime = %g", cooling->units.time_units);
  message("\tScale Factor = %g (units: %g)", cooling->units.a_value,
          cooling->units.a_units);

  message("Grackle parameters:");
  message("grackle_chemistry_data.use_grackle = %d",
          cooling->chemistry_data.use_grackle);
  message("grackle_chemistry_data.with_radiative_cooling %d",
          cooling->chemistry_data.with_radiative_cooling);
  message("grackle_chemistry_data.primordial_chemistry = %d",
          cooling->chemistry_data.primordial_chemistry);
  message("grackle_chemistry_data.three_body_rate = %d",
          cooling->chemistry_data.three_body_rate);
  message("grackle_chemistry_data.cmb_temperature_floor = %d",
          cooling->chemistry_data.cmb_temperature_floor);
  message("grackle_chemistry_data.cie_cooling = %d",
          cooling->chemistry_data.cie_cooling);
  message("grackle_chemistry_data.dust_chemistry = %d",
          cooling->chemistry_data.dust_chemistry);
  message("grackle_chemistry_data.metal_cooling = %d",
          cooling->chemistry_data.metal_cooling);
  message("grackle_chemistry_data.UVbackground = %d",
          cooling->chemistry_data.UVbackground);
  message("grackle_chemistry_data.CaseBRecombination = %d",
          cooling->chemistry_data.CaseBRecombination);
  message("grackle_chemistry_data.grackle_data_file = %s",
          cooling->chemistry_data.grackle_data_file);
  message("grackle_chemistry_data.use_radiative_transfer = %d",
          cooling->chemistry_data.use_radiative_transfer);
  message("grackle_chemistry_data.use_volumetric_heating_rate = %d",
          cooling->chemistry_data.use_volumetric_heating_rate);
  message("grackle_chemistry_data.use_specific_heating_rate = %d",
          cooling->chemistry_data.use_specific_heating_rate);
  message("grackle_chemistry_data.self_shielding_method = %d",
          cooling->chemistry_data.self_shielding_method);
  message("grackle_chemistry_data.HydrogenFractionByMass = %.3g",
          cooling->chemistry_data.HydrogenFractionByMass);
  message("grackle_chemistry_data.Gamma = %.6g", cooling->chemistry_data.Gamma);
  message("grackle_chemistry_data.cie_cooling = %d",
          cooling->chemistry_data.cie_cooling);
  message("grackle_chemistry_data.three_body_rate = %d",
          cooling->chemistry_data.three_body_rate);
  message("grackle_chemistry_data.h2_on_dust = %d",
          cooling->chemistry_data.h2_on_dust);
  message("grackle_chemistry_data.use_dust_density_field = %d",
          cooling->chemistry_data.use_dust_density_field);
  message("grackle_chemistry_data.local_dust_to_gas_ratio = %.3g",
          cooling->chemistry_data.local_dust_to_gas_ratio);
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
void cooling_copy_to_grackle(grackle_field_data* data, const struct part* p,
                             struct xpart* xp, gr_float rho,
                             gr_float species_densities[12],
                             const struct cooling_function_data* cooling,
                             const struct phys_const* phys_const) {

  const float time_units = cooling->units.time_units;

  cooling_copy_to_grackle1(data, p, xp, rho, species_densities);
  cooling_copy_to_grackle2(data, p, xp, rho, species_densities);
  cooling_copy_to_grackle3(data, p, xp, rho, species_densities);

  if (cooling->chemistry_data.use_volumetric_heating_rate) {
    gr_float* volumetric_heating_rate = (gr_float*)malloc(sizeof(gr_float));
    *volumetric_heating_rate = cooling->volumetric_heating_rates;
    data->volumetric_heating_rate = volumetric_heating_rate;
  }

  if (cooling->chemistry_data.use_specific_heating_rate) {
    gr_float* specific_heating_rate = (gr_float*)malloc(sizeof(gr_float));
    *specific_heating_rate = cooling->specific_heating_rates;
    data->specific_heating_rate = specific_heating_rate;
  }

  if (cooling->chemistry_data.use_radiative_transfer) {

    /* heating rate */
    gr_float* RT_heating_rate = (gr_float*)malloc(sizeof(gr_float));
    *RT_heating_rate = cooling->RT_heating_rate;
    /* Note to self:
     * If cooling->RT_heating_rate is computed properly, i.e. using
     * the HI density, and then being HI density dependent, we need
     * to divide it as follow. If it is assumed to be already normed
     * as it is so when providing it via some parameters, we keep it
     * unchanged.
     */
    /* Grackle wants heating rate in units of / nHI_cgs */
    // const double nHI_cgs = species_densities[0]
    //                      / phys_const->const_proton_mass
    //                      / pow(length_units,3);
    //*RT_heating_rate /= nHI_cgs;
    data->RT_heating_rate = RT_heating_rate;

    /* HI ionization rate */
    gr_float* RT_HI_ionization_rate = (gr_float*)malloc(sizeof(gr_float));
    *RT_HI_ionization_rate = cooling->RT_HI_ionization_rate;
    /* Grackle wants it in 1/internal_time_units */
    *RT_HI_ionization_rate /= (1. / time_units);
    data->RT_HI_ionization_rate = RT_HI_ionization_rate;

    /* HeI ionization rate */
    gr_float* RT_HeI_ionization_rate = (gr_float*)malloc(sizeof(gr_float));
    *RT_HeI_ionization_rate = cooling->RT_HeI_ionization_rate;
    /* Grackle wants it in 1/internal_time_units */
    *RT_HeI_ionization_rate /= (1. / time_units);
    data->RT_HeI_ionization_rate = RT_HeI_ionization_rate;

    /* HeII ionization rate */
    gr_float* RT_HeII_ionization_rate = (gr_float*)malloc(sizeof(gr_float));
    *RT_HeII_ionization_rate = cooling->RT_HeII_ionization_rate;
    /* Grackle wants it in 1/internal_time_units */
    *RT_HeII_ionization_rate /= (1. / time_units);
    data->RT_HeII_ionization_rate = RT_HeII_ionization_rate;

    /* H2 ionization rate */
    gr_float* RT_H2_dissociation_rate = (gr_float*)malloc(sizeof(gr_float));
    *RT_H2_dissociation_rate = cooling->RT_H2_dissociation_rate;
    /* Grackle wants it in 1/internal_time_units */
    *RT_H2_dissociation_rate /= (1. / time_units);
    data->RT_H2_dissociation_rate = RT_H2_dissociation_rate;

  } else {
    data->volumetric_heating_rate = NULL;
    data->specific_heating_rate = NULL;
    data->RT_heating_rate = NULL;
    data->RT_HI_ionization_rate = NULL;
    data->RT_HeI_ionization_rate = NULL;
    data->RT_HeII_ionization_rate = NULL;
    data->RT_H2_dissociation_rate = NULL;
  }

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
                               struct xpart* xp, gr_float rho,
                               const struct cooling_function_data* cooling) {
  cooling_copy_from_grackle1(data, p, xp, rho);
  cooling_copy_from_grackle2(data, p, xp, rho);
  cooling_copy_from_grackle3(data, p, xp, rho);

  if (cooling->chemistry_data.use_volumetric_heating_rate)
    free(data->volumetric_heating_rate);
  if (cooling->chemistry_data.use_specific_heating_rate)
    free(data->specific_heating_rate);

  if (cooling->chemistry_data.use_radiative_transfer) {
    free(data->RT_heating_rate);
    free(data->RT_HI_ionization_rate);
    free(data->RT_HeI_ionization_rate);
    free(data->RT_HeII_ionization_rate);
    free(data->RT_H2_dissociation_rate);
  }

  free(data->metal_density);
}

/**
 * @brief Apply the self shielding (if needed) by turning on/off the UV
 * background.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param chemistry The chemistry_data structure from grackle.
 * @param p Pointer to the particle data.
 * @param cosmo The #cosmology.
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
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 *
 * @return du / dt
 */
gr_float cooling_new_energy(const struct phys_const* phys_const,
                            const struct unit_system* us,
                            const struct cosmology* cosmo,
                            const struct hydro_props* hydro_props,
                            const struct cooling_function_data* cooling,
                            const struct part* p, struct xpart* xp, double dt,
                            double dt_therm) {

  /* set current time */
  code_units units = cooling->units;
  chemistry_data chemistry_grackle = cooling->chemistry_data;
  chemistry_data_storage rates_grackle = cooling->chemistry_rates;

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
                    dt_therm * hydro_get_physical_internal_energy_dt(p, cosmo);
  energy = max(energy, hydro_props->minimal_internal_energy);
  /* keep this array here so you can point to and from it in
   * copy to/from grackle */
  gr_float species_densities[12];

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  /* copy to grackle structure */
  cooling_copy_to_grackle(&data, p, xp, density, species_densities, cooling,
                          phys_const);

  /* Apply the self shielding if requested */
  cooling_apply_self_shielding(cooling, &chemistry_grackle, p, cosmo);

  /* solve chemistry */
  if (local_solve_chemistry(&chemistry_grackle, &rates_grackle, &units, &data,
                            dt) == 0) {
    error("Error in solve_chemistry.");
  }

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(&data, p, xp, density, cooling);

  return energy;
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
gr_float cooling_time(const struct phys_const* phys_const,
                      const struct unit_system* us,
                      const struct hydro_props* hydro_props,
                      const struct cosmology* cosmo,
                      const struct cooling_function_data* cooling,
                      const struct part* p, struct xpart* xp) {

  /* set current time */
  code_units units = cooling->units;

  /* initialize data */
  grackle_field_data data;
  chemistry_data chemistry_grackle = cooling->chemistry_data;
  chemistry_data_storage rates_grackle = cooling->chemistry_rates;

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
  energy = max(energy, hydro_props->minimal_internal_energy);

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  gr_float species_densities[12];
  /* copy data from particle to grackle data */
  cooling_copy_to_grackle(&data, p, xp, density, species_densities, cooling,
                          phys_const);

  /* Apply the self shielding if requested */
  cooling_apply_self_shielding(cooling, &chemistry_grackle, p, cosmo);

  /* Compute cooling time */
  gr_float cooling_time;
  if (local_calculate_cooling_time(&chemistry_grackle, &rates_grackle, &units,
                                   &data, &cooling_time) == 0) {
    error("Error in calculate_cooling_time.");
  }

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(&data, p, xp, density, cooling);

  /* compute rate */
  return cooling_time;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param floor_props Properties of the entropy floor.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
void cooling_cool_part(const struct phys_const* phys_const,
                       const struct unit_system* us,
                       const struct cosmology* cosmo,
                       const struct hydro_props* hydro_props,
                       const struct entropy_floor_properties* floor_props,
                       const struct pressure_floor_props* pressure_floor,
                       const struct cooling_function_data* cooling,
                       struct part* p, struct xpart* xp, const double dt,
                       const double dt_therm, const double time) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Energy after the adiabatic cooling */
  float u_ad_before =
      u_old + dt_therm * hydro_get_physical_internal_energy_dt(p, cosmo);

  /* We now need to check that we are not going to go below any of the limits */
  const double u_minimal = hydro_props->minimal_internal_energy;
  if (u_ad_before < u_minimal) {
    u_ad_before = u_minimal;
    const float du_dt = (u_ad_before - u_old) / dt_therm;
    hydro_set_physical_internal_energy_dt(p, cosmo, du_dt);
  }

  /* Calculate energy after dt */
  gr_float u_new = 0;

  /* Is the cooling turn off */
  if (time - xp->cooling_data.time_last_event < cooling->thermal_time) {
    u_new = u_ad_before;
  } else {
    u_new = cooling_new_energy(phys_const, us, cosmo, hydro_props, cooling, p,
                               xp, dt, dt_therm);
  }

  /* Get the change in internal energy due to hydro forces */
  float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* We now need to check that we are not going to go below any of the limits */
  u_new = max(u_new, u_minimal);

  /* Calculate the cooling rate */
  float cool_du_dt = (u_new - u_ad_before) / dt_therm;
  float du_dt = cool_du_dt + hydro_du_dt;

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cool_du_dt * dt_therm;
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
float cooling_get_temperature(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp) {
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
 * We return FLT_MAX so as to impose no limit on the time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The #cosmology.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 */
float cooling_timestep(const struct cooling_function_data* cooling,
                       const struct phys_const* phys_const,
                       const struct cosmology* cosmo,
                       const struct unit_system* us,
                       const struct hydro_props* hydro_props,
                       const struct part* p, const struct xpart* xp) {

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

  chemistry_data* chemistry = &cooling->chemistry_data;

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
  chemistry->three_body_rate = cooling->H2_three_body_rate;
  chemistry->cmb_temperature_floor = cooling->cmb_temperature_floor;
  chemistry->cie_cooling = cooling->H2_cie_cooling;
  chemistry->h2_on_dust = cooling->H2_on_dust;
  chemistry->grackle_data_file = cooling->cloudy_table;

  if (cooling->local_dust_to_gas_ratio > 0)
    chemistry->local_dust_to_gas_ratio = cooling->local_dust_to_gas_ratio;

    /* radiative transfer */
#if COOLING_GRACKLE_MODE == 0
  if (cooling->use_radiative_transfer)
    error(
        "The parameter use_radiative_transfer cannot be set to 1 in Grackle "
        "mode 0 !");
#endif

  chemistry->use_radiative_transfer = cooling->use_radiative_transfer;

  if (cooling->volumetric_heating_rates > 0)
    chemistry->use_volumetric_heating_rate = 1;

  if (cooling->specific_heating_rates > 0)
    chemistry->use_specific_heating_rate = 1;

  /* hydrogen fraction by mass */
  chemistry->HydrogenFractionByMass = cooling->HydrogenFractionByMass;

  /* use the Case B recombination rates */
  chemistry->CaseBRecombination = 1;

  if (cooling->specific_heating_rates > 0 &&
      cooling->volumetric_heating_rates > 0)
    error(
        "You should specified either the specific or the volumetric "
        "heating rates, not both");

  /* self shielding */
  if (cooling->self_shielding_method <= 0)
    chemistry->self_shielding_method = 0;
  else
    chemistry->self_shielding_method = cooling->self_shielding_method;

  if (local_initialize_chemistry_data(&cooling->chemistry_data,
                                      &cooling->chemistry_rates,
                                      &cooling->units) == 0) {
    error("Error in initialize_chemistry_data");
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
  /* Clean up grackle data. This is a call to a grackle function */
  local_free_chemistry_data(&cooling->chemistry_data,
                            &cooling->chemistry_rates);
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
