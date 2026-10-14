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
 * @file src/cooling/geajay/cooling.c
 * @brief Cooling using the geajay library.
 */

#include <config.h>

/* Include header */
#include "cooling.h"

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* The geajay library itself */
#include <geajay_message.h>
#include <geajay_particle.h>
#include <grackle.h>

/* Local includes. */
#include "atomic.h"
#include "chemistry.h"
#include "cooling_io.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* TODO: Move these to cooling.h */
gr_float cooling_new_energy(const struct phys_const *phys_const,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct hydro_props *hydro_properties,
                            const struct pressure_floor_props *pressure_floor,
                            const struct cooling_function_data *cooling,
                            struct part *p, struct xpart *xp, double dt,
                            double dt_therm, double time);

double cooling_get_physical_density(
    const struct part *p, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling);

/*! geajay solves on this rank that failed after every retry since the last
 * report. Each kept the particle's previous state. */
static long long cooling_geajay_failed_solves = 0;

/*! geajay solves on this rank that failed once and then succeeded in a
 * retry since the last report. */
static long long cooling_geajay_recovered_solves = 0;

/*! Set by the first particle above the density range of the Cloudy table. */
static int cooling_geajay_table_edge_warned = 0;

/**
 * @brief Return the numbers of failed and of recovered geajay solves on this
 * rank since the last call, and reset both to zero.
 *
 * @param failed (return) Solves that kept the particle's previous state.
 * @param recovered (return) Solves that succeeded in a retry.
 */
void cooling_pop_geajay_solve_failure_counts(long long *failed,
                                             long long *recovered) {
  *failed = atomic_swap(&cooling_geajay_failed_solves, 0);
  *recovered = atomic_swap(&cooling_geajay_recovered_solves, 0);
}

/**
 * @brief Test that a stored value is finite and, optionally, strictly
 * positive, from its IEEE-754 bits.
 *
 * -ffast-math marks a float parameter passed by value nofpclass(nan inf),
 * making the bit test legally foldable even in a noinline function; a
 * pointer to the stored field carries no such annotation.
 *
 * @param x Pointer to the stored value.
 * @param positive Also require *x > 0.
 * @return 1 if the test passes, 0 otherwise.
 */
__attribute__((noinline)) static int cooling_stored_value_is_valid(
    const float *x, const int positive) {
  uint32_t bits;
  memcpy(&bits, x, sizeof(bits));
  if (((bits >> 23) & 0xffU) == 0xffU) return 0;
  if (!positive) return 1;
  return (bits >> 31) == 0 && bits != 0;
}

/**
 * @brief Record the time when cooling was switched off for a particle.
 *
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @param time The time when the cooling was switched off.
 */
INLINE void cooling_set_part_time_cooling_off(struct part *p, struct xpart *xp,
                                              const double time) {
  xp->cooling_data.time_last_event = time;
}

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
void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time) {
  /* set current time */
  if (cooling->redshift == -1)
    cooling->units.a_value = cosmo->a;
  else
    cooling->units.a_value = 1. / (1. + cooling->redshift);

  if (geajay_context_update(cooling->context, &cooling->units) == GR_FAIL)
    error("Error in geajay_context_update.");

  /* Report the geajay failures of the previous step */
  long long failed, recovered;
  cooling_pop_geajay_solve_failure_counts(&failed, &recovered);
  if (failed > 0 || recovered > 0)
    message(
        "geajay solves in the previous step: %lld failed (particles kept "
        "their previous state), %lld recovered by retrying.",
        failed, recovered);
}

/**
 * @brief Print the chemical network
 *
 * @param xp The #xpart to print
 */
void cooling_print_fractions(const struct xpart *restrict xp) {

  const struct cooling_xpart_data tmp = xp->cooling_data;
#if COOLING_GEAJAY_MODE > 0
  message("HI %g, HII %g, HeI %g, HeII %g, HeIII %g, e %g", tmp.HI_frac,
          tmp.HII_frac, tmp.HeI_frac, tmp.HeII_frac, tmp.HeIII_frac,
          tmp.e_frac);
#endif

#if COOLING_GEAJAY_MODE > 1
  message("HM %g, H2I %g, H2II %g", tmp.HM_frac, tmp.H2I_frac, tmp.H2II_frac);
#endif

#if COOLING_GEAJAY_MODE > 2
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
int cooling_converged(const struct xpart *restrict xp,
                      const struct xpart *restrict old, const float limit) {

#if COOLING_GEAJAY_MODE > 0
  cooling_check_field(xp, old, HI_frac, limit);
  cooling_check_field(xp, old, HII_frac, limit);
  cooling_check_field(xp, old, HeI_frac, limit);
  cooling_check_field(xp, old, HeII_frac, limit);
  cooling_check_field(xp, old, HeIII_frac, limit);
  cooling_check_field(xp, old, e_frac, limit);
#endif
#if COOLING_GEAJAY_MODE > 1
  cooling_check_field(xp, old, HM_frac, limit);
  cooling_check_field(xp, old, H2I_frac, limit);
  cooling_check_field(xp, old, H2II_frac, limit);
#endif

#if COOLING_GEAJAY_MODE > 2
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
void cooling_compute_equilibrium(const struct phys_const *phys_const,
                                 const struct unit_system *us,
                                 const struct hydro_props *hydro_properties,
                                 const struct cosmology *cosmo,
                                 const struct cooling_function_data *cooling,
                                 const struct part *p, struct xpart *xp) {

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
void cooling_first_init_part(const struct phys_const *phys_const,
                             const struct unit_system *us,
                             const struct hydro_props *hydro_props,
                             const struct cosmology *cosmo,
                             const struct cooling_function_data *cooling,
                             const struct part *p, struct xpart *xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.time_last_event = -cooling->thermal_time;

#if COOLING_GEAJAY_MODE >= 1
  gr_float zero = 1.e-20;

  /* NOTE: if the ratio with respect to hydrogen is given, use it. Instead,
   * we assume neutral gas.
   * A better determination of the abundances can be done in
   * cooling_post_init_part */

  /* Compute nH (formally divided by the gas density and assuming the proton
   * mass to be one) */
  const int has_context = (cooling->context != NULL);
  const double nH =
      has_context ? cooling->chemistry_data.HydrogenFractionByMass : 1.;

  /* Compute nHe (formally divided by the gas density and assuming the proton
   * mass to be one) */
  const double nHe =
      has_context ? (1. - cooling->chemistry_data.HydrogenFractionByMass) / 4.
                  : 0.;

  /* Electron density */
  double ne = zero;

  /* Sum of mass fraction (to test its consistency) */
  double Xtot = 0.f;

  /* primordial chemistry >= 1 */

  /* Hydrogen I */
  if (cooling->initial_nHII_to_nH_ratio >= 0.f)
    xp->cooling_data.HI_frac = nH * (1.f - cooling->initial_nHII_to_nH_ratio);
  else
    xp->cooling_data.HI_frac = nH;

  Xtot += xp->cooling_data.HI_frac;

  /* Hydrogen II */
  if (cooling->initial_nHII_to_nH_ratio >= 0.f) {
    double nHII = nH * cooling->initial_nHII_to_nH_ratio;
    xp->cooling_data.HII_frac = nHII;
    ne += nHII;
  } else
    xp->cooling_data.HII_frac = zero;

  Xtot += xp->cooling_data.HII_frac;

  /* Helium I */
  if (cooling->initial_nHeI_to_nH_ratio >= 0.f) {
    double nHeI = nH * cooling->initial_nHeI_to_nH_ratio;
    xp->cooling_data.HeI_frac = nHeI * 4.f;
  } else
    xp->cooling_data.HeI_frac = nHe * 4.f;

  Xtot += xp->cooling_data.HeI_frac;

  /* Helium II */
  if (cooling->initial_nHeII_to_nH_ratio >= 0.f) {
    double nHeII = nH * cooling->initial_nHeII_to_nH_ratio;
    xp->cooling_data.HeII_frac = nHeII * 4.f;
    ne += nHeII;
  } else
    xp->cooling_data.HeII_frac = zero;

  Xtot += xp->cooling_data.HeII_frac;

  /* Helium III */
  if (cooling->initial_nHeIII_to_nH_ratio >= 0.f) {
    double nHeIII = nH * cooling->initial_nHeIII_to_nH_ratio;
    xp->cooling_data.HeIII_frac = nHeIII * 4.f;
    ne += 2.f * nHeIII;
  } else
    xp->cooling_data.HeIII_frac = zero;

  Xtot += xp->cooling_data.HeIII_frac;

  /* electron mass fraction (multiplied by the proton mass (Grackle convention)
   */
  xp->cooling_data.e_frac = ne;

#endif  // MODE >= 1

#if COOLING_GEAJAY_MODE >= 2
  /* primordial chemistry >= 2 */

  /* Hydrogen- */
  if (cooling->initial_nHM_to_nH_ratio >= 0.f) {
    double nHM = nH * cooling->initial_nHM_to_nH_ratio;
    xp->cooling_data.HM_frac = nHM;
    ne -= nHM;
  } else
    xp->cooling_data.HM_frac = zero;

  Xtot += xp->cooling_data.HM_frac;

  /* H2I */
  if (cooling->initial_nH2I_to_nH_ratio >= 0.f) {
    double nH2I = nH * cooling->initial_nH2I_to_nH_ratio;
    xp->cooling_data.H2I_frac = nH2I * 2.f;
  } else
    xp->cooling_data.H2I_frac = zero;

  Xtot += xp->cooling_data.H2I_frac;

  /* H2II */
  if (cooling->initial_nH2II_to_nH_ratio >= 0.f) {
    double nH2II = nH * cooling->initial_nH2II_to_nH_ratio;
    xp->cooling_data.H2II_frac = nH2II * 2.f;
    ne += nH2II;
  } else
    xp->cooling_data.H2II_frac = zero;

  Xtot += xp->cooling_data.H2II_frac;

  /* electron mass fraction (multiplied by the proton mass (Grackle convention)
   */
  xp->cooling_data.e_frac = ne;

#endif  // MODE >= 2

#if COOLING_GEAJAY_MODE >= 3
  /* primordial chemistry >= 3 */

  /* Deuterium I */
  if (cooling->initial_nDI_to_nH_ratio >= 0.f) {
    double nDI = nH * cooling->initial_nDI_to_nH_ratio;
    xp->cooling_data.DI_frac = nDI * 2.f;
  } else {
    xp->cooling_data.DI_frac =
        cooling->chemistry_data.DeuteriumToHydrogenRatio *
        cooling->chemistry_data.HydrogenFractionByMass;
  }

  Xtot += xp->cooling_data.DI_frac;

  /* Deuterium II */
  if (cooling->initial_nDII_to_nH_ratio >= 0.f) {
    double nDII = nH * cooling->initial_nDII_to_nH_ratio;
    xp->cooling_data.DII_frac = nDII * 2.f;
    ne += nDII;
  } else
    xp->cooling_data.DII_frac = zero;

  Xtot += xp->cooling_data.DII_frac;

  /* HD I */
  if (cooling->initial_nHDI_to_nH_ratio >= 0.f) {
    double nHDI = nH * cooling->initial_nHDI_to_nH_ratio;
    xp->cooling_data.HDI_frac = nHDI * 3.f;
  } else
    xp->cooling_data.HDI_frac = zero;

  Xtot += xp->cooling_data.HDI_frac;

  /* electron mass fraction (multiplied by the proton mass (Grackle convention)
   */
  xp->cooling_data.e_frac = ne;

#endif  // MODE >= 3

#if COOLING_GEAJAY_MODE >= 1
  /* Consistency check: the species mass fractions accumulated above
     (whichever tiers are compiled in for this COOLING_GEAJAY_MODE) must
     sum to the total gas mass fraction. */
  if (fabs(Xtot - 1.0) > 1e-3)
    error("Got total mass fraction of gas = %.6g", Xtot);
#endif  // MODE >= 1
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
void cooling_post_init_part(const struct phys_const *phys_const,
                            const struct unit_system *us,
                            const struct hydro_props *hydro_props,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            const struct part *p, struct xpart *xp) {

#if COOLING_GEAJAY_MODE > 0
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
float cooling_get_radiated_energy(const struct xpart *xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  if (engine_rank != 0) {
    return;
  }

  message("Cooling function is 'geajay'.");
  message("Using geajay = %i", cooling->chemistry_data.use_grackle);
  message("Chemical network = %i",
          cooling->chemistry_data.primordial_chemistry);
  message("CloudyTable = %s", cooling->cloudy_table);
  message("Redshift = %g", cooling->redshift);
  message("UV background = %d", cooling->with_uv_background);
  message("Metal cooling = %i", cooling->chemistry_data.metal_cooling);
  message("Self Shielding = %i", cooling->self_shielding_method);
  message(
      "Cloudy table density range = [%g, %g] cm^-3, beyond it the lookup %s",
      cooling->table_density_min_Hpcm3, cooling->table_density_max_Hpcm3,
      cooling->table_density_edge == GEAJAY_TABLE_EDGE_HOLD
          ? "holds the edge value (table_density_edge: edge)"
          : "extrapolates (table_density_edge: extrapolate)");
  message("AGORA CMB temperature floor = %d",
          cooling->agora_cmb_temperature_floor);
  if (cooling->self_shielding_method == -1) {
    message("Self Shelding density = %g", cooling->self_shielding_threshold);
  }

  message("Thermal time = %g", cooling->thermal_time);
  message("Max geajay sub-steps = %d", cooling->max_step);
  message("geajay retries on failure = %d", cooling->subcycle_on_failure);
  if (cooling->disable_cooling_for_debugging)
    warning(
        "GeajayCooling:disable_cooling_for_debugging is set: the "
        "per-particle cooling/heating update is skipped entirely. Never "
        "use this in a production run.");
  message("Specific Heating Rates = %g", cooling->specific_heating_rates);
  message("Volumetric Heating Rates = %g", cooling->volumetric_heating_rates);

  message("chemistry_data.RT_heating_rate = %g", cooling->RT_heating_rate);
  message("chemistry_data.RT_HI_ionization_rate = %g",
          cooling->RT_HI_ionization_rate);
  message("chemistry_data.RT_HeI_ionization_rate = %g",
          cooling->RT_HeI_ionization_rate);
  message("chemistry_data.RT_HeII_ionization_rate = %g",
          cooling->RT_HeII_ionization_rate);
  message("chemistry_data.RT_H2_dissociation_rate = %g",
          cooling->RT_H2_dissociation_rate);
  message("cooling.initial_nHII_to_nH_ratio= %g",
          cooling->initial_nHII_to_nH_ratio);
  message("cooling.initial_nHeI_to_nH_ratio= %g",
          cooling->initial_nHeI_to_nH_ratio);
  message("cooling.initial_nHeII_to_nH_ratio= %g",
          cooling->initial_nHeII_to_nH_ratio);
  message("cooling.initial_nHeIII_to_nH_ratio= %g",
          cooling->initial_nHeIII_to_nH_ratio);
  message("cooling.initial_nDI_to_nH_ratio= %g",
          cooling->initial_nDI_to_nH_ratio);
  message("cooling.initial_nDII_to_nH_ratio= %g",
          cooling->initial_nDII_to_nH_ratio);
  message("cooling.initial_nHM_to_nH_ratio= %g",
          cooling->initial_nHM_to_nH_ratio);
  message("cooling.initial_nH2I_to_nH_ratio= %g",
          cooling->initial_nH2I_to_nH_ratio);
  message("cooling.initial_nH2II_to_nH_ratio= %g",
          cooling->initial_nH2II_to_nH_ratio);
  message("cooling.initial_nHDI_to_nH_ratio= %g",
          cooling->initial_nHDI_to_nH_ratio);
  message("Units:");
  message("\tComoving = %i", cooling->units.comoving_coordinates);
  message("\tLength = %g", cooling->units.length_units);
  message("\tDensity = %g", cooling->units.density_units);
  message("\tTime = %g", cooling->units.time_units);
  message("\tScale Factor = %g (units: %g)", cooling->units.a_value,
          cooling->units.a_units);

  message("geajay parameters:");
  message("chemistry_data.use_grackle = %d",
          cooling->chemistry_data.use_grackle);
  message("chemistry_data.with_radiative_cooling %d",
          cooling->chemistry_data.with_radiative_cooling);
  message("chemistry_data.primordial_chemistry = %d",
          cooling->chemistry_data.primordial_chemistry);
  message("chemistry_data.three_body_rate = %d",
          cooling->chemistry_data.three_body_rate);
  message("chemistry_data.cmb_temperature_floor = %d",
          cooling->chemistry_data.cmb_temperature_floor);
  message("chemistry_data.cie_cooling = %d",
          cooling->chemistry_data.cie_cooling);
  message("chemistry_data.dust_chemistry = %d",
          cooling->chemistry_data.dust_chemistry);
  message("chemistry_data.metal_cooling = %d",
          cooling->chemistry_data.metal_cooling);
  message("chemistry_data.UVbackground = %d",
          cooling->chemistry_data.UVbackground);
  message("chemistry_data.CaseBRecombination = %d",
          cooling->chemistry_data.CaseBRecombination);
  message("chemistry_data.grackle_data_file = %s",
          cooling->chemistry_data.grackle_data_file);
  message("chemistry_data.use_radiative_transfer = %d",
          cooling->chemistry_data.use_radiative_transfer);
  message("chemistry_data.use_volumetric_heating_rate = %d",
          cooling->chemistry_data.use_volumetric_heating_rate);
  message("chemistry_data.use_specific_heating_rate = %d",
          cooling->chemistry_data.use_specific_heating_rate);
  message("chemistry_data.self_shielding_method = %d",
          cooling->chemistry_data.self_shielding_method);
  message("chemistry_data.H2_self_shielding = %d",
          cooling->chemistry_data.H2_self_shielding);
  message("chemistry_data.HydrogenFractionByMass = %.3g",
          cooling->chemistry_data.HydrogenFractionByMass);
  message("chemistry_data.Gamma = %.6g", cooling->chemistry_data.Gamma);
  message("chemistry_data.h2_on_dust = %d", cooling->chemistry_data.h2_on_dust);
  message("chemistry_data.use_dust_density_field = %d",
          cooling->chemistry_data.use_dust_density_field);
  message("chemistry_data.local_dust_to_gas_ratio = %.3g",
          cooling->chemistry_data.local_dust_to_gas_ratio);
}

/**
 * @brief Fill the inputs of one geajay call for a particle.
 *
 * The species are not part of the inputs: geajay reads them in place in the
 * #xpart, through the layout of the cooling.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The #part.
 * @param cosmo The #cosmology.
 * @param density The physical density passed to geajay.
 * @param energy The physical specific internal energy passed to geajay.
 * @param inputs (return) The inputs.
 */
static void cooling_fill_inputs(const struct cooling_function_data *cooling,
                                const struct part *p,
                                const struct cosmology *cosmo,
                                const double density, const double energy,
                                geajay_particle_inputs *inputs) {

  const float time_units = cooling->units.time_units;

  geajay_particle_inputs_init(cooling->context, inputs);

  inputs->density = density;
  inputs->internal_energy = energy;

  inputs->has_metal = 1;
  inputs->metal =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) * density;

  if (cooling->chemistry_data.use_volumetric_heating_rate)
    inputs->volumetric_heating_rate = cooling->volumetric_heating_rates;

  if (cooling->chemistry_data.use_specific_heating_rate)
    inputs->specific_heating_rate = cooling->specific_heating_rates;

  /* The rates are given in 1/s: the library wants 1/internal_time_units */
  if (cooling->chemistry_data.use_radiative_transfer) {
    inputs->RT_heating_rate = cooling->RT_heating_rate;
    inputs->RT_HI_ionization_rate =
        cooling->RT_HI_ionization_rate / (1. / time_units);
    inputs->RT_HeI_ionization_rate =
        cooling->RT_HeI_ionization_rate / (1. / time_units);
    inputs->RT_HeII_ionization_rate =
        cooling->RT_HeII_ionization_rate / (1. / time_units);
    inputs->RT_H2_dissociation_rate =
        cooling->RT_H2_dissociation_rate / (1. / time_units);
  }

  /* The H2 column is N_H2 = 2 n_H2 l, so l = path / 2 gives the column
     through the chosen path of H2_self_shielding_path_in_kernel_radii kernel
     support radii H = kernel_gamma h. With comoving_coordinates = 0, a
     physical length in length_units is expected. */
  if (cooling->chemistry_data.H2_self_shielding == 2)
    inputs->H2_self_shielding_length =
        0.5f * cooling->H2_self_shielding_path_in_kernel_radii * kernel_gamma *
        p->h * cosmo->a;

  /* Self shielding by a density threshold: switch the UV background off */
  if (cooling->with_uv_background && cooling->self_shielding_method < 0 &&
      hydro_get_physical_density(p, cosmo) > cooling->self_shielding_threshold)
    inputs->uvb_off = 1;
}

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data, or keep the previous state if geajay fails.
 *
 * A failed solve (e.g. more than max_substeps sub-cycle iterations) is
 * retried subcycle_on_failure times by geajay, retry i with 2^i sub-steps of
 * dt/2^i. If every retry fails, the species are not changed and this step's
 * cooling and heating are skipped, not deferred: the next step integrates
 * only its own dt. The failure is counted. A non-positive energy after the
 * minimal-energy floor (a legitimate state, not a geajay failure) takes the
 * same kept-previous-state path without ever calling geajay.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current simulation time.
 * @param kept_previous_state (return) 1 if the particle kept its previous
 * state (every solve failed, or the solve was skipped), 0 otherwise.
 *
 * @return The new internal energy, or the energy passed to geajay if
 * kept_previous_state is 1.
 */
static gr_float cooling_new_energy_or_keep_previous(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double dt, double dt_therm, double time,
    int *kept_previous_state) {

  *kept_previous_state = 0;

  /* A non-finite stored value, or a non-positive stored density, cannot be
     integrated: geajay would freeze the cell at its energy floor. The
     stored energy is not required positive here: the minimal-energy floor
     just below can legitimately clamp it to exactly 0 (the default
     hydro_props_default_min_temp is 0 K), which is not itself an error. */
  const float rho_stored = hydro_get_comoving_density(p);
  const float u_stored = hydro_get_comoving_internal_energy(p, xp);
  const float u_dt_stored = hydro_get_comoving_internal_energy_dt(p);
  if (!cooling_stored_value_is_valid(&rho_stored, /*positive=*/1) ||
      !cooling_stored_value_is_valid(&u_stored, /*positive=*/0) ||
      !cooling_stored_value_is_valid(&u_dt_stored, /*positive=*/0))
    error(
        "Particle %lld has an invalid state before the geajay solve: "
        "comoving density = %e, comoving internal energy = %e, its time "
        "derivative = %e.",
        p->id, rho_stored, u_stored, u_dt_stored);

  /* general particle data */
  const double density = cooling_get_physical_density(p, cosmo, cooling);
  gr_float energy = hydro_get_physical_internal_energy(p, xp, cosmo) +
                    dt_therm * hydro_get_physical_internal_energy_dt(p, cosmo);
  energy = max(energy, hydro_props->minimal_internal_energy);

  /* The stored density already passed the finite/positive guard above, so a
     non-positive physical density here is a real bug (e.g. a bad comoving
     conversion): stop the run. */
  if (density <= 0.)
    error(
        "Particle %lld has a non-positive physical density before the "
        "geajay solve: physical density = %e.",
        p->id, density);

  /* A non-positive energy after the minimal-energy floor is not a bug:
     geajay cannot integrate from it, so skip the solve and keep the
     particle's previous state instead of calling error(). */
  if (energy <= 0.) {
    atomic_inc(&cooling_geajay_failed_solves);
    *kept_previous_state = 1;
    return energy;
  }

  geajay_particle_inputs inputs;
  cooling_fill_inputs(cooling, p, cosmo, density, energy, &inputs);

  /* The particle keeps its previous state on a final failure */
  inputs.store_on_failure = 0;

  geajay_particle_output output = {.size = sizeof(output)};
  geajay_particle_status status = {.size = sizeof(status)};

  const int solved =
      geajay_particle_solve(cooling->context, &cooling->layout, xp, &inputs, dt,
                            /*n_substeps=*/1, &output, &status) != GR_FAIL;

  if (status.reason == GEAJAY_PARTICLE_BAD_ARGUMENT)
    error("Particle %lld: geajay refused the arguments of the solve.", p->id);

  if (!solved) {
    atomic_inc(&cooling_geajay_failed_solves);
    *kept_previous_state = 1;
    return energy;
  }

  if (status.retries > 0) atomic_inc(&cooling_geajay_recovered_solves);

  return output.internal_energy;
}

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data
 *
 * See #cooling_new_energy_or_keep_previous for the handling of a failed
 * geajay solve.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current simulation time.
 *
 * @return The new internal energy.
 */
gr_float cooling_new_energy(const struct phys_const *phys_const,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct hydro_props *hydro_props,
                            const struct pressure_floor_props *pressure_floor,
                            const struct cooling_function_data *cooling,
                            struct part *p, struct xpart *xp, double dt,
                            double dt_therm, double time) {
  int kept_previous_state;
  return cooling_new_energy_or_keep_previous(
      phys_const, us, cosmo, hydro_props, pressure_floor, cooling, p, xp, dt,
      dt_therm, time, &kept_previous_state);
}

/**
 * @brief A quiet NaN made from its bits, for a value that is only stored.
 *
 * SWIFT builds with -ffinite-math-only, under which a NaN constant is
 * undefined: the volatile keeps the compiler from seeing one.
 *
 * @return A quiet NaN.
 */
static double cooling_stored_nan(void) {
  volatile uint64_t bits = 0x7ff8000000000000ULL;
  const uint64_t value = bits;
  double nan_value;
  memcpy(&nan_value, &value, sizeof(nan_value));
  return nan_value;
}

/**
 * @brief Compute the cooling time of a particle at its stored state.
 *
 * Signed as in Grackle: negative for net cooling, +inf at a net rate of 0.
 * geajay reads the particle and writes nothing to it.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 *
 * @return The physical cooling time in internal units, NaN if geajay refuses
 * the particle.
 */
double cooling_time(const struct phys_const *phys_const,
                    const struct unit_system *us,
                    const struct hydro_props *hydro_props,
                    const struct cosmology *cosmo,
                    const struct cooling_function_data *cooling,
                    const struct part *p, const struct xpart *xp) {

  /* general particle data */
  const double density = cooling_get_physical_density(p, cosmo, cooling);
  gr_float energy = hydro_get_physical_internal_energy(p, xp, cosmo);
  energy = max(energy, hydro_props->minimal_internal_energy);

  geajay_particle_inputs inputs;
  cooling_fill_inputs(cooling, p, cosmo, density, energy, &inputs);

  double cooling_time;
  if (geajay_particle_cooling_time(cooling->context, &cooling->layout, xp,
                                   &inputs, &cooling_time) != GR_SUCCESS)
    return cooling_stored_nan();

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
void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_props,
                       const struct entropy_floor_properties *floor_props,
                       const struct pressure_floor_props *pressure_floor,
                       const struct cooling_function_data *cooling,
                       struct part *p, struct xpart *xp, const double dt,
                       const double dt_therm, const double time) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Debug/test-only: skip the per-particle energy update entirely,
     leaving hydro's own internal-energy derivative untouched. See
     cooling_function_data.disable_cooling_for_debugging's own doxygen. */
  if (cooling->disable_cooling_for_debugging) return;

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Energy after the adiabatic cooling */
  float u_ad_before =
      u_old + dt_therm * hydro_get_physical_internal_energy_dt(p, cosmo);

  double u_CMB_agora = 0.0;
  if (cooling->agora_cmb_temperature_floor) {
    u_CMB_agora = cooling_agora_cmb_floor_internal_energy(
        phys_const, us, cosmo, hydro_props, cooling, p, xp);

    if (u_ad_before < u_CMB_agora) {
      u_ad_before = u_CMB_agora;
      const float du_dt = (u_ad_before - u_old) / dt_therm;
      hydro_set_physical_internal_energy_dt(p, cosmo, du_dt);
    }
  }

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
    int kept_previous_state;
    u_new = cooling_new_energy_or_keep_previous(
        phys_const, us, cosmo, hydro_props, pressure_floor, cooling, p, xp, dt,
        dt_therm, time, &kept_previous_state);

    /* No cooling this step: exactly the adiabatic energy */
    if (kept_previous_state) u_new = u_ad_before;

    /* The solve evolves the ion/molecule fractions that set mu;
       recompute the floor from the post-solve composition. */
    if (cooling->agora_cmb_temperature_floor) {
      u_CMB_agora = cooling_agora_cmb_floor_internal_energy(
          phys_const, us, cosmo, hydro_props, cooling, p, xp);
    }
  }

  /* Get the change in internal energy due to hydro forces */
  float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* We now need to check that we are not going to go below any of the limits.
     Any energy this adds is a floor, not radiative cooling/heating: fold it
     into hydro_du_dt (excluded from radiated_energy below) instead of
     cool_du_dt, matching the pre-solve floor clamp above, so both book the
     floor the same way regardless of when it triggers. */
  const gr_float u_new_solved = u_new;
  u_new = max3(u_new, u_minimal, u_CMB_agora);
  const float floor_injection = u_new - u_new_solved;
  if (floor_injection > 0.f) hydro_du_dt += floor_injection / dt_therm;

  /* Calculate the cooling rate */
  float cool_du_dt = (u_new_solved - u_ad_before) / dt_therm;
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
float cooling_get_temperature(const struct phys_const *phys_const,
                              const struct hydro_props *hydro_props,
                              const struct unit_system *us,
                              const struct cosmology *cosmo,
                              const struct cooling_function_data *cooling,
                              const struct part *p, const struct xpart *xp) {
  const double u = hydro_get_drifted_physical_internal_energy(p, cosmo);

#if COOLING_GEAJAY_MODE == 0
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;
  const double mu_neutral = hydro_props->mu_neutral;
  const double mu_ionised = hydro_props->mu_ionised;
  const double T_transition = hydro_props->hydrogen_ionization_temperature;
  const double T_over_mu = hydro_gamma_minus_one * u * m_H / k_B;

  /* Are we above or below the HII -> HI transition? */
  if (T_over_mu > (T_transition + 1.) / mu_ionised)
    return T_over_mu * mu_ionised;
  else if (T_over_mu < (T_transition - 1.) / mu_neutral)
    return T_over_mu * mu_neutral;
  else
    return T_transition;

#elif COOLING_GEAJAY_MODE == 1
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);
  return cooling_temperature_from_internal_energy(u, mu, k_B, m_H);

#else /* COOLING_GEAJAY_MODE >= 2 */
  return cooling_get_temperature_h2_gamma_corrected(phys_const, cosmo, p, xp,
                                                    u);
#endif
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
double cooling_get_ycompton(const struct phys_const *phys_const,
                            const struct hydro_props *hydro_props,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            const struct part *p, const struct xpart *xp) {

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
float cooling_timestep(const struct cooling_function_data *cooling,
                       const struct phys_const *phys_const,
                       const struct cosmology *cosmo,
                       const struct unit_system *us,
                       const struct hydro_props *hydro_props,
                       const struct part *p, const struct xpart *xp) {

  return FLT_MAX;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part *p, struct xpart *xp, double n) {

  error("Loic: to be implemented");
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param phys_const The #phys_const.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_units(const struct unit_system *us,
                        const struct phys_const *phys_const,
                        struct cooling_function_data *cooling) {

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
 * @brief Describe where the species live in the #xpart.
 *
 * The species are stored as mass fractions of the particle density, as floats.
 *
 * @param layout (return) The layout.
 */
static void cooling_init_species_layout(geajay_species_layout *layout) {

  geajay_species_layout_init(layout, GEAJAY_FLOAT, GEAJAY_MASS_FRACTION);

#if COOLING_GEAJAY_MODE >= 1
  layout->offset[GEAJAY_HI] = offsetof(struct xpart, cooling_data.HI_frac);
  layout->offset[GEAJAY_HII] = offsetof(struct xpart, cooling_data.HII_frac);
  layout->offset[GEAJAY_HeI] = offsetof(struct xpart, cooling_data.HeI_frac);
  layout->offset[GEAJAY_HeII] = offsetof(struct xpart, cooling_data.HeII_frac);
  layout->offset[GEAJAY_HeIII] =
      offsetof(struct xpart, cooling_data.HeIII_frac);
  layout->offset[GEAJAY_E] = offsetof(struct xpart, cooling_data.e_frac);
#endif

#if COOLING_GEAJAY_MODE >= 2
  layout->offset[GEAJAY_HM] = offsetof(struct xpart, cooling_data.HM_frac);
  layout->offset[GEAJAY_H2I] = offsetof(struct xpart, cooling_data.H2I_frac);
  layout->offset[GEAJAY_H2II] = offsetof(struct xpart, cooling_data.H2II_frac);
#endif

#if COOLING_GEAJAY_MODE >= 3
  layout->offset[GEAJAY_DI] = offsetof(struct xpart, cooling_data.DI_frac);
  layout->offset[GEAJAY_DII] = offsetof(struct xpart, cooling_data.DII_frac);
  layout->offset[GEAJAY_HDI] = offsetof(struct xpart, cooling_data.HDI_frac);
#endif
}

/**
 * @brief Initialises geajay.
 *
 * @param cooling The cooling properties to initialize
 */
void cooling_init_geajay(struct cooling_function_data *cooling) {

  chemistry_data *chemistry = &cooling->chemistry_data;

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
  chemistry->H2_photodissociation_heating =
      cooling->H2_photodissociation_heating;
  chemistry->grackle_data_file = cooling->cloudy_table;

  if (cooling->local_dust_to_gas_ratio > 0)
    chemistry->local_dust_to_gas_ratio = cooling->local_dust_to_gas_ratio;

  /* radiative transfer */
#if COOLING_GEAJAY_MODE == 0
  if (cooling->use_radiative_transfer)
    error(
        "The parameter GeajayCooling:use_radiative_transfer_fields cannot be "
        "set "
        "to 1 in geajay mode 0 !");
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

  chemistry->H2_self_shielding = cooling->H2_self_shielding;

  /* Sub-cycle limit. Exceeding it is a failure that geajay retries, instead of
     a state integrated over only part of the time-step. */
  chemistry->max_iterations = cooling->max_step;
  chemistry->exit_after_iterations_exceeded = 1;

  if (local_initialize_chemistry_data(&cooling->chemistry_data,
                                      &cooling->chemistry_rates,
                                      &cooling->units) == 0) {
    error("Error in initialize_chemistry_data");
  }

  /* geajay resolves a negative hydrogen mass fraction to its own value */
  cooling->HydrogenFractionByMass =
      cooling->chemistry_data.HydrogenFractionByMass;

  /* Here, so that a restart pads the reloaded tables again */
  if (geajay_set_table_density_edge(&cooling->chemistry_rates,
                                    cooling->table_density_edge) == GR_FAIL)
    error("geajay refused GeajayCooling:table_density_edge: %s",
          geajay_last_error_message());

  if (geajay_table_density_range(&cooling->chemistry_rates,
                                 &cooling->table_density_min_Hpcm3,
                                 &cooling->table_density_max_Hpcm3) == GR_FAIL)
    error("geajay could not give the density range of the Cloudy table.");
  cooling->table_edge_density =
      cooling->table_density_max_Hpcm3 * cooling->proton_mass_cgs /
      (cooling->HydrogenFractionByMass * cooling->units.density_units);

  geajay_context_error context_error;
  cooling->context =
      geajay_context_create(&cooling->chemistry_data, &cooling->chemistry_rates,
                            &cooling->units, &context_error);
  if (cooling->context == NULL)
    error("geajay refused the chemistry parameters (error %d).",
          (int)context_error);

  if (cooling->redshift != -1 &&
      geajay_context_set_redshift(cooling->context, cooling->redshift) ==
          GR_FAIL)
    error("Invalid redshift %g for geajay.", cooling->redshift);

  geajay_context_set_subcycle_on_failure(cooling->context,
                                         cooling->subcycle_on_failure);

  cooling_init_species_layout(&cooling->layout);
  if (geajay_context_check_layout(cooling->context, &cooling->layout,
                                  sizeof(struct xpart)) == GR_FAIL)
    error("The species layout of the xpart does not fit the geajay context.");
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
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling) {

  /* read parameters */
  cooling_read_parameters(parameter_file, cooling, phys_const, us);

  /* Set up the units system. */
  cooling_init_units(us, phys_const, cooling);

  /* Set up geajay */
  cooling_init_geajay(cooling);
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data *cooling) {
  geajay_context_free(cooling->context);
  cooling->context = NULL;

  /* Clean up the rate data. This is a call to a geajay function */
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
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream) {
  restart_write_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
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
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo) {
  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  /* Set up geajay: the context pointer read above is stale */
  cooling_init_geajay(cooling);
}

/**
 * @brief Get the physical density of the #part passed to geajay.
 *
 * The first particle of a run above the density range of the Cloudy table
 * prints a warning. Its n_H is X_H rho / m_p, close to the n_H of geajay,
 * which sums the hydrogen species.
 *
 * Note: This function is called in cooling_time() and cooling_new_energy().
 *
 * @param p #part data.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 */
double cooling_get_physical_density(
    const struct part *p, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling) {

  const double part_density = hydro_get_physical_density(p, cosmo);

  if (part_density > cooling->table_edge_density &&
      !cooling_geajay_table_edge_warned &&
      atomic_cas(&cooling_geajay_table_edge_warned, 0, 1) == 0)
    warning(
        "Gas particle %lld has n_H = %e cm^-3, above the Cloudy table range "
        "(%e cm^-3): the table lookup %s. Printed once per rank, and again "
        "after a restart.",
        p->id,
        cooling->table_density_max_Hpcm3 * part_density /
            cooling->table_edge_density,
        cooling->table_density_max_Hpcm3,
        cooling->table_density_edge == GEAJAY_TABLE_EDGE_HOLD
            ? "holds the edge value"
            : "extrapolates");

  return part_density;
}
