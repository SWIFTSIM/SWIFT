/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
 * @file src/cooling/QLA_EAGLE/cooling.c
 * @brief QLA_EAGLE cooling functions
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "active.h"
#include "chemistry.h"
#include "cooling.h"
#include "cooling_rates.h"
#include "cooling_struct.h"
#include "cooling_tables.h"
#include "entropy_floor.h"
#include "error.h"
#include "exp10.h"
#include "hydro.h"
#include "interpolate.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/* Maximum number of iterations for bisection scheme */
static const int bisection_max_iterations = 150;

/* Tolerances for termination criteria. */
static const float explicit_tolerance = 0.05;
static const float bisection_tolerance = 1.0e-6;
static const double bracket_factor = 1.5;

/**
 * @brief Find the index of the current redshift along the redshift dimension
 * of the cooling tables.
 *
 * Since the redshift table is not evenly spaced, compare z with each
 * table value in decreasing order starting with the previous redshift index
 *
 * The returned difference is expressed in units of the table separation. This
 * means dx = (x - table[i]) / (table[i+1] - table[i]). It is always between
 * 0 and 1.
 *
 * @param z Redshift we are searching for.
 * @param z_index (return) Index of the redshift in the table.
 * @param dz (return) Difference in redshift between z and table[z_index].
 * @param cooling #cooling_function_data structure containing redshift table.
 */
__attribute__((always_inline)) INLINE void get_redshift_index(
    const float z, int *z_index, float *dz,
    struct cooling_function_data *restrict cooling) {

  /* Before the earliest redshift or before hydrogen reionization, flag for
   * collisional cooling */
  if (z > cooling->H_reion_z) {
    *z_index = qla_eagle_cooling_N_redshifts;
    *dz = 0.0;
  }

  /* From reionization use the cooling tables */
  else if (z > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1] &&
           z <= cooling->H_reion_z) {
    *z_index = qla_eagle_cooling_N_redshifts + 1;
    *dz = 0.0;
  }

  /* At the end, just use the last value */
  else if (z <= cooling->Redshifts[0]) {
    *z_index = 0;
    *dz = 0.0;
  }

  /* Normal case: search... */
  else {

    /* start at the previous index and search */
    for (int i = cooling->previous_z_index; i >= 0; i--) {
      if (z > cooling->Redshifts[i]) {

        *z_index = i;
        cooling->previous_z_index = i;

        *dz = (z - cooling->Redshifts[i]) /
              (cooling->Redshifts[i + 1] - cooling->Redshifts[i]);
        break;
      }
    }
  }
}

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * Also calls the additional H reionisation energy injection if need be.
 *
 * @param cosmo The current cosmological model.
 * @param pressure_floor The properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The space data, including a pointer to array of particles
 * @param time The current system time
 */
void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time) {

  /* Current redshift */
  const float redshift = cosmo->z;

  /* What is the current table index along the redshift axis? */
  int z_index = -1;
  float dz = 0.f;
  get_redshift_index(redshift, &z_index, &dz, cooling);
  cooling->dz = dz;

  /* Extra energy for reionization? */
  if (!cooling->H_reion_done) {

    /* Does this timestep straddle Hydrogen reionization? If so, we need to
     * input extra heat */
    if (cosmo->z <= cooling->H_reion_z && cosmo->z_old > cooling->H_reion_z) {

      if (s == NULL) error("Trying to do H reionization on an empty space!");

      /* Inject energy to all particles */
      cooling_Hydrogen_reionization(cooling, cosmo, pressure_floor, s);

      /* Flag that reionization happened */
      cooling->H_reion_done = 1;
    }
  }

  /* Do we already have the correct tables loaded? */
  if (cooling->z_index == z_index) return;

  /* Which table should we load ? */
  if (z_index >= qla_eagle_cooling_N_redshifts) {

    if (z_index == qla_eagle_cooling_N_redshifts + 1) {

      /* Bewtween re-ionization and first table */
      get_redshift_invariant_table(cooling, /* photodis=*/0);

    } else {

      /* Above re-ionization */
      get_redshift_invariant_table(cooling, /* photodis=*/1);
    }

  } else {

    /* Normal case: two tables bracketing the current z */
    const int low_z_index = z_index;
    const int high_z_index = z_index + 1;

    get_cooling_table(cooling, low_z_index, high_z_index);
  }

  /* Store the currently loaded index */
  cooling->z_index = z_index;
}

/**
 * @brief Bisection integration scheme
 *
 * @param u_ini_cgs Internal energy at beginning of hydro step in CGS.
 * @param n_H_cgs Hydrogen number density in CGS.
 * @param redshift Current redshift.
 * @param n_H_index Particle hydrogen number density index.
 * @param d_n_H Particle hydrogen number density offset.
 * @param He_index Particle helium fraction index.
 * @param d_He Particle helium fraction offset.
 * @param Lambda_He_reion_cgs Cooling rate coming from He reionization.
 * @param ratefact_cgs Multiplication factor to get a cooling rate.
 * @param cooling #cooling_function_data structure.
 * @param abundance_ratio Array of ratios of metal abundance to solar.
 * @param dt_cgs timestep in CGS.
 * @param ID ID of the particle (for debugging).
 */
INLINE static double bisection_iter(
    const double u_ini_cgs, const double n_H_cgs, const double redshift,
    const int n_H_index, const float d_n_H, const int He_index,
    const float d_He, const double Lambda_He_reion_cgs,
    const double ratefact_cgs,
    const struct cooling_function_data *restrict cooling,
    const float abundance_ratio[qla_eagle_cooling_N_abundances],
    const double dt_cgs, const long long ID) {

  /* Bracketing */
  double u_lower_cgs = u_ini_cgs;
  double u_upper_cgs = u_ini_cgs;

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  double LambdaNet_cgs =
      Lambda_He_reion_cgs + qla_eagle_cooling_rate(log10(u_ini_cgs), redshift,
                                                   n_H_cgs, abundance_ratio,
                                                   n_H_index, d_n_H, He_index,
                                                   d_He, cooling);

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  if (LambdaNet_cgs < 0) {

    /* we're cooling! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs = Lambda_He_reion_cgs +
                    qla_eagle_cooling_rate(log10(u_lower_cgs), redshift,
                                           n_H_cgs, abundance_ratio, n_H_index,
                                           d_n_H, He_index, d_He, cooling);

    int i = 0;
    while (u_lower_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs >
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs /= bracket_factor;
      u_upper_cgs /= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          qla_eagle_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs,
                                 abundance_ratio, n_H_index, d_n_H, He_index,
                                 d_He, cooling);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling, u_ini_cgs %.5e n_H_cgs %.5e",
          ID, u_ini_cgs, n_H_cgs);
    }
  } else {

    /* we are heating! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs = Lambda_He_reion_cgs +
                    qla_eagle_cooling_rate(log10(u_upper_cgs), redshift,
                                           n_H_cgs, abundance_ratio, n_H_index,
                                           d_n_H, He_index, d_He, cooling);

    int i = 0;
    while (u_upper_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs <
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs *= bracket_factor;
      u_upper_cgs *= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          qla_eagle_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs,
                                 abundance_ratio, n_H_index, d_n_H, He_index,
                                 d_He, cooling);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "heating, u_ini_cgs %.5e n_H_cgs %.5e",
          ID, u_ini_cgs, n_H_cgs);
    }
  }

  /********************************************/
  /* We now have an upper and lower bound.    */
  /* Let's iterate by reducing the bracketing */
  /********************************************/

  /* bisection iteration */
  int i = 0;
  double u_next_cgs;

  do {

    /* New guess */
    u_next_cgs = 0.5 * (u_lower_cgs + u_upper_cgs);

    /* New rate */
    LambdaNet_cgs = Lambda_He_reion_cgs +
                    qla_eagle_cooling_rate(log10(u_next_cgs), redshift, n_H_cgs,
                                           abundance_ratio, n_H_index, d_n_H,
                                           He_index, d_He, cooling);
#ifdef SWIFT_DEBUG_CHECKS
    if (u_next_cgs <= 0)
      error(
          "Got negative energy! u_next_cgs=%.5e u_upper=%.5e u_lower=%.5e "
          "Lambda=%.5e",
          u_next_cgs, u_upper_cgs, u_lower_cgs, LambdaNet_cgs);
#endif

    /* Where do we go next? */
    if (u_next_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs > 0.0) {
      u_upper_cgs = u_next_cgs;
    } else {
      u_lower_cgs = u_next_cgs;
    }

    i++;
  } while (fabs(u_upper_cgs - u_lower_cgs) / u_next_cgs > bisection_tolerance &&
           i < bisection_max_iterations);

  if (i >= bisection_max_iterations)
    error("Particle id %llu failed to converge", ID);

  return u_upper_cgs;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * We want to compute u_new such that u_new = u_old + dt * du/dt(u_new, X),
 * where X stands for the metallicity, density and redshift. These are
 * kept constant.
 *
 * We first compute du/dt(u_old). If dt * du/dt(u_old) is small enough, we
 * use an explicit integration and use this as our solution.
 *
 * Otherwise, we try to find a solution to the implicit time-integration
 * problem. This leads to the root-finding problem:
 *
 * f(u_new) = u_new - u_old - dt * du/dt(u_new, X) = 0
 *
 * We first try a few Newton-Raphson iteration if it does not converge, we
 * revert to a bisection scheme.
 *
 * This is done by first bracketing the solution and then iterating
 * towards the solution by reducing the window down to a certain tolerance.
 * Note there is always at least one solution since
 * f(+inf) is < 0 and f(-inf) is > 0.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_properties the hydro_props struct
 * @param floor_props Properties of the entropy floor.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The cooling time-step of this particle.
 * @param dt_therm The hydro time-step of this particle.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct pressure_floor_props *pressure_floor,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm,
                       const double time) {

  /* No cooling happens over zero time */
  if (dt == 0.) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--cooling runtime flag?");
#endif

  /* Get internal energy at the last kick step */
  const float u_start = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Get the change in internal energy due to hydro forces */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Get internal energy at the end of the step (assuming dt does not
   * increase) */
  double u_0 = (u_start + hydro_du_dt * dt_therm);

  /* Check for minimal energy */
  u_0 = max(u_0, hydro_properties->minimal_internal_energy);

  /* Convert to CGS units */
  const double u_0_cgs = u_0 * cooling->internal_energy_to_cgs;
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Change in redshift over the course of this time-step
     (See cosmology theory document for the derivation) */
  const double delta_redshift = -dt * cosmo->H * cosmo->a_inv;

  /* Get this particle's abundance ratios compared to solar
   * Note that we need to add S and Ca that are in the tables but not tracked
   * by the particles themselves.
   * The order is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe] */
  float abundance_ratio[qla_eagle_cooling_N_abundances];
  abundance_ratio_to_solar(p, cooling, phys_const, abundance_ratio);

  /* Get the Hydrogen and Helium mass fractions */
  const float XH = 1. - phys_const->const_primordial_He_fraction;
  const float XHe = phys_const->const_primordial_He_fraction;

  /* Get the Helium mass fraction. Note that this is He / (H + He), i.e. a
   * metal-free Helium mass fraction as per the Wiersma+08 definition */
  const float HeFrac = XHe / (XH + XHe);

  /* convert Hydrogen mass fraction into physical Hydrogen number density */
  const double n_H =
      hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* ratefact = n_H * n_H / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  const double ratefact_cgs = n_H_cgs * (XH * cooling->inv_proton_mass_cgs);

  /* compute hydrogen number density and helium fraction table indices and
   * offsets (These are fixed for any value of u, so no need to recompute them)
   */
  int He_index, n_H_index;
  float d_He, d_n_H;
  get_index_1d(cooling->HeFrac, qla_eagle_cooling_N_He_frac, HeFrac, &He_index,
               &d_He);
  get_index_1d(cooling->nH, qla_eagle_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Start by computing the cooling (heating actually) rate from Helium
     re-ionization as this needs to be added on no matter what */

  /* Get helium and hydrogen reheating term */
  const double Helium_reion_heat_cgs = qla_eagle_helium_reionization_extraheat(
      cosmo->z, delta_redshift, cooling);

  /* Convert this into a rate */
  const double Lambda_He_reion_cgs =
      Helium_reion_heat_cgs / (dt_cgs * ratefact_cgs);

  /* Let's compute the internal energy at the end of the step */
  /* Initialise to the initial energy to appease compiler; this will never not
     be overwritten. */
  double u_final_cgs = u_0_cgs;

  /* First try an explicit integration (note we ignore the derivative) */
  const double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      qla_eagle_cooling_rate(log10(u_0_cgs), cosmo->z, n_H_cgs, abundance_ratio,
                             n_H_index, d_n_H, He_index, d_He, cooling);

  /* if cooling rate is small, take the explicit solution */
  if (fabs(ratefact_cgs * LambdaNet_cgs * dt_cgs) <
      explicit_tolerance * u_0_cgs) {

    u_final_cgs = u_0_cgs + ratefact_cgs * LambdaNet_cgs * dt_cgs;

  } else {

    /* Otherwise, go the bisection route. */
    u_final_cgs =
        bisection_iter(u_0_cgs, n_H_cgs, cosmo->z, n_H_index, d_n_H, He_index,
                       d_He, Lambda_He_reion_cgs, ratefact_cgs, cooling,
                       abundance_ratio, dt_cgs, p->id);
  }

  /* Convert back to internal units */
  double u_final = u_final_cgs * cooling->internal_energy_from_cgs;

  /* We now need to check that we are not going to go below any of the limits */

  /* Absolute minimum */
  const double u_minimal = hydro_properties->minimal_internal_energy;
  u_final = max(u_final, u_minimal);

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho_physical = hydro_get_physical_density(p, cosmo);
  const double u_floor =
      gas_internal_energy_from_entropy(rho_physical, A_floor);
  u_final = max(u_final, u_floor);

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u = u_final - max(u_start, u_floor);

  /* Turn this into a rate of change (including cosmology term) */
  const float cooling_du_dt = delta_u / dt_therm;

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);
}

/**
 * @brief Computes the cooling time-step.
 *
 * The time-step is not set by the properties of cooling.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const #phys_const data struct.
 * @param us The internal system of units.
 * @param cosmo #cosmology struct.
 * @param hydro_props the properties of the hydro scheme.
 * @param p #part data.
 * @param xp extended particle data.
 */
__attribute__((always_inline)) INLINE float cooling_timestep(
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props, const struct part *restrict p,
    const struct xpart *restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const #phys_const data structure.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
__attribute__((always_inline)) INLINE void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp) {}

/**
 * @brief Perform additional init on the cooling properties of the
 * (x-)particles that requires the density to be known.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE void cooling_post_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp) {}

/**
 * @brief Compute the temperature based on gas properties.
 *
 * Dummy function!
 *
 * @param phys_const #phys_const data structure.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param rho_phys Density of the gas in internal physical units.
 * @param logZZsol Logarithm base 10 of the gas' metallicity in units of solar
 * metallicity.
 * @param XH The Hydrogen abundance of the gas.
 * @param u_phys Internal energy of the gas in internal physical units.
 * @param HII_region Is this patch of gas an HII region? (Not implemented in
 * QLA_EAGLE)
 */
float cooling_get_temperature_from_gas(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const float rho_phys,
    const float XH, const float logZZsol, const float u_phys,
    const int HII_region) {

  error("Do not call this function");
  return -1.f;
}

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * We use the Temperature table of the Wiersma+08 set. This computes the
 * equilibirum temperature of a gas for a given redshift, Hydrogen density,
 * internal energy per unit mass and Helium fraction.
 *
 * The temperature returned is consistent with the cooling rates.
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
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--temperature runtime flag?");
#endif

  /* Get physical internal energy */
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const double u_cgs = u * cooling->internal_energy_to_cgs;

  /* Get the Hydrogen and Helium mass fractions */
  const float XH = 1. - phys_const->const_primordial_He_fraction;
  const float XHe = phys_const->const_primordial_He_fraction;

  /* Get the Helium mass fraction. Note that this is He / (H + He), i.e. a
   * metal-free Helium mass fraction as per the Wiersma+08 definition */
  const float HeFrac = XHe / (XH + XHe);

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float rho = hydro_get_physical_density(p, cosmo);
  const double n_H = rho * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* compute hydrogen number density and helium fraction table indices and
   * offsets */
  int He_index, n_H_index;
  float d_He, d_n_H;
  get_index_1d(cooling->HeFrac, qla_eagle_cooling_N_He_frac, HeFrac, &He_index,
               &d_He);
  get_index_1d(cooling->nH, qla_eagle_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Compute the log10 of the temperature by interpolating the table */
  const double log_10_T = qla_eagle_convert_u_to_temp(
      log10(u_cgs), cosmo->z, n_H_index, He_index, d_n_H, d_He, cooling);

  /* Undo the log! */
  return exp10(log_10_T);
}

/**
 * @brief Compute the electron number density of a #part based on the cooling
 * function.
 *
 * Returns -1 in this model.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_electron_density(const struct phys_const *phys_const,
                                   const struct hydro_props *hydro_props,
                                   const struct unit_system *us,
                                   const struct cosmology *cosmo,
                                   const struct cooling_function_data *cooling,
                                   const struct part *p,
                                   const struct xpart *xp) {

  return -1.;
}

/**
 * @brief Compute the electron pressure of a #part based on the cooling
 * function.
 *
 * There are no subgrid properties in this model, we return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
double cooling_get_electron_pressure(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  return 0.;
}

/**
 * @brief Compute the y-Compton contribution of a #part based on the cooling
 * function.
 *
 * Returns -1 in this model.
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

  return -1.f;
}

/**
 * @brief Compute the HI fraction of a #part based on the cooling function.
 *
 * There are no subgrid properties in this model, we return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_particle_subgrid_HI_fraction(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  return 0.f;
}

/**
 * @brief Compute the HI fraction of a #part based on the cooling function.
 *
 * There are no subgrid properties in this model, we return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_particle_subgrid_HII_fraction(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  return 0.f;
}

/**
 * @brief Compute the H2 fraction of a #part based on the cooling function.
 *
 * There are no subgrid properties in this model, we return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_particle_subgrid_H2_fraction(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  return 0.f;
}

double compute_subgrid_property(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const float rho_phys, const float logZZsol,
    const float XH, const float P_phys, const float log10_T,
    const float log10_T_EOS_max, const int HII_region,
    const float *abundance_ratio, const double log_u_cgs,
    const enum cooling_subgrid_properties isub) {

  error("Do not call this function");
  return -1.f;
}

/**
 * @brief Compute the physical subgrid density of the gas.
 *
 * There is no subgrid density in this model so we just
 * return the regular density.
 *
 * Note that we return the density in physical coordinates.
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param rho_phys Density of the gas in internal physical units.
 * @param logZZsol Logarithm base 10 of the gas' metallicity in units of solar
 * metallicity.
 * @param XH The Hydrogen abundance of the gas.
 * @param P_phys Pressure of the gas in internal physical units.
 * @param log10_T The logarithm base 10 of the temperature of the gas.
 * @param log10_T_EOS_max The logarithm base 10 of the maximal temperature
 * to be considered on the EOS at the density of the gas.
 */
double compute_subgrid_density(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const float rho_phys, const float logZZsol,
    const float XH, const float P_phys, const float log10_T,
    const float log10_T_EOS_max) {

  return rho_phys;
}

/**
 * @brief Returns the subgrid temperature of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

/**
 * @brief Returns the subgrid density of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
float cooling_get_subgrid_density(const struct part *p,
                                  const struct xpart *xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp #xpart data struct
 */
__attribute__((always_inline)) INLINE float cooling_get_radiated_energy(
    const struct xpart *restrict xp) {

  return -1.f;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part *p, struct xpart *xp, double n) {}

/**
 * @brief Inject a fixed amount of energy to each particle in the simulation
 * to mimic Hydrogen reionization.
 *
 * @param cooling The properties of the cooling model.
 * @param cosmo The cosmological model.
 * @param pressure_floor Properties of the pressure floor.
 * @param s The #space containing the particles.
 */
void cooling_Hydrogen_reionization(
    const struct cooling_function_data *cooling, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, struct space *s) {

  struct part *parts = s->parts;
  struct xpart *xparts = s->xparts;

  /* Energy to inject in internal units */
  const float extra_heat =
      cooling->H_reion_heat_cgs * cooling->internal_energy_from_cgs;

  message("Applying extra energy for H reionization!");

  /* Loop through particles and set new heat */
  for (size_t i = 0; i < s->nr_parts; i++) {

    struct part *p = &parts[i];
    struct xpart *xp = &xparts[i];

    if (part_is_inhibited(p, s->e)) continue;

    const float old_u = hydro_get_physical_internal_energy(p, xp, cosmo);
    const float new_u = old_u + extra_heat;

    hydro_set_physical_internal_energy(p, xp, cosmo, new_u);
    hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, new_u);
  }
}

/**
 * @brief Initialises properties stored in the cooling_function_data struct
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const #phys_const data structure
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling #cooling_function_data struct to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling) {

  /* Read model parameters */

  /* Directory for cooling tables */
  parser_get_param_string(parameter_file, "QLACooling:dir_name",
                          cooling->cooling_table_path);

  /* Despite the names, the values of H_reion_heat_cgs and He_reion_heat_cgs
   * that are read in are actually in units of electron volts per proton mass.
   * We later convert to units just below */

  cooling->H_reion_done = 0;
  cooling->H_reion_z =
      parser_get_param_float(parameter_file, "QLACooling:H_reion_z");
  cooling->H_reion_heat_cgs =
      parser_get_param_float(parameter_file, "QLACooling:H_reion_eV_p_H");
  cooling->He_reion_z_centre =
      parser_get_param_float(parameter_file, "QLACooling:He_reion_z_centre");
  cooling->He_reion_z_sigma =
      parser_get_param_float(parameter_file, "QLACooling:He_reion_z_sigma");
  cooling->He_reion_heat_cgs =
      parser_get_param_float(parameter_file, "QLACooling:He_reion_eV_p_H");

  /* Convert H_reion_heat_cgs and He_reion_heat_cgs to cgs
   * (units used internally by the cooling routines). This is done by
   * multiplying by 'eV/m_H' in internal units, then converting to cgs units.
   * Note that the dimensions of these quantities are energy/mass = velocity^2
   */

  cooling->H_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  cooling->He_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Read in the list of redshifts */
  get_cooling_redshifts(cooling);

  /* Read in cooling table header */
  char fname[qla_eagle_table_path_name_length + 12];
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  read_cooling_header(fname, cooling);

  /* Allocate space for cooling tables */
  allocate_cooling_tables(cooling);

  /* Compute conversion factors */
  cooling->internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  cooling->internal_energy_from_cgs = 1. / cooling->internal_energy_to_cgs;
  cooling->number_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Store some constants in CGS units */
  const double proton_mass_cgs =
      phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->inv_proton_mass_cgs = 1. / proton_mass_cgs;
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Compute the coefficient at the front of the Compton cooling expression */
  const double radiation_constant =
      4. * phys_const->const_stefan_boltzmann / phys_const->const_speed_light_c;
  const double compton_coefficient =
      4. * radiation_constant * phys_const->const_thomson_cross_section *
      phys_const->const_boltzmann_k /
      (phys_const->const_electron_mass * phys_const->const_speed_light_c);
  const float dimension_coefficient[5] = {1, 2, -3, 0, -5};

  /* This should be ~1.0178085e-37 g cm^2 s^-3 K^-5 */
  const double compton_coefficient_cgs =
      compton_coefficient *
      units_general_cgs_conversion_factor(us, dimension_coefficient);

#ifdef SWIFT_DEBUG_CHECKS
  const double expected_compton_coefficient_cgs = 1.0178085e-37;
  if (fabs(compton_coefficient_cgs - expected_compton_coefficient_cgs) /
          expected_compton_coefficient_cgs >
      0.01)
    error("compton coefficient incorrect.");
#endif

  /* And now the Compton rate */
  cooling->compton_rate_cgs = compton_coefficient_cgs * cooling->T_CMB_0 *
                              cooling->T_CMB_0 * cooling->T_CMB_0 *
                              cooling->T_CMB_0;

  /* Set the redshift indices to invalid values */
  cooling->z_index = -10;

  /* set previous_z_index and to last value of redshift table*/
  cooling->previous_z_index = qla_eagle_cooling_N_redshifts - 2;
}

/**
 * @brief Restore cooling tables (if applicable) after
 * restart
 *
 * @param cooling the #cooling_function_data structure
 * @param cosmo #cosmology structure
 */
void cooling_restore_tables(struct cooling_function_data *cooling,
                            const struct cosmology *cosmo) {

  /* Read redshifts */
  get_cooling_redshifts(cooling);

  /* Read cooling header */
  char fname[qla_eagle_table_path_name_length + 12];
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  read_cooling_header(fname, cooling);

  /* Allocate memory for the tables */
  allocate_cooling_tables(cooling);

  /* Force a re-read of the cooling tables */
  cooling->z_index = -10;
  cooling->previous_z_index = qla_eagle_cooling_N_redshifts - 2;
  cooling_update(/*phys_const=*/NULL, cosmo, /*pfloor=*/NULL, cooling,
                 /*space=*/NULL, /*time=*/0);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling #cooling_function_data struct.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  message(
      "Cooling function is 'Quick Lyman-alpha (Wiersma+09 tables with "
      "primordial Z only)'.");
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * We simply free all the arrays.
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data *cooling) {

  /* Free the side arrays */
  swift_free("cooling", cooling->Redshifts);
  swift_free("cooling", cooling->nH);
  swift_free("cooling", cooling->Temp);
  swift_free("cooling", cooling->HeFrac);
  swift_free("cooling", cooling->Therm);
  swift_free("cooling", cooling->SolarAbundances);
  swift_free("cooling", cooling->SolarAbundances_inv);

  /* Free the tables */
  swift_free("cooling-tables", cooling->table.metal_heating);
  swift_free("cooling-tables", cooling->table.electron_abundance);
  swift_free("cooling-tables", cooling->table.temperature);
  swift_free("cooling-tables", cooling->table.H_plus_He_heating);
  swift_free("cooling-tables", cooling->table.H_plus_He_electron_abundance);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the cooling routines. Helps debugging. */
  struct cooling_function_data cooling_copy = *cooling;
  cooling_copy.Redshifts = NULL;
  cooling_copy.nH = NULL;
  cooling_copy.Temp = NULL;
  cooling_copy.Therm = NULL;
  cooling_copy.SolarAbundances = NULL;
  cooling_copy.SolarAbundances_inv = NULL;
  cooling_copy.table.metal_heating = NULL;
  cooling_copy.table.H_plus_He_heating = NULL;
  cooling_copy.table.H_plus_He_electron_abundance = NULL;
  cooling_copy.table.temperature = NULL;
  cooling_copy.table.electron_abundance = NULL;

  restart_write_blocks((void *)&cooling_copy,
                       sizeof(struct cooling_function_data), 1, stream,
                       "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the cooling tables by
 * re-reading them.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo) {
  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  cooling_restore_tables(cooling, cosmo);
}
