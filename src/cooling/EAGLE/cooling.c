/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 * @file src/cooling/EAGLE/cooling.c
 * @brief EAGLE cooling functions
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
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
#include "units.h"

/* Maximum number of iterations for newton
 * and bisection integration schemes */
static const int newton_max_iterations = 15;
static const int bisection_max_iterations = 150;

/* Tolerances for termination criteria. */
static const float explicit_tolerance = 0.05;
static const float newton_tolerance = 1.0e-4;
static const float bisection_tolerance = 1.0e-6;
static const float rounding_tolerance = 1.0e-4;
static const double bracket_factor = 1.5;              /* sqrt(1.1) */
static const double newton_log_u_guess_cgs = 12.30103; /* log10(2e12) */

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
    *z_index = eagle_cooling_N_redshifts;
    *dz = 0.0;
  }

  /* From reionization use the cooling tables */
  else if (z > cooling->Redshifts[eagle_cooling_N_redshifts - 1] &&
           z <= cooling->H_reion_z) {
    *z_index = eagle_cooling_N_redshifts + 1;
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
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 */
void cooling_update(const struct cosmology *cosmo,
                    struct cooling_function_data *cooling) {

  /* Current redshift */
  const float redshift = cosmo->z;

  /* What is the current table index along the redshift axis? */
  int z_index = -1;
  float dz = 0.f;
  get_redshift_index(redshift, &z_index, &dz, cooling);
  cooling->dz = dz;

  /* Do we already have the correct tables loaded? */
  if (cooling->z_index == z_index) return;

  /* Which table should we load ? */
  if (z_index >= eagle_cooling_N_redshifts) {

    if (z_index == eagle_cooling_N_redshifts + 1) {

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
 * @brief Newton Raphson integration scheme to calculate particle cooling over
 * timestep. This replaces bisection scheme used in EAGLE to minimize the
 * number of array accesses. Integration defaults to bisection scheme (see
 * function bisection_iter) if this function does not converge within a
 * specified number of steps
 *
 * @param logu_init Initial guess for log(internal energy)
 * @param u_ini Internal energy at beginning of hydro step
 * @param n_H_index Particle hydrogen number density index
 * @param d_n_H Particle hydrogen number density offset
 * @param He_index Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param He_reion_heat Heating due to helium reionization
 * (only depends on redshift, so passed as parameter)
 * @param p #part structure
 * @param cosmo #cosmology structure
 * @param cooling #cooling_function_data structure
 * @param phys_const #phys_const data structure
 * @param abundance_ratio Array of ratios of metal abundance to solar
 * @param dt timestep
 * @param bisection_flag Flag to identify if scheme failed to converge
 */
INLINE static float newton_iter(
    float logu_init, double u_ini, int n_H_index, float d_n_H, int He_index,
    float d_He, float He_reion_heat, struct part *restrict p,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const float abundance_ratio[chemistry_element_count + 2], float dt,
    int *bisection_flag) {

  double logu, logu_old;
  double dLambdaNet_du = 0.0, LambdaNet;

  /* table bounds */
  const float log_table_bound_high =
      (cooling->Therm[eagle_cooling_N_temperature - 1] - 0.05) / M_LOG10E;
  const float log_table_bound_low = (cooling->Therm[0] + 0.05) / M_LOG10E;

  /* convert Hydrogen mass fraction in Hydrogen number density */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const double n_H =
      hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* compute ratefact = n_H * n_H / rho; Might lead to round-off error:
   * replaced by equivalent expression below */
  const double ratefact_cgs = n_H_cgs * XH * cooling->inv_proton_mass_cgs;

  logu_old = logu_init;
  logu = logu_old;
  int i = 0;

  float LambdaNet_old = 0;
  LambdaNet = 0;
  do /* iterate to convergence */
  {
    logu_old = logu;
    LambdaNet_old = LambdaNet;
    LambdaNet = (He_reion_heat / (dt * ratefact_cgs)) +
                eagle_cooling_rate(logu_old, cosmo->z, n_H_cgs, abundance_ratio,
                                   n_H_index, d_n_H, He_index, d_He, cooling,
                                   &dLambdaNet_du);

    /* Newton iteration. For details on how the cooling equation is integrated
     * see documentation in theory/Cooling/ */
    logu = logu_old - (1.0 - u_ini * exp(-logu_old) -
                       LambdaNet * ratefact_cgs * dt * exp(-logu_old)) /
                          (1.0 - dLambdaNet_du * ratefact_cgs * dt);
    /* Check if first step passes over equilibrium solution, if it does adjust
     * next guess */
    if (i == 1 && LambdaNet_old * LambdaNet < 0) logu = newton_log_u_guess_cgs;

    /* check whether iterations go within about 10% of the table bounds,
     * if they do default to bisection method */
    if (logu > log_table_bound_high) {
      i = newton_max_iterations;
      break;
    } else if (logu < log_table_bound_low) {
      i = newton_max_iterations;
      break;
    }

    i++;
  } while (fabs(logu - logu_old) > newton_tolerance &&
           i < newton_max_iterations);
  if (i >= newton_max_iterations) {
    /* flag to trigger bisection scheme */
    *bisection_flag = 1;
  }

  return logu;
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
    int n_H_index, float d_n_H, int He_index, float d_He,
    double Lambda_He_reion_cgs, double ratefact_cgs,
    const struct cooling_function_data *restrict cooling,
    const float abundance_ratio[chemistry_element_count + 2], double dt_cgs,
    long long ID) {

  /* Bracketing */
  double u_lower_cgs = u_ini_cgs;
  double u_upper_cgs = u_ini_cgs;

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      eagle_cooling_rate(log(u_ini_cgs), redshift, n_H_cgs, abundance_ratio,
                         n_H_index, d_n_H, He_index, d_He, cooling,
                         /*dLambdaNet_du=*/NULL);

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  if (LambdaNet_cgs < 0) {

    /* we're cooling! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        eagle_cooling_rate(log(u_lower_cgs), redshift, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, He_index, d_He, cooling,
                           /*dLambdaNet_du=*/NULL);

    int i = 0;
    while (u_lower_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs >
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs /= bracket_factor;
      u_upper_cgs /= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs = Lambda_He_reion_cgs +
                      eagle_cooling_rate(log(u_lower_cgs), redshift, n_H_cgs,
                                         abundance_ratio, n_H_index, d_n_H,
                                         He_index, d_He, cooling,
                                         /*dLambdaNet_du=*/NULL);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling",
          ID);
    }
  } else {

    /* we are heating! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        eagle_cooling_rate(log(u_upper_cgs), redshift, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, He_index, d_He, cooling,
                           /*dLambdaNet_du=*/NULL);

    int i = 0;
    while (u_upper_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs <
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs *= bracket_factor;
      u_upper_cgs *= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs = Lambda_He_reion_cgs +
                      eagle_cooling_rate(log(u_upper_cgs), redshift, n_H_cgs,
                                         abundance_ratio, n_H_index, d_n_H,
                                         He_index, d_He, cooling,
                                         /*dLambdaNet_du=*/NULL);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "heating",
          ID);
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
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        eagle_cooling_rate(log(u_next_cgs), redshift, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, He_index, d_He, cooling,
                           /*dLambdaNet_du=*/NULL);

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
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The cooling time-step of this particle.
 * @param dt_therm The hydro time-step of this particle.
 */
void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm) {

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

  /* Get internal energy at the end of the next kick step (assuming dt does not
   * increase) */
  double u_0 = (u_start + hydro_du_dt * dt_therm);

  /* Check for minimal energy */
  u_0 = max(u_0, hydro_properties->minimal_internal_energy);

  /* Convert to CGS units */
  const double u_start_cgs = u_start * cooling->internal_energy_to_cgs;
  const double u_0_cgs = u_0 * cooling->internal_energy_to_cgs;
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Change in redshift over the course of this time-step
     (See cosmology theory document for the derivation) */
  const double delta_redshift = -dt * cosmo->H * cosmo->a_inv;

  /* Get this particle's abundance ratios compared to solar
   * Note that we need to add S and Ca that are in the tables but not tracked
   * by the particles themselves.
   * The order is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe] */
  float abundance_ratio[chemistry_element_count + 2];
  abundance_ratio_to_solar(p, cooling, abundance_ratio);

  /* Get the Hydrogen and Helium mass fractions */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float XHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He];

  /* Get the Helium mass fraction. Note that this is He / (H + He), i.e. a
   * metal-free Helium mass fraction as per the Wiersma+08 definition */
  const float HeFrac = XHe / (XH + XHe);

  /* convert Hydrogen mass fraction into Hydrogen number density */
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
  get_index_1d(cooling->HeFrac, eagle_cooling_N_He_frac, HeFrac, &He_index,
               &d_He);
  get_index_1d(cooling->nH, eagle_cooling_N_density, log10(n_H_cgs), &n_H_index,
               &d_n_H);

  /* Start by computing the cooling (heating actually) rate from Helium
     re-ionization as this needs to be added on no matter what */

  /* Get helium and hydrogen reheating term */
  const double Helium_reion_heat_cgs = eagle_helium_reionization_extraheat(
      cooling->z_index, delta_redshift, cooling);

  /* Convert this into a rate */
  const double Lambda_He_reion_cgs =
      Helium_reion_heat_cgs / (dt_cgs * ratefact_cgs);

  /* Let's compute the internal energy at the end of the step */
  double u_final_cgs;

  /* First try an explicit integration (note we ignore the derivative) */
  const double LambdaNet_cgs =
      Lambda_He_reion_cgs + eagle_cooling_rate(log(u_0_cgs), cosmo->z, n_H_cgs,
                                               abundance_ratio, n_H_index,
                                               d_n_H, He_index, d_He, cooling,
                                               /*dLambdaNet_du=*/NULL);

  /* if cooling rate is small, take the explicit solution */
  if (fabs(ratefact_cgs * LambdaNet_cgs * dt_cgs) <
      explicit_tolerance * u_0_cgs) {

    u_final_cgs = u_0_cgs + ratefact_cgs * LambdaNet_cgs * dt_cgs;

  } else {

    int bisection_flag = 1;

    // MATTHIEU: TO DO restore the Newton-Raphson scheme
    if (0 && cooling->newton_flag) {

      /* Ok, try a Newton-Raphson scheme instead */
      double log_u_final_cgs =
          newton_iter(log(u_0_cgs), u_0_cgs, n_H_index, d_n_H, He_index, d_He,
                      Lambda_He_reion_cgs, p, cosmo, cooling, phys_const,
                      abundance_ratio, dt_cgs, &bisection_flag);

      /* Check if newton scheme sent us to a higher energy despite being in
         a  cooling regime If it did try newton scheme with a better guess.
         (Guess internal energy near equilibrium solution).  */
      if (LambdaNet_cgs < 0 && log_u_final_cgs > log(u_0_cgs)) {
        bisection_flag = 0;
        log_u_final_cgs =
            newton_iter(newton_log_u_guess_cgs, u_0_cgs, n_H_index, d_n_H,
                        He_index, d_He, Lambda_He_reion_cgs, p, cosmo, cooling,
                        phys_const, abundance_ratio, dt_cgs, &bisection_flag);
      }

      u_final_cgs = exp(log_u_final_cgs);
    }

    /* Alright, all else failed, let's bisect */
    if (bisection_flag || !(cooling->newton_flag)) {
      u_final_cgs =
          bisection_iter(u_0_cgs, n_H_cgs, cosmo->z, n_H_index, d_n_H, He_index,
                         d_He, Lambda_He_reion_cgs, ratefact_cgs, cooling,
                         abundance_ratio, dt_cgs, p->id);
    }
  }

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u_cgs = u_final_cgs - u_start_cgs;

  /* Convert back to internal units */
  double delta_u = delta_u_cgs * cooling->internal_energy_from_cgs;

  /* We now need to check that we are not going to go below any of the limits */

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_floor = gas_internal_energy_from_entropy(rho, A_floor);

  /* Absolute minimum */
  const double u_minimal = hydro_properties->minimal_internal_energy;

  /* Largest of both limits */
  const double u_limit = max(u_minimal, u_floor);

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_u. */
  if (u_start + 1.5 * delta_u < u_limit) {
    delta_u = (u_limit - u_start) / 1.5;
  }

  /* Second, check whether the energy used in the prediction could get negative.
   * We need to check for the 1/2 dt kick followed by a full time-step drift
   * that could potentially be for a time-step twice as big. We hence check
   * for 2.5 delta_u but this time against 0 energy not the minimum.
   * To avoid numerical rounding bringing us below 0., we add a tiny tolerance.
   */
  if (u_start + 2.5 * delta_u < 0.) {
    delta_u = -u_start / (2.5 + rounding_tolerance);
  }

  /* Turn this into a rate of change (including cosmology term) */
  const float cooling_du_dt = delta_u / dt_therm;

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cooling_du_dt * dt;
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
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
__attribute__((always_inline)) INLINE void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
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
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float XHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He];

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
  get_index_1d(cooling->HeFrac, eagle_cooling_N_He_frac, HeFrac, &He_index,
               &d_He);
  get_index_1d(cooling->nH, eagle_cooling_N_density, log10(n_H_cgs), &n_H_index,
               &d_n_H);

  /* Compute the log10 of the temperature by interpolating the table */
  const double log_10_T = eagle_convert_u_to_temp(
      log10(u_cgs), cosmo->z, /*compute_dT_du=*/0, /*dT_du=*/NULL, n_H_index,
      He_index, d_n_H, d_He, cooling);

  /* Undo the log! */
  return exp10(log_10_T);
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp #xpart data struct
 */
__attribute__((always_inline)) INLINE float cooling_get_radiated_energy(
    const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Initialises properties stored in the cooling_function_data struct
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const #phys_const data structure
 * @param cooling #cooling_function_data struct to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          struct cooling_function_data *cooling) {

  /* read some parameters */
  parser_get_param_string(parameter_file, "EAGLECooling:dir_name",
                          cooling->cooling_table_path);
  cooling->H_reion_z =
      parser_get_param_float(parameter_file, "EAGLECooling:H_reion_z");
  cooling->He_reion_z_centre =
      parser_get_param_float(parameter_file, "EAGLECooling:He_reion_z_centre");
  cooling->He_reion_z_sigma =
      parser_get_param_float(parameter_file, "EAGLECooling:He_reion_z_sigma");
  cooling->He_reion_heat_cgs =
      parser_get_param_float(parameter_file, "EAGLECooling:He_reion_eV_p_H");

  /* Optional parameters to correct the abundances */
  cooling->Ca_over_Si_ratio_in_solar = parser_get_opt_param_float(
      parameter_file, "EAGLECooling:Ca_over_Si_in_solar", 1.f);
  cooling->S_over_Si_ratio_in_solar = parser_get_opt_param_float(
      parameter_file, "EAGLECooling:S_over_Si_in_solar", 1.f);

  /* Convert to cgs (units used internally by the cooling routines) */
  cooling->He_reion_heat_cgs *=
      phys_const->const_electron_volt *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Read in the list of redshifts */
  get_cooling_redshifts(cooling);

  /* Read in cooling table header */
  char fname[eagle_table_path_name_length + 12];
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
  cooling->previous_z_index = eagle_cooling_N_redshifts - 2;

  /* Check if we are running with the newton scheme */
  cooling->newton_flag = parser_get_opt_param_int(
      parameter_file, "EAGLECooling:newton_integration", 0);
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
  char fname[eagle_table_path_name_length + 12];
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  read_cooling_header(fname, cooling);

  /* Allocate memory for the tables */
  allocate_cooling_tables(cooling);

  /* Force a re-read of the cooling tables */
  cooling->z_index = -10;
  cooling->previous_z_index = eagle_cooling_N_redshifts - 2;
  cooling_update(cosmo, cooling);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling #cooling_function_data struct.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  message("Cooling function is 'EAGLE'.");
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
  free(cooling->Redshifts);
  free(cooling->nH);
  free(cooling->Temp);
  free(cooling->HeFrac);
  free(cooling->Therm);
  free(cooling->SolarAbundances);
  free(cooling->SolarAbundances_inv);

  /* Free the tables */
  free(cooling->table.metal_heating);
  free(cooling->table.electron_abundance);
  free(cooling->table.temperature);
  free(cooling->table.H_plus_He_heating);
  free(cooling->table.H_plus_He_electron_abundance);
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
