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
#include "eagle_cool_tables.h"
#include "error.h"
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
static const double bracket_factor = sqrt(1.1); /* Taken from EAGLE */
static const double newton_log_u_guess_cgs = log(2.e12);

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param restart_flag Flag indicating restarted run.
 */
void cooling_update(const struct cosmology *cosmo,
                    struct cooling_function_data *cooling,
                    const int restart_flag) {

  /* Current redshift */
  const float redshift = cosmo->z;

  /* Get index along the redshift index of the tables */
  int z_index = -1;
  float dz = 0.f;
  if (redshift > cooling->reionisation_redshift) {
    z_index = -2;
  } else if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    z_index = -1;
  } else {
    get_redshift_index(redshift, &z_index, &dz, cooling);
  }
  cooling->z_index = z_index;
  cooling->dz = dz;

  eagle_check_cooling_tables(cooling, restart_flag);
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
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
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
    float logu_init, double u_ini, int n_h_i, float d_n_h, int He_i, float d_He,
    float He_reion_heat, struct part *restrict p,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const, float *abundance_ratio,
    float dt, int *bisection_flag) {

  double logu, logu_old;
  double dLambdaNet_du = 0.0, LambdaNet;

  /* table bounds */
  const float log_table_bound_high =
      (cooling->Therm[cooling->N_Temp - 1] - 0.05) / M_LOG10E;
  const float log_table_bound_low = (cooling->Therm[0] + 0.05) / M_LOG10E;

  /* convert Hydrogen mass fraction in Hydrogen number density */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const double n_h = hydro_get_physical_density(p, cosmo) * XH /
                     phys_const->const_proton_mass *
                     cooling->number_density_scale;

  /* compute ratefact = n_h * n_h / rho; Might lead to round-off error:
   * replaced by equivalent expression below */
  const double ratefact = n_h * (XH / cooling->proton_mass_cgs);

  logu_old = logu_init;
  logu = logu_old;
  int i = 0;

  float LambdaNet_old = 0;
  LambdaNet = 0;
  do /* iterate to convergence */
  {
    logu_old = logu;
    LambdaNet_old = LambdaNet;
    LambdaNet =
        (He_reion_heat / (dt * ratefact)) +
        eagle_cooling_rate(logu_old, &dLambdaNet_du, n_h_i, d_n_h, He_i, d_He,
                           p, cooling, cosmo, phys_const, abundance_ratio);

    /* Newton iteration. For details on how the cooling equation is integrated
     * see documentation in theory/Cooling/ */
    logu = logu_old - (1.0 - u_ini * exp(-logu_old) -
                       LambdaNet * ratefact * dt * exp(-logu_old)) /
                          (1.0 - dLambdaNet_du * ratefact * dt);
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
 * @brief Bisection integration scheme used when Newton Raphson method fails to
 * converge
 *
 * @param logu_init Initial guess for log(internal energy)
 * @param u_ini Internal energy at beginning of hydro step
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p #part structure
 * @param cosmo #cosmology structure
 * @param cooling #cooling_function_data structure
 * @param phys_const #phys_const data structure
 * @param abundance ratio array of ratios of metal abundance to solar
 * @param dt timestep
 */
INLINE static float bisection_iter(
    float logu_init, double u_ini, int n_h_i, float d_n_h, int He_i, float d_He,
    float He_reion_heat, struct part *restrict p,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const, float *abundance_ratio,
    float dt) {
  double u_upper, u_lower, u_next, LambdaNet;
  double *dLambdaNet_du = NULL;
  int i = 0;
  double u_init = exp(logu_init);

  /* convert Hydrogen mass fraction in Hydrogen number density */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const double n_h = hydro_get_physical_density(p, cosmo) * XH /
                     phys_const->const_proton_mass *
                     cooling->number_density_scale;

  /* compute ratefact = n_h * n_h / rho; Might lead to round-off error:
   * replaced by equivalent expression  below */
  const double ratefact = n_h * (XH / cooling->proton_mass_cgs);

  /* Bracketing */
  u_lower = u_init;
  u_upper = u_init;
  LambdaNet =
      (He_reion_heat / (dt * ratefact)) +
      eagle_cooling_rate(log(u_init), dLambdaNet_du, n_h_i, d_n_h, He_i, d_He,
                         p, cooling, cosmo, phys_const, abundance_ratio);

  i = 0;
  if (LambdaNet < 0) {
    /* we're cooling */
    u_lower /= bracket_factor;
    u_upper *= bracket_factor;

    LambdaNet = (He_reion_heat / (dt * ratefact)) +
                eagle_cooling_rate(log(u_lower), dLambdaNet_du, n_h_i, d_n_h,
                                   He_i, d_He, p, cooling, cosmo, phys_const,
                                   abundance_ratio);
    while (u_lower - u_ini - LambdaNet * ratefact * dt > 0 &&
           i < bisection_max_iterations) {
      u_lower /= bracket_factor;
      u_upper /= bracket_factor;
      LambdaNet = (He_reion_heat / (dt * ratefact)) +
                  eagle_cooling_rate(log(u_lower), dLambdaNet_du, n_h_i, d_n_h,
                                     He_i, d_He, p, cooling, cosmo, phys_const,
                                     abundance_ratio);
      i++;
    }
    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling",
          p->id);
    }
  } else {
    /* heating */
    u_lower /= bracket_factor;
    u_upper *= bracket_factor;

    LambdaNet = (He_reion_heat / (dt * ratefact)) +
                eagle_cooling_rate(log(u_upper), dLambdaNet_du, n_h_i, d_n_h,
                                   He_i, d_He, p, cooling, cosmo, phys_const,
                                   abundance_ratio);
    while (u_upper - u_ini - LambdaNet * ratefact * dt < 0 &&
           i < bisection_max_iterations) {
      u_lower *= bracket_factor;
      u_upper *= bracket_factor;
      LambdaNet = (He_reion_heat / (dt * ratefact)) +
                  eagle_cooling_rate(log(u_upper), dLambdaNet_du, n_h_i, d_n_h,
                                     He_i, d_He, p, cooling, cosmo, phys_const,
                                     abundance_ratio);
      i++;
    }
    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "heating",
          p->id);
    }
  }

  /* bisection iteration */
  i = 0;
  do {
    u_next = 0.5 * (u_lower + u_upper);
    LambdaNet =
        (He_reion_heat / (dt * ratefact)) +
        eagle_cooling_rate(log(u_next), dLambdaNet_du, n_h_i, d_n_h, He_i, d_He,
                           p, cooling, cosmo, phys_const, abundance_ratio);
    if (u_next - u_ini - LambdaNet * ratefact * dt > 0.0) {
      u_upper = u_next;
    } else {
      u_lower = u_next;
    }

    i++;
  } while (fabs(u_upper - u_lower) / u_next > bisection_tolerance &&
           i < bisection_max_iterations);

  if (i >= bisection_max_iterations) {
    error("Particle id %llu failed to converge", p->id);
  }

  return log(u_upper);
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_properties the hydro_props struct, used for
 * getting the minimal internal energy allowed in by SWIFT.
 * Read from yml file into engine struct.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt_cool The cooling time-step of this particle.
 * @param dt_therm The hydro time-step of this particle.
 */
void cooling_cool_part(const struct phys_const *restrict phys_const,
                       const struct unit_system *restrict us,
                       const struct cosmology *restrict cosmo,
                       const struct hydro_props *restrict hydro_properties,
                       const struct cooling_function_data *restrict cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm) {

  /* No cooling happens over zero time */
  if (dt == 0.) return;

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
  const double u_start_cgs = u_start * cooling->internal_energy_scale;
  const double u_0_cgs = u_0 * cooling->internal_energy_scale;
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Get this particle's abundance ratios */
  float abundance_ratio[chemistry_element_count + 2];
  abundance_ratio_to_solar(p, cooling, abundance_ratio);

  /* Get the H and He mass fractions */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float HeFrac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
      (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  /* convert Hydrogen mass fraction in Hydrogen number density */
  const double n_h = hydro_get_physical_density(p, cosmo) * XH /
                     phys_const->const_proton_mass *
                     cooling->number_density_scale;

  /* ratefact = n_h * n_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  const double ratefact = n_h * (XH / cooling->proton_mass_cgs);

  /* Get helium and hydrogen reheating term */
  const double LambdaTune = eagle_helium_reionization_extraheat(
      cooling->z_index, -dt * cosmo->H * cosmo->a_inv, cooling);

  /* compute hydrogen number density and helium fraction table indices and
   * offsets (These are fixed for of u, so no need to recompute them) */
  int He_i, n_h_i;
  float d_He, d_n_h;
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(n_h), &n_h_i, &d_n_h);

  /* Let's compute the internal energy at the end of the step */
  double u_final_cgs;

  /* First try an explicit integration (note we ignore the derivative) */
  const double LambdaNet =
      LambdaTune / (dt_cgs * ratefact) +
      eagle_cooling_rate(log(u_0_cgs), /*dLambdaNet_du=*/NULL, n_h_i, d_n_h,
                         He_i, d_He, p, cooling, cosmo, phys_const,
                         abundance_ratio);

  /* if cooling rate is small, take the explicit solution */
  if (fabs(ratefact * LambdaNet * dt_cgs) < explicit_tolerance * u_0_cgs) {

    u_final_cgs = u_0_cgs + ratefact * LambdaNet * dt_cgs;

  } else {

    int bisection_flag = 0;
    double log_u_final_cgs = 0.0;

    if (cooling->newton_flag) {
      /* Ok, try a Newton-Raphson scheme instead */
      log_u_final_cgs = newton_iter(
          log(u_0_cgs), u_0_cgs, n_h_i, d_n_h, He_i, d_He, LambdaTune, p, cosmo,
          cooling, phys_const, abundance_ratio, dt_cgs, &bisection_flag);

      /* Check if newton scheme sent us to a higher energy despite being in a
       * cooling regime If it did try newton scheme with a better guess. (Guess
       * internal energy near equilibrium solution).  */
      if (LambdaNet < 0 && log_u_final_cgs > log(u_0_cgs)) {
        bisection_flag = 0;
        log_u_final_cgs =
            newton_iter(newton_log_u_guess_cgs, u_0_cgs, n_h_i, d_n_h, He_i,
                        d_He, LambdaTune, p, cosmo, cooling, phys_const,
                        abundance_ratio, dt_cgs, &bisection_flag);
      }
    }

    /* Alright, all else failed, let's bisect */
    if (bisection_flag || !(cooling->newton_flag)) {
      log_u_final_cgs = bisection_iter(
          log(u_0_cgs), u_0_cgs, n_h_i, d_n_h, He_i, d_He, LambdaTune, p, cosmo,
          cooling, phys_const, abundance_ratio, dt_cgs);
    }

    u_final_cgs = exp(log_u_final_cgs);
  }

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u_cgs = u_final_cgs - u_start_cgs;

  /* Convert back to internal units */
  double delta_u = delta_u_cgs / cooling->internal_energy_scale;

  /* We now need to check that we are not going to go below any of the limits */

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_u. */
  if (u_start + 1.5 * delta_u < hydro_properties->minimal_internal_energy) {
    delta_u = (hydro_properties->minimal_internal_energy - u_start) / 1.5;
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
 * @param p #part data.
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
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @param cooling #cooling_function_data struct.
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
  parser_get_param_string(parameter_file, "EagleCooling:filename",
                          cooling->cooling_table_path);
  cooling->reionisation_redshift = parser_get_param_float(
      parameter_file, "EagleCooling:reionisation_redshift");
  cooling->calcium_over_silicon_ratio = parser_get_param_float(
      parameter_file, "EAGLEChemistry:CalciumOverSilicon");
  cooling->sulphur_over_silicon_ratio = parser_get_param_float(
      parameter_file, "EAGLEChemistry:SulphurOverSilicon");
  cooling->He_reion_z_centre =
      parser_get_param_float(parameter_file, "EagleCooling:He_reion_z_centre");
  cooling->He_reion_z_sigma =
      parser_get_param_float(parameter_file, "EagleCooling:He_reion_z_sigma");
  cooling->He_reion_ev_pH =
      parser_get_param_float(parameter_file, "EagleCooling:He_reion_ev_pH");

  /* convert to cgs */
  cooling->He_reion_ev_pH *= phys_const->const_electron_volt *
                             units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* read in cooling table header */
  get_cooling_redshifts(cooling);
  char fname[eagle_table_path_name_length + 12];
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  read_cooling_header(fname, cooling);

  /* Allocate space for cooling tables */
  allocate_cooling_tables(cooling);

  /* compute conversion factors */
  cooling->internal_energy_scale =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) /
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->number_density_scale =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) /
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);

  cooling->proton_mass_cgs = phys_const->const_proton_mass *
                             units_cgs_conversion_factor(us, UNIT_CONV_MASS);
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

  /* set low_z_index to -10 to indicate we haven't read any tables yet */
  cooling->low_z_index = -10;
  /* set previous_z_index and to last value of redshift table*/
  cooling->previous_z_index = cooling->N_Redshifts - 2;

  /* Check if we are running with the newton scheme */
  cooling->newton_flag = parser_get_opt_param_int(
      parameter_file, "EagleCooling:newton_integration", 0);
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

  /* Read relevant cooling tables.
   * Third variable in cooling_update flag to mark restart*/
  allocate_cooling_tables(cooling);
  cooling_update(cosmo, cooling, /*restart=*/1);
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

  /* Free the tables */
  free(cooling->table.metal_heating);
  free(cooling->table.electron_abundance);
  free(cooling->table.temperature);
  free(cooling->table.H_plus_He_heating);
  free(cooling->table.H_plus_He_electron_abundance);
}
