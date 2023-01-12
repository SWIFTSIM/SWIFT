//
// Created by yuyttenh on 10/01/23.
//

#ifndef SWIFTSIM_COOLING_DE_RIJCKE_H
#define SWIFTSIM_COOLING_DE_RIJCKE_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "cooling_properties.h"
#include "cooling_tables.h"
#include "error.h"
#include "hydro.h"

/* Maximum number of iterations for bisection scheme */
static const int bisection_max_iterations = 150;

/* Tolerances for termination criteria. */
static const float explicit_tolerance = 0.05;
static const float bisection_tolerance = 1.0e-6;
static const double bracket_factor = 1.5;

INLINE static float cooling_get_temperature_from_u(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling, const double u);

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
  // Add content if required.
}

/**
 * @brief Finds the index of a value in a table and compute delta to nearest
 * element.
 *
 * This function assumes the table is monotonically increasing with a constant
 * difference between adjacent values.
 *
 * The returned difference is expressed in units of the table separation. This
 * means dx = (x - table[i]) / (table[i+1] - table[i]). It is always between
 * 0 and 1.
 *
 * We use a small epsilon of 1e-4 to avoid out-of-range accesses due to
 * rounding errors.
 *
 * @param table The table to search in.
 * @param size The number of elements in the table.
 * @param x The value to search for.
 * @param i (return) The index in the table of the element.
 * @param *dx (return) The difference between x and table[i]
 */
__attribute__((always_inline)) INLINE void get_index_1d(
    const float* restrict table, const int size, const float x, int* i,
    float* restrict dx) {

  /* Small epsilon to avoid rounding issues leading to out-of-bound
   * access when using the indices later to read data from the tables. */
  const float epsilon = 1e-4f;

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Distance between elements in the array */
  const float delta = (size - 1) / (table[size - 1] - table[0]);

  if (x < table[0] + epsilon) {
    /* We are below the first element */
    *i = 0;
    *dx = 0.f;
  } else if (x < table[size - 1] - epsilon) {
    /* Normal case */
    *i = (x - table[0]) * delta;

#ifdef SWIFT_DEBUG_CHECKS
    if (*i > size || *i < 0) {
      error(
          "trying to get index for value outside table range. Table size: %d, "
          "calculated index: %d, value: %.5e, table[0]: %.5e, grid size: %.5e",
          size, *i, x, table[0], delta);
    }
#endif

    *dx = (x - table[*i]) * delta;
  } else {
    /* We are after the last element */
    *i = size - 2;
    *dx = 1.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (*dx < -0.001f || *dx > 1.001f) error("Invalid distance found dx=%e", *dx);
#endif
}

__attribute__((always_inline)) INLINE static float interp_1d(
    const float* restrict table, int i, float dx) {
  return table[i] * (1.f - dx) + table[i + 1] * dx;
}

INLINE static double bisection_iter(
    const double u_ini_cgs, const double cooling_rate_ini,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props,
    const struct cooling_function_data* restrict cooling, const double dt_cgs,
    const double rate_fact, const long long ID) {

  /* Bracketing */
  double u_lower_cgs = u_ini_cgs / bracket_factor;
  double u_upper_cgs = u_ini_cgs * bracket_factor;

  double cooling_rate = cooling_rate_ini;
  int T_index;
  float d_T;

  /* Try to bracket the solution */
  if (cooling_rate < 0.) {
    /* Compute a new rate */
    float T_lower = cooling_get_temperature_from_u(
        phys_const, hydro_props, us, cosmo, cooling,
        u_lower_cgs * cooling->internal_energy_from_cgs);
    get_index_1d(cooling->table.temperature, de_rijcke_cooling_N_temperatures,
                 T_lower, &T_index, &d_T);

    /* New cooling rate */
    cooling_rate =
        rate_fact * interp_1d(cooling->table.cooling_rate, T_index, d_T);

    int i = 0;
    while (u_lower_cgs - u_ini_cgs - cooling_rate * dt_cgs > 0 &&
           i < bisection_max_iterations) {
      u_lower_cgs /= bracket_factor;
      u_upper_cgs /= bracket_factor;

      /* Compute a new rate */
      T_lower = cooling_get_temperature_from_u(
          phys_const, hydro_props, us, cosmo, cooling,
          u_lower_cgs * cooling->internal_energy_from_cgs);
      get_index_1d(cooling->table.temperature, de_rijcke_cooling_N_temperatures,
                   T_lower, &T_index, &d_T);

      /* New cooling rate */
      cooling_rate =
          rate_fact * interp_1d(cooling->table.cooling_rate, T_index, d_T);
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling, u_ini_cgs %.5e",
          ID, u_ini_cgs);
    }
  } else {
    error("We only implemented cooling at this point!");
  }

  /* bisection iteration */
  int i = 0;
  double u_next_cgs;
  float T_next;

  do {

    /* New guess */
    u_next_cgs = 0.5 * (u_lower_cgs + u_upper_cgs);
    T_next = cooling_get_temperature_from_u(
        phys_const, hydro_props, us, cosmo, cooling,
        u_next_cgs * cooling->internal_energy_from_cgs);

    get_index_1d(cooling->table.temperature, de_rijcke_cooling_N_temperatures,
                 T_next, &T_index, &d_T);

    /* New cooling rate */
    cooling_rate =
        rate_fact * interp_1d(cooling->table.cooling_rate, T_index, d_T);
#ifdef SWIFT_DEBUG_CHECKS
    if (u_next_cgs <= 0)
      error(
          "Got negative energy! u_next_cgs=%.5e u_upper=%.5e u_lower=%.5e "
          "Lambda=%.5e",
          u_next_cgs, u_upper_cgs, u_lower_cgs, LambdaNet_cgs);
#endif

    /* Where do we go next? */
    if (u_next_cgs - u_ini_cgs - cooling_rate * dt_cgs > 0.0) {
      u_upper_cgs = u_next_cgs;
    } else {
      u_lower_cgs = u_next_cgs;
    }

    i++;
  } while (fabs(u_upper_cgs - u_lower_cgs) / u_next_cgs > bisection_tolerance &&
           i < bisection_max_iterations);

  if (i >= bisection_max_iterations)
    error("Particle id %llu failed to converge", ID);

  return cooling_rate;
}

/**
 * @brief Calculates du/dt in CGS units for a particle.
 *
 * The cooling rate is \f$\frac{du}{dt} = -\frac{\Lambda}{n_H^2}
 * \frac{n_H^2}{\rho} \f$, where \f$ \frac{\Lambda}{n_H^2} \f$ is a constant in
 * this model (lambda_nH2_cgs in #cooling_function_data).
 * The returned value is in physical [erg * g^-1 * s^-1].
 *
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @return The change in energy per unit mass due to cooling for this particle
 * in cgs units [erg * g^-1 * s^-1].
 */
__attribute__((always_inline)) INLINE static double cooling_rate_cgs(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const struct cosmology* cosmo,
    const struct hydro_props* hydro_properties,
    const struct entropy_floor_properties* floor_props,
    const struct cooling_function_data* cooling, struct part* p,
    struct xpart* xp, const float dt, const float dt_therm) {

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

  /* Get the temperature */
  const float T = cooling_get_temperature_from_u(phys_const, hydro_properties,
                                                 us, cosmo, cooling, u_0);

  /* Get correction factor (table is for n_H = 1 cm^-3) */
  const double n_H =
      hydro_get_physical_density(p, cosmo) / phys_const->const_proton_mass;
  const double cm_3 =
      us->UnitLength_in_cgs * us->UnitLength_in_cgs * us->UnitLength_in_cgs;
  const double one_over_cm_6 = 1. / (cm_3 * cm_3);
  const double rate_fact = n_H * n_H * one_over_cm_6;

  int T_index;
  float d_T;
  get_index_1d(cooling->table.temperature, de_rijcke_cooling_N_temperatures, T,
               &T_index, &d_T);

  /* Compute the cooling rate */
  double cooling_rate =
      rate_fact * interp_1d(cooling->table.cooling_rate, T_index, d_T);

  /* if cooling rate is small, take the explicit solution */
  if (fabs(cooling_rate * dt_cgs) > explicit_tolerance * u_0_cgs) {
    /* Otherwise, go the bisection route. */
    cooling_rate =
        bisection_iter(u_0_cgs, cooling_rate, phys_const, us, cosmo,
                       hydro_properties, cooling, dt_cgs, rate_fact, p->id);
  }

  return cooling_rate;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time Time since Big Bang (or start of the simulation) in internal
 * units.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, const float dt,
    const float dt_therm, const double time) {

  /* Nothing to do here? */
  if (dt == 0.) return;

  /* Current energy (in internal units) */
  const float u_old_com = hydro_get_comoving_internal_energy(p, xp);

  /* Y' = RHS of the comoving equation for du/dt that will be integrated
     forward in time using dt_therm */
  const float hydro_du_dt_com = hydro_get_comoving_internal_energy_dt(p);

  /* Calculate cooling du_dt (in cgs units) */
  const double cooling_du_dt_cgs =
      cooling_rate_cgs(phys_const, us, cosmo, hydro_props, floor_props, cooling,
                       p, xp, dt, dt_therm);

  /* Convert to internal units */
  const float cooling_du_dt_physical =
      cooling_du_dt_cgs * cooling->conv_factor_energy_rate_from_cgs;

  /* Add cosmological term to get Y_cooling' */
  const float cooling_du_dt = cooling_du_dt_physical * cosmo->a * cosmo->a /
                              cosmo->a_factor_internal_energy;

  /* Y_total' */
  float total_du_dt = hydro_du_dt_com + cooling_du_dt;

  /* We now need to check that we are not going to go below any of the limits */

  /* Limit imposed by the entropy floor (comoving)
   * (Recall entropy is the same in physical and comoving frames) */
  const float A_floor_com = entropy_floor(p, cosmo, floor_props);
  const float rho_com = hydro_get_comoving_density(p);
  const float u_floor_com =
      gas_internal_energy_from_entropy(rho_com, A_floor_com);

  /* Absolute minimum */
  const float u_minimal_com =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Largest of both limits */
  const float u_limit_com = max(u_minimal_com, u_floor_com);

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_t. */
  if (u_old_com + total_du_dt * 1.5 * dt_therm < u_limit_com) {
    total_du_dt = (u_limit_com - u_old_com) / (1.5f * dt_therm);
  }

  /* Second, check whether the energy used in the prediction could get negative.
   * We need to check for the 1/2 dt kick followed by a full time-step drift
   * that could potentially be for a time-step twice as big. We hence check
   * for 2.5 delta_t but this time against 0 energy not the minimum */
  if (u_old_com + total_du_dt * 2.5 * dt_therm < 0.) {
    total_du_dt = -u_old_com / ((2.5f + 0.0001f) * dt_therm);
  }

  if (cooling->rapid_cooling) {
    const float u_new_com = u_old_com + total_du_dt * dt_therm;
    const float u_new_phys = u_new_com * cosmo->a_factor_internal_energy;
    hydro_set_physical_internal_energy(p, xp, cosmo, u_new_phys);
    hydro_set_drifted_physical_internal_energy(p, cosmo, u_new_phys);
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.);
  } else {
    /* Update the internal energy time derivative */
    hydro_set_comoving_internal_energy_dt(p, total_du_dt);
  }

  const float actual_cooling_du_dt = total_du_dt - hydro_du_dt_com;
  const float actual_cooling_du_dt_physical = actual_cooling_du_dt / cosmo->a /
                                              cosmo->a *
                                              cosmo->a_factor_internal_energy;
  /* Store the radiated energy (assuming dt will not change) */
  xp->cooling_data.radiated_energy +=
      -hydro_get_mass(p) * actual_cooling_du_dt_physical * dt;
}

/**
 * @brief Computes the time-step due to cooling for this particle.
 *
 * The time-step is not set by the properties of cooling.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended data of the particle.
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
 * @brief Compute the electron pressure of a #part based on the cooling
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
__attribute__((always_inline)) INLINE static double
cooling_get_electron_pressure(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp) {
  return 0.;
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
__attribute__((always_inline)) INLINE static double cooling_get_ycompton(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp) {
  error("This cooling model does not compute Compton Y!");
  return 0.;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here. Just set the radiated energy counter to 0.
 *
 * @param phys_const The physical constants in internal units.
 * @param cooling The properties of the cooling function.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
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
INLINE static float cooling_get_temperature_from_u(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling, const double u) {

  /* Physical constants */
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  /* Gas properties */
  const double T_transition = hydro_props->hydrogen_ionization_temperature;
  const double mu_neutral = hydro_props->mu_neutral;
  const double mu_ionised = hydro_props->mu_ionised;

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
INLINE static float cooling_get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {
  return cooling_get_temperature_from_u(
      phys_const, hydro_props, us, cosmo, cooling,
      hydro_get_physical_internal_energy(p, xp, cosmo));
}

/**
 * @brief Returns the subgrid temperature of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_temperature(const struct part* p,
                                                    const struct xpart* xp) {
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
INLINE static float cooling_get_subgrid_density(const struct part* p,
                                                const struct xpart* xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
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
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
static INLINE void cooling_split_part(struct part* p, struct xpart* xp,
                                      double n) {

  xp->cooling_data.radiated_energy /= n;
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
static INLINE void cooling_init_backend(struct swift_params* parameter_file,
                                        const struct unit_system* us,
                                        const struct phys_const* phys_const,
                                        const struct hydro_props* hydro_props,
                                        struct cooling_function_data* cooling) {
  /* Read in parameters */
  cooling->rapid_cooling = parser_get_opt_param_int(
      parameter_file, "DeRijckeCooling:rapid_cooling", 0);

  /* Directory for cooling tables */
  parser_get_param_string(parameter_file, "DeRijckeCooling:dir_name",
                          cooling->cooling_table_path);

  /* Read the actual tables */
  get_cooling_tables(cooling);

  /* Some useful conversion values */
  cooling->conv_factor_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->conv_factor_energy_rate_from_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Useful constants */
  cooling->proton_mass_cgs_inv =
      1. / (phys_const->const_proton_mass *
            units_cgs_conversion_factor(us, UNIT_CONV_MASS));

  /* Compute conversion factors */
  cooling->internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  cooling->internal_energy_from_cgs = 1. / cooling->internal_energy_to_cgs;
  cooling->number_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
}

/**
 * @brief Restore cooling tables (if applicable) after
 * restart
 *
 * @param cooling the cooling_function_data structure
 * @param cosmo cosmology structure
 */
static INLINE void cooling_restore_tables(struct cooling_function_data* cooling,
                                          const struct cosmology* cosmo) {
  get_cooling_tables(cooling);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'De Rijcke et al. (2013) Cooling tables'");

  if (cooling->rapid_cooling) {
    message("Using rapid cooling");
  } else {
    message("Using normal cooling");
  }
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
static INLINE void cooling_clean(struct cooling_function_data* cooling) {
  swift_free("cooling", cooling->table.temperature);
  swift_free("cooling", cooling->table.cooling_rate);
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

#endif  // SWIFTSIM_COOLING_DE_RIJCKE_H
