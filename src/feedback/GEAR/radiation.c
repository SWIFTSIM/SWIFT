/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
 * @file src/feedback/GEAR/radiation.c
 * @brief Subgrid radiation feedback for GEAR. This files contains functions to
 * compute quantities for the radiation feedback.
 */

/* Include header */
#include "radiation.h"

#include "interpolation.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"
#include "kernel_hydro.h"
#include "units.h"

/**
 * Get the gas number of hydrogen atoms.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Number of hydrogen atoms.
 */
__attribute__((always_inline)) INLINE
double radiation_get_part_number_hydrogen_atoms(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp) {

  const float m = hydro_get_mass(p);
  const double m_p = phys_const->const_proton_mass;
  const float X_H = cooling_get_hydrogen_mass_fraction(cooling, p, xp);
  const float mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);

  /* Number of hydrogen atoms in b */
  const double N_H = (X_H * m) / (mu * m_p);

  return N_H;
}

/**
 * Get the gas ionizing rate needed to fully ionize the #part.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Ionizing photon rate to ionize this #part (physical units).
 */
__attribute__((always_inline)) INLINE
double radiation_get_part_rate_to_fully_ionize(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp) {

  /* const float m = hydro_get_mass(p); */
  const float rho = hydro_get_physical_density(p, cosmo);
  const double beta = phys_const->const_caseb_recomb;
  const double m_e = phys_const->const_electron_mass;
  const float X_H = cooling_get_hydrogen_mass_fraction(cooling, p, xp);
  const float mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);

  /* Number of hydrogen atoms in b */
  const double N_H = radiation_get_part_number_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);

  /* Electron density assuming full ionization */
  const double n_e = (X_H * rho) / (mu * m_e);

  /* Required ionizing rate in [photons / internal time unit] */
  const double Delta_N_dot = N_H * beta * n_e;

  return Delta_N_dot;
}

/**
 * Get the #spart ionization photon emission rate.
 *
 * @param sp The star.
 * @return Ionizing photon rate.
 */
__attribute__((always_inline)) INLINE
double radiation_get_star_ionization_rate(const struct spart* sp) {
  return sp->feedback_data.radiation.dot_N_ion;
}

/**
 * Consume the #spart ionizing photon budget.
 *
 * @param sp The star.
 * @param Delta_dot_N_ion The ionizing photon rate to remove.
 */
__attribute__((always_inline)) INLINE
void radiation_consume_ionizing_photons(struct spart* sp,
                                        double Delta_dot_N_ion) {
  sp->feedback_data.radiation.dot_N_ion -= Delta_dot_N_ion;
  return;
}

/**
 * Tag the #part as ionized to be ionized in feedback_update_part().
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE
void radiation_tag_part_as_ionized(struct part* p, struct xpart* xp) {
  xp->feedback_data.radiation.is_ionized = 1;
  return;
}

/**
 * Reset the #part ionization tag.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE
void radiation_reset_part_ionized_tag(struct part* p, struct xpart* xp) {
  xp->feedback_data.radiation.is_ionized = 0;
  return;
}

/**
 * Is this #part *tagged* as ionized ?
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Is the particle *tagged* ionized?
 */
__attribute__((always_inline)) INLINE
int radiation_is_part_tagged_as_ionized(struct part* p, struct xpart* xp) {
  return xp->feedback_data.radiation.is_ionized;
}

/**
 * Determines whether a gas #part is ionized or not based on its
 * thermodynamical properties.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param hydro_properties The #hydro_props.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Is the particle ionized?
 */
__attribute__((always_inline)) INLINE
int radiation_is_part_ionized(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp) {

  /* Is T > 10^4 K ? */
  const float T = cooling_get_temperature(phys_const, hydro_props, us, cosmo,
                                          cooling, p, xp);
  const float ten_to_four_kelvin =
      1e4 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Is the particle ionized ? */
  return (T > ten_to_four_kelvin || xp->feedback_data.radiation.is_ionized);
}

/**
 * Compute the gas comoving column density at the star's location using the
 * Sobolev approximation.
 *
 * @param sp The #spart.
 * @return Comoving gas column density at the star's location.
 */
__attribute__((always_inline)) INLINE
float radiation_get_comoving_gas_column_density_at_star(const struct spart* sp) {
  const float rho_gas = sp->feedback_data.rho_star;
  const float grad_rho[3] = {sp->feedback_data.grad_rho_star[0],
                             sp->feedback_data.grad_rho_star[1],
                             sp->feedback_data.grad_rho_star[2]};
  const float norm_grad_rho =
      sqrtf(grad_rho[0] * grad_rho[0] + grad_rho[1] * grad_rho[1] +
            grad_rho[2] * grad_rho[2]);

  const float length_gas = sp->h * kernel_gamma + rho_gas / norm_grad_rho;
  return length_gas * rho_gas;
}

/**
 * Compute the physical infrared opacity around a star.
 *
 * @param sp The #spart.
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @return Infrared gas opacity around the star.
 */
__attribute__((always_inline)) INLINE
float radiation_get_physical_IR_opacity(const struct spart* sp,
                               const struct unit_system* us,
                               const struct phys_const* phys_const,
			       const struct cosmology* cosmo) {
  const float Z_gas = sp->feedback_data.Z_star;
  const float Z_sun = 0.02;
  const float value = 10.0 * units_cgs_conversion_factor(us, UNIT_CONV_AREA) /
                      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  return value * Z_gas / Z_sun;
}

/**
 * Compute the physical infrared optical depth around a star.
 *
 * @param sp The #spart.
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @return Infrared gas optical depth around the star.
 */
__attribute__((always_inline)) INLINE
float radiation_get_physical_IR_optical_depth(const struct spart* sp,
					      const struct unit_system* us,
					      const struct phys_const* phys_const,
					      const struct cosmology* cosmo) {
  const float Sigma_gas_c = radiation_get_comoving_gas_column_density_at_star(sp);
  const float Sigma_gas_p = Sigma_gas_c * cosmo->a2_inv;
  const float kappa_IR = radiation_get_physical_IR_opacity(sp, us, phys_const, cosmo);
  return kappa_IR * Sigma_gas_p;
}

/**
 * Compute the physical radiation pressure emitted by the star.
 *
 * @param sp The #spart.
 * @param Delta_t The current #spart timestep.
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @return Radiation pressure emittied by the star.
 */
__attribute__((always_inline)) INLINE
float radiation_get_star_physical_radiation_pressure(
    const struct spart* sp, const float Delta_t, const struct phys_const* phys_const,
    const struct unit_system* us, const struct cosmology* cosmo) {

  const float tau_IR = radiation_get_physical_IR_optical_depth(sp, us, phys_const, cosmo);
  const float L_bol = sp->feedback_data.radiation.L_bol; /* In physical units */
  const float c = phys_const->const_speed_light_c;

  return Delta_t * L_bol / c * (1 + tau_IR);
}

/**
 * Compute the radius of a single star from empirical mass-radius relations.
 *
 * This function get the value for an individual star. For a SSP, this function
 * is used to compute and IMF-average.
 *
 * @param sp Pointer to star particle.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Radius in code units.
 */
float radiation_get_individual_star_radius(
    const float mass, const struct unit_system* us,
    const struct phys_const* phys_const) {

  /* Perform some units conversions */
  const float R_sun = phys_const->const_solar_radius;
  const float M_solar = phys_const->const_solar_mass;
  const float M_in_solar = mass / M_solar;

  if (M_in_solar < 1.f) {
    return R_sun * powf(M_in_solar, 0.8f);
  } else if (M_in_solar < 8.f) {
    return R_sun * powf(M_in_solar, 0.57f);
  } else {
    return R_sun * powf(M_in_solar, 0.5f);
  }
}

/**
 * Compute the temperature of a single star from empirical mass-temperature
 * relations.
 *
 * This function get the value for an individual star. For a SSP, this function
 * is used to compute and IMF-average.
 *
 * @param sp Pointer to star particle.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Temperature in code units.
 */
float radiation_get_individual_star_temperature(
    const float mass, const struct unit_system* us,
    const struct phys_const* phys_const) {

  const float M_solar = phys_const->const_solar_mass;
  const float M_in_solar = mass / M_solar; /* In solar masses */

  float T_K = 0.0;

  if (M_in_solar < 1.f) {
    T_K = 3500.f * powf(M_in_solar, 0.5f);
  } else if (M_in_solar < 8.f) {
    T_K = 5800.f * powf(M_in_solar, 0.5f);
  } else {
    T_K = 25000.f * powf(M_in_solar / 20.f, 0.1f);
  }

  /* Convert from Kelvin to internal units using unit_system_temperature_in_cgs
   */
  const float T_internal =
      T_K / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  return T_internal;
}

/**
 * Computes the bolometric luminosity of a single star from empirical
 * mass-luminosity relations.
 *
 * This function get the value for an individual star. For a SSP, this function
 * is used to compute and IMF-average.
 *
 * @param sp Pointer to star particle.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Luminosity in code units.
 */
float radiation_get_individual_star_luminosity(
    const float mass, const struct unit_system* us,
    const struct phys_const* phys_const) {

  /* Convert mass to solar masses */
  const float M_in_solar = mass / phys_const->const_solar_mass;

  /* Piecewise empirical mass-luminosity relation */
  float lum_sol;
  if (M_in_solar < 0.43f) {
    lum_sol = 0.185f * M_in_solar * M_in_solar;
  } else if (M_in_solar < 2.0f) {
    lum_sol = M_in_solar * M_in_solar * M_in_solar * M_in_solar;
  } else if (M_in_solar < 54.0f) {
    lum_sol = 1.5f * M_in_solar * M_in_solar * M_in_solar * sqrtf(M_in_solar);
  } else {
    lum_sol = 32000.0f * M_in_solar;
  }

  /* Convert from solar luminosities to code units */
  const float luminosity = lum_sol * phys_const->const_solar_luminosity;
  return luminosity;
}

/**
 * @brief Get the #spart ionizing photon emission rate using a fit on the
 * Blackbody spectrum (from Gizmo).
 *
 * This function get the value for an individual star. For a SSP, this function
 * is used to compute and IMF-average.
 *
 * @param sp The #spart to consider.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @return N_dot_ion The ionizing photon emission rate in code units
 * [photons/U_T].
 */
double radiation_get_individual_star_ionizing_photon_emission_rate_fit(
    const float mass, const struct unit_system* us,
    const struct phys_const* phys_const) {

  /* Get star properties in internal units */
  /* const float T = radiation_get_individual_star_temperature(sp, us,
   * phys_const); */
  const float R = radiation_get_individual_star_radius(mass, us, phys_const);
  const float L = radiation_get_individual_star_luminosity(mass, us, phys_const);

  const float R_in_R_sun = R / phys_const->const_solar_radius;
  const float L_in_L_sun = L / phys_const->const_solar_luminosity;

  if (R <= 0.f || L <= 0.f) {
    // TODO: print a warning or an error
    return 0.f;
  }

  /* Get the Blackbody effective temperature in K */
  const double T_K = 5780 *
                     powf((L_in_L_sun) / (R_in_R_sun * R_in_R_sun), 0.25) /
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Compute dimensionless photon cutoff x_0 = h*nu / kT for 13.6 eV. Note :
     13.6 eV = 2.17872e-11 erg to be converted to internal unit */
  const double x_0 =
      (2.17872e-11 / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY)) /
      (phys_const->const_boltzmann_k * T_K);

  /* Fit for fraction of ionizing luminosity */
  if (x_0 < 30.f) {
    const double q =
        18.0 / (x_0 * x_0) + 1.0 / (8.0 + x_0 + 20.0 * expf(-x_0 / 10.0));
    const double f_ion = expf(-1.0 / q);

    /* Ionizing luminosity in internal units */
    const double L_ion =
        f_ion * L;  //* units_cgs_conversion_factor(us, UNIT_CONV_POWER);

    /* Assume average ionizing photon energy ~ 20 eV = 20 * 1.60218e-12 erg */
    const double photon_energy_cgs = 20.0 * 1.60218e-12;
    const double photon_energy_internal =
        photon_energy_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

    const double N_dot_ion = L_ion / photon_energy_internal;

    return N_dot_ion;
  } else {
    return 0.0;
  }
}


/******************************************************************************/
/* Functions to deal with integrated data over an IMF. These functions read,
   interpolate and integrate. */
/******************************************************************************/

/**
 * @brief Print the radiation model.
 *
 * @param rad The #radiation.
 */
void radiation_print(const struct radiation *rad) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  /* message("Mass range for RAD = [%g, %g]", rad->mass_min, rad->mass_max); */
}

/**
 * @brief Initialize the #radiation structure.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_init(struct radiation *rad, struct swift_params *params,
		    const struct stellar_model *sm,
		    const struct unit_system *us,
		    const struct phys_const* phys_const) {

  /* Read the data */
  radiation_read_data(rad, params, sm, us, phys_const, /* restart */ 0);
}

/**
 * @brief Write a radiation struct to the given FILE as a stream of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void radiation_dump(const struct radiation *rad, FILE *stream,
                        const struct stellar_model *sm) {

  restart_write_blocks((void*)rad,
                       sizeof(struct radiation), 1, stream,
                       "radiation", "radiation");
  message("Dumping GEAR radiation...");
}

/**
 * @brief Restore a radiation struct from the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void radiation_restore(struct radiation *rad, FILE *stream,
		       const struct stellar_model *sm) {

  restart_read_blocks((void*)rad,
                      sizeof(struct radiation), 1, stream, NULL,
                      "radiation");
  message("Restoring GEAR radiation struct...");
}

/**
 * @brief Clean the allocated memory.
 *
 * @param rad the #radiation.
 */
void radiation_clean(struct radiation *rad) {

  interpolate_1d_free(&rad->integrated.luminosities);
  interpolate_1d_free(&rad->raw.luminosities);
  interpolate_1d_double_free(&rad->integrated.dot_N_ion);
  interpolate_1d_double_free(&rad->raw.dot_N_ion);
}

/**
 * @brief Get the IMF-averaged nolometric luminosity per mass.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The bolometric luminosity.
 */
float radiation_get_luminosities_from_integral(const struct radiation *rad,
					      float log_m1, float log_m2) {

    float luminosity_1 = interpolate_1d(&rad->integrated.luminosities, log_m1);
    float luminosity_2 = interpolate_1d(&rad->integrated.luminosities, log_m2);
    return luminosity_2 - luminosity_1;
};

/**
 * @brief Get the IMF-averaged bolometric luminosity per mass.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param The bolometric luminosity.
 */
float radiation_get_luminosities_from_raw(const struct radiation *rad,
                                       float log_m) {
    return interpolate_1d(&rad->raw.luminosities, log_m);
};


/**
 * @brief Get the IMF-averaged ionization rate per mass.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The ionization rate;
 */
double radiation_get_ionization_rate_from_integral(const struct radiation *rad,
					      float log_m1, float log_m2) {

    double dot_N_ion_1 = interpolate_1d_double(&rad->integrated.dot_N_ion, log_m1);
    double dot_N_ion_2 = interpolate_1d_double(&rad->integrated.dot_N_ion, log_m2);
    return dot_N_ion_2 - dot_N_ion_1;
};

/**
 * @brief Get the non-IMF-integrated ionization rate per mass.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param The ionization rate;
 */
double radiation_get_ionization_rate_from_raw(const struct radiation *rad,
					      float log_m) {
    return interpolate_1d_double(&rad->raw.dot_N_ion, log_m);
};


/**
 * @brief Read an array of luminosities data from the table.
 *
 * @param rad The #radiation model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm * The #stellar_model.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size Number of element to keep in the interpolation
 * data.
 */
void radiation_read_luminosities_array(
    struct radiation *rad, struct interpolation_1d *interp_raw,
    struct interpolation_1d *interp_int, const struct stellar_model *sm,
    int interpolation_size, const struct unit_system *us,
    const struct phys_const* phys_const) {

  /* Allocate the memory */
  const int count = 500;
  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for luminosities.");

  const float mass_min = sm->imf.mass_min;
  const float mass_max = sm->imf.mass_max;
  const float log_mass_min = log10f(mass_min);
  const float log_mass_max = log10f(mass_max);
  const float step_size = (log_mass_max - log_mass_min) / (count - 1);

  /* Fill the table */
  for (size_t j = 0; j < count; j++) {
    /* Compute the log-mass and mass */
    const float log_mass = log_mass_min + j * step_size;
    const float mass = exp10(log_mass) * phys_const->const_solar_mass;

    /* Get bolometric luminosity for this mass, in internal units */
    data[j] = radiation_get_individual_star_luminosity(mass, us, phys_const);
  }

  /* Initialize the raw interpolation */
  interpolate_1d_init(interp_raw, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_error);

  initial_mass_function_integrate(&sm->imf, data, count, log_mass_min,
                                  step_size);
  // TODO: decrease count in order to keep the same distance between points

  /* Initialize the integrated interpolation */
  interpolate_1d_init(interp_int, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_const);

  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read an array of ionizing emission rates data from the table.
 *
 * @param rad The #radiation model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm * The #stellar_model.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size Number of element to keep in the interpolation
 * data.
 */
void radiation_read_ionization_rate_array(
    struct radiation *rad, struct interpolation_1d_double *interp_raw,
    struct interpolation_1d_double *interp_int, const struct stellar_model *sm,
    int interpolation_size, const struct unit_system *us, const struct phys_const* phys_const) {

  /* Allocate the memory */
  const int count = 500;
  double *data = (double *)malloc(sizeof(double) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for luminosities.");

  const float mass_min = sm->imf.mass_min;
  const float mass_max = sm->imf.mass_max;
  const float log_mass_min = log10f(mass_min);
  const float log_mass_max = log10f(mass_max);
  const float step_size = (log_mass_max - log_mass_min) / (count - 1);

  /* Fill the table */
  for (size_t j = 0; j < count; j++) {
    /* Compute the log-mass and mass */
    const float log_mass = log_mass_min + j * step_size;
    const float mass = exp10(log_mass) * phys_const->const_solar_mass;

    /* Get bolometric luminosity for this mass, in internal units */
    data[j] = radiation_get_individual_star_ionizing_photon_emission_rate_fit(mass, us, phys_const);
  }

  /* Initialize the raw interpolation */
  interpolate_1d_double_init(interp_raw, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_error);

  initial_mass_function_integrate_double(&sm->imf, data, count, log_mass_min,
                                  step_size);
  // TODO: decrease count in order to keep the same distance between points

  /* Initialize the integrated interpolation */
  interpolate_1d_double_init(interp_int, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_const);

  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read the RAD yields from the table.
 *
 * The tables are in internal units at the end of this function.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void radiation_read_data(struct radiation *rad, struct swift_params *params,
			 const struct stellar_model *sm,
			 const struct unit_system *us,
			 const struct phys_const* phys_const,
			 const int restart) {

  if (!restart) {
    /* TODO: Maybe update this */
    rad->interpolation_size = parser_get_opt_param_int(
        params, "GEARSupernovaeII:interpolation_size", 200);
  }

  /* Read the luminosities */
  radiation_read_luminosities_array(rad, &rad->raw.luminosities,
				    &rad->integrated.luminosities, sm,
				    rad->interpolation_size, us, phys_const);

  /* Read the ionization emission rates */
  radiation_read_ionization_rate_array(rad, &rad->raw.dot_N_ion,
				       &rad->integrated.dot_N_ion,
				       sm, rad->interpolation_size, us, phys_const);
};
