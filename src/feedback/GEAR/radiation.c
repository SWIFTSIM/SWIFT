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

#include "kernel_hydro.h"
#include "units.h"

/* TODO: Check unit... similaru in radiation_blackbody_etc */
float radiation_get_blackbody_luminosity_band(const float nu_min,
                                              const float nu_max, const float T,
                                              const float R, const float kB,
                                              const float h, const float c) {

  const float dnu = (nu_max - nu_min) / BB_NU_INTEGRATION_STEPS;
  float luminosity = 0.f;

  for (int i = 0; i < BB_NU_INTEGRATION_STEPS; ++i) {
    const float nu = nu_min + (i + 0.5f) * dnu;
    const float intensity =
        radiation_blackbody_spectrum_intensity(nu, T, kB, h, c);

    // Calculate the luminosity contribution from the band
    const float dL = 4.f * M_PI * R * R * intensity * dnu;

    luminosity += dL;
  }
  return luminosity;
}

/* TODO: Check unit... similaru in radiation_blackbody_etc */
float radiation_get_ionizing_photon_emission_rate(const float nu_min,
                                                  const float nu_max,
                                                  const float T, const float R,
                                                  const float kB, const float h,
                                                  const float c) {

  const float dnu = (nu_max - nu_min) / RADIATION_N_IONIZATION_STEPS;
  float integral = 0.f;

  for (int i = 0; i < RADIATION_N_IONIZATION_STEPS; i++) {
    const float nu1 = nu_min + i * dnu;
    const float nu2 = nu1 + dnu;

    const float B1 = radiation_blackbody_spectrum_intensity(nu1, T, kB, h, c);
    const float B2 = radiation_blackbody_spectrum_intensity(nu2, T, kB, h, c);

    const float integrand1 = B1 / (h * nu1);
    const float integrand2 = B2 / (h * nu2);

    integral += 0.5f * (integrand1 + integrand2) * dnu;
  }

  const float surface_area = 4.f * M_PI * R * R;

  return surface_area * integral;  // [photons / second]
}

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
double radiation_get_star_ionization_rate(const struct spart* sp) {
  return sp->feedback_data.radiation.dot_N_ion;
}

/**
 * Consume the #spart ionizing photon budget.
 *
 * @param sp The star.
 * @param Delta_dot_N_ion The ionizing photon rate to remove.
 */
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
int radiation_is_part_tagged_as_ionized(struct part* p, struct xpart* xp) {
  return xp->feedback_data.radiation.is_ionized;
}

/**
 * Compute the temperature of a single star from empirical mass-temperature
 * relations.
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

// TODO: DO converstion to physical? Or to comoving?
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

float radiation_get_physical_IR_optical_depth(const struct spart* sp,
					      const struct unit_system* us,
					      const struct phys_const* phys_const,
					      const struct cosmology* cosmo) {
  const float Sigma_gas_c = radiation_get_comoving_gas_column_density_at_star(sp);
  const float Sigma_gas_p = Sigma_gas_c * cosmo->a2_inv;
  const float kappa_IR = radiation_get_physical_IR_opacity(sp, us, phys_const, cosmo);
  return kappa_IR * Sigma_gas_p;
}

float radiation_get_star_physical_radiation_pressure(
    const struct spart* sp, const float Delta_t, const struct unit_system* us,
    const struct phys_const* phys_const, const struct cosmology* cosmo) {

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
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const) {

  /* Perform some units conversions */
  const float R_sun = phys_const->const_solar_radius;
  const float M_solar = phys_const->const_solar_mass;
  const float M = sp->mass;
  const float M_in_solar = M / M_solar;

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
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const) {

  const float M_solar = phys_const->const_solar_mass;
  const float M = sp->mass;             /* In internal units */
  const float M_in_solar = M / M_solar; /* In solar masses */

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
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const) {

  /* Convert mass to solar masses */
  const float M_in_solar = sp->mass / phys_const->const_solar_mass;

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
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const) {

  /* Get star properties in internal units */
  /* const float T = radiation_get_individual_star_temperature(sp, us,
   * phys_const); */
  const float R = radiation_get_individual_star_radius(sp, us, phys_const);
  const float L = radiation_get_individual_star_luminosity(sp, us, phys_const);

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
