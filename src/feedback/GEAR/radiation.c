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

/* Include header */
#include "radiation.h"


/* TODO: Check unit... similaru in radiation_blackbody_etc */
float radiation_get_blackbody_luminosity_band(
    const float nu_min, const float nu_max, const float T, const float R,
    const float kB, const float h, const float c) {

  const float dnu = (nu_max - nu_min) / BB_NU_INTEGRATION_STEPS;
  float luminosity = 0.f;

  for (int i = 0; i < BB_NU_INTEGRATION_STEPS; ++i) {
    const float nu = nu_min + (i + 0.5f) * dnu;
    const float intensity = radiation_blackbody_spectrum_intensity(nu, T, kB, h, c);

    // Calculate the luminosity contribution from the band
    const float dL = 4.f * M_PI * R * R * intensity * dnu;

    luminosity += dL;
  }
  return luminosity;
}

/* TODO: Check unit... similaru in radiation_blackbody_etc */
float radiation_get_ionizing_photon_emission_rate(const float nu_min, const float nu_max,
                                  const float T, const float R,
                                  const float kB, const float h, const float c) {

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

  return surface_area * integral; // [photons / second]
}


float radiation_get_star_ionisation_rate(const struct spart* sp) {
  return sp->feedback_data.radiation.dot_N_ion ;
}

float radiation_get_part_rate_to_fully_ionize(const struct part* p, const struct xpart* xp) {

  return 0.0;
}

void radiation_consume_ionizing_photons(struct spart* sp, float Delta_dot_N_ion) {
  sp->feedback_data.radiation.dot_N_ion -= Delta_dot_N_ion;
  return;
}

void radiation_tag_part_as_ionized(struct part* p, struct xpart* xpj) {
  xpj->feedback_data.radiation.is_ionized = 1;
  return;
}

int radiation_is_part_ionized(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp) {

  /* Is T > 10^4 K ? */
  const float T =  cooling_get_temperature(phys_const, hydro_props, us,
					   cosmo, cooling, p, xp);
  const float ten_to_four_kelvin = 1e4 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Is the particle ionized ? */
  return (T > ten_to_four_kelvin || xp->feedback_data.radiation.is_ionized);
}

float radiation_get_individual_star_radius(const struct spart* sp,
					   const struct unit_system* us,
					   const struct phys_const* phys_const) {

  const float R_sun = phys_const->const_solar_radius; /* In internal units */
  const float M_solar = phys_const->const_solar_mass;

  const float M = sp->mass; /* In internal units */
  const float M_in_solar = M / M_solar; /* In solar masses */

  /* TODO: Correct unit conversion */
  if (M_in_solar < 1.f) {
    return R_sun * powf(M, 0.8f);
  } else if (M_in_solar < 8.f) {
    return R_sun * powf(M, 0.57f);
  } else {
    return R_sun * powf(M, 0.5f);
  }
}

float radiation_get_individual_star_temperature(const struct spart* sp,
						const struct unit_system* us,
						const struct phys_const* phys_const) {

  const float M_solar = phys_const->const_solar_mass;
  const float M = sp->mass; /* In internal units */
  const float M_in_solar = M / M_solar; /* In solar masses */

  float T_K = 0.0;

  if (M_in_solar < 1.f) {
    T_K = 3500.f * powf(M_in_solar, 0.5f);
  } else if (M_in_solar < 8.f) {
    T_K = 5800.f * powf(M_in_solar, 0.5f);
  } else {
    T_K = 25000.f * powf(M_in_solar / 20.f, 0.1f);
  }

  /* Convert from Kelvin to internal units using unit_system_temperature_in_cgs */
  const float T_internal = T_K /  units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  return T_internal;
}

/**
 * Return the bolometric luminosity of a single star from empirical
 * mass-luminosity relations.
 *
 * Uses a piecewise power-law approximation:
 * - L ∝ M^2 for M < 0.43 Msun
 * - L ∝ M^4 for 0.43 ≤ M < 2.0 Msun
 * - L ∝ M^3.5 for 2.0 ≤ M < 54 Msun
 * - L ∝ 32000 M for M ≥ 54 Msun
 *
 * Returns luminosity in code units.
 *
 * @param sp Pointer to star particle.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Luminosity in code units.
 */
float radiation_get_individual_star_luminosity(
					       const struct spart* sp,
					       const struct unit_system* us,
					       const struct phys_const* phys_const) {

  /* Convert mass to solar masses */
  const float M_in_solar = sp->mass / phys_const->const_solar_mass;

  /* Piecewise empirical mass-luminosity relation */
  float lum_sol;
  if (M_in_solar < 0.43f) {
    lum_sol = 0.185f *  M_in_solar *  M_in_solar;
  } else if ( M_in_solar < 2.0f) {
    lum_sol =  M_in_solar *  M_in_solar *  M_in_solar *  M_in_solar;
  } else if ( M_in_solar < 54.0f) {
    lum_sol = 1.5f *  M_in_solar *  M_in_solar *  M_in_solar * sqrtf( M_in_solar);
  } else {
    lum_sol = 32000.0f *  M_in_solar;
  }

  /* Convert from solar luminosities to code units */
  const float luminosity = lum_sol * phys_const->const_solar_luminosity;
  return luminosity;
}
