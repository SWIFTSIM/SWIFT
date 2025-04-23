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
#ifndef SWIFT_RADIATION_GEAR_H
#define SWIFT_RADIATION_GEAR_H

#include "cooling.h"
#include "hdf5_functions.h"
#include "hydro.h"
#include "interpolation.h"
#include "part.h"
#include "physical_constants.h"
#include "stellar_evolution_struct.h"
#include "units.h"

#define BB_NU_INTEGRATION_STEPS 500
#define RADIATION_N_IONIZATION_STEPS 1000

/**
 * @brief Return the specific intensity of the blackbody spectrum
 *
 * @param nu frequency at which to compute specific intensity
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 * @param c speed of light
 */
__attribute__((always_inline)) INLINE double
radiation_blackbody_spectrum_intensity(const double nu, const double T,
                                       const double kB, const double h_planck,
                                       const double c) {

  const double hnu = h_planck * nu;
  const double kT = kB * T;
  const double nu2 = nu * nu;
  double temp;
  if (hnu / kT < 1e-6) {
    /* prevent division by zero, use Taylor approximation */
    temp = kT;
  } else if (hnu / kT > 700.) {
    /* prevent infs */
    temp = 0.;
  } else {
    temp = 1. / (exp(hnu / kT) - 1.);
  }
  return 2. * hnu * nu2 / (c * c) * temp;
}

/**
 * Return the blackbody spectrum energy density
 *
 * @param nu frequency at which to compute specific intensity
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 * @param c speed of light
 */
__attribute__((always_inline)) INLINE double
radiation_blackbody_spectrum_energy_density(const double nu, const double T,
                                            const double kB,
                                            const double h_planck,
                                            const double c) {
  return 4. * M_PI / c *
         radiation_blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
}

float radiation_get_blackbody_luminosity_band(const float nu_min,
                                              const float nu_max, const float T,
                                              const float R, const float kB,
                                              const float h, const float c);

float radiation_get_ionizing_photon_emission_rate(const float nu_min,
                                                  const float nu_max,
                                                  const float T, const float R,
                                                  const float kB, const float h,
                                                  const float c);

double radiation_get_star_ionization_rate(const struct spart* sp);

double radiation_get_part_number_hydrogen_atoms(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp);

double radiation_get_part_rate_to_fully_ionize(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp);

void radiation_tag_part_as_ionized(struct part* p, struct xpart* xpj);
void radiation_reset_part_ionized_tag(struct part* p, struct xpart* xpj);
int radiation_is_part_tagged_as_ionized(struct part* p, struct xpart* xpj);
void radiation_consume_ionizing_photons(struct spart* sp,
                                        double Delta_dot_N_ion);
float radiation_get_star_gas_column_density(const struct spart* sp);

float radiation_get_IR_opacity(const struct spart* sp,
                               const struct unit_system* us,
                               const struct phys_const* phys_const);

float radiation_get_IR_optical_depth(const struct spart* sp,
                                     const struct unit_system* us,
                                     const struct phys_const* phys_const);

float radiation_get_star_radiation_pressure(
    const struct spart* sp, const float Delta_t, const struct unit_system* us,
    const struct phys_const* phys_const);

int radiation_is_part_ionized(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp);

float radiation_get_individual_star_radius(const struct spart* sp,
                                           const struct unit_system* us,
                                           const struct phys_const* phys_const);
float radiation_get_individual_star_temperature(
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const);

float radiation_get_individual_star_luminosity(
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const);

double radiation_get_individual_star_ionizing_photon_emission_rate_fit(
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const);

#endif /* SWIFT_RADIATION_GEAR_H */
