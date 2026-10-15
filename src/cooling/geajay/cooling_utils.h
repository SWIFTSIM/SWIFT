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
#ifndef SWIFT_COOLING_GEAJAY_COOLING_UTILS_H
#define SWIFT_COOLING_GEAJAY_COOLING_UTILS_H
/**
 * @file src/cooling/geajay/cooling_utils.h
 * @brief Cooling utilities functions for geajay.
 */

#include "chemistry.h"
#include "cooling_properties.h"
#include "error.h"
#include "hydro.h"
#include "units.h"

#include <math.h>

/**
 * Compute the Hydrogen mass fraction.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static float
cooling_get_hydrogen_mass_fraction(const struct cooling_function_data *cooling,
                                   const struct part *p,
                                   const struct xpart *xp) {

#if COOLING_GEAJAY_MODE == 0
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Approximation that takes into account metals */
  return cooling->HydrogenFractionByMass - Z;
#else

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const float X_H = cool_data->HI_frac + cool_data->HII_frac;
  return X_H;
#endif
}

/**
 * Compute the Helium mass fraction.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static float
cooling_get_helium_mass_fraction(const struct cooling_function_data *cooling,
                                 const struct part *p, const struct xpart *xp) {

#if COOLING_GEAJAY_MODE == 0
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Approximation that takes into account metals */
  const float X_H = cooling->HydrogenFractionByMass - Z;
  return 1 - X_H;
#else

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const float X_He =
      cool_data->HeI_frac + cool_data->HeII_frac + cool_data->HeIII_frac;
  return X_He;
#endif
}

/**
 * @brief Compute the mean molecular weight of a particle, as geajay uses it.
 *
 * Calls geajay_particle_mean_molecular_weight() at the drifted state. Stops the
 * run if geajay refuses the particle.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Mean molecular weight.
 */
double cooling_get_mean_molecular_weight(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

/**
 * @brief Compute the specific internal energy at which a particle has the
 * temperature T_target: the inverse of cooling_get_temperature().
 *
 * Calls geajay_particle_internal_energy_from_temperature(), which uses the
 * species of the particle. Stops the run if geajay refuses the particle or
 * returns an energy that is not positive.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param T_target Desired temperature.
 * @return Specific internal energy that gives T_target under
 *   cooling_get_temperature().
 */
double cooling_get_internal_energy_from_temperature(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp, const double T_target);

/**
 * @brief Compute the mean molecular weight of gas in ionisation equilibrium
 * from its density, energy and metal content.
 *
 * Calls geajay_particle_equilibrium_mean_molecular_weight(). Use it where the
 * species of the gas are not known.
 *
 * @param u Internal energy in physical units.
 * @param rho_physical Physical gas density.
 * @param metal_mass_fraction Metal mass fraction of the gas.
 * @param phys_const Physical constants.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @return Mean molecular weight.
 */
double cooling_get_equilibrium_mean_molecular_weight(
    const float u, const double rho_physical, const double metal_mass_fraction,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling);

/**
 * @brief compute the temperature of an ideal gas for a given specific internal
 * energy and mean molecular weight.
 *
 * @param u Specific internal energy of the gas.
 * @param mu Mean molecular weight of the gas.
 * @param kB Boltzmann constant.
 * @param mp Proton mass.
 * */
__attribute__((always_inline)) INLINE static float
cooling_temperature_from_internal_energy(const double u, const double mu,
                                         const double kB, const double mp) {
  return u * hydro_gamma_minus_one * mu * mp / kB;
}

/**
 * @brief Compute the temperature of a particle at a given energy.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param u Physical specific internal energy.
 * @return The temperature in K.
 */
double cooling_get_temperature_at_energy(
    const struct cosmology *cosmo, const struct cooling_function_data *cooling,
    const struct part *p, const struct xpart *xp, const double u);

/**
 * @brief compute the AGORA redshift-dependent CMB-floor specific internal
 * energy for a particle, from its current composition.
 *
 * Call again after the composition changes (e.g. a Grackle chemistry solve)
 * to keep it consistent with cooling_get_temperature()'s later mu.
 *
 * The inverse is costly in mode 0, so it runs only if the particle at u_now is
 * colder than the floor. Otherwise the floor does not bind at u_now and the
 * function returns 0.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param u_now Physical specific internal energy that the floor is applied to.
 * @return Specific internal energy corresponding to T_CMB,0 * (1+z), or 0 if
 *   the temperature at u_now is not below it.
 */
__attribute__((always_inline)) INLINE static double
cooling_agora_cmb_floor_internal_energy(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp, const double u_now) {

  const double z = (cooling->redshift == -1) ? cosmo->z : cooling->redshift;
  const double T_CMB_agora = phys_const->const_T_CMB_0 * (z + 1.0);

  /* geajay takes no non-positive energy: that case goes to the inverse. */
  if (u_now > 0. && cooling_get_temperature_at_energy(cosmo, cooling, p, xp,
                                                      u_now) >= T_CMB_agora)
    return 0.;

  return cooling_get_internal_energy_from_temperature(
      phys_const, us, cosmo, hydro_props, cooling, p, xp, T_CMB_agora);
}
#endif /* SWIFT_COOLING_GEAJAY_COOLING_UTILS_H */
