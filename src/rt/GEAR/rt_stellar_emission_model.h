/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STELLAR_EMISSION_MODEL_GEAR_H
#define SWIFT_RT_STELLAR_EMISSION_MODEL_GEAR_H

/**
 * @file src/rt/GEAR/rt_stellar_emission_model.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * stellar radiation emission models.
 */

enum rt_stellar_emission_models {
  rt_stellar_emission_model_none = 0,
  rt_stellar_emission_model_const,
  rt_stellar_emission_model_IlievTest,
  rt_stellar_emission_model_count
};

/**
 * @brief Compute the energy emitted from a star during the time step dt.
 * This is for the constant emission rate model.
 *
 * @param emission_this_step (return) the emitted radiation energy of a star
 * during the time interval dt
 * @param const_stellar_emission_rates the constant emission rates used in this
 * run
 * @param dt time step size (in internal units)
 */
__attribute__((always_inline)) INLINE static void
rt_get_emission_this_step_const(
    double emission_this_step[RT_NGROUPS],
    const double const_stellar_emission_rates[RT_NGROUPS], double dt) {

  /* The read-in constant stellar emisison rates are in units of L_sol.
   * But they have been read in assuming they are in cgs. Convert this
   * only now to proper internal units to avoid float overflows. We only
   * store the energy that is to be distributed from this spart to its
   * neighbours in this step in internal units.*/
  const double solar_luminosity = 3.828e33; /* erg/s */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const double emission_rate_internal_units =
        const_stellar_emission_rates[g] * solar_luminosity;
    emission_this_step[g] = emission_rate_internal_units * dt;
  }
}

/**
 * @brief Compute the energy emitted from a star during the time step dt.
 * This is for the Iliev+2006 Test 4.
 *
 * @param emission_this_step (return) the emitted radiation energy of a star
 * during the time interval dt
 * @param M the star mass (in internal units)
 * @param dt time step size (in internal units)
 * @param photon_number_integral Integrated photon numbers over frequency
 * interval
 * @param average_photon_energy average photon energy in each frequency bin, in
 * erg
 * @param phys_const struct holding physical constants
 * @param internal_units units struct containing internal units
 */
__attribute__((always_inline)) INLINE static void
rt_get_emission_this_step_IlievTest(
    double emission_this_step[RT_NGROUPS], float M, const double dt,
    const double photon_number_integral[RT_NGROUPS],
    const double average_photon_energy[RT_NGROUPS],
    const struct phys_const* phys_const,
    const struct unit_system* internal_units) {

  /* Note that this model uses the halo mass to determine the luminosity
   * of a source. I'm cheating the system here by storing the required halo
   * mass as the star mass. This is only ok because the test is supposed to
   * run with all dynamics and gravity turned off. */

  const double Omega_b = 0.043;
  const double Omega_0 = 0.27;
  const double m_p = phys_const->const_proton_mass;
  const double t_s = 3e6 * phys_const->const_year;
  const double f_gamma = 250.;
  const double Ndot_gamma = (f_gamma * M * Omega_b) / (Omega_0 * m_p * t_s);

  double Nsum = 0.;
  for (int g = 0; g < RT_NGROUPS; g++) Nsum += photon_number_integral[g];

  if (Nsum <= 0.) error("No photons in spectrum...???");

  const double energy_units =
      units_cgs_conversion_factor(internal_units, UNIT_CONV_ENERGY);
  for (int g = 0; g < RT_NGROUPS; g++) {
    const double fi = photon_number_integral[g] / Nsum;
    const double Ndot_i = fi * Ndot_gamma;
    /* average photon densities are in cgs! */
    const double Edot_i = average_photon_energy[g] * Ndot_i / energy_units;
    emission_this_step[g] = Edot_i * dt;
  }
}

#endif /* SWIFT_RT_STELLAR_EMISSION_MODEL_GEAR_H */
