/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IONIZATION_EQUILIBRIUM_H
#define SWIFT_RT_IONIZATION_EQUILIBRIUM_H

#include "rt_thermochemistry_utils.h"

/* some local definitions */
#define RT_ION_EQUIL_ITER_MAX 50
#define RT_ION_EQUIL_EPSILON 1e-4

/**
 * @file src/rt/GEAR/rt_ionization_equilibrium.h
 * @brief Function to get the number densities for primordial species as a
 * function of the temperature, assuming ionization equilibrium (see Katz et
 * al. 1996; ui.adsabs.harvard.edu/abs/1996ApJS..105...19K)
 */

/**
 * @brief compute the ionization equilibrium mass fractions
 * for a given temperature
 * @param T temperature in K
 * @param X total mass fraction of all hydrogen species
 * @param Y total mass fraction of all helium species
 * @param XHI (return) mass fraction of HI
 * @param XHII (return) mass fraction of HII
 * @param XHeI (return) mass fraction of HeI
 * @param XHeII (return) mass fraction of HeII
 * @param XHeIII (return) mass fraction of HeII
 */
__attribute__((always_inline)) INLINE static void
rt_ion_equil_mass_fractions_from_T(double T, float X, float Y, float* XHI,
                                   float* XHII, float* XHeI, float* XHeII,
                                   float* XHeIII) {

  if (fabsf(X + Y - 1.f) > 1e-4)
    error("mass fractions don't add up to one: X=%.3e, Y=%.3e", X, Y);

  /* Total number densities for all H and He species, respectively.
   * We don't really care about the actual number densities of the
   * species, but only want the mass fractions. So all quantities
   * will be per unit volume, specifically per rho_gas / m_proton. */
  double nH, nHe;

  /* TODO: this assumes no metals. */
  if (X <= RT_GEAR_TINY_MASS_FRACTION) {
    nH = RT_GEAR_TINY_MASS_FRACTION;
    nHe = 1.;
  } else {
    nH = X;
    nHe = 0.25 * Y;
  }

  /* Number densities of the actual species */
  float nHI, nHII, nHeI, nHeII, nHeIII;

  if (T < 5e3) {
    /* Below 5000K, we just go with fully neutral gas.
     * Not only is this the case, but some divisions
     * would otherwise not be safe and lead to problems
     * if we don't exception handle it here. */
    nHI = nH;
    nHII = 0.f;
    nHeI = nHe;
    nHeII = 0.f;
    nHeIII = 0.f;

  } else {

    const double T_inv = 1. / T;
    const double sqrtT = sqrt(T);
    const double temp_pow07 = 1. / (1. + pow(T * 1e-6, 0.7));
    const double temp_sqrtT5 = 1. / (1. + sqrt(T * 1e-5));
    const double temp_pow02 = pow(T * 1e-3, -0.2);
    /* Recombination rate for H+ in units of cm^3 s^-1  */
    const double A_Hp = 8.40e-11 / sqrtT * temp_pow02 * temp_pow07;
    /* Dielectronic recombination rate for He+ in units of cm^3 s^-1 */
    const double A_d = 1.9e-3 / pow(T, 1.5) * exp(-470000.0 * T_inv) *
                       (1.0 + 0.3 * exp(-94000.0 * T_inv));
    /* Recombination rate for He+ in units of cm^3 s^-1 */
    const double A_Hep = 1.5e-10 / pow(T, 0.6353);
    /* Recombination rate for He++ in units of cm^3 s^-1 */
    const double A_Hepp = 3.36e-10 / sqrtT * temp_pow02 * temp_pow07;
    /* collisional ionization rate for H0 in units of cm^3 s^-1 */
    const double G_H0 = 5.85e-11 * sqrtT * exp(-157809.1 * T_inv) * temp_sqrtT5;
    /* collisional ionization rate for He0 in units of cm^3 s^-1 */
    const double G_He0 =
        2.38e-11 * sqrtT * exp(-285335.4 * T_inv) * temp_sqrtT5;
    /* collisional ionization rate for He+ in units of cm^3 s^-1 */
    const double G_Hep =
        5.68e-12 * sqrtT * exp(-631515.0 * T_inv) * temp_sqrtT5;

    /* Katz et al. 1996 eq. 33 - 38 */
    /* Note: We assume all photoionization rates to be zero. */
    nHI = nH * A_Hp / (A_Hp + G_H0);
    nHII = nH - nHI;
    const double div1 = (A_Hep + A_d) / G_He0;
    const double div2 = G_Hep / A_Hepp;
    nHeII = nHe / (1. + div1 + div2);
    nHeI = nHeII * div1;
    nHeIII = nHeII * div2;
    /* ne = nHII + nHeII + 2.f * nHeIII; */
  }

  /* With the number densities per unit volume given,
   * let's compute the mass fractions now. */
  const double mHI = nHI;
  const double mHII = nHII;
  const double mHeI = 4.0f * nHeI;
  const double mHeII = 4.0f * nHeII;
  const double mHeIII = 4.0f * nHeIII;
  /* we ignore the electron mass fraction. */

  const double mtot_inv = 1. / (mHI + mHII + mHeI + mHeII + mHeIII);

  *XHI = (float)(mHI * mtot_inv);
  *XHI = max(*XHI, RT_GEAR_TINY_MASS_FRACTION);
  *XHII = (float)(mHII * mtot_inv);
  *XHII = max(*XHII, RT_GEAR_TINY_MASS_FRACTION);
  *XHeI = (float)(mHeI * mtot_inv);
  *XHeI = max(*XHeI, RT_GEAR_TINY_MASS_FRACTION);
  *XHeII = (float)(mHeII * mtot_inv);
  *XHeII = max(*XHeII, RT_GEAR_TINY_MASS_FRACTION);
  *XHeIII = (float)(mHeIII * mtot_inv);
  *XHeIII = max(*XHeIII, RT_GEAR_TINY_MASS_FRACTION);
}

/**
 * @brief compute the ionization equilibrium mass fractions
 * for a given particle.
 * @param XHI (return) mass fraction of HI
 * @param XHII (return) mass fraction of HII
 * @param XHeI (return) mass fraction of HeI
 * @param XHeII (return) mass fraction of HeII
 * @param XHeIII (return) mass fraction of HeII
 * @param p part to work with
 * @param rt_props rt_properties struct
 * @param hydro_props hydro properties struct
 * @param phys_const physical constants struct
 * @param us unit system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void
rt_ion_equil_get_mass_fractions(float* XHI, float* XHII, float* XHeI,
                                float* XHeII, float* XHeIII,
                                struct part* restrict p,
                                const struct rt_props* rt_props,
                                const struct hydro_props* hydro_props,
                                const struct phys_const* restrict phys_const,
                                const struct unit_system* restrict us,
                                const struct cosmology* restrict cosmo) {

  /* get conversions and constants */
  const double internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  const float dimension_kB[5] = {1, 2, -2, 0, -1}; /* [g cm^2 s^-2 K^-1] */
  const double kB_to_cgs =
      units_general_cgs_conversion_factor(us, dimension_kB);
  const double kB_cgs = phys_const->const_boltzmann_k * kB_to_cgs;
  const double mp_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double mp_cgs = phys_const->const_proton_mass * mp_to_cgs;

  /* Get the specific internal energy of the gas */
  const float u_minimal = hydro_props->minimal_internal_energy;
  /* Using 'drifted' version here because I'm lazy and don't want to pass
   * the xpart down to use in this function. */
  const float u_part = hydro_get_drifted_physical_internal_energy(p, cosmo);
  const double u_expect =
      ((double)max(u_minimal, u_part)) * internal_energy_to_cgs;
  double mu_guess, T_guess;

  /* Get a first estimate for gas temperature. */
  const float X = rt_props->hydrogen_mass_fraction;
  const float Y = rt_props->helium_mass_fraction;
  *XHI = X;
  *XHII = 0.f;
  *XHeI = Y;
  *XHeII = 0.f;
  *XHeIII = 0.f;
  mu_guess =
      rt_tchem_get_mean_molecular_weight(*XHI, *XHII, *XHeI, *XHeII, *XHeIII);
  T_guess = rt_tchem_temperature_from_internal_energy(u_expect, mu_guess,
                                                      kB_cgs, mp_cgs);

  if (T_guess > 1e4) {
    /* T_guess should be in K, so it's fine to hardcode
     * the ionization limit of ~1e4 K for the first guess here
     * If we're above the temperature threshold with this guess,
     * assume we're fully ionized as first guess instead. */
    *XHI = 0.f;
    *XHII = X;
    *XHeI = 0.f;
    *XHeII = 0.f;
    *XHeIII = Y;
    mu_guess =
        rt_tchem_get_mean_molecular_weight(*XHI, *XHII, *XHeI, *XHeII, *XHeIII);
    T_guess = rt_tchem_temperature_from_internal_energy(u_expect, mu_guess,
                                                        kB_cgs, mp_cgs);
  }

  /* Now given the first temperature guess, update
   * the mass fractions and mean molecular weight */
  rt_ion_equil_mass_fractions_from_T(T_guess, X, Y, XHI, XHII, XHeI, XHeII,
                                     XHeIII);
  mu_guess =
      rt_tchem_get_mean_molecular_weight(*XHI, *XHII, *XHeI, *XHeII, *XHeIII);

  /* Get first guess for internal energy */
  double u_guess =
      rt_tchem_internal_energy_from_T(T_guess, mu_guess, kB_cgs, mp_cgs);

  int iter = 0;
  double du = u_expect - u_guess;
  double du_old = du; /* initialize as same value */

  while (fabs(du) >= RT_ION_EQUIL_EPSILON * u_expect) {
    iter += 1;
    if (iter > RT_ION_EQUIL_ITER_MAX) {
      message(
          "Warning: Ionization Equilibrium iteration didn't converge; "
          "T=%.6g, 1 - u/u_correct = %.6g",
          T_guess, 1. - u_guess / u_expect);
      break;
    }

    if (T_guess < 0.) {
      message("Warning: Got negative temperature, resetting");
      T_guess = 0.1; /* Note: T_guess should be in K here */
    }

    /* find next temperature guess by solving linear equation
     * m * T_next + n = u_expect - u_guess_new ~ 0
     * NOTE: we pretend that the function that we're looking the
     * root of is f(T) = u_expect - u_guess(T), with u_expect = const.
     * so df/dT = - du_guess/dT, therefore we add a minus sign here. */
    double m = -rt_tchem_internal_energy_dT(mu_guess, kB_cgs, mp_cgs);
    double n = u_expect - u_guess - m * T_guess;
    double T_next = -n / m;

    /* Given the new temperature guess, compute the
     * expected mean molecular weight */
    rt_ion_equil_mass_fractions_from_T(T_next, X, Y, XHI, XHII, XHeI, XHeII,
                                       XHeIII);
    double mu_next =
        rt_tchem_get_mean_molecular_weight(*XHI, *XHII, *XHeI, *XHeII, *XHeIII);

    /* now given the new temperature and mass fraction guess,
     * update the expected gas internal energy */
    double u_next =
        rt_tchem_internal_energy_from_T(T_next, mu_next, kB_cgs, mp_cgs);

    /* save the old internal energy, and get the new one */
    du_old = du;
    du = u_expect - u_next;

    /* if we're oscillating between positive and negative
     * values, try a bisection to help out */
    if (du_old * du < 0.0) {
      T_next = 0.5 * (T_guess + T_next);
      rt_ion_equil_mass_fractions_from_T(T_next, X, Y, XHI, XHII, XHeI, XHeII,
                                         XHeIII);
      mu_next = rt_tchem_get_mean_molecular_weight(*XHI, *XHII, *XHeI, *XHeII,
                                                   *XHeIII);
      u_next = rt_tchem_internal_energy_from_T(T_next, mu_next, kB_cgs, mp_cgs);
    }

    /* prepare for next iteration */
    T_guess = T_next;
    mu_guess = mu_next;
    u_guess = u_next;
    du = u_expect - u_next;
  }
}

#undef RT_ION_EQUIL_ITER_MAX
#undef RT_ION_EQUIL_EPSILON

#endif /* SWIFT_RT_IONIZATION_EQUILIBRIUM_H */
