/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_RATES_TREECOOL_H
#define SWIFT_COOLING_RATES_TREECOOL_H

/**
 * @file src/cooling/TREECOOL/cooling_rates.h
 * @brief Rate coefficients, ionization equilibrium and cooling rates of the
 * TREECOOL cooling function, i.e. the primordial cooling function of Katz,
 * Weinberg & Hernquist (1996).
 *
 * All the rate coefficients are the fits collected in Table 2 of
 * Katz, Weinberg & Hernquist (1996), ApJS, 105, 19 (hereafter KWH96), which
 * are themselves mostly taken from Cen (1992), ApJS, 78, 341. The cooling and
 * heating terms are the ones listed in Table 1 of the same paper.
 *
 * The gas is assumed to be primordial (Hydrogen and Helium only) and to be in
 * ionization equilibrium with the collisional processes and with an optically
 * thin, uniform UV background (KWH96, eq. 33-38).
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "adiabatic_index.h"
#include "cooling_properties.h"
#include "error.h"
#include "exp10.h"
#include "inline.h"
#include "minmax.h"

/*! Maximal number of iterations of the ionization equilibrium solver. */
#define treecool_max_iterations 150

/*! Number of iterations of the damped fixed-point temperature solver after
 * which we give up and fall back to bisection. The scheme converges in fewer
 * than 15 iterations in the overwhelming majority of cases. */
#define treecool_temperature_max_fixed_point_iterations 30

/*! Maximal number of iterations of the bisection fallback of the temperature
 * solver. The bracket only spans the ratio of the neutral to the fully ionized
 * mean molecular weight, so ~12 iterations are enough in practice. */
#define treecool_temperature_max_bisection_iterations 100

/*! Convergence criterion of the ionization equilibrium solver, expressed as an
 * absolute change in the electron fraction n_e / n_H. */
#define treecool_ionization_tolerance 1.e-4

/*! Convergence criterion of the temperature solver, expressed as a relative
 * change in temperature. */
#define treecool_temperature_tolerance 1.e-3

/*! Number below which an electron density is considered to be zero [cm^-3] */
#define treecool_min_electron_density 1.e-25

/*! Generic tiny number used to protect divisions */
#define treecool_small_number 1.e-60

/*! Largest exponent we evaluate in the Boltzmann factors. Beyond this value
 * the corresponding rate is set to zero. */
#define treecool_max_exponent 70.

/**
 * @brief The state of the gas at a given temperature and density.
 *
 * All the number densities are expressed in units of the Hydrogen number
 * density n_H. All the rate coefficients are in physical cgs units.
 */
struct treecool_gas_state {

  /*! Number density of HI (in units of n_H) */
  double n_H0;

  /*! Number density of HII (in units of n_H) */
  double n_Hp;

  /*! Number density of HeI (in units of n_H) */
  double n_He0;

  /*! Number density of HeII (in units of n_H) */
  double n_Hep;

  /*! Number density of HeIII (in units of n_H) */
  double n_Hepp;

  /*! Number density of electrons (in units of n_H) */
  double n_e;

  /*! Recombination rate of HII [cm^3 * s^-1] */
  double alpha_Hp_cgs;

  /*! Recombination rate of HeII [cm^3 * s^-1] */
  double alpha_Hep_cgs;

  /*! Recombination rate of HeIII [cm^3 * s^-1] */
  double alpha_Hepp_cgs;

  /*! Dielectronic recombination rate of HeII [cm^3 * s^-1] */
  double alpha_d_cgs;

  /*! Collisional ionization rate of HI [cm^3 * s^-1] */
  double gamma_eH0_cgs;

  /*! Collisional ionization rate of HeI [cm^3 * s^-1] */
  double gamma_eHe0_cgs;

  /*! Collisional ionization rate of HeII [cm^3 * s^-1] */
  double gamma_eHep_cgs;

  /*! Collisional excitation cooling coefficient of HI [erg * cm^3 * s^-1] */
  double beta_H0_cgs;

  /*! Collisional excitation cooling coefficient of HeII
   * [erg * cm^3 * s^-1] */
  double beta_Hep_cgs;

  /*! Free-free cooling coefficient [erg * cm^3 * s^-1] */
  double beta_ff_cgs;
};

/**
 * @brief Free-free (Bremsstrahlung) cooling coefficient (Cen 1992).
 *
 * This is used both to build the table and above its upper end. Note that we
 * use the same 1.43e-27 normalisation in both places, whereas Arepo switches
 * to 1.42e-27 above T_max; keeping a single value avoids a jump in the cooling
 * rate at T_max.
 *
 * @param T The temperature [K].
 * @param log10_T The log10 of the temperature (in K).
 *
 * @return The free-free cooling coefficient [erg * cm^3 * s^-1].
 */
__attribute__((always_inline)) INLINE static double
treecool_free_free_coefficient(const double T, const double log10_T) {

  return 1.43e-27 * sqrt(T) *
         (1.1 + 0.34 * exp(-(5.5 - log10_T) * (5.5 - log10_T) / 3.));
}

/**
 * @brief Construct the tables of rate coefficients.
 *
 * The tables are built on a regular grid in log10(T) running from
 * cooling->log10_T_min_cgs to cooling->log10_T_max_cgs. The fits are the ones
 * of KWH96, Table 2.
 *
 * @param cooling The #cooling_function_data to fill in.
 */
__attribute__((always_inline)) INLINE static void treecool_make_rate_table(
    struct cooling_function_data *cooling) {

  cooling->delta_log10_T =
      (cooling->log10_T_max_cgs - cooling->log10_T_min_cgs) /
      (double)(treecool_cooling_N_temperature - 1);
  cooling->inv_delta_log10_T = 1. / cooling->delta_log10_T;

  for (int i = 0; i < treecool_cooling_N_temperature; ++i) {

    const double T =
        exp10(cooling->log10_T_min_cgs + cooling->delta_log10_T * (double)i);
    const double inv_T = 1. / T;
    const double sqrt_T = sqrt(T);
    const double log10_T = log10(T);

    /* Common high-temperature suppression factor of the collisional terms */
    const double T_fact = 1. / (1. + sqrt(T * 1.e-5));

    /* Collisional excitation cooling (Cen 1992) */
    cooling->table_Beta_H0_cgs[i] =
        (118348. * inv_T < treecool_max_exponent)
            ? 7.50e-19 * exp(-118348. * inv_T) * T_fact
            : 0.;
    cooling->table_Beta_Hep_cgs[i] =
        (473638. * inv_T < treecool_max_exponent)
            ? 5.54e-17 * pow(T, -0.397) * exp(-473638. * inv_T) * T_fact
            : 0.;

    /* Free-free (Bremsstrahlung) cooling (Cen 1992) */
    cooling->table_Beta_ff_cgs[i] = treecool_free_free_coefficient(T, log10_T);

    /* Radiative recombination (Cen 1992) */
    cooling->table_Alpha_Hp_cgs[i] =
        8.40e-11 * pow(T * 1.e-3, -0.2) / (1. + pow(T * 1.e-6, 0.7)) / sqrt_T;
    cooling->table_Alpha_Hep_cgs[i] = 1.50e-10 * pow(T, -0.6353);
    cooling->table_Alpha_Hepp_cgs[i] = 4. * cooling->table_Alpha_Hp_cgs[i];

    /* Dielectronic recombination of HeII (Cen 1992) */
    cooling->table_Alpha_d_cgs[i] = (470000. * inv_T < treecool_max_exponent)
                                        ? 1.90e-3 * pow(T, -1.5) *
                                              exp(-470000. * inv_T) *
                                              (1. + 0.3 * exp(-94000. * inv_T))
                                        : 0.;

    /* Collisional ionization (Cen 1992) */
    cooling->table_Gamma_eH0_cgs[i] =
        (157809.1 * inv_T < treecool_max_exponent)
            ? 5.85e-11 * sqrt_T * exp(-157809.1 * inv_T) * T_fact
            : 0.;
    cooling->table_Gamma_eHe0_cgs[i] =
        (285335.4 * inv_T < treecool_max_exponent)
            ? 2.38e-11 * sqrt_T * exp(-285335.4 * inv_T) * T_fact
            : 0.;
    cooling->table_Gamma_eHep_cgs[i] =
        (631515.0 * inv_T < treecool_max_exponent)
            ? 5.68e-12 * sqrt_T * exp(-631515.0 * inv_T) * T_fact
            : 0.;
  }
}

/**
 * @brief Set all the rate coefficients of a gas state to zero.
 *
 * This is used outside the range covered by the tables, where the collisional
 * rates are either negligible (low temperatures) or where the fits are no
 * longer valid (high temperatures).
 *
 * @param gas (return) The #treecool_gas_state to reset.
 */
__attribute__((always_inline)) INLINE static void
treecool_zero_rate_coefficients(struct treecool_gas_state *gas) {

  gas->alpha_Hp_cgs = 0.;
  gas->alpha_Hep_cgs = 0.;
  gas->alpha_Hepp_cgs = 0.;
  gas->alpha_d_cgs = 0.;
  gas->gamma_eH0_cgs = 0.;
  gas->gamma_eHe0_cgs = 0.;
  gas->gamma_eHep_cgs = 0.;
  gas->beta_H0_cgs = 0.;
  gas->beta_Hep_cgs = 0.;
  gas->beta_ff_cgs = 0.;
}

/**
 * @brief Interpolate the tables of rate coefficients at a given temperature.
 *
 * The caller must guarantee that log10_T lies inside the range covered by the
 * tables.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param log10_T The log10 of the temperature (in K).
 * @param gas (return) The #treecool_gas_state whose coefficients we set.
 */
__attribute__((always_inline)) INLINE static void
treecool_interpolate_rate_table(const struct cooling_function_data *cooling,
                                const double log10_T,
                                struct treecool_gas_state *gas) {

  const double t =
      (log10_T - cooling->log10_T_min_cgs) * cooling->inv_delta_log10_T;

  int index = (int)t;
  index = max(index, 0);
  index = min(index, treecool_cooling_N_temperature - 2);

  const double f_hi = t - (double)index;
  const double f_lo = 1. - f_hi;

#define treecool_interpolate(table) \
  (f_lo * cooling->table[index] + f_hi * cooling->table[index + 1])

  gas->alpha_Hp_cgs = treecool_interpolate(table_Alpha_Hp_cgs);
  gas->alpha_Hep_cgs = treecool_interpolate(table_Alpha_Hep_cgs);
  gas->alpha_Hepp_cgs = treecool_interpolate(table_Alpha_Hepp_cgs);
  gas->alpha_d_cgs = treecool_interpolate(table_Alpha_d_cgs);
  gas->gamma_eH0_cgs = treecool_interpolate(table_Gamma_eH0_cgs);
  gas->gamma_eHe0_cgs = treecool_interpolate(table_Gamma_eHe0_cgs);
  gas->gamma_eHep_cgs = treecool_interpolate(table_Gamma_eHep_cgs);
  gas->beta_H0_cgs = treecool_interpolate(table_Beta_H0_cgs);
  gas->beta_Hep_cgs = treecool_interpolate(table_Beta_Hep_cgs);
  gas->beta_ff_cgs = treecool_interpolate(table_Beta_ff_cgs);

#undef treecool_interpolate
}

/**
 * @brief Compute the equilibrium abundances of the primordial species.
 *
 * This solves eq. 33-38 of KWH96 for the number densities of HI, HII, HeI,
 * HeII, HeIII and of the electrons, all expressed in units of the Hydrogen
 * number density. Since the photo-ionization terms depend on the electron
 * density, the system is solved iteratively, starting from the guess passed
 * in via the gas state.
 *
 * Below the range covered by the tables the gas is assumed to be entirely
 * neutral; above it, entirely ionized.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param log10_T The log10 of the temperature (in K).
 * @param n_H_cgs The Hydrogen number density in physical cgs units [cm^-3].
 * @param gas (return) The #treecool_gas_state to fill in. On entry, its n_e
 * field is used as the starting guess of the iteration.
 */
__attribute__((always_inline)) INLINE static void treecool_abundances(
    const struct cooling_function_data *cooling, const double log10_T,
    const double n_H_cgs, struct treecool_gas_state *gas) {

  /* Everything neutral below the table */
  if (log10_T <= cooling->log10_T_min_cgs) {

    gas->n_H0 = 1.;
    gas->n_Hp = 0.;
    gas->n_He0 = cooling->y_He;
    gas->n_Hep = 0.;
    gas->n_Hepp = 0.;
    gas->n_e = 0.;
    treecool_zero_rate_coefficients(gas);
    return;
  }

  /* Everything ionized above the table */
  if (log10_T >= cooling->log10_T_max_cgs) {

    gas->n_H0 = 0.;
    gas->n_Hp = 1.;
    gas->n_He0 = 0.;
    gas->n_Hep = 0.;
    gas->n_Hepp = cooling->y_He;
    gas->n_e = gas->n_Hp + 2. * gas->n_Hepp;
    treecool_zero_rate_coefficients(gas);
    return;
  }

  /* Collect the rate coefficients at this temperature */
  treecool_interpolate_rate_table(cooling, log10_T, gas);

  /* Make sure we start from a sensible guess */
  if (gas->n_e <= 0.) gas->n_e = 1.;

  int iter = 0;
  while (1) {

    const double n_e_old = gas->n_e;
    const double n_e_cgs = gas->n_e * n_H_cgs;

    /* Photo-ionization rates divided by the electron number density. This is
     * the combination entering eq. 33-37 of KWH96. */
    double gamma_H0_over_ne = 0.;
    double gamma_He0_over_ne = 0.;
    double gamma_Hep_over_ne = 0.;

    if (cooling->UV_background_on && n_e_cgs > treecool_min_electron_density) {

      gamma_H0_over_ne = cooling->gamma_H0_cgs / n_e_cgs;
      gamma_He0_over_ne = cooling->gamma_He0_cgs / n_e_cgs;
      gamma_Hep_over_ne = cooling->gamma_Hep_cgs / n_e_cgs;
    }

    /* Hydrogen (KWH96, eq. 33 and 34) */
    gas->n_H0 = gas->alpha_Hp_cgs /
                (gas->alpha_Hp_cgs + gas->gamma_eH0_cgs + gamma_H0_over_ne);
    gas->n_Hp = 1. - gas->n_H0;

    /* Helium (KWH96, eq. 35, 36 and 37) */
    const double ionization_He0 = gas->gamma_eHe0_cgs + gamma_He0_over_ne;

    if (ionization_He0 <= treecool_small_number) {

      /* No ionization of HeI at all */
      gas->n_He0 = cooling->y_He;
      gas->n_Hep = 0.;
      gas->n_Hepp = 0.;

    } else {

      const double recombination_Hep = gas->alpha_Hep_cgs + gas->alpha_d_cgs;
      const double ionization_Hep = gas->gamma_eHep_cgs + gamma_Hep_over_ne;

      gas->n_Hep = cooling->y_He / (1. + recombination_Hep / ionization_He0 +
                                    ionization_Hep / gas->alpha_Hepp_cgs);
      gas->n_He0 = gas->n_Hep * recombination_Hep / ionization_He0;
      gas->n_Hepp = gas->n_Hep * ionization_Hep / gas->alpha_Hepp_cgs;
    }

    /* Electrons (KWH96, eq. 38) */
    gas->n_e = gas->n_Hp + gas->n_Hep + 2. * gas->n_Hepp;

    /* Without a UV background the solution does not depend on n_e, so the
     * first pass is already the exact answer. */
    if (!cooling->UV_background_on) break;

    /* Damp the iteration to avoid oscillating around the solution */
    gas->n_e = 0.5 * (gas->n_e + n_e_old);

    if (fabs(gas->n_e - n_e_old) < treecool_ionization_tolerance) break;

    ++iter;
    if (iter >= treecool_max_iterations)
      error(
          "Ionization equilibrium failed to converge: log10(T)=%e n_H=%e "
          "n_e=%e",
          log10_T, n_H_cgs, gas->n_e);
  }
}

/**
 * @brief Compute the mean molecular weight of the gas.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param n_e The electron number density in units of n_H.
 */
__attribute__((always_inline)) INLINE static double
treecool_mean_molecular_weight(const struct cooling_function_data *cooling,
                               const double n_e) {

  return (1. + 4. * cooling->y_He) / (1. + cooling->y_He + n_e);
}

/**
 * @brief Compute the temperature of the gas from its internal energy.
 *
 * The temperature depends on the mean molecular weight, which itself depends
 * on the electron density, which in turn depends on the temperature. We
 * therefore iterate, using the damping scheme of Katz et al. (1996), until the
 * temperature and the abundances are mutually consistent.
 *
 * The relation u(T) is discontinuous at T_min: the gas is assumed to be
 * entirely neutral below it but can be highly ionized (by the UV background)
 * just above it, so mu jumps by a factor of ~2 there. For internal energies
 * falling in the gap there is no self-consistent temperature and the damped
 * iteration oscillates for ever around T_min. If the iteration has not
 * converged after a fixed number of steps we therefore fall back to a
 * bisection of g(T) = T - T(mu(T)), which is bracketed by the temperatures
 * obtained with the fully ionized and with the neutral mean molecular
 * weights. Where a solution exists, g is monotonic and the bisection finds the
 * same root as the fixed-point scheme; in the gap it returns T_min.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param u_cgs The internal energy per unit mass in physical cgs units
 * [erg * g^-1].
 * @param n_H_cgs The Hydrogen number density in physical cgs units [cm^-3].
 * @param gas (return) The #treecool_gas_state to fill in. On entry, its n_e
 * field is used as the starting guess of the iteration.
 *
 * @return The temperature of the gas [K].
 */
__attribute__((always_inline)) INLINE static double treecool_temperature_from_u(
    const struct cooling_function_data *cooling, const double u_cgs,
    const double n_H_cgs, struct treecool_gas_state *gas) {

  /* (gamma - 1) * u * m_H / k_B, i.e. the temperature divided by mu */
  const double T_over_mu = hydro_gamma_minus_one * u_cgs *
                           cooling->proton_mass_cgs / cooling->boltzmann_k_cgs;

  double mu = treecool_mean_molecular_weight(cooling, gas->n_e);
  double T = T_over_mu * mu;

  double T_old, damping = 0.;
  int iter = 0;

  do {

    const double n_e_old = gas->n_e;

    treecool_abundances(cooling, log10(T), n_H_cgs, gas);

    T_old = T;

    mu = treecool_mean_molecular_weight(cooling, gas->n_e);
    const double T_new = T_over_mu * mu;

    /* Estimate how strongly the temperature reacts to a change in n_e and damp
     * the update accordingly (Katz et al. 1996). */
    damping =
        max(damping, T_new / (1. + cooling->y_He + gas->n_e) *
                         fabs((gas->n_e - n_e_old) / (T_new - T_old + 1.)));

    T = T_old + (T_new - T_old) / (1. + damping);

    ++iter;

  } while (fabs(T - T_old) > treecool_temperature_tolerance * T &&
           iter < treecool_temperature_max_fixed_point_iterations);

  if (iter < treecool_temperature_max_fixed_point_iterations) return T;

  /* The fixed-point iteration did not converge: bisect instead. The
   * temperature is bracketed by the fully ionized (smallest mu) and the
   * neutral (largest mu) solutions. */
  double T_lo = T_over_mu * treecool_mean_molecular_weight(
                                cooling, /*n_e=*/1. + 2. * cooling->y_He);
  double T_hi = T_over_mu * treecool_mean_molecular_weight(cooling, /*n_e=*/0.);

  iter = 0;
  do {

    T = 0.5 * (T_lo + T_hi);

    treecool_abundances(cooling, log10(T), n_H_cgs, gas);
    const double T_new =
        T_over_mu * treecool_mean_molecular_weight(cooling, gas->n_e);

    /* g(T) = T - T_new is increasing in T: a negative value means that the
     * root lies above the current guess */
    if (T_new > T)
      T_lo = T;
    else
      T_hi = T;

    ++iter;

  } while (T_hi - T_lo > treecool_temperature_tolerance * T_lo &&
           iter < treecool_temperature_max_bisection_iterations);

  if (iter >= treecool_temperature_max_bisection_iterations)
    error("Temperature failed to converge: u=%e n_H=%e T_lo=%e T_hi=%e", u_cgs,
          n_H_cgs, T_lo, T_hi);

  /* Leave the abundances consistent with the temperature we return */
  T = 0.5 * (T_lo + T_hi);
  treecool_abundances(cooling, log10(T), n_H_cgs, gas);

  return T;
}

/**
 * @brief Compute the net cooling rate of the gas at a given temperature.
 *
 * This is the quantity (Heating - Cooling) / n_H^2 of KWH96, Table 1, in
 * physical cgs units [erg * cm^3 * s^-1]. A negative value means that the gas
 * is cooling.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param log10_T The log10 of the temperature (in K).
 * @param n_H_cgs The Hydrogen number density in physical cgs units [cm^-3].
 * @param gas (return) The #treecool_gas_state to fill in. On entry, its n_e
 * field is used as the starting guess of the iteration.
 */
__attribute__((always_inline)) INLINE static double treecool_cooling_rate(
    const struct cooling_function_data *cooling, double log10_T,
    const double n_H_cgs, struct treecool_gas_state *gas) {

  /* Never evaluate the rates below the table: the gas would be entirely
   * neutral and the cooling rate exactly zero. Instead, as in KWH96's
   * implementation, evaluate them in the middle of the first bin. */
  if (log10_T <= cooling->log10_T_min_cgs)
    log10_T = cooling->log10_T_min_cgs + 0.5 * cooling->delta_log10_T;

  const double T = exp10(log10_T);

  double Lambda_Compton = 0.;
  double Lambda, Heat;

  if (log10_T < cooling->log10_T_max_cgs) {

    treecool_abundances(cooling, log10_T, n_H_cgs, gas);

    /* Collisional excitation */
    const double Lambda_exc_H0 = gas->beta_H0_cgs * gas->n_e * gas->n_H0;
    const double Lambda_exc_Hep = gas->beta_Hep_cgs * gas->n_e * gas->n_Hep;

    /* Collisional ionization. Each ionization removes the ionization
     * potential of the species: 13.6 eV (HI), 24.6 eV (HeI) and 54.4 eV
     * (HeII), i.e. 2.18e-11, 3.94e-11 and 8.72e-11 erg. */
    const double Lambda_ion_H0 =
        2.18e-11 * gas->gamma_eH0_cgs * gas->n_e * gas->n_H0;
    const double Lambda_ion_He0 =
        3.94e-11 * gas->gamma_eHe0_cgs * gas->n_e * gas->n_He0;
    const double Lambda_ion_Hep =
        8.72e-11 * gas->gamma_eHep_cgs * gas->n_e * gas->n_Hep;

    /* Recombination. Each radiative recombination removes ~0.75 k_B T
     * (1.036e-16 erg/K * T). The dielectronic term uses the ratio of the
     * cooling and rate fits of Cen (1992), 1.24e-13 / 1.90e-3 = 6.526e-11 erg,
     * i.e. ~40.7 eV per recombination. */
    const double Lambda_rec_Hp =
        1.036e-16 * T * gas->n_e * gas->alpha_Hp_cgs * gas->n_Hp;
    const double Lambda_rec_Hep =
        1.036e-16 * T * gas->n_e * gas->alpha_Hep_cgs * gas->n_Hep;
    const double Lambda_rec_Hepp =
        1.036e-16 * T * gas->n_e * gas->alpha_Hepp_cgs * gas->n_Hepp;
    const double Lambda_rec_Hep_d =
        6.526e-11 * gas->alpha_d_cgs * gas->n_e * gas->n_Hep;

    /* Free-free (Bremsstrahlung) */
    const double Lambda_ff = gas->beta_ff_cgs * gas->n_e *
                             (gas->n_Hp + gas->n_Hep + 4. * gas->n_Hepp);

    Lambda = Lambda_exc_H0 + Lambda_exc_Hep + Lambda_ion_H0 + Lambda_ion_He0 +
             Lambda_ion_Hep + Lambda_rec_Hp + Lambda_rec_Hep + Lambda_rec_Hepp +
             Lambda_rec_Hep_d + Lambda_ff;

    /* Photo-heating by the UV background */
    if (cooling->UV_background_on) {

      Heat = (gas->n_H0 * cooling->epsilon_H0_cgs +
              gas->n_He0 * cooling->epsilon_He0_cgs +
              gas->n_Hep * cooling->epsilon_Hep_cgs) /
             n_H_cgs;
    } else {

      Heat = 0.;
    }

  } else {

    /* Above the table the gas is fully ionized. Only free-free emission and
     * the inverse Compton cooling survive and there is no heating. */

    gas->n_H0 = 0.;
    gas->n_Hp = 1.;
    gas->n_He0 = 0.;
    gas->n_Hep = 0.;
    gas->n_Hepp = cooling->y_He;
    gas->n_e = gas->n_Hp + 2. * gas->n_Hepp;

    treecool_zero_rate_coefficients(gas);

    gas->beta_ff_cgs = treecool_free_free_coefficient(T, log10_T);

    Lambda = gas->beta_ff_cgs * gas->n_e * (gas->n_Hp + 4. * gas->n_Hepp);
    Heat = 0.;
  }

  /* Inverse Compton cooling off the CMB. The coefficient is
   * 4 sigma_T a_rad k_B T_CMB,0^4 / (m_e c) [erg * s^-1 * K^-1] evaluated
   * with T_CMB,0 = 2.73 K, as in Arepo; the (1 + z)^4 factor then
   * scales the CMB energy density to the current redshift. */
  if (cooling->with_Compton_cooling) {

    Lambda_Compton = 5.65e-36 * gas->n_e * (T - cooling->T_CMB_cgs) *
                     cooling->one_plus_z_to_the_4 / n_H_cgs;
  }

  return Heat - Lambda - Lambda_Compton;
}

/**
 * @brief Compute the net cooling rate of the gas from its internal energy.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param u_cgs The internal energy per unit mass in physical cgs units
 * [erg * g^-1].
 * @param n_H_cgs The Hydrogen number density in physical cgs units [cm^-3].
 * @param gas (return) The #treecool_gas_state to fill in. On entry, its n_e
 * field is used as the starting guess of the iteration.
 *
 * @return (Heating - Cooling) / n_H^2 in physical cgs units
 * [erg * cm^3 * s^-1].
 */
__attribute__((always_inline)) INLINE static double
treecool_cooling_rate_from_u(const struct cooling_function_data *cooling,
                             const double u_cgs, const double n_H_cgs,
                             struct treecool_gas_state *gas) {

  const double T = treecool_temperature_from_u(cooling, u_cgs, n_H_cgs, gas);

  return treecool_cooling_rate(cooling, log10(T), n_H_cgs, gas);
}

#endif /* SWIFT_COOLING_RATES_TREECOOL_H */
