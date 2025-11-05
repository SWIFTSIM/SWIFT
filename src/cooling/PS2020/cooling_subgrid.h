/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PS2020_COOLING_SUBGRID_H
#define SWIFT_PS2020_COOLING_SUBGRID_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "cooling.h"
#include "equation_of_state.h"

/**
 * @brief Computes the gas pressure at thermal equilibrium (cooling = heating)
 * for gas with at a given redshift, abundance ratio, metallicity, and density
 *
 * @return The logarithm base 10 of the thermal pressure at thermal equilibrium
 * in cgs
 *
 * @param cooling The properties of the cooling scheme.
 * @param abundance_ratio element abundance ratio of gas particle
 * @param log_u_cgs The logarithm base 10 of the internal energy in cgs
 * @param log10nH_local The logarithm base 10 of the hydrogen particle density
 * in cgs
 * @param rho_cgs Physical density of the gas in cgs
 * @param redshift The redshift
 * @param ired Index of redshift within the redshift dimension of the cooling
 * tables
 * @param imet Index of metallicity within the metallicity dimension of the
 * cooling tables
 * @param dred delta redshift between redshift bins, for interpolation
 * @param dmet delta metallicity between metallicity, for interpolation
 */

double get_thermal_equilibrium_pressure(
    const struct cooling_function_data *cooling,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const double log_u_cgs, const float log10nH_local, const double rho_cgs,
    const double redshift, const int ired, const int imet, const float dred,
    const float dmet) {

  /* The thermal equilibrium might have multiple solutions for a given density
   * (WNM, CNM), so it is dangerous to do a bisection; instead the internal
   * energy is decreased in steps of dlogu until heating > cooling (lambda > 0),
   * then interpolate between u_previous and u_current*/
  const double du = exp10(0.3f);

  /* Get H number density */
  double nH_cgs = exp10(log10nH_local);

  /* Get index along the density index of the table */
  int iden;
  float dden;
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10nH_local, &iden,
               &dden);

  /* Prepare for iterations */
  double u_current = exp10(log_u_cgs);
  double u_prev = u_current;

  double Lambda_current = colibre_cooling_rate(
      log10f(u_current), redshift, nH_cgs, abundance_ratio, iden, dden, imet,
      dmet, ired, dred, cooling, /*onlyicool=*/0,
      /*onlyiheat=*/0, /*icool=*/0, /*iheat=*/0);

  /* Are we in the case where we actually heat?
   * This can happen but should not */
  if (Lambda_current > 0.f) {

    /* Increase u_current by a factor of 10 */
    u_current *= 10.f;

    Lambda_current = colibre_cooling_rate(
        log10f(u_current), redshift, nH_cgs, abundance_ratio, iden, dden, imet,
        dmet, ired, dred, cooling, /*onlyicool=*/0,
        /*onlyiheat=*/0, /*icool=*/0, /*iheat=*/0);

    /* If we are still heating, print a message */
    if (Lambda_current > 0.f) message("Should not be here!");
  }

  /* The equilibirum energy we are looking for */
  double u_eq_cgs = -1.;

  /* Main iteration */
  do {

    Lambda_current = colibre_cooling_rate(
        log10f(u_current), redshift, nH_cgs, abundance_ratio, iden, dden, imet,
        dmet, ired, dred, cooling, 0, 0, 0, 0);

    /* still cooling, need to decrease u */
    if (Lambda_current < 0) {

      u_current = u_current / du;

      /* heating: interpolate to find eq */
    } else {

      u_prev = u_current * du;
      double Lambda_prev = colibre_cooling_rate(
          log10f(u_prev), redshift, nH_cgs, abundance_ratio, iden, dden, imet,
          dmet, ired, dred, cooling, /*onlyicool=*/0,
          /*onlyiheat=*/0, /*icool=*/0, /*iheat=*/0);

      if (Lambda_prev * Lambda_current > 0.) {
        message("Something wrong! %.4e, %.4e", Lambda_prev, Lambda_current);
      }

      /* Interpolate to get the equilibrium internal energy */
      u_eq_cgs = u_prev - Lambda_prev / (Lambda_current - Lambda_prev) *
                              (u_current - u_prev);

      /* And we are done here */
      break;
    }

  } while (u_current >= cooling->umin_cgs);

  /* Equilibrium internal energy is below the minimum internal energy
   * set u_eq_cgs to minimum internal energy */
  if (u_eq_cgs < 0.f) {
    u_eq_cgs = cooling->umin_cgs;
  }

  /* Finish by computing the pressure from the density
   * and equilibirum internal energy */
  const float P_eq_cgs = gas_pressure_from_internal_energy(rho_cgs, u_eq_cgs);
  return (double)log10f(P_eq_cgs);
}

/**
 * @brief Compute the subgrid properties of the gas.
 * For particles on the entropy floor, we use pressure equilibrium to
 * infer the properties of the particle.
 * Note that we return the density in physical coordinates.
 *
 * @return depends on the input for isub:
 *   isub = colibre_compute_subgrid_density: the physical subgrid density in
 * internal units isub = colibre_compute_subgrid_temperature: the subgrid
 * temperature temperature isub = colibre_compute_subgrid_HI_fraction: the
 * subgrid HI fraction (n_HI / n_H) isub =
 * colibre_compute_subgrid_HII_fraction: the subgrid HII fraction (n_HII /
 * n_H) isub = colibre_compute_subgrid_H2_fraction: the subgrid HII fraction
 * (n_H2 / n_H)
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param rho_phys Density of the gas in internal physical units.
 * @param logZZsol Logarithm base 10 of the gas' metallicity in units of solar
 * metallicity.
 * @param XH The Hydrogen abundance of the gas.
 * @param P_phys Pressure of the gas in internal physical units.
 * @param log10_T The logarithm base 10 of the temperature of the gas.
 * @param log10_u_EOS_max_cgs The logarithm base 10 of the maximal energy in cgs
 * to be considered on the EOS at the density of the gas.
 * @param HII_region Is this patch of gas in an HII region?
 * @param abundance_ratio element abundance ratio of gas particle
 * @param log_u_cgs The logarithm base 10 of the internal energy in cgs
 * @param isub which subgrid property to calculate
 */
double compute_subgrid_property(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const float rho_phys, const float logZZsol,
    const float XH, const float P_phys, const float log10_T,
    const float log10_u_EOS_max_cgs, const int HII_region,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const double log_u_cgs, const enum cooling_subgrid_properties isub) {

  if (HII_region)
    error("HII regions are not implemented in the EAGLE-XL flavour");

  const float weights[3] = {1.0, 1.0, 1.0};

  /* Set the value for returning the values based on the hydro properties,
   * e.g. for particles above the EOS, or if no EOS is present */

  double standard_return;

  /* Start by computing the subgrid property as if it were not subgrid */

  if (isub == cooling_compute_subgrid_density) {
    standard_return = rho_phys;
  } else if (isub == cooling_compute_subgrid_temperature) {
    standard_return = exp10(log10_T);
  } else {

    /* Get the Hydrogen number density */
    const double n_H = XH * rho_phys / phys_const->const_proton_mass;
    const double n_H_cgs = n_H * cooling->number_density_to_cgs;

    /* Get index along the different table axis for this particle */
    int ired, imet, iden, item;
    float dred, dmet, dden, dtem;
    get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
                 &ired, &dred);
    get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
                 &imet, &dmet);
    get_index_1d(cooling->nH, colibre_cooling_N_density, log10f(n_H_cgs), &iden,
                 &dden);
    get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
                 &dtem);

    if (isub == cooling_compute_subgrid_HI_fraction) {

      const float nHI_over_nH =
          interpolation4d_plus_summation(cooling->table.logHfracs_all,  /* */
                                         weights, neutral, neutral,     /* */
                                         ired, item, imet, iden,        /* */
                                         dred, dtem, dmet, dden,        /* */
                                         colibre_cooling_N_redshifts,   /* */
                                         colibre_cooling_N_temperature, /* */
                                         colibre_cooling_N_metallicity, /* */
                                         colibre_cooling_N_density,     /* */
                                         3);

      standard_return = nHI_over_nH;

    } else if (isub == cooling_compute_subgrid_HII_fraction) {

      const float nHII_over_nH =
          interpolation4d_plus_summation(cooling->table.logHfracs_all,  /* */
                                         weights, ionized, ionized,     /* */
                                         ired, item, imet, iden,        /* */
                                         dred, dtem, dmet, dden,        /* */
                                         colibre_cooling_N_redshifts,   /* */
                                         colibre_cooling_N_temperature, /* */
                                         colibre_cooling_N_metallicity, /* */
                                         colibre_cooling_N_density,     /* */
                                         3);

      standard_return = nHII_over_nH;

    } else if (isub == cooling_compute_subgrid_H2_fraction) {

      const float mH2_over_mH =
          interpolation4d_plus_summation(cooling->table.logHfracs_all,  /* */
                                         weights, molecular, molecular, /* */
                                         ired, item, imet, iden,        /* */
                                         dred, dtem, dmet, dden,        /* */
                                         colibre_cooling_N_redshifts,   /* */
                                         colibre_cooling_N_temperature, /* */
                                         colibre_cooling_N_metallicity, /* */
                                         colibre_cooling_N_density,     /* */
                                         3);

      /* molecular fraction are stored as mass fractions in the table (2nH2/nH)
       * but we want particle fractions (nH2/nH) as output */
      standard_return = 0.5 * mH2_over_mH;

    } else {

      /* This should never happen! */
      standard_return = -1.;
      error("Unknown subgrid property to calculate");
    }
  }

  /* We can now look at whether we need to correct that first guess */

  /* Return case 1:
   * Particle is above the equation of state, nothing to be done */
  if (log_u_cgs >= log10_u_EOS_max_cgs) {
    return standard_return;
  }

  /* Return case 2:
   * Particle is in an HII region
   * get density by mapping particle onto T_HII and mu_HII at equal pressure */

  /* NOTE: HII regions DO NOT exit in the EAGLE-XL flavour
   * --> Nothing to do here, move straight to case 3 */

  double Lambda;
  /* Get the Hydrogen number density */
  const double n_H = XH * rho_phys / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis for this particle */
  int ired, imet, iden, item;
  float dred, dmet, dden, dtem;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z, &ired,
               &dred);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10f(n_H_cgs), &iden,
               &dden);
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
               &dtem);

  /* Get the cooling rate for the particle properties */
  Lambda =
      colibre_cooling_rate(log_u_cgs, cosmo->z, n_H_cgs, abundance_ratio, iden,
                           dden, imet, dmet, ired, dred, cooling, 0, 0, 0, 0);

  /* Return case 3:
   * Particle is on the EOS, but heating dominates cooling, so the
   * particle is below or at its thermal equilibrium temperature;
   * This might happen for a high EOS and low densities, where
   * radiative cooling does not dominate the total (net) cooling rate;
   * In this case we return the hydro properties */
  if (Lambda >= 0.) {
    return standard_return;
  }

  /* Particle is on the EOS and cooling dominates heating, so the
   * subgrid property is calculated by projecting the hydro properties
   * on the thermal equilibrium temperature where the pressure equals
   * the hydro pressure */

  if (Lambda < 0) {
    float log10_Teq = interpolation_3d(cooling->table.logTeq,         /* */
                                       ired, imet, iden,              /* */
                                       dred, dmet, dden,              /* */
                                       colibre_cooling_N_redshifts,   /* */
                                       colibre_cooling_N_metallicity, /* */
                                       colibre_cooling_N_density);

    const float P_cgs = P_phys * cooling->pressure_to_cgs;
    const float log10_P_cgs = log10f(P_cgs);

    /* Maximal equilibrium pressure (n = n_max, and Teq = Teq(n_max)) from the
     * table at the current redshift and metallicity */
    const float log10_Peq_max_cgs =
        interpolation_3d_no_z(cooling->table.logPeq,                     /* */
                              ired, imet, colibre_cooling_N_density - 1, /* */
                              dred, dmet, 0.,                            /* */
                              colibre_cooling_N_redshifts,               /* */
                              colibre_cooling_N_metallicity,             /* */
                              colibre_cooling_N_density);

    /* Return case 4:
     * EOS pressure (logP) is larger than maximum Peq (possible for very
     * steep EOS); use the equilibrium temperature from the highest density
     * bin */
    if (log10_P_cgs > log10_Peq_max_cgs) {

      if (isub == cooling_compute_subgrid_density) {

        const double rho_cgs =
            exp10(cooling->nH[colibre_cooling_N_density - 1]) *
            cooling->proton_mass_cgs / XH;

        return rho_cgs * cooling->density_from_cgs;

      } else if (isub == cooling_compute_subgrid_temperature) {

        const float log10_T_at_Peq = interpolation_3d_no_z(
            cooling->table.logTeq, ired, imet, colibre_cooling_N_density - 1,
            dred, dmet, 0., colibre_cooling_N_redshifts,
            colibre_cooling_N_metallicity, colibre_cooling_N_density);

        return exp10(log10_T_at_Peq);

      } else if (isub == cooling_compute_subgrid_HI_fraction) {

        const float log10_HI_over_nH = interpolation_4d_no_z_no_w(
            cooling->table.logHfracs_Teq,                       /* */
            ired, imet, colibre_cooling_N_density - 1, neutral, /* */
            dred, dmet, 0., 0.,                                 /* */
            colibre_cooling_N_redshifts,                        /* */
            colibre_cooling_N_metallicity,                      /* */
            colibre_cooling_N_density,                          /* */
            3);

        return exp10(log10_HI_over_nH);

      } else if (isub == cooling_compute_subgrid_HII_fraction) {

        const float log10_HII_over_nH = interpolation_4d_no_z_no_w(
            cooling->table.logHfracs_Teq,                       /* */
            ired, imet, colibre_cooling_N_density - 1, ionized, /* */
            dred, dmet, 0., 0.,                                 /* */
            colibre_cooling_N_redshifts,                        /* */
            colibre_cooling_N_metallicity,                      /* */
            colibre_cooling_N_density,                          /* */
            3);

        return exp10(log10_HII_over_nH);

      } else if (isub == cooling_compute_subgrid_H2_fraction) {

        const float log10_H2_over_nH = interpolation_4d_no_z_no_w(
            cooling->table.logHfracs_Teq,                         /* */
            ired, imet, colibre_cooling_N_density - 1, molecular, /* */
            dred, dmet, 0., 0.,                                   /* */
            colibre_cooling_N_redshifts,                          /* */
            colibre_cooling_N_metallicity,                        /* */
            colibre_cooling_N_density,                            /* */
            3);

        return 0.5 * exp10(log10_H2_over_nH);

      } else {
        error("Unknown subgrid property to calculate");
        return 0.;
      }
    }

    /* Find equilibrium property by increasing the density
     * Note that the logPeq table is neither equally spaced nor
     * necessarilly monotically increasing.
     * We hence loop over densities and pick the first one where
     * log10_P < log10_Peq. We start with the next resolved density (index
     * iden+1), as we require the subgrid density to be larger than the
     * physical density */

    float log10_Peq_prev = log10_P_cgs;

    for (int i = iden + 1; i < colibre_cooling_N_density; i++) {

      float log10_Peq_current, log10_n_at_Peq;

      if ((log10_Teq <= log10_T) &&
          (logZZsol <=
           cooling->Metallicity[colibre_cooling_N_metallicity - 1])) {

        /* Standard case: the thermal equilibrium temperature is below the
         * EOS temperature and the metallicity is within the table range */
        log10_Peq_current = interpolation_3d_no_z(
            cooling->table.logPeq, ired, imet, i, dred, dmet, 0.,
            colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
            colibre_cooling_N_density);

      } else {
        /* The thermal equilibrium tables cannot be used, need to calculate
         * the thermal equilibrium temperatures; this only occurs when
         * the metallicity is higher than the maximum metallicity in the
         * tables */

        const double rho_cgs =
            exp10(cooling->nH[i]) * cooling->proton_mass_cgs / XH;
        log10_Peq_current = get_thermal_equilibrium_pressure(
            cooling, abundance_ratio, log_u_cgs, cooling->nH[i], rho_cgs,
            cosmo->z, ired, imet, dred, dmet);
      }

      if (log10_Peq_current >= log10_P_cgs) {

        if (log10_Peq_current == log10_Peq_prev) {

          /* the interpolated pressure exactly equals the input pressure
           * on the first density iteration (where log10_Peq_prev =
           * log10_P_cgs), no interpolation necessary */

          log10_n_at_Peq = cooling->nH[i];

        } else {

          /* first density bin where log P_eq >= log P_SPH
           * the solution is therefore between log10_Peq_prev and
           * log10_Peq_current */

          const float delta_P_eq = (log10_P_cgs - log10_Peq_prev) /
                                   (log10_Peq_current - log10_Peq_prev);

          /* Interpolate to get the density at equilibrium */
          log10_n_at_Peq = cooling->nH[i - 1] +
                           delta_P_eq * (cooling->nH[i] - cooling->nH[i - 1]);
        }

        /* Return case 5:
         * SPH density and temperatures are projected onto the thermal
         * equilibrium temperature function Teq(density,
         * metallicity/abundances, redshift) for equal pressure */

        if (isub == cooling_compute_subgrid_density) {
          const double rho_cgs =
              exp10(log10_n_at_Peq) * cooling->proton_mass_cgs / XH;
          return rho_cgs * cooling->density_from_cgs;

        } else {

          /* We found a valid density: Get the index along the density
           * axis of the tables. */
          int iden_eq;
          float dden_eq;

          get_index_1d(cooling->nH, colibre_cooling_N_density, log10_n_at_Peq,
                       &iden_eq, &dden_eq);

          /* Finish by interpolating the tables to get the temperature */
          const float log10_T_at_Peq =
              interpolation_3d(cooling->table.logTeq,         /* */
                               ired, imet, iden_eq,           /* */
                               dred, dmet, dden_eq,           /* */
                               colibre_cooling_N_redshifts,   /* */
                               colibre_cooling_N_metallicity, /* */
                               colibre_cooling_N_density);

          if (isub == cooling_compute_subgrid_temperature) {
            return exp10(log10_T_at_Peq);

          } else {

            int item_eq;
            float dtem_eq;

            get_index_1d(cooling->Temp, colibre_cooling_N_temperature,
                         log10_T_at_Peq, &item_eq, &dtem_eq);

            if (isub == cooling_compute_subgrid_HI_fraction) {

              const float nHI_over_nH_eq = interpolation4d_plus_summation(
                  cooling->table.logHfracs_all,  /* */
                  weights, neutral, neutral,     /* */
                  ired, item_eq, imet, iden_eq,  /* */
                  dred, dtem_eq, dmet, dden_eq,  /* */
                  colibre_cooling_N_redshifts,   /* */
                  colibre_cooling_N_temperature, /* */
                  colibre_cooling_N_metallicity, /* */
                  colibre_cooling_N_density,     /* */
                  3);

              return nHI_over_nH_eq;

            } else if (isub == cooling_compute_subgrid_HII_fraction) {

              const float nHII_over_nH_eq = interpolation4d_plus_summation(
                  cooling->table.logHfracs_all,  /* */
                  weights, ionized, ionized,     /* */
                  ired, item_eq, imet, iden_eq,  /* */
                  dred, dtem_eq, dmet, dden_eq,  /* */
                  colibre_cooling_N_redshifts,   /* */
                  colibre_cooling_N_temperature, /* */
                  colibre_cooling_N_metallicity, /* */
                  colibre_cooling_N_density,     /* */
                  3);

              return nHII_over_nH_eq;

            } else if (isub == cooling_compute_subgrid_H2_fraction) {

              const float mH2_over_mH_eq = interpolation4d_plus_summation(
                  cooling->table.logHfracs_all,  /* */
                  weights, molecular, molecular, /* */
                  ired, item_eq, imet, iden_eq,  /* */
                  dred, dtem_eq, dmet, dden_eq,  /* */
                  colibre_cooling_N_redshifts,   /* */
                  colibre_cooling_N_temperature, /* */
                  colibre_cooling_N_metallicity, /* */
                  colibre_cooling_N_density,     /* */
                  3);

              return 0.5 * mH2_over_mH_eq;

            } else {

              error("Unknown subgrid property to calculate");
              return 0.;
            }
          }
        }
      }

      /* Move to the next iteration */
      log10_Peq_prev = log10_Peq_current;
    }
  }

  /* Return case 6:
   * Nothing worked, return the property at the SPH density and temperature
   * The code should not reach this point during the simulation, but still
   * does before and at the first timestep */
  return standard_return;
}

#endif /* SWIFT_PS2020_COOLING_SUBGRID_H */
