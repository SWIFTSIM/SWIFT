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
#ifndef SWIFT_COLIBRE_COOLING_SUBGRID_H
#define SWIFT_COLIBRE_COOLING_SUBGRID_H

/* Config parameters. */
#include "../config.h"

/**
 * @brief Compute the fraction of Hydrogen that is in HI based
 * on the pressure of the gas.
 *
 * For particles on the entropy floor, we use pressure equilibrium to
 * infer the properties of the particle.
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param p The #part.
 * @param xp The #xpart.
 * @param log10_T The logarithm base 10 of the temperature of the particle.
 * @param log10_T_EOS_max The logarithm base 10 of the maximal temperature
 * to be considered on the EOS at the density of this particle.
 */
double compute_subgrid_HI_fraction(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const struct part *p, const struct xpart *xp,
    const float log10_T, const float log10_T_EOS_max) {

  /* Physical density of this particle */
  const double rho_phys = hydro_get_physical_density(p, cosmo);

  /* Get the total metallicity */
  float dummy[colibre_cooling_N_elementtypes];
  const float logZZsol = abundance_ratio_to_solar(p, cooling, dummy);

  /* Get the Hydrogen number density */
  const float *const metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];
  const double n_H = XH * rho_phys / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis for this particle */
  int ired, imet, iden, item;
  float dred, dmet, dden, dtem;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z, &ired,
               &dred);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &iden,
               &dden);
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
               &dtem);

  /* Are we above or below the EoS? */
  if (log10_T < log10_T_EOS_max) {

    /* We are below the EoS. Use subgrid properties assuming P equilibrium */

    const float P = hydro_get_physical_pressure(p, cosmo);
    const float P_cgs = P * cooling->pressure_to_cgs;
    const float log10_P_cgs = log10(P_cgs);

    /* Recover the maximal equilibrium pressure from the table at the current
     * redshift and metallicity */
    const float log10_Peq_max_cgs =
        interpolation_3d_no_z(cooling->table.logPeq,                     /* */
                              ired, imet, colibre_cooling_N_density - 1, /* */
                              dred, dmet, 0.,                            /* */
                              colibre_cooling_N_redshifts,               /* */
                              colibre_cooling_N_metallicity,             /* */
                              colibre_cooling_N_density);

    /* Are we above or below the max pressure equilibrium? */
    if (log10_P_cgs > log10_Peq_max_cgs) {

      /* EOS pressure (logP) is larger than maximum Peq (can happen for very
       * steep EOS) HI fraction from the highest density bin */

      const float log10_HI_over_nH = interpolation_4d_no_z_no_w(
          cooling->table.logHfracs_Teq,                       /* */
          ired, imet, colibre_cooling_N_density - 1, neutral, /* */
          dred, dmet, 0., 0.,                                 /* */
          colibre_cooling_N_redshifts,                        /* */
          colibre_cooling_N_metallicity,                      /* */
          colibre_cooling_N_density,                          /* */
          3);

      return exp10f(log10_HI_over_nH);

    } else {

      /* Normal case: We are not beyond the table range */

      /* Need to find thermal equilibrium state with the same pressure
       * by interpolating the equilibrium table
       *
       * Note that the logPeq table is neither equally spaced nor
       * necessarilly monotically increasing.
       * We hance loop over densities and pick the first one where
       * log10_P < log10_Peq. We start with the resolved density (index iden),
       * as we expect the subgrid density to be larger */

      /* If the solution can't be found, we revert to the non-subgrid
       * density */
      int iden_eq = iden;
      float dden_eq = dden;

      /* Loop over the density bins */
      for (int i = iden; i < colibre_cooling_N_density; i++) {

        /* Equilibirum pressure at this density */
        const float log10_Peq_interp =
            interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                  ired, imet, i,                 /* */
                                  dred, dmet, 0.,                /* */
                                  colibre_cooling_N_redshifts,   /* */
                                  colibre_cooling_N_metallicity, /* */
                                  colibre_cooling_N_density);

        /* Did we find a solution? */
        if (log10_Peq_interp > log10_P_cgs) {

          /* Equilibrium pressure at the previous density point */
          const float log10_Peq_prev =
              interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                    ired, imet, i - 1,             /* */
                                    dred, dmet, 0.,                /* */
                                    colibre_cooling_N_redshifts,   /* */
                                    colibre_cooling_N_metallicity, /* */
                                    colibre_cooling_N_density);

          /* How far from the equilibrium point are we? */
          const float delta_P_eq = (log10_P_cgs - log10_Peq_prev) /
                                   (log10_Peq_interp - log10_Peq_prev);

          /* Interpolate to get the density at equilibrium */
          const float log10_n_at_Peq =
              cooling->nH[i - 1] +
              delta_P_eq * (cooling->nH[i] - cooling->nH[i - 1]);

          /* We found a valid density: Get the index along the density
           * axis of the tables. */
          get_index_1d(cooling->nH, colibre_cooling_N_density, log10_n_at_Peq,
                       &iden_eq, &dden_eq);

          break;
        }
      }

      /* Finish by interpolating the tables to get the HI density */
      const float log10_HI_over_nH =
          interpolation_4d_no_w(cooling->table.logHfracs_Teq,  /* */
                                ired, imet, iden_eq, neutral,  /* */
                                dred, dmet, dden_eq, 0.,       /* */
                                colibre_cooling_N_redshifts,   /* */
                                colibre_cooling_N_metallicity, /* */
                                colibre_cooling_N_density,     /* */
                                3);

      return exp10f(log10_HI_over_nH);
    }

  } else {

    /* We are above the EoS. */
    const float weights[3] = {1.0, 1.0, 1.0};
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

    return nHI_over_nH;
  }
}

/**
 * @brief Compute the fraction of Hydrogen that is in HII based
 * on the pressure of the gas.
 *
 * For particles on the entropy floor, we use pressure equilibrium to
 * infer the properties of the particle.
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param p The #part.
 * @param xp The #xpart.
 * @param log10_T The logarithm base 10 of the temperature of the particle.
 * @param log10_T_EOS_max The logarithm base 10 of the maximal temperature
 * to be considered on the EOS at the density of this particle.
 */
double compute_subgrid_HII_fraction(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const struct part *p, const struct xpart *xp,
    const float log10_T, const float log10_T_EOS_max) {

  /* Physical density of this particle */
  const double rho_phys = hydro_get_physical_density(p, cosmo);

  /* Get the total metallicity */
  float dummy[colibre_cooling_N_elementtypes];
  const float logZZsol = abundance_ratio_to_solar(p, cooling, dummy);

  /* Get the Hydrogen number density */
  const float *const metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];
  const double n_H = XH * rho_phys / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis for this particle */
  int ired, imet, iden, item;
  float dred, dmet, dden, dtem;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z, &ired,
               &dred);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &iden,
               &dden);
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
               &dtem);

  /* Are we above or below the EoS? */
  if (log10_T < log10_T_EOS_max) {

    /* We are below the EoS. Use subgrid properties assuming P equilibrium */

    const float P = hydro_get_physical_pressure(p, cosmo);
    const float P_cgs = P * cooling->pressure_to_cgs;
    const float log10_P_cgs = log10(P_cgs);

    /* Recover the maximal equilibrium pressure from the table at the current
     * redshift and metallicity */
    const float log10_Peq_max_cgs =
        interpolation_3d_no_z(cooling->table.logPeq,                     /* */
                              ired, imet, colibre_cooling_N_density - 1, /* */
                              dred, dmet, 0.,                            /* */
                              colibre_cooling_N_redshifts,               /* */
                              colibre_cooling_N_metallicity,             /* */
                              colibre_cooling_N_density);

    /* Are we above or below the max pressure equilibrium? */
    if (log10_P_cgs > log10_Peq_max_cgs) {

      /* EOS pressure (logP) is larger than maximum Peq (can happen for very
       * steep EOS) HII fraction from the highest density bin */

      const float log10_HII_over_nH = interpolation_4d_no_z_no_w(
          cooling->table.logHfracs_Teq,                       /* */
          ired, imet, colibre_cooling_N_density - 1, ionized, /* */
          dred, dmet, 0., 0.,                                 /* */
          colibre_cooling_N_redshifts,                        /* */
          colibre_cooling_N_metallicity,                      /* */
          colibre_cooling_N_density,                          /* */
          3);

      return exp10f(log10_HII_over_nH);

    } else {

      /* Normal case: We are not beyond the table range */

      /* Need to find thermal equilibrium state with the same pressure
       * by interpolating the equilibrium table
       *
       * Note that the logPeq table is neither equally spaced nor
       * necessarilly monotically increasing.
       * We hance loop over densities and pick the first one where
       * log10_P < log10_Peq. We start with the resolved density (index iden),
       * as we expect the subgrid density to be larger */

      /* If the solution can't be found, we revert to the non-subgrid
       * density */
      int iden_eq = iden;
      float dden_eq = dden;

      /* Loop over the density bins */
      for (int i = iden; i < colibre_cooling_N_density; i++) {

        /* Equilibirum pressure at this density */
        const float log10_Peq_interp =
            interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                  ired, imet, i,                 /* */
                                  dred, dmet, 0.,                /* */
                                  colibre_cooling_N_redshifts,   /* */
                                  colibre_cooling_N_metallicity, /* */
                                  colibre_cooling_N_density);

        /* Did we find a solution? */
        if (log10_Peq_interp > log10_P_cgs) {

          /* Equilibrium pressure at the previous density point */
          const float log10_Peq_prev =
              interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                    ired, imet, i - 1,             /* */
                                    dred, dmet, 0.,                /* */
                                    colibre_cooling_N_redshifts,   /* */
                                    colibre_cooling_N_metallicity, /* */
                                    colibre_cooling_N_density);

          /* How far from the equilibrium point are we? */
          const float delta_P_eq = (log10_P_cgs - log10_Peq_prev) /
                                   (log10_Peq_interp - log10_Peq_prev);

          /* Interpolate to get the density at equilibrium */
          const float log10_n_at_Peq =
              cooling->nH[i - 1] +
              delta_P_eq * (cooling->nH[i] - cooling->nH[i - 1]);

          /* We found a valid density: Get the index along the density
           * axis of the tables. */
          get_index_1d(cooling->nH, colibre_cooling_N_density, log10_n_at_Peq,
                       &iden_eq, &dden_eq);

          break;
        }
      }

      /* Finish by interpolating the tables to get the HII density */
      const float log10_HII_over_nH =
          interpolation_4d_no_w(cooling->table.logHfracs_Teq,  /* */
                                ired, imet, iden_eq, ionized,  /* */
                                dred, dmet, dden_eq, 0.,       /* */
                                colibre_cooling_N_redshifts,   /* */
                                colibre_cooling_N_metallicity, /* */
                                colibre_cooling_N_density,     /* */
                                3);

      return exp10f(log10_HII_over_nH);
    }

  } else {

    /* We are above the EoS. */
    const float weights[3] = {1.0, 1.0, 1.0};
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

    return nHII_over_nH;
  }
}

/**
 * @brief Compute the fraction of Hydrogen that is in H_2 based
 * on the pressure of the gas.
 *
 * For particles on the entropy floor, we use pressure equilibrium to
 * infer the properties of the particle.
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param p The #part.
 * @param xp The #xpart.
 * @param log10_T The logarithm base 10 of the temperature of the particle.
 * @param log10_T_EOS_max The logarithm base 10 of the maximal temperature
 * to be considered on the EOS at the density of this particle.
 */
double compute_subgrid_H2_fraction(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const struct part *p, const struct xpart *xp,
    const float log10_T, const float log10_T_EOS_max) {

  /* Physical density of this particle */
  const double rho_phys = hydro_get_physical_density(p, cosmo);

  /* Get the total metallicity */
  float dummy[colibre_cooling_N_elementtypes];
  const float logZZsol = abundance_ratio_to_solar(p, cooling, dummy);

  /* Get the Hydrogen number density */
  const float *const metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];
  const double n_H = XH * rho_phys / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis for this particle */
  int ired, imet, iden, item;
  float dred, dmet, dden, dtem;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z, &ired,
               &dred);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &iden,
               &dden);
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
               &dtem);

  /* Are we above or below the EoS? */
  if (log10_T < log10_T_EOS_max) {

    /* We are below the EoS. Use subgrid properties assuming P equilibrium */

    const float P = hydro_get_physical_pressure(p, cosmo);
    const float P_cgs = P * cooling->pressure_to_cgs;
    const float log10_P_cgs = log10(P_cgs);

    /* Recover the maximal equilibrium pressure from the table at the current
     * redshift and metallicity */
    const float log10_Peq_max_cgs =
        interpolation_3d_no_z(cooling->table.logPeq,                     /* */
                              ired, imet, colibre_cooling_N_density - 1, /* */
                              dred, dmet, 0.,                            /* */
                              colibre_cooling_N_redshifts,               /* */
                              colibre_cooling_N_metallicity,             /* */
                              colibre_cooling_N_density);

    /* Are we above or below the max pressure equilibrium? */
    if (log10_P_cgs > log10_Peq_max_cgs) {

      /* EOS pressure (logP) is larger than maximum Peq (can happen for very
       * steep EOS) H2 fraction from the highest density bin */

      const float log10_H2_over_nH = interpolation_4d_no_z_no_w(
          cooling->table.logHfracs_Teq,                         /* */
          ired, imet, colibre_cooling_N_density - 1, molecular, /* */
          dred, dmet, 0., 0.,                                   /* */
          colibre_cooling_N_redshifts,                          /* */
          colibre_cooling_N_metallicity,                        /* */
          colibre_cooling_N_density,                            /* */
          3);

      return 0.5f * exp10f(log10_H2_over_nH);

    } else {

      /* Normal case: We are not beyond the table range */

      /* Need to find thermal equilibrium state with the same pressure
       * by interpolating the equilibrium table
       *
       * Note that the logPeq table is neither equally spaced nor
       * necessarilly monotically increasing.
       * We hance loop over densities and pick the first one where
       * log10_P < log10_Peq. We start with the resolved density (index iden),
       * as we expect the subgrid density to be larger */

      /* If the solution can't be found, we revert to the non-subgrid
       * density */
      int iden_eq = iden;
      float dden_eq = dden;

      /* Loop over the density bins */
      for (int i = iden; i < colibre_cooling_N_density; i++) {

        /* Equilibirum pressure at this density */
        const float log10_Peq_interp =
            interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                  ired, imet, i,                 /* */
                                  dred, dmet, 0.,                /* */
                                  colibre_cooling_N_redshifts,   /* */
                                  colibre_cooling_N_metallicity, /* */
                                  colibre_cooling_N_density);

        /* Did we find a solution? */
        if (log10_Peq_interp > log10_P_cgs) {

          /* Equilibrium pressure at the previous density point */
          const float log10_Peq_prev =
              interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                    ired, imet, i - 1,             /* */
                                    dred, dmet, 0.,                /* */
                                    colibre_cooling_N_redshifts,   /* */
                                    colibre_cooling_N_metallicity, /* */
                                    colibre_cooling_N_density);

          /* How far from the equilibrium point are we? */
          const float delta_P_eq = (log10_P_cgs - log10_Peq_prev) /
                                   (log10_Peq_interp - log10_Peq_prev);

          /* Interpolate to get the density at equilibrium */
          const float log10_n_at_Peq =
              cooling->nH[i - 1] +
              delta_P_eq * (cooling->nH[i] - cooling->nH[i - 1]);

          /* We found a valid density: Get the index along the density
           * axis of the tables. */
          get_index_1d(cooling->nH, colibre_cooling_N_density, log10_n_at_Peq,
                       &iden_eq, &dden_eq);

          break;
        }
      }

      /* Finish by interpolating the tables to get the H2 density */
      const float log10_H2_over_nH =
          interpolation_4d_no_w(cooling->table.logHfracs_Teq,   /* */
                                ired, imet, iden_eq, molecular, /* */
                                dred, dmet, dden_eq, 0.,        /* */
                                colibre_cooling_N_redshifts,    /* */
                                colibre_cooling_N_metallicity,  /* */
                                colibre_cooling_N_density,      /* */
                                3);

      return 0.5 * exp10f(log10_H2_over_nH);
    }

  } else {

    /* We are above the EoS. */
    const float weights[3] = {1.0, 1.0, 1.0};
    const float nHII_over_nH =
        interpolation4d_plus_summation(cooling->table.logHfracs_all,  /* */
                                       weights, molecular, molecular, /* */
                                       ired, item, imet, iden,        /* */
                                       dred, dtem, dmet, dden,        /* */
                                       colibre_cooling_N_redshifts,   /* */
                                       colibre_cooling_N_temperature, /* */
                                       colibre_cooling_N_metallicity, /* */
                                       colibre_cooling_N_density,     /* */
                                       3);

    return 0.5 * nHII_over_nH;
  }
}

/**
 * @brief Compute the subgrid temperature of the gas.
 *
 * For particles on the entropy floor, we use pressure equilibrium to
 * infer the properties of the particle.
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param p The #part.
 * @param xp The #xpart.
 * @param log10_T The logarithm base 10 of the temperature of the particle.
 * @param log10_T_EOS_max The logarithm base 10 of the maximal temperature
 * to be considered on the EOS at the density of this particle.
 */
double compute_subgrid_temperature(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const struct part *p, const struct xpart *xp,
    const float log10_T, const float log10_T_EOS_max) {

  /* Physical density of this particle */
  const double rho_phys = hydro_get_physical_density(p, cosmo);

  /* Get the total metallicity */
  float dummy[colibre_cooling_N_elementtypes];
  const float logZZsol = abundance_ratio_to_solar(p, cooling, dummy);

  /* Get the Hydrogen number density */
  const float *const metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];
  const double n_H = XH * rho_phys / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis for this particle */
  int ired, imet, iden, item;
  float dred, dmet, dden, dtem;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z, &ired,
               &dred);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &iden,
               &dden);
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
               &dtem);

  /* Are we above or below the EoS? */
  if (log10_T < log10_T_EOS_max) {

    /* We are below the EoS. Use subgrid properties assuming P equilibrium */

    const float P = hydro_get_physical_pressure(p, cosmo);
    const float P_cgs = P * cooling->pressure_to_cgs;
    const float log10_P_cgs = log10(P_cgs);

    /* Recover the maximal equilibrium pressure from the table at the current
     * redshift and metallicity */
    const float log10_Peq_max_cgs =
        interpolation_3d_no_z(cooling->table.logPeq,                     /* */
                              ired, imet, colibre_cooling_N_density - 1, /* */
                              dred, dmet, 0.,                            /* */
                              colibre_cooling_N_redshifts,               /* */
                              colibre_cooling_N_metallicity,             /* */
                              colibre_cooling_N_density);

    /* Are we above or below the max pressure equilibrium? */
    if (log10_P_cgs > log10_Peq_max_cgs) {

      /* EOS pressure (logP) is larger than maximum Peq (can happen for very
       * steep EOS) HII fraction from the highest density bin */

      const float log10_T_at_Peq =
          interpolation_3d_no_z(cooling->table.logTeq,                     /* */
                                ired, imet, colibre_cooling_N_density - 1, /* */
                                dred, dmet, 0.,                            /* */
                                colibre_cooling_N_redshifts,               /* */
                                colibre_cooling_N_metallicity,             /* */
                                colibre_cooling_N_density);

      const float T_at_Peq = exp10f(log10_T_at_Peq);

      return T_at_Peq;

    } else {

      /* Normal case: We are not beyond the table range */

      /* Need to find thermal equilibrium state with the same pressure
       * by interpolating the equilibrium table
       *
       * Note that the logPeq table is neither equally spaced nor
       * necessarilly monotically increasing.
       * We hance loop over densities and pick the first one where
       * log10_P < log10_Peq. We start with the resolved density (index iden),
       * as we expect the subgrid density to be larger */

      /* If the solution can't be found, we revert to the non-subgrid
       * density */
      int iden_eq = iden;
      float dden_eq = dden;

      /* Loop over the density bins */
      for (int i = iden; i < colibre_cooling_N_density; i++) {

        /* Equilibirum pressure at this density */
        const float log10_Peq_interp =
            interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                  ired, imet, i,                 /* */
                                  dred, dmet, 0.,                /* */
                                  colibre_cooling_N_redshifts,   /* */
                                  colibre_cooling_N_metallicity, /* */
                                  colibre_cooling_N_density);

        /* Did we find a solution? */
        if (log10_Peq_interp > log10_P_cgs) {

          /* Equilibrium pressure at the previous density point */
          const float log10_Peq_prev =
              interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                    ired, imet, i - 1,             /* */
                                    dred, dmet, 0.,                /* */
                                    colibre_cooling_N_redshifts,   /* */
                                    colibre_cooling_N_metallicity, /* */
                                    colibre_cooling_N_density);

          /* How far from the equilibrium point are we? */
          const float delta_P_eq = (log10_P_cgs - log10_Peq_prev) /
                                   (log10_Peq_interp - log10_Peq_prev);

          /* Interpolate to get the density at equilibrium */
          const float log10_n_at_Peq =
              cooling->nH[i - 1] +
              delta_P_eq * (cooling->nH[i] - cooling->nH[i - 1]);

          /* We found a valid density: Get the index along the density
           * axis of the tables. */
          get_index_1d(cooling->nH, colibre_cooling_N_density, log10_n_at_Peq,
                       &iden_eq, &dden_eq);

          break;
        }
      }

      /* Finish by interpolating the tables to get the temperature */
      const float log10_T_at_Peq =
          interpolation_3d(cooling->table.logTeq,         /* */
                           ired, imet, iden_eq,           /* */
                           dred, dmet, dden_eq,           /* */
                           colibre_cooling_N_redshifts,   /* */
                           colibre_cooling_N_metallicity, /* */
                           colibre_cooling_N_density);

      const float T_at_Peq = exp10f(log10_T_at_Peq);

      return T_at_Peq;
    }
  } else {

    /* We are above the EoS. */
    /* Use the normal temperature */
    return exp10(log10_T);
  }
}

/**
 * @brief Compute the physical subgrid density of the gas.
 *
 * For particles on the entropy floor, we use pressure equilibrium to
 * infer the properties of the particle.
 *
 * Note that we return the density in physical coordinates.
 *
 * @param cooling The properties of the cooling scheme.
 * @param phys_const The physical constants.
 * @param floor_props The properties of the entropy floor.
 * @param cosmo The cosmological model.
 * @param p The #part.
 * @param xp The #xpart.
 * @param log10_T The logarithm base 10 of the temperature of the particle.
 * @param log10_T_EOS_max The logarithm base 10 of the maximal temperature
 * to be considered on the EOS at the density of this particle.
 */
double compute_subgrid_density(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const struct part *p, const struct xpart *xp,
    const float log10_T, const float log10_T_EOS_max) {

  /* Physical density of this particle */
  const double rho_phys = hydro_get_physical_density(p, cosmo);

  /* Get the total metallicity */
  float dummy[colibre_cooling_N_elementtypes];
  const float logZZsol = abundance_ratio_to_solar(p, cooling, dummy);

  /* Get the Hydrogen number density */
  const float *const metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];
  const double n_H = XH * rho_phys / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis for this particle */
  int ired, imet, iden, item;
  float dred, dmet, dden, dtem;
  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z, &ired,
               &dred);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &iden,
               &dden);
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10_T, &item,
               &dtem);

  /* Are we above or below the EoS? */
  if (log10_T < log10_T_EOS_max) {

    /* We are below the EoS. Use subgrid properties assuming P equilibrium */

    const float P = hydro_get_physical_pressure(p, cosmo);
    const float P_cgs = P * cooling->pressure_to_cgs;
    const float log10_P_cgs = log10(P_cgs);
    float log10_n_at_Peq;

    /* Recover the maximal equilibrium pressure from the table at the current
     * redshift and metallicity */
    const float log10_Peq_max_cgs =
        interpolation_3d_no_z(cooling->table.logPeq,                     /* */
                              ired, imet, colibre_cooling_N_density - 1, /* */
                              dred, dmet, 0.,                            /* */
                              colibre_cooling_N_redshifts,               /* */
                              colibre_cooling_N_metallicity,             /* */
                              colibre_cooling_N_density);

    /* Are we above or below the max pressure equilibrium? */
    if (log10_P_cgs > log10_Peq_max_cgs) {

      /* EOS pressure (logP) is larger than maximum Peq (can happen for very
       * steep EOS); use the highest density bin */

      const float log10_T_at_Peq = interpolation_3d_no_z(
          cooling->table.logTeq, ired, imet, colibre_cooling_N_density - 1,
          dred, dmet, 0., colibre_cooling_N_redshifts,
          colibre_cooling_N_metallicity, colibre_cooling_N_density);

      const float mu_at_Peq = interpolation_3d_no_z(
          cooling->table.meanpartmass_Teq, ired, imet,
          colibre_cooling_N_density - 1, dred, dmet, 0.,
          colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
          colibre_cooling_N_density);

      const double log10_kB = cooling->log10_kB_cgs;

      log10_n_at_Peq = log10_P_cgs - log10_T_at_Peq + log10(XH) +
                       log10(mu_at_Peq) - log10_kB;

    } else {

      /* Normal case: We are not beyond the table range */

      /* Need to find thermal equilibrium state with the same pressure
       * by interpolating the equilibrium table
       *
       * Note that the logPeq table is neither equally spaced nor
       * necessarilly monotically increasing.
       * We hance loop over densities and pick the first one where
       * log10_P < log10_Peq. We start with the resolved density (index iden),
       * as we expect the subgrid density to be larger */

      /* If the solution can't be found, we revert to the non-subgrid
       * density */
      log10_n_at_Peq = log10(n_H_cgs);

      /* Loop over the density bins */
      for (int i = iden; i < colibre_cooling_N_density; i++) {

        /* Equilibirum pressure at this density */
        const float log10_Peq_interp =
            interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                  ired, imet, i,                 /* */
                                  dred, dmet, 0.,                /* */
                                  colibre_cooling_N_redshifts,   /* */
                                  colibre_cooling_N_metallicity, /* */
                                  colibre_cooling_N_density);

        /* Did we find a solution? */
        if (log10_Peq_interp > log10_P_cgs) {

          /* Equilibrium pressure at the previous density point */
          const float log10_Peq_prev =
              interpolation_3d_no_z(cooling->table.logPeq,         /* */
                                    ired, imet, i - 1,             /* */
                                    dred, dmet, 0.,                /* */
                                    colibre_cooling_N_redshifts,   /* */
                                    colibre_cooling_N_metallicity, /* */
                                    colibre_cooling_N_density);

          /* How far from the equilibrium point are we? */
          const float delta_P_eq = (log10_P_cgs - log10_Peq_prev) /
                                   (log10_Peq_interp - log10_Peq_prev);

          /* Interpolate to get the density at equilibrium */
          log10_n_at_Peq = cooling->nH[i - 1] +
                           delta_P_eq * (cooling->nH[i] - cooling->nH[i - 1]);

          break;
        }
      }
    }

    const float n_at_Peq = exp10(log10_n_at_Peq);
    return n_at_Peq * cooling->number_density_from_cgs *
           phys_const->const_proton_mass / XH;

  } else {

    /* We are above the EoS. */

    return rho_phys;
  }
}

#endif /* SWIFT_COLIBRE_COOLING_SUBGRID_H */
