/**
 * @file Backup file containing the now-defunct Newton integration scheme
 * for the cooling.
 */

/**
 * @brief Computes the log_10 of the temperature corresponding to a given
 * internal energy, hydrogen number density, Helium fraction and redshift.
 *
 * Note that the redshift is implicitly passed in via the currently loaded
 * tables in the #cooling_function_data.
 *
 * For the low-z case, we interpolate the flattened 4D table 'u_to_temp' that
 * is arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * For the high-z case, we interpolate the flattened 3D table 'u_to_temp' that
 * is arranged in the following way:
 * - 1st dim: Hydrogen density, length = eagle_cooling_N_density
 * - 2nd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 3rd dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * @param log_10_u_cgs Log base 10 of internal energy in cgs.
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param He_index Index along the Helium fraction dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param d_He Offset between helium fraction and table[He_index].
 * @param cooling #cooling_function_data structure.
 *
 * @param compute_dT_du Do we want to compute dT/du ?
 * @param dT_du (return) The value of dT/du
 *
 * @return log_10 of the temperature.
 */
__attribute__((always_inline)) INLINE double eagle_convert_u_to_temp(
    const double log_10_u_cgs, const float redshift, const int compute_dT_du,
    float *dT_du, int n_H_index, int He_index, float d_n_H, float d_He,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the internal energy axis */
  int u_index;
  float d_u;
  get_index_1d(cooling->Therm, eagle_cooling_N_temperature, log_10_u_cgs,
               &u_index, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  float log_10_T;
  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    log_10_T = interpolation_3d(cooling->table.temperature,   /* */
                                n_H_index, He_index, u_index, /* */
                                d_n_H, d_He, d_u,             /* */
                                eagle_cooling_N_density,      /* */
                                eagle_cooling_N_He_frac,      /* */
                                eagle_cooling_N_temperature); /* */
  } else {

    log_10_T =
        interpolation_4d(cooling->table.temperature,                  /* */
                         /*z_index=*/0, n_H_index, He_index, u_index, /* */
                         cooling->dz, d_n_H, d_He, d_u,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */
  }

  if (compute_dT_du) {

    float log_10_T_high, log_10_T_low;

    /* Interpolate temperature table to return temperature for internal energy
     * at grid point above current internal energy for computing dT_du used for
     * calculation of dlambda_du in cooling.c (use 3D interpolation for high
     * redshift table, otherwise 4D) */
    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      log_10_T_high = interpolation_3d(cooling->table.temperature,   /* */
                                       n_H_index, He_index, u_index, /* */
                                       d_n_H, d_He, /*delta_u=*/1.f, /* */
                                       eagle_cooling_N_density,      /* */
                                       eagle_cooling_N_He_frac,      /* */
                                       eagle_cooling_N_temperature); /* */

    } else {

      log_10_T_high =
          interpolation_4d(cooling->table.temperature,                  /* */
                           /*z_index=*/0, n_H_index, He_index, u_index, /* */
                           cooling->dz, d_n_H, d_He, /*delta_u=*/1.f,   /* */
                           eagle_cooling_N_loaded_redshifts,            /* */
                           eagle_cooling_N_density,                     /* */
                           eagle_cooling_N_He_frac,                     /* */
                           eagle_cooling_N_temperature);                /* */
    }

    /* Interpolate temperature table to return temperature for internal energy
     * at grid point below current internal energy for computing dT_du used for
     * calculation of dlambda_du in cooling.c (use 3D interpolation for high
     * redshift table, otherwise 4D) */
    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      log_10_T_low = interpolation_3d(cooling->table.temperature,   /* */
                                      n_H_index, He_index, u_index, /* */
                                      d_n_H, d_He, /*delta_u=*/0.f, /* */
                                      eagle_cooling_N_density,      /* */
                                      eagle_cooling_N_He_frac,      /* */
                                      eagle_cooling_N_temperature); /* */

    } else {

      log_10_T_low =
          interpolation_4d(cooling->table.temperature,                  /* */
                           /*z_index=*/0, n_H_index, He_index, u_index, /* */
                           cooling->dz, d_n_H, d_He, /*delta_u=*/0.f,   /* */
                           eagle_cooling_N_loaded_redshifts,            /* */
                           eagle_cooling_N_density,                     /* */
                           eagle_cooling_N_He_frac,                     /* */
                           eagle_cooling_N_temperature);                /* */
    }

    /* Calculate dT/du */
    const float delta_u = exp(cooling->Therm[u_index + 1] * M_LN10) -
                          exp(cooling->Therm[u_index] * M_LN10);
    *dT_du =
        (exp(M_LN10 * log_10_T_high) - exp(M_LN10 * log_10_T_low)) / delta_u;
  }

  /* Special case for temperatures below the start of the table */
  if (u_index == 0 && d_u == 0.f) {

    /* The temperature is multiplied by u / 10^T[0]
     * where T[0] is the first entry in the table */
    log_10_T += log_10_u_cgs - cooling->Temp[0];
  }

  return log_10_T;
}

/**
 * @brief Computes the cooling rate corresponding to a given internal energy,
 * hydrogen number density, Helium fraction, redshift and metallicity from
 * all the possible channels.
 *
 * 1) Metal-free cooling:
 * We interpolate the flattened 4D table 'H_and_He_net_heating' that is
 * arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 2) Electron abundance
 * We compute the electron abundance by interpolating the flattened 4d table
 * 'H_and_He_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 3) Compton cooling is applied via the analytic formula.
 *
 * 4) Solar electron abudance
 * We compute the solar electron abundance by interpolating the flattened 3d
 * table 'solar_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 5) Metal-line cooling
 * For each tracked element we interpolate the flattened 4D table
 * 'table_metals_net_heating' that is arrange in the following way:
 * - 1st dim: element, length = eagle_cooling_N_metal
 * - 2nd dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * Note that this is a fake 4D interpolation as we do not interpolate
 * along the 1st dimension. We just do this once per element.
 *
 * Since only the temperature changes when cooling a given particle,
 * the redshift, hydrogen number density and helium fraction indices
 * and offsets passed in.
 *
 * If the arguement dlambda_du is non-NULL, the routine also
 * calculates derivative of cooling rate with respect to internal
 * energy.
 *
 * If the argument element_lambda is non-NULL, the routine also
 * returns the cooling rate per element in the array.
 *
 * @param log10_u_cgs Log base 10 of internal energy per unit mass in CGS units.
 * @param redshift The current redshift
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param solar_ratio Array of ratios of particle metal abundances
 * to solar metal abundances
 *
 * @param n_H_index Particle hydrogen number density index
 * @param d_n_H Particle hydrogen number density offset
 * @param He_index Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling Cooling data structure
 *
 * @param dlambda_du (return) Derivative of the cooling rate with respect to u.
 * @param element_lambda (return) Cooling rate from each element
 *
 * @return The cooling rate
 */
INLINE static double eagle_metal_cooling_rate(
    double log10_u_cgs, double redshift, double n_H_cgs,
    const float solar_ratio[eagle_cooling_N_abundances], int n_H_index,
    float d_n_H, int He_index, float d_He,
    const struct cooling_function_data *restrict cooling, double *dlambda_du,
    double *element_lambda) {

  /* used for calculating dlambda_du */
  double temp_lambda_high = 0, temp_lambda_low = 0;
  double h_plus_he_electron_abundance_high = 0;
  double h_plus_he_electron_abundance_low = 0;
  double solar_electron_abundance_high = 0;
  double solar_electron_abundance_low = 0;
  double elem_cool_low = 0, elem_cool_high = 0;

  /* We only need dT_du if dLambda_du is non-NULL */
  const int compute_dT_du = (dlambda_du != NULL) ? 1 : 0;

  /* Temperature */
  float dT_du = -1.f;
  const double log_10_T =
      eagle_convert_u_to_temp(log10_u_cgs, redshift, compute_dT_du, &dT_du,
                              n_H_index, He_index, d_n_H, d_He, cooling);

  /* Get index along temperature dimension of the tables */
  int T_index;
  float d_T;
  get_index_1d(cooling->Temp, eagle_cooling_N_temperature, log_10_T, &T_index,
               &d_T);

  /* Difference between entries on the temperature table around u */
  const float delta_T = exp(M_LN10 * cooling->Temp[T_index + 1]) -
                        exp(M_LN10 * cooling->Temp[T_index]);

  /**********************/
  /* Metal-free cooling */
  /**********************/

  double Lambda_free;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    Lambda_free = interpolation_3d(cooling->table.H_plus_He_heating, /* */
                                   n_H_index, He_index, T_index,     /* */
                                   d_n_H, d_He, d_T,                 /* */
                                   eagle_cooling_N_density,          /* */
                                   eagle_cooling_N_He_frac,          /* */
                                   eagle_cooling_N_temperature);     /* */

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du. Pass in NULL pointer for
     * dlambda_du in order to skip */
    if (dlambda_du != NULL) {
      temp_lambda_high = interpolation_3d(
          cooling->table.H_plus_He_heating, n_H_index, He_index, T_index, d_n_h,
          d_He, 1.f, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low = interpolation_3d(
          cooling->table.H_plus_He_heating, n_H_index, He_index, T_index, d_n_h,
          d_He, 0.f, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }

  } else {

    /* Using normal tables, have to interpolate in redshift */
    Lambda_free =
        interpolation_4d(cooling->table.H_plus_He_heating,            /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_H, d_He, d_T,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      temp_lambda_high =
          interpolation_4d(cooling->table.H_plus_He_heating, 0, n_H_index,
                           He_index, T_index, cooling->dz, d_n_h, d_He, 1.f, 2,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low =
          interpolation_4d(cooling->table.H_plus_He_heating, 0, n_H_index,
                           He_index, T_index, cooling->dz, d_n_h, d_He, 0.f, 2,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }
  }

  if (dlambda_du != NULL) {
    *dlambda_du += (temp_lambda_high - temp_lambda_low) / delta_T * dT_du;
  }

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[0] = Lambda_free;
  }

  /**********************/
  /* Electron abundance */
  /**********************/

  double H_plus_He_electron_abundance;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    H_plus_He_electron_abundance =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, /* */
                         n_H_index, He_index, T_index,                /* */
                         d_n_H, d_He, d_T,                            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du. Pass in NULL pointer for
     * dlambda_du in order to skip */

    h_plus_he_electron_abundance_high =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, n_H_index,
                         He_index, T_index, d_n_h, d_He, 1.f, cooling->N_nH,
                         cooling->N_He, cooling->N_Temp);
    h_plus_he_electron_abundance_low =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, n_H_index,
                         He_index, T_index, d_n_h, d_He, 0.f, cooling->N_nH,
                         cooling->N_He, cooling->N_Temp);

  } else {

    H_plus_He_electron_abundance =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_H, d_He, d_T,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    h_plus_he_electron_abundance_high =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, 0,
                         n_H_index, He_index, T_index, cooling->dz, d_n_h, d_He,
                         1.f, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    h_plus_he_electron_abundance_low =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, 0,
                         n_H_index, He_index, T_index, cooling->dz, d_n_h, d_He,
                         0.f, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
  }

  /**********************/
  /* Compton cooling    */
  /**********************/

  double Lambda_Compton = 0.;

  /* Do we need to add the inverse Compton cooling? */
  /* It is *not* stored in the tables before re-ionisation */
  if ((redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) ||
      (redshift > cooling->H_reion_z)) {

    const double T = exp10(log_10_T);

    /* Note the minus sign */
    Lambda_Compton -= eagle_Compton_cooling_rate(cooling, redshift, n_H_cgs, T,
                                                 H_plus_He_electron_abundance);
  }

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[1] = Lambda_Compton;
  }

  /*******************************/
  /* Solar electron abundance    */
  /*******************************/

  double solar_electron_abundance;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    solar_electron_abundance =
        interpolation_2d(cooling->table.electron_abundance, /* */
                         n_H_index, T_index,                /* */
                         d_n_H, d_T,                        /* */
                         eagle_cooling_N_density,           /* */
                         eagle_cooling_N_temperature);      /* */

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high =
          interpolation_2d(cooling->table.electron_abundance, n_H_index,
                           T_index, d_n_h, 1.f, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low =
          interpolation_2d(cooling->table.electron_abundance, n_H_index,
                           T_index, d_n_h, 0.f, cooling->N_nH, cooling->N_Temp);
    }

  } else {

    /* Using normal tables, have to interpolate in redshift */
    solar_electron_abundance =
        interpolation_3d(cooling->table.electron_abundance, /* */
                         /*z_index=*/0, n_H_index, T_index, /* */
                         cooling->dz, d_n_H, d_T,           /* */
                         eagle_cooling_N_loaded_redshifts,  /* */
                         eagle_cooling_N_density,           /* */
                         eagle_cooling_N_temperature);      /* */

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high = interpolation_3d(
          cooling->table.electron_abundance, 0, n_H_index, T_index, cooling->dz,
          d_n_h, 1.f, 2, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low = interpolation_3d(
          cooling->table.electron_abundance, 0, n_H_index, T_index, cooling->dz,
          d_n_h, 0.f, 2, cooling->N_nH, cooling->N_Temp);
    }

    const double electron_abundance_ratio =
        H_plus_He_electron_abundance / solar_electron_abundance;

    /**********************/
    /* Metal-line cooling */
    /**********************/

    /* for each element the cooling rate is multiplied by the ratio of H, He
     * electron abundance to solar electron abundance then by the ratio of the
     * particle metal abundance to solar metal abundance. */

    double lambda_metal[eagle_cooling_N_metal + 2] = {0.};

    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      /* Loop over the metals (ignore H and He) */
      for (int elem = 2; elem < eagle_cooling_N_metal + 2; elem++) {

        if (solar_ratio[elem] > 0.) {

          /* Note that we do not interpolate along the x-axis
           * (element dimension) */
          lambda_metal[elem] =
              interpolation_3d_no_x(cooling->table.metal_heating,   /* */
                                    elem - 2, n_H_index, T_index,   /* */
                                    /*delta_elem=*/0.f, d_n_H, d_T, /* */
                                    eagle_cooling_N_metal,          /* */
                                    eagle_cooling_N_density,        /* */
                                    eagle_cooling_N_temperature);   /* */

          lambda_metal[elem] *= electron_abundance_ratio;
          lambda_metal[elem] *= solar_ratio[elem];
        }

        /* compute values at temperature gridpoints above and below input
         * temperature for calculation of dlambda_du */
        if (dlambda_du != NULL) {
          elem_cool_high = interpolation_3d_no_x(
              cooling->table.metal_heating, elem, n_H_index, T_index, 0.f,
              d_n_h, 1.f, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);

          elem_cool_low = interpolation_3d_no_x(
              cooling->table.metal_heating, elem, n_H_index, T_index, 0.f,
              d_n_h, 0.f, cooling->N_nH, cooling->N_Temp, cooling->N_Elements);

          *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                              solar_electron_abundance_high -
                          elem_cool_low * h_plus_he_electron_abundance_low /
                              solar_electron_abundance_low) /
                         delta_T * dT_du * solar_ratio[elem + 2];
        }
      }

    } else {

      /* Loop over the metals (ignore H and He) */
      for (int elem = 2; elem < eagle_cooling_N_metal + 2; elem++) {

        if (solar_ratio[elem] > 0.) {

          /* Note that we do not interpolate along the x-axis
           * (element dimension) */
          lambda_metal[elem] = interpolation_4d_no_x(
              cooling->table.metal_heating,                /* */
              elem - 2, /*z_index=*/0, n_H_index, T_index, /* */
              /*delta_elem=*/0.f, cooling->dz, d_n_H, d_T, /* */
              eagle_cooling_N_metal,                       /* */
              eagle_cooling_N_loaded_redshifts,            /* */
              eagle_cooling_N_density,                     /* */
              eagle_cooling_N_temperature);                /* */

          lambda_metal[elem] *= electron_abundance_ratio;
          lambda_metal[elem] *= solar_ratio[elem];
        }

        /* compute values at temperature gridpoints above and below input
         * temperature for calculation of dlambda_du */
        if (dlambda_du != NULL) {
          elem_cool_high = interpolation_4d_no_x(
              cooling->table.metal_heating, elem, 0, n_H_index, T_index, 0.,
              cooling->dz, d_n_h, 1.f, cooling->N_Elements, 2, cooling->N_nH,
              cooling->N_Temp);

          elem_cool_low = interpolation_4d_no_x(
              cooling->table.metal_heating, elem, 0, n_H_index, T_index, 0.,
              cooling->dz, d_n_h, 0.f, cooling->N_Elements, 2, cooling->N_nH,
              cooling->N_Temp);

          *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                              solar_electron_abundance_high -
                          elem_cool_low * h_plus_he_electron_abundance_low /
                              solar_electron_abundance_low) /
                         delta_T * dT_du * solar_ratio[elem + 2];
        }
      }
    }

    if (element_lambda != NULL) {
      for (int elem = 2; elem < eagle_cooling_N_metal + 2; ++elem) {
        element_lambda[elem] = lambda_metal[elem];
      }
    }

    /* Sum up all the contributions */
    double Lambda_net = Lambda_free + Lambda_Compton;
    for (int elem = 2; elem < eagle_cooling_N_metal + 2; ++elem) {
      Lambda_net += lambda_metal[elem];
    }

    return Lambda_net;
  }

  /**
   * @brief Wrapper function used to calculate cooling rate and dLambda_du.
   * Table indices and offsets for redshift, hydrogen number density and
   * helium fraction are passed it so as to compute them only once per particle.
   *
   * @param log_u_cgs Natural log of internal energy per unit mass in CGS units.
   * @param redshift The current redshift.
   * @param n_H_cgs Hydrogen number density in CGS units.
   * @param abundance_ratio Ratio of element abundance to solar.
   *
   * @param n_H_index Particle hydrogen number density index
   * @param d_n_H Particle hydrogen number density offset
   * @param He_index Particle helium fraction index
   * @param d_He Particle helium fraction offset
   * @param cooling #cooling_function_data structure
   *
   * @param dLambdaNet_du (return) Derivative of the cooling rate with respect
   * to u.
   *
   * @return The cooling rate
   */
  INLINE static double eagle_cooling_rate(
      double log_u_cgs, double redshift, double n_H_cgs,
      const float abundance_ratio[eagle_cooling_N_abundances], int n_H_index,
      float d_n_H, int He_index, float d_He,
      const struct cooling_function_data *restrict cooling,
      double *dLambdaNet_du) {

    return eagle_metal_cooling_rate(log_u_cgs / M_LN10, redshift, n_H_cgs,
                                    abundance_ratio, n_H_index, d_n_H, He_index,
                                    d_He, cooling, dLambdaNet_du,
                                    /*element_lambda=*/NULL);
  }

  /**
   * @brief Newton Raphson integration scheme to calculate particle cooling over
   * timestep. This replaces bisection scheme used in EAGLE to minimize the
   * number of array accesses. Integration defaults to bisection scheme (see
   * function bisection_iter) if this function does not converge within a
   * specified number of steps
   *
   * @param logu_init Initial guess for log(internal energy)
   * @param u_ini Internal energy at beginning of hydro step
   * @param n_H_index Particle hydrogen number density index
   * @param d_n_H Particle hydrogen number density offset
   * @param He_index Particle helium fraction index
   * @param d_He Particle helium fraction offset
   * @param He_reion_heat Heating due to helium reionization
   * (only depends on redshift, so passed as parameter)
   * @param p #part structure
   * @param cosmo #cosmology structure
   * @param cooling #cooling_function_data structure
   * @param phys_const #phys_const data structure
   * @param abundance_ratio Array of ratios of metal abundance to solar
   * @param dt timestep
   * @param bisection_flag Flag to identify if scheme failed to converge
   */
  float newton_iter(float logu_init, double u_ini, int n_H_index, float d_n_H,
                    int He_index, float d_He, float He_reion_heat,
                    struct part *restrict p,
                    const struct cosmology *restrict cosmo,
                    const struct cooling_function_data *restrict cooling,
                    const struct phys_const *restrict phys_const,
                    const float abundance_ratio[eagle_cooling_N_abundances],
                    float dt, int *bisection_flag) {

    double logu, logu_old;
    double dLambdaNet_du = 0.0, LambdaNet;

    /* table bounds */
    const float log_table_bound_high =
        (cooling->Therm[eagle_cooling_N_temperature - 1] - 0.05) / M_LOG10E;
    const float log_table_bound_low = (cooling->Therm[0] + 0.05) / M_LOG10E;

    /* convert Hydrogen mass fraction in Hydrogen number density */
    const float XH =
        p->chemistry_data.smoothed_metal_mass_fraction[chemistry_element_H];
    const double n_H = hydro_get_physical_density(p, cosmo) * XH /
                       phys_const->const_proton_mass;
    const double n_H_cgs = n_H * cooling->number_density_to_cgs;

    /* compute ratefact = n_H * n_H / rho; Might lead to round-off error:
     * replaced by equivalent expression below */
    const double ratefact_cgs = n_H_cgs * XH * cooling->inv_proton_mass_cgs;

    logu_old = logu_init;
    logu = logu_old;
    int i = 0;

    float LambdaNet_old = 0;
    LambdaNet = 0;
    do /* iterate to convergence */
    {
      logu_old = logu;
      LambdaNet_old = LambdaNet;
      LambdaNet = (He_reion_heat / (dt * ratefact_cgs)) +
                  eagle_cooling_rate(logu_old, cosmo->z, n_H_cgs,
                                     abundance_ratio, n_H_index, d_n_H,
                                     He_index, d_He, cooling, &dLambdaNet_du);

      /* Newton iteration. For details on how the cooling equation is integrated
       * see documentation in theory/Cooling/ */
      logu = logu_old - (1.0 - u_ini * exp(-logu_old) -
                         LambdaNet * ratefact_cgs * dt * exp(-logu_old)) /
                            (1.0 - dLambdaNet_du * ratefact_cgs * dt);
      /* Check if first step passes over equilibrium solution, if it does adjust
       * next guess */
      if (i == 1 && LambdaNet_old * LambdaNet < 0)
        logu = newton_log_u_guess_cgs;

      /* check whether iterations go within about 10% of the table bounds,
       * if they do default to bisection method */
      if (logu > log_table_bound_high) {
        i = newton_max_iterations;
        break;
      } else if (logu < log_table_bound_low) {
        i = newton_max_iterations;
        break;
      }

      i++;
    } while (fabs(logu - logu_old) > newton_tolerance &&
             i < newton_max_iterations);
    if (i >= newton_max_iterations) {
      /* flag to trigger bisection scheme */
      *bisection_flag = 1;
    }

    return logu;
  }
