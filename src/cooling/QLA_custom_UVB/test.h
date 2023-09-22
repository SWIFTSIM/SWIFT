/*! \brief Compute gas temperature from internal energy per unit mass.
 *
 *   This function determines the electron fraction, and hence the mean
 *   molecular weight. With it arrives at a self-consistent temperature.
 *   Element abundances and the rates for the emission are also computed.
 *
 *  \param[in] u   internal energy per unit mass.
 *  \param[in] rho gas density.
 *  \param[in, out] ne_guess electron number density relative to hydrogen
 *                  number density
 *
 *  \return The gas temperature.
 */
double convert_u_to_temp(double u, double rho, double *ne_guess,
                         const struct cooling_function_data *cooling,
                         GasState *gs) {
  double temp, temp_old, temp_new, max = 0, ne_old;
  double mu;
  int iter = 0;

  double u_input, rho_input, ne_input;

  u_input = u;
  rho_input = rho;
  ne_input = *ne_guess;

  mu = (1 + 4 * gs->yhelium) / (1 + gs->yhelium + *ne_guess);
  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  do {
    ne_old = *ne_guess;

    find_abundances_and_rates(log10(temp), rho, ne_guess, cooling, &gs);
    temp_old = temp;

    mu = (1 + 4 * gs->yhelium) / (1 + gs->yhelium + *ne_guess);

    temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

    max =
        dmax(max, temp_new / (1 + gs->yhelium + *ne_guess) *
                      fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

    temp = temp_old + (temp_new - temp_old) / (1 + max);
    iter++;

    if (iter > (MAXITER - 10)) printf("-> temp= %g ne=%g\n", temp, *ne_guess);
  } while (fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

  if (iter >= MAXITER) {
    printf("failed to converge in convert_u_to_temp()\n");
    printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input,
           ne_input);
    printf("convergence failure");
    exit(1);
  }

  // gs->mu = mu;

  return temp;
}

/*! \brief Computes the actual abundance ratios.
 *
 *  The chemical composition of the gas is primordial (no metals are present).
 *
 *  \param[in] logT log10 of gas temperature.
 *  \param[in] rho Gas density.
 *  \param[in, out] ne_guess Electron number density relative to hydrogen
 *                  number density.
 *
 *  \return void
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess,
                               const struct cooling_function_data *cooling,
                               GasState *gs) {
  double neold, nenew;
  int j, niter;
  double flow, fhi, t;

  double logT_input, rho_input, ne_input;

  logT_input = logT;
  rho_input = rho;
  ne_input = *ne_guess;

  // if(!gsl_finite(logT))
  //   terminate("logT=%g\n", logT);

  if (logT <= Tmin) /* everything neutral */
  {
    gs->nH0 = 1.0;
    gs->nHe0 = gs->yhelium;
    gs->nHp = 0;
    gs->nHep = 0;
    gs->nHepp = 0;
    gs->ne = 0;
    *ne_guess = 0;
    return;
  }

  if (logT >= Tmax) /* everything is ionized */
  {
    gs->nH0 = 0;
    gs->nHe0 = 0;
    gs->nHp = 1.0;
    gs->nHep = 0;
    gs->nHepp = gs->yhelium;
    gs->ne = gs->nHp + 2.0 * gs->nHepp;
    *ne_guess = gs->ne; /* note: in units of the hydrogen number density */
    return;
  }

  t = (logT - Tmin) / deltaT;
  j = (int)t;
  fhi = t - j;
  flow = 1 - fhi;

  if (*ne_guess == 0) *ne_guess = 1.0;

  gs->nHcgs = gs->XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */

  gs->ne = *ne_guess;
  neold = gs->ne;
  niter = 0;
  gs->necgs = gs->ne * gs->nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do {
    niter++;

    gs->aHp =
        flow * cooling->RateT[j].AlphaHp + fhi * cooling->RateT[j + 1].AlphaHp;
    gs->aHep = flow * cooling->RateT[j].AlphaHep +
               fhi * cooling->RateT[j + 1].AlphaHep;
    gs->aHepp = flow * cooling->RateT[j].AlphaHepp +
                fhi * cooling->RateT[j + 1].AlphaHepp;
    gs->ad =
        flow * cooling->RateT[j].Alphad + fhi * cooling->RateT[j + 1].Alphad;
    gs->geH0 = flow * cooling->RateT[j].GammaeH0 +
               fhi * cooling->RateT[j + 1].GammaeH0;
    gs->geHe0 = flow * cooling->RateT[j].GammaeHe0 +
                fhi * cooling->RateT[j + 1].GammaeHe0;
    gs->geHep = flow * cooling->RateT[j].GammaeHep +
                fhi * cooling->RateT[j + 1].GammaeHep;

    if (gs->necgs <= 1.e-25 || cooling->pc.J_UV == 0) {
      gs->gJH0ne = gs->gJHe0ne = gs->gJHepne = 0;
    } else {
      gs->gJH0ne = cooling->pc.gJH0 / gs->necgs;
      gs->gJHe0ne = cooling->pc.gJHe0 / gs->necgs;
      gs->gJHepne = cooling->pc.gJHep / gs->necgs;
    }

    gs->nH0 = gs->aHp / (gs->aHp + gs->geH0 + gs->gJH0ne); /* eqn (33) */
    gs->nHp = 1.0 - gs->nH0;                               /* eqn (34) */

    if ((gs->gJHe0ne + gs->geHe0) <= SMALLNUM) /* no ionization at all */
    {
      gs->nHep = 0.0;
      gs->nHepp = 0.0;
      gs->nHe0 = gs->yhelium;
    } else {
      gs->nHep =
          gs->yhelium / (1.0 + (gs->aHep + gs->ad) / (gs->geHe0 + gs->gJHe0ne) +
                         (gs->geHep + gs->gJHepne) / gs->aHepp); /* eqn (35) */
      gs->nHe0 = gs->nHep * (gs->aHep + gs->ad) /
                 (gs->geHe0 + gs->gJHe0ne); /* eqn (36) */
      gs->nHepp =
          gs->nHep * (gs->geHep + gs->gJHepne) / gs->aHepp; /* eqn (37) */
    }

    neold = gs->ne;

    gs->ne = gs->nHp + gs->nHep + 2 * gs->nHepp; /* eqn (38) */
    gs->necgs = gs->ne * gs->nHcgs;

    if (cooling->pc.J_UV == 0) break;

    nenew = 0.5 * (gs->ne + neold);
    gs->ne = nenew;
    gs->necgs = gs->ne * gs->nHcgs;

    if (fabs(gs->ne - neold) < 1.0e-4) break;

    if (niter > (MAXITER - 10)) printf("ne= %g  niter=%d\n", gs->ne, niter);
  } while (niter < MAXITER);

  if (niter >= MAXITER) {
    printf("gs->aHp = %le\n", gs->aHp);
    char buff[1000];
    // sprintf(buff, "%s/cooling_task%d.dat", All.OutputDir, ThisTask);
    FILE *fp = fopen(buff, "w");
    // fwrite(&All.Time, sizeof(double), 1, fp);
    fwrite(&logT_input, sizeof(double), 1, fp);
    fwrite(&rho_input, sizeof(double), 1, fp);
    fwrite(&ne_input, sizeof(double), 1, fp);
    fclose(fp);
    printf("no convergence reached in find_abundances_and_rates");
    exit(1);
  }
  gs->bH0 =
      flow * cooling->RateT[j].BetaH0 + fhi * cooling->RateT[j + 1].BetaH0;
  gs->bHep =
      flow * cooling->RateT[j].BetaHep + fhi * cooling->RateT[j + 1].BetaHep;
  gs->bff =
      flow * cooling->RateT[j].Betaff + fhi * cooling->RateT[j + 1].Betaff;

  *ne_guess = gs->ne;
}

/*! \brief  Calculate (heating rate-cooling rate)/n_h^2 in cgs units.
 *
 *  \param[in] logT log10 of gas temperature.
 *  \param[in] rho Gas density.
 *  \param[in, out] nelec Electron number density relative to hydrogen number
 *                  density.
 *
 *  \return (heating rate-cooling rate)/n_h^2.
 */
double CoolingRate(double logT, double rho, double redz, double *nelec,
                   const struct cooling_function_data *cooling, GasState *gs) {
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;
  double LambdaPrim = 0.0, LambdaMet = 0.0, LambdaDust = 0.0, LambdaMol = 0.0;

  if (logT <= Tmin) logT = Tmin + 0.5 * deltaT; /* floor at Tmin */

  gs->nHcgs = gs->XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */

  if (logT < Tmax) {
    find_abundances_and_rates(logT, rho, nelec, cooling, &gs);

    /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
    T = pow(10.0, logT);

    LambdaExcH0 = gs->bH0 * gs->ne * gs->nH0;
    LambdaExcHep = gs->bHep * gs->ne * gs->nHep;
    LambdaExc = LambdaExcH0 + LambdaExcHep; /* excitation */
    LambdaIonH0 = 2.18e-11 * gs->geH0 * gs->ne * gs->nH0;
    LambdaIonHe0 = 3.94e-11 * gs->geHe0 * gs->ne * gs->nHe0;
    LambdaIonHep = 8.72e-11 * gs->geHep * gs->ne * gs->nHep;
    LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep; /* ionization */
    LambdaRecHp = 1.036e-16 * T * gs->ne * (gs->aHp * gs->nHp);
    LambdaRecHep = 1.036e-16 * T * gs->ne * (gs->aHep * gs->nHep);
    LambdaRecHepp = 1.036e-16 * T * gs->ne * (gs->aHepp * gs->nHepp);
    LambdaRecHepd = 6.526e-11 * gs->ad * gs->ne * gs->nHep;
    LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;
    LambdaFF = gs->bff * (gs->nHp + gs->nHep + 4 * gs->nHepp) * gs->ne;
    LambdaPrim = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

    redshift = redz;
    LambdaCmptn = 5.65e-36 * gs->ne * (T - 2.73 * (1. + redshift)) *
                  pow(1. + redshift, 4.) / gs->nHcgs;
    // LambdaCmptn = 0.0;
    Lambda = LambdaPrim + LambdaMet + LambdaDust + LambdaCmptn + LambdaMol;

    Heat = 0;
    if (cooling->pc.J_UV != 0)
      Heat += (gs->nH0 * cooling->pc.epsH0 + gs->nHe0 * cooling->pc.epsHe0 +
               gs->nHep * cooling->pc.epsHep) /
              gs->nHcgs;
  } else /* here we're outside of tabulated rates, T>Tmax K */
  {
    /* at high T (fully ionized); only free-free and Compton cooling are
     * present. Assumes no heating. */
    Heat = 0;

    LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
        LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

    /* very hot: H and He both fully ionized */
    gs->nHp = 1.0;
    gs->nHep = 0;
    gs->nHepp = gs->yhelium;
    gs->ne = gs->nHp + 2.0 * gs->nHepp;
    *nelec = gs->ne; /* note: in units of the hydrogen number density */

    T = pow(10.0, logT);
    LambdaFF = 1.42e-27 * sqrt(T) *
               (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) *
               (gs->nHp + 4 * gs->nHepp) * gs->ne;

    redshift = redz;
    /* add inverse Compton cooling off the microwave background */
    LambdaCmptn = 5.65e-36 * gs->ne * (T - 2.73 * (1. + redshift)) *
                  pow(1. + redshift, 4.) / gs->nHcgs;
    // LambdaCmptn = 0.0;
    Lambda = LambdaFF + LambdaCmptn;
  }

  return (Heat - Lambda);
}
