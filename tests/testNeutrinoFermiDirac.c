/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

/* Some standard headers */
#include <fenv.h>
#include <math.h>

/* Includes. */
#include "../config.h"
#include "neutrino/Default/fermi_dirac.h"
#include "swift.h"

/* Riemann function zeta(3) and zeta(5) */
#define M_ZETA_3 1.2020569031595942853997
#define M_ZETA_5 1.0369277551433699263314

int main(int argc, char *argv[]) {
  /* Exact integrals of x^n / (exp(x) + 1) on (0, infinity) */
  double integral2 = M_ZETA_3 * 1.5;
  double integral3 = M_PI * M_PI * M_PI * M_PI * 7.0 / 120.0;
  double integral4 = M_ZETA_5 * 22.5;
  double integral5 = M_PI * M_PI * M_PI * M_PI * M_PI * M_PI * 31. / 252.0;
  /* Exact moments */
  double mu = integral3 / integral2;
  double mu2 = mu * mu;
  double sigma2 = integral4 / integral2 - mu2;
  double sigma = sqrt(sigma2);
  double sigma3 = sigma2 * sigma;
  double skewns = (integral5 / integral2 - 3 * mu * sigma2 - mu2 * mu) / sigma3;

  /* Print FermiDirac(n) for n = 1, 2, 3 for external testing */
  message("fermi_dirac(1) = %e", neutrino_seed_to_fermi_dirac(1));
  message("fermi_dirac(2) = %e", neutrino_seed_to_fermi_dirac(2));
  message("fermi_dirac(3) = %e\n", neutrino_seed_to_fermi_dirac(3));

  /* Generate Fermi-Dirac numbers and compute the sum, min, and max */
  int N = 1e7;
  double sum = 0;
  double min = FLT_MAX, max = 0;
  long long seed = 290009001901;
  for (int i = 0; i < N; i++) {
    double x = neutrino_seed_to_fermi_dirac(seed + i);
    sum += x;
    if (x < min) min = x;
    if (x > max) max = x;
  }

  /* Do a second pass for the sample variance and skewness */
  double mean = sum / N;
  double ss_sum = 0;
  double sss_sum = 0;

  /* We also construct a histogram */
  int bins = 1000;
  int *histogram1 = calloc(bins, sizeof(int));

  /* Generate the same numbers again and compute statistics and histogram */
  for (int i = 0; i < N; i++) {
    double x = neutrino_seed_to_fermi_dirac(seed + i);
    ss_sum += (x - mean) * (x - mean);
    sss_sum += (x - mean) * (x - mean) * (x - mean);

    int bin = (int)((x - min) / (max - min) * bins);
    histogram1[bin]++;
  }

  /* Sample statistics */
  double var = ss_sum / (N - 1);
  double sdev = sqrt(var);
  double mu3 = sss_sum / (N - 2) / (var * sdev);

  /* Relative errors with the exact moments */
  double err_mu = mean / mu - 1.;
  double err_sig = var / sigma2 - 1.;
  double err_skew = mu3 / skewns - 1.;

  message("Sample mean is %f, exact: %f (rel. %e).", mean, mu, err_mu);
  message("Sample variance is %f, exact: %f (rel. %e).", var, sigma2, err_sig);
  message("Sample skewness is %f, exact: %f (rel. %e).", mu3, skewns, err_skew);

  assert(fabs(err_mu) < 1e-3);
  assert(fabs(err_sig) < 1e-2);
  assert(fabs(err_skew) < 1e-2);

  /* Construct Kolmogorov-Smirnov statistic by integrating the histogram */
  double sum_empirical = 0.;
  double sum_analytical = 0.;
  double KS_statistic = 0.;
  double dx = (max - min) / bins;
  for (int bin = 0; bin < bins; bin++) {
    double x = min + (bin + 0.5) * dx;
    sum_empirical += histogram1[bin] * 1. / N;
    sum_analytical += fermi_dirac_density(x) * x * x * dx / integral2;

    double delta = fabs(sum_empirical - sum_analytical);
    if (delta > KS_statistic) {
      KS_statistic = delta;
    }
  }

  /* Can we reject the hypothesis that these numbers are FD distributed? */
  double crit_val = 5.146991e-05 * sqrt(1e9 / N);  // 99% confidence level
  message("KS statistic = %e (99%% that KS < %e)\n", KS_statistic, crit_val);

  /* We should not reject this */
  assert(KS_statistic < crit_val);

  free(histogram1);

  message("Success.");

  return 0;
}
