/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2019 Folkert Nobels    (nobels@strw.leidenuniv.nl)
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

/* Config parameters. */
#include "../config.h"

#include <fenv.h>

/* Local headers. */
#include "swift.h"

double pearsonfunc(double mean1, double mean2, double total12, double var1, double var2, int counter) {
  
  const double mean12 = total12 / (double)counter;
  const double correlation = (mean12 - mean1 * mean2)/ pow(var1 * var2, .5f);
  return fabs(correlation); 
}

int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* Time-step size */
  const int time_bin = 29;

  /* Try a few different values for the ID */
  for (int i = 0; i < 20; ++i) {

    const long long id = rand() * (1LL << 31) + rand();
    const integertime_t increment = (1LL << time_bin);
    const long long idoffset = id + 2;

    message("Testing id=%lld time_bin=%d", id, time_bin);

    double total = 0., total2 = 0.;
    int count = 0;

    /* Pearson correlation variables for different times */
    double sum_previous_current = 0.;
    double previous = 0.;

    /* Pearson correlation for two different IDs */
    double pearsonIDs = 0.;
    double totalID = 0.;
    double total2ID = 0.;

    /* Pearson correlation for different processes */
    double pearson_star_sf = 0.;
    double pearson_star_se = 0.;
    double pearson_star_bh = 0.;
    double pearson_sf_se = 0.;
    double pearson_sf_bh = 0.;
    double pearson_se_bh = 0.;

    /* Calculate the mean and <x^2> for these processes */
    double total_sf = 0.;
    double total_se = 0.;
    double total_bh = 0.;

    double total2_sf = 0.;
    double total2_se = 0.;
    double total2_bh = 0.;

    /* Check that the numbers are uniform over the full-range of useful
     * time-steps */
    for (integertime_t ti_current = 0LL; ti_current < max_nr_timesteps;
         ti_current += increment) {

      ti_current += increment;

      const double r =
          random_unit_interval(id, ti_current, random_number_star_formation);

      total += r;
      total2 += r * r;
      count++;

      /* Calculate for correlation between time.
       * For this we use the pearson correlation of time i and i-1 */
      sum_previous_current += r * previous;
      previous = r;

      /* Calculate if there is a correlation between different ids */
      const double r_2ndid = random_unit_interval(idoffset, ti_current,
                                                  random_number_star_formation);

      /* Pearson correlation for small different IDs */
      pearsonIDs += r * r_2ndid;
      totalID += r_2ndid;
      total2ID += r_2ndid * r_2ndid;

      /* Calculate random numbers for the different processes and check
       * that they are uncorrelated */
      
      const double r_sf = 
          random_unit_interval(id, ti_current, random_number_stellar_feedback);

      const double r_se = 
          random_unit_interval(id, ti_current, random_number_stellar_enrichment);

      const double r_bh = 
          random_unit_interval(id, ti_current, random_number_BH_feedback);

      /* Calculate the correlation between the different processes */
      total_sf += r_sf;
      total_se += r_se;
      total_bh += r_bh;

      total2_sf += r_sf * r_sf;
      total2_se += r_se * r_se;
      total2_bh += r_bh * r_bh;

      pearson_star_sf += r * r_sf;
      pearson_star_se += r * r_se;
      pearson_star_bh += r * r_bh;
      pearson_sf_se += r_sf * r_se;
      pearson_sf_bh += r_sf * r_bh;
      pearson_se_bh += r_se * r_bh;
    }

    const double mean = total / (double)count;
    const double var = total2 / (double)count - mean * mean;

    /* Pearson correlation calculation for different times */
    //const double mean_xy = sum_previous_current / ((double)count - 1.f);
    //const double correlation = (mean_xy - mean * mean) / var;
    const double correlation = pearsonfunc(mean,mean, sum_previous_current, var, var, count-1);

    /* Mean for different IDs */
    const double meanID = totalID / (double)count;
    const double varID = total2ID / (double)count - meanID * meanID;

    /* Pearson correlation between different IDs*/
    const double correlationID = pearsonfunc(mean, meanID, pearsonIDs, var, varID, count);

    /* Mean and <x^2> for different processes */
    const double mean_sf = total_sf / (double)count;
    const double mean_se = total_se / (double)count;
    const double mean_bh = total_bh / (double)count;
    
    const double var_sf = total2_sf / (double)count - mean_sf * mean_sf;
    const double var_se = total2_se / (double)count - mean_se * mean_se;
    const double var_bh = total2_bh / (double)count - mean_bh * mean_bh;

    /* Correlation between different processes */
    const double corr_star_sf = pearsonfunc(mean,mean_sf,pearson_star_sf, var, var_sf, count);
    const double corr_star_se = pearsonfunc(mean,mean_se,pearson_star_se, var, var_se, count);
    const double corr_star_bh = pearsonfunc(mean,mean_bh,pearson_star_bh, var, var_bh, count);
    const double corr_sf_se = pearsonfunc(mean_sf,mean_se,pearson_sf_se, var_sf, var_se, count);
    const double corr_sf_bh = pearsonfunc(mean_sf,mean_bh,pearson_sf_bh, var_sf, var_bh, count);
    const double corr_se_bh = pearsonfunc(mean_se,mean_bh,pearson_se_bh, var_se, var_bh, count);
    message("%e %e %e %e %e %e",corr_star_sf, corr_star_se, corr_star_bh, corr_sf_se, corr_sf_bh, corr_se_bh);
    
    /* Verify that the mean and variance match the expected values for a uniform
     * distribution */
    const double tolmean = 2e-4;
    const double tolvar = 1e-3;
    const double tolcorr = 3e-4;

    if ((fabs(mean - 0.5) / 0.5 > tolmean) ||
        (fabs(var - 1. / 12.) / (1. / 12.) > tolvar) ||
        (correlation > tolcorr) || (correlationID > tolcorr) ||
        (fabs(meanID - 0.5) / 0.5 > tolmean) ||
        (fabs(varID - 1. / 12.) / (1. / 12.) > tolvar) || 
        (corr_star_sf > tolcorr) || (corr_star_se > tolcorr) ||
        (corr_star_bh > tolcorr) || (corr_sf_se > tolcorr) ||
        (corr_sf_bh > tolcorr) || (corr_se_bh > tolcorr) || 
        (fabs(mean_sf - 0.5) / 0.5 > tolmean) ||
        (fabs(mean_se - 0.5) / 0.5 > tolmean) ||
        (fabs(mean_bh - 0.5) / 0.5 > tolmean) ||
        (fabs(var_sf - 1. / 12.) / (1. / 12.) > tolvar) || 
        (fabs(var_se - 1. / 12.) / (1. / 12.) > tolvar) || 
        (fabs(var_bh - 1. / 12.) / (1. / 12.) > tolvar)) {
      message("Test failed!");
      message("Global result:");
      message(
          "Result:    count=%d mean=%f var=%f, correlation=%f",
          count, mean, var, correlation);
      message(
          "Expected:  count=%d mean=%f var=%f, correlation=%f",
          count, 0.5f, 1. / 12., 0.);
      message("ID part");
      message("Result:     count%d mean=%f var=%f"
          " correlation=%f", count, meanID, varID, correlationID);
      message("Expected:   count%d mean=%f var=%f"
          " correlation=%f", count, .5f, 1. / 12., 0.);
      message("Different physical processes:");
      message("Means:    stars=%f stellar feedback=%f stellar "
          " enrichement=%f black holes=%f", mean, mean_sf, mean_se,
          mean_bh);
      message("Expected: stars=%f stellar feedback=%f stellar "
          " enrichement=%f black holes=%f", .5f, .5f, .5f, .5f);
      message("Var:      stars=%f stellar feedback=%f stellar "
          " enrichement=%f black holes=%f", var, var_sf, var_se,
          var_bh);
      message("Expected: stars=%f stellar feedback=%f stellar "
          " enrichement=%f black holes=%f", 1./12., 1./12., 1/12.,
          1./12.);
      message("Correlation: stars-sf=%f stars-se=%f stars-bh=%f"
          "sf-se=%f sf-bh=%f se-bh=%f", corr_star_sf, corr_star_se,
          corr_star_bh, corr_sf_se, corr_sf_bh, corr_se_bh);
      message("Expected:    stars-sf=%f stars-se=%f stars-bh=%f"
          "sf-se=%f sf-bh=%f se-bh=%f", 0., 0., 0., 0., 0., 0.);
      return 1;
    }
  }

  return 0;
}
