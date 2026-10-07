/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
/**
 * @file src/feedback/GEAR/radiation_star.c
 * @brief Star-side state of the GEAR subgrid radiation: the ionizing photon
 * rate and budget of a #spart.
 */

/* Config parameters. */
#include <config.h>

/* Local headers */
#include "error.h"
#include "inline.h"
#include "minmax.h"
#include "radiation.h"

#ifdef GEAR_SUBGRID_RADIATION_HII
/**
 * @brief Set the #spart's ionizing photon rate, split evenly across the active
 * angular pixels.
 *
 * @param sp The star.
 * @param dot_N_ion_total The total ionizing photon rate for this star.
 * @param n_HII_pixels Number of active angular pixels (from
 * GEARFeedback:HII_angular_nside via #radiation.n_HII_pixels).
 */
__attribute__((always_inline)) INLINE void radiation_set_ionizing_photon_rate(
    struct spart *sp, double dot_N_ion_total, int n_HII_pixels) {

  sp->feedback_data.radiation.n_HII_pixels = n_HII_pixels;

  const double dot_N_ion_per_pixel = dot_N_ion_total / n_HII_pixels;
  for (int p = 0; p < n_HII_pixels; p++) {
    sp->feedback_data.radiation.dot_N_ion_pix[p] = dot_N_ion_per_pixel;
  }
}
#endif /* GEAR_SUBGRID_RADIATION_HII */

/**
 * @brief Zero a #spart's radiation output, for when no radiation table is
 * loaded (#radiation.is_active = 0).
 *
 * The zero-pointered tables must not be read, and n_HII_pixels=0 would divide
 * by zero in radiation_set_ionizing_photon_rate().
 *
 * @param sp The star to zero.
 */
__attribute__((always_inline)) INLINE void radiation_zero_spart_output(
    struct spart *sp) {
  radiation_set_star_bolometric_luminosity(sp, 0.f);
  radiation_set_star_mean_excess_photon_energy_HI(sp, 0.f);
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
    radiation_set_star_band_luminosity(sp, (enum radiation_isrf_moment)m, 0.);
  radiation_set_star_teff(sp, 0.f);
  radiation_set_ionizing_photon_rate(sp, 0.0, 1);
}

#ifdef GEAR_SUBGRID_RADIATION_HII
/**
 * @brief Open this #spart's ionizing photon budget for one HII rebuild pass:
 * photons emitted over dt_back plus any overdraft carried from the last pass.
 *
 * The rate is integrated with the trapezoid rule, because a declining SSP rate
 * makes the rectangle rule under-issue photons.
 *
 * @param sp The star.
 * @param dt_back Time elapsed since this star's last HII rebuild pass.
 */
__attribute__((always_inline)) INLINE void
radiation_open_ionizing_photon_budget(struct spart *sp, double dt_back) {

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  /* All pixels share one rate, so log once per star. */
  double issued_total = 0.;
  double rate_prev_dbg = 0., rate_now_dbg = 0., rate_used_dbg = 0.;
#endif

  for (int p = 0; p < sp->feedback_data.radiation.n_HII_pixels; p++) {
    /* A pass overdraws its pixel by up to one particle's cost. Carry the debt
       forward: forgiving it would over-issue photons as 1/dt_back. Unspent
       positive budget escaped and is dropped. */
    const double debt =
        min(sp->feedback_data.radiation.N_ion_budget_pix[p], 0.);

    const double rate_now = sp->feedback_data.radiation.dot_N_ion_pix[p];
    const double rate_prev = sp->feedback_data.radiation.dot_N_ion_pix_prev[p];
    /* rate_prev < 0: first pass, no previous sample, so use rate_now. */
    const double rate_used =
        rate_prev < 0. ? rate_now : 0.5 * (rate_prev + rate_now);
    const double issued = rate_used * dt_back;

    sp->feedback_data.radiation.N_ion_budget_pix[p] = debt + issued;
    sp->feedback_data.radiation.dot_N_ion_pix_prev[p] = rate_now;

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
    issued_total += issued;
    rate_prev_dbg = rate_prev;
    rate_now_dbg = rate_now;
    rate_used_dbg = rate_used;
#endif
  }

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message(
      "HII budget open: star %lld dt_back=%e rate_prev=%e rate_now=%e "
      "rate_used=%e issued_total=%e",
      sp->id, dt_back, rate_prev_dbg, rate_now_dbg, rate_used_dbg,
      issued_total);
#endif
}

/**
 * @brief Resync the cached previous rate to the rate now, without opening a
 * budget.
 *
 * Call it when a gas-free cell skips the budget but advances
 * HII_region_last_attempt. Otherwise the next trapezoid averages against a
 * stale, too-high rate and over-issues photons.
 *
 * @param sp The star.
 */
__attribute__((always_inline)) INLINE void
radiation_resync_ionizing_photon_rate_cache(struct spart *sp) {

  for (int p = 0; p < sp->feedback_data.radiation.n_HII_pixels; p++) {
    sp->feedback_data.radiation.dot_N_ion_pix_prev[p] =
        sp->feedback_data.radiation.dot_N_ion_pix[p];
  }
}

/**
 * @brief Consume the #spart ionizing photon budget.
 *
 * @param sp The star.
 * @param pixel The angular pixel to consume from.
 * @param Delta_N_ion The ionizing photon count to remove.
 */
__attribute__((always_inline)) INLINE void radiation_consume_ionizing_photons(
    struct spart *sp, int pixel, double Delta_N_ion) {
  sp->feedback_data.radiation.N_ion_budget_pix[pixel] -= Delta_N_ion;
  return;
}
#endif /* GEAR_SUBGRID_RADIATION_HII */
