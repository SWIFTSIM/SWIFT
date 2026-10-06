/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
 * @file src/feedback/GEAR/radiation_gas.c
 * @brief Per-gas-particle radiation feedback physics for GEAR: hydrogen
 * content, ionized-state temperature and recombination rate, the ionizing
 * photon budget, and the ionization tag carried on each #part.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "atomic.h"
#include "chemistry.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "inline.h"
#include "minmax.h"
#include "radiation.h"
#include "units.h"

/**
 * @brief Total hydrogen mass fraction of this #part, from its composition.
 *
 * Not cooling_get_hydrogen_mass_fraction(): at COOLING_GRACKLE_MODE >= 2 that
 * omits the hydrogen bound in H2/H-.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @return Total hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static double
radiation_get_part_total_hydrogen_mass_fraction(
    const struct cooling_function_data *cooling, const struct part *p) {

  const double Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Clamped: Z > HydrogenFractionByMass would give a negative photon cost. */
  return max(cooling->HydrogenFractionByMass - Z, 0.);
}

/**
 * @brief Get the gas number of hydrogen atoms.
 *
 * @param phys_const Physical constants.
 * @param hydro_props The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Number of hydrogen atoms.
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_number_hydrogen_atoms(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  const float m = hydro_get_mass(p);
  const double m_p = phys_const->const_proton_mass;
  const double X_H =
      radiation_get_part_total_hydrogen_mass_fraction(cooling, p);

  /* Number of hydrogen atoms in b (Hu et al. 2017; Smith et al. 2021). */
  const double N_H = (X_H * m) / m_p;

  return N_H;
}

/**
 * @brief Get the gas number of NEUTRAL hydrogen atoms, from the tracked
 * species fractions.
 *
 * Prices only the one-off cost of claiming a fresh candidate: an already
 * (partly) ionized particle need not pay again. The maintenance cost uses the
 * total. At COOLING_GRACKLE_MODE == 0 it falls back to the total N_H.
 *
 * @param phys_const Physical constants.
 * @param hydro_props The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Number of neutral hydrogen atoms.
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_number_neutral_hydrogen_atoms(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

#if COOLING_GRACKLE_MODE >= 1
  const float m = hydro_get_mass(p);
  const double m_p = phys_const->const_proton_mass;
  const struct cooling_xpart_data *cool_data = &xp->cooling_data;

  double X_HI = cool_data->HI_frac;
#if COOLING_GRACKLE_MODE >= 2
  /* H2 and H- hydrogen is neutral; H2II is already ionized, so excluded. */
  X_HI += cool_data->H2I_frac + cool_data->HM_frac;
#endif
  const double N_HI = (X_HI * m) / m_p;

  /* Capped at the composition total: species can transiently exceed it. */
  const double N_H = radiation_get_part_number_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);
  return min(N_HI, N_H);
#else
  return radiation_get_part_number_hydrogen_atoms(phys_const, hydro_props, us,
                                                  cosmo, cooling, p, xp);
#endif
}

/**
 * @brief Collisional-equilibrium temperature floor as a function of Z
 * (Hopkins 2023 fit; theory/GEAR/Radiation/01_algorithm.tex, Eq. tcollisional).
 *
 * Compile with -DIONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K=<value>
 * to force a fixed value regardless of Z.
 *
 * @param Z Metal mass fraction.
 * @return Collisional-equilibrium temperature (Kelvin).
 */
__attribute__((always_inline)) INLINE double radiation_get_T_collisional_K(
    const double Z) {

#ifdef IONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K
  return IONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K;
#else
  const double Z_sun = 0.02;
  const double ten_to_four_K = 1e4;

  /* Guard against Z << Z_sun, where the fit below would give T < 0. */
  if (Z >= Z_sun * 1e-3) {
    /* Hopkins (2023)'s fit is in log10(Z/Z_sun), not ln. */
    const double tmp = 0.86 / (1 + 0.22 * log10(Z / Z_sun));
    return ten_to_four_K * min(6.62, tmp);
  } else {
    return 6.62 * ten_to_four_K; /* High-temperature asymptote */
  }
#endif
}

/**
 * @brief Specific internal energy of this #part once ionized: the minimum of
 * the energy to fully ionize it and the collisional-equilibrium energy.
 *
 * Shared by cooling_ionize_part_subgrid() and
 * radiation_get_part_rate_to_fully_ionize().
 *
 * @param phys_const Physical constants.
 * @param hydro_props The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Specific internal energy (physical, code units).
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_ionized_internal_energy(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  const double m_p = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  const double N_H = radiation_get_part_number_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);
  const double E_ion =
      2.17872e-11 / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  const double Delta_u_ionized = N_H * E_ion / hydro_get_mass(p);

  const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(p);
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);

  const double T_collisional =
      radiation_get_T_collisional_K(Z) /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  const double u_collisional =
      cooling_internal_energy_from_T(T_collisional, mu, k_B, m_p);

  return min(Delta_u_ionized, u_collisional);
}

/**
 * @brief Case-B hydrogen recombination coefficient, temperature-dependent
 * (Hui & Gnedin 1997, MNRAS 292, 27, Appendix A; their fit to Ferland et
 * al. 1992, accurate to 0.7% from 1 K to 1e9 K).
 *
 * @param T Temperature in Kelvin.
 * @return alpha_B in cm^3/s (CGS).
 */
__attribute__((always_inline)) INLINE double
radiation_get_case_b_recombination_coefficient_cgs(const double T) {
  /* Floor avoids lambda=inf and inf*0 = NaN for T=0. */
  const double lambda = 315614.0 / max(T, 1.0);
  return 2.753e-14 * pow(lambda, 1.5) *
         pow(1.0 + pow(lambda / 2.740, 0.407), -2.242);
}

/**
 * @brief Get the gas ionizing rate needed to fully ionize the #part.
 *
 * @param phys_const Physical constants.
 * @param hydro_props The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Ionizing photon rate to ionize this #part (physical units).
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_rate_to_fully_ionize(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  const float rho = hydro_get_physical_density(p, cosmo);
  const double m_p = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;
  const double X_H =
      radiation_get_part_total_hydrogen_mass_fraction(cooling, p);

  /* Number of hydrogen atoms in b */
  const double N_H = radiation_get_part_number_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);

  /* Z >= HydrogenFractionByMass gives N_H = 0: nothing to ionize, and T=0
     would give a NaN in the recombination fit. */
  if (N_H <= 0.) return 0.;

  /* Electron density assuming full ionization (n_e ~= n_H). */
  const double n_e = (X_H * rho) / m_p;

  /* Case-B coefficient at the temperature the gas is held at once ionized. */
  const double u_ionized = radiation_get_part_ionized_internal_energy(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);
  const double T_ionized_K =
      cooling_temperature_from_internal_energy(u_ionized, mu, k_B, m_p) *
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  const double beta_cgs =
      radiation_get_case_b_recombination_coefficient_cgs(T_ionized_K);
  const float dimension_alphaB[5] = {0, 3, -1, 0, 0}; /* [cm^3 s^-1] */
  const double beta =
      beta_cgs / units_general_cgs_conversion_factor(us, dimension_alphaB);

  /* Required ionizing rate in [photons / internal time unit] */
  const double Delta_N_dot = N_H * beta * n_e;

  return Delta_N_dot;
}

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
  sp->feedback_data.radiation.L_bol = 0.f;
  sp->feedback_data.radiation.mean_excess_photon_energy_HI = 0.f;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
    sp->feedback_data.radiation.L_band[m] = 0.;
  sp->feedback_data.radiation.teff = 0.f;
  radiation_set_ionizing_photon_rate(sp, 0.0, 1);
}

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

/**
 * @brief Tag the #part as ionized to be ionized in feedback_update_part().
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param star_id The id of the star that ionized this particle.
 * @param end_time The simulation time until which the particle stays tagged
 * (the star's next HII rebuild).
 * @param excess_photon_energy_HI Mean photon energy above the 13.6 eV HI
 * threshold, erg; 0 unless GEARFeedback:HII_couple_ionization_rate is on.
 * @param photoionization_rate_HI Gamma_HI at the particle, internal 1/time; 0
 * unless GEARFeedback:HII_couple_ionization_rate is on.
 */
__attribute__((always_inline)) INLINE void radiation_tag_part_as_ionized(
    struct part *p, struct xpart *xp, long long star_id, double end_time,
    float excess_photon_energy_HI, float photoionization_rate_HI) {
  p->feedback_data.is_ionized = 1;
  p->feedback_data.star_id = star_id;
  p->feedback_data.end_time = end_time;
  xp->feedback_data.HII_region.excess_photon_energy_HI =
      excess_photon_energy_HI;
  xp->feedback_data.HII_region.photoionization_rate_HI =
      photoionization_rate_HI;
  return;
}

/**
 * @brief Reset the #part ionization tag.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE void radiation_reset_part_ionized_tag(
    struct part *p, struct xpart *xp) {
  p->feedback_data.is_ionized = 0;
  return;
}

/**
 * @brief Is this #part *tagged* as ionized ?
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Is the particle *tagged* ionized?
 */
__attribute__((always_inline)) INLINE char radiation_is_part_tagged_as_ionized(
    const struct part *p, const struct xpart *xp) {
  return p->feedback_data.is_ionized;
}

/**
 * @brief The simulation time at which the ionization tag expires. Valid only
 * while tagged ionized.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return The simulation time at which the ionization tag expires.
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_ionized_end_time(const struct part *p,
                                    const struct xpart *xp) {
  return p->feedback_data.end_time;
}

/**
 * @brief Clear #part::feedback_data.is_illuminated_ISRF once its illumination
 * window has lapsed. Called once per step per particle (feedback_reset_part).
 *
 * With LW/PE propagation off, the band u values are also zeroed: nothing else
 * decays them. With propagation on, its own per-step update does.
 *
 * @param p The particle.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE void
radiation_reset_part_ISRF_illumination_tag(struct part *p,
                                           const struct engine *e) {
  if (!p->feedback_data.is_illuminated_ISRF) return;
  if (e->ti_current < p->feedback_data.ISRF_illumination_end_ti) return;

  p->feedback_data.is_illuminated_ISRF = 0;

  if (!e->feedback_props->ISRF_propagation) {
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
      p->feedback_data.isrf_moment[m].u = 0.f;
  }
}

/**
 * @brief Id of the star that ionized this #part. Valid only while tagged
 * ionized.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return The id of the ionizing star.
 */
__attribute__((always_inline)) INLINE long long
radiation_get_part_ionized_star_id(const struct part *p,
                                   const struct xpart *xp) {
  return p->feedback_data.star_id;
}

/**
 * @brief Mean photon energy above the 13.6 eV HI threshold of the tagging star,
 * in cgs (erg), frozen at tag time. Valid only while tagged ionized.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Mean excess photon energy above the HI threshold, cgs erg.
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_excess_photon_energy_HI(const struct part *p,
                                           const struct xpart *xp) {
  return xp->feedback_data.HII_region.excess_photon_energy_HI;
}

/**
 * @brief Photoionization rate coefficient Gamma_HI frozen at tag time (internal
 * 1/time). Valid only while tagged ionized.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Photoionization rate coefficient Gamma_HI, internal 1/time.
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_photoionization_rate_coefficient(const struct part *p,
                                                    const struct xpart *xp) {
  return xp->feedback_data.HII_region.photoionization_rate_HI;
}

/**
 * @brief Photoionization rate coefficient Gamma_HI from an HI-ionizing photon
 * flux (internal units), with sigma_HI = 6.3e-18 cm^2 at the Lyman limit
 * (Osterbrock & Ferland 2006). Computed at tag time: the raw flux overflows
 * float32 but the product with the cross section does not.
 *
 * @param us Unit system.
 * @param ionizing_flux_HI HI-ionizing photon flux (internal units).
 * @return Gamma_HI (internal units).
 */
__attribute__((always_inline)) INLINE double
radiation_get_photoionization_rate_coefficient_from_flux_HI(
    const struct unit_system *us, const double ionizing_flux_HI) {

  const double sigma_HI_cgs = 6.3e-18; /* [cm^2], Osterbrock & Ferland 2006 */
  const float dimension_area[5] = {0, 2, 0, 0, 0}; /* [cm^2] */
  const double sigma_HI =
      sigma_HI_cgs / units_general_cgs_conversion_factor(us, dimension_area);

  return sigma_HI * ionizing_flux_HI;
}

/* Warning throttle: the first N clamp events are reported, then one summary
   per interval. */
#define RADIATION_NEGATIVE_CLAMP_WARN_LIMIT 10
#define RADIATION_NEGATIVE_CLAMP_SUMMARY_INTERVAL 10000

/**
 * @brief Clamp a quantity about to reach Grackle to be non-negative, with a
 * throttled warning.
 *
 * The propagated PE/LW energy can undershoot below zero, and a negative flux
 * would act as spurious cooling or dissociation in Grackle. This masks the
 * undershoot at the interface, it does not fix it.
 *
 * @param name Human-readable name of the quantity, for the warning message.
 * @param value The value about to be sent to Grackle.
 * @param count Running clamp-event count for this quantity (updated).
 * @param worst Most negative value seen for this quantity so far (updated).
 * @return value, or 0 if value was negative.
 */
static double radiation_clamp_nonnegative_for_grackle(const char *name,
                                                      double value,
                                                      volatile long long *count,
                                                      volatile double *worst) {

  if (value >= 0.) return value;

  atomic_min_d(worst, value);
  const long long n = atomic_add(count, 1LL) + 1LL;

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message("Clamped negative %s = %g to 0 before passing it to Grackle.", name,
          value);
#endif

  if (n <= RADIATION_NEGATIVE_CLAMP_WARN_LIMIT ||
      n % RADIATION_NEGATIVE_CLAMP_SUMMARY_INTERVAL == 0) {
    warning(
        "Clamped %lld negative %s value(s) reaching Grackle to zero so far "
        "this run (this occurrence: %g, worst seen: %g). A negative flux is "
        "unphysical; it indicates an undershoot in the LW/PE propagation "
        "scheme that this clamp only masks at the Grackle interface.",
        n, name, value, *worst);
  }

  return 0.;
}

/*! Moment names for the clamp warnings, indexed by #radiation_isrf_moment. */
static const char *const radiation_isrf_moment_clamp_name[] = {
    "PE-band specific energy", "LW-band specific energy",
    "LW-band photon-number specific energy"};
_Static_assert(sizeof(radiation_isrf_moment_clamp_name) /
                       sizeof(radiation_isrf_moment_clamp_name[0]) ==
                   ISRF_MOMENT_COUNT,
               "radiation_isrf_moment_clamp_name needs one initialiser per "
               "ISRF_MOMENT_COUNT entry.");

/**
 * @brief Fetch one ISRF band's specific energy, clamped to be non-negative.
 *
 * Each band is clamped before any sum, so a negative LW energy cannot cancel
 * a positive PE energy unnoticed. The count is per band read, not per
 * particle. The stored value is unchanged.
 *
 * @param p The particle.
 * @param m The moment to read (#radiation_isrf_moment).
 * @return The band's specific energy, or 0 if it was negative.
 */
static double radiation_get_band_u_nonnegative(const struct part *p,
                                               const int m) {

  static volatile long long band_clamp_count[ISRF_MOMENT_COUNT] = {0};
  static volatile double band_clamp_worst[ISRF_MOMENT_COUNT] = {0.};

  return radiation_clamp_nonnegative_for_grackle(
      radiation_isrf_moment_clamp_name[m],
      (double)p->feedback_data.isrf_moment[m].u, &band_clamp_count[m],
      &band_clamp_worst[m]);
}

/**
 * @brief Local ISRF strength in Habing units: G0 = c*rho*u /
 * #RADIATION_HABING_FLUX_CGS, with u the PE+LW specific energy. Feeds
 * Grackle's isrf_habing array. Zero for a particle never illuminated.
 *
 * Each band is clamped to be non-negative before the sum, so an undershot band
 * cannot cancel the other (#radiation_get_band_u_nonnegative).
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param p The particle.
 * @return G0, dimensionless (Habing units), never negative.
 */
double radiation_get_part_isrf_habing(const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct cosmology *cosmo,
                                      const struct part *p) {

  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_sum = radiation_get_band_u_nonnegative(p, ISRF_MOMENT_PE) +
                       radiation_get_band_u_nonnegative(p, ISRF_MOMENT_LW);
  const double flux = phys_const->const_speed_light_c * rho * u_sum;
  const double flux_cgs =
      flux *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_FLUX_PER_UNIT_SURFACE);

  static volatile long long isrf_habing_clamp_count = 0;
  static volatile double isrf_habing_clamp_worst = 0.;

  return radiation_clamp_nonnegative_for_grackle(
      "ISRF Habing flux", flux_cgs / RADIATION_HABING_FLUX_CGS,
      &isrf_habing_clamp_count, &isrf_habing_clamp_worst);
}

/**
 * @brief H2 Lyman-Werner photodissociation rate from this #part's LW-band
 * specific energy: k_diss = (sigma_H2 / E_LW) * energy flux.
 *
 * Only the quotient is constrained (Sternberg et al. 2014), so it is read from
 * #RADIATION_SIGMA_H2_OVER_E_LW_CGS. #radiation_lw_photon_energy_cgs is a
 * diagnostic and has no effect here. Feeds Grackle's RT_H2_dissociation_rate
 * (COOLING_GRACKLE_MODE > 1 only). The band and the rate are both clamped to
 * be non-negative.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param p The particle.
 * @return H2 photodissociation rate, internal 1/time, never negative.
 */
double radiation_get_part_LW_dissociation_rate_internal(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct part *p) {

  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_LW = radiation_get_band_u_nonnegative(p, ISRF_MOMENT_LW);
  const double flux_LW = phys_const->const_speed_light_c * rho * u_LW;
  const double flux_LW_cgs =
      flux_LW *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_FLUX_PER_UNIT_SURFACE);

  const double k_diss_cgs = RADIATION_SIGMA_H2_OVER_E_LW_CGS * flux_LW_cgs;
  const double k_diss =
      k_diss_cgs / units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);

  static volatile long long lw_dissociation_clamp_count = 0;
  static volatile double lw_dissociation_clamp_worst = 0.;

  return radiation_clamp_nonnegative_for_grackle(
      "LW photodissociation rate", k_diss, &lw_dissociation_clamp_count,
      &lw_dissociation_clamp_worst);
}
