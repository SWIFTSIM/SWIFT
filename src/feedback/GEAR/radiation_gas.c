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

/* Include header */
#include "radiation.h"

#include "chemistry.h"
#include "cooling.h"
#include "minmax.h"
#include "units.h"

/**
 * Total hydrogen mass fraction of this #part, from its composition alone.
 *
 * Not cooling_get_hydrogen_mass_fraction(): that returns HI_frac + HII_frac at
 * COOLING_GRACKLE_MODE >= 1, i.e. the *atomic* hydrogen only, so at mode >= 2
 * it omits the hydrogen bound in H2/H- and under-reports what this model has to
 * ionize and then keep ionized. Composition is also the one source available in
 * every cooling mode.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @return Total hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static double
radiation_get_part_total_hydrogen_mass_fraction(
    const struct cooling_function_data *cooling, const struct part *p) {

  const double Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Clamped: a pathologically enriched particle (Z > HydrogenFractionByMass)
     would otherwise give a negative N_H, hence a negative photon cost, and
     radiation_consume_ionizing_photons() would manufacture photons rather
     than spend them. The star-level formula below guards this the same way. */
  return max(cooling->HydrogenFractionByMass - Z, 0.);
}

/**
 * Get the gas number of hydrogen atoms.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
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
 * Get the gas number of NEUTRAL hydrogen atoms, from the tracked species
 * fractions rather than total composition. Used only to price the one-off
 * cost of claiming a fresh candidate (feedback_iact_HII_ionization): a
 * particle whose species are already partly or fully ionized -- re-tagged
 * after its previous tag lapsed on a marginal budget shortfall, or
 * pre-ionized by a UV background -- does not need to pay to strip
 * electrons it has already lost. The maintenance cost
 * (radiation_get_part_rate_to_fully_ionize) keeps using the total, since
 * upkeep of a fully-ionized particle scales with its full electron/proton
 * content, not with what was neutral before this pass.
 *
 * At COOLING_GRACKLE_MODE == 0 (no species tracking) there is nothing to
 * distinguish neutral from ionized, so this falls back to the total N_H
 * (radiation_get_part_number_hydrogen_atoms) -- the pre-existing,
 * conservative behaviour.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
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
  /* Hydrogen locked in H2/H- is also neutral (not yet stripped); H2II is
     already singly-ionized, so it is excluded -- the same H2II
     approximation radiation_get_part_total_hydrogen_mass_fraction's own
     doxygen already notes for the total-hydrogen accounting. */
  X_HI += cool_data->H2I_frac + cool_data->HM_frac;
#endif
  const double N_HI = (X_HI * m) / m_p;

  /* Capped at the composition total: species can transiently drift above
     it (e.g. mid-way through a Grackle sub-step), and a neutral count
     above the particle's total hydrogen content would make this exceed
     radiation_get_part_number_hydrogen_atoms itself, defeating the point
     of pricing on neutral content specifically. */
  const double N_H = radiation_get_part_number_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);
  return min(N_HI, N_H);
#else
  return radiation_get_part_number_hydrogen_atoms(phys_const, hydro_props, us,
                                                  cosmo, cooling, p, xp);
#endif
}

/**
 * Metallicity-dependent collisional-equilibrium temperature floor
 * (Hopkins 2023's photoionization temperature fit; see
 * theory/GEAR/Radiation/01_algorithm.tex, Eq. tcollisional). Depends
 * only on Z (used by radiation_get_part_ionized_internal_energy).
 *
 * Compile with -DIONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K=<value>
 * to force this to a fixed value regardless of Z -- e.g. to reproduce a paper's
 * own flat T_i=1e4 K convention at Z=0 (pure hydrogen, no metal-line
 * cooling), decoupling the ionized-gas temperature from the metallicity
 * this fit would otherwise require to hit that value (Z/Zsun~0.231 for
 * 1e4 K, which brings real metal cooling along with it).
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
 * Get the specific internal energy this #part would be held at once
 * ionized: the minimum of the energy needed to fully ionize it and the
 * metallicity-dependent collisional-equilibrium energy (see
 * cooling_ionize_part_subgrid in cooling_gear_subgrid.h). Pure
 * computation, no side effects -- shared by cooling_ionize_part_subgrid
 * (which actually floors the particle's temperature) and
 * radiation_get_part_rate_to_fully_ionize (which evaluates the case-B
 * recombination coefficient at the temperature the gas is actually held
 * at, instead of a fixed 1e4 K), so the two stay consistent with each
 * other.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
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
 * Case-B hydrogen recombination coefficient, temperature-dependent
 * (Hui & Gnedin 1997, MNRAS 292, 27, Appendix A; their fit to Ferland et
 * al. 1992, accurate to 0.7% from 1 K to 1e9 K).
 *
 * @param T Temperature in Kelvin.
 * @return alpha_B in cm^3/s (CGS).
 */
__attribute__((always_inline)) INLINE double
radiation_get_case_b_recombination_coefficient_cgs(const double T) {
  /* Floored: T=0 (e.g. a fully-metal particle's N_H==0 chain, or any other
     caller's edge case) would otherwise give lambda=inf and the fit below
     inf*0 = NaN. 1 K is far outside the fit's accurate range (1-1e9 K per
     Hui & Gnedin) and only ever reached through a guard like this one. */
  const double lambda = 315614.0 / max(T, 1.0);
  return 2.753e-14 * pow(lambda, 1.5) *
         pow(1.0 + pow(lambda / 2.740, 0.407), -2.242);
}

/**
 * Get the gas ionizing rate needed to fully ionize the #part.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
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

  /* Z >= HydrogenFractionByMass clamps N_H to 0 (this file's X_H clamp);
     nothing to ionize, and the temperature machinery below would give
     u_ionized=0 -> T=0 -> NaN in the recombination fit for zero cost. */
  if (N_H <= 0.) return 0.;

  /* Electron density assuming full ionization (n_e ~= n_H). */
  const double n_e = (X_H * rho) / m_p;

  /* Case-B recombination coefficient, evaluated at the temperature this
     particle is actually held at once ionized (Hui & Gnedin 1997),
     instead of the fixed 1e4 K convention. */
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
 * Set the #spart's ionizing photon rate, split evenly across the active
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
 * Open this #spart's ionizing photon budget for one HII rebuild pass:
 * convert its emission rate into the photon count emitted over dt_back, the
 * time elapsed since the previous pass, plus any overdraft carried from it.
 *
 * @param sp The star.
 * @param dt_back Time elapsed since this star's last HII rebuild pass.
 */
__attribute__((always_inline)) INLINE void
radiation_open_ionizing_photon_budget(struct spart *sp, double dt_back) {

  for (int p = 0; p < sp->feedback_data.radiation.n_HII_pixels; p++) {
    /* A pass overdraws its pixel by up to one particle's cost, since the
       boundary particle is claimed in full. Carry that debt forward instead of
       forgiving it: forgiven once per pass, it would over-issue photons in
       proportion to the number of passes, i.e. as 1/dt_back -- a cadence
       dependence. Unspent *positive* budget is not carried, those photons
       reached no gas and escaped. */
    const double debt =
        min(sp->feedback_data.radiation.N_ion_budget_pix[p], 0.);
    sp->feedback_data.radiation.N_ion_budget_pix[p] =
        debt + sp->feedback_data.radiation.dot_N_ion_pix[p] * dt_back;
  }
}

/**
 * Consume the #spart ionizing photon budget.
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
 * Tag the #part as ionized to be ionized in feedback_update_part().
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param star_id The id of the star that ionized this particle.
 * @param end_time The simulation time until which this particle should
 * stay flagged as ionized (the ionizing star's next HII rebuild) -- cooling
 * keeps re-flooring its temperature until then instead of undoing the
 * ionization on the very next step.
 */
__attribute__((always_inline)) INLINE void radiation_tag_part_as_ionized(
    struct part *p, struct xpart *xp, long long star_id, double end_time,
    float excess_photon_energy_HI, float photoionization_rate_HI) {
  p->feedback_data.is_ionized = 1;
  p->feedback_data.star_id = star_id;
  p->feedback_data.end_time = end_time;
  xp->tracers_data.HII_region.excess_photon_energy_HI = excess_photon_energy_HI;
  xp->tracers_data.HII_region.photoionization_rate_HI = photoionization_rate_HI;
  return;
}

/**
 * Reset the #part ionization tag.
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
 * Is this #part *tagged* as ionized ?
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
 * The simulation time until which this #part should stay flagged as
 * ionized. Only meaningful while radiation_is_part_tagged_as_ionized()
 * is true.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_ionized_end_time(const struct part *p,
                                    const struct xpart *xp) {
  return p->feedback_data.end_time;
}

/**
 * Id of the star that ionized this #part. Only meaningful while
 * radiation_is_part_tagged_as_ionized() is true.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE long long
radiation_get_part_ionized_star_id(const struct part *p,
                                   const struct xpart *xp) {
  return p->feedback_data.star_id;
}

/**
 * Mean photon energy above the 13.6 eV HI ionization threshold of the
 * star that tagged this #part, in cgs (erg), frozen at tag time. Only
 * meaningful while radiation_is_part_tagged_as_ionized() is true.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_excess_photon_energy_HI(const struct part *p,
                                           const struct xpart *xp) {
  return xp->tracers_data.HII_region.excess_photon_energy_HI;
}

/**
 * Photoionization rate coefficient Gamma_HI frozen on this #part at tag
 * time (internal 1/time). Only meaningful while
 * radiation_is_part_tagged_as_ionized() is true.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_photoionization_rate_coefficient(const struct part *p,
                                                    const struct xpart *xp) {
  return xp->tracers_data.HII_region.photoionization_rate_HI;
}

/**
 * Photoionization rate coefficient Gamma_HI from an HI-ionizing photon
 * flux (photons / area / time, internal units), via the standard hydrogen
 * photoionization cross-section at the Lyman limit (sigma_HI = 6.3e-18
 * cm^2, Osterbrock & Ferland 2006 -- a physical constant, not a tunable
 * parameter). Called once at tag time (feedback_iact_HII_ionization): the
 * raw flux is too large for float32 in this unit system, but the product
 * with the tiny cross-section is safely representable, so only that
 * product is stored.
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
