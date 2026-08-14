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
 * @file src/feedback/GEAR/radiation.c
 * @brief Subgrid radiation feedback for GEAR. This files contains functions to
 * compute quantities for the radiation feedback.
 */

/* Include header */
#include "radiation.h"

#include "chemistry.h"
#include "cooling.h"
#include "equation_of_state.h"
#include "interpolation.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"
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
  double X_HI;
#ifdef WITH_MPI
  if (xp == NULL) {
    /* Foreign copy: struct part's feedback_data.neutral_H_frac already
       folds in the H2/H- arms at mode >= 2 (see cooling_cache_neutral_H_
       fraction_subgrid), so this is the same quantity the branch below
       would read from a real xpart. */
    X_HI = p->feedback_data.neutral_H_frac;
  } else
#endif
  {
    const struct cooling_xpart_data *cool_data = &xp->cooling_data;
    X_HI = cool_data->HI_frac;
#if COOLING_GRACKLE_MODE >= 2
    /* Hydrogen locked in H2/H- is also neutral (not yet stripped); H2II is
       already singly-ionized, so it is excluded -- the same H2II
       approximation radiation_get_part_total_hydrogen_mass_fraction's own
       doxygen already notes for the total-hydrogen accounting. */
    X_HI += cool_data->H2I_frac + cool_data->HM_frac;
#endif
  }

  const float m = hydro_get_mass(p);
  const double m_p = phys_const->const_proton_mass;
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
 * Mean molecular weight of a #part, foreign-copy safe.
 *
 * cooling_get_mean_molecular_weight reads xp->cooling_data species
 * fractions at Grackle modes >= 1, which segfaults for a foreign copy
 * with no struct xpart. This accessor owns that xp == NULL contract:
 * a real xpart computes exactly as today. xp == NULL inverts the same
 * relation cooling_cache_eligibility_temperature_subgrid used to write
 * struct part's feedback_data.T_eligibility, mu = k_B T / ((gamma-1) u
 * m_p), using the particle's current internal energy (rides the normal
 * xv/rho exchange). T_eligibility and u are always written from the same
 * cooling_new_energy call and shipped together whenever the owner-side
 * state actually gets read by a foreign gather, so this recovers the
 * owner's exact mu rather than approximating it, exactly like the
 * temperature cache itself.
 *
 * @param phys_const Physical constants.
 * @param hydro_properties The #hydro_props.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle, or NULL for a foreign
 * copy with no local xpart.
 * @return Mean molecular weight.
 */
__attribute__((always_inline)) INLINE double
radiation_get_part_mean_molecular_weight(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

#ifdef WITH_MPI
  if (xp == NULL) {
    const double u = hydro_get_drifted_physical_internal_energy(p, cosmo);
    const double T = p->feedback_data.T_eligibility;
    return phys_const->const_boltzmann_k * T /
           (hydro_gamma_minus_one * u * phys_const->const_proton_mass);
  }
#endif

  return cooling_get_mean_molecular_weight(phys_const, us, cosmo, hydro_props,
                                           cooling, p, xp);
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
  const double mu = radiation_get_part_mean_molecular_weight(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);

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

  /* Z >= HydrogenFractionByMass clamps N_H to 0 (radiation.c's X_H clamp);
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
  const double mu = radiation_get_part_mean_molecular_weight(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);
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

/**
 * Compute the gas comoving column density at the star's location using the
 * Sobolev approximation.
 *
 * @param sp The #spart.
 * @return Comoving gas column density at the star's location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_star(const struct spart *sp) {
  /* enrichment_weight is the star's SPH-averaged local gas density. */
  const float rho_gas = sp->feedback_data.enrichment_weight;
  const float grad_rho[3] = {sp->feedback_data.grad_rho_star[0],
                             sp->feedback_data.grad_rho_star[1],
                             sp->feedback_data.grad_rho_star[2]};
  const float norm_grad_rho =
      sqrtf(grad_rho[0] * grad_rho[0] + grad_rho[1] * grad_rho[1] +
            grad_rho[2] * grad_rho[2]);

  /* A locally uniform density field (zero gradient, e.g. unperturbed
     glass/grid ICs) makes the Sobolev length rho/|grad rho| undefined;
     fall back to just the kernel support radius in that case. */
  const float sobolev_length =
      norm_grad_rho > 0.0f ? rho_gas / norm_grad_rho : 0.0f;
  const float length_gas = sp->h * kernel_gamma + sobolev_length;
  return length_gas * rho_gas;
}

/**
 * Compute the physical infrared opacity around a star.
 *
 * @param sp The #spart.
 * @param us Unit system.
 * @return Infrared gas opacity around the star.
 */
__attribute__((always_inline)) INLINE float radiation_get_physical_IR_opacity(
    const struct spart *sp, const struct unit_system *us) {
  const float Z_gas = sp->feedback_data.Z_star;
  const float Z_sun = 0.02;
  const float value = 10.0 * units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
                      units_cgs_conversion_factor(us, UNIT_CONV_AREA);
  return value * Z_gas / Z_sun;
}

/**
 * Compute the physical infrared optical depth around a star.
 *
 * @param sp The #spart.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @return Infrared gas optical depth around the star.
 */
__attribute__((always_inline)) INLINE float
radiation_get_physical_IR_optical_depth(const struct spart *sp,
                                        const struct unit_system *us,
                                        const struct cosmology *cosmo) {
  const float Sigma_gas_c =
      radiation_get_comoving_gas_column_density_at_star(sp);
  const float Sigma_gas_p = Sigma_gas_c * cosmo->a2_inv;
  const float kappa_IR = radiation_get_physical_IR_opacity(sp, us);
  return kappa_IR * Sigma_gas_p;
}

/**
 * Compute the physical radiation pressure emitted by the star.
 *
 * @param sp The #spart.
 * @param Delta_t The current #spart timestep.
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @return Radiation pressure emittied by the star.
 */
__attribute__((always_inline)) INLINE float
radiation_get_star_physical_radiation_pressure(
    const struct spart *sp, const float Delta_t,
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo) {

  const float tau_IR = radiation_get_physical_IR_optical_depth(sp, us, cosmo);
  const float L_bol = sp->feedback_data.radiation.L_bol;
  const float c = phys_const->const_speed_light_c;

  return Delta_t * L_bol / c * (1 + tau_IR);
}

/**
 * Compute the radius of a single star from empirical mass-radius relations.
 *
 * This function gets the value for an individual star. For a SSP, this
 * function is used to compute an IMF-average.
 *
 * @param mass Mass of the star.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Radius in code units.
 */
float radiation_get_individual_star_radius(
    const float mass, const struct unit_system *us,
    const struct phys_const *phys_const) {

  /* Perform some units conversions */
  const float R_sun = phys_const->const_solar_radius;
  const float M_solar = phys_const->const_solar_mass;
  const float M_in_solar = mass / M_solar;

  if (M_in_solar < 1.f) {
    return R_sun * powf(M_in_solar, 0.8f);
  } else if (M_in_solar < 8.f) {
    return R_sun * powf(M_in_solar, 0.57f);
  } else {
    return R_sun * powf(M_in_solar, 0.5f);
  }
}

/**
 * Compute the temperature of a single star from empirical mass-temperature
 * relations.
 *
 * This function gets the value for an individual star. For a SSP, this
 * function is used to compute an IMF-average.
 *
 * @param mass Mass of the star.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Temperature in code units.
 */
float radiation_get_individual_star_temperature(
    const float mass, const struct unit_system *us,
    const struct phys_const *phys_const) {

  const float M_solar = phys_const->const_solar_mass;
  const float M_in_solar = mass / M_solar; /* In solar masses */

  float T_K = 0.0;

  if (M_in_solar < 1.f) {
    T_K = 3500.f * powf(M_in_solar, 0.5f);
  } else if (M_in_solar < 8.f) {
    T_K = 5800.f * powf(M_in_solar, 0.5f);
  } else {
    T_K = 25000.f * powf(M_in_solar / 20.f, 0.1f);
  }

  /* Convert from Kelvin to internal units using unit_system_temperature_in_cgs
   */
  const float T_internal =
      T_K / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  return T_internal;
}

/**
 * Computes the bolometric luminosity of a single star from empirical
 * mass-luminosity relations.
 *
 * This function gets the value for an individual star. For a SSP, this
 * function is used to compute an IMF-average.
 *
 * @param mass Mass of the star.
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @return Luminosity in code units.
 */
float radiation_get_individual_star_luminosity(
    const float mass, const struct unit_system *us,
    const struct phys_const *phys_const) {

  /* Convert mass to solar masses */
  const float M_in_solar = mass / phys_const->const_solar_mass;

  /* Piecewise empirical mass-luminosity relation */
  float lum_sol;
  if (M_in_solar < 0.43f) {
    lum_sol = 0.185f * M_in_solar * M_in_solar;
  } else if (M_in_solar < 2.0f) {
    lum_sol = M_in_solar * M_in_solar * M_in_solar * M_in_solar;
  } else if (M_in_solar < 54.0f) {
    lum_sol = 1.5f * M_in_solar * M_in_solar * M_in_solar * sqrtf(M_in_solar);
  } else {
    lum_sol = 32000.0f * M_in_solar;
  }

  /* Convert from solar luminosities to code units */
  const float luminosity = lum_sol * phys_const->const_solar_luminosity;
  return luminosity;
}

/**
 * @brief Get the #spart ionizing photon emission rate using an analytical
 * series expansion of the Blackbody spectrum.
 *
 * This provides a physically justified, highly accurate formulation for an
 * individual stellar source without relying on crude empirical fitting curves
 * or fixed average photon energy assumptions.
 *
 * @param mass Mass of the star particle.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @return N_dot_ion The ionizing photon emission rate in code units
 * [photons/U_T].
 */
double radiation_get_individual_star_ionizing_photon_emission_rate_fit(
    const float mass, const struct unit_system *us,
    const struct phys_const *phys_const) {

  /* Get star properties in internal units */
  const float R = radiation_get_individual_star_radius(mass, us, phys_const);
  const float L =
      radiation_get_individual_star_luminosity(mass, us, phys_const);

  if (R <= 0.f || L <= 0.f) {
    return 0.0;
  }

  const float R_in_R_sun = R / phys_const->const_solar_radius;
  const float L_in_L_sun = L / phys_const->const_solar_luminosity;

  /* Get the Blackbody effective temperature in K */
  const double T_K =
      5780.0 * pow((double)(L_in_L_sun / (R_in_R_sun * R_in_R_sun)), 0.25) /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Compute dimensionless photon cutoff x_0 = h*nu_0 / k_B T for 13.605 eV */
  const double E_threshold_internal = 13.605 * phys_const->const_electron_volt;
  const double x_0 =
      E_threshold_internal / (phys_const->const_boltzmann_k * T_K);

  /* If x_0 is highly elevated, the stellar temperature is too low to produce
     any significant UV-ionizing radiation. (e.g., x_0 > 45 means exp(-x_0) <
     1e-20) */
  if (x_0 > 45.0) {
    return 0.0;
  }

  /* Evaluate the integral using a fast-converging series expansion
     Integral(x^2 / (e^x - 1)) =
                  Sum_{n=1}^inf [ e^(-n*x_0)/n * (x_0^2 + 2x_0/n + 2/n^2) ] */
  double photon_integral_sum = 0.0;
  const int max_terms = 5;

  for (int n = 1; n <= max_terms; ++n) {
    const double exp_term = exp(-((double)n) * x_0);

    /* Break early if subsequent terms underflow our interest bounds */
    if (exp_term < 1e-10) {
      break;
    }

    const double n_double = (double)n;
    const double term =
        (exp_term / n_double) *
        (x_0 * x_0 + (2.0 * x_0) / n_double + 2.0 / (n_double * n_double));
    photon_integral_sum += term;
  }

  /* Prefactor for the blackbody photon number density fraction. Normalized
   * via 15 / pi^4 */
  const double prefactor = 15.0 / (M_PI * M_PI * M_PI * M_PI);

  /* Total photon production rate above the ionization edge:
                 N_dot = (L / (k_B * T)) * (15 / pi^4) * photon_integral_sum
  */
  const double N_dot_ion = (L / (phys_const->const_boltzmann_k * T_K)) *
                           prefactor * photon_integral_sum;

  return N_dot_ion;
}

/**
 * @brief Get the star's mean excess photon energy above the 13.6 eV
 * hydrogen ionization threshold, from the same blackbody spectrum used by
 * radiation_get_individual_star_ionizing_photon_emission_rate_fit.
 *
 * Only used when GEARFeedback:HII_couple_ionization_rate is on, to derive
 * Grackle's RT_heating_rate for a rate-coupled particle
 * (RT_heating_rate = Gamma_HI * mean_excess_photon_energy_HI).
 *
 * @param mass Mass of the star particle.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @return Mean excess photon energy above 13.6 eV, in cgs (erg) -- not
 * internal units, since this project's internal mass unit makes the
 * absolute per-particle value underflow float precision (see caller).
 */
double radiation_get_individual_star_mean_excess_photon_energy_HI(
    const float mass, const struct unit_system *us,
    const struct phys_const *phys_const) {

  const float R = radiation_get_individual_star_radius(mass, us, phys_const);
  const float L =
      radiation_get_individual_star_luminosity(mass, us, phys_const);

  if (R <= 0.f || L <= 0.f) {
    return 0.0;
  }

  const float R_in_R_sun = R / phys_const->const_solar_radius;
  const float L_in_L_sun = L / phys_const->const_solar_luminosity;

  const double T_K =
      5780.0 * pow((double)(L_in_L_sun / (R_in_R_sun * R_in_R_sun)), 0.25) /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  const double E_threshold_internal = 13.605 * phys_const->const_electron_volt;
  const double x_0 =
      E_threshold_internal / (phys_const->const_boltzmann_k * T_K);

  if (x_0 > 45.0) {
    return 0.0;
  }

  const int max_terms = 5;

  /* Photon NUMBER integral, Integral(x^2/(e^x-1))dx (same series as the
     ionizing photon rate fit above). */
  double number_integral_sum = 0.0;
  /* Photon ENERGY integral, Integral(x^3/(e^x-1))dx -- one moment higher,
     by the identical by-parts derivation. */
  double energy_integral_sum = 0.0;

  for (int n = 1; n <= max_terms; ++n) {
    const double exp_term = exp(-((double)n) * x_0);
    if (exp_term < 1e-10) {
      break;
    }
    const double n_double = (double)n;

    number_integral_sum +=
        (exp_term / n_double) *
        (x_0 * x_0 + (2.0 * x_0) / n_double + 2.0 / (n_double * n_double));

    energy_integral_sum += (exp_term / n_double) *
                           (x_0 * x_0 * x_0 + (3.0 * x_0 * x_0) / n_double +
                            (6.0 * x_0) / (n_double * n_double) +
                            6.0 / (n_double * n_double * n_double));
  }

  if (number_integral_sum <= 0.0) {
    return 0.0;
  }

  /* Mean photon energy above threshold, in units of k_B*T: ratio of the
     energy-weighted to the number-weighted integral. */
  const double mean_hnu_over_kT = energy_integral_sum / number_integral_sum;
  const double E_excess_internal =
      phys_const->const_boltzmann_k * T_K * (mean_hnu_over_kT - x_0);

  /* Return in cgs, not internal units: the internal mass unit (1e10 Msun)
     makes this absolute per-particle energy ~1e-65 internally, underflowing
     to exactly 0 once narrowed to float at the caching site. The cgs value
     (~1e-11 erg) is safely representable. */
  return E_excess_internal * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
}

/******************************************************************************/
/* Functions to deal with integrated data over an IMF. These functions read,
   interpolate and integrate. */
/******************************************************************************/

/**
 * @brief Print the radiation model.
 *
 * @param rad The #radiation.
 */
void radiation_print(const struct radiation *rad) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  message("Angular pixels for HII ionization = %d", rad->n_HII_pixels);
  message("Interpolation table size = %d", rad->interpolation_size);
}

/**
 * @brief Initialize the #radiation structure.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_init(struct radiation *rad, struct swift_params *params,
                    const struct stellar_model *sm,
                    const struct unit_system *us,
                    const struct phys_const *phys_const) {

  /* Read the data */
  radiation_read_data(rad, params, sm, us, phys_const, /* restart */ 0);

  /* Angular (HEALPix) splitting of the HII ionization budget. nside=0 means
     spherical (HEALPix disabled, today's behaviour, n_HII_pixels=1); any
     nside>=1 means the standard HEALPix RING-scheme tessellation
     (n_HII_pixels=12*nside^2) -- RING, unlike NEST, has no power-of-2
     restriction on nside, so any positive integer is mathematically valid
     here (see /usr/include/chealpix.h). The practical ceiling is memory,
     not geometry: every star carries a fixed-size
     dot_N_ion_pix[HII_MAX_ANGULAR_PIXELS] array
     (src/feedback/GEAR_thermal/feedback_struct.h) sized by
     ./configure --with-number-of-hii-angular-pixels (default 12, i.e.
     nside<=1); a run requesting more pixels than that build was
     configured for errors clearly below rather than overflowing the
     array. */
  const int nside =
      parser_get_opt_param_int(params, "GEARFeedback:HII_angular_nside", 0);
  if (nside < 0) {
    error("GEARFeedback:HII_angular_nside must be >= 0; got %d.", nside);
  }
  const int n_HII_pixels_requested = (nside == 0) ? 1 : 12 * nside * nside;
  if (n_HII_pixels_requested > HII_MAX_ANGULAR_PIXELS) {
    error(
        "GEARFeedback:HII_angular_nside=%d requires %d HealPix pixels, but "
        "this build only supports up to HII_MAX_ANGULAR_PIXELS=%d. "
        "Reconfigure with "
        "--with-number-of-hii-angular-pixels=%d (or higher) and rebuild, "
        "or lower HII_angular_nside.",
        nside, n_HII_pixels_requested, HII_MAX_ANGULAR_PIXELS,
        n_HII_pixels_requested);
  }
#ifndef HAVE_CHEALPIX
  if (nside != 0) {
    error(
        "GEARFeedback:HII_angular_nside > 0 requires the HEALPix C API "
        "(chealpix). Reconfigure with --with-chealpix, or set nside=0.");
  }
#endif
  rad->n_HII_pixels = n_HII_pixels_requested;
}

/**
 * @brief Write a radiation struct to the given FILE as a stream of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void radiation_dump(const struct radiation *rad, FILE *stream,
                    const struct stellar_model *sm) {

  restart_write_blocks((void *)rad, sizeof(struct radiation), 1, stream,
                       "radiation", "radiation");
  message("Dumping GEAR radiation...");
}

/**
 * @brief Restore a radiation struct from the given FILE as a stream of
 * bytes.
 *
 * The flat restore below copies the interpolation tables' internal data
 * pointers as raw bytes -- meaningless in the new process, since they held
 * the old process's heap addresses. radiation_read_data() re-derives those
 * tables from scratch instead of trying to serialize them (they are
 * computed from mass/Z, not read from a file, so re-deriving is exact and
 * avoids ever leaving a dangling pointer for radiation_clean() to free().
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 */
void radiation_restore(struct radiation *rad, FILE *stream,
                       const struct stellar_model *sm,
                       const struct unit_system *us,
                       const struct phys_const *phys_const) {

  restart_read_blocks((void *)rad, sizeof(struct radiation), 1, stream, NULL,
                      "radiation");
  radiation_read_data(rad, NULL, sm, us, phys_const, /*restart=*/1);
  message("Restoring GEAR radiation struct...");
}

/**
 * @brief Clean the allocated memory.
 *
 * @param rad the #radiation.
 */
void radiation_clean(struct radiation *rad) {

  interpolate_1d_free(&rad->integrated.luminosities);
  interpolate_1d_free(&rad->raw.luminosities);
  interpolate_1d_free(&rad->integrated.dot_N_ion);
  interpolate_1d_free(&rad->raw.dot_N_ion);
  interpolate_1d_free(&rad->integrated.dot_E_excess);
  interpolate_1d_free(&rad->raw.dot_E_excess);
}

/**
 * @brief Get the IMF-averaged nolometric luminosity per mass.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The bolometric luminosity.
 */
float radiation_get_luminosities_from_integral(const struct radiation *rad,
                                               float log_m1, float log_m2) {

  float luminosity_1 = interpolate_1d(&rad->integrated.luminosities, log_m1);
  float luminosity_2 = interpolate_1d(&rad->integrated.luminosities, log_m2);
  return luminosity_2 - luminosity_1;
};

/**
 * @brief Get the IMF-averaged bolometric luminosity per mass.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param The bolometric luminosity.
 */
float radiation_get_luminosities_from_raw(const struct radiation *rad,
                                          float log_m) {
  return interpolate_1d(&rad->raw.luminosities, log_m);
};

/**
 * @brief Get the IMF-averaged ionization rate per mass.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The ionization rate;
 */
double radiation_get_ionization_rate_from_integral(const struct radiation *rad,
                                                   float log_m1, float log_m2) {

  double dot_N_ion_1 = interpolate_1d(&rad->integrated.dot_N_ion, log_m1) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  double dot_N_ion_2 = interpolate_1d(&rad->integrated.dot_N_ion, log_m2) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  return dot_N_ion_2 - dot_N_ion_1;
};

/**
 * @brief Get the non-IMF-integrated ionization rate per mass.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param The ionization rate;
 */
double radiation_get_ionization_rate_from_raw(const struct radiation *rad,
                                              float log_m) {
  return interpolate_1d(&rad->raw.dot_N_ion, log_m) *
         RADIATION_DOT_N_ION_TABLE_SCALING;
};

/**
 * @brief Get the IMF-averaged, Q-weighted mean excess photon energy above
 * the 13.6 eV HI ionization threshold, for a population over a mass window.
 *
 * Ratio of the integrated dot_E_excess and dot_N_ion tables, taken directly
 * on the raw (still /RADIATION_DOT_N_ION_TABLE_SCALING) interpolated
 * values rather than through their public accessors: the scaling constant
 * multiplies both tables identically, so it cancels in the ratio without
 * ever needing to be undone, leaving a result in cgs erg (see
 * #radiation_get_individual_star_mean_excess_photon_energy_HI).
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return Q-weighted mean excess photon energy in cgs erg, or 0 if no
 * ionizing photons are produced over the window (dot_N_ion difference is
 * 0 -- e.g. no alive ionizing stars).
 */
double radiation_get_mean_excess_photon_energy_HI_from_integral(
    const struct radiation *rad, float log_m1, float log_m2) {

  const double dot_N_ion_1 = interpolate_1d(&rad->integrated.dot_N_ion, log_m1);
  const double dot_N_ion_2 = interpolate_1d(&rad->integrated.dot_N_ion, log_m2);
  const double delta_dot_N_ion = dot_N_ion_2 - dot_N_ion_1;

  /* The cumulative table is monotonically non-decreasing in mass, so any
     non-positive difference (no alive ionizing stars in this window, a
     zero-width window, or roundoff noise between two nearly-equal table
     entries) is degenerate -- guard against dividing by it rather than
     testing for exact 0, which a near-cancellation could slip past and
     amplify into a meaningless huge or negative result. */
  if (delta_dot_N_ion <= 0.) return 0.;

  const double dot_E_excess_1 =
      interpolate_1d(&rad->integrated.dot_E_excess, log_m1);
  const double dot_E_excess_2 =
      interpolate_1d(&rad->integrated.dot_E_excess, log_m2);
  const double delta_dot_E_excess = dot_E_excess_2 - dot_E_excess_1;

  return delta_dot_E_excess / delta_dot_N_ion;
};

/**
 * @brief Read an array of luminosities data from the table.
 *
 * @param rad The #radiation model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm * The #stellar_model.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size Number of element to keep in the interpolation
 * data.
 */
void radiation_read_luminosities_array(struct radiation *rad,
                                       struct interpolation_1d *interp_raw,
                                       struct interpolation_1d *interp_int,
                                       const struct stellar_model *sm,
                                       int interpolation_size,
                                       const struct unit_system *us,
                                       const struct phys_const *phys_const) {

  /* Allocate the memory */
  const int count = 500;
  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for luminosities.");

  const float mass_min = sm->imf.mass_min;
  const float mass_max = sm->imf.mass_max;
  const float log_mass_min = log10f(mass_min);
  const float log_mass_max = log10f(mass_max);
  const float step_size = (log_mass_max - log_mass_min) / (count - 1);

  /* Fill the table */
  for (size_t j = 0; j < count; j++) {
    /* Compute the log-mass and mass */
    const float log_mass = log_mass_min + j * step_size;
    const float mass = exp10(log_mass) * phys_const->const_solar_mass;

    /* Get bolometric luminosity for this mass, in internal units */
    data[j] = radiation_get_individual_star_luminosity(mass, us, phys_const);
  }

  /* Initialize the raw interpolation */
  interpolate_1d_init(interp_raw, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_error);

  initial_mass_function_integrate(&sm->imf, data, count, log_mass_min,
                                  step_size);
  // TODO: decrease count in order to keep the same distance between points

  /* Initialize the integrated interpolation */
  interpolate_1d_init(interp_int, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_const);

  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read an array of ionizing emission rates data from the table.
 *
 * @param rad The #radiation model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm * The #stellar_model.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size Number of element to keep in the interpolation
 * data.
 */
void radiation_read_ionization_rate_array(struct radiation *rad,
                                          struct interpolation_1d *interp_raw,
                                          struct interpolation_1d *interp_int,
                                          const struct stellar_model *sm,
                                          int interpolation_size,
                                          const struct unit_system *us,
                                          const struct phys_const *phys_const) {

  /* Allocate the memory */
  const int count = 500;
  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for luminosities.");

  const float mass_min = sm->imf.mass_min;
  const float mass_max = sm->imf.mass_max;
  const float log_mass_min = log10f(mass_min);
  const float log_mass_max = log10f(mass_max);
  const float step_size = (log_mass_max - log_mass_min) / (count - 1);

  /* Fill the table */
  for (size_t j = 0; j < count; j++) {
    /* Compute the log-mass and mass */
    const float log_mass = log_mass_min + j * step_size;
    const float mass = exp10(log_mass) * phys_const->const_solar_mass;

    /* Get bolometric luminosity for this mass, in internal units */
    data[j] = radiation_get_individual_star_ionizing_photon_emission_rate_fit(
                  mass, us, phys_const) /
              RADIATION_DOT_N_ION_TABLE_SCALING;
  }

  /* Initialize the raw interpolation */
  interpolate_1d_init(interp_raw, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_error);

  initial_mass_function_integrate(&sm->imf, data, count, log_mass_min,
                                  step_size);
  // TODO: decrease count in order to keep the same distance between points

  /* Initialize the integrated interpolation */
  interpolate_1d_init(interp_int, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_const);

  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read an array of excess-photon-energy emission rate data from the
 * table: dot_E_excess(m) = dot_N_ion(m) * mean_excess_photon_energy_HI(m).
 *
 * IMF-averaging this product (rather than mean_excess_photon_energy_HI(m)
 * alone) is what makes the eventual ratio of integrated tables a
 * Q-weighted mean: a star that contributes more ionizing photons should
 * weigh more in the population's mean excess energy. Divided by
 * RADIATION_DOT_N_ION_TABLE_SCALING for the same reason as the dot_N_ion
 * table (dot_N_ion(m) alone already needs it; the product would otherwise
 * overflow float storage) -- the same constant multiplies both raw tables,
 * so it cancels exactly in
 * #radiation_get_mean_excess_photon_energy_HI_from_integral's ratio.
 *
 * @param rad The #radiation model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm The #stellar_model.
 * @param interpolation_size Number of element to keep in the interpolation
 * data.
 */
void radiation_read_mean_excess_photon_energy_array(
    struct radiation *rad, struct interpolation_1d *interp_raw,
    struct interpolation_1d *interp_int, const struct stellar_model *sm,
    int interpolation_size, const struct unit_system *us,
    const struct phys_const *phys_const) {

  /* Allocate the memory */
  const int count = 500;
  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for excess photon energy.");

  const float mass_min = sm->imf.mass_min;
  const float mass_max = sm->imf.mass_max;
  const float log_mass_min = log10f(mass_min);
  const float log_mass_max = log10f(mass_max);
  const float step_size = (log_mass_max - log_mass_min) / (count - 1);

  /* Fill the table */
  for (size_t j = 0; j < count; j++) {
    /* Compute the log-mass and mass */
    const float log_mass = log_mass_min + j * step_size;
    const float mass = exp10(log_mass) * phys_const->const_solar_mass;

    const double dot_N_ion =
        radiation_get_individual_star_ionizing_photon_emission_rate_fit(
            mass, us, phys_const);
    const double E_excess =
        radiation_get_individual_star_mean_excess_photon_energy_HI(mass, us,
                                                                   phys_const);
    data[j] = (float)(dot_N_ion * E_excess / RADIATION_DOT_N_ION_TABLE_SCALING);
  }

  /* Initialize the raw interpolation */
  interpolate_1d_init(interp_raw, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_error);

  initial_mass_function_integrate(&sm->imf, data, count, log_mass_min,
                                  step_size);

  /* Initialize the integrated interpolation */
  interpolate_1d_init(interp_int, log_mass_min, log_mass_max,
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_const);

  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read the RAD yields from the table.
 *
 * The tables are in internal units at the end of this function.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void radiation_read_data(struct radiation *rad, struct swift_params *params,
                         const struct stellar_model *sm,
                         const struct unit_system *us,
                         const struct phys_const *phys_const,
                         const int restart) {

  if (!restart) {
    /* TODO: Maybe update this */
    rad->interpolation_size = parser_get_opt_param_int(
        params, "GEARSupernovaeII:interpolation_size", 200);
  }

  /* Read the luminosities */
  radiation_read_luminosities_array(rad, &rad->raw.luminosities,
                                    &rad->integrated.luminosities, sm,
                                    rad->interpolation_size, us, phys_const);

  /* Read the ionization emission rates */
  radiation_read_ionization_rate_array(rad, &rad->raw.dot_N_ion,
                                       &rad->integrated.dot_N_ion, sm,
                                       rad->interpolation_size, us, phys_const);

  /* Read the excess-photon-energy emission rates */
  radiation_read_mean_excess_photon_energy_array(
      rad, &rad->raw.dot_E_excess, &rad->integrated.dot_E_excess, sm,
      rad->interpolation_size, us, phys_const);
};
