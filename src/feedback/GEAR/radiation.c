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
#include "engine.h"
#include "equation_of_state.h"
#include "hdf5_functions.h"
#include "interpolation.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"
#include "units.h"

#include <float.h>
#include <string.h>

/*! Floor applied before taking log10() of a Data/Radiation CGS value, so a
    genuinely-zero table entry (Q_H/DotEExcess below the source table's own
    ionization-threshold mass, where pychem defines them to be exactly 0 --
    see PyChemInitTable/libradiation.py) does not produce log10(0) = -inf.
    Matches pychem's own floor bit-for-bit (PyChemInitTable/
    libparsec_radiation.py's `_LOG_FLOOR`), per Darwin's explicit ruling that
    SWIFT's radiation interpolation must match pychem's own
    log10(mass)-vs-log10(value) interpolation scheme, not just approximate
    it. Applied to the CGS value itself (before the internal-unit
    conversion), not the converted value: this is the number pychem's own
    floor actually clamps, so a query that is "native-zero" in SWIFT means
    exactly what it means in pychem, independently of SWIFT's unit
    system. */
#define RADIATION_LOG_FLOOR_CGS 1e-300

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
  message("Interpolation table size (mass) = %d", rad->interpolation_size);
  if (rad->is_2d) {
    message("Interpolation table size (metallicity) = %d",
            rad->interpolation_size_metallicity);
  }
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
 * tables from scratch instead of trying to serialize them, avoiding ever
 * leaving a dangling pointer for radiation_clean() to free(). Re-derivation
 * reads sm->yields_table again (Data/Radiation) rather than recomputing
 * from mass/Z alone, so it is exact only if that path still resolves and
 * the file is unchanged since the run started -- the same uncanonicalized-
 * path caveat already noted for GEARFeedback:yields_table in general
 * (feedback_properties.h); a restart resubmitted from a different working
 * directory with a relative path can fail here.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param with_radiation Are we restoring with photoionization and/or
 * radiation pressure? The raw struct bytes are always read back
 * (radiation_dump() always writes them, unlike e.g. stellar_wind_dump()),
 * but the tables -- and #radiation.is_active -- are only re-derived, and
 * sm->yields_table only re-opened, when this is set; otherwise
 * #radiation_zero_pointers overwrites whatever stale value the raw restore
 * above just wrote into is_active.
 */
void radiation_restore(struct radiation *rad, FILE *stream,
                       const struct stellar_model *sm,
                       const struct unit_system *us,
                       const struct phys_const *phys_const,
                       const char with_radiation) {

  restart_read_blocks((void *)rad, sizeof(struct radiation), 1, stream, NULL,
                      "radiation");

  if (!with_radiation) {
    /* The bytes just restored are another process's heap addresses (see the
       function's own doxygen above); never dereference or free them. */
    radiation_zero_pointers(rad);
    return;
  }

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

  interpolate_2d_free(&rad->raw.luminosities_2d);
  interpolate_2d_free(&rad->raw.dot_N_ion_2d);
  interpolate_2d_free(&rad->raw.dot_E_excess_2d);
}

/**
 * @brief Zero a #radiation struct -- pointers, dimensions and #is_active --
 * so a struct that was never (or is not yet) initialized can be safely
 * passed to #radiation_clean, printed, or read by any of the getters (which
 * must check #is_active first; the getters themselves do not, and would
 * dereference the zeroed pointers below). Also called internally by
 * #radiation_read_data before it (re)builds the tables; see the scalar
 * save/restore around that call for why is_2d is the only field here it can
 * leave zeroed going in.
 *
 * @param rad The #radiation.
 */
void radiation_zero_pointers(struct radiation *rad) {

  rad->is_active = 0;
  rad->is_2d = 0;
  rad->interpolation_size = 0;
  rad->interpolation_size_metallicity = 0;
  rad->n_HII_pixels = 0;

  interpolate_1d_zero_pointers(&rad->raw.luminosities);
  interpolate_1d_zero_pointers(&rad->raw.dot_N_ion);
  interpolate_1d_zero_pointers(&rad->raw.dot_E_excess);
  interpolate_2d_zero_pointers(&rad->raw.luminosities_2d);
  interpolate_2d_zero_pointers(&rad->raw.dot_N_ion_2d);
  interpolate_2d_zero_pointers(&rad->raw.dot_E_excess_2d);
  interpolate_1d_zero_pointers(&rad->integrated.luminosities);
  interpolate_1d_zero_pointers(&rad->integrated.dot_N_ion);
  interpolate_1d_zero_pointers(&rad->integrated.dot_E_excess);
}

/**
 * @brief Abort with a clear message if #rad holds a 2D (mass x
 * metallicity) table: the getters below are 1D-only for now, since no
 * caller passes a metallicity yet.
 *
 * @param rad The #radiation model.
 * @param caller Name of the calling getter, for the error message.
 */
__attribute__((always_inline)) INLINE static void radiation_check_is_1d(
    const struct radiation *rad, const char *caller) {
  if (rad->is_2d) {
    error(
        "%s does not support a mass x metallicity (2D) radiation table yet: "
        "no caller passes a metallicity.",
        caller);
  }
}

/**
 * @brief Get the IMF-averaged nolometric luminosity per mass.
 *
 * Reads #rad->integrated.luminosities directly (no exponentiation): unlike
 * the raw table, the IMF-integrated table stays in linear value space --
 * see radiation_build_tables()'s doxygen for why pychem's log-log
 * convention does not apply to this SWIFT-side cumulative-integral
 * quantity.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The bolometric luminosity.
 */
float radiation_get_luminosities_from_integral(const struct radiation *rad,
                                               float log_m1, float log_m2) {

  radiation_check_is_1d(rad, __func__);
  float luminosity_1 = interpolate_1d(&rad->integrated.luminosities, log_m1);
  float luminosity_2 = interpolate_1d(&rad->integrated.luminosities, log_m2);
  return luminosity_2 - luminosity_1;
};

/**
 * @brief Get the non-IMF-integrated bolometric luminosity at a given mass.
 *
 * #rad->raw.luminosities holds log10(luminosity), pychem's own
 * log10(mass)-vs-log10(value) convention (see radiation_build_tables()'s
 * doxygen); the interpolated log-value is exponentiated back here. The
 * narrowing to float happens implicitly on return (exp10() itself returns
 * double); Luminosity is never 0 in this table (unlike Q_H/DotEExcess), so
 * there is no floor-underflow case to worry about here.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return The bolometric luminosity, internal units.
 */
float radiation_get_luminosities_from_raw(const struct radiation *rad,
                                          float log_m) {
  radiation_check_is_1d(rad, __func__);
  return (float)exp10(interpolate_1d(&rad->raw.luminosities, log_m));
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

  radiation_check_is_1d(rad, __func__);
  double dot_N_ion_1 = interpolate_1d(&rad->integrated.dot_N_ion, log_m1) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  double dot_N_ion_2 = interpolate_1d(&rad->integrated.dot_N_ion, log_m2) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  return dot_N_ion_2 - dot_N_ion_1;
};

/**
 * @brief Get the non-IMF-integrated ionization rate at a given mass.
 *
 * #rad->raw.dot_N_ion holds log10(dot_N_ion /
 * #RADIATION_DOT_N_ION_TABLE_SCALING in internal units) (see
 * radiation_build_tables()'s doxygen for the log-log storage convention);
 * exp10() undoes the log, exactly recovering the pre-log-transform
 * (already-scaled-down) internal value, then the existing
 * *RADIATION_DOT_N_ION_TABLE_SCALING undoes the scaling as before. The
 * narrowing to `float` below is load-bearing, not cosmetic: near
 * #RADIATION_LOG_FLOOR_CGS (1e-300, applied before the internal-unit
 * conversion divides it further down), float32's underflow floor
 * (~1e-45) is reached well before float64's, so this narrowing is what
 * makes a below-ionization-threshold query reliably return exactly
 * 0.0f -- the exact query mass at which that happens depends on the run's
 * own unit system. See radiation_read_cgs_array()'s own doxygen for the
 * reasoning behind #RADIATION_LOG_FLOOR_CGS's specific value; this getter
 * is why it needs to be that extreme.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return The ionization rate, internal units.
 */
double radiation_get_ionization_rate_from_raw(const struct radiation *rad,
                                              float log_m) {
  radiation_check_is_1d(rad, __func__);
  const float dot_N_ion_scaled =
      (float)exp10(interpolate_1d(&rad->raw.dot_N_ion, log_m));
  return (double)dot_N_ion_scaled * RADIATION_DOT_N_ION_TABLE_SCALING;
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
 * #radiation_read_mean_excess_photon_energy_array's own doxygen for why
 * that unit convention holds).
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

  radiation_check_is_1d(rad, __func__);
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
 * @brief Get the non-IMF-integrated mean excess photon energy above the
 * 13.6 eV HI ionization threshold, for a single star of a given mass.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_integral on the
 * raw (single-mass) tables: the ratio of the raw dot_E_excess and dot_N_ion
 * tables, taken directly on their still-/RADIATION_DOT_N_ION_TABLE_SCALING
 * values, for the same reason that ratio is exact and comes out in cgs erg
 * there. Both tables now hold log10(value) (see radiation_build_tables()'s
 * doxygen for the log-log storage convention); each is exponentiated
 * back, narrowed to `float`, before the ratio -- the same load-bearing
 * float32-underflow narrowing #radiation_get_ionization_rate_from_raw
 * uses to return exactly 0.0f below the table's native ionization
 * threshold, which is what makes the `dot_N_ion <= 0.` guard below
 * actually trigger there instead of dividing by a near-zero double.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return Mean excess photon energy in cgs erg, or 0 if this mass produces
 * no ionizing photons (dot_N_ion(log_m) <= 0).
 */
double radiation_get_mean_excess_photon_energy_HI_from_raw(
    const struct radiation *rad, float log_m) {

  radiation_check_is_1d(rad, __func__);
  const double dot_N_ion =
      (float)exp10(interpolate_1d(&rad->raw.dot_N_ion, log_m));

  /* Same degenerate-ratio guard as the integrated getter above. */
  if (dot_N_ion <= 0.) return 0.;

  const double dot_E_excess =
      (float)exp10(interpolate_1d(&rad->raw.dot_E_excess, log_m));
  return dot_E_excess / dot_N_ion;
};

/**
 * @brief Read a scalar HDF5 string attribute (fixed- or variable-length)
 * into a NUL-terminated buffer.
 *
 * SWIFT's generic io_read_attribute() (common_io.h) only supports the
 * numeric/char/bool #IO_DATA_TYPE variants; there is no read counterpart to
 * common_io.c's write-only io_writeStringAttribute() for a python/h5py
 * variable-length UTF-8 string, which is what pychem writes for
 * Data/Radiation's "dimensionality" attribute. This is the minimal reader
 * needed to dispatch on it.
 *
 * @param group_id Open HDF5 group id.
 * @param name Attribute name.
 * @param out Output buffer, NUL-terminated on return.
 * @param out_size Size of out, including the terminating NUL.
 */
static void radiation_read_string_attribute(hid_t group_id, const char *name,
                                            char *out, size_t out_size) {

  const hid_t h_attr = H5Aopen(group_id, name, H5P_DEFAULT);
  if (h_attr < 0) error("Error while opening attribute '%s'", name);

  const hid_t h_type = H5Aget_type(h_attr);
  if (h_type < 0) error("Error while getting the type of attribute '%s'", name);

  if (H5Tis_variable_str(h_type) > 0) {
    /* Read into the attribute's own native type (variable-length,
       H5T_CSET_UTF8, as pychem/h5py writes it) rather than a freshly
       crafted H5T_C_S1/H5T_VARIABLE memory type: the two differ in
       character set (ASCII vs UTF-8), and this HDF5 build has no
       registered ASCII<->UTF-8 conversion path, so H5Aread() into the
       mismatched type fails ("no appropriate function for conversion
       path") -- caught by actually running this against a real pychem
       table, not just compiling it. */
    char *tmp = NULL;
    if (H5Aread(h_attr, h_type, &tmp) < 0)
      error("Error while reading string attribute '%s'", name);

    const size_t len = strlen(tmp);
    if (len >= out_size) {
      error(
          "String attribute '%s' (%zu bytes) does not fit in the %zu-byte "
          "buffer.",
          name, len, out_size);
    }

    strncpy(out, tmp, out_size - 1);
    out[out_size - 1] = '\0';

    const hid_t h_space = H5Aget_space(h_attr);
    H5Dvlen_reclaim(h_type, h_space, H5P_DEFAULT, &tmp);
    H5Sclose(h_space);
  } else {
    const size_t fixed_size = H5Tget_size(h_type);
    if (fixed_size >= out_size)
      error(
          "String attribute '%s' (%zu bytes) does not fit in the %zu-byte "
          "buffer.",
          name, fixed_size, out_size);

    char *tmp = (char *)calloc(fixed_size + 1, sizeof(char));
    if (tmp == NULL) error("Failed to allocate string attribute buffer.");
    if (H5Aread(h_attr, h_type, tmp) < 0)
      error("Error while reading string attribute '%s'", name);
    memcpy(out, tmp, fixed_size);
    out[fixed_size] = '\0';
    free(tmp);
  }

  H5Tclose(h_type);
  H5Aclose(h_attr);
}

/**
 * @brief Assert that a Data/Radiation dataset's own "units" attribute
 * matches what this reader is about to assume.
 *
 * pychem writes an explicit "units" attribute on every dataset in
 * Data/Radiation specifically so a unit mismatch is never silently
 * implicit (PyChemInitTable/libradiation.py's own write_h5_table()
 * docstring: "Grackle's unit conventions are a known source of silent
 * errors, so units are never left implicit here"). SWIFT's reader
 * (radiation_read_cgs_array()) converts every dataset with a fixed,
 * hardcoded physical-dimension assumption (CGS erg/s for Luminosity, CGS
 * 1/s for Q_H, CGS erg/s for DotEExcess); this check makes that assumption
 * self-verifying against the table itself instead of trusting it blindly,
 * the same fail-loud-on-mismatch policy this file already applies to
 * float overflow (radiation_read_cgs_array()) and table-coverage (the
 * inline check in radiation_read_data()).
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset whose "units" attribute to check.
 * @param expected_units The unit string this reader is about to assume
 * (e.g. "erg/s"), compared verbatim (case-sensitive) against the table's
 * own attribute.
 */
static void radiation_check_dataset_units(hid_t group_id,
                                          const char *dataset_name,
                                          const char *expected_units) {

  const hid_t h_dataset = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
  if (h_dataset < 0)
    error("Error while opening dataset '%s' to check its units.", dataset_name);

  char actual_units[64];
  radiation_read_string_attribute(h_dataset, "units", actual_units,
                                  sizeof(actual_units));

  H5Dclose(h_dataset);

  if (strcmp(actual_units, expected_units) != 0) {
    error(
        "Data/Radiation/%s declares units='%s', but SWIFT's reader assumes "
        "'%s'. This table's unit convention no longer matches what this "
        "code converts -- aborting rather than silently misinterpreting "
        "the physics.",
        dataset_name, actual_units, expected_units);
  }
}

/**
 * @brief Turn one field's edge_policy_<field>_below/above attribute pair
 * into the #interpolate_boundary_condition SWIFT's 2D interpolator can
 * apply on the mass axis.
 *
 * pychem's schema allows "constant", "zero" or "linear" independently per
 * side (PyChemInitTable/libradiation.py's ExtrapolationPolicy); SWIFT has
 * no "linear" extrapolation mode, and #interpolate_boundary_condition can
 * only express "zero both sides", "zero below/constant above" or "constant
 * both sides" -- not "constant below/zero above". Every real shipped
 * table so far only uses the three representable combinations (see
 * radiation_parsec_popIII.hdf5's edge_policy_* attributes); this errors
 * loudly on anything else rather than silently misapplying a policy.
 *
 * @param below The field's edge_policy_<field>_below value.
 * @param above The field's edge_policy_<field>_above value.
 * @param field_name Name of the field these came from, for the error
 * message only.
 * @return The matching #interpolate_boundary_condition.
 */
static enum interpolate_boundary_condition radiation_parse_edge_policy(
    const char *below, const char *above, const char *field_name) {

  int zero_below = 0;
  if (strcmp(below, "zero") == 0) {
    zero_below = 1;
  } else if (strcmp(below, "constant") != 0) {
    error(
        "Data/Radiation field '%s': edge_policy_%s_below='%s' has no SWIFT "
        "#interpolate_boundary_condition equivalent (only \"zero\" and "
        "\"constant\" are supported).",
        field_name, field_name, below);
  }

  int zero_above = 0;
  if (strcmp(above, "zero") == 0) {
    zero_above = 1;
  } else if (strcmp(above, "constant") != 0) {
    error(
        "Data/Radiation field '%s': edge_policy_%s_above='%s' has no SWIFT "
        "#interpolate_boundary_condition equivalent (only \"zero\" and "
        "\"constant\" are supported).",
        field_name, field_name, above);
  }

  if (zero_below && zero_above) return boundary_condition_zero;
  if (zero_below && !zero_above) return boundary_condition_zero_const;
  if (!zero_below && !zero_above) return boundary_condition_const;

  error(
      "Data/Radiation field '%s': edge policy 'constant below / zero "
      "above' has no SWIFT #interpolate_boundary_condition equivalent.",
      field_name);
  return boundary_condition_error;
}

/**
 * @brief Read the Data/Radiation group's own grid metadata: the
 * "dimensionality" attribute ("M" or "M,Z") and the group-level
 * "m0"/"dm"/"nm" mass-grid attributes shared by every dataset in the
 * group, plus, for a 2D ("M,Z") table, "nz", the "Metallicity" dataset,
 * and the mass-axis edge_policy_* attributes (dispatched through "source"
 * for Q_H; see radiation_parse_edge_policy()).
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid (output) The #radiation_grid_metadata to fill in.
 */
static void radiation_read_grid_metadata(hid_t group_id,
                                         struct radiation_grid_metadata *grid) {

  radiation_read_string_attribute(group_id, "dimensionality",
                                  grid->dimensionality,
                                  sizeof(grid->dimensionality));

  io_read_attribute(group_id, "m0", FLOAT, &grid->log_mass_min);
  io_read_attribute(group_id, "dm", FLOAT, &grid->mass_step);
  io_read_attribute(group_id, "nm", INT, &grid->n_mass);

  if (!(grid->mass_step > 0.f))
    error(
        "Data/Radiation's 'dm' attribute is %.4g; it must be a strictly "
        "positive mass-grid step.",
        (double)grid->mass_step);

  if (grid->n_mass < 2)
    error(
        "Data/Radiation's 'nm' attribute is %d; at least 2 mass points are "
        "needed to interpolate.",
        grid->n_mass);

  grid->n_metallicity = 0;
  grid->metallicity = NULL;

  if (strcmp(grid->dimensionality, "M,Z") == 0) {
    grid->is_2d = 1;

    io_read_attribute(group_id, "nz", INT, &grid->n_metallicity);

    if (grid->n_metallicity < 2)
      error(
          "Data/Radiation's 'nz' attribute is %d; at least 2 metallicity "
          "points are needed to interpolate the log10(Z) axis.",
          grid->n_metallicity);

    grid->metallicity = (float *)malloc(sizeof(float) * grid->n_metallicity);
    if (grid->metallicity == NULL)
      error("Failed to allocate the RAD metallicity grid.");

    io_read_array_dataset(group_id, "Metallicity", FLOAT, grid->metallicity,
                          grid->n_metallicity);

    for (int i = 0; i < grid->n_metallicity; i++) {
      if (grid->metallicity[i] <= 0.f)
        error(
            "Data/Radiation's 'Metallicity' dataset entry %d is %.4g "
            "(<= 0): the metallicity axis is interpolated in log10(Z), "
            "which requires every entry to be strictly positive.",
            i, (double)grid->metallicity[i]);
    }

    char source[32];
    radiation_read_string_attribute(group_id, "source", source, sizeof(source));

    char luminosity_below[16], luminosity_above[16];
    char q_h_blackbody_below[16], q_h_blackbody_above[16];
    char q_h_parsec_below[16], q_h_parsec_above[16];
    radiation_read_string_attribute(group_id, "edge_policy_luminosity_below",
                                    luminosity_below, sizeof(luminosity_below));
    radiation_read_string_attribute(group_id, "edge_policy_luminosity_above",
                                    luminosity_above, sizeof(luminosity_above));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_blackbody_below",
                                    q_h_blackbody_below,
                                    sizeof(q_h_blackbody_below));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_blackbody_above",
                                    q_h_blackbody_above,
                                    sizeof(q_h_blackbody_above));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_parsec_below",
                                    q_h_parsec_below, sizeof(q_h_parsec_below));
    radiation_read_string_attribute(group_id, "edge_policy_q_h_parsec_above",
                                    q_h_parsec_above, sizeof(q_h_parsec_above));

    grid->edge_policy_luminosity = radiation_parse_edge_policy(
        luminosity_below, luminosity_above, "luminosity");
    const enum interpolate_boundary_condition edge_policy_q_h_blackbody =
        radiation_parse_edge_policy(q_h_blackbody_below, q_h_blackbody_above,
                                    "q_h_blackbody");
    const enum interpolate_boundary_condition edge_policy_q_h_parsec =
        radiation_parse_edge_policy(q_h_parsec_below, q_h_parsec_above,
                                    "q_h_parsec");
    grid->edge_policy_dot_e_excess = edge_policy_q_h_blackbody;

    if (strcmp(source, "parsec_blackbody") == 0) {
      grid->edge_policy_q_h = edge_policy_q_h_blackbody;
    } else if (strcmp(source, "parsec_qtable") == 0) {
      grid->edge_policy_q_h = edge_policy_q_h_parsec;
    } else {
      error(
          "Data/Radiation's 'source' attribute is '%s' (expected "
          "'parsec_blackbody' or 'parsec_qtable'): cannot tell which "
          "Q_H edge policy applies to the table's primary 'Q_H' dataset.",
          source);
    }
  } else if (strcmp(grid->dimensionality, "M") == 0) {
    grid->is_2d = 0;
    grid->edge_policy_luminosity = boundary_condition_error;
    grid->edge_policy_q_h = boundary_condition_error;
    grid->edge_policy_dot_e_excess = boundary_condition_error;
  } else {
    error(
        "Data/Radiation has an unrecognised 'dimensionality' attribute "
        "'%s' (expected 'M' or 'M,Z').",
        grid->dimensionality);
  }
}

/**
 * @brief Read one CGS-valued dataset from an open Data/Radiation group,
 * convert it to internal units and (for Q_H/DotEExcess)
 * #RADIATION_DOT_N_ION_TABLE_SCALING, narrow it to float, and (optionally)
 * also compute its log10, pychem-style, for the caller's raw (log-log)
 * interpolation table.
 *
 * Shared by radiation_read_luminosities_array(),
 * radiation_read_ionization_rate_array() and
 * radiation_read_mean_excess_photon_energy_array() to avoid tripling the
 * read/convert/guard boilerplate.
 *
 * Guards against float overflow the way stellar_evolution.c:676-695 does:
 * error() aborts (MPI_Abort/swift_abort, src/error.h) rather than capping
 * -- a units/scaling bug should stop the run, not silently corrupt the
 * physics. Also flags (debug-checks only) an implausible collapse to
 * exactly zero for a CGS input that was not itself zero: pychem bakes a
 * literal 0 into Q_H/DotEExcess below its own ionization threshold, so an
 * exact-zero result is only suspicious when the source value was nonzero.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to read.
 * @param count Number of elements to read (the group's "nm", or "nm" *
 * "nz" for a 2D table).
 * @param conversion_factor units_cgs_conversion_factor() for this
 * dataset's physical dimension; CGS values are divided by this to reach
 * internal units (SWIFT's convention).
 * @param extra_scaling Additional SWIFT-side-only divisor applied after
 * unit conversion (#RADIATION_DOT_N_ION_TABLE_SCALING for Q_H/DotEExcess,
 * 1 for Luminosity).
 * @param expected_units The dataset's own "units" attribute is asserted to
 * equal this string before any conversion happens (see
 * radiation_check_dataset_units()).
 * @param log_data_internal (output, optional) If not NULL, a caller-owned
 * float array of length @p count filled with log10 of the same
 * internal-unit value #RADIATION_LOG_FLOOR_CGS-floored on the CGS side
 * before conversion (see that macro's own doxygen) -- exactly pychem's
 * log-log convention, used to build a raw table's #interpolation_1d /
 * #interpolation_2d in log-value space instead of raw-value space. Left
 * untouched if NULL (the integrated-table caller has no use for it: see
 * radiation_build_tables()'s own doxygen for why the IMF-integrated table
 * stays in linear space).
 * @return Newly malloc'd float array of length count, in internal
 * (optionally rescaled) units, NOT logged -- this is the value the
 * IMF integration in radiation_build_tables() needs. Caller must free().
 */
static float *radiation_read_cgs_array(hid_t group_id, const char *dataset_name,
                                       hsize_t count, double conversion_factor,
                                       double extra_scaling,
                                       const char *expected_units,
                                       float *log_data_internal) {

  radiation_check_dataset_units(group_id, dataset_name, expected_units);

  double *data_cgs = (double *)malloc(sizeof(double) * count);
  if (data_cgs == NULL)
    error("Failed to allocate the RAD yields for %s.", dataset_name);

  io_read_array_dataset(group_id, dataset_name, DOUBLE, data_cgs, count);

  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the RAD yields for %s.", dataset_name);

  /* log10(internal value) = log10(cgs value) - log10(conversion_factor *
     extra_scaling); computed once here rather than per-entry below. */
  const double log_conversion = log10(conversion_factor) + log10(extra_scaling);

  for (hsize_t j = 0; j < count; j++) {
    const double value_internal =
        data_cgs[j] / conversion_factor / extra_scaling;

    if (fabs(value_internal) > (double)FLT_MAX) {
      error(
          "Radiation table '%s' entry %llu (%e cgs) converts to %e in "
          "internal units, exceeding FLT_MAX. This is a units/scaling bug; "
          "aborting rather than silently corrupting the physics.",
          dataset_name, (unsigned long long)j, data_cgs[j], value_internal);
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (data_cgs[j] != 0. && (float)value_internal == 0.0f) {
      message(
          "WARNING: radiation table '%s' entry %llu (%e cgs, nonzero) "
          "collapsed to exactly 0 in internal units after conversion -- "
          "check RADIATION_DOT_N_ION_TABLE_SCALING and the unit system.",
          dataset_name, (unsigned long long)j, data_cgs[j]);
    }
#endif

    data[j] = (float)value_internal;

    if (log_data_internal != NULL) {
      const double floored_cgs = max(data_cgs[j], RADIATION_LOG_FLOOR_CGS);
      const double log_value_internal = log10(floored_cgs) - log_conversion;
      log_data_internal[j] = (float)log_value_internal;

#ifdef SWIFT_DEBUG_CHECKS
      /* Self-check the log10/exp10 round-trip against the independently
         computed linear value above, away from the floor (where the two
         are expected to diverge by construction): a sign error or a
         dropped extra_scaling term in log_conversion would silently make
         every raw getter wrong by many orders of magnitude while still
         running to completion, so this is checked at load time on every
         run rather than trusted from inspection alone. */
      if (data_cgs[j] > RADIATION_LOG_FLOOR_CGS * 1e10) {
        const double round_trip = exp10(log_value_internal);
        const double rel_diff =
            fabs(round_trip - value_internal) / fabs(value_internal);
        if (rel_diff > 1e-4) {
          error(
              "Radiation table '%s' entry %llu: log10/exp10 round-trip "
              "mismatch (internal=%e, round-trip=%e, rel_diff=%e). This "
              "indicates a bug in the log-log conversion, not the data.",
              dataset_name, (unsigned long long)j, value_internal, round_trip,
              rel_diff);
        }
      }
#endif
    }
  }

  free(data_cgs);
  return data;
}

/**
 * @brief Build the raw (and, for a 1D table, IMF-integrated) interpolation
 * table for one Data/Radiation quantity, dispatching on the table's
 * dimensionality.
 *
 * The raw table (1D or 2D) is built in log10(value) space, floored and
 * exponentiated pychem-style (see #RADIATION_LOG_FLOOR_CGS and
 * radiation_read_cgs_array()): #interpolate_1d_init()/#interpolate_2d_init()
 * are otherwise-unmodified generic linear interpolators, so feeding them
 * already-logged data makes their existing linear interpolation a log-log
 * interpolation for free. Every raw getter (radiation_get_*_from_raw())
 * must exponentiate the result back; see their own doxygen.
 *
 * The IMF-integrated table (1D only) is deliberately left in linear
 * (un-logged) value space, unlike the raw table above: it is not the same
 * kind of quantity pychem's log-log scheme governs. Pychem's convention
 * interpolates a single tabulated property value as a function of mass;
 * #initial_mass_function_integrate() instead turns `data` in place into a
 * cumulative IMF-weighted integral (a running sum starting at exactly 0 at
 * the IMF's own mass_min) -- a SWIFT-side population-synthesis technique
 * with no pychem analogue to match (pychem never computes or interpolates
 * such a quantity). Integrating log10(value) instead of value itself would
 * not even compute the right integral, so this ordering (raw table built
 * from the logged copy, then the *un-logged* `data` handed to
 * initial_mass_function_integrate()) is required, not just a style choice.
 *
 * The 2D ("M,Z") branch only builds the raw table: no caller integrates a
 * 2D table over the IMF yet (see #radiation's own doxygen), so
 * @p integrated_1d is left untouched -- radiation_read_data() zeroes it
 * beforehand so radiation_clean() stays safe either way. It also
 * approximates the metallicity axis as log-uniformly spaced from the
 * "Metallicity" dataset's first/last values: pychem does not guarantee
 * this (it is the curated set of PARSEC metallicities actually collapsed
 * into the table, not a synthetic grid), but interpolate_2d_init()
 * requires a uniform grid. Acceptable only because this 2D path has zero
 * test coverage until the PARSEC validation track (plan Phase 6)
 * exercises it.
 *
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param dataset_name Name of the dataset to read.
 * @param grid The group's own grid metadata (see
 * radiation_read_grid_metadata()).
 * @param sm The #stellar_model (for the output mass-grid bounds).
 * @param interpolation_size_mass Number of points in the mass
 * interpolation output grid.
 * @param interpolation_size_metallicity Number of points in the
 * metallicity interpolation output grid (2D tables only).
 * @param conversion_factor See radiation_read_cgs_array().
 * @param extra_scaling See radiation_read_cgs_array().
 * @param expected_units See radiation_read_cgs_array().
 * @param raw_1d (output) Raw 1D interpolation table (1D tables only),
 * holding log10(value in internal units), pychem-floored.
 * @param integrated_1d (output) IMF-integrated 1D interpolation table (1D
 * tables only), holding the linear (un-logged) cumulative value -- see
 * this function's own doxygen for why.
 * @param raw_2d (output) Raw 2D interpolation table (2D tables only),
 * holding log10(value in internal units), pychem-floored.
 * @param boundary_condition_mass Mass-axis #interpolate_boundary_condition
 * for @p dataset_name (2D tables only; ignored for a 1D table, which always
 * clamps -- see radiation_read_grid_metadata()'s edge_policy_* fields for
 * where this comes from). The metallicity axis always clamps
 * (boundary_condition_const), matching pychem's "clamp to nearest grid Z,
 * never extrapolated" convention.
 */
static void radiation_build_tables(
    hid_t group_id, const char *dataset_name,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    int interpolation_size_mass, int interpolation_size_metallicity,
    double conversion_factor, double extra_scaling, const char *expected_units,
    struct interpolation_1d *raw_1d, struct interpolation_1d *integrated_1d,
    struct interpolation_2d *raw_2d,
    enum interpolate_boundary_condition boundary_condition_mass) {

  const float log_mass_min_out = log10f(sm->imf.mass_min);
  const float log_mass_max_out = log10f(sm->imf.mass_max);

  if (grid->is_2d) {

    const hsize_t count = (hsize_t)grid->n_mass * (hsize_t)grid->n_metallicity;
    float *log_data = (float *)malloc(sizeof(float) * count);
    if (log_data == NULL)
      error("Failed to allocate the RAD 2D log-value yields for %s.",
            dataset_name);
    float *data = radiation_read_cgs_array(group_id, dataset_name, count,
                                           conversion_factor, extra_scaling,
                                           expected_units, log_data);

    const float log_z_min = log10f(grid->metallicity[0]);
    const float log_z_max = log10f(grid->metallicity[grid->n_metallicity - 1]);
    const float log_z_step =
        grid->n_metallicity > 1
            ? (log_z_max - log_z_min) / (grid->n_metallicity - 1)
            : 0.f;

    /* interpolate_2d_init() takes a double source array (its internal
       storage is float; see interpolation.h); re-widen the already
       guarded/narrowed/logged float data rather than duplicating the guard
       for a double codepath. */
    double *log_data_double = (double *)malloc(sizeof(double) * count);
    if (log_data_double == NULL)
      error("Failed to allocate the RAD 2D log-value yields for %s.",
            dataset_name);
    for (hsize_t i = 0; i < count; i++)
      log_data_double[i] = (double)log_data[i];

    interpolate_2d_init(raw_2d, log_z_min, log_z_max,
                        interpolation_size_metallicity, log_mass_min_out,
                        log_mass_max_out, interpolation_size_mass, log_z_min,
                        grid->log_mass_min, log_z_step, grid->mass_step,
                        grid->n_metallicity, grid->n_mass, log_data_double,
                        boundary_condition_const, boundary_condition_mass);

    free(log_data_double);
    free(log_data);
    free(data);
    return;
  }

  float *log_data = (float *)malloc(sizeof(float) * grid->n_mass);
  if (log_data == NULL)
    error("Failed to allocate the RAD log-value yields for %s.", dataset_name);
  float *data = radiation_read_cgs_array(
      group_id, dataset_name, (hsize_t)grid->n_mass, conversion_factor,
      extra_scaling, expected_units, log_data);

  interpolate_1d_init(raw_1d, log_mass_min_out, log_mass_max_out,
                      interpolation_size_mass, grid->log_mass_min,
                      grid->mass_step, grid->n_mass, log_data,
                      boundary_condition_const);

  /* initial_mass_function_integrate() mutates data (the LINEAR, un-logged
     copy -- see this function's own doxygen) in place into its cumulative
     IMF integral; must run after the raw table above is built from the
     un-integrated values. */
  initial_mass_function_integrate(&sm->imf, data, grid->n_mass,
                                  grid->log_mass_min, grid->mass_step);

  interpolate_1d_init(integrated_1d, log_mass_min_out, log_mass_max_out,
                      interpolation_size_mass, grid->log_mass_min,
                      grid->mass_step, grid->n_mass, data,
                      boundary_condition_const);

  free(data);
  free(log_data);
}

/**
 * @brief Read an array of luminosities data from the table.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_luminosities_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(
      group_id, "Luminosity", grid, sm, rad->interpolation_size,
      rad->interpolation_size_metallicity,
      units_cgs_conversion_factor(us, UNIT_CONV_POWER), 1., "erg/s",
      &rad->raw.luminosities, &rad->integrated.luminosities,
      &rad->raw.luminosities_2d, grid->edge_policy_luminosity);
}

/**
 * @brief Read an array of ionizing emission rates data from the table.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_ionization_rate_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(
      group_id, "Q_H", grid, sm, rad->interpolation_size,
      rad->interpolation_size_metallicity,
      units_cgs_conversion_factor(us, UNIT_CONV_PHOTONS_PER_TIME),
      RADIATION_DOT_N_ION_TABLE_SCALING, "1/s", &rad->raw.dot_N_ion,
      &rad->integrated.dot_N_ion, &rad->raw.dot_N_ion_2d,
      grid->edge_policy_q_h);
}

/**
 * @brief Read an array of excess-photon-energy emission rate data from the
 * table: DotEExcess(m) = Q_H(m) * MeanExcessPhotonEnergyHI(m).
 *
 * Converted with the SAME (rate-only) #UNIT_CONV_PHOTONS_PER_TIME factor
 * used for Q_H, not #UNIT_CONV_POWER -- deliberately, even though the file
 * stores DotEExcess in erg/s (a power): dividing only the rate part by
 * #RADIATION_DOT_N_ION_TABLE_SCALING and the unit conversion, while
 * leaving the "erg" part of the product in cgs, reproduces the mixed-unit
 * convention #radiation_get_mean_excess_photon_energy_HI_from_integral (an
 * existing, unmodified function/call site) already relies on: both raw
 * tables share the same rate-only scaling, so it cancels exactly in that
 * ratio, and the result comes out in cgs erg -- matching
 * feedback_struct.h's documented cgs-erg convention for
 * mean_excess_photon_energy_HI, without needing a unit_system argument on
 * the getters. Using #UNIT_CONV_POWER here instead would leave that ratio
 * in internal energy units, silently changing the existing getter's
 * output. The units check below still asserts "erg/s" (not the rate-only
 * factor's implied "1/s"): it verifies the table's stored unit, which the
 * deliberate mismatch above depends on staying exactly "erg/s" for the
 * cancellation to hold.
 *
 * @param rad The #radiation model.
 * @param group_id Open HDF5 "Data/Radiation" group id.
 * @param grid The group's own grid metadata.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void radiation_read_mean_excess_photon_energy_array(
    struct radiation *rad, hid_t group_id,
    const struct radiation_grid_metadata *grid, const struct stellar_model *sm,
    const struct unit_system *us) {

  radiation_build_tables(
      group_id, "DotEExcess", grid, sm, rad->interpolation_size,
      rad->interpolation_size_metallicity,
      units_cgs_conversion_factor(us, UNIT_CONV_PHOTONS_PER_TIME),
      RADIATION_DOT_N_ION_TABLE_SCALING, "erg/s", &rad->raw.dot_E_excess,
      &rad->integrated.dot_E_excess, &rad->raw.dot_E_excess_2d,
      grid->edge_policy_dot_e_excess);
}

/**
 * @brief Open the "Data/Radiation" group of a yields table, with an
 * actionable error naming the file and the regeneration step if the group
 * is missing.
 *
 * Every yields table generated before this migration (the ones currently
 * committed/fetched for the 8 radiation examples included) lacks this
 * group, and #h5_open_group's own generic "unable to open group" message
 * gives no hint that regenerating the table, not fixing a typo, is the
 * actual fix -- mirrors the actionable-message convention this file
 * already uses for GEARFeedback:HII_angular_nside (radiation_init(),
 * above).
 *
 * @param filename The yields table filename (sm->yields_table or
 * sm->yields_table for the first-stars model -- same check either way).
 * @param file_id (output) The opened HDF5 file id.
 * @param group_id (output) The opened "Data/Radiation" group id.
 */
static void radiation_open_data_group(const char *filename, hid_t *file_id,
                                      hid_t *group_id) {

  *file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (*file_id < 0) error("unable to open file %s.\n", filename);

  const htri_t exists = H5Lexists(*file_id, "Data/Radiation", H5P_DEFAULT);
  if (exists <= 0) {
    error(
        "'%s' has no 'Data/Radiation' group. This yields table was "
        "generated before the radiation-table migration and needs "
        "regenerating: run pychem's pychem_generate_hdf5_parameters on "
        "this table's own chimieparam file, then point "
        "GEARFeedback:yields_table (or yields_table_first_stars, for the "
        "PopIII model) at the regenerated file.",
        filename);
  }

  *group_id = H5Gopen(*file_id, "Data/Radiation", H5P_DEFAULT);
  if (*group_id < 0)
    error("unable to open group 'Data/Radiation' in %s.\n", filename);
}

/**
 * @brief Read the RAD yields from the table.
 *
 * The tables are in internal units at the end of this function.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void radiation_read_data(struct radiation *rad, struct swift_params *params,
                         const struct stellar_model *sm,
                         const struct unit_system *us,
                         const struct phys_const *phys_const,
                         const int restart) {

  if (!restart) {
    rad->interpolation_size = parser_get_opt_param_int(
        params, "GEARRadiation:interpolation_size_mass", 500);
    rad->interpolation_size_metallicity = parser_get_opt_param_int(
        params, "GEARRadiation:interpolation_size_metallicity", 110);
    if (rad->interpolation_size < 2) {
      error(
          "GEARRadiation:interpolation_size_mass must be >= 2; got "
          "%d.",
          rad->interpolation_size);
    }
    if (rad->interpolation_size_metallicity < 2) {
      error(
          "GEARRadiation:interpolation_size_metallicity must be >= "
          "2; got %d.",
          rad->interpolation_size_metallicity);
    }
  }

  /* radiation_zero_pointers() below also clears interpolation_size(_
     metallicity) and n_HII_pixels, so callers with radiation disabled
     don't inherit uninitialized garbage in them -- but on this call path
     those three either were just set above (!restart) or hold the values
     radiation_restore() flat-restored moments ago (restart), and
     radiation_build_tables() below needs them either way. Round-trip them
     around the call. */
  const int interpolation_size_before = rad->interpolation_size;
  const int interpolation_size_metallicity_before =
      rad->interpolation_size_metallicity;
  const int n_HII_pixels_before = rad->n_HII_pixels;

  /* Zero every table up front: radiation_build_tables() only populates the
     _1d or _2d variant matching this table's dimensionality, and only the
     2D path's raw tables at that (see its own doxygen) -- the rest must be
     safe no-ops for radiation_clean()'s interpolate_1d_free()/
     interpolate_2d_free() calls regardless of which branch ran. */
  radiation_zero_pointers(rad);

  rad->interpolation_size = interpolation_size_before;
  rad->interpolation_size_metallicity = interpolation_size_metallicity_before;
  rad->n_HII_pixels = n_HII_pixels_before;

  hid_t file_id, group_id;
  radiation_open_data_group(sm->yields_table, &file_id, &group_id);

  struct radiation_grid_metadata grid;
  radiation_read_grid_metadata(group_id, &grid);
  rad->is_2d = grid.is_2d;

  /* Table-coverage check: the raw table's boundary condition is
     boundary_condition_const (radiation_build_tables()), so a star outside
     the table's own native mass grid silently gets the nearest edge's
     value instead of its own -- e.g. pychem's real PopIII table has a
     native floor of 13 Msun, well above a typical IMF's own mass_min.
     Fail loudly instead, matching this file's existing convention
     (the FLT_MAX guard, the HII_angular_nside checks above).

     Compared with a half-grid-cell tolerance, not exact equality: both
     sides are independently accumulated (log_mass_min_imf/log_mass_max_imf
     from log10f() of the IMF's own bounds, log_mass_max_table from the
     grid's own log_mass_min/mass_step/n_mass), so a star whose mass range
     was generated to sit exactly on the table's own edge can differ from
     it by a few ULP of float rounding. A tolerance-free comparison rejects
     the overwhelming majority of otherwise-valid (mass_min, mass_max, nm)
     combinations for no physical reason; a genuine shortfall (e.g. the
     documented 13 Msun PopIII floor case) is orders of magnitude outside
     this tolerance and still aborts. */
  const float log_mass_min_imf = log10f(sm->imf.mass_min);
  const float log_mass_max_imf = log10f(sm->imf.mass_max);
  const double log_mass_max_table =
      (double)grid.log_mass_min + (grid.n_mass - 1) * (double)grid.mass_step;
  const double tol = 0.5 * (double)grid.mass_step;
  if ((double)log_mass_min_imf < (double)grid.log_mass_min - tol ||
      (double)log_mass_max_imf > log_mass_max_table + tol) {
    error(
        "'%s': Data/Radiation's native mass grid [%.4g, %.4g] Msun does not "
        "cover the IMF's mass range [%.4g, %.4g] Msun. A star outside the "
        "table's own range would silently receive the nearest edge's "
        "photon budget (boundary_condition_const) instead of its own "
        "value. Regenerate the table over (at least) the IMF's mass range "
        "with pychem, or adjust the IMF's own mass_min/mass_max to fit "
        "inside the table.",
        sm->yields_table, (double)exp10(grid.log_mass_min),
        (double)exp10(log_mass_max_table), (double)sm->imf.mass_min,
        (double)sm->imf.mass_max);
  }

  /* Fail at load time, not at the first star's feedback computation: every
     radiation_get_*_from_raw()/_from_integral() getter aborts on a 2D
     table (see radiation_check_is_1d()), since no caller passes a
     metallicity yet. */
  if (rad->is_2d && engine_rank == 0) {
    message(
        "WARNING: '%s' is a mass x metallicity (2D) radiation table; no "
        "individual-star or population feedback getter supports one yet -- "
        "any star reaching feedback will abort the run.",
        sm->yields_table);
  }

  /* Read the luminosities */
  radiation_read_luminosities_array(rad, group_id, &grid, sm, us);

  /* Read the ionization emission rates */
  radiation_read_ionization_rate_array(rad, group_id, &grid, sm, us);

  /* Read the excess-photon-energy emission rates */
  radiation_read_mean_excess_photon_energy_array(rad, group_id, &grid, sm, us);

  free(grid.metallicity);
  h5_close_group(file_id, group_id);

  /* The tables above are now valid: mark this #radiation active so callers
     use them instead of skipping to their zeroed defaults. */
  rad->is_active = 1;
};
