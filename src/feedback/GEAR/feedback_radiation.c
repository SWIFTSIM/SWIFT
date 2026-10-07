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
/* Include header */
#include "feedback_radiation.h"

#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "feedback_common.h"
#include "hydro_properties.h"
#include "minmax.h"
#include "part.h"
#include "radiation.h"
#include "radiation_isrf.h"
#include "timeline.h"
#include "timestep_sync_part.h"
#include "units.h"

#include <float.h>

/**
 * @brief Determines whether a gas #part can be ionized.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param e The #engine.
 * @return Is the particle ionized?
 */
__attribute__((always_inline)) INLINE char feedback_part_can_be_ionized(
    const struct part *p, const struct xpart *xp, const struct engine *e) {

  const struct phys_const *phys_const = e->physical_constants;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct feedback_props *feedback_props = e->feedback_props;
  const struct unit_system *us = e->internal_units;
  const struct cosmology *cosmo = e->cosmology;
  const struct cooling_function_data *cooling = e->cooling_func;

  /* Is T > 10^4 K ? */
  const float T = cooling_get_temperature(phys_const, hydro_props, us, cosmo,
                                          cooling, p, xp);
  const float ten_to_four_kelvin =
      1e4 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* The 1.01 factor is here for safety margin and numerical stability */
  const char is_cold = (T <= 1.01 * ten_to_four_kelvin);

  /* Density threshold criterion */
  const float rho = hydro_get_physical_density(p, cosmo);
  const float rho_threshold = feedback_props->HII_min_density;
  const char is_dense = rho >= rho_threshold;

  /* Can the particle be ionized? */
  return (is_cold && is_dense && !radiation_is_part_tagged_as_ionized(p, xp));
}

/**
 * @brief Photoionization rate coefficient this star delivers at a gas
 * particle's location, frozen at tag time.
 *
 * The cooling task cannot recompute it later (it has no neighbour search),
 * so it is stored on the particle. The intermediate photon flux would
 * overflow float32 in this unit system, so only the final coefficient,
 * computed in double up to that point, is returned. Zero unless
 * GEARFeedback:HII_couple_ionization_rate is on.
 *
 * @param si The #spart (star) providing photons.
 * @param r2 Squared comoving distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param us Internal unit system.
 * @param cosmo The current cosmological model.
 * @param cooling Cooling function data.
 * @return Photoionization rate coefficient Gamma_HI, internal 1/time, or 0
 * if GEARFeedback:HII_couple_ionization_rate is off.
 */
__attribute__((always_inline)) INLINE static float
feedback_hii_photoionization_rate_HI(
    const struct spart *restrict si, float r2, int pixel,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling) {

  if (!cooling->HII_couple_ionization_rate) return 0.f;

  const float Omega_pixel =
      4.0f * (float)M_PI / si->feedback_data.radiation.n_HII_pixels;

  /* r2 is comoving, so the physical 1/r^2 dilution carries an a^-2. Floored at
     a fraction of the star's own smoothing length: a gas particle sitting on
     top of the star would otherwise divide by zero and hand Grackle an
     infinite heating rate. */
  const float r2_min = 1e-6f * si->h * si->h;
  const double r2_phys = max(r2, r2_min) * cosmo->a * cosmo->a;
  const double ionizing_flux_HI =
      feedback_get_star_ionization_rate(si, pixel) / (Omega_pixel * r2_phys);
  return (float)radiation_get_photoionization_rate_coefficient_from_flux_HI(
      us, ionizing_flux_HI);
}

/**
 * @brief Absolute time until which an ionization tag stamped now stays valid.
 *
 * Sized on the interval the pass just covered, so a tag outlives the gap to
 * its star's next rebuild and cooling cannot expire it first.
 *
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 * @return Absolute simulation time until which the tag stays valid.
 */
__attribute__((always_inline)) INLINE static double feedback_hii_tag_end_time(
    const double time, const double dt_back) {
  return time + RADIATION_TAG_LIFETIME_INTERVALS * dt_back;
}

/**
 * @brief Tag a gas particle as ionized by this star and charge its photons.
 *
 * Shared by the three acceptance branches of #feedback_iact_HII_ionization.
 * No atomics: task_type_stars_hii_ionization_feedback's cell locking
 * (src/task.c) already serializes every writer of these fields.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being ionized.
 * @param xpj The #xpart (gas) being ionized.
 * @param r2 Squared distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param us Internal unit system.
 * @param cosmo The current cosmological model.
 * @param cooling Cooling function data.
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 * @param cost Ionizing photon count to charge for this particle.
 */
__attribute__((always_inline)) INLINE static void feedback_hii_claim_part(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const double time,
    const double dt_back, const double cost) {

  if (radiation_is_part_tagged_as_ionized(pj, xpj)) return;

  radiation_tag_part_as_ionized(
      pj, xpj, si->id, feedback_hii_tag_end_time(time, dt_back),
      si->feedback_data.radiation.mean_excess_photon_energy_HI,
      feedback_hii_photoionization_rate_HI(si, r2, pixel, us, cosmo, cooling));
  timestep_sync_part(pj);

  radiation_consume_ionizing_photons(si, pixel, cost);

  /* Update HII region properties */
  si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);
  si->h_hii = max(si->h_hii, sqrtf(r2) * kernel_gamma_inv);
}

/**
 * @brief Charge an already-ionized gas particle's recombination losses and
 * renew its tag.
 *
 * Recombinations must be replaced continuously to hold gas ionized, so a
 * particle already held ionized costs photons every pass: charging only
 * newly-claimed ones is what let the ionized volume grow without bound as the
 * rebuild cadence was refined. There is no one-off N_H term here: those
 * electrons are already stripped.
 *
 * Renewal is budget-gated: once a pixel is spent, its remaining gas is not
 * renewed and cooling expires those tags at end_time. The caller
 * (runner_dosub_stars_hii_ionization_feedback) invokes this in ascending
 * (r2, id) order over the star's whole held region, so a shortfall lapses
 * the outermost shell (largest r2, the recession front) first, rather than
 * a set determined by cell-traversal order.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being maintained.
 * @param xpj The #xpart (gas) being maintained.
 * @param r2 Squared distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param phys_const Physics constants.
 * @param hydro_props Hydrodynamics properties.
 * @param us Internal unit system.
 * @param cosmo Cosmology.
 * @param cooling Cooling function data.
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 * @return Photons charged for this particle, 0 if the budget was exhausted.
 */
__attribute__((always_inline)) INLINE double
feedback_iact_HII_maintain_ionized_part(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const double time,
    const double dt_back) {

  /* Counted whether or not its photons can be afforded: this field is the
     mass the star HOLDS ionized, matching h_hii's extent (see stars_io.h). */
  si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);

  if (feedback_get_star_ionization_budget(si, pixel) <= 0.0) return 0.0;

  const double cost =
      dt_back * radiation_get_part_rate_to_fully_ionize(
                    phys_const, hydro_props, us, cosmo, cooling, pj, xpj);
  radiation_consume_ionizing_photons(si, pixel, cost);

  /* Renew in place: the tag's other fields would otherwise keep values
     frozen at first claim, growing staler every pass the particle survives.
     Deliberately no timestep_sync_part() here, unlike a fresh claim: waking
     every held particle every pass would drag the whole region down to the
     shortest time bin. Expiry is therefore checked on the gas particle's own
     next cooling task, so recession can lag by up to one of its timesteps. */
  radiation_tag_part_as_ionized(
      pj, xpj, si->id, feedback_hii_tag_end_time(time, dt_back),
      si->feedback_data.radiation.mean_excess_photon_energy_HI,
      feedback_hii_photoionization_rate_HI(si, r2, pixel, us, cosmo, cooling));

  return cost;
}

/**
 * @brief Perform the ionization of a gas particle by a star.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being ionized.
 * @param xpj The #xpart (gas) being ionized.
 * @param r2 Squared distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param phys_const Physics constants.
 * @param hydro_props Hydrodynamics properties.
 * @param us Internal unit system.
 * @param cosmo Cosmology.
 * @param cooling Cooling function data.
 * @param feedback_props The #feedback_props.
 * @param ti_begin Integer time at the start of the step (for RNG).
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 */
__attribute__((always_inline)) INLINE void feedback_iact_HII_ionization(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling,
    const struct feedback_props *feedback_props, const integertime_t ti_begin,
    const double time, const double dt_back) {

  const int deterministic_boundary =
      feedback_props->HII_deterministic_boundary_ionization;

  /* If already ionized (by another thread, or by an earlier pass), just move
     on: recombination losses for gas already held ionized are charged by
     feedback_iact_HII_maintain_ionized_part, not here. */
  if (radiation_is_part_tagged_as_ionized(pj, xpj)) return;

  /* Photons this candidate costs: a one-off payment to strip its remaining
     NEUTRAL hydrogen (not its total hydrogen content: a particle whose
     tag lapsed on a marginal budget shortfall, or one pre-ionized by a UV
     background, is already partway or fully stripped, and re-paying full
     N_H on reclaim would be a cadence-coupled photon sink), plus a
     maintenance reserve sized by the *elapsed* interval (the next
     interval's recombinations are charged by
     feedback_iact_HII_maintain_ionized_part, not here). That reserve is what
     makes region growth implicit (dS/dt = (Q-S)/t_rec integrates as
     dS = (Q-S)*dt/(t_rec+dt)), and so unconditionally stable at the
     dt >> t_rec the default HII_rebuild_time_Myr produces. Dropping it would
     give explicit Euler, which overshoots and then churns. */
  const double N_HI = radiation_get_part_number_neutral_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, pj, xpj);
  const double Delta_dot_N_ion_maintenance_rate =
      radiation_get_part_rate_to_fully_ionize(phys_const, hydro_props, us,
                                              cosmo, cooling, pj, xpj);
  const double cost = N_HI + dt_back * Delta_dot_N_ion_maintenance_rate;

  const double budget = feedback_get_star_ionization_budget(si, pixel);

  /* Case 1: Ionization is guaranteed */
  if (cost <= budget) {
    feedback_hii_claim_part(si, pj, xpj, r2, pixel, us, cosmo, cooling, time,
                            dt_back, cost);
  } else if (deterministic_boundary) {
    /* Deterministic mode: always ionize the boundary particle, letting the
       pixel's budget go slightly negative (bounded by one particle's
       ionization cost). */
    feedback_hii_claim_part(si, pj, xpj, r2, pixel, us, cosmo, cooling, time,
                            dt_back, cost);
  } else {
    /* Probabilistic mode: a weighted coin flip decides whether to fully
       ionize pj. On a win, the full cost is consumed (more than what was
       left, going slightly negative); on a loss, nothing is consumed, so
       the budget survives intact to try the next, farther candidate. This
       keeps the expected photon consumption equal to what's actually
       available: proba*cost + (1-proba)*0 = budget. */
    const float proba = budget / cost;
    /* Keyed on the star and the pixel, deliberately *not* on pj: this must be
       one trial per pixel per pass for the identity below to hold. The same
       number is drawn for every candidate offered to this pixel, so a loss
       rejects the whole remaining shell, which is exactly the (1 - proba)
       branch. Rolling per candidate instead would give each of the up-to
       HII_max_retry_full_buffer x max_ngbs candidates an independent shot in a
       loop that stops at the first win, so a pass would claim one particle
       almost surely and spend `cost` rather than `budget`. */
    const float random_number = random_unit_interval_part_ID_and_index(
        si->id, pixel, ti_begin, random_number_HII_regions);

    if (random_number <= proba) {
      /* We won the roll! Claim the particle. */
      feedback_hii_claim_part(si, pj, xpj, r2, pixel, us, cosmo, cooling, time,
                              dt_back, cost);
    }
    /* Lost the roll: consume nothing, budget carries over. */
  } /* End of probability handling */
}

/**
 * @brief Is this gas particle currently tagged as HII-ionized?
 *
 * Thin dispatch wrapper so callers outside this feedback model (e.g.
 * star_formation/GEAR, sink/GEAR, both of which are selectable
 * independently of the feedback model) can query ionization state without
 * depending on this model being the one actually compiled in: every
 * feedback model provides this function, matching #radiation_is_part_
 * tagged_as_ionized() here for GEAR and unconditionally returning false
 * everywhere else.
 *
 * @param p The #part to query.
 * @param xp The #part's extended data.
 * @return 1 if the particle is tagged as HII-ionized, 0 otherwise.
 */
char feedback_is_part_tagged_as_ionized(const struct part *p,
                                        const struct xpart *xp) {
  return radiation_is_part_tagged_as_ionized(p, xp);
}

/**
 * @brief Id of the star that tagged this gas particle as HII-ionized.
 *
 * Thin dispatch wrapper, same reasoning as
 * #feedback_is_part_tagged_as_ionized. Only meaningful while that function
 * returns true.
 *
 * @param p The #part to query.
 * @param xp The #part's extended data.
 * @return The id of the star that tagged this particle as HII-ionized.
 */
long long feedback_get_part_ionized_star_id(const struct part *p,
                                            const struct xpart *xp) {
  return radiation_get_part_ionized_star_id(p, xp);
}

/**
 * @brief Local specific PE-band radiation field, see
 * #feedback_part_data.isrf_moment[ISRF_MOMENT_PE].u. Thin dispatch wrapper,
 * same reasoning as #feedback_is_part_tagged_as_ionized: every feedback model
 * provides this function, returning 0 everywhere except here for GEAR.
 *
 * @param p The #part to query.
 * @return Local specific PE-band radiation field.
 */
double feedback_get_part_u_PE(const struct part *p) {
  return p->feedback_data.isrf_moment[ISRF_MOMENT_PE].u;
}

/**
 * @brief Local specific Lyman-Werner-band radiation field, see
 * #feedback_get_part_u_PE.
 *
 * @param p The #part to query.
 * @return Local specific Lyman-Werner-band radiation field.
 */
double feedback_get_part_u_LW(const struct part *p) {
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW].u;
}

/**
 * @brief Local Lyman-Werner-band photon-number moment, energy-equivalent
 * at a fixed reference photon energy, see #feedback_get_part_u_PE.
 *
 * @param p The #part to query.
 * @return Local Lyman-Werner-band photon-number moment.
 */
double feedback_get_part_u_LW_PHOTON(const struct part *p) {
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u;
}

/**
 * @brief Negativity-triggered artificial-dissipation coefficient, see
 * #feedback_part_data.isrf_operator[ISRF_OPERATOR_PE].dissipation_alpha_trigger
 * and
 * #feedback_part_data.isrf_operator[ISRF_OPERATOR_PE].dissipation_alpha_floor.
 * Thin dispatch wrapper, same reasoning as #feedback_get_part_u_PE.
 *
 * Per-particle SUMMARY for I/O only: the coefficient the force loop uses is
 * the per-pair `alpha_ij = max(trigger_i, trigger_j, floor_i, floor_j)`,
 * which has no single-particle representation. This returns the value the
 * particle would contribute against an identical partner.
 *
 * @param p The #part to query.
 * @return Negativity-triggered artificial-dissipation coefficient, PE band.
 */
float feedback_get_part_dissipation_alpha_PE(const struct part *p) {
  return max(
      p->feedback_data.isrf_operator[ISRF_OPERATOR_PE]
          .dissipation_alpha_trigger,
      p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].dissipation_alpha_floor);
}

/**
 * @brief See #feedback_get_part_dissipation_alpha_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @return Negativity-triggered artificial-dissipation coefficient,
 * Lyman-Werner band.
 */
float feedback_get_part_dissipation_alpha_LW(const struct part *p) {
  return max(
      p->feedback_data.isrf_operator[ISRF_OPERATOR_LW]
          .dissipation_alpha_trigger,
      p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].dissipation_alpha_floor);
}

/**
 * @brief See #feedback_get_part_dissipation_alpha_PE, Lyman-Werner-band
 * photon-number moment. Decorative: this moment shares #ISRF_OPERATOR_LW
 * with #ISRF_MOMENT_LW, so it reads the identical operator state as
 * #feedback_get_part_dissipation_alpha_LW.
 *
 * @param p The #part to query.
 * @return Negativity-triggered artificial-dissipation coefficient,
 * Lyman-Werner-band photon-number moment.
 */
float feedback_get_part_dissipation_alpha_LW_PHOTON(const struct part *p) {
  return max(
      p->feedback_data.isrf_operator[ISRF_OPERATOR_LW]
          .dissipation_alpha_trigger,
      p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].dissipation_alpha_floor);
}

/**
 * @brief `(1/rho) div(rho F)` accumulator, see
 * #feedback_part_data.isrf_moment[ISRF_MOMENT_PE].div_specific_flux. Thin
 * dispatch wrapper, same reasoning as #feedback_get_part_u_PE.
 *
 * @param p The #part to query.
 * @return `(1/rho) div(rho F)` accumulator, PE band.
 */
float feedback_get_part_div_specific_flux_PE(const struct part *p) {
  return p->feedback_data.isrf_moment[ISRF_MOMENT_PE].div_specific_flux;
}

/**
 * @brief See #feedback_get_part_div_specific_flux_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @return `(1/rho) div(rho F)` accumulator, Lyman-Werner band.
 */
float feedback_get_part_div_specific_flux_LW(const struct part *p) {
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW].div_specific_flux;
}

/**
 * @brief See #feedback_get_part_div_specific_flux_PE, Lyman-Werner-band
 * photon-number moment.
 *
 * @param p The #part to query.
 * @return `(1/rho) div(rho F)` accumulator, Lyman-Werner-band photon-number
 * moment.
 */
float feedback_get_part_div_specific_flux_LW_PHOTON(const struct part *p) {
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].div_specific_flux;
}

/**
 * @brief Tracked specific flux moment, see
 * #feedback_part_data.isrf_moment[ISRF_MOMENT_PE].specific_flux, returned as
 * the true physical flux `F_true`. The stored field is the reduced flux
 * `Ft = F_true/c_hyp` (radiation_propagation_iact.h's file header), so this
 * getter rescales it back (`F_true = c_hyp*Ft`) before returning: the io
 * field's meaning (PESpecificFluxes, feedback_io.h) and every downstream
 * consumer (the ISRFHyperbolicPropagation check scripts included) do not
 * depend on this internal representation. A multiply, never a division: no
 * zero-denominator hazard.
 *
 * @param p The #part to query.
 * @param ret (return) The three components.
 */
void feedback_get_part_specific_flux_PE(const struct part *p, float *ret) {
  const float rescale = p->feedback_data.c_hyp;
  ret[0] =
      rescale * p->feedback_data.isrf_moment[ISRF_MOMENT_PE].specific_flux[0];
  ret[1] =
      rescale * p->feedback_data.isrf_moment[ISRF_MOMENT_PE].specific_flux[1];
  ret[2] =
      rescale * p->feedback_data.isrf_moment[ISRF_MOMENT_PE].specific_flux[2];
}

/**
 * @brief See #feedback_get_part_specific_flux_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @param ret (return) The three components.
 */
void feedback_get_part_specific_flux_LW(const struct part *p, float *ret) {
  const float rescale = p->feedback_data.c_hyp;
  ret[0] =
      rescale * p->feedback_data.isrf_moment[ISRF_MOMENT_LW].specific_flux[0];
  ret[1] =
      rescale * p->feedback_data.isrf_moment[ISRF_MOMENT_LW].specific_flux[1];
  ret[2] =
      rescale * p->feedback_data.isrf_moment[ISRF_MOMENT_LW].specific_flux[2];
}

/**
 * @brief See #feedback_get_part_specific_flux_PE, Lyman-Werner-band
 * photon-number moment.
 *
 * @param p The #part to query.
 * @param ret (return) The three components.
 */
void feedback_get_part_specific_flux_LW_PHOTON(const struct part *p,
                                               float *ret) {
  const float rescale = p->feedback_data.c_hyp;
  ret[0] = rescale *
           p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].specific_flux[0];
  ret[1] = rescale *
           p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].specific_flux[1];
  ret[2] = rescale *
           p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].specific_flux[2];
}

/**
 * @brief Most negative PE-band specific energy written since the previous
 * snapshot, see
 * #feedback_part_data.isrf_moment[ISRF_MOMENT_PE].u_min_since_snapshot.
 *
 * Values stamped with an older snapshot index belong to an interval that
 * saw no update of this particle, so they read as 0. Always 0 without
 * SWIFT_DEBUG_CHECKS.
 *
 * @param p The #part to query.
 * @param e The #engine.
 * @return Most negative PE-band specific energy since the previous
 * snapshot, or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_u_min_since_snapshot_PE(const struct part *p,
                                                const struct engine *e) {
#ifdef SWIFT_DEBUG_CHECKS
  if (p->feedback_data.u_min_snapshot_index == e->snapshot_output_count)
    return p->feedback_data.isrf_moment[ISRF_MOMENT_PE].u_min_since_snapshot;
#endif
  return 0.f;
}

/**
 * @brief See #feedback_get_part_u_min_since_snapshot_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @param e The #engine.
 * @return Most negative Lyman-Werner-band specific energy since the
 * previous snapshot, or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_u_min_since_snapshot_LW(const struct part *p,
                                                const struct engine *e) {
#ifdef SWIFT_DEBUG_CHECKS
  if (p->feedback_data.u_min_snapshot_index == e->snapshot_output_count)
    return p->feedback_data.isrf_moment[ISRF_MOMENT_LW].u_min_since_snapshot;
#endif
  return 0.f;
}

/**
 * @brief See #feedback_get_part_u_min_since_snapshot_PE, Lyman-Werner-band
 * photon-number moment.
 *
 * @param p The #part to query.
 * @param e The #engine.
 * @return Most negative Lyman-Werner-band photon-number moment since the
 * previous snapshot, or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_u_min_since_snapshot_LW_PHOTON(const struct part *p,
                                                       const struct engine *e) {
#ifdef SWIFT_DEBUG_CHECKS
  if (p->feedback_data.u_min_snapshot_index == e->snapshot_output_count)
    return p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON]
        .u_min_since_snapshot;
#endif
  return 0.f;
}

/**
 * @brief Cumulative PE-band raw injected dose since first init, see
 * #feedback_part_data.isrf_moment[ISRF_MOMENT_PE].cumulative_injected. Always 0
 * without SWIFT_DEBUG_CHECKS.
 *
 * @param p The #part to query.
 * @return Cumulative PE-band raw injected dose since first init, or 0
 * without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_cumulative_injected_PE(const struct part *p) {
#ifdef SWIFT_DEBUG_CHECKS
  return p->feedback_data.isrf_moment[ISRF_MOMENT_PE].cumulative_injected;
#else
  return 0.f;
#endif
}

/**
 * @brief See #feedback_get_part_cumulative_injected_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @return Cumulative Lyman-Werner-band raw injected dose since first init,
 * or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_cumulative_injected_LW(const struct part *p) {
#ifdef SWIFT_DEBUG_CHECKS
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW].cumulative_injected;
#else
  return 0.f;
#endif
}

/**
 * @brief See #feedback_get_part_cumulative_injected_PE, Lyman-Werner-band
 * photon-number moment.
 *
 * @param p The #part to query.
 * @return Cumulative Lyman-Werner-band photon-number moment raw injected
 * dose since first init, or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_cumulative_injected_LW_PHOTON(const struct part *p) {
#ifdef SWIFT_DEBUG_CHECKS
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON]
      .cumulative_injected;
#else
  return 0.f;
#endif
}

/**
 * @brief Cumulative PE-band absorbed/transport-and-dissipation-attributed
 * specific energy since first init, see
 * #feedback_part_data.isrf_moment[ISRF_MOMENT_PE].cumulative_absorbed. Always 0
 * without SWIFT_DEBUG_CHECKS.
 *
 * @param p The #part to query.
 * @return Cumulative PE-band absorbed/dissipation-attributed specific
 * energy since first init, or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_cumulative_absorbed_PE(const struct part *p) {
#ifdef SWIFT_DEBUG_CHECKS
  return p->feedback_data.isrf_moment[ISRF_MOMENT_PE].cumulative_absorbed;
#else
  return 0.f;
#endif
}

/**
 * @brief See #feedback_get_part_cumulative_absorbed_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @return Cumulative Lyman-Werner-band absorbed/dissipation-attributed
 * specific energy since first init, or 0 without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_cumulative_absorbed_LW(const struct part *p) {
#ifdef SWIFT_DEBUG_CHECKS
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW].cumulative_absorbed;
#else
  return 0.f;
#endif
}

/**
 * @brief See #feedback_get_part_cumulative_absorbed_PE, Lyman-Werner-band
 * photon-number moment.
 *
 * @param p The #part to query.
 * @return Cumulative Lyman-Werner-band photon-number moment
 * absorbed/dissipation-attributed specific energy since first init, or 0
 * without SWIFT_DEBUG_CHECKS.
 */
float feedback_get_part_cumulative_absorbed_LW_PHOTON(const struct part *p) {
#ifdef SWIFT_DEBUG_CHECKS
  return p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON]
      .cumulative_absorbed;
#else
  return 0.f;
#endif
}

/**
 * @brief Hyperbolic propagation speed, see #feedback_part_data.c_hyp. Thin
 * dispatch wrapper, same reasoning as #feedback_get_part_u_PE.
 *
 * Shared by both bands, and physical: a physical speed, whichever scheme
 * (#isrf_c_hyp_scheme) set it.
 *
 * @param p The #part to query.
 * @return Hyperbolic propagation speed.
 */
float feedback_get_part_c_hyp(const struct part *p) {
  return p->feedback_data.c_hyp;
}

/**
 * @brief Specific energy owed to a moment by finer neighbours, added by the
 * particle's next update: `c_hyp` times the two pending amounts.
 *
 * @param p The #part to query.
 * @param m The #radiation_isrf_moment.
 * @return The pending specific energy.
 */
float feedback_get_part_pending_specific_energy(const struct part *p, int m) {
  const struct feedback_isrf_moment_data *mo = &p->feedback_data.isrf_moment[m];
  return p->feedback_data.c_hyp *
         (mo->pending_transport_u + mo->pending_dissipation_u);
}

/**
 * @brief Radiation timestep bound of a particle, see
 * #radiation_isrf_part_timestep.
 *
 * Stops the run if the bound is below TimeIntegration:dt_min.
 *
 * @param p The particle to consider.
 * @param e The #engine.
 * @return The radiation timestep bound (before the cosmology factor), or
 *     FLT_MAX if none applies.
 */
float feedback_radiation_compute_part_timestep(const struct part *restrict p,
                                               const struct engine *e) {
  const float dt_isrf = radiation_isrf_part_timestep(p, e);
  /* No bound: FLT_MAX times the cosmology factor must not overflow. */
  if (dt_isrf == FLT_MAX) return dt_isrf;
  /* Compared to dt_min after the cosmology factor, like the other
   * candidates in get_part_timestep(). */
  const float dt_isrf_scaled = dt_isrf * e->cosmology->time_step_factor;
  if (dt_isrf_scaled < e->dt_min)
    error(
        "part (id=%lld) wants an ISRF radiation time-step (%e, %e after "
        "the cosmology factor) below TimeIntegration:dt_min (%e): "
        "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c=%g forces dt_rad = "
        "C_hyp*h/(f*c) below dt_min for this particle's h. Lower the "
        "fraction (dt_rad grows as 1/f), or lower dt_min.",
        p->id, dt_isrf, dt_isrf_scaled, e->dt_min,
        e->feedback_props->ISRF_c_hyp_fixed_fraction_of_c);
  return dt_isrf;
}
