/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_FEEDBACK_GEAR_MECHANICAL_H
#define SWIFT_FEEDBACK_GEAR_MECHANICAL_H

#include "../GEAR/feedback_common.h"
#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <float.h>
#include <strings.h>

void feedback_update_part(struct part *p, struct xpart *xp,
                          const struct engine *e);
void feedback_end_density(struct part *p, struct xpart *xp,
                          const struct engine *e);
void feedback_reset_part(struct part *p, struct xpart *xp,
                         const struct engine *e);
int feedback_is_active(const struct spart *sp, const struct engine *e);
int feedback_should_inject_SN_feedback(const struct spart *sp);
int feedback_should_inject_wind_feedback(const struct spart *sp);
int feedback_should_inject_feedback(const struct spart *sp);
void feedback_init_spart(struct spart *sp);
void feedback_reset_feedback(struct spart *sp,
                             const struct feedback_props *feedback_props);
void feedback_prepare_spart(struct spart *sp,
                            const struct feedback_props *feedback_props);
void feedback_prepare_feedback(struct spart *restrict sp,
                               const struct feedback_props *feedback_props,
                               const struct cosmology *cosmo,
                               const struct unit_system *us,
                               const struct phys_const *phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology);

/**
 * @brief Sets all particle fields to sensible values when the #part has 0
 * neighbours. Nothing to do here.
 *
 * @param p The particle to act upon.
 * @param xp The extra particle to act upon.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
feedback_part_has_no_neighbours(struct part *p, struct xpart *xp,
                                const struct engine *e) {}

/**
 * @brief Finishes the #part gradient calculation. Nothing to do here:
 * this feedback model does not track a propagated radiation flux.
 *
 * @param p The particle.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void feedback_end_gradient(
    struct part *p, const struct engine *e) {}

/**
 * @brief Finishes the #part force calculation. Nothing to do here.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void feedback_end_force(
    struct part *p, const struct engine *e) {}

/**
 * @brief Radiation timestep contribution. This module tracks no propagated
 * radiation flux, so this imposes no timestep limit.
 *
 * @param p The particle to consider.
 * @param e The #engine.
 * @return FLT_MAX, always.
 */
__attribute__((always_inline)) INLINE static float
feedback_compute_part_timestep(const struct part *restrict p,
                               const struct engine *e) {
  return FLT_MAX;
}

/**
 * @brief Re-initialise the gas particle-carried fields related to
 * feedback at the start of each density h-iteration. Nothing to do here.
 *
 * @param p The particle.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void feedback_init_part(
    struct part *p, const struct engine *e) {}

/**
 * @brief First-init of a #part's feedback-model state. Nothing to do here.
 *
 * @param p The particle.
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_part(
    struct part *restrict p) {}

/**
 * @brief Is this star particle done evolving, i.e. finished with its
 * feedback-relevant lifetime?
 *
 * @param sp The #spart to query.
 * @return sp->feedback_data.is_dead.
 */
__attribute__((always_inline)) INLINE static int feedback_is_star_dead(
    const struct spart *sp) {

  return sp->feedback_data.is_dead;
}

/**
 * @brief Is this gas particle currently tagged as HII-ionized?
 *
 * Nothing to do here: this module has no subgrid radiation.
 *
 * @param p The #part to query.
 * @param xp The #part's extended data.
 */
__attribute__((always_inline)) INLINE static char
feedback_is_part_tagged_as_ionized(const struct part *p,
                                   const struct xpart *xp) {
  return 0;
}

/**
 * @brief Id of the star that tagged this gas particle as HII-ionized.
 *
 * Nothing to do here: this module has no subgrid radiation.
 *
 * @param p The #part to query.
 * @param xp The #part's extended data.
 */
__attribute__((always_inline)) INLINE static long long
feedback_get_part_ionized_star_id(const struct part *p,
                                  const struct xpart *xp) {
  return 0;
}
/**
 * @brief Local specific PE-band radiation field. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float feedback_get_part_u_PE(
    const struct part *p) {
  return 0.f;
}

/**
 * @brief Local specific Lyman-Werner-band radiation field. Nothing to do
 * here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float feedback_get_part_u_LW(
    const struct part *p) {
  return 0.f;
}

/**
 * @brief Local Lyman-Werner-band photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_u_LW_PHOTON(const struct part *p) {
  return 0.f;
}

/**
 * @brief Negativity-triggered artificial-dissipation coefficient. Nothing to do
 * here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_dissipation_alpha_PE(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_dissipation_alpha_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_dissipation_alpha_LW(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_dissipation_alpha_PE, Lyman-Werner-band
 * photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_dissipation_alpha_LW_PHOTON(const struct part *p) {
  return 0.f;
}

/**
 * @brief `(1/rho) div(rho F)` accumulator. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_div_specific_flux_PE(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_div_specific_flux_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_div_specific_flux_LW(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_div_specific_flux_PE, Lyman-Werner-band
 * photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_div_specific_flux_LW_PHOTON(const struct part *p) {
  return 0.f;
}

/**
 * @brief Tracked specific flux moment. Nothing to do here.
 *
 * @param p The #part to query.
 * @param ret (return) The three components, zeroed.
 */
__attribute__((always_inline)) INLINE static void
feedback_get_part_specific_flux_PE(const struct part *p, float *ret) {
  ret[0] = 0.f;
  ret[1] = 0.f;
  ret[2] = 0.f;
}

/**
 * @brief See #feedback_get_part_specific_flux_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @param ret (return) The three components, zeroed.
 */
__attribute__((always_inline)) INLINE static void
feedback_get_part_specific_flux_LW(const struct part *p, float *ret) {
  ret[0] = 0.f;
  ret[1] = 0.f;
  ret[2] = 0.f;
}

/**
 * @brief See #feedback_get_part_specific_flux_PE, Lyman-Werner-band
 * photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 * @param ret (return) The three components, zeroed.
 */
__attribute__((always_inline)) INLINE static void
feedback_get_part_specific_flux_LW_PHOTON(const struct part *p, float *ret) {
  ret[0] = 0.f;
  ret[1] = 0.f;
  ret[2] = 0.f;
}

struct engine;

/**
 * @brief Most negative PE-band specific energy since the previous snapshot.
 * Nothing to do here.
 *
 * @param p The #part to query.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_u_min_since_snapshot_PE(const struct part *p,
                                          const struct engine *e) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_u_min_since_snapshot_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_u_min_since_snapshot_LW(const struct part *p,
                                          const struct engine *e) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_u_min_since_snapshot_PE, Lyman-Werner-band
 * photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_u_min_since_snapshot_LW_PHOTON(const struct part *p,
                                                 const struct engine *e) {
  return 0.f;
}

/**
 * @brief Cumulative PE-band raw injected dose since first init. Nothing to
 * do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_cumulative_injected_PE(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_cumulative_injected_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_cumulative_injected_LW(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_cumulative_injected_PE, Lyman-Werner-band
 * photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_cumulative_injected_LW_PHOTON(const struct part *p) {
  return 0.f;
}

/**
 * @brief Cumulative PE-band absorbed/transport-and-dissipation-attributed
 * specific energy since first init. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_cumulative_absorbed_PE(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_cumulative_absorbed_PE, Lyman-Werner band.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_cumulative_absorbed_LW(const struct part *p) {
  return 0.f;
}

/**
 * @brief See #feedback_get_part_cumulative_absorbed_PE, Lyman-Werner-band
 * photon-number moment. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float
feedback_get_part_cumulative_absorbed_LW_PHOTON(const struct part *p) {
  return 0.f;
}

/**
 * @brief Hyperbolic propagation speed of the ISRF. Nothing to do here.
 *
 * @param p The #part to query.
 */
__attribute__((always_inline)) INLINE static float feedback_get_part_c_hyp(
    const struct part *p) {
  return 0.f;
}

/**
 * @brief Get the comoving gas density around the star, averaged with the
 * norms of the
 * isotropic vector weights |w_j|.
 *
 * Only valid after runner_iact_nonsym_feedback_prep3().
 *
 * @param sp The #spart.
 */
INLINE static float feedback_get_weighted_gas_density(const struct spart *sp) {
  if (sp->feedback_data.enrichment_weight <= 0.f) return 0.f;
  return sp->feedback_data.weighted_gas_density /
         sp->feedback_data.enrichment_weight;
}

/**
 * @brief Get the gas metal mass fraction around the star, averaged with the
 * norms of the
 * isotropic vector weights |w_j|.
 *
 * Only valid after runner_iact_nonsym_feedback_prep3().
 *
 * @param sp The #spart.
 */
INLINE static double feedback_get_weighted_gas_metallicity(
    const struct spart *sp) {
  if (sp->feedback_data.enrichment_weight <= 0.f) return 0.;
  return sp->feedback_data.weighted_gas_metallicity /
         sp->feedback_data.enrichment_weight;
}

/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The #feedback_props.
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour(struct feedback_props *feedback,
                                          hid_t h_grp) {
#if FEEDBACK_GEAR_MECHANICAL_MODE == 1
  io_write_attribute_s(h_grp, "Feedback Model", "GEAR-mechanical_1");
#elif FEEDBACK_GEAR_MECHANICAL_MODE == 2
  io_write_attribute_s(h_grp, "Feedback Model", "GEAR-mechanical_2");
#else
  error(
      "This function should be called only with one of the GEAR-mechanical "
      "feedback mode.");
#endif
};

void feedback_compute_scalar_weight(const float r2, const float *dx,
                                    const float hi, const float hj,
                                    const struct spart *restrict si,
                                    const struct part *restrict pj,
                                    double dx_ij_plus[3], double dx_ij_minus[3],
                                    double *scalar_weight_j);

void feedback_compute_vector_weight_non_normalized(
    const float r2, const float *dx, const float hi, const float hj,
    const struct spart *restrict si, const struct part *restrict pj,
    double f_plus_i[3], double f_minus_i[3], double w_j[3]);

void feedback_compute_vector_weight_normalized(const float r2, const float *dx,
                                               const float hi, const float hj,
                                               const struct spart *restrict si,
                                               const struct part *restrict pj,
                                               double w_j_bar[3]);

double feedback_get_physical_SN_terminal_momentum(
    const struct spart *restrict sp, const struct part *restrict p,
    const struct xpart *restrict xp, const struct phys_const *phys_const,
    const struct unit_system *us, const struct feedback_props *feedback_props,
    const struct cosmology *cosmo);

float feedback_get_physical_SN_cooling_radius(const struct spart *restrict sp,
                                              float p_SN_initial,
                                              float p_terminal,
                                              const struct cosmology *cosmo);

float feedback_compute_momentum_correction_factor_for_multiple_sn_events(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo);

void feedback_accumulate_kinetic_energy_for_multiple_sn_events(
    struct xpart *xp, const float mj, const float new_mass,
    const float v_pec[3], const float v_hubble[3], const double dp[3],
    const double dp_ejecta[3]);

float feedback_compute_residual_internal_energy_for_multiple_sn_events(
    const struct part *p, const struct xpart *xp, const struct cosmology *cosmo,
    const float old_mass, const float new_mass, const float f_corr);
#endif /* SWIFT_FEEDBACK_GEAR_MECHANICAL_H */
