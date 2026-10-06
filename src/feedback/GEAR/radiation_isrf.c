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
 * @file src/feedback/GEAR/radiation_isrf.c
 * @brief Receiver-side LW/PE dust extinction and hyperbolic M1
 * propagation physics for GEAR: closed by a variable Eddington tensor
 * (not a fixed E/3) that adapts between the free-streaming and diffusive
 * limits.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "active.h"
#include "chemistry.h"
#include "cooling.h"
#include "cosmology.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "radiation.h"
#include "radiation_isrf.h"
#include "radiation_propagation_iact.h"
#include "timeline.h"

#include <float.h>
#include <math.h>

/**
 * @brief First-init of a #part's LW/PE radiation-field state.
 *
 * Does not zero #feedback_isrf_moment_data.u, which an IC may supply through
 * the optional "PESpecificEnergy"/"LWSpecificEnergy" fields. u_prev is seeded
 * from u, because the initial pass rebuilds `u` from `u_prev` before the first
 * snapshot. The specific flux is always zeroed.
 *
 * @param p The #part to initialise.
 */
void radiation_first_init_part(struct part *restrict p) {
  struct feedback_part_data *fd = &p->feedback_data;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    moment->u_prev = moment->u;
    moment->specific_flux[0] = 0.f;
    moment->specific_flux[1] = 0.f;
    moment->specific_flux[2] = 0.f;
    moment->u_dose_reservoir = 0.f;
    moment->u_source_rate = 0.f;
    moment->dissipation_u = 0.f;
    moment->div_specific_flux = 0.f;
#ifdef SWIFT_DEBUG_CHECKS
    moment->u_min_since_snapshot = 0.f;
    moment->cumulative_injected = 0.f;
    moment->cumulative_absorbed = 0.f;
#endif
  }
  /* The LW and LW photon IC fields are optional and independent. A zero one
   * is set equal to the other, in either direction, so the two stay
   * nonzero together: #radiation_isrf_part_timestep relies on this. */
  if (fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u == 0.f &&
      fd->isrf_moment[ISRF_MOMENT_LW].u != 0.f) {
    fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u =
        fd->isrf_moment[ISRF_MOMENT_LW].u;
    fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u_prev =
        fd->isrf_moment[ISRF_MOMENT_LW].u;
  } else if (fd->isrf_moment[ISRF_MOMENT_LW].u == 0.f &&
             fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u != 0.f) {
    fd->isrf_moment[ISRF_MOMENT_LW].u =
        fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u;
    fd->isrf_moment[ISRF_MOMENT_LW].u_prev =
        fd->isrf_moment[ISRF_MOMENT_LW_PHOTON].u;
  }
  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++) {
    struct feedback_isrf_operator_data *op = &fd->isrf_operator[o];
    op->kappa = 0.f;
    op->dissipation_alpha_trigger = 0.f;
    op->dissipation_alpha_floor = 0.f;
  }
#ifdef SWIFT_DEBUG_CHECKS
  fd->u_min_snapshot_index = 0;
#endif
  fd->ISRF_last_touch_ti = -1;
  /* MPI foreign particles are not bzero'd: keep them from reading as
     illuminated. */
  fd->is_illuminated_ISRF = 0;
  fd->ISRF_illumination_end_ti = -1;
  /* Placeholder density: a 0 seed makes `1/rho_i` in grad(u) infinite
   * during the initial pass, before the first snapshot. Safe since `u` and
   * `F` are still 0. */
  fd->rho_prev = 1.0f;
  fd->c_hyp = 0.f;
  fd->dt_prev = 0.f;
  fd->ISRF_reservoir_end_ti = -1;
  radiation_init_part_propagation(p);
  radiation_cache_m1_closure_part(p);
}

/**
 * @brief Once per step, before the density loop: snapshot `u`, cache the
 * absorption rates, the comoving density, `dt_prev` and (fixed_fraction
 * scheme) `c_hyp`, zero the gradient accumulators and, for active particles,
 * draw down `u_source_rate` from the dose reservoir.
 *
 * Runs at drift time, the last point where the metal fraction and `p->rho`
 * still hold the previous step's converged values.
 *
 * @param p The #part to reset.
 * @param e The #engine.
 */
void radiation_snapshot_part_propagation(struct part *p,
                                         const struct engine *e) {
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &p->feedback_data.isrf_moment[m];
    moment->u_prev = moment->u;

    moment->grad_u[0] = 0.f;
    moment->grad_u[1] = 0.f;
    moment->grad_u[2] = 0.f;

    /* The force loop also runs once per step, with no h-iteration redo. */
    moment->dissipation_u = 0.f;
  }

  /* Cached even without ISRF_propagation: the gradient loop always reads
   * it, and it must never be 0. */
  const float rho_comoving = hydro_get_comoving_density(p);
  p->feedback_data.rho_prev = rho_comoving > 0.f ? rho_comoving : 1.0f;

  if (!e->feedback_props->ISRF_propagation) {
    p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa = 0.f;
    p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa = 0.f;
    return;
  }

  const float rho_phys = hydro_get_physical_density(p, e->cosmology);
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  /* Resolved value, not the raw `-1` sentinel. */
  const float local_dust_to_gas_ratio =
      (float)e->cooling_func->chemistry_data.local_dust_to_gas_ratio;
  p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa =
      radiation_get_part_linear_absorption_rate(
          e->internal_units, e->physical_constants, Z, rho_phys,
          RADIATION_SIGMA_D_PE_CGS, local_dust_to_gas_ratio);
  p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa =
      radiation_get_part_linear_absorption_rate(
          e->internal_units, e->physical_constants, Z, rho_phys,
          RADIATION_SIGMA_D_LW_CGS, local_dust_to_gas_ratio);

  /* Physical timestep of this particle. Can be exactly 0 before its first
   * step: cached unfloored so readers can detect that case. */
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_step = get_integer_timestep(p->time_bin);
  /* Start of the step ending now. Using `ti_current` instead would
   * over-drain the reservoir one sub-step early. */
  const integertime_t ti_begin =
      get_integer_time_begin(e->ti_current, p->time_bin);
  float dt_phys;
  if (with_cosmology) {
    dt_phys = (float)cosmology_get_delta_time(e->cosmology, ti_begin,
                                              ti_begin + ti_step);
  } else {
    dt_phys = (float)get_timestep(p->time_bin, e->time_base);
  }

  /* The kernel-local scheme sets c_hyp in radiation_end_density_propagation.
   * The fixed-fraction scheme sets it here. Any other value is invalid. */
  const int c_hyp_scheme = e->feedback_props->ISRF_c_hyp_scheme;
  if (c_hyp_scheme == isrf_c_hyp_scheme_fixed_fraction) {
    /* Uniform c_hyp. The CFL is enforced by radiation_isrf_part_timestep(). */
    p->feedback_data.c_hyp = e->feedback_props->ISRF_c_hyp_fixed_fraction_of_c *
                             (float)e->physical_constants->const_speed_light_c;
  } else if (c_hyp_scheme != isrf_c_hyp_scheme_kernel_local_reduced_flux) {
    error(
        "Invalid GEARFeedback:ISRF_c_hyp_scheme (internal value %d): must be "
        "fixed_fraction (2) or kernel_local (4).",
        c_hyp_scheme);
  }
  p->feedback_data.dt_prev = dt_phys;

  /* Drawdown for active particles only. #radiation_part_has_no_neighbours
   * refunds `u_source_rate * dt_prev`, so the rate must use the same raw
   * dt_phys. At dt_phys == 0 nothing is drawn and the dose is kept. */
  struct feedback_part_data *fd = &p->feedback_data;
  if (part_is_active(p, e) && dt_phys > 0.f &&
      (fd->isrf_moment[ISRF_MOMENT_PE].u_dose_reservoir > 0.f ||
       fd->isrf_moment[ISRF_MOMENT_LW].u_dose_reservoir > 0.f)) {
    double t_rem;
    if (fd->ISRF_reservoir_end_ti <= ti_begin) {
      t_rem = 0.0;
    } else if (with_cosmology) {
      t_rem = cosmology_get_delta_time(e->cosmology, ti_begin,
                                       fd->ISRF_reservoir_end_ti);
    } else {
      t_rem = (double)(fd->ISRF_reservoir_end_ti - ti_begin) * e->time_base;
    }
    const float f =
        (t_rem <= (double)dt_phys) ? 1.f : (float)((double)dt_phys / t_rem);
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
      struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
      moment->u_source_rate = f * moment->u_dose_reservoir / dt_phys;
      moment->u_dose_reservoir -= f * moment->u_dose_reservoir;
    }
  } else if (part_is_active(p, e)) {
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
      fd->isrf_moment[m].u_source_rate = 0.f;
  }
}

/**
 * @brief Radiation timestep bound `C_hyp*h_i/(f*c)` for the fixed_fraction
 * c_hyp scheme (#feedback_props.ISRF_c_hyp_fixed_fraction_of_c = f), which
 * has no built-in guarantee that `f*c*dt_i <= C_hyp*h_i`.
 *
 * FLT_MAX (no constraint) if f is 0, ISRF_propagation is off, the debug
 * off-switch is set, or the particle is not near a field. Near a field means
 * it is illuminated, any band's u != 0, or any band's `ngb_mean_abs_u_V > 0`.
 *
 * @param p The #part to consider.
 * @param e The #engine.
 * @return The radiation timestep bound, or FLT_MAX if none applies.
 */
float radiation_isrf_part_timestep(const struct part *restrict p,
                                   const struct engine *e) {
  const float f = e->feedback_props->ISRF_c_hyp_fixed_fraction_of_c;
  /* Gating on f alone is safe because feedback_props_check_c_hyp_scheme()
   * makes f > 0 imply the fixed_fraction scheme. */
  if (f <= 0.f) return FLT_MAX;
  if (!e->feedback_props->ISRF_propagation) return FLT_MAX;
  if (e->feedback_props->ISRF_c_hyp_timestep_term_off_for_debugging)
    return FLT_MAX;

  const struct feedback_part_data *fd = &p->feedback_data;
  /* LW_PHOTON.u is nonzero iff LW.u is (see #radiation_first_init_part),
   * so the LW test covers it. */
  const int near_field =
      fd->is_illuminated_ISRF || fd->isrf_moment[ISRF_MOMENT_PE].u != 0.f ||
      fd->isrf_moment[ISRF_MOMENT_LW].u != 0.f ||
      fd->isrf_operator[ISRF_OPERATOR_PE].ngb_mean_abs_u_V > 0.f ||
      fd->isrf_operator[ISRF_OPERATOR_LW].ngb_mean_abs_u_V > 0.f;
  if (!near_field) return FLT_MAX;

  const float h_phys = (float)e->cosmology->a * p->h;
  /* Explicit branch: under -ffast-math a clamp would not absorb a NaN. */
  if (h_phys <= 0.f) return FLT_MAX;

  const float c_M = f * (float)e->physical_constants->const_speed_light_c;
  return e->feedback_props->ISRF_c_hyp_margin * h_phys / c_M;
}

/**
 * @brief Reset the per-h-iteration accumulators: the kernel-mean field and
 * #feedback_part_data.max_ngb_time_bin (set to the particle's own time bin).
 *
 * Safe to call several times per step. The force-loop accumulators are reset
 * elsewhere (`dissipation_u` in #radiation_snapshot_part_propagation,
 * `div_specific_flux` in #radiation_end_gradient_propagation).
 *
 * @param p The #part to reset.
 */
void radiation_init_part_propagation(struct part *p) {
  p->feedback_data.max_ngb_time_bin = p->time_bin;
  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
    p->feedback_data.isrf_operator[o].ngb_mean_abs_u_V = 0.f;
}

/**
 * @brief Set the kernel-local `c_hyp` once the density loop has converged,
 * then refresh the M1 closure cache.
 *
 * `c_hyp_i = min(C_hyp*h_i/dt_max(i), c)`, with `dt_max(i)` the physical
 * duration of the step at #feedback_part_data.max_ngb_time_bin that contains
 * `e->ti_current`. Acts only for the kernel_local scheme; any invalid scheme
 * stops the run. Neighbours with `H_i <= r < H_j` are not covered, see
 * #feedback_part_data.c_hyp.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_density_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;
  const int c_hyp_scheme = e->feedback_props->ISRF_c_hyp_scheme;
  if (c_hyp_scheme == isrf_c_hyp_scheme_fixed_fraction) return;
  if (c_hyp_scheme != isrf_c_hyp_scheme_kernel_local_reduced_flux)
    error(
        "Invalid GEARFeedback:ISRF_c_hyp_scheme (internal value %d): must be "
        "fixed_fraction (2) or kernel_local (4).",
        c_hyp_scheme);

  struct feedback_part_data *fd = &p->feedback_data;

  float dt_max;
  if (fd->max_ngb_time_bin == p->time_bin) {
    dt_max = fd->dt_prev;
  } else {
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const integertime_t ti_step_max =
        get_integer_timestep(fd->max_ngb_time_bin);
    const integertime_t ti_begin_max =
        get_integer_time_begin(e->ti_current, fd->max_ngb_time_bin);
    if (with_cosmology) {
      dt_max = (float)cosmology_get_delta_time(e->cosmology, ti_begin_max,
                                               ti_begin_max + ti_step_max);
    } else {
      dt_max = (float)get_timestep(fd->max_ngb_time_bin, e->time_base);
    }
  }

  const float h_phys = (float)e->cosmology->a * p->h;
  float c_hyp;
  if (dt_max <= 0.f) {
    /* dt_prev is 0 before the first drift. Set c directly: min(+inf, c) is
     * not guaranteed under -ffast-math. */
    c_hyp = (float)e->physical_constants->const_speed_light_c;
  } else {
    c_hyp = e->feedback_props->ISRF_c_hyp_margin * h_phys / dt_max;
    c_hyp = min(c_hyp, (float)e->physical_constants->const_speed_light_c);
  }

  fd->c_hyp = c_hyp;
  radiation_cache_m1_closure_part(p);
}

/**
 * @brief Exact-relaxation integrating factor `phi(a) = (1 - exp(-a))/a`,
 * used by the `u`- and `F`-updates. A Taylor series is used near `a = 0`.
 *
 * @param a Dimensionless relaxation depth `c_hyp*(kappa + H/c)*dt` (or
 * `c_hyp*kappa*Delta_t_star` at injection). Always `>= 0`.
 * @return phi(a).
 */
float radiation_relaxation_phi_factor(float a) {
  if (a < 1e-6f) return 1.0f - 0.5f * a + (1.0f / 6.0f) * a * a;
  return -expm1f(-a) / a;
}

/**
 * @brief Double-precision #radiation_relaxation_phi_factor for the
 * `u`-update, whose depth `a` can be ~1e-8, below float32 resolution.
 *
 * @param a Dimensionless relaxation depth, see
 * #radiation_relaxation_phi_factor. Always `>= 0`.
 * @return phi(a), in double precision.
 */
static double radiation_relaxation_phi_factor_double(double a) {
  if (a < 1e-6) return 1.0 - 0.5 * a + (1.0 / 6.0) * a * a;
  return -expm1(-a) / a;
}

/**
 * @brief The three outcomes of the M1 flux limiter, so that moments sharing
 * an operator apply the same branch.
 */
enum radiation_isrf_flux_limiter_state {
  ISRF_LIMITER_ZERO, /*!< `u <= 0`: every flux component is zeroed. */
  ISRF_LIMITER_SKIP, /*!< `|F|^2 <= 0`: left untouched, no multiply. */
  ISRF_LIMITER_SCALE /*!< Scaled by #scale below. */
};

/**
 * @brief M1 flux limiter decision for one particle, one band:
 * `F <- F*min(1, u/|F|)` for `u > 0`, `F <- 0` for `u <= 0`, on the reduced
 * flux `Ft = F_true/c_hyp`. The scale is applied by
 * #radiation_apply_flux_limiter_band, possibly to another band.
 *
 * `F.F` and the ratio are in double: in float32 `F.F` underflows for
 * `|F| < 1.1e-19` and the limiter would be skipped for a nonzero flux.
 *
 * @param u This band's specific field `u^n`.
 * @param F This particle's tracked reduced flux (this band), unmodified.
 * @param scale (return) The multiplier to apply, valid only when the
 * return value is #ISRF_LIMITER_SCALE.
 * @return Which of the three outcomes applies.
 */
__attribute__((
    always_inline)) INLINE static enum radiation_isrf_flux_limiter_state
radiation_compute_flux_limiter_scale_band(float u, const float F[3],
                                          float *scale) {

  if (u <= 0.f) return ISRF_LIMITER_ZERO;

  const double F2 = (double)F[0] * (double)F[0] + (double)F[1] * (double)F[1] +
                    (double)F[2] * (double)F[2];
  if (F2 <= 0.) return ISRF_LIMITER_SKIP;

  *scale = (float)min(1., (double)u / sqrt(F2));
  return ISRF_LIMITER_SCALE;
}

/**
 * @brief Apply a flux-limiter decision
 * (#radiation_compute_flux_limiter_scale_band) to one band's flux. The SKIP
 * case does no multiply, since "scale by 1.0" may differ under `-ffast-math`.
 *
 * @param state The decision, from #radiation_compute_flux_limiter_scale_band.
 * @param scale The multiplier, meaningful only under #ISRF_LIMITER_SCALE.
 * @param F (in/out) This particle's tracked flux (this band).
 */
__attribute__((always_inline)) INLINE static void
radiation_apply_flux_limiter_band(enum radiation_isrf_flux_limiter_state state,
                                  float scale, float F[3]) {

  switch (state) {
    case ISRF_LIMITER_ZERO:
      F[0] = 0.f;
      F[1] = 0.f;
      F[2] = 0.f;
      break;
    case ISRF_LIMITER_SKIP:
      break;
    case ISRF_LIMITER_SCALE:
      F[0] *= scale;
      F[1] *= scale;
      F[2] *= scale;
      break;
  }
}

/**
 * @brief Exact-relaxation update of #feedback_isrf_moment_data.u from the
 * `div(F)` and dissipation accumulators and this step's `u_source_rate`.
 *
 * `u^{n+1} = e*(u_prev + dt*phi*diss) + dt*phi*((c_hyp/c)*source_rate -
 * div_F)`, with `e = exp(-a)`, `phi = (1-e)/a` and
 * `a = (c_hyp*kappa + lambda(m)*(c_hyp/c)*H)*dt`. The Hubble term is dilated
 * by `c_hyp/c` like absorption and injection. `div_F` comes from the flux
 * already relaxed in #radiation_end_gradient_propagation. Derivation:
 * theory/GEAR/Radiation/02_fuv_isrf.tex, "From the Eulerian moments to the
 * system we integrate" and "Time integration: exact relaxation, staggered".
 *
 * `lambda(m)` is the moment's band-edge weight (photons redshift across the
 * PE/LW edges); the LW loss at 11.2 eV is transferred to PE (`transfer`).
 * Injection deposits the raw dose, so the `c_hyp/c` rescale is applied only
 * here.
 *
 * Idempotent: `u` is rebuilt from `u_prev`. The debug-only ledger counters
 * are the exception and need exactly one call per active particle per step.
 * Reads `dt_prev`, `c_hyp` and `kappa`, and takes no `dt` so that it stays
 * consistent with the flux update. No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_force_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  const double dt = (double)fd->dt_prev;
  const double c_hyp = (double)fd->c_hyp;
  const double rescale = c_hyp / e->physical_constants->const_speed_light_c;
  const double H = e->cosmology->H;
  /* Double, since the relaxation depth can be ~1e-8. H = 0 without
   * cosmology, which makes this exactly 0. */
  const double H_dilated = rescale * H;

#ifdef SWIFT_DEBUG_CHECKS
  if (fd->u_min_snapshot_index != e->snapshot_output_count) {
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
      fd->isrf_moment[m].u_min_since_snapshot = 0.f;
    fd->u_min_snapshot_index = e->snapshot_output_count;
  }
#endif

  /* Per moment, not per operator: LW and LW_PHOTON share `kappa` but have
   * different weights. */
  const double lambda[ISRF_MOMENT_COUNT] = {
      e->feedback_props->band_edge_weight_pe,
      e->feedback_props->band_edge_weight_lw,
      e->feedback_props->band_edge_photon_weight_lw};

  /* LW-to-PE band-edge transfer, computed before the moment loop from LW's
   * own relaxation depth `a_lw`, so it does not depend on moment order.
   * `f_edge` reuses `a_lw`: independently written expressions can disagree
   * under -ffast-math and unbalance the LW loss and the PE gain. Skipped
   * for H = 0. */
  double transfer = 0.;
  if (H_dilated != 0.) {
    const struct feedback_isrf_operator_data *op_lw =
        &fd->isrf_operator[radiation_isrf_moment_to_operator[ISRF_MOMENT_LW]];
    const struct feedback_isrf_moment_data *moment_lw =
        &fd->isrf_moment[ISRF_MOMENT_LW];
    const double a_lw =
        (c_hyp * (double)op_lw->kappa + lambda[ISRF_MOMENT_LW] * H_dilated) *
        dt;
    /* expm1 avoids the cancellation in `1 - exp(-a)` at small `a`. */
    const double one_minus_decay_lw = -expm1(-a_lw);
    const double phi_lw = radiation_relaxation_phi_factor_double(a_lw);
    /* Fraction of LW's `u_prev` and source/transport terms absorbed this
     * step. Not debug-only: it feeds the transfer. */
    const double absorbed_lw =
        (moment_lw->u_prev + dt * phi_lw * (double)moment_lw->dissipation_u) *
            one_minus_decay_lw +
        (rescale * (double)moment_lw->u_source_rate -
         (double)moment_lw->div_specific_flux) *
            dt * (1. - phi_lw);
    /* Avoid `0/0` at `a_lw = 0` (metal-free particle): a NaN cannot be
     * caught under -ffast-math. */
    if (a_lw > 0.) {
      const double f_edge =
          (lambda[ISRF_MOMENT_LW] - 1.) * H_dilated * dt / a_lw;
      /* Signed: clamping would break the 11.2 eV cancellation in G0 and the
       * ledger. */
      transfer = f_edge * absorbed_lw;
    }
  }

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    const struct feedback_isrf_operator_data *op =
        &fd->isrf_operator[radiation_isrf_moment_to_operator[m]];
    const double a = (c_hyp * (double)op->kappa + lambda[m] * H_dilated) * dt;
    const double decay = exp(-a);
    const double phi = radiation_relaxation_phi_factor_double(a);

#ifdef SWIFT_DEBUG_CHECKS
    /* Energy ledger: raw dose attempted and fraction absorbed, from the same
     * `decay` and `phi` as the update. Stored as float. */
    const double one_minus_decay = -expm1(-a);
    moment->cumulative_injected +=
        (float)(dt * rescale * (double)moment->u_source_rate);
    moment->cumulative_absorbed +=
        (float)((moment->u_prev + dt * phi * (double)moment->dissipation_u) *
                    one_minus_decay +
                (rescale * (double)moment->u_source_rate -
                 (double)moment->div_specific_flux) *
                    dt * (1. - phi));
    /* LW's absorbed amount already contains the transfer, so book it as PE
     * injection to keep `E + Abs - Inj = 0`. */
    if (m == ISRF_MOMENT_PE) moment->cumulative_injected += (float)transfer;
#endif

    moment->u =
        decay * (moment->u_prev + dt * phi * (double)moment->dissipation_u) +
        dt * phi *
            (rescale * (double)moment->u_source_rate -
             (double)moment->div_specific_flux);

    /* Raw add, not subject to PE's absorption this step. A separate statement
     * so -ffast-math cannot merge it into a form that is not 0 at H = 0. */
    if (m == ISRF_MOMENT_PE) moment->u += transfer;

#ifdef SWIFT_DEBUG_CHECKS
    if (moment->u < (double)moment->u_min_since_snapshot)
      moment->u_min_since_snapshot = (float)moment->u;
#endif
  }
}

/**
 * @brief Refund this step's dose-reservoir drawdown (`source_rate*dt_prev`)
 * for a #part whose density h-iteration gives up with no neighbours, and zero
 * the rates. No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_part_has_no_neighbours(struct part *p, const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    moment->u_dose_reservoir += moment->u_source_rate * fd->dt_prev;
    moment->u_source_rate = 0.f;
  }
}

/**
 * @brief One band's negativity-triggered dissipation coefficient: raised
 * instantly to the target, or decayed toward it. Compares `u_V = rho_prev*u^n`
 * with the neighbours' kernel mean. The result is used by this step's force
 * loop, so an undershoot is corrected one step later.
 *
 * @param u_V This band's volumetric field, rho_prev*u^n.
 * @param ngb_mean_abs_u_V This band's kernel-mean |rho_prev*u_prev| scratch
 * accumulator (radiation_propagation_iact.h).
 * @param alpha_prev This band's
 * #feedback_isrf_operator_data.dissipation_alpha_trigger from the previous
 * step (not the floor-combined value).
 * @param alpha_max #feedback_props.ISRF_dissipation_alpha_max.
 * @param eps_1 #feedback_props.ISRF_dissipation_negativity_threshold.
 * @param c_hyp The particle's own #c_hyp.
 * @param kappa This band's #feedback_isrf_operator_data.kappa.
 * @param dt The particle's own #dt_prev.
 * @param h_phys The particle's own physical smoothing length.
 * @return This step's updated dissipation coefficient for this band.
 */
__attribute__((always_inline)) INLINE static float
radiation_update_dissipation_alpha_band(float u_V, float ngb_mean_abs_u_V,
                                        float alpha_prev, float alpha_max,
                                        float eps_1, float c_hyp, float kappa,
                                        float dt, float h_phys) {

  /* Branch on the reference: for a faint band `u_V` underflows to `-0.f`
   * while the sign test sees a negative double, giving `0/0`. With no
   * resolvable field scale the trigger stays zero. */
  const float u_V_ref = max(ngb_mean_abs_u_V, -u_V);
  const float eps = (u_V < 0.f && u_V_ref > 0.f) ? -u_V / u_V_ref : 0.f;
  const float x = min(eps / eps_1, 1.f);
  const float alpha_aim = alpha_max * x * x * (3.f - 2.f * x);

  if (alpha_aim >= alpha_prev) return alpha_aim;

  const float a_kappa = c_hyp * kappa * dt;
  const float decay =
      expf(-c_hyp * dt / (RADIATION_ISRF_DISSIPATION_DECAY_LENGTH * h_phys) -
           a_kappa);
  return alpha_aim + (alpha_prev - alpha_aim) * decay;
}

/**
 * @brief The `h/lambda`-gated dissipation floor under the trigger of
 * #radiation_update_dissipation_alpha_band. It damps the dispersive wake of
 * an optically-thin front, where the trigger is zero, and rolls off as
 * `alpha_floor/(1 + (h*kappa/eps_lambda)^4)`.
 *
 * @param kappa This band's #feedback_isrf_operator_data.kappa.
 * @param h_phys The particle's own physical smoothing length.
 * @param alpha_floor #feedback_props.ISRF_dissipation_alpha_floor.
 * @param eps_lambda #feedback_props.ISRF_dissipation_floor_h_over_lambda.
 * @return This band's floor value for this step, stored on its own
 * #feedback_part_data field for the force loop to combine.
 */
__attribute__((always_inline)) INLINE static float
radiation_dissipation_alpha_floor_band(float kappa, float h_phys,
                                       float alpha_floor, float eps_lambda) {

  const float x = h_phys * kappa / eps_lambda;
  const float x2 = x * x;
  return alpha_floor / (1.f + x2 * x2);
}

/**
 * @brief Flux-relaxation residual gate on the dissipation floor. It lowers
 * the floor on a particle whose flux is at the fixed point `w*Ft + grad_u = 0`
 * of the unlimited flux update (`w = kappa + lambda*H/c`), and keeps the full
 * floor (`s = 1`) on a fresh front. It can only lower the floor. See
 * theory/GEAR/Radiation/02_fuv_isrf.tex, "The flux-relaxation-residual gate".
 *
 * `R = |w*Ft + grad_u| / (w*|Ft| + |grad_u|)`, `s = min(1, (R/eps_R)^2)`.
 * Returns 1 for `eps_R <= 0`, `c_hyp <= 0` or `w <= 0`, and 0 for a quiescent
 * particle (zero denominator). `R` and `s` are in double so the squared
 * denominator cannot underflow under fast-math.
 *
 * @param F This band's #feedback_isrf_moment_data.specific_flux (the reduced
 * flux), from BEFORE this step's own update (already post-limiter).
 * @param grad_u This band's #feedback_isrf_moment_data.grad_u accumulator.
 * @param c_hyp The particle's own #c_hyp. Only its sign is read.
 * @param kappa This band's #feedback_isrf_operator_data.kappa.
 * @param H The Hubble rate, #cosmology.H.
 * @param c The true speed of light, #phys_const.const_speed_light_c.
 * @param lambda The owning moment's band-edge weight.
 * @param eps_R #feedback_props.ISRF_dissipation_floor_relaxation_residual.
 * @return The floor-aim multiplier `s`, in `[0, 1]`.
 */
__attribute__((always_inline)) INLINE static float
radiation_dissipation_floor_relaxation_gate(const float F[3],
                                            const float grad_u[3], float c_hyp,
                                            float kappa, float H, float c,
                                            float lambda, float eps_R) {

  if (eps_R <= 0.f) return 1.f;
  if (c_hyp <= 0.f) return 1.f;

  /* Multiplied through by w to avoid float32 overflow of `grad_u/w` at
   * small kappa. Uses the true speed c. */
  const float w = kappa + lambda * H / c;

  /* No relaxation timescale (kappa = 0 and H = 0): full floor. Also stops F
   * dropping out of R. */
  if (w <= 0.f) return 1.f;

  /* Double: float32 squares underflow below 1.1e-19 and would send a front
   * to the quiescent branch. */
  const double w_d = w;
  const double F_d[3] = {F[0], F[1], F[2]};
  const double G_d[3] = {grad_u[0], grad_u[1], grad_u[2]};
  const double wx = w_d * F_d[0] + G_d[0];
  const double wy = w_d * F_d[1] + G_d[1];
  const double wz = w_d * F_d[2] + G_d[2];
  const double num = sqrt(wx * wx + wy * wy + wz * wz);
  const double F_norm =
      sqrt(F_d[0] * F_d[0] + F_d[1] * F_d[1] + F_d[2] * F_d[2]);
  const double G_norm =
      sqrt(G_d[0] * G_d[0] + G_d[1] * G_d[1] + G_d[2] * G_d[2]);
  const double den = w_d * F_norm + G_norm;

  /* Quiescent particle: at the fixed point, `R = 0`. A denominator epsilon
   * would underflow once squared and give `0/0`. */
  if (den <= 0.) return 0.f;

  const double R = num / den;
  const double ratio2 = (R / eps_R) * (R / eps_R);
  return (float)min(1.0, ratio2);
}

/**
 * @brief Exact-relaxation update of the reduced flux
 * #feedback_isrf_moment_data.specific_flux from `grad(u)`, followed by the M1
 * flux limiter against `u^n`. Also updates the two dissipation coefficients
 * (#radiation_update_dissipation_alpha_band,
 * #radiation_dissipation_alpha_floor_band) and zeroes `div_specific_flux` for
 * the force loop.
 *
 * Runs once per step in the extra ghost. The limiter scale is decided once per
 * operator from its owning moment and applied to every moment sharing it.
 * The stored flux is `Ft = F_true/c_hyp` and obeys
 * `Ft_new = e*Ft - c_hyp*dt*phi*grad(u)` with `a = c_hyp*(kappa + H/c)*dt`;
 * the limiter bound is `|Ft| <= u`. No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_gradient_propagation(struct part *p,
                                        const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  const float dt = fd->dt_prev;
  const float c_hyp = fd->c_hyp;
  const float H = (float)e->cosmology->H;
  const float c = (float)e->physical_constants->const_speed_light_c;
  const float H_dilated = (c_hyp / c) * H;
  const float h_phys = (float)e->cosmology->a * p->h;

  const float alpha_pin =
      e->feedback_props->ISRF_dissipation_alpha_pin_for_debugging;
  const float alpha_max = e->feedback_props->ISRF_dissipation_alpha_max;
  const float eps_1 = e->feedback_props->ISRF_dissipation_negativity_threshold;
  const float alpha_floor = e->feedback_props->ISRF_dissipation_alpha_floor;
  const float eps_lambda =
      e->feedback_props->ISRF_dissipation_floor_h_over_lambda;
  const float eps_R =
      e->feedback_props->ISRF_dissipation_floor_relaxation_residual;

  /* Per-moment weight, as in #radiation_end_force_propagation: `F` and `u`
   * of one moment must relax at the same rate. */
  const float lambda[ISRF_MOMENT_COUNT] = {
      (float)e->feedback_props->band_edge_weight_pe,
      (float)e->feedback_props->band_edge_weight_lw,
      (float)e->feedback_props->band_edge_photon_weight_lw};

  /* Pre-update flux per operator, written by the owning moment only. */
  float F_old_stash[ISRF_OPERATOR_COUNT][3];

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    const int o = radiation_isrf_moment_to_operator[m];
    const struct feedback_isrf_operator_data *op = &fd->isrf_operator[o];

    const float a = (c_hyp * op->kappa + lambda[m] * H_dilated) * dt;
    const float decay = expf(-a);
    const float phi = radiation_relaxation_phi_factor(a);
    const float coeff = c_hyp * dt * phi;

    /* For the relaxation-residual gate: `u^n` came from this flux. */
    if (m == (int)radiation_isrf_operator_owner[o]) {
      F_old_stash[o][0] = moment->specific_flux[0];
      F_old_stash[o][1] = moment->specific_flux[1];
      F_old_stash[o][2] = moment->specific_flux[2];
    }

    for (int k = 0; k < 3; k++) {
      moment->specific_flux[k] =
          decay * moment->specific_flux[k] - coeff * moment->grad_u[k];
    }

    /* Not zeroed at drift, which would blank the output field at snapshots. */
    moment->div_specific_flux = 0.f;
  }

  /* The limiter decision and dissipation coefficients are operator state:
   * computed once per operator from its owning moment, applied to all
   * moments sharing it. */
  enum radiation_isrf_flux_limiter_state limiter_state[ISRF_OPERATOR_COUNT];
  float limiter_scale[ISRF_OPERATOR_COUNT];

  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++) {
    const int m = radiation_isrf_operator_owner[o];
    const struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    struct feedback_isrf_operator_data *op = &fd->isrf_operator[o];

    limiter_state[o] = radiation_compute_flux_limiter_scale_band(
        moment->u, moment->specific_flux, &limiter_scale[o]);

    const float u_V = fd->rho_prev * moment->u;

    if (alpha_pin > 0.f) {
      /* Debug pin: held in the trigger component, floor zeroed, so it stays
       * uniform and is not rolled off by `h/lambda`. */
      op->dissipation_alpha_trigger = alpha_pin;
      op->dissipation_alpha_floor = 0.f;
    } else {
      op->dissipation_alpha_trigger = radiation_update_dissipation_alpha_band(
          u_V, op->ngb_mean_abs_u_V, op->dissipation_alpha_trigger, alpha_max,
          eps_1, c_hyp, op->kappa, dt, h_phys);

      /* Kept apart from the trigger: the force loop combines them. */
      const float s = radiation_dissipation_floor_relaxation_gate(
          F_old_stash[o], moment->grad_u, c_hyp, op->kappa, H, c,
          lambda[radiation_isrf_operator_owner[o]], eps_R);
      op->dissipation_alpha_floor =
          s * radiation_dissipation_alpha_floor_band(op->kappa, h_phys,
                                                     alpha_floor, eps_lambda);
    }
  }

  /* Apply each operator's decision to all its moments, owner or not. */
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    const int o = radiation_isrf_moment_to_operator[m];
    radiation_apply_flux_limiter_band(limiter_state[o], limiter_scale[o],
                                      fd->isrf_moment[m].specific_flux);
  }
}

/**
 * @brief Comoving path length of the receiver-side LW/PE dust extinction
 * column, chosen by #isrf_extinction_path_mechanism.
 *
 * The arguments are the union over all mechanisms; each case uses a subset.
 * The call is inside the star-gas pair loop. The temperature_capped_jeans
 * mechanism is costly there (one temperature and a square root per pair).
 *
 * @param fb_props Properties of the feedback scheme.
 * @param p The receiving #part.
 * @param xp The receiving #xpart (tracked species, for the temperature).
 * @param r Comoving separation of the illuminating star and the receiver.
 * @param cosmo The current cosmological model.
 * @param phys_const The physical constants.
 * @param hydro_props The hydro scheme properties.
 * @param us Unit system.
 * @param cooling The cooling function properties.
 * @return Comoving extinction path length.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_extinction_path(
    const struct feedback_props *fb_props, const struct part *p,
    const struct xpart *xp, const float r, const struct cosmology *cosmo,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling) {

  /* One support radius: the longest admissible path, and the fallback for a
   * degenerate gas state. */
  const float h_gas = p->h * kernel_gamma;

  switch ((enum isrf_extinction_path_mechanism)
              fb_props->ISRF_extinction_path_mechanism) {

    case isrf_extinction_path_constant_kernel_path:
      return fb_props->ISRF_extinction_path_in_kernel_radii * h_gas;

    case isrf_extinction_path_pair_separation:
      return r;

    case isrf_extinction_path_temperature_capped_jeans: {
      const float T = cooling_get_temperature(phys_const, hydro_props, us,
                                              cosmo, cooling, p, xp);
      const float rho_phys = hydro_get_physical_density(p, cosmo);

      /* Return the cap: a clamped division is unsafe under -ffast-math. */
      if (T <= 0.f || rho_phys <= 0.f) return h_gas;

      /* c_s^2 scales with T, so the cap multiplies it by T_cap/T. */
      const float cs = hydro_get_physical_soundspeed(p, cosmo);
      const float T_ratio =
          min(1.f, fb_props->ISRF_extinction_jeans_temperature_cap_K / T);
      const float cs2_capped = cs * cs * T_ratio;

      /* lambda_J = sqrt(pi c_s^2 / (G rho)), physical, then comoving. */
      const float lambda_J_phys = sqrtf(
          (float)(M_PI * cs2_capped / (phys_const->const_newton_G * rho_phys)));
      return min(lambda_J_phys * cosmo->a_inv, h_gas);
    }
  }

  error("Unknown GEARFeedback:ISRF_extinction_path mechanism %d.",
        (int)fb_props->ISRF_extinction_path_mechanism);
  return 0.f;
}

/**
 * @brief Comoving gas column density at a gas particle: its own density times
 * the extinction path.
 *
 * @param p The #part.
 * @param extinction_path Comoving path length, from
 * #radiation_get_comoving_extinction_path.
 * @return Comoving gas column density at the particle's own location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_part(const struct part *p,
                                                  const float extinction_path) {
  return extinction_path * p->rho;
}

/**
 * @brief Dust-to-gas ratio relative to the Milky Way, following Grackle:
 * `(fgr/fgr_default) * Z/z_solar`, see
 * #RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO.
 *
 * @param Z Gas metal mass fraction.
 * @param local_dust_to_gas_ratio Resolved chemistry_data value, not the raw
 * `-1` sentinel.
 * @return Dust-to-gas ratio relative to the Milky Way.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_to_gas_ratio_relative_to_MW(float Z,
                                               float local_dust_to_gas_ratio) {
  const float D_relative_Z =
      max(Z, 0.f) / RADIATION_GRACKLE_SOLAR_METAL_FRACTION;
  return D_relative_Z * (local_dust_to_gas_ratio /
                         (float)RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO);
}

/**
 * @brief Band dust mass opacity in internal units (area/mass):
 * `sigma_d_band * D(Z) / (mu_H * m_H)`.
 *
 * @param us Unit system.
 * @param phys_const Physical constants (for the proton mass).
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Dust mass opacity, internal units.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_mass_opacity(const struct unit_system *us,
                                const struct phys_const *phys_const, float Z,
                                float sigma_d_band_cgs,
                                float local_dust_to_gas_ratio) {

  const float D_relative = radiation_get_dust_to_gas_ratio_relative_to_MW(
      Z, local_dust_to_gas_ratio);
  return sigma_d_band_cgs * D_relative /
         (RADIATION_MU_H * phys_const->const_proton_mass *
          units_cgs_conversion_factor(us, UNIT_CONV_AREA));
}

/**
 * @brief Band dust extinction factor `exp(-kappa_eff * Sigma_gas_p)`.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param Sigma_gas_p Physical gas column density, internal units.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Dust extinction factor, in (0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_extinction_factor(const struct unit_system *us,
                                     const struct phys_const *phys_const,
                                     float Z, float sigma_d_band_cgs,
                                     float Sigma_gas_p,
                                     float local_dust_to_gas_ratio) {

  const float kappa_eff = radiation_get_dust_mass_opacity(
      us, phys_const, Z, sigma_d_band_cgs, local_dust_to_gas_ratio);
  return expf(-kappa_eff * Sigma_gas_p);
}

/**
 * @brief Local linear dust absorption rate `kappa_eff(Z) * rho` in 1/length.
 * It is local: no column density is involved.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param Z Gas metal mass fraction.
 * @param rho_p Physical gas density, internal units.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Local linear dust absorption rate, internal units (1/length).
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_linear_absorption_rate(const struct unit_system *us,
                                          const struct phys_const *phys_const,
                                          float Z, float rho_p,
                                          float sigma_d_band_cgs,
                                          float local_dust_to_gas_ratio) {
  return radiation_get_dust_mass_opacity(us, phys_const, Z, sigma_d_band_cgs,
                                         local_dust_to_gas_ratio) *
         rho_p;
}

/**
 * @brief Receiver-side LW/PE dust extinction factors for a gas particle.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The receiving #part.
 * @param Z The receiving particle's own metal mass fraction.
 * @param cooling The cooling function properties (for the resolved
 * chemistry_data.local_dust_to_gas_ratio).
 * @param extinction_path Comoving extinction path length, from
 * #radiation_get_comoving_extinction_path.
 * @param extinction (return) Extinction factor of each operator, indexed by
 * #radiation_isrf_operator.
 */
__attribute__((always_inline)) INLINE void
radiation_get_part_ISRF_extinction_factors(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct part *p, float Z,
    const struct cooling_function_data *cooling, const float extinction_path,
    float extinction[ISRF_OPERATOR_COUNT]) {

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(p, extinction_path) *
      cosmo->a2_inv;
  /* Resolved value, not the raw `-1` sentinel. */
  const float local_dust_to_gas_ratio =
      (float)cooling->chemistry_data.local_dust_to_gas_ratio;

  extinction[ISRF_OPERATOR_PE] = radiation_get_dust_extinction_factor(
      us, phys_const, Z, RADIATION_SIGMA_D_PE_CGS, Sigma_gas_p,
      local_dust_to_gas_ratio);
  extinction[ISRF_OPERATOR_LW] = radiation_get_dust_extinction_factor(
      us, phys_const, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p,
      local_dust_to_gas_ratio);
}
