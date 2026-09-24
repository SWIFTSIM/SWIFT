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
 * @brief Receiver-side LW/PE dust extinction and hyperbolic
 * P1-relaxation propagation physics for GEAR.
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

/* See this global's own doxygen, radiation_isrf.h. */
int isrf_c_hyp_consistent_variable_c = 0;

/**
 * @brief First-init of a #part's LW/PE radiation-field state. Shared
 * across GEAR feedback variants: independent of the injection mechanism.
 *
 * Deliberately does NOT zero #feedback_isrf_moment_data.u: this
 * runs (space_first_init.c) after the IC file has been read into #part
 * (single_io.c/parallel_io.c/serial_io.c), and an IC may supply them via
 * the optional "PESpecificEnergy"/"LWSpecificEnergy" fields (see
 * src/feedback/GEAR_thermal/feedback_io.h) for a validation setup that
 * bypasses star injection entirely; zeroing here would silently stomp
 * that value back to 0.f. An IC that does not supply them is unaffected:
 * every #part is bzero'd before the IC read runs, so every band's u is
 * already 0.f by the time this function is reached, identical to the
 * previous unconditional assignment.
 *
 * #feedback_isrf_moment_data.u_prev is seeded from #feedback_isrf_moment_data.u
 * rather than left at 0.f, for the same reason: with `ISRF_propagation` on, the
 * engine's initial pass (before the first real step's
 * #radiation_snapshot_part_propagation call has ever run) reaches
 * #radiation_end_force_propagation directly off this first-init state, and
 * that update rebuilds `u` from `u_prev`. A 0.f seed there would wipe an
 * IC-supplied field to (near) 0 before it is ever seen. Seeding from
 * #feedback_isrf_moment_data.u reduces to today's behaviour when the IC does
 * not supply them (both already 0.f from the bzero above).
 *
 * #feedback_isrf_moment_data.specific_flux are zeroed unconditionally: no IC
 * field is proposed for them, so the "don't stomp an IC value" concern above
 * does not apply.
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
  /* -1 so an MPI foreign particle's uninitialized memory (not covered by
     the IC-read bzero above) can never read as "still illuminated". */
  fd->is_illuminated_ISRF = 0;
  fd->ISRF_illumination_end_ti = -1;
  /* Placeholder, not a real density: #part.rho is not yet computed at this
   * point (the "Densities" IC field, if supplied, is output-only; SPH
   * density is always computed from scratch by the first density loop),
   * and the engine's initial density/gradient pass calls the propagation
   * ghosts directly off this first-init state, before
   * #radiation_snapshot_part_propagation has ever run. A 0.f seed here
   * would turn `grad(u)`'s `1/rho_i` into +inf, poisoning every particle
   * with NaN before the run's first real step even begins (see this
   * field's own doxygen, feedback_struct.h). 1.0f is safe regardless of
   * this particle's real density, since `u`/`F` are also still 0 at this
   * point. */
  fd->rho_prev = 1.0f;
  fd->c_hyp = 0.f;
  fd->dt_prev = 0.f;
  fd->ISRF_reservoir_end_ti = -1;
  radiation_init_part_propagation(p);
  radiation_cache_m1_closure_part(p);
}

/**
 * @brief Snapshot #feedback_isrf_moment_data.u once per step, before the
 * density loop's h-iterations begin (see #feedback_isrf_moment_data.u_prev),
 * cache this step's per-band absorption rate, cache a stable comoving-density
 * snapshot the propagation loops need (see #feedback_part_data.rho_prev's own
 * doxygen for why), cache this step's own physical timestep
 * #feedback_part_data.dt_prev and, for #feedback_props.ISRF_c_hyp_scheme 0
 * (shipped) or 2 (fixed-fraction), #feedback_part_data.c_hyp itself (scheme 1,
 * kernel-local, defers c_hyp to radiation_end_density_propagation, once the
 * density loop's neighbour-bin maximum is known), zero every per-step
 * gradient-loop accumulator, and, for active particles only, draw down this
 * step's
 * #feedback_isrf_moment_data.u_source_rate from
 * #feedback_isrf_moment_data.u_dose_reservoir.
 *
 * Must run here, not in #radiation_init_part_propagation: this call site
 * (cell_drift.c) precedes chemistry_init_part's per-step reset of
 * smoothed_metal_mass_fraction and hydro_init_part's per-step reset of
 * #part.rho, so it is the last point where those fields still hold the
 * previous step's converged value.
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

    /* Force-loop accumulator, zeroed here for the same reason as grad_u
     * above: that loop, like the gradient loop, runs exactly once per step,
     * with no h-iteration redo. */
    moment->dissipation_u = 0.f;
  }

  /* Stable comoving density snapshot, cached unconditionally (not gated on
   * ISRF_propagation below): the gradient loop's `grad(u)` accumulation
   * (radiation_propagation_iact.h) always runs, even in injection-only
   * mode, and reads this snapshot. See this field's own doxygen
   * (feedback_struct.h) for why the density loop cannot use a live `p->rho`
   * instead, and why it must never be left at 0.f. */
  const float rho_comoving = hydro_get_comoving_density(p);
  p->feedback_data.rho_prev = rho_comoving > 0.f ? rho_comoving : 1.0f;

  if (!e->feedback_props->ISRF_propagation) {
    p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa = 0.f;
    p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa = 0.f;
    return;
  }

  const float rho_phys = hydro_get_physical_density(p, e->cosmology);
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  /* Resolved value, never the raw `-1`-sentinel
   * cooling->local_dust_to_gas_ratio. */
  const float local_dust_to_gas_ratio =
      (float)e->cooling_func->chemistry_data.local_dust_to_gas_ratio;
  p->feedback_data.isrf_operator[ISRF_OPERATOR_PE].kappa =
      radiation_get_part_linear_absorption_rate(e->internal_units, Z, rho_phys,
                                                RADIATION_SIGMA_D_PE_CGS,
                                                local_dust_to_gas_ratio);
  p->feedback_data.isrf_operator[ISRF_OPERATOR_LW].kappa =
      radiation_get_part_linear_absorption_rate(e->internal_units, Z, rho_phys,
                                                RADIATION_SIGMA_D_LW_CGS,
                                                local_dust_to_gas_ratio);

  /* This particle's own physical timestep, using its already-decided
   * integer timestep, not a new timestep-computation hook. dt_i is floored
   * at FLT_MIN so a not-yet-assigned time_bin (only possible before this
   * particle's very first real step) cannot divide by an exact zero. */
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_step = get_integer_timestep(p->time_bin);
  /* One-step lookback (`ti_current - ti_step`, not `ti_current` itself) is
   * load-bearing for the dose-reservoir drawdown below: it keeps the drain
   * rate constant across a star step's sub-steps. Using `ti_current` would
   * over-drain and empty the reservoir one sub-step early. */
  const integertime_t ti_begin =
      get_integer_time_begin(e->ti_current, p->time_bin);
  float dt_phys;
  if (with_cosmology) {
    dt_phys = (float)cosmology_get_delta_time(e->cosmology, ti_begin,
                                              ti_begin + ti_step);
  } else {
    dt_phys = (float)get_timestep(p->time_bin, e->time_base);
  }
  dt_phys = max(dt_phys, FLT_MIN);

  /* Schemes "kernel-local" and "kernel-local + variable-c" both defer
   * c_hyp entirely to radiation_end_density_propagation, once the density
   * loop's neighbour-bin maximum (dt_max(i)) is known; drift only caches
   * dt_i here (below). The other schemes decide c_hyp now, from dt_i, and
   * radiation_end_density_propagation is a no-op for them, so this branch
   * is the ENTIRE definition of c_hyp for schemes 0, 2 and 3, bit-identical
   * to the pre-comparison-branch shipped formula when scheme is 0. */
  if (e->feedback_props->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_kernel_local &&
      e->feedback_props->ISRF_c_hyp_scheme !=
          isrf_c_hyp_scheme_kernel_local_plus_variable_c) {
    const float h_phys = (float)e->cosmology->a * p->h;
    float c_hyp;
    if (e->feedback_props->ISRF_c_hyp_scheme ==
        isrf_c_hyp_scheme_fixed_fraction) {
      /* Uniform reduced light-speed candidate: every particle gets the
       * same c_hyp, independent of h_phys/dt_phys above. The receiver-side
       * CFL this removes coverage for is instead enforced by a dedicated
       * timestep term, see radiation_isrf_part_timestep(). */
      c_hyp = e->feedback_props->ISRF_c_hyp_fixed_fraction_of_c *
              (float)e->physical_constants->const_speed_light_c;
    } else {
      c_hyp = e->feedback_props->ISRF_c_hyp_margin * h_phys / dt_phys;
      c_hyp = min(c_hyp, (float)e->physical_constants->const_speed_light_c);
    }
    /* The debug pin is applied after the light-speed clamp above and is
     * not itself clamped: a pin value above c gives a superluminal
     * propagation speed on purpose, for isolating dispersion behaviour at
     * chosen values of the Courant number. Never set it above c outside
     * of that use. */
    if (e->feedback_props->ISRF_c_hyp_pin_for_debugging > 0.f)
      c_hyp = e->feedback_props->ISRF_c_hyp_pin_for_debugging;

    p->feedback_data.c_hyp = c_hyp;
  }
  p->feedback_data.dt_prev = dt_phys;

  /* Dose-reservoir drawdown, for active particles only: a cell drifted for an
   * inactive particle must not draw down a dose it will not integrate this
   * step. An inactive particle's u_*_source_rate is simply left at last
   * step's value; it is never read again before this function next runs
   * for it (as active) and overwrites it. */
  struct feedback_part_data *fd = &p->feedback_data;
  if (part_is_active(p, e) &&
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
 * @brief Radiation timestep term for the uniform reduced light-speed
 * candidate (#feedback_props.ISRF_c_hyp_scheme ==
 * #isrf_c_hyp_scheme_fixed_fraction, magnitude
 * #feedback_props.ISRF_c_hyp_fixed_fraction_of_c): every particle's c_hyp
 * is fixed at `f*c` there, independent of h/dt, so unlike the other two
 * schemes' `c_hyp_i` (each derived from a timestep) this speed carries
 * no built-in guarantee that `f*c*dt_i <= C_hyp*h_i`. This returns that
 * bound directly, `C_hyp*h_i/(f*c)`, following the same #ISRF_c_hyp_margin
 * used by the shipped formula.
 *
 * FLT_MAX (no constraint) whenever: the fixed fraction is off (0, the
 * default -- the other two schemes need no such term, since their own
 * c_hyp is already derived from a timestep); ISRF_propagation is off (no
 * flux transport,
 * so no receiver-side CFL to protect); the debug off-switch is set
 * (diagnostic-only, see that parameter's own doxygen); or this particle is
 * outside the narrow eligible set below.
 *
 * Eligible set (narrowest defensible, not "every particle"): this
 * particle's own field is live (#feedback_part_data.is_illuminated_ISRF,
 * or, since that tag lapses while #feedback_isrf_moment_data.u itself is
 * held indefinitely -- "never cleared by cooling", see that field's own
 * doxygen -- any band's u != 0), OR a neighbour inside this particle's own
 * kernel carries field this step (any band's
 * #feedback_isrf_operator_data.ngb_mean_abs_u_V > 0, the same kernel-mean the
 * negativity trigger reads, giving one kernel of margin before the front
 * itself arrives). A particle with neither can only start receiving flux
 * next step via a neighbour that is itself constrained, or via direct
 * star injection, which unconditionally calls timestep_sync_part on first
 * touch (radiation_iact.h) and so picks up this constraint on the step it
 * needs it, without waiting for its own next unforced timestep
 * recomputation.
 *
 * @param p The #part to consider.
 * @param e The #engine.
 * @return The radiation timestep bound, or FLT_MAX if none applies.
 */
float radiation_isrf_part_timestep(const struct part *restrict p,
                                   const struct engine *e) {
  const float f = e->feedback_props->ISRF_c_hyp_fixed_fraction_of_c;
  /* Gating on f alone (not e->feedback_props->ISRF_c_hyp_scheme itself) is
   * safe only because feedback_props_check_c_hyp_scheme() forces the two
   * to agree at parse time: f > 0 implies isrf_c_hyp_scheme_fixed_fraction.
   * If that pairing check is ever relaxed, this gate must switch to the
   * scheme directly. */
  if (f <= 0.f) return FLT_MAX;
  if (!e->feedback_props->ISRF_propagation) return FLT_MAX;
  if (e->feedback_props->ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging)
    return FLT_MAX;

  const struct feedback_part_data *fd = &p->feedback_data;
  const int near_field =
      fd->is_illuminated_ISRF || fd->isrf_moment[ISRF_MOMENT_PE].u != 0.f ||
      fd->isrf_moment[ISRF_MOMENT_LW].u != 0.f ||
      fd->isrf_operator[ISRF_OPERATOR_PE].ngb_mean_abs_u_V > 0.f ||
      fd->isrf_operator[ISRF_OPERATOR_LW].ngb_mean_abs_u_V > 0.f;
  if (!near_field) return FLT_MAX;

  const float h_phys = (float)e->cosmology->a * p->h;
  /* Explicit branch, not a clamp: under -ffast-math a clamp does not
   * shield a NaN/inf that a 0-numerator division could otherwise produce
   * downstream (see swift-knowledge.md's floating-point-hazards section).
   * h_phys <= 0 cannot happen for a real particle; guard it anyway rather
   * than trust a clamp to absorb it. */
  if (h_phys <= 0.f) return FLT_MAX;

  const float c_M = f * (float)e->physical_constants->const_speed_light_c;
  return e->feedback_props->ISRF_c_hyp_margin * h_phys / c_M;
}

/**
 * @brief Zero the kernel-mean per-h-iteration accumulator and the ISRF
 * density-loop neighbour-bin maximum #feedback_part_data.max_ngb_time_bin
 * (reset to this particle's own #part.time_bin, so a particle with no
 * neighbours this iteration still yields `dt_max(i) = dt_i`; read only by
 * the kernel-local scheme, but accumulated unconditionally -- one byte
 * compare per pair -- so switching #feedback_props.ISRF_c_hyp_scheme at
 * runtime needs no separate code path here). Mirrors chemistry_init_part's
 * own per-iteration reset (called from the same sites: part_init.h and the
 * ghost h-iteration redo path), so it is safe to call once or several times
 * per step. The force-loop accumulators are zeroed once per step elsewhere:
 * #feedback_isrf_moment_data.dissipation_u by
 * #radiation_snapshot_part_propagation, and
 * #feedback_isrf_moment_data.div_specific_flux by
 * #radiation_end_gradient_propagation.
 *
 * @param p The #part to reset.
 */
void radiation_init_part_propagation(struct part *p) {
  p->feedback_data.max_ngb_time_bin = p->time_bin;
  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++)
    p->feedback_data.isrf_operator[o].ngb_mean_abs_u_V = 0.f;
}

/**
 * @brief Cache this active particle's kernel-local hyperbolic propagation
 * speed #feedback_part_data.c_hyp, once the density loop's h-iteration has
 * converged and #feedback_part_data.max_ngb_time_bin therefore holds the
 * true neighbour-bin maximum for this step's kernel, then rebuild the M1
 * closure cache from it (#radiation_cache_m1_closure_part), so the
 * gradient loop that follows reads a closure built from THIS step's speed
 * rather than the stale value the drift-time call left behind.
 *
 * `c_hyp_i = min(C_hyp*h_i/dt_max(i), c)`, `dt_max(i)` the physical
 * duration of a step at #max_ngb_time_bin, computed with the same
 * get_integer_timestep/get_integer_time_begin/cosmology_get_delta_time
 * calls #radiation_snapshot_part_propagation uses for this particle's own
 * `dt_i`, just evaluated at the neighbour-maximum bin instead. Same-bin
 * case (#max_ngb_time_bin equal to #part.time_bin, i.e. every neighbour on
 * this particle's own clock): reuses #dt_prev, already this step's `dt_i`
 * from the drift, rather than a second call with the same bin, so the
 * result is bit-identical to the shipped per-particle scheme there,
 * independent of codegen. `dt_max` is floored at FLT_MIN in both branches:
 * before this particle's first drift has ever run (the initial,
 * pre-any-step gradient pass), #dt_prev is still its first-init 0.f, which
 * would otherwise divide by an exact zero.
 *
 * No-op unless #feedback_props.ISRF_c_hyp_scheme is
 * #isrf_c_hyp_scheme_kernel_local or
 * #isrf_c_hyp_scheme_kernel_local_plus_variable_c: for the other schemes,
 * drift-time #radiation_snapshot_part_propagation already decided #c_hyp
 * (and #feedback_reset_part already cached the M1 closure built from it),
 * and this function must leave that alone, bit-identical to the
 * pre-comparison-branch behaviour for the shipped scheme. Under
 * #isrf_c_hyp_consistent_variable_c the rebuilt M1 closure's own `c_M` is
 * pinned to 1 regardless of `c_hyp` (see #radiation_cache_m1_closure_part's
 * own doxygen), so for scheme 4 this function still updates `c_hyp` itself
 * (read by the pairwise operators' receiver-side multiply), even though the
 * closure rebuild that follows does not change value because of it.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_density_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;
  if (e->feedback_props->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_kernel_local &&
      e->feedback_props->ISRF_c_hyp_scheme !=
          isrf_c_hyp_scheme_kernel_local_plus_variable_c)
    return;

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
  dt_max = max(dt_max, FLT_MIN);

  const float h_phys = (float)e->cosmology->a * p->h;
  float c_hyp = e->feedback_props->ISRF_c_hyp_margin * h_phys / dt_max;
  c_hyp = min(c_hyp, (float)e->physical_constants->const_speed_light_c);
  /* The debug pin is applied after the light-speed clamp above and is not
   * itself clamped: a pin value above c gives a superluminal propagation
   * speed on purpose, for isolating dispersion behaviour at chosen values
   * of the Courant number. Never set it above c outside of that use. */
  if (e->feedback_props->ISRF_c_hyp_pin_for_debugging > 0.f)
    c_hyp = e->feedback_props->ISRF_c_hyp_pin_for_debugging;

  fd->c_hyp = c_hyp;

  radiation_cache_m1_closure_part(p);
}

/**
 * @brief The exact-relaxation integrating factor `phi(a) = (1 -
 * exp(-a))/a`, shared by both the `u`- and `F`-updates below:
 * `phi -> 1` as `a -> 0` (no absorption within the step: the update
 * reduces to plain forward transport) and `phi -> 1/a` as `a -> inf` (the
 * stiff limit: the update relaxes exactly to its local target every
 * step). Implemented via `-expm1(-a)/a`, which is well-conditioned at
 * both ends except very close to `a=0`, where a short Taylor series is
 * used instead to avoid a `0/0` cancellation.
 *
 * @param a Dimensionless relaxation depth for this step: `c_hyp*(kappa +
 * H/c)*dt` at the two ghost sites, which relax an absorption rate and the
 * cosmological redshift rate together, both dilated by the same `c_hyp/c`
 * factor (see #radiation_end_force_propagation), or the injection-site
 * equivalent, `c_hyp*kappa*Delta_t_star`. Always `>= 0`.
 * @return phi(a).
 */
float radiation_relaxation_phi_factor(float a) {
  if (a < 1e-6f) return 1.0f - 0.5f * a + (1.0f / 6.0f) * a * a;
  return -expm1f(-a) / a;
}

/**
 * @brief The three mutually exclusive outcomes of the M1 flux limiter,
 * split out so the operator that shares a limiter decision across several
 * moments (#radiation_end_gradient_propagation) applies the exact same
 * branch to each of them, rather than a value that happens to be
 * algebraically equivalent.
 */
enum radiation_isrf_flux_limiter_state {
  ISRF_LIMITER_ZERO, /*!< `u <= 0`: every flux component is zeroed. */
  ISRF_LIMITER_SKIP, /*!< `|F|^2 <= 0`: left untouched, no multiply. */
  ISRF_LIMITER_SCALE /*!< Scaled by #scale below. */
};

/**
 * @brief M1 flux limiter decision for one particle, one band:
 * `F <- F*min(1, c_M*u/|F|)` for `u > 0`, `F <- 0` for `u <= 0`. Enforces
 * the reduced-flux closure's guarantee (the interior field is `|F|/c_M`,
 * not more) against whatever `u` this call is given. Split from the apply
 * half (#radiation_apply_flux_limiter_band) so a scale computed once, from
 * one band, can be applied identically to another band that shares this
 * band's operator.
 *
 * Guarded rather than relying on algebraic cancellation: `F = 0` under
 * `u > 0` needs no division at all (scaling the zero vector is still
 * zero), so that case is skipped outright instead of computing
 * `c_M*u/|F|` unguarded. `F.F` and the limiter ratio are formed in double:
 * in float32, `F.F` underflows to zero once `|F| < sqrt(FLT_MIN) ~ 1.1e-19`
 * (internal units), which would skip the limiter for a nonzero flux.
 *
 * @param u This band's specific field `u^n`.
 * @param c_M This particle's own #feedback_part_data.c_hyp, or 1 under
 * #isrf_c_hyp_consistent_variable_c (`F` is then already the reduced flux
 * `Ft = F_true/c_hyp`, whose own bound is `|Ft| <= u`).
 * @param F This particle's tracked flux (this band), unmodified.
 * @param scale (return) The multiplier to apply, valid only when the
 * return value is #ISRF_LIMITER_SCALE.
 * @return Which of the three outcomes applies.
 */
__attribute__((
    always_inline)) INLINE static enum radiation_isrf_flux_limiter_state
radiation_compute_flux_limiter_scale_band(float u, float c_M, const float F[3],
                                          float *scale) {

  if (u <= 0.f) return ISRF_LIMITER_ZERO;

  const double F2 = (double)F[0] * (double)F[0] + (double)F[1] * (double)F[1] +
                    (double)F[2] * (double)F[2];
  if (F2 <= 0.) return ISRF_LIMITER_SKIP;

  *scale = (float)min(1., (double)c_M * (double)u / sqrt(F2));
  return ISRF_LIMITER_SCALE;
}

/**
 * @brief Apply a flux-limiter decision
 * (#radiation_compute_flux_limiter_scale_band) to one band's flux. A `switch`
 * on the three states rather than a single scale value, so the
 * `ISRF_LIMITER_SKIP` case takes no multiply at all, exactly like the unsplit
 * function this replaces: "scale by 1.0" is the same value but not necessarily
 * the same instruction sequence under
 * `-ffast-math`.
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
 * @brief Exact-relaxation update of #feedback_isrf_moment_data.u, from the
 * `div(F)` and negativity-triggered dissipation accumulators
 * radiation_propagation_iact.h filled during the force loop, and this step's
 * own #feedback_isrf_moment_data.u_source_rate (drawn down from the dose
 * reservoir by #radiation_snapshot_part_propagation).
 *
 * `u^{n+1} = e*(u_prev + dt*phi*diss) + dt*phi*((c_hyp/c)*source_rate -
 * div_F)`, with `e = exp(-a)`, `phi = (1-e)/a`, `a = c_hyp*(kappa + H/c)*dt`
 * (#radiation_relaxation_phi_factor). Without the dissipation this is the
 * exact solution, UNDER THE SAME CHANGE OF VARIABLE the reduced-speed method
 * applies to every other rate (see #isrf_c_hyp_consistent_variable_c's own
 * doxygen, radiation_propagation_iact.h), of `du/dt = -u/tau - H*u +
 * (c_hyp/c)*source_rate - div(F)` over one step with `source_rate` and
 * `div(F)` frozen, `tau = 1/(c_hyp*kappa)`. The Hubble term is dilated by the
 * SAME `c_hyp/c` factor as the absorption and injection terms: the
 * reduced-speed-of-light method is only correct if every rate in the
 * equation carries that factor, and `H`, unlike `kappa`, does not itself
 * scale with `c_hyp` (it is a property of the expanding background, not of
 * the radiation transport), so it must be dilated explicitly here rather
 * than picking it up "for free" the way `c_hyp*kappa` does. Leaving `H`
 * undilated (as this file did until this fix) makes the fixed point of the
 * homogeneous (`div_F = 0`) equation `u* = (source_rate/c)/(kappa +
 * H/c_hyp)` instead of the true-speed `u*_true = (source_rate/c)/(kappa +
 * H/c)`: since `c_hyp << c`, `H/c_hyp >> H/c`, suppressing `u*` by the
 * factor `x/(1+x)`, `x = c_hyp*kappa/H`, a spurious sink strongest exactly
 * where `c_hyp*kappa` is smallest relative to `H` (low density, low
 * metallicity, i.e. an early ultra-faint-dwarf's ISM). With `H` dilated,
 * `c_hyp` cancels out of the fixed point identically to every other term,
 * for any per-particle speed field, restoring `u*_true`. `div_F` is the
 * divergence of this step's already-relaxed
 * flux `F^{n+1}` (#radiation_end_gradient_propagation): the flux is advanced
 * first, from `grad(u^n)`, then `u` from the new flux, a staggered order with
 * the same linear stability as advancing `u` first.
 *
 * The dissipation, built from `u^n`, is decayed together with `u_prev`
 * rather than added undamped: that keeps the per-step amplification matrix
 * equal to the one of a dissipation correction applied to an already-relaxed
 * `u`, whose stability margin under the parameter guard is the larger of the
 * two placements (theory/GEAR/Radiation/verify_isrf_dissipation.py, Part I).
 * The discrete steady state is therefore
 * `(1-e)*u = dt*phi*((c_hyp/c)*source_rate - div_F + e*diss)`.
 *
 * `-H*u` is the cosmological expansion term, ONE power of the Hubble rate,
 * not three: `u` is MASS-SPECIFIC (energy per unit gas mass), so the volume
 * dilution is already carried by the physical gas density it is measured
 * against. Writing `E` for the physical volumetric band energy density,
 * `E ~ a^-4` under pure expansion (`a^-3` volume, `a^-1` redshift) while
 * `rho ~ a^-3`, so `u = E/rho` loses exactly the redshift residual:
 * `du/dt = (dE/dt)/rho - u*(drho/dt)/rho = -4H*u + 3H*u = -H*u`. This is the
 * TRUE-SPEED rate; what this update actually relaxes is its `c_hyp/c`-dilated
 * counterpart `-(c_hyp/c)*H*u`, for the fixed-point reason given above. The
 * same single power, and the same dilation, applies to the specific flux
 * (#radiation_end_gradient_propagation), whose volumetric counterpart
 * free-streams and therefore dilutes like `E`.
 *
 * It is folded into the relaxation depth rather than added as a separate
 * explicit decrement because it is a linear decay of the SAME state variable
 * the absorption term relaxes: `exp(-(1/tau + (c_hyp/c)*H)*dt)` is then exact
 * for the homogeneous problem at any `H*dt`, and cannot drive `u` negative
 * the way an explicit `-(c_hyp/c)*H*dt*u` can at high redshift with a long
 * step. A consequence of the dilation, not a defect: with no absorption at
 * all (`kappa = 0`, e.g. the `ISRFCosmology` `free_field` fixture), the
 * pure-expansion transient itself now decays at the SLOWED rate `(c_hyp/c)*H`
 * rather than the true `H`, exactly like every other transient the reduced
 * speed of light slows down; `examples/SubgridTests/StellarFeedback/ISRF/
 * ISRFCosmology/isrf_cosmology_check.py`'s `free_field` reference is
 * re-derived for this. Ungated: SWIFT
 * sets `cosmo->H = 0` for a non-cosmological run (`cosmology_init_no_cosmo`),
 * so the term vanishes there by construction (multiplying it by `c_hyp/c`
 * first does not change this: `(c_hyp/c)*0 = 0` exactly), exactly as for
 * `hydro.h`'s own `div_v + hydro_dimension*cosmo->H`. No gate on the spectrum
 * shape, unlike `src/rt/GEAR/rt.h`'s own redshift term: photons also
 * redshift ACROSS these two narrow band edges, a loss `-H*u` does not model,
 * so `-H*u` is a lower bound on the true band loss rather than an
 * overestimate to be suppressed. `c_hyp` here plays the role of the M1 reduced
 * light speed `c_M`: the `c_M/c` rescale (replacing the old,
 * P1-Yukawa-tuned `3*c_hyp/c`) is applied exclusively
 * here; injection (`radiation_iact.h`) deposits the raw, unrescaled dose.
 *
 * Runs in the `end_force` task, after the force loop and before cooling
 * (engine_maketasks.c). The negativity trigger that set this step's
 * dissipation coefficient (the extra ghost) read `u^n`, so an undershoot
 * created by this update is seen by the next step's trigger and corrected
 * by the next step's call; cooling reads it uncorrected for one step, through
 * its own non-negative clamp.
 *
 * Idempotent: `u` is rebuilt from the stable `u_prev` snapshot and this
 * step's accumulators, never incremented, so a repeated call gives the same
 * state. The M1 flux limiter is not applied here: the extra ghost already
 * limited this step's flux against `u^n`, the `u` its closure was built from.
 * The debug-only energy-ledger counters below (#feedback_isrf_moment_data.
 * cumulative_injected/cumulative_absorbed) are the one exception: they ARE
 * incremented, relying on the task graph calling this exactly once per
 * active particle per step (no h-iteration-style redo exists for the force
 * ghost, unlike the density loop).
 *
 * Reads #dt_prev (cached earlier this step by
 * #radiation_snapshot_part_propagation), #c_hyp (cached later, once the
 * density loop's neighbour-bin maximum is known, by
 * #radiation_end_density_propagation) and #feedback_isrf_operator_data.kappa,
 * and deliberately takes no `dt` of its own: the call site computes its local
 * `dt` from a different timestep-begin convention (`ti_current - 1`), and using
 * it here would make this update's `dt*phi` inconsistent with the flux
 * update's. The thin `(p, e)` signature exists to make that mistake
 * structurally impossible. No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_force_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->ISRF_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  const float dt = fd->dt_prev;
  const float c_hyp = fd->c_hyp;
  const float rescale =
      c_hyp / (float)e->physical_constants->const_speed_light_c;
  const float H = (float)e->cosmology->H;
  /* Dilated by the same c_hyp/c factor as the absorption term: see this
   * function's own doxygen. Bit-identical to the plain H when H = 0.f
   * (SWIFT's non-cosmological cosmology_init_no_cosmo sets cosmo->H = 0):
   * rescale * 0.f is exactly 0.f for any finite rescale, no rounding. */
  const float H_dilated = rescale * H;

#ifdef SWIFT_DEBUG_CHECKS
  if (fd->u_min_snapshot_index != e->snapshot_output_count) {
    for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
      fd->isrf_moment[m].u_min_since_snapshot = 0.f;
    fd->u_min_snapshot_index = e->snapshot_output_count;
  }
#endif

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    const struct feedback_isrf_operator_data *op =
        &fd->isrf_operator[radiation_isrf_moment_to_operator[m]];
    /* H_dilated is this particle's own scalar, not looked up per operator:
     * every moment applies its own residual against it here, so two
     * moments sharing one operator can still carry different residuals. */
    const float a = (c_hyp * op->kappa + H_dilated) * dt;
    const float decay = expf(-a);
    const float phi = radiation_relaxation_phi_factor(a);

#ifdef SWIFT_DEBUG_CHECKS
    /* Energy-ledger accumulation, read from this step's own inputs before
     * #u below overwrites #u_prev's role: the raw (unrelaxed) dose this
     * step attempted, and the exact-relaxation update's own split of
     * `u_prev` and the frozen source/transport terms into the fraction
     * that decayed/never-arrived this step, using the SAME `decay`/`phi`
     * #u's update uses. */
    moment->cumulative_injected += dt * rescale * moment->u_source_rate;
    moment->cumulative_absorbed +=
        (moment->u_prev + dt * phi * moment->dissipation_u) * (1.f - decay) +
        (rescale * moment->u_source_rate - moment->div_specific_flux) * dt *
            (1.f - phi);
#endif

    moment->u =
        decay * (moment->u_prev + dt * phi * moment->dissipation_u) +
        dt * phi *
            (rescale * moment->u_source_rate - moment->div_specific_flux);

#ifdef SWIFT_DEBUG_CHECKS
    if (moment->u < moment->u_min_since_snapshot)
      moment->u_min_since_snapshot = moment->u;
#endif
  }
}

/**
 * @brief Return this step's dose-reservoir drawdown to the reservoir for a
 * #part whose density h-iteration gives up with no neighbours found
 * (`runner_ghost.c`'s `has_no_neighbours` give-up path). Such a particle still
 * reaches #radiation_end_force_propagation, so the drawn rate would otherwise
 * be applied on a step whose kernel sums are meaningless; the dose is kept for
 * a later step instead. Restores it into
 * #feedback_isrf_moment_data.u_dose_reservoir exactly (the drawn amount is
 * `source_rate*dt_prev`) and zeroes the rates, so that update then carries
 * only the decay of `u_prev` and whatever the force loop accumulated.
 * No-op when propagation is off.
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
 * @brief One band's negativity-triggered artificial-dissipation coefficient
 * update: raised instantly to a negativity-triggered target, or decayed
 * toward it otherwise. Reads `u_V = rho_prev*u^n`, the state the previous
 * step's #radiation_end_force_propagation produced, against the kernel mean
 * of the neighbours' `u^n`: both at the same time level. The coefficient this
 * returns is consumed by THIS step's force loop, which dissipates `u^n`, so an
 * undershoot is corrected one step after it appears.
 *
 * @param u_V This band's volumetric field, rho_prev*u^n.
 * @param ngb_mean_abs_u_V This band's kernel-mean |rho_prev*u_prev| scratch
 * accumulator (radiation_propagation_iact.h).
 * @param alpha_prev This band's
 * #feedback_isrf_operator_data.dissipation_alpha_trigger from the previous
 * step, i.e. the trigger's own previous output, never the floor-combined
 * coefficient.
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

  /* max(ngb_mean_abs_u_V, -u_V) >= -u_V > 0 in this branch, so the
   * division below is always well-defined. */
  const float eps = (u_V < 0.f) ? -u_V / max(ngb_mean_abs_u_V, -u_V) : 0.f;
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
 * @brief The `h/lambda`-gated floor under
 * #radiation_update_dissipation_alpha_band's trigger: the trigger
 * fires only on negativity and is exactly zero on the positive delta-shell
 * front of an optically-thin P1 pulse, so a purely reactive coefficient
 * cannot damp the resulting dispersive wake there. This floor supplies
 * dissipation the trigger structurally cannot, rolling off as
 * `(eps_lambda/(h*kappa))^4` once `h/lambda` exceeds #ISRF_dissipation_
 * floor_h_over_lambda, which bounds its steady-state cost by construction.
 *
 * The quartic roll-off separates two regimes the quadratic one could not:
 * a moderately-resolved front (`h/lambda ~ 0.25`) still needs most of the
 * floor, while an optically-thick region (`h/lambda ~ 6-10`) needs it
 * essentially absent, since there the floor buys no negativity protection
 * and costs pure accuracy. A quadratic tail leaves a percent-level
 * coefficient in the thick regime; the quartic one leaves ~1e-5.
 *
 * This is the floor's only gating: the force loop applies it to a pair
 * unconditionally, as `max(trigger_i, trigger_j, floor_i, floor_j)`.
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
 * @brief Flux-relaxation residual gate on the floor's aim: a particle whose
 * flux is already in Fickian balance with this step's own gradient (`F ~=
 * -C*grad_u`, `C = c_hyp/(kappa+H/c)` the fixed point of
 * #radiation_end_gradient_propagation's own UNLIMITED flux-update
 * recurrence, i.e. before #radiation_apply_flux_limiter_band clamps it) is
 * at the discrete steady state the floor's cost formula assumes; a
 * particle on a genuine front, or with `tau = 1/(c_hyp*(kappa+H/c)) >> dt`
 * so the flux has not relaxed yet, is not. A particle whose flux is instead
 * pinned by the M1 limiter (`|F| = c_M*u`, the free-streaming branch)
 * generally never reaches that fixed point either, so `R` stays finite
 * there too: a conservative false positive that keeps part of the floor
 * where the limiter is active, never removes protection where a front is
 * present. `R` measures the mismatch (0 at the fixed point, ~1 away from
 * it); the floor's aim is multiplied by `min(1, (R/eps_R)^2)`, so it can
 * only ever be lowered, never raised: `s=1` whenever exactly one of `F`,
 * `grad_u` is zero (`R=1`), so a fresh front or a limiter-zeroed flux
 * keeps the full floor, provided the relaxation weight `w = kappa +
 * H/c` is nonzero (see below). The exception is a quiescent particle
 * with both `F` and `grad_u` zero: that is trivially at the fixed point,
 * so `R=0` and `s=0` there instead. `eps_R = 0` disables the gate
 * (returns 1 identically); `c_hyp <= 0` likewise (the relaxation has no
 * timescale to be settled against). `w <= 0` (`kappa=0` and `H=0`) also
 * returns 1 unconditionally, rather than falling through to the `R`
 * formula below: at `w=0`, `F` drops out of both that formula's numerator
 * and denominator, which would otherwise break the `s=1` guarantee above
 * whenever `F` alone is nonzero. Rescaled by `(kappa+H/c)` relative to
 * the `|F+C*grad_u|` form (the two are algebraically identical; this one
 * avoids computing `C` as its own value, which can overflow float32 at
 * near-primordial `kappa`). `w` divides `H` by the TRUE speed of light
 * `c`, not `c_hyp`: it is `a/(c_hyp*dt)` for this function's own fixed-point
 * `a = c_hyp*(kappa+H/c)*dt` (#radiation_end_gradient_propagation), so
 * `c_hyp` cancels out of `w` itself, unlike `a`. `R` and `(R/eps_R)^2` are
 * formed in double so that the squared denominator cannot underflow to zero
 * under this build's fast-math folding of the ratio and its square into one
 * division.
 *
 * @param F This band's #feedback_isrf_moment_data.specific_flux, from BEFORE
 * this step's own update (the flux `u^n` was produced from, as left by the
 * previous step's #radiation_apply_flux_limiter_band, already post-limiter).
 * @param grad_u This band's #feedback_isrf_moment_data.grad_u accumulator.
 * @param c_hyp The particle's own #c_hyp.
 * @param kappa This band's #feedback_isrf_operator_data.kappa.
 * @param H The Hubble rate, #cosmology.H.
 * @param c The TRUE speed of light, #phys_const.const_speed_light_c (not
 * `c_hyp`: see this function's own doxygen for why `w` uses the true speed).
 * @param eps_R #feedback_props.ISRF_dissipation_floor_relaxation_residual.
 * @return The floor-aim multiplier `s`, in `[0, 1]`.
 */
__attribute__((always_inline)) INLINE static float
radiation_dissipation_floor_relaxation_gate(const float F[3],
                                            const float grad_u[3], float c_hyp,
                                            float kappa, float H, float c,
                                            float eps_R) {

  if (eps_R <= 0.f) return 1.f;
  if (c_hyp <= 0.f) return 1.f;

  /* Rescaled by (kappa + H/c) relative to the doxygen's |F + C*grad_u|
   * form: algebraically identical (this factor cancels top and bottom),
   * but every term here stays O(1)-to-O(1e10) on production fixtures,
   * where computing C = c_hyp/(kappa+H/c) as its own value first
   * can overflow float32 at near-primordial kappa. Divides by the TRUE
   * speed c, not c_hyp: see this function's own doxygen. */
  const float w = kappa + H / c;

  /* No relaxation timescale to settle against (kappa = 0 and H = 0): keep
   * the floor at full strength. Also avoids F dropping out of both the R
   * numerator and denominator below, which would otherwise make R = 1
   * regardless of grad_u and break the s=1 guarantee documented above. */
  if (w <= 0.f) return 1.f;

  /* The norms are formed in double: in float32 each squared component
   * underflows to zero once its magnitude drops below sqrt(FLT_MIN) ~
   * 1.1e-19 (internal units), which would report a small nonzero flux or
   * gradient as exactly zero and send a front to the quiescent branch. */
  const double w_d = w;
  const double c_d = c_hyp;
  const double F_d[3] = {F[0], F[1], F[2]};
  const double G_d[3] = {grad_u[0], grad_u[1], grad_u[2]};
  const double wx = w_d * F_d[0] + c_d * G_d[0];
  const double wy = w_d * F_d[1] + c_d * G_d[1];
  const double wz = w_d * F_d[2] + c_d * G_d[2];
  const double num = sqrt(wx * wx + wy * wy + wz * wz);
  const double F_norm =
      sqrt(F_d[0] * F_d[0] + F_d[1] * F_d[1] + F_d[2] * F_d[2]);
  const double G_norm =
      sqrt(G_d[0] * G_d[0] + G_d[1] * G_d[1] + G_d[2] * G_d[2]);
  const double den = w_d * F_norm + c_d * G_norm;

  /* Quiescent particle (`F` and `grad_u` both zero): trivially at the fixed
   * point, so `R = 0`. Branched rather than kept finite by an epsilon added
   * to the denominator: `R/eps_R` is squared below, and the optimizer folds
   * that into `num^2/(den/eps_R)^2`, where a denominator epsilon small
   * enough not to perturb a real `den` underflows to zero once squared,
   * turning this case into `0/0`. */
  if (den <= 0.) return 0.f;

  /* R is formed and squared in double: under -ffast-math this expression
   * gets folded into a single division by (den*eps_R)^2, which underflows
   * float32 to zero for small enough den and yields an unfiltered NaN;
   * double's exponent range keeps that denominator representable. */
  const double R = num / den;
  const double ratio2 = (R / eps_R) * (R / eps_R);
  return (float)min(1.0, ratio2);
}

/**
 * @brief Exact-relaxation update of #feedback_isrf_moment_data.specific_flux,
 * from the `grad(u)` accumulators radiation_propagation_iact.h filled during
 * the gradient loop from `u^n`, followed by the M1 flux limiter
 * (#radiation_compute_flux_limiter_scale_band /
 * #radiation_apply_flux_limiter_band) against `u^n`, the same `u` the
 * gradient loop's closure tensor was built from. The limiter's scale is
 * decided once per operator, from its one owning moment
 * (#radiation_isrf_operator_owner), and applied identically to every moment
 * that shares that operator, so two such moments cannot be limited by
 * different factors. Runs once per step in the extra ghost, never re-run:
 * the gradient loop itself only runs once per step. Also updates the two
 * negativity-triggered dissipation components,
 * #feedback_isrf_operator_data.dissipation_alpha_trigger
 * (see #radiation_update_dissipation_alpha_band) and
 * #feedback_isrf_operator_data.dissipation_alpha_floor (see
 * #radiation_dissipation_alpha_floor_band), once per step and separately,
 * for THIS step's force loop to combine, and zeroes
 * #feedback_isrf_moment_data.div_specific_flux for that loop to accumulate
 * `div(F^{n+1})` into.
 *
 * `F_new = e*F - c_hyp^2*dt*phi*grad(u)`: the exact solution of
 * `dF/dt = -F/tau - H*F - (D/tau)*grad(u)` over one step with the source
 * term frozen at this step's value, `D/tau = c_hyp^2`. `-H*F` is the flux
 * counterpart of the `-H*u` derived at #radiation_end_force_propagation,
 * one power of the Hubble rate for the same mass-specific reason, dilated by
 * the same `c_hyp/c` factor for the same fixed-point reason, and carried the
 * same way, inside the relaxation depth `a = c_hyp*(kappa + H/c)*dt`.
 * No-op when propagation is off.
 *
 * Under #isrf_c_hyp_consistent_variable_c, #feedback_isrf_moment_data.
 * specific_flux stores the reduced flux `Ft = F_true/c_hyp` instead of
 * `F_true` (see radiation_propagation_iact.h's file header): substituting
 * `F = c_hyp*Ft` into the recurrence above and dividing through by the
 * (this-step-constant) `c_hyp` removes exactly one power of it, giving
 * `Ft_new = e*Ft - c_hyp*dt*phi*grad(u)`, and the M1 limiter's own bound
 * becomes `|Ft| <= u` (#radiation_compute_flux_limiter_scale_band with
 * `c_M = 1`, matching #radiation_cache_m1_closure_part's own selection for
 * the same scheme).
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
  /* Under the consistent-variable-c scheme #specific_flux is already the
   * reduced flux Ft = F/c_hyp, whose own limiter bound is |Ft| <= u: see
   * this function's own doxygen. */
  const float c_M = isrf_c_hyp_consistent_variable_c ? 1.f : c_hyp;
  const float H = (float)e->cosmology->H;
  const float c = (float)e->physical_constants->const_speed_light_c;
  /* Dilated by the same c_hyp/c factor as the absorption term: see this
   * function's own doxygen and #radiation_end_force_propagation's. Bit-
   * identical to the plain H when H = 0.f (non-cosmological runs). */
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

  /* Pre-update flux, kept per OPERATOR (not per moment): the floor's
   * relaxation-residual gate below is evaluated once per operator, from
   * its one owning moment, so only the owner's snapshot is ever read back.
   * Written from the moment loop below, guarded to the owner so a moment
   * that merely shares an operator cannot overwrite its owner's value. */
  float F_old_stash[ISRF_OPERATOR_COUNT][3];

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    const int o = radiation_isrf_moment_to_operator[m];
    const struct feedback_isrf_operator_data *op = &fd->isrf_operator[o];

    /* H_dilated is this particle's own scalar, not looked up per operator:
     * see #radiation_end_force_propagation's matching comment. */
    const float a = (c_hyp * op->kappa + H_dilated) * dt;
    const float decay = expf(-a);
    const float phi = radiation_relaxation_phi_factor(a);
    /* One power of c_hyp under the consistent-variable-c scheme: see this
     * function's own doxygen for the substitution F = c_hyp*Ft. */
    const float coeff = isrf_c_hyp_consistent_variable_c
                            ? c_hyp * dt * phi
                            : c_hyp * c_hyp * dt * phi;

    /* Snapshot for the floor's relaxation-residual gate below: `u^n` was
     * produced from THIS flux, not the one about to be computed. */
    if (m == (int)radiation_isrf_operator_owner[o]) {
      F_old_stash[o][0] = moment->specific_flux[0];
      F_old_stash[o][1] = moment->specific_flux[1];
      F_old_stash[o][2] = moment->specific_flux[2];
    }

    for (int k = 0; k < 3; k++) {
      moment->specific_flux[k] =
          decay * moment->specific_flux[k] - coeff * moment->grad_u[k];
    }

    /* Zeroed here rather than in the drift snapshot, unlike dissipation_u:
     * every drift, including the one before a snapshot dump, would otherwise
     * blank the PE/LWSpecificFluxDivergences output field. */
    moment->div_specific_flux = 0.f;
  }

  /* The flux-limiter decision and the two dissipation coefficients are
   * OPERATOR state: computed once per operator, from its one owning
   * moment, then applied to every moment that shares it (the moment loop
   * below). An operator-bounded loop can only run #ISRF_OPERATOR_COUNT
   * times per particle, so this cannot turn into a last-writer-wins or a
   * self-referential re-entry the way a moment-bounded loop over the same
   * writes would once a second moment shares an operator. */
  enum radiation_isrf_flux_limiter_state limiter_state[ISRF_OPERATOR_COUNT];
  float limiter_scale[ISRF_OPERATOR_COUNT];

  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++) {
    const int m = radiation_isrf_operator_owner[o];
    const struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    struct feedback_isrf_operator_data *op = &fd->isrf_operator[o];

    /* The just-updated, still-unlimited flux: the moment loop above has
     * already run to completion, so #radiation_apply_flux_limiter_band has
     * not yet touched this moment's flux. */
    limiter_state[o] = radiation_compute_flux_limiter_scale_band(
        moment->u, c_M, moment->specific_flux, &limiter_scale[o]);

    const float u_V = fd->rho_prev * moment->u;

    if (alpha_pin > 0.f) {
      /* Bypass the trigger entirely: every particle's coefficient is held at
       * the pinned value (see this parameter's own doxygen,
       * feedback_properties.h). The pinned value goes into the trigger
       * component and the floor component is zeroed, so that the pin stays
       * a spatially uniform coefficient: routing it through the floor would
       * subject it to the floor's own `h/lambda` roll-off and make it a
       * different quantity from the one the A/B reference arms measure. */
      op->dissipation_alpha_trigger = alpha_pin;
      op->dissipation_alpha_floor = 0.f;
    } else {
      /* The trigger's decay memory is its OWN previous value, not the
       * previous combined coefficient: a high floor must not hold up the
       * trigger's decay tail. */
      op->dissipation_alpha_trigger = radiation_update_dissipation_alpha_band(
          u_V, op->ngb_mean_abs_u_V, op->dissipation_alpha_trigger, alpha_max,
          eps_1, c_hyp, op->kappa, dt, h_phys);

      /* Stored separately from the trigger rather than combined here, so
       * the two mechanisms stay separately readable; the force loop
       * combines them. The relaxation-residual gate only ever lowers this
       * value (`s <= 1`), computed from the incoming flux snapshotted
       * above. */
      const float s = radiation_dissipation_floor_relaxation_gate(
          F_old_stash[o], moment->grad_u, c_hyp, op->kappa, H, c, eps_R);
      op->dissipation_alpha_floor =
          s * radiation_dissipation_alpha_floor_band(op->kappa, h_phys,
                                                     alpha_floor, eps_lambda);
    }
  }

  /* Apply each operator's shared limiter decision to every moment that
   * reads it, including a non-owner moment: the decision was computed
   * once above, from the owner alone. */
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    const int o = radiation_isrf_moment_to_operator[m];
    radiation_apply_flux_limiter_band(limiter_state[o], limiter_scale[o],
                                      fd->isrf_moment[m].specific_flux);
  }
}

/**
 * Comoving path length of the receiver-side LW/PE dust extinction column.
 * This is the ONLY place a path is chosen: every mechanism of
 * #isrf_extinction_path_mechanism is one case here, and the column itself
 * is formed once, in #radiation_get_comoving_gas_column_density_at_part,
 * with no knowledge of which mechanism produced the length.
 *
 * The argument list is the union over every mechanism, so each case uses a
 * subset and ignores the rest. Adding a mechanism is one enum value, one
 * case here, one string in feedback_props_init() and one entry in
 * examples/parameter_example.yml, with no call site re-cut.
 *
 * The call stays inside the star-gas pair loop for all mechanisms:
 * #isrf_extinction_path_pair_separation depends on the pair, so hoisting it
 * would fork the code path. For that mechanism and for
 * #isrf_extinction_path_constant_kernel_path the cost is a handful of flops
 * per neighbour. #isrf_extinction_path_temperature_capped_jeans is far more
 * expensive: it calls cooling_get_temperature and a square root on every
 * star-gas pair, and recomputes the same value for each star illuminating
 * the same particle.
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

  /* One support radius: the longest path the geometry admits, since the
   * illuminating star sits inside the receiver's own kernel. Also the
   * fallback for a degenerate gas state below. */
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

      /* Return the cap outright rather than clamping a division: a guarded
       * denominator is not safe under -ffast-math, where the reciprocal is
       * reassociated and can underflow before the clamp sees it. */
      if (T <= 0.f || rho_phys <= 0.f) return h_gas;

      /* c_s^2 = gamma k_B T / (mu m_H), so capping the temperature scales
       * the squared sound speed by T_cap/T with the scheme's own mu, which
       * never has to be recovered here. */
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
 * Comoving gas column density at a gas particle's own location: the
 * receiver-side analogue of the star-side Sobolev column
 * (#radiation_get_comoving_gas_column_density_at_star), used for LW/PE
 * extinction: the receiver's own local density times the extinction path
 * (no resolved density gradient on the gas side).
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
 * Dust-to-gas ratio relative to the Milky Way, matching Grackle's own
 * default convention (dust2gas = fgr * Z/z_solar): D(Z)/D(Z_sun) =
 * (fgr/fgr_default) * Z/z_solar. The fgr factor is 1 only when a run
 * leaves GrackleCooling:local_dust_to_gas_ratio at Grackle's own
 * compiled default; see #RADIATION_GRACKLE_DEFAULT_DUST_TO_GAS_RATIO.
 *
 * @param Z Gas metal mass fraction.
 * @param local_dust_to_gas_ratio Resolved
 * chemistry_data.local_dust_to_gas_ratio (never the raw `-1`-sentinel
 * cooling->local_dust_to_gas_ratio field).
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
 * Band-specific dust mass opacity (area/mass, internal units): kappa_eff
 * = sigma_d_band * D(Z) / (mu_H * m_H).
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Dust mass opacity, internal units.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_mass_opacity(const struct unit_system *us, float Z,
                                float sigma_d_band_cgs,
                                float local_dust_to_gas_ratio) {

  const float D_relative = radiation_get_dust_to_gas_ratio_relative_to_MW(
      Z, local_dust_to_gas_ratio);
  const float kappa_eff_cgs = sigma_d_band_cgs * D_relative /
                              (RADIATION_MU_H * RADIATION_HYDROGEN_MASS_CGS);
  return kappa_eff_cgs * units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
         units_cgs_conversion_factor(us, UNIT_CONV_AREA);
}

/**
 * Band-specific dust extinction factor: exp(-kappa_eff * Sigma_gas_p).
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param Sigma_gas_p Physical gas column density, internal units.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Dust extinction factor, in (0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_extinction_factor(const struct unit_system *us, float Z,
                                     float sigma_d_band_cgs, float Sigma_gas_p,
                                     float local_dust_to_gas_ratio) {

  const float kappa_eff = radiation_get_dust_mass_opacity(
      us, Z, sigma_d_band_cgs, local_dust_to_gas_ratio);
  return expf(-kappa_eff * Sigma_gas_p);
}

/**
 * Local linear dust absorption rate (1/length): kappa_eff(Z) * rho. The
 * The `lambda = 1/kappa` screening length and hyperbolic relaxation
 * time are both built from this. Purely local: no column density
 * involved, unlike the injection extinction above. Returns the raw
 * physical rate.
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param rho_p Physical gas density, internal units.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param local_dust_to_gas_ratio See
 * #radiation_get_dust_to_gas_ratio_relative_to_MW.
 * @return Local linear dust absorption rate, internal units (1/length).
 */
__attribute__((always_inline)) INLINE float
radiation_get_part_linear_absorption_rate(const struct unit_system *us, float Z,
                                          float rho_p, float sigma_d_band_cgs,
                                          float local_dust_to_gas_ratio) {
  return radiation_get_dust_mass_opacity(us, Z, sigma_d_band_cgs,
                                         local_dust_to_gas_ratio) *
         rho_p;
}

/**
 * @brief Receiver-side LW/PE dust extinction factors for a gas particle.
 *
 * @param us Unit system.
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
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct part *p, float Z, const struct cooling_function_data *cooling,
    const float extinction_path, float extinction[ISRF_OPERATOR_COUNT]) {

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(p, extinction_path) *
      cosmo->a2_inv;
  /* Resolved value, never the raw `-1`-sentinel
   * cooling->local_dust_to_gas_ratio. */
  const float local_dust_to_gas_ratio =
      (float)cooling->chemistry_data.local_dust_to_gas_ratio;

  extinction[ISRF_OPERATOR_PE] = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_PE_CGS, Sigma_gas_p, local_dust_to_gas_ratio);
  extinction[ISRF_OPERATOR_LW] = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p, local_dust_to_gas_ratio);
}
