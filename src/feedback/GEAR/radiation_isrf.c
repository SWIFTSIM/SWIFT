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
 * @brief Receiver-side LW/FUV dust extinction and hyperbolic
 * P1-relaxation propagation physics for GEAR.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "active.h"
#include "chemistry.h"
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
#include "timeline.h"

#include <float.h>
#include <math.h>

/**
 * @brief First-init of a #part's LW/FUV radiation-field state. Shared
 * across GEAR feedback variants: independent of the injection mechanism.
 *
 * Deliberately does NOT zero #feedback_part_data.u_FUV/u_LW: this runs
 * (space_first_init.c) after the IC file has been read into #part
 * (single_io.c/parallel_io.c/serial_io.c), and an IC may supply them via
 * the optional "FUVSpecificEnergy"/"LWSpecificEnergy" fields (see
 * src/feedback/GEAR_thermal/feedback_io.h) for a validation setup that
 * bypasses star injection entirely; zeroing here would silently stomp
 * that value back to 0.f. An IC that does not supply them is unaffected:
 * every #part is bzero'd before the IC read runs, so u_FUV/u_LW are
 * already 0.f by the time this function is reached, identical to the
 * previous unconditional assignment.
 *
 * #u_FUV_prev/#u_LW_prev are seeded from #u_FUV/#u_LW rather than left at
 * 0.f, for the same reason: with `LW_FUV_propagation` on, the engine's
 * initial density computation (before the first real step's
 * #radiation_snapshot_part_propagation call has ever run) calls
 * #radiation_end_density_propagation directly off this first-init state.
 * A 0.f seed there would make that very first propagation update read
 * `u_FUV_prev=0` while `u_FUV` holds the IC value, wiping an IC-supplied
 * field to (near) 0 before it is ever seen. Seeding from #u_FUV/#u_LW
 * reduces to today's behaviour when the IC does not supply them (both
 * already 0.f from the bzero above).
 *
 * #specific_flux_FUV/#specific_flux_LW are zeroed unconditionally: no IC field
 * is proposed for them, so the "don't stomp an IC value" concern above does not
 * apply.
 *
 * @param p The #part to initialise.
 */
void radiation_first_init_part(struct part *restrict p) {
  p->feedback_data.u_FUV_prev = p->feedback_data.u_FUV;
  p->feedback_data.u_LW_prev = p->feedback_data.u_LW;
  p->feedback_data.kappa_FUV = 0.f;
  p->feedback_data.kappa_LW = 0.f;
  p->feedback_data.LW_FUV_last_touch_ti = -1;
  /* -1 so an MPI foreign particle's uninitialized memory (not covered by
     the IC-read bzero above) can never read as "still illuminated". */
  p->feedback_data.is_illuminated_LW_FUV = 0;
  p->feedback_data.LW_FUV_illumination_end_ti = -1;
  p->feedback_data.specific_flux_FUV[0] = 0.f;
  p->feedback_data.specific_flux_FUV[1] = 0.f;
  p->feedback_data.specific_flux_FUV[2] = 0.f;
  p->feedback_data.specific_flux_LW[0] = 0.f;
  p->feedback_data.specific_flux_LW[1] = 0.f;
  p->feedback_data.specific_flux_LW[2] = 0.f;
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
  p->feedback_data.rho_prev = 1.0f;
  p->feedback_data.c_hyp = 0.f;
  p->feedback_data.dt_prev = 0.f;
  p->feedback_data.u_FUV_dose_reservoir = 0.f;
  p->feedback_data.u_LW_dose_reservoir = 0.f;
  p->feedback_data.u_FUV_source_rate = 0.f;
  p->feedback_data.u_LW_source_rate = 0.f;
  p->feedback_data.LW_FUV_reservoir_end_ti = -1;
  p->feedback_data.dissipation_alpha_FUV = 0.f;
  p->feedback_data.dissipation_alpha_LW = 0.f;
  p->feedback_data.dissipation_u_FUV = 0.f;
  p->feedback_data.dissipation_u_LW = 0.f;
#ifdef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
  for (int k = 0; k < 3; k++) {
    p->feedback_data.grad_u_FUV_prev[k] = 0.f;
    p->feedback_data.grad_u_LW_prev[k] = 0.f;
  }
#endif
#ifdef RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX
  p->feedback_data.dissipation_alpha_flux_FUV = 0.f;
  p->feedback_data.dissipation_alpha_flux_LW = 0.f;
  p->feedback_data.div_specific_flux_FUV_prev = 0.f;
  p->feedback_data.div_specific_flux_LW_prev = 0.f;
  for (int k = 0; k < 3; k++) {
    p->feedback_data.dissipation_F_FUV[k] = 0.f;
    p->feedback_data.dissipation_F_LW[k] = 0.f;
  }
#endif
  radiation_init_part_propagation(p);
}

/**
 * @brief Snapshot #u_FUV/#u_LW once per step, before the density loop's
 * h-iterations begin (see #feedback_part_data.u_FUV_prev), cache this
 * step's per-band absorption rate, cache a stable comoving-density
 * snapshot the propagation loops need (see
 * #feedback_part_data.rho_prev's own doxygen for why), cache
 * this step's hyperbolic propagation speed and physical timestep, zero
 * every per-step gradient-loop accumulator, and, for active particles only,
 * draw down this step's #u_FUV_source_rate/#u_LW_source_rate from
 * #u_FUV_dose_reservoir/#u_LW_dose_reservoir.
 *
 * #grad_u_FUV_prev/#grad_u_LW_prev are NOT written here: this function runs
 * at drift time for every particle regardless of activity, so copying
 * #grad_u_FUV/LW here would overwrite an inactive particle's last real
 * gradient with whatever this same unconditional zeroing already left there
 * on a previous drift. They are instead written at the end of
 * #radiation_end_gradient_propagation, which runs only for active
 * particles, from that step's own just-finalised gradient -- the same
 * pattern #div_specific_flux_FUV_prev/LW_prev already uses, and the one
 * MAGMA2's `hydro_prepare_force` uses for the analogous hydro gradient.
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
  p->feedback_data.u_FUV_prev = p->feedback_data.u_FUV;
  p->feedback_data.u_LW_prev = p->feedback_data.u_LW;

  p->feedback_data.grad_u_FUV[0] = 0.f;
  p->feedback_data.grad_u_FUV[1] = 0.f;
  p->feedback_data.grad_u_FUV[2] = 0.f;
  p->feedback_data.grad_u_LW[0] = 0.f;
  p->feedback_data.grad_u_LW[1] = 0.f;
  p->feedback_data.grad_u_LW[2] = 0.f;

  /* Force-loop accumulator, zeroed here for the same reason as grad_u
   * above: that loop, like the gradient loop, runs exactly once per step,
   * with no h-iteration redo. */
  p->feedback_data.dissipation_u_FUV = 0.f;
  p->feedback_data.dissipation_u_LW = 0.f;

  /* The other gradient-loop accumulators, zeroed here for the same reason as
   * grad_u above: that loop runs exactly once per step. */
#ifdef RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX
  for (int k = 0; k < 3; k++) {
    p->feedback_data.dissipation_F_FUV[k] = 0.f;
    p->feedback_data.dissipation_F_LW[k] = 0.f;
  }
#endif
  /* Stable comoving density snapshot, cached unconditionally (not gated on
   * LW_FUV_propagation below): the gradient loop's `grad(u)` accumulation
   * (radiation_propagation_iact.h) always runs, even in injection-only
   * mode (mirroring Design A's own always-accumulate pattern), and reads
   * this snapshot. See this field's own doxygen (feedback_struct.h) for
   * why the density loop cannot use a live `p->rho` instead, and why it
   * must never be left at 0.f. */
  const float rho_comoving = hydro_get_comoving_density(p);
  p->feedback_data.rho_prev = rho_comoving > 0.f ? rho_comoving : 1.0f;

  if (!e->feedback_props->LW_FUV_propagation) {
    p->feedback_data.kappa_FUV = 0.f;
    p->feedback_data.kappa_LW = 0.f;
    return;
  }

  const float rho_phys = hydro_get_physical_density(p, e->cosmology);
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  /* Resolved value, never the raw `-1`-sentinel
   * cooling->local_dust_to_gas_ratio. */
  const float local_dust_to_gas_ratio =
      (float)e->cooling_func->chemistry_data.local_dust_to_gas_ratio;
  p->feedback_data.kappa_FUV = radiation_get_part_linear_absorption_rate(
      e->internal_units, Z, rho_phys, RADIATION_SIGMA_D_FUV_CGS,
      local_dust_to_gas_ratio);
  p->feedback_data.kappa_LW = radiation_get_part_linear_absorption_rate(
      e->internal_units, Z, rho_phys, RADIATION_SIGMA_D_LW_CGS,
      local_dust_to_gas_ratio);

  /* Hyperbolic propagation speed closure: c_hyp_i = min(C_hyp*h_i/dt_i, c),
   * using this particle's own already-decided integer timestep -- not a
   * new timestep-computation hook. dt_i is floored at FLT_MIN so a
   * not-yet-assigned time_bin (only possible before this particle's very
   * first real step) cannot divide by an exact zero. */
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

  const float h_phys = (float)e->cosmology->a * p->h;
  float c_hyp = e->feedback_props->LW_FUV_c_hyp_margin * h_phys / dt_phys;
  c_hyp = min(c_hyp, (float)e->physical_constants->const_speed_light_c);
  /* The debug pin is applied after the light-speed clamp above and is not
   * itself clamped: a pin value above c gives a superluminal propagation
   * speed on purpose, for isolating dispersion behaviour at chosen values
   * of the Courant number. Never set it above c outside of that use. */
  if (e->feedback_props->LW_FUV_c_hyp_pin_for_debugging > 0.f)
    c_hyp = e->feedback_props->LW_FUV_c_hyp_pin_for_debugging;

  p->feedback_data.c_hyp = c_hyp;
  p->feedback_data.dt_prev = dt_phys;

  /* Dose-reservoir drawdown (design-lw-fuv-design-b-dissipation.md
   * Section 4.6.5), for active particles only: a cell drifted for an
   * inactive particle must not draw down a dose it will not integrate this
   * step. An inactive particle's u_*_source_rate is simply left at last
   * step's value; it is never read again before this function next runs
   * for it (as active) and overwrites it. */
  struct feedback_part_data *fd = &p->feedback_data;
  if (part_is_active(p, e) &&
      (fd->u_FUV_dose_reservoir > 0.f || fd->u_LW_dose_reservoir > 0.f)) {
    double t_rem;
    if (fd->LW_FUV_reservoir_end_ti <= ti_begin) {
      t_rem = 0.0;
    } else if (with_cosmology) {
      t_rem = cosmology_get_delta_time(e->cosmology, ti_begin,
                                       fd->LW_FUV_reservoir_end_ti);
    } else {
      t_rem = (double)(fd->LW_FUV_reservoir_end_ti - ti_begin) * e->time_base;
    }
    const float f =
        (t_rem <= (double)dt_phys) ? 1.f : (float)((double)dt_phys / t_rem);
    fd->u_FUV_source_rate = f * fd->u_FUV_dose_reservoir / dt_phys;
    fd->u_LW_source_rate = f * fd->u_LW_dose_reservoir / dt_phys;
    fd->u_FUV_dose_reservoir -= f * fd->u_FUV_dose_reservoir;
    fd->u_LW_dose_reservoir -= f * fd->u_LW_dose_reservoir;
  } else if (part_is_active(p, e)) {
    fd->u_FUV_source_rate = 0.f;
    fd->u_LW_source_rate = 0.f;
  }
}

/**
 * @brief Zero the `div(F)` and kernel-mean per-h-iteration accumulators.
 * Mirrors chemistry_init_part's own per-iteration reset (called from the
 * same sites: part_init.h and the ghost h-iteration redo path), so it is
 * safe to call once or several times per step. Pure scratch space (no
 * restart I/O); the density snapshot and c_hyp/dt are cached separately,
 * once per step, by #radiation_snapshot_part_propagation,
 * #dissipation_u_FUV/LW is a force-loop accumulator zeroed there too, and
 * #dissipation_alpha_FUV/LW are persistent and updated once per step by
 * #radiation_end_gradient_propagation, not here.
 *
 * @param p The #part to reset.
 */
void radiation_init_part_propagation(struct part *p) {
  p->feedback_data.div_specific_flux_FUV = 0.f;
  p->feedback_data.div_specific_flux_LW = 0.f;
  p->feedback_data.ngb_mean_abs_u_V_FUV = 0.f;
  p->feedback_data.ngb_mean_abs_u_V_LW = 0.f;
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
 * @param a Dimensionless absorption depth for this step, `c_hyp*kappa*dt`
 * (or the injection-site equivalent, `c_hyp*kappa*Delta_t_star`). Always
 * `>= 0`.
 * @return phi(a).
 */
float radiation_relaxation_phi_factor(float a) {
  if (a < 1e-6f) return 1.0f - 0.5f * a + (1.0f / 6.0f) * a * a;
  return -expm1f(-a) / a;
}

/**
 * @brief Exact-relaxation update of #u_FUV/#u_LW, from the `div(F)`
 * accumulators radiation_propagation_iact.h filled during the density loop
 * and this step's own #u_FUV_source_rate/#u_LW_source_rate (drawn down from
 * the dose reservoir by #radiation_snapshot_part_propagation,
 * design-lw-fuv-design-b-dissipation.md Section 4.6.5). Idempotent: always
 * recomputed from the stable #u_FUV_prev snapshot and this h-iteration's
 * `div(F)` accumulator, so repeated calls across h-iterations converge to
 * the same answer regardless of how many there are.
 *
 * `u_star = e*u_prev + dt*phi*((3*c_hyp/c)*source_rate - div_F)`, with
 * `e = exp(-a)`, `phi = (1-e)/a`, `a = c_hyp*kappa*dt`
 * (#radiation_relaxation_phi_factor): the exact solution of
 * `du/dt = -u/tau + (3*c_hyp/c)*source_rate - div(F)` over one step with
 * `source_rate` and `div(F)` frozen at this h-iteration's value, `tau =
 * 1/(c_hyp*kappa)`. The `3*c_hyp/c` rescale is applied exclusively here;
 * injection (`radiation_iact.h`) deposits the raw, unrescaled dose.
 *
 * The result is the INTERMEDIATE state `u*`, not this step's final `u`:
 * the Stage-1 artificial-dissipation correction is added on top of it by
 * #radiation_end_force_propagation, after the force loop has accumulated
 * the mirrored pairwise term. The gradient loop and the negativity trigger
 * therefore both see `u*`, which is what closes the trigger's one-step lag
 * (design-lw-fuv-design-b-dissipation.md Section 4.3). No-op when
 * propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_density_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->LW_FUV_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  const float dt = fd->dt_prev;
  const float c_hyp = fd->c_hyp;
  const float rescale =
      3.0f * c_hyp / (float)e->physical_constants->const_speed_light_c;

  const float a_FUV = c_hyp * fd->kappa_FUV * dt;
  const float a_LW = c_hyp * fd->kappa_LW * dt;
  const float decay_FUV = expf(-a_FUV);
  const float decay_LW = expf(-a_LW);
  const float phi_FUV = radiation_relaxation_phi_factor(a_FUV);
  const float phi_LW = radiation_relaxation_phi_factor(a_LW);

  fd->u_FUV = decay_FUV * fd->u_FUV_prev +
              dt * phi_FUV *
                  (rescale * fd->u_FUV_source_rate - fd->div_specific_flux_FUV);
  fd->u_LW =
      decay_LW * fd->u_LW_prev +
      dt * phi_LW * (rescale * fd->u_LW_source_rate - fd->div_specific_flux_LW);
}

/**
 * @brief Stage-1 artificial-dissipation correction of #u_FUV/#u_LW, from
 * the accumulators radiation_propagation_iact.h filled during the force
 * loop: `u = u_star + dt*phi*dissipation_u`, closing the exact-relaxation
 * update #radiation_end_density_propagation left at its intermediate state
 * `u_star`.
 *
 * Runs in the `end_force` task, which SWIFT places after the force loop
 * and before cooling (engine_maketasks.c), so a coefficient raised by this
 * step's negativity trigger (the extra ghost, which precedes the force
 * loop) acts on this step's own `u` rather than the next step's.
 *
 * HARD INVARIANT: `end_force` runs exactly once per active particle per
 * step, with no h-iteration redo of the kind the density ghost has. The
 * `+=` below is NOT idempotent; a second call would double-correct.
 *
 * Reads #dt_prev/#c_hyp/#kappa_FUV/LW, all cached earlier in this same step
 * by #radiation_snapshot_part_propagation, and deliberately takes no `dt`
 * of its own: the call site computes its local `dt` from a different
 * timestep-begin convention (`ti_current - 1`), and using it here would
 * make this correction's `dt*phi` inconsistent with the rest of the step's
 * radiation update. The thin `(p, e)` signature exists to make that
 * mistake structurally impossible.
 *
 * `phi` is recomputed rather than cached: this is an additive correction,
 * not a decay-weighted blend, so no `exp(-a)` memory term is needed.
 * No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_force_propagation(struct part *p, const struct engine *e) {

  if (!e->feedback_props->LW_FUV_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  const float dt = fd->dt_prev;
  const float c_hyp = fd->c_hyp;

  const float a_FUV = c_hyp * fd->kappa_FUV * dt;
  const float a_LW = c_hyp * fd->kappa_LW * dt;
  const float phi_FUV = radiation_relaxation_phi_factor(a_FUV);
  const float phi_LW = radiation_relaxation_phi_factor(a_LW);

  fd->u_FUV += dt * phi_FUV * fd->dissipation_u_FUV;
  fd->u_LW += dt * phi_LW * fd->dissipation_u_LW;
}

/**
 * @brief Undo this step's dose-reservoir drawdown for a #part whose density
 * h-iteration gives up with no neighbours found:
 * #radiation_end_density_propagation, the only consumer of
 * #u_FUV_source_rate/#u_LW_source_rate, is never reached in that case
 * (`runner_ghost.c`'s `has_no_neighbours` give-up path), so the rate drawn down
 * by #radiation_snapshot_part_propagation would otherwise be discarded rather
 * than applied. Restores it into #u_FUV_dose_reservoir/#u_LW_dose_reservoir
 * exactly (the drawn amount is `source_rate*dt_prev`) and zeroes the rates.
 * No-op when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_part_has_no_neighbours(struct part *p, const struct engine *e) {

  if (!e->feedback_props->LW_FUV_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  fd->u_FUV_dose_reservoir += fd->u_FUV_source_rate * fd->dt_prev;
  fd->u_LW_dose_reservoir += fd->u_LW_source_rate * fd->dt_prev;
  fd->u_FUV_source_rate = 0.f;
  fd->u_LW_source_rate = 0.f;
}

/**
 * @brief One band's Stage-3 anisotropic flux-dissipation coefficient update
 * (design-lw-fuv-design-b-dissipation.md Section 5.2): Chan et al. 2021
 * Eq. 36-37, raised instantly on a steepening `div(F)` in compression and
 * decayed at the same rate as the Stage-1 coefficient otherwise.
 *
 * A particle whose `u_V` has already gone non-positive gets the ceiling,
 * mirroring the `urad == 0` branch of `src/rt/SPHM1RT/rt.h`: the switch's
 * denominator is the local radiation energy the term is meant to protect,
 * so where there is none left, no smallness argument applies.
 *
 * @param div_F This band's own finalized `div(F)`.
 * @param div_F_prev This band's #div_specific_flux_FUV_prev/LW_prev.
 * @param u_V This band's live volumetric field, rho_prev*u.
 * @param alpha_prev This band's #dissipation_alpha_flux_FUV/LW from the
 * previous step.
 * @param c_hyp The particle's own #c_hyp.
 * @param kappa This band's #kappa_FUV/LW.
 * @param dt The particle's own #dt_prev.
 * @param h_phys The particle's own physical smoothing length.
 * @return This step's updated flux-dissipation coefficient for this band.
 */
__attribute__((always_inline)) INLINE static float
radiation_update_dissipation_alpha_flux_band(float div_F, float div_F_prev,
                                             float u_V, float alpha_prev,
                                             float c_hyp, float kappa, float dt,
                                             float h_phys) {

  float alpha_aim = 0.f;

  /* The flux term only acts in compression, as Chan et al. Eq. 36 does. */
  if (div_F < 0.f) {
    float shock_estimate = 1.f;
    if (u_V > 0.f && c_hyp > 0.f) {
      const float div_F_rate = (div_F - div_F_prev) / dt;
      shock_estimate = -RADIATION_LW_FUV_DISSIPATION_FLUX_SWITCH_AMPLITUDE *
                       h_phys * h_phys * div_F_rate / (u_V * c_hyp * c_hyp);
    }
    const float shock_capped = min(shock_estimate, 1.f);
    alpha_aim = max(shock_capped, 0.f);
  }

  if (alpha_aim >= alpha_prev) return alpha_aim;

  const float a_kappa = c_hyp * kappa * dt;
  const float decay =
      expf(-c_hyp * dt / (RADIATION_LW_FUV_DISSIPATION_DECAY_LENGTH * h_phys) -
           a_kappa);
  return alpha_aim + (alpha_prev - alpha_aim) * decay;
}

/**
 * @brief One band's Stage-1 artificial-dissipation coefficient update
 * (design-lw-fuv-design-b-dissipation.md Section 4.3): raised instantly to
 * a negativity-triggered target, or decayed toward it otherwise. Reads the
 * particle's LIVE, this-step `u_V = rho_prev*u_star` rather than the
 * `u_*_prev` snapshot the design document's own text specifies: with the
 * dose-reservoir injection form (Section 4.6.5) already landed, injection
 * no longer writes `u` at all, so the value seen here (after
 * #radiation_end_density_propagation has already run in the density ghost,
 * before any star touches this step's `u` again) is this step's own
 * transported state. The coefficient this returns is consumed by THIS
 * step's force loop, so the trigger has no lag left: an undershoot is
 * corrected in the step it appears, before cooling reads `u`.
 *
 * @param u_V This band's live volumetric field, rho_prev*u_star.
 * @param ngb_mean_abs_u_V This band's kernel-mean |rho_prev*u_prev| scratch
 * accumulator (radiation_propagation_iact.h).
 * @param alpha_prev This band's #dissipation_alpha_FUV/LW from the
 * previous step.
 * @param alpha_max #feedback_props.LW_FUV_dissipation_alpha_max.
 * @param eps_1 #feedback_props.LW_FUV_dissipation_negativity_threshold.
 * @param c_hyp The particle's own #c_hyp.
 * @param kappa This band's #kappa_FUV/LW.
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
      expf(-c_hyp * dt / (RADIATION_LW_FUV_DISSIPATION_DECAY_LENGTH * h_phys) -
           a_kappa);
  return alpha_aim + (alpha_prev - alpha_aim) * decay;
}

/**
 * @brief The `h/lambda`-gated floor under
 * #radiation_update_dissipation_alpha_band's trigger
 * (PHASE5B_diffuse_phase_fable_review_2026-09-11.md Section 4): the trigger
 * fires only on negativity and is exactly zero on the positive delta-shell
 * front of an optically-thin P1 pulse, so a purely reactive coefficient
 * cannot damp the resulting dispersive wake there. This floor supplies
 * dissipation the trigger structurally cannot, rolling off as
 * `(eps_lambda/(h*kappa))^2` once `h/lambda` exceeds #LW_FUV_dissipation_
 * floor_h_over_lambda, which bounds its steady-state cost by construction.
 *
 * @param kappa This band's #kappa_FUV/LW.
 * @param h_phys The particle's own physical smoothing length.
 * @param alpha_floor #feedback_props.LW_FUV_dissipation_alpha_floor.
 * @param eps_lambda #feedback_props.LW_FUV_dissipation_floor_h_over_lambda.
 * @return This band's floor value for this step, to be combined with the
 * trigger's own output via max().
 */
__attribute__((always_inline)) INLINE static float
radiation_dissipation_alpha_floor_band(float kappa, float h_phys,
                                       float alpha_floor, float eps_lambda) {

  const float x = h_phys * kappa / eps_lambda;
  return alpha_floor / (1.f + x * x);
}

/**
 * @brief Exact-relaxation update of #specific_flux_FUV/#specific_flux_LW, from
 * the `grad(u)` accumulators radiation_propagation_iact.h filled during the
 * gradient loop, which reads this step's already-relaxed `u`
 * (#radiation_end_density_propagation having already run in the density ghost).
 * Runs once per step in the extra ghost, never re-run: the gradient loop itself
 * only runs once per step. Also updates #dissipation_alpha_FUV/LW once per
 * step (see #radiation_update_dissipation_alpha_band), for THIS step's
 * force loop to consume: the extra ghost precedes the force loop.
 *
 * `F_new = e*F - c_hyp^2*dt*phi*grad(u) + dt*phi*dissipation_F`: the exact
 * solution of `dF/dt = -F/tau - (D/tau)*grad(u) + dissipation_F` over one
 * step with both source terms frozen at this step's value, `D/tau =
 * c_hyp^2`. `dissipation_F` is the Stage-3 anisotropic term
 * radiation_propagation_iact.h accumulates alongside `grad(u)` and is zero
 * unless #RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX is defined. No-op
 * when propagation is off.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void radiation_end_gradient_propagation(struct part *p,
                                        const struct engine *e) {

  if (!e->feedback_props->LW_FUV_propagation) return;

  struct feedback_part_data *fd = &p->feedback_data;
  const float dt = fd->dt_prev;
  const float c_hyp = fd->c_hyp;

  const float a_FUV = c_hyp * fd->kappa_FUV * dt;
  const float a_LW = c_hyp * fd->kappa_LW * dt;
  const float decay_FUV = expf(-a_FUV);
  const float decay_LW = expf(-a_LW);

  const float phi_FUV = radiation_relaxation_phi_factor(a_FUV);
  const float phi_LW = radiation_relaxation_phi_factor(a_LW);
  const float coeff_FUV = c_hyp * c_hyp * dt * phi_FUV;
  const float coeff_LW = c_hyp * c_hyp * dt * phi_LW;

  for (int k = 0; k < 3; k++) {
    fd->specific_flux_FUV[k] =
        decay_FUV * fd->specific_flux_FUV[k] - coeff_FUV * fd->grad_u_FUV[k];
    fd->specific_flux_LW[k] =
        decay_LW * fd->specific_flux_LW[k] - coeff_LW * fd->grad_u_LW[k];
  }

  const float h_phys = (float)e->cosmology->a * p->h;
  const float u_V_FUV = fd->rho_prev * fd->u_FUV;
  const float u_V_LW = fd->rho_prev * fd->u_LW;

  /* Stage 3's own fields only exist when the stage is built, so they are
   * reached through pointers selected here; the block below stays compiled
   * in either state and is removed by the optimizer when the flag is 0. */
#ifdef RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX
  const int use_anisotropic_flux = 1;
  float *const diss_F_FUV = fd->dissipation_F_FUV;
  float *const diss_F_LW = fd->dissipation_F_LW;
  float *const alpha_f_FUV = &fd->dissipation_alpha_flux_FUV;
  float *const alpha_f_LW = &fd->dissipation_alpha_flux_LW;
  float *const div_F_FUV_prev = &fd->div_specific_flux_FUV_prev;
  float *const div_F_LW_prev = &fd->div_specific_flux_LW_prev;
#else
  const int use_anisotropic_flux = 0;
  float absent_diss_F_FUV[3] = {0.f, 0.f, 0.f};
  float absent_diss_F_LW[3] = {0.f, 0.f, 0.f};
  float absent_alpha_f_FUV = 0.f;
  float absent_alpha_f_LW = 0.f;
  float absent_div_F_FUV_prev = 0.f;
  float absent_div_F_LW_prev = 0.f;
  float *const diss_F_FUV = absent_diss_F_FUV;
  float *const diss_F_LW = absent_diss_F_LW;
  float *const alpha_f_FUV = &absent_alpha_f_FUV;
  float *const alpha_f_LW = &absent_alpha_f_LW;
  float *const div_F_FUV_prev = &absent_div_F_FUV_prev;
  float *const div_F_LW_prev = &absent_div_F_LW_prev;
#endif

  if (use_anisotropic_flux) {
    /* Added separately rather than as a third term of the relaxation above,
     * so a build without Stage 3 evaluates the same expression it did
     * before the stage existed. */
    for (int k = 0; k < 3; k++) {
      fd->specific_flux_FUV[k] += dt * phi_FUV * diss_F_FUV[k];
      fd->specific_flux_LW[k] += dt * phi_LW * diss_F_LW[k];
    }

    *alpha_f_FUV = radiation_update_dissipation_alpha_flux_band(
        fd->div_specific_flux_FUV, *div_F_FUV_prev, u_V_FUV, *alpha_f_FUV,
        c_hyp, fd->kappa_FUV, dt, h_phys);
    *alpha_f_LW = radiation_update_dissipation_alpha_flux_band(
        fd->div_specific_flux_LW, *div_F_LW_prev, u_V_LW, *alpha_f_LW, c_hyp,
        fd->kappa_LW, dt, h_phys);
    *div_F_FUV_prev = fd->div_specific_flux_FUV;
    *div_F_LW_prev = fd->div_specific_flux_LW;
  }

  const float alpha_pin =
      e->feedback_props->LW_FUV_dissipation_alpha_pin_for_debugging;
  if (alpha_pin > 0.f) {
    /* Bypass the trigger entirely: every particle's coefficient is held at
     * the pinned value (see this parameter's own doxygen,
     * feedback_properties.h). */
    fd->dissipation_alpha_FUV = alpha_pin;
    fd->dissipation_alpha_LW = alpha_pin;
  } else {
    const float alpha_max = e->feedback_props->LW_FUV_dissipation_alpha_max;
    const float eps_1 =
        e->feedback_props->LW_FUV_dissipation_negativity_threshold;
    const float alpha_floor = e->feedback_props->LW_FUV_dissipation_alpha_floor;
    const float eps_lambda =
        e->feedback_props->LW_FUV_dissipation_floor_h_over_lambda;

    const float alpha_trigger_FUV = radiation_update_dissipation_alpha_band(
        u_V_FUV, fd->ngb_mean_abs_u_V_FUV, fd->dissipation_alpha_FUV, alpha_max,
        eps_1, c_hyp, fd->kappa_FUV, dt, h_phys);
    const float alpha_trigger_LW = radiation_update_dissipation_alpha_band(
        u_V_LW, fd->ngb_mean_abs_u_V_LW, fd->dissipation_alpha_LW, alpha_max,
        eps_1, c_hyp, fd->kappa_LW, dt, h_phys);

    /* Floor the trigger cannot suppress: applies unconditionally, not only
     * on negativity (see radiation_dissipation_alpha_floor_band's doxygen). */
    const float alpha_floor_FUV = radiation_dissipation_alpha_floor_band(
        fd->kappa_FUV, h_phys, alpha_floor, eps_lambda);
    const float alpha_floor_LW = radiation_dissipation_alpha_floor_band(
        fd->kappa_LW, h_phys, alpha_floor, eps_lambda);

    fd->dissipation_alpha_FUV = max(alpha_trigger_FUV, alpha_floor_FUV);
    fd->dissipation_alpha_LW = max(alpha_trigger_LW, alpha_floor_LW);
  }

  /* Written here, at the end of this active-gated ghost, from the gradient
   * this same active step just finalised above (#grad_u_FUV/LW), so the
   * value read by a neighbour is always this particle's own last real
   * gradient regardless of how many inactive steps it takes in between --
   * never a value zeroed by an unrelated drift. Stage 2 reads it in the
   * FORCE loop, which this ghost precedes, so an active particle's
   * reconstruction now uses this step's own finalised gradient rather than
   * the previous step's. */
#ifdef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
  for (int k = 0; k < 3; k++) {
    fd->grad_u_FUV_prev[k] = fd->grad_u_FUV[k];
    fd->grad_u_LW_prev[k] = fd->grad_u_LW[k];
  }
#endif
}

/**
 * Comoving gas column density at a gas particle's own location: the
 * receiver-side analogue of the star-side Sobolev column
 * (#radiation_get_comoving_gas_column_density_at_star), used for LW/FUV
 * extinction. Simplified to the kernel-radius fallback (no resolved
 * density gradient on the gas side).
 *
 * @param p The #part.
 * @return Comoving gas column density at the particle's own location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_part(const struct part *p) {
  const float h_gas = p->h * kernel_gamma;
  return 2.0f * h_gas * p->rho;
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
 * @brief Receiver-side LW/FUV dust extinction factors for a gas particle.
 *
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param p The receiving #part.
 * @param Z The receiving particle's own metal mass fraction.
 * @param cooling The cooling function properties (for the resolved
 * chemistry_data.local_dust_to_gas_ratio).
 * @param extinction_FUV (return) FUV-band extinction factor.
 * @param extinction_LW (return) Lyman-Werner-band extinction factor.
 */
__attribute__((always_inline)) INLINE void
radiation_get_part_LW_FUV_extinction_factors(
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct part *p, float Z, const struct cooling_function_data *cooling,
    float *extinction_FUV, float *extinction_LW) {

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(p) * cosmo->a2_inv;
  /* Resolved value, never the raw `-1`-sentinel
   * cooling->local_dust_to_gas_ratio. */
  const float local_dust_to_gas_ratio =
      (float)cooling->chemistry_data.local_dust_to_gas_ratio;

  *extinction_FUV = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_FUV_CGS, Sigma_gas_p, local_dust_to_gas_ratio);
  *extinction_LW = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p, local_dust_to_gas_ratio);
}
