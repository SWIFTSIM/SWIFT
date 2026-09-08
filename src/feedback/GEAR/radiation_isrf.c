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
  radiation_init_part_propagation(p);
}

/**
 * @brief Snapshot #u_FUV/#u_LW once per step, before the density loop's
 * h-iterations begin (see #feedback_part_data.u_FUV_prev), cache this
 * step's per-band absorption rate, cache a stable comoving-density
 * snapshot the propagation loops need (see
 * #feedback_part_data.rho_prev's own doxygen for why), cache
 * this step's hyperbolic propagation speed and physical timestep, and zero
 * the per-step `grad(u)` accumulators.
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
  float dt_phys;
  if (with_cosmology) {
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current, p->time_bin);
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
}

/**
 * @brief Zero the `div(F)` per-h-iteration accumulators. Mirrors
 * chemistry_init_part's own per-iteration reset (called from the same
 * sites: part_init.h and the ghost h-iteration redo path), so it is safe
 * to call once or several times per step. Pure scratch space (no restart
 * I/O); the density snapshot and c_hyp/dt are cached separately, once per
 * step, by #radiation_snapshot_part_propagation.
 *
 * @param p The #part to reset.
 */
void radiation_init_part_propagation(struct part *p) {
  p->feedback_data.div_specific_flux_FUV = 0.f;
  p->feedback_data.div_specific_flux_LW = 0.f;
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
 * accumulators
 * radiation_propagation_iact.h filled during the density loop, which runs
 * before star feedback: this always stamps LW_FUV_last_touch_ti, so
 * injection adds on top instead of resetting. Idempotent: always
 * recomputed from the stable #u_FUV_prev snapshot and this h-iteration's
 * `div(F)` accumulator, so repeated calls across h-iterations converge to
 * the same answer regardless of how many there are.
 *
 * `u_new = e*u_prev - dt*phi*div_F`, with `e = exp(-a)`,
 * `phi = (1-e)/a`, `a = c_hyp*kappa*dt` (#radiation_relaxation_phi_factor): the
 * exact solution of `du/dt = -u/tau - div(F)` over one step with `div(F)`
 * frozen at this h-iteration's value, `tau = 1/(c_hyp*kappa)`. No-op when
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

  const float a_FUV = c_hyp * fd->kappa_FUV * dt;
  const float a_LW = c_hyp * fd->kappa_LW * dt;
  const float decay_FUV = expf(-a_FUV);
  const float decay_LW = expf(-a_LW);
  const float phi_FUV = radiation_relaxation_phi_factor(a_FUV);
  const float phi_LW = radiation_relaxation_phi_factor(a_LW);

  fd->u_FUV =
      decay_FUV * fd->u_FUV_prev - dt * phi_FUV * fd->div_specific_flux_FUV;
  fd->u_LW = decay_LW * fd->u_LW_prev - dt * phi_LW * fd->div_specific_flux_LW;
  fd->LW_FUV_last_touch_ti = e->ti_current;
}

/**
 * @brief Exact-relaxation update of #specific_flux_FUV/#specific_flux_LW, from
 * the `grad(u)` accumulators radiation_propagation_iact.h filled during the
 * gradient loop, which reads this step's already-relaxed `u`
 * (#radiation_end_density_propagation having already run in the density ghost).
 * Runs once per step in the extra ghost, never re-run: the gradient loop itself
 * only runs once per step.
 *
 * `F_new = e*F - c_hyp^2*dt*phi*grad(u)`: the exact solution of
 * `dF/dt = -F/tau - (D/tau)*grad(u)` over one step with `grad(u)` frozen
 * at this step's value, `D/tau = c_hyp^2`. No-op when propagation is off.
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
  const float coeff_FUV =
      c_hyp * c_hyp * dt * radiation_relaxation_phi_factor(a_FUV);
  const float coeff_LW =
      c_hyp * c_hyp * dt * radiation_relaxation_phi_factor(a_LW);

  for (int k = 0; k < 3; k++) {
    fd->specific_flux_FUV[k] =
        decay_FUV * fd->specific_flux_FUV[k] - coeff_FUV * fd->grad_u_FUV[k];
    fd->specific_flux_LW[k] =
        decay_LW * fd->specific_flux_LW[k] - coeff_LW * fd->grad_u_LW[k];
  }
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
