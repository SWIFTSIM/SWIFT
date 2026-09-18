/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_FEEDBACK_PROPERTIES_H
#define SWIFT_GEAR_FEEDBACK_PROPERTIES_H

#include "../GEAR/radiation_isrf.h"
#include "../GEAR/stellar_evolution.h"
#include "../GEAR/stellar_evolution_struct.h"
#include "chemistry.h"
#include "hydro_properties.h"
#include "minmax.h"

#include <math.h>
#include <string.h>

#define default_HII_min_density_Hpcm3 1.0
#define default_HII_max_age_Myr 50.0
#define default_HII_rebuild_time_Myr 0.5
#define default_HII_deterministic_boundary_ionization 0
#define default_HII_rebuild_floor_Myr 1e-4
#define default_dt_evolution_factor_max 300.0
#define default_event_dt_floor_Myr 1e-4

/**
 * @brief The different subgrid radiation feedback processes GEAR models.
 */
enum radiation_policy {
  radiation_policy_none = 0,
  /*! Do we want the ionization effect (Strömgren sphere)? */
  radiation_policy_photoionization = (1 << 0),
  /*! Radiation pressure from the stars' bolometric luminosity */
  radiation_policy_radiation_pressure = (1 << 1),
  /* Local Lyman-Werner/FUV feedback: photoelectric (PE) heating by FUV
     radiation on dust, and H2 photodissociation by the Lyman-Werner band.
     One switch for both, since they share the same two band luminosities
     and injected fields. */
  radiation_policy_photoelectric_heating = (1 << 2),
};

/**
 * @brief Which ISRF hyperbolic-propagation scheme is active.
 *
 * Two independent choices sit behind this one selector:
 *  - the FORMULA for the propagation speed `c_hyp_i` (radiation_isrf.c):
 *    #isrf_c_hyp_scheme_shipped's per-particle `dt_i` (the default) or
 *    #isrf_c_hyp_scheme_kernel_local's per-kernel `dt_max(i)`;
 *    #isrf_c_hyp_scheme_fixed_fraction's uniform `f*c` is a third, standalone
 *    option that does not combine with either;
 *  - which PAIRWISE OPERATORS (radiation_propagation_iact.h) consume that
 *    speed: the shipped operators, or the consistent-variable-c change of
 *    variable (#isrf_c_hyp_consistent_variable_c) that rebuilds them for a
 *    genuinely per-particle `c_hyp_i`.
 *
 * #isrf_c_hyp_scheme_kernel_local_plus_variable_c selects BOTH non-default
 * choices at once: the kernel-local speed formula feeding the
 * change-of-variable operators. The two axes are independent by
 * construction (the operator rewrite never reads how `c_hyp_i` was computed,
 * only its value), so this is not a new formula, just the selected
 * combination; see radiation_propagation_iact.h's file header for why the
 * two do not double-throttle the same term the way a kernel-local +
 * pair-weight stack would. Every value is enforced to be one of these five
 * at parse time by feedback_props_init(); see
 * #feedback_props.ISRF_c_hyp_scheme's own doxygen for the specifics of each.
 */
enum isrf_c_hyp_scheme {
  /*! Speed: `c_hyp_i = min(C_hyp*h_i/dt_i, c)`, `dt_i` this particle's own
   * timestep. Operators: shipped. Bit-identical to the scheme this
   * comparison branch was built from; no longer the default (see
   * #isrf_c_hyp_scheme_kernel_local_plus_variable_c), still reachable by
   * setting #feedback_props.ISRF_c_hyp_scheme explicitly. */
  isrf_c_hyp_scheme_shipped = 0,
  /*! Speed: `c_hyp_i = min(C_hyp*h_i/dt_max(i), c)`, `dt_max(i)` the longest
   * timestep among this particle and every neighbour in its kernel.
   * Operators: shipped. */
  isrf_c_hyp_scheme_kernel_local = 1,
  /*! Speed: `c_hyp_i = ISRF_c_hyp_fixed_fraction_of_c * c` for every
   * particle, independent of `h`/timestep. Operators: shipped. Does not
   * combine with either other speed formula. */
  isrf_c_hyp_scheme_fixed_fraction = 2,
  /*! Speed: #isrf_c_hyp_scheme_shipped's own `dt_i` formula, unchanged.
   * Operators: every pairwise transport/dissipation operator
   * (radiation_propagation_iact.h) is rebuilt as the consistent
   * generalisation of the single-uniform-speed reduced-speed-of-light
   * method to a per-particle `c_hyp_i`: a change of variable (every
   * operator at particle i becomes `c_hyp_i/c` times the true-speed
   * equation, and the state stored becomes the reduced flux `Ft =
   * F_true/c_hyp`, so a time-bin change rescales nothing), not a new pair
   * weight. See radiation_propagation_iact.h's file header and
   * #isrf_c_hyp_consistent_variable_c, the global flag that carries this
   * selection into that file's pairwise dispatch. Reduces bit-for-bit to
   * #isrf_c_hyp_scheme_shipped whenever `c_hyp_i` is spatially uniform. */
  isrf_c_hyp_scheme_consistent_variable_c = 3,
  /*! Speed: #isrf_c_hyp_scheme_kernel_local's own `dt_max(i)` formula.
   * Operators: #isrf_c_hyp_scheme_consistent_variable_c's own change of
   * variable. The kernel-local speed reduces the speed contrast between
   * neighbours (and so the negativity); the change of variable fixes the
   * pairwise operators' amplitude error; neither touches the other's
   * mechanism, so the two compose without a combined re-derivation.
   * Default. */
  isrf_c_hyp_scheme_kernel_local_plus_variable_c = 4,
};

/**
 * @brief Properties of the GEAR feedback model.
 */
struct feedback_props {

  /*! Whether sinks are configured; set by sink_props_init(), 0 otherwise. */
  int with_sinks;

  /*! Supernovae energy effectively deposited */
  float supernovae_efficiency;

  /* ------------- Stellar model properties ------------- */

  /*! The stellar model */
  struct stellar_model stellar_model;

  /*! The stellar model for first stars */
  struct stellar_model stellar_model_first_stars;

  /*! Metallicity limits for the first stars */
  float metallicity_max_first_stars;

  /*! Metallicity [Fe/H] transition for the first stars */
  float imf_transition_metallicity;

  /* ------------- Star evolution timestep properties ------------- */

  /*! Timestep refinement factor as lifetime_myr -> 0, used by
   * feedback_compute_spart_timestep()'s logistic transition (SSP particles
   * only; single_star particles use the exact, zero-tuning-parameter
   * dt_event instead, see event_dt_floor_Myr below). The logistic's
   * midpoint and steepness are fixed internal constants
   * (GEAR_dt_evolution_lifetime_myr_0/GEAR_dt_evolution_steepness in
   * feedback_common.c), not parsed from params.yml. */
  float dt_evolution_factor_max;

  /*! Floors the event-anchored timestep terms (single_star's dt_event,
   * and dt_HII_safe for every star type) in internal units, before they are
   * combined with the coarser min_star_timestep-floored non-event terms
   * (stars_compute_timestep(), src/stars/GEAR/stars.h). Runs for every
   * star regardless of radiation, so parsed unconditionally. Unlike
   * HII_rebuild_floor_Myr, which stays 0.0 (no floor at all) whenever
   * with_photoionization is off, this constant must never be zero, since
   * dt_event applies to every single_star particle unconditionally. */
  float event_dt_floor_Myr;

  /* ------------- Subgrid Radiation properties ------------- */

  /* The radiation processes enabled */
  int radiation_policy;

  /*! Radiation pressure momentum effectively injected */
  float radiation_pressure_efficiency;

  /*! Run the hyperbolic P1-relaxation propagation update on top of
   * injection + receiver-side extinction? Only meaningful when
   * radiation_policy_photoelectric_heating is set. */
  char ISRF_propagation;

  /*! Path of the receiver-side LW/FUV dust extinction column, in kernel
   * support radii kernel_gamma * h: 2 for "kernel_diameter", 1 for
   * "kernel_radius" (GEARFeedback:ISRF_extinction_path). */
  float ISRF_extinction_path_in_kernel_radii;

  /*! Stability-margin coefficient in the `c_hyp_i = C_hyp*h_i/dt_max(i)`
   * closure, `dt_max(i)` the longest time step among particle i and every
   * neighbour in its kernel, not just i's own step: an independently-
   * tunable multiple of the hydro CFL margin, rather than silently
   * inheriting whatever SPH:CFL_condition happens to be. Documented valid
   * range (0, sqrt(2/0.70)] (~1.6903); the staggered exact-relaxation
   * scheme is stable for every lambda/h only below that bound (see
   * radiation_isrf.c). That ~1.6903 ceiling only applies at
   * max(#ISRF_dissipation_alpha_max, #ISRF_dissipation_alpha_floor) =
   * 0: the joint stability bound checked in feedback_props_init() is
   * tighter whenever either dissipation coefficient is nonzero (e.g. the
   * shipped alpha_floor=0.5 caps this margin at 0.571), so the effective
   * range depends on both dissipation parameters, not just this one. */
  float ISRF_c_hyp_margin;

  /*! Selects the active ISRF scheme (#isrf_c_hyp_scheme): 0 (shipped),
   * 1 (kernel-local speed), 2 (fixed fraction of c, magnitude
   * #ISRF_c_hyp_fixed_fraction_of_c), 3 (consistent variable-c operators,
   * shipped speed formula) or 4 (kernel-local speed feeding the
   * consistent-variable-c operators, default). See that enum's own doxygen for
   * the specifics of each value. #isrf_c_hyp_scheme_fixed_fraction is an
   * alternative to the other four, not a layer: feedback_props_init()
   * errors if #ISRF_c_hyp_fixed_fraction_of_c is positive with this not set
   * to it, or this is set to it with #ISRF_c_hyp_fixed_fraction_of_c left at
   * 0. */
  int ISRF_c_hyp_scheme;

  /*! Debug/test-only: pin every particle's own `c_hyp_i` (radiation_isrf.c)
   * to this fixed physical value instead of computing it from whichever
   * formula #ISRF_c_hyp_scheme selects, whenever positive. Needed by the
   * causal-reach validation leg
   * (a single, unambiguous wavefront speed to check the field against) and
   * by the steady-state amplitude leg's two-`c_hyp` cross-check (confirming
   * the source-rescaling cancellation empirically). 0 (default): disabled,
   * use the formula. Never set in a production run. */
  float ISRF_c_hyp_pin_for_debugging;

  /*! Uniform reduced light-speed candidate: fraction of the true speed of
   * light `c` every particle's `c_hyp_i` is set to, `c_hyp_i = f*c`, with
   * no dependence on `h_i` or on the particle's own timestep. Only read
   * when #ISRF_c_hyp_scheme is #isrf_c_hyp_scheme_fixed_fraction (see that
   * parameter's own mutual-exclusion check). 0 (default): disabled. A
   * positive value removes the per-particle speed spread that
   * #ISRF_c_hyp_margin alone cannot (two neighbours on different time bins or
   * with different `h` otherwise get different `c_hyp_i`); the resulting
   * receiver-side CFL violation is instead prevented up front by a dedicated
   * timestep term (see #ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging).
   * Mutually exclusive with #ISRF_c_hyp_pin_for_debugging
   * (feedback_props_init() errors if both are set): the pin exists to probe a
   * single chosen `c_hyp`, this exists to ship a uniform one with its own
   * stability term, and stacking them would silently let one override the
   * other. */
  float ISRF_c_hyp_fixed_fraction_of_c;

  /*! Debug/test-only: with #ISRF_c_hyp_fixed_fraction_of_c positive, skip
   * the radiation timestep term (`C_hyp*h_i/(f*c)`) that keeps that
   * particle's own receiver-side CFL condition satisfied. 0 (default): the
   * term is applied, which is the only supported configuration whenever
   * the fixed fraction is on. 1: the term is skipped, so
   * #ISRF_c_hyp_fixed_fraction_of_c is exactly as unstable at a seam as the
   * shipped per-particle formula -- this exists solely to measure, by A/B run,
   * how much of the fixed fraction's step-count cost the timestep term itself
   * is responsible for. Never set in a production run. No effect when
   * #ISRF_c_hyp_fixed_fraction_of_c is 0. */
  char ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging;

  /*! Ceiling of the triggered artificial-conductivity coefficient. On by
   * default at the calibrated ceiling; 0 disables the term (for A/B runs).
   * The
   * enforced range depends on #ISRF_c_hyp_margin (see
   * feedback_properties_init()'s own range check), so raising
   * #ISRF_c_hyp_margin can require lowering this. */
  float ISRF_dissipation_alpha_max;

  /*! Undershoot of a particle's own `rho_prev*u` below the neighbours'
   * kernel-mean `|rho_prev*u_prev|`, relative, at which the
   * negativity-triggered dissipation coefficient reaches
   * #ISRF_dissipation_alpha_max. */
  float ISRF_dissipation_negativity_threshold;

  /*! Floor under the negativity trigger, `h/lambda`-gated: the trigger fires
   * only on negativity and is exactly zero on the positive delta-shell
   * front of an optically-thin P1 pulse, so a purely reactive coefficient
   * cannot damp the resulting dispersive wake there. This floor supplies
   * dissipation the trigger structurally cannot. Combined
   * with #ISRF_dissipation_floor_h_over_lambda as
   * `alpha_floor/(1+(h*kappa/eps_lambda)^4)`, then taken as a max against
   * the trigger's own output. 0 disables the floor and recovers the
   * trigger-only behaviour exactly. */
  float ISRF_dissipation_alpha_floor;

  /*! Screening-length error budget (`eps_lambda`) gating where the floor
   * (#ISRF_dissipation_alpha_floor) applies: the floor rolls off as
   * `(eps_lambda/(h*kappa))^4` once `h/lambda` exceeds this value, since
   * the joint stability bound's own steady-state distortion,
   * `lambda_eff/lambda <= sqrt(1+0.27*alpha_floor*eps_lambda)`, is bounded
   * only near `h/lambda ~ eps_lambda`; away from it the roll-off keeps the
   * floor negligible where physical absorption or the trigger's own decay
   * memory already dominates. */
  float ISRF_dissipation_floor_h_over_lambda;

  /*! Flux-relaxation residual threshold gating the floor's own aim:
   * `R = |F + C*grad_u| / (|F| + C*|grad_u|)`, `C =
   * c_hyp^2/(c_hyp*kappa+H)` the exact fixed point of the UNLIMITED
   * flux-update recurrence, is identically 0 at the code's own discrete
   * steady state and ~1 away from it; the triangle inequality bounds R to
   * [0, 1], so this parameter's valid range is [0, 1] too (checked in
   * feedback_props_init(): a value above 1 would scale the floor down on
   * every particle, fronts included). A particle whose flux is instead
   * pinned by the M1 limiter (#radiation_apply_flux_limiter_band, `|F| =
   * c_M*u`, the free-streaming branch) generally never reaches that fixed
   * point either, so R stays finite there too: a conservative false
   * positive that keeps part of the floor where the limiter is active,
   * never removes protection where a front is present. The floor's aim is
   * multiplied by `min(1, (R/eps_R)^2)`, so a particle already at the
   * flux-relaxation fixed point (a resolved, settled profile) gets a
   * reduced floor, while a genuine front (R large) keeps the floor at
   * full strength. Can only lower the floor relative to
   * #ISRF_dissipation_floor_h_over_lambda's own roll-off, never raise
   * it: `s=1` whenever exactly one of `F`, `grad_u` is zero (`R=1`), so a
   * fresh front or a limiter-zeroed flux keeps the full floor, provided
   * the relaxation weight `w = kappa + H/c_hyp` is nonzero. The exception
   * is a quiescent particle with both `F` and `grad_u` zero and `w`
   * nonzero: that is trivially at the fixed point, giving `R=0` and
   * `s=0`. At `w=0` (`kappa=0` and `H=0`: no relaxation timescale to
   * settle against), the gate returns `s=1` unconditionally instead. 0
   * disables this gate (R treated as always saturating, i.e. `s=1`
   * everywhere) and recovers the `h/lambda`-only floor exactly. */
  float ISRF_dissipation_floor_relaxation_residual;

  /*! Debug/test-only: bypass the negativity trigger and hold every
   * particle's dissipation coefficient (both bands) at this fixed value,
   * whenever positive. Not itself subject to the joint
   * (ISRF_dissipation_alpha_max, ISRF_c_hyp_margin) bound; only a
   * warning fires if it exceeds that bound. 0 (default): disabled, use the
   * trigger. Never set in a production run. */
  float ISRF_dissipation_alpha_pin_for_debugging;

  /*! Minimal density to consider a particle eligible for HII ionization */
  float HII_min_density;

  /*! HII region rebuild frequency */
  float HII_rebuild_time;

  /*! Maximun age of star particle to trigger the HII region algorithm */
  float HII_max_age;

  /*! How to treat the boundary gas particle a star's remaining photon
   * budget cannot fully ionize: 0 = probabilistic (weighted coin flip,
   * unbiased in expectation), 1 = deterministic (always ionize it,
   * letting the budget go slightly negative). */
  char HII_deterministic_boundary_ionization;

  /*! Floors the elapsed interval the per-pass ionizing photon budget is
   * integrated over, in every cadence mode. */
  float HII_rebuild_floor_Myr;

  /* ------------- Stellar winds properties ------------- */

  /*! Pre-supernova feedback energy effectively deposited */
  float winds_efficiency;

  /*! Do stellar wind feedback? */
  char with_stellar_wind_feedback;
};

/**
 * @brief Does this run need Grackle's chemistry_data actually resolved
 * (cooling_init() having run, via --cooling or --temperature), for
 * GEAR's own local Lyman-Werner/FUV photoelectric-heating channel to
 * read a real (not silently zero) local_dust_to_gas_ratio? See
 * radiation_isrf.c's dust-opacity helpers.
 *
 * @param feedback_props The #feedback_props.
 * @return True if with_photoelectric_heating is enabled.
 */
__attribute__((always_inline)) INLINE static int
feedback_props_needs_cooling_initialized(
    const struct feedback_props *feedback_props) {
  return (feedback_props->radiation_policy &
          radiation_policy_photoelectric_heating) != 0;
}

/**
 * @brief Print the feedback model.
 *
 * @param feedback_props The #feedback_props
 */
__attribute__((always_inline)) INLINE static void feedback_props_print(
    const struct feedback_props *feedback_props) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  /* Print the name of the elements */
  char txt[GEAR_CHEMISTRY_ELEMENT_COUNT * (GEAR_LABELS_SIZE + 2)] = "";
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (i != 0) {
      strcat(txt, ", ");
    }
    strcat(txt, stellar_evolution_get_element_name(
                    &feedback_props->stellar_model, i));
  }

  if (engine_rank == 0) {
    message("Chemistry elements: %s", txt);
  }

  /* Print the feedback properties, grouped by mechanism to match
   * GEARFeedback's layout in parameter_example.yml: stellar evolution, SN,
   * stellar winds, radiation pressure, HII regions, ISRF. */

  /* Stellar evolution */
  message("Yields table                                               = %s",
          feedback_props->stellar_model.yields_table);
  message("dt_evolution factor_max                                    = %g",
          feedback_props->dt_evolution_factor_max);
  message("event_dt_floor (internal units)                            = %g",
          feedback_props->event_dt_floor_Myr);

  /* Print the stellar model */
  stellar_model_print(&feedback_props->stellar_model);

  /* Print the first stars */
  if (feedback_props->metallicity_max_first_stars != -1) {
    message("Yields table first stars                                 = %s",
            feedback_props->stellar_model_first_stars.yields_table);
    stellar_model_print(&feedback_props->stellar_model_first_stars);
    message("Metallicity max for the first stars (in abundance)       = %g",
            feedback_props->imf_transition_metallicity);
    message("Metallicity max for the first stars (in mass fraction)   = %g",
            feedback_props->metallicity_max_first_stars);
  }

  /* Supernovae */
  message("Supernovae efficiency                                      = %.2g",
          feedback_props->supernovae_efficiency);

  /* Stellar winds */
  message("Stellar wind feedback                                      = %s",
          feedback_props->with_stellar_wind_feedback ? "ON" : "OFF");
  message("Stellar winds efficiency                                   = %.2g",
          feedback_props->winds_efficiency);

  /* Radiation pressure */
  message(
      "Radiation pressure                                         = %i",
      feedback_props->radiation_policy & radiation_policy_radiation_pressure);
  message("Radiation pressure efficiency                              = %.2g",
          feedback_props->radiation_pressure_efficiency);

  /* HII regions */
  const char do_photoionization =
      feedback_props->radiation_policy & radiation_policy_photoionization;
  message("Photoionization                                            = %i",
          do_photoionization);

  if (do_photoionization) {
    message("HII region minimal gas density (internal units)            = %g",
            feedback_props->HII_min_density);
    message("HII boundary ionization mode                               = %s",
            feedback_props->HII_deterministic_boundary_ionization
                ? "deterministic"
                : "probabilistic");
    message("HII max age (internal units)                               = %g",
            feedback_props->HII_max_age);
    message("HII rebuild time (internal units)                          = %g",
            feedback_props->HII_rebuild_time);
    message("HII rebuild floor (internal units)                         = %g",
            feedback_props->HII_rebuild_floor_Myr);
  }

  /* ISRF */
  const char do_photoelectric_heating =
      feedback_props->radiation_policy & radiation_policy_photoelectric_heating;
  message("Photo-electric heating / H2 photodissociation (ISRF)       = %i",
          do_photoelectric_heating);
  if (do_photoelectric_heating) {
    message("ISRF propagation                                           = %s",
            feedback_props->ISRF_propagation ? "ON" : "OFF (injection only)");
    if (feedback_props->ISRF_propagation) {
      message("ISRF propagation speed margin (C_hyp)                      = %g",
              feedback_props->ISRF_c_hyp_margin);
      const char *isrf_c_hyp_scheme_name = "shipped (per-particle timestep)";
      if (feedback_props->ISRF_c_hyp_scheme == isrf_c_hyp_scheme_kernel_local)
        isrf_c_hyp_scheme_name = "kernel-local (per-kernel slowest timestep)";
      else if (feedback_props->ISRF_c_hyp_scheme ==
               isrf_c_hyp_scheme_fixed_fraction)
        isrf_c_hyp_scheme_name = "fixed fraction of c";
      else if (feedback_props->ISRF_c_hyp_scheme ==
               isrf_c_hyp_scheme_consistent_variable_c)
        isrf_c_hyp_scheme_name =
            "consistent variable-c operators (shipped c_hyp formula)";
      else if (feedback_props->ISRF_c_hyp_scheme ==
               isrf_c_hyp_scheme_kernel_local_plus_variable_c)
        isrf_c_hyp_scheme_name =
            "consistent variable-c operators (kernel-local c_hyp formula)";
      message("ISRF c_hyp scheme                                          = %s",
              isrf_c_hyp_scheme_name);
      if (feedback_props->ISRF_c_hyp_pin_for_debugging > 0.f)
        message(
            "ISRF c_hyp pinned for debugging (physical units)          = %g",
            feedback_props->ISRF_c_hyp_pin_for_debugging);
      if (feedback_props->ISRF_c_hyp_fixed_fraction_of_c > 0.f) {
        message(
            "ISRF c_hyp fixed fraction of c                             = %g",
            feedback_props->ISRF_c_hyp_fixed_fraction_of_c);
        message(
            "ISRF c_hyp fixed fraction timestep term                    = %s",
            feedback_props->ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging
                ? "OFF for debugging (diagnostic only, not a supported "
                  "configuration)"
                : "ON");
      }
      message("ISRF dissipation alpha_max                                 = %g",
              feedback_props->ISRF_dissipation_alpha_max);
      message("ISRF dissipation negativity threshold                      = %g",
              feedback_props->ISRF_dissipation_negativity_threshold);
      message("ISRF dissipation alpha_floor                               = %g",
              feedback_props->ISRF_dissipation_alpha_floor);
      message("ISRF dissipation floor h/lambda budget (eps_lambda)        = %g",
              feedback_props->ISRF_dissipation_floor_h_over_lambda);
      message("ISRF dissipation floor relaxation-residual gate (eps_R)    = %g",
              feedback_props->ISRF_dissipation_floor_relaxation_residual);
      if (feedback_props->ISRF_dissipation_alpha_pin_for_debugging > 0.f)
        message(
            "ISRF dissipation alpha pinned for debugging                = %g",
            feedback_props->ISRF_dissipation_alpha_pin_for_debugging);
    }
  }
}

/**
 * @brief Enforce that #feedback_props.ISRF_c_hyp_scheme and
 * #feedback_props.ISRF_c_hyp_fixed_fraction_of_c stay a matched pair: the
 * kernel-local and fixed-fraction c_hyp schemes are alternatives, not
 * layers, so a magnitude set without its scheme selected (or a scheme
 * selected without its magnitude) would otherwise silently do nothing or
 * silently pick up a stale value. A standalone function (not inlined into
 * feedback_props_init()'s own body) so a unit test can call it directly,
 * with neither a #swift_params nor the stellar-evolution tables
 * feedback_props_init() also reads.
 *
 * @param scheme #feedback_props.ISRF_c_hyp_scheme's parsed value.
 * @param fixed_fraction #feedback_props.ISRF_c_hyp_fixed_fraction_of_c's
 * parsed value.
 */
__attribute__((always_inline)) INLINE static void
feedback_props_check_c_hyp_scheme(int scheme, float fixed_fraction) {
  if (scheme == isrf_c_hyp_scheme_fixed_fraction && fixed_fraction <= 0.f)
    error(
        "GEARFeedback:ISRF_c_hyp_scheme is set to 2 (fixed fraction of c) "
        "but GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c is 0: there is "
        "nothing to select. Set a positive fraction.");
  if (scheme != isrf_c_hyp_scheme_fixed_fraction && fixed_fraction > 0.f)
    error(
        "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c is set (%g) but "
        "GEARFeedback:ISRF_c_hyp_scheme is %d, not 2 (fixed fraction of c): "
        "the two speed schemes are alternatives, not layers. Set "
        "ISRF_c_hyp_scheme to 2 to use this fraction, or leave it at 0 if "
        "it was set by mistake.",
        fixed_fraction, scheme);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void feedback_props_init(
    struct feedback_props *fp, const struct phys_const *phys_const,
    const struct unit_system *us, struct swift_params *params,
    const struct hydro_props *hydro_props, const struct cosmology *cosmo) {

  /* Several fields below (ISRF_*, HII_*, stellar_model_first_stars) are only
   * assigned inside a conditional branch; zero first so the disabled case
   * reads a defined 0, not uninitialised stack memory, since those fields
   * are read downstream unconditionally. sink_props_init() writes
   * with_sinks after this function returns, not before. */
  bzero(fp, sizeof(struct feedback_props));

  /* Supernovae energy efficiency */
  double e_efficiency =
      parser_get_param_double(params, "GEARFeedback:supernovae_efficiency");
  fp->supernovae_efficiency = e_efficiency;

  /* Activate the stellar wind feedback */
  char with_stellar_wind_feedback = (char)parser_get_param_int(
      params, "GEARFeedback:with_stellar_wind_feedback");
  fp->with_stellar_wind_feedback = with_stellar_wind_feedback;

  /* Are we running with photoionization? Read early, like
   * with_stellar_wind_feedback above, so stellar_evolution_props_init() can
   * skip opening the radiation table for non-radiation runs. */
  const char with_photoionization = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_photoionization", 0);

  /* Radiation pressure efficiency, read early for the same reason: its sign
   * decides whether the radiation table is needed too, independently of
   * with_photoionization. */
  const float radiation_pressure_efficiency = parser_get_opt_param_float(
      params, "GEARFeedback:radiation_pressure_efficiency", 0.0);

  /* Are we running with the local Lyman-Werner/FUV feedback (photoelectric
   * heating + H2 photodissociation)? Read early, for the same reason as
   * with_photoionization: it needs both the radiation table (for L_bol,
   * from which L_FUV/L_LW are split via Teff) and the Teff table itself. */
  const char with_photoelectric_heating = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_photoelectric_heating", 0);

  /* The radiation table backs the HII photoionization band, the bolometric
   * radiation-pressure band, and (since Teff is stored on the same table)
   * the Lyman-Werner/FUV band split (see
   * stellar_evolution_compute_preSN_feedback_individual_star()/_spart()).
   * Previously omitted radiation_policy_photoelectric_heating here because
   * "photoelectric heating has no downstream consumer yet". Now that it
   * does, leaving it out would silently skip opening the radiation table
   * (and its Teff dataset) for a with_photoelectric_heating-only run,
   * and desync from feedback_struct_restore()'s own copy of this same
   * condition on restart (see that function's matching comment). */
  const char with_radiation = with_photoionization ||
                              (radiation_pressure_efficiency > 0.0f) ||
                              with_photoelectric_heating;

  /* Pre-Supernovae energy efficiency */
  double w_efficiency = 0.0;
  if (with_stellar_wind_feedback) {
    w_efficiency = parser_get_param_double(
        params, "GEARFeedback:stellar_winds_efficiency");
  }

  fp->winds_efficiency = w_efficiency;

  /* filename of the chemistry tables. */
  parser_get_param_string(params, "GEARFeedback:yields_table",
                          fp->stellar_model.yields_table);

  /* Initialize the stellar models. */
  stellar_evolution_props_init(&fp->stellar_model, phys_const, us, params,
                               cosmo, fp->with_stellar_wind_feedback,
                               with_radiation);

  /* Read the metallicity threashold */
  fp->imf_transition_metallicity = parser_get_opt_param_float(
      params, "GEARFeedback:imf_transition_metallicity", 0);

  /* Read and get the solar abundances */
  struct chemistry_global_data data;
  bzero(&data, sizeof(struct chemistry_global_data));
  chemistry_read_solar_abundances(params, &data);

  const int iFe = stellar_evolution_get_element_index(&fp->stellar_model, "Fe");
  const float XFe = data.solar_abundances[iFe];

  if (fp->imf_transition_metallicity == 0)
    fp->metallicity_max_first_stars = -1;
  else
    fp->metallicity_max_first_stars =
        exp10(fp->imf_transition_metallicity) * XFe;

  /* Now initialize the first stars. */
  if (fp->metallicity_max_first_stars == -1) {
    message("First stars are disabled.");
  } else {
    if (fp->metallicity_max_first_stars < 0) {
      error(
          "The metallicity threshold for the first stars is in mass fraction. "
          "It cannot be lower than 0.");
    }
    if (engine_rank == 0) {
      message("Reading the stellar model for the first stars");
    }
    parser_get_param_string(params, "GEARFeedback:yields_table_first_stars",
                            fp->stellar_model_first_stars.yields_table);
    stellar_evolution_props_init(&fp->stellar_model_first_stars, phys_const, us,
                                 params, cosmo, fp->with_stellar_wind_feedback,
                                 with_radiation);
  }

  /* ------------- Star evolution timestep properties ------------- */
  /* Runs for every star regardless of radiation, so parsed unconditionally
   * (feedback_compute_spart_timestep() uses this factor for all stars). */

  /* Needed unconditionally below (event_dt_floor_Myr's conversion), unlike
   * HII_max_age/HII_rebuild_time/HII_rebuild_floor_Myr further down, which
   * stay gated behind with_photoionization. */
  const double Myr_internal_units = 1e6 * phys_const->const_year;

  fp->dt_evolution_factor_max =
      parser_get_opt_param_float(params, "GEARFeedback:dt_evolution_factor_max",
                                 default_dt_evolution_factor_max);

  if (fp->dt_evolution_factor_max < 1.f)
    error("GEARFeedback:dt_evolution_factor_max must be >= 1 (got %g).",
          fp->dt_evolution_factor_max);

  fp->event_dt_floor_Myr = parser_get_opt_param_float(
      params, "GEARFeedback:event_dt_floor_Myr", default_event_dt_floor_Myr);

  if (fp->event_dt_floor_Myr <= 0.f)
    error(
        "GEARFeedback:event_dt_floor_Myr must be > 0 (got %g): it floors "
        "single_star's event-anchored death timestep (dt_event) in every "
        "run, not just radiation ones -- <= 0 reopens the get_spart_timestep "
        "crash this parameter exists to prevent.",
        fp->event_dt_floor_Myr);

  fp->event_dt_floor_Myr *= Myr_internal_units;

  /* Startup cross-check: event_dt_floor_Myr is the sole crash guard for
   * every single_star death (get_spart_timestep()'s dt_min error()), so it
   * must itself clear dt_min. TimeIntegration:dt_min is a mandatory
   * parameter already fully parsed into `params` at this point (engine_init
   * has not run yet, but that only matters for e->dt_min the struct field --
   * the parsed value is available directly from `params` regardless of
   * init order). */
  const double dt_min =
      parser_get_param_double(params, "TimeIntegration:dt_min");
  if (fp->event_dt_floor_Myr <= dt_min)
    error(
        "GEARFeedback:event_dt_floor_Myr (%g, internal units) must exceed "
        "TimeIntegration:dt_min (%g, internal units): otherwise the floor "
        "itself can still trip get_spart_timestep()'s dt_min error() at a "
        "single_star's death.",
        fp->event_dt_floor_Myr, dt_min);

  /* ------------- Subgrid Radiation properties ------------- */
  fp->radiation_policy = 0;

  /* Reject the pre-rename GEARFeedback:LW_FUV_* keys: an unrecognised key
   * is otherwise only reported as unused, so an old parameter file would
   * silently lose these settings instead of failing loudly. A restart
   * bypasses this check by design: feedback_struct_restore() re-reads the
   * saved struct, not the parameter file. */
  const char *const deprecated_ISRF_keys[9] = {
      "GEARFeedback:LW_FUV_propagation",
      "GEARFeedback:LW_FUV_c_hyp_margin",
      "GEARFeedback:LW_FUV_c_hyp_pin_for_debugging",
      "GEARFeedback:LW_FUV_dissipation_alpha_max",
      "GEARFeedback:LW_FUV_dissipation_negativity_threshold",
      "GEARFeedback:LW_FUV_dissipation_alpha_floor",
      "GEARFeedback:LW_FUV_dissipation_floor_h_over_lambda",
      "GEARFeedback:LW_FUV_dissipation_floor_relaxation_residual",
      "GEARFeedback:LW_FUV_dissipation_alpha_pin_for_debugging"};
  const char *const renamed_ISRF_keys[9] = {
      "GEARFeedback:ISRF_propagation",
      "GEARFeedback:ISRF_c_hyp_margin",
      "GEARFeedback:ISRF_c_hyp_pin_for_debugging",
      "GEARFeedback:ISRF_dissipation_alpha_max",
      "GEARFeedback:ISRF_dissipation_negativity_threshold",
      "GEARFeedback:ISRF_dissipation_alpha_floor",
      "GEARFeedback:ISRF_dissipation_floor_h_over_lambda",
      "GEARFeedback:ISRF_dissipation_floor_relaxation_residual",
      "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging"};
  for (int i = 0; i < 9; ++i) {
    if (parser_does_param_exist(params, deprecated_ISRF_keys[i]))
      error("%s has been renamed to %s. Update the parameter file.",
            deprecated_ISRF_keys[i], renamed_ISRF_keys[i]);
  }

  /* TODO: For the future, enforce these to have a non-zero value */

  /* Radiation pressure */
  fp->radiation_pressure_efficiency = radiation_pressure_efficiency;

  if (fp->radiation_pressure_efficiency > 0.0) {
    fp->radiation_policy |= radiation_policy_radiation_pressure;
  }

  /* Parsed unconditionally, so the LW/FUV injection never reads an unset
   * path. */
  char extinction_path[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(params, "GEARFeedback:ISRF_extinction_path",
                              extinction_path, "kernel_diameter");
  if (strcmp(extinction_path, "kernel_diameter") == 0)
    fp->ISRF_extinction_path_in_kernel_radii = 2.0f;
  else if (strcmp(extinction_path, "kernel_radius") == 0)
    fp->ISRF_extinction_path_in_kernel_radii = 1.0f;
  else
    error(
        "GEARFeedback:ISRF_extinction_path must be kernel_diameter or "
        "kernel_radius, got '%s'.",
        extinction_path);

  if (with_photoelectric_heating) {
    fp->radiation_policy |= radiation_policy_photoelectric_heating;

    fp->ISRF_propagation = (char)parser_get_opt_param_int(
        params, "GEARFeedback:ISRF_propagation", 0);

    /* Which c_hyp scheme runs; see #isrf_c_hyp_scheme's own doxygen.
     * Parsed unconditionally, like the pin/fraction below, so a validation
     * run can set it even with ISRF_propagation off in the base config.
     * Default is scheme 4: best measured negativity and amplitude of the
     * five; the previously shipped scheme 0 stays reachable by setting this
     * parameter explicitly. */
    fp->ISRF_c_hyp_scheme = parser_get_opt_param_int(
        params, "GEARFeedback:ISRF_c_hyp_scheme",
        isrf_c_hyp_scheme_kernel_local_plus_variable_c);
    if (fp->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_shipped &&
        fp->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_kernel_local &&
        fp->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_fixed_fraction &&
        fp->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_consistent_variable_c &&
        fp->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_kernel_local_plus_variable_c)
      error(
          "GEARFeedback:ISRF_c_hyp_scheme must be 0 (shipped), 1 "
          "(kernel-local), 2 (fixed fraction of c), 3 (consistent "
          "variable-c operators) or 4 (kernel-local speed with the "
          "consistent variable-c operators) (got %d).",
          fp->ISRF_c_hyp_scheme);
    /* Carries the selection into radiation_propagation_iact.h's pairwise
     * dispatch, which has no #engine pointer to read #feedback_props from;
     * see #isrf_c_hyp_consistent_variable_c's own doxygen. Scheme 4 selects
     * the same operator rewrite as scheme 3, just fed by the kernel-local
     * speed instead of the shipped one -- the two axes are independent. */
    isrf_c_hyp_consistent_variable_c =
        (fp->ISRF_c_hyp_scheme == isrf_c_hyp_scheme_consistent_variable_c ||
         fp->ISRF_c_hyp_scheme ==
             isrf_c_hyp_scheme_kernel_local_plus_variable_c);

    /* Debug/test-only: see ISRF_c_hyp_pin_for_debugging's own doxygen.
     * Parsed unconditionally (like the stability margin and dissipation
     * parameters below) so a validation run can set it even with
     * ISRF_propagation off in the base config and toggled on
     * separately; only meaningful when it is. */
    fp->ISRF_c_hyp_pin_for_debugging = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_c_hyp_pin_for_debugging", 0.0f);

    /* Uniform reduced light-speed candidate: 0 disables it. Only takes
     * effect under ISRF_c_hyp_scheme == isrf_c_hyp_scheme_fixed_fraction
     * (enforced below). Parsed and validated unconditionally, like the pin
     * above, so a validation run can set it even with ISRF_propagation off
     * in the base config. */
    fp->ISRF_c_hyp_fixed_fraction_of_c = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c", 0.0f);
    fp->ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging =
        (char)parser_get_opt_param_int(
            params,
            "GEARFeedback:ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging",
            0);

    if (fp->ISRF_c_hyp_fixed_fraction_of_c < 0.f ||
        fp->ISRF_c_hyp_fixed_fraction_of_c > 1.f)
      error(
          "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c must lie in [0, 1] "
          "(got %g): it is a fraction of the true speed of light c. 0 "
          "disables it.",
          fp->ISRF_c_hyp_fixed_fraction_of_c);

    /* See feedback_props_check_c_hyp_scheme()'s own doxygen: the two speed
     * schemes are alternatives, not layers. */
    feedback_props_check_c_hyp_scheme(fp->ISRF_c_hyp_scheme,
                                      fp->ISRF_c_hyp_fixed_fraction_of_c);

    if (fp->ISRF_c_hyp_fixed_fraction_of_c > 0.f &&
        fp->ISRF_c_hyp_pin_for_debugging > 0.f)
      error(
          "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c and "
          "GEARFeedback:ISRF_c_hyp_pin_for_debugging cannot both be set: "
          "they are two different ways to override c_hyp_i and stacking "
          "them leaves it ambiguous which one actually took effect.");

    /* Parsed and validated unconditionally, like the pin above: a
     * validation run can set and check the stability margin (and the
     * dissipation coefficients below) even with ISRF_propagation off. */
    fp->ISRF_c_hyp_margin = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_c_hyp_margin", 0.5f);

    /* Absolute static bound: even with dissipation fully disabled
     * (alpha_max = alpha_floor = 0), the joint stability bound below
     * (6.2*alpha*C_hyp + 0.70*C_hyp^2 <= 2) still applies at alpha = 0,
     * giving C_hyp <= sqrt(2/0.70) ~ 1.6903. Increasing alpha only
     * tightens the bound further (enforced separately below, once
     * alpha_max/alpha_floor are known), so alpha = 0 is the loosest case
     * and no configuration can ever exceed this value. The previous 1.7
     * ceiling admitted C_hyp values in (1.6903, 1.7] that are unstable
     * even with dissipation off, since the joint check below only fires
     * when alpha_max or alpha_floor is nonzero. */
    const float ISRF_c_hyp_absolute_bound = sqrtf(2.f / 0.70f);

    if (fp->ISRF_c_hyp_margin <= 0.f ||
        fp->ISRF_c_hyp_margin > ISRF_c_hyp_absolute_bound)
      error(
          "GEARFeedback:ISRF_c_hyp_margin must lie in "
          "(0, %g] (got %g): above this bound the staggered "
          "exact-relaxation scheme is no longer stable for every "
          "lambda/h at this project's kernel/eta_neighbours choice, even "
          "with dissipation (ISRF_dissipation_alpha_max/alpha_floor) "
          "fully disabled (6.2*alpha*C_hyp + 0.70*C_hyp^2 <= 2 at "
          "alpha = 0).",
          ISRF_c_hyp_absolute_bound, fp->ISRF_c_hyp_margin);

    /* Negativity-triggered artificial dissipation. Shipped default 0.5,
     * the same value as the floor's own ceiling, so the joint stability
     * bound checked below
     * is unchanged. Parsed and validated unconditionally,
     * like the pin and the stability margin above, so a validation run can
     * exercise these even with ISRF_propagation off. */
    fp->ISRF_dissipation_alpha_max = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_alpha_max", 0.5f);
    fp->ISRF_dissipation_negativity_threshold = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_negativity_threshold", 0.01f);
    fp->ISRF_dissipation_alpha_pin_for_debugging = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging", 0.0f);

    if (fp->ISRF_dissipation_alpha_max < 0.f)
      error("GEARFeedback:ISRF_dissipation_alpha_max must be >= 0 (got %g).",
            fp->ISRF_dissipation_alpha_max);

    if (fp->ISRF_dissipation_negativity_threshold <= 0.f ||
        fp->ISRF_dissipation_negativity_threshold > 1.f)
      error(
          "GEARFeedback:ISRF_dissipation_negativity_threshold must lie "
          "in (0, 1] (got %g).",
          fp->ISRF_dissipation_negativity_threshold);

    /* Diffuse-phase floor under the trigger: shipped defaults 0.5/0.5. The
     * eps_lambda default moved 0.05 -> 0.5 together with the roll-off
     * exponent 2 -> 4: the quartic tail is what keeps the thick regime
     * negligible, so the knee itself no longer has to sit an order of
     * magnitude below the regimes that need the floor. */
    fp->ISRF_dissipation_alpha_floor = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_alpha_floor", 0.5f);
    fp->ISRF_dissipation_floor_h_over_lambda = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_floor_h_over_lambda", 0.5f);

    if (fp->ISRF_dissipation_alpha_floor < 0.f)
      error(
          "GEARFeedback:ISRF_dissipation_alpha_floor must be >= 0 (got "
          "%g).",
          fp->ISRF_dissipation_alpha_floor);

    if (fp->ISRF_dissipation_floor_h_over_lambda <= 0.f)
      error(
          "GEARFeedback:ISRF_dissipation_floor_h_over_lambda must be > 0 "
          "(got %g).",
          fp->ISRF_dissipation_floor_h_over_lambda);

    /* Flux-relaxation residual gate: 0 disables it (recovers the
     * h/lambda-only floor). */
    fp->ISRF_dissipation_floor_relaxation_residual = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_floor_relaxation_residual",
        0.40f);

    if (fp->ISRF_dissipation_floor_relaxation_residual < 0.f ||
        fp->ISRF_dissipation_floor_relaxation_residual > 1.f)
      error(
          "GEARFeedback:ISRF_dissipation_floor_relaxation_residual must "
          "lie in [0, 1] (got %g): the residual R it gates is itself bounded "
          "to [0, 1] by the triangle inequality, so a value above 1 would "
          "scale the floor down on every particle, fronts included. 0 "
          "disables the gate.",
          fp->ISRF_dissipation_floor_relaxation_residual);

    /* Joint (alpha_max, C_hyp) stability bound: 6.2 = 2*I_W and
     * 0.70 = nu_max_coeff^2/2, the kernel's own Wendland-C2 lattice
     * constants (I_W = 3.10,
     * nu_max_coeff = 1.18), reproduced by
     * theory/GEAR/Radiation/verify_isrf_dissipation.py's Part C. The
     * floor can dissipate even where the trigger never fires (it is not
     * gated on negativity), so it must satisfy the same bound as the
     * trigger's own ceiling: check max(alpha_max, alpha_floor). */
    const float C_hyp = fp->ISRF_c_hyp_margin;
    const float alpha_bound = (2.f - 0.70f * C_hyp * C_hyp) / (6.2f * C_hyp);
    const float alpha_joint_ceiling =
        max(fp->ISRF_dissipation_alpha_max, fp->ISRF_dissipation_alpha_floor);
    if (alpha_joint_ceiling > 0.f && alpha_joint_ceiling > alpha_bound)
      error(
          "max(GEARFeedback:ISRF_dissipation_alpha_max, "
          "GEARFeedback:ISRF_dissipation_alpha_floor) = %g exceeds the "
          "stability bound %g at GEARFeedback:ISRF_c_hyp_margin = %g "
          "(6.2*alpha*C_hyp + 0.70*C_hyp^2 <= 2). alpha_max = %g, "
          "alpha_floor = %g. Lower GEARFeedback:ISRF_c_hyp_margin to "
          "raise the bound, or lower alpha_max/alpha_floor to fit it.",
          alpha_joint_ceiling, alpha_bound, C_hyp,
          fp->ISRF_dissipation_alpha_max, fp->ISRF_dissipation_alpha_floor);

    if (fp->ISRF_dissipation_alpha_pin_for_debugging < 0.f)
      error(
          "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging must be "
          ">= 0 (got %g).",
          fp->ISRF_dissipation_alpha_pin_for_debugging);

    if (fp->ISRF_dissipation_alpha_pin_for_debugging > 0.f &&
        fp->ISRF_dissipation_alpha_pin_for_debugging > alpha_bound)
      warning(
          "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging (%g) "
          "exceeds the stability bound %g at the run's C_hyp: not itself "
          "clamped (debug/test only). Never use this in a production run.",
          fp->ISRF_dissipation_alpha_pin_for_debugging, alpha_bound);

    if (fp->ISRF_propagation) {
      if (fp->ISRF_c_hyp_pin_for_debugging > 0.f)
        warning(
            "GEARFeedback:ISRF_c_hyp_pin_for_debugging is set: every "
            "particle's own hyperbolic propagation speed is pinned to %g "
            "(physical units) instead of C_hyp*h/dt_max. Never use this in "
            "a production run.",
            fp->ISRF_c_hyp_pin_for_debugging);

      /* Tripwire, not a fix (see radiation_get_dust_mass_opacity() in
       * radiation_isrf.c): IC metallicity is per-particle HDF5 data, not
       * visible here, so this warns unconditionally rather than gating on
       * a metallicity value that cannot bound the risk (kappa -> 0
       * continuously as Z -> 0, with no floor on the opacity itself). */
      warning(
          "GEARFeedback:ISRF_propagation is on together with "
          "GEARFeedback:with_photoelectric_heating. The propagation's only "
          "loss channel is dust absorption, whose rate is proportional to "
          "the gas metallicity. Gas at or near zero metallicity has no "
          "loss channel. Check GEARChemistry:initial_metallicity and any "
          "MetalMassFraction field in the initial conditions before a long "
          "run.");
    }
  }

  if (with_photoionization) {
    fp->radiation_policy |= radiation_policy_photoionization;

    /* Read the minimal density */
    fp->HII_min_density =
        parser_get_opt_param_float(params, "GEARFeedback:HII_min_density_Hpcm3",
                                   default_HII_min_density_Hpcm3);

    /* Read the HII region maximal age */
    fp->HII_max_age = parser_get_opt_param_float(
        params, "GEARFeedback:HII_max_age_Myr", default_HII_max_age_Myr);

    /* Read the HII region rebuild frequency */
    fp->HII_rebuild_time =
        parser_get_opt_param_float(params, "GEARFeedback:HII_rebuild_time_Myr",
                                   default_HII_rebuild_time_Myr);

    /* Read the boundary-particle ionization mode */
    fp->HII_deterministic_boundary_ionization = parser_get_opt_param_int(
        params, "GEARFeedback:HII_deterministic_boundary_ionization",
        default_HII_deterministic_boundary_ionization);

    fp->HII_rebuild_floor_Myr =
        parser_get_opt_param_float(params, "GEARFeedback:HII_rebuild_floor_Myr",
                                   default_HII_rebuild_floor_Myr);
    if (fp->HII_rebuild_floor_Myr <= 0.f)
      error(
          "GEARFeedback:HII_rebuild_floor_Myr must be > 0 (got %g): it "
          "floors the interval every per-pass ionizing photon budget is "
          "integrated over, in every cadence mode -- <= 0 silently zeroes "
          "every star's first-pass budget.",
          fp->HII_rebuild_floor_Myr);

    /* Convert to internal units */
    const double m_p_cgs = phys_const->const_proton_mass *
                           units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    fp->HII_min_density *=
        m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

    /* Myr_internal_units already computed above, unconditionally. */
    fp->HII_max_age *= Myr_internal_units;
    fp->HII_rebuild_time *= Myr_internal_units;
    fp->HII_rebuild_floor_Myr *= Myr_internal_units;
  }

  /* -------------------------------------------- */
  /* Print the stellar properties */
  feedback_props_print(fp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Stellar feedback initialized");
  }
}

#endif /* SWIFT_GEAR_FEEDBACK_PROPERTIES_H */
