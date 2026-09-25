/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_STRUCT_GEAR_H
#define SWIFT_FEEDBACK_STRUCT_GEAR_H

#include "chemistry_struct.h"
#include "timeline.h"

/*! Maximum number of HEALPix angular pixels the HII ionization budget can
    be split across (12*nside_max^2). Every star carries a fixed-size
    array of this length regardless of the run's actual
    GEARFeedback:HII_angular_nside, so this is a memory/generality
    trade-off, not a physics one. Set via
    ./configure --with-number-of-hii-angular-pixels=N (default 12, i.e.
    nside<=1) rather than hardcoded here, so builds that only ever need
    nside<=1 don't pay for a finer split they'll never request.
    radiation_init() (src/feedback/GEAR/radiation.c) errors clearly at
    startup if GEARFeedback:HII_angular_nside implies more pixels than
    this build was configured for. */
#ifndef HII_MAX_ANGULAR_PIXELS
#error "HII_MAX_ANGULAR_PIXELS should be defined by configure (config.h)"
#endif

/**
 * @brief The three moments the ISRF module transports: PE (6-11.2 eV) and
 * Lyman-Werner (11.2-13.6 eV) specific energy, and the Lyman-Werner band's
 * photon-number moment. Indexes #feedback_part_data.isrf_moment and
 * #feedback_spart_data.radiation.L_band.
 *
 * #ISRF_MOMENT_LW_PHOTON shares #ISRF_OPERATOR_LW with #ISRF_MOMENT_LW
 * (#radiation_isrf_moment_to_operator below): it carries the same field as
 * #ISRF_MOMENT_LW, at a fixed reference photon energy, so its own opacity,
 * M1 closure and dissipation state are never independent of LW's.
 */
enum radiation_isrf_moment {
  ISRF_MOMENT_PE = 0,
  ISRF_MOMENT_LW,
  ISRF_MOMENT_LW_PHOTON,
  ISRF_MOMENT_COUNT
};

/**
 * @brief The physical band operators (dust opacity, M1 closure,
 * dissipation coefficients) a moment is evaluated against. Indexes
 * #feedback_part_data.isrf_operator; a moment reaches its own operator
 * through #radiation_isrf_moment_to_operator.
 */
enum radiation_isrf_operator {
  ISRF_OPERATOR_PE = 0,
  ISRF_OPERATOR_LW,
  ISRF_OPERATOR_COUNT
};

/**
 * @brief Compile-time moment -> operator map: which #radiation_isrf_operator
 * each #radiation_isrf_moment reads its opacity/closure/dissipation
 * coefficients from. A moment can point at an operator an earlier moment
 * already owns instead of getting its own. The cosmological redshift
 * residual is NOT looked up through this map: every consumer evaluates it
 * once per particle, not per operator, at each moment's own update site,
 * so two moments sharing one operator can still carry different residuals.
 *
 * Declared unsized: for a sized array, sizeof(a)/sizeof(a[0]) equals the
 * declared count regardless of how many initialisers were written, so a
 * length assert against a sized array is vacuous. Unsized, the array's
 * size comes from the initialiser list itself, so a missing entry is a
 * compile error rather than a silent zero-fill.
 */
static const enum radiation_isrf_operator radiation_isrf_moment_to_operator[] =
    {ISRF_OPERATOR_PE, ISRF_OPERATOR_LW, ISRF_OPERATOR_LW};

_Static_assert(sizeof(radiation_isrf_moment_to_operator) /
                       sizeof(radiation_isrf_moment_to_operator[0]) ==
                   ISRF_MOMENT_COUNT,
               "radiation_isrf_moment_to_operator needs one entry per "
               "ISRF_MOMENT_COUNT.");

/**
 * @brief Compile-time operator -> owning-moment map, the reverse of
 * #radiation_isrf_moment_to_operator: which moment's state a WRITER to an
 * operator field must read, when several moments share that operator.
 *
 * Every site that writes an #feedback_isrf_operator_data field must loop
 * over operators and fetch its one owning moment's data through this map,
 * rather than looping over moments and writing through the forward map: a
 * loop bounded by #ISRF_OPERATOR_COUNT can physically only touch each
 * operator once per particle, so a later moment sharing an operator cannot
 * turn an assignment into last-writer-wins or an accumulation into a
 * double count.
 *
 * Each entry MUST be the FIRST (lowest-index) moment that maps to that
 * operator: a bare round trip through the forward map accepts any sharer,
 * not only the first, so it cannot tell a correct entry from a wrong one
 * when two moments share an operator. #feedback_check_isrf_operator_owner_map()
 * checks the stronger property this map actually needs.
 *
 * Declared unsized for the same reason as the forward map above.
 */
static const enum radiation_isrf_moment radiation_isrf_operator_owner[] = {
    ISRF_MOMENT_PE, ISRF_MOMENT_LW};

_Static_assert(sizeof(radiation_isrf_operator_owner) /
                       sizeof(radiation_isrf_operator_owner[0]) ==
                   ISRF_OPERATOR_COUNT,
               "radiation_isrf_operator_owner needs one entry per "
               "ISRF_OPERATOR_COUNT.");

/**
 * @brief Per-moment ISRF transport state carried by each hydro particle,
 * one instance per #radiation_isrf_moment in #feedback_part_data.isrf_moment.
 */
struct feedback_isrf_moment_data {

  /*! Local specific radiation field of each band, internal
      specific-energy units (per-unit-mass, like this codebase's own
      hydro `u`; NOT cgs, unlike
      #feedback_spart_data.radiation.mean_excess_photon_energy_HI). The
      LW band feeds Grackle's RT_H2_dissociation_rate (COOLING_GRACKLE_MODE
      > 1 only) separately from the PE band, since the two bands carry
      different dust opacities. An
      instantaneous field strength, not an accumulated dose: holds the
      illuminating star(s)' most recently computed contribution, summed
      across every star that touched this particle in the same step
      (#feedback_part_data.ISRF_last_touch_ti), then held unchanged until the
      next step any star touches it again. Never cleared by cooling: a reader
      gets whatever was last written, however long ago that was.

      `a`-SCALING: PHYSICAL and mass-specific, with no scale-factor exponent
      of its own beyond what UNIT_CONV_ENERGY_PER_UNIT_MASS already implies
      (snapshot output declares `0.f`, tracers_io.h, and is correct). Being
      per-unit-mass, it carries the volume part of cosmological dilution
      automatically through the physical gas density it is measured against;
      only the redshift residual `-H*u` is an explicit term, applied in
      radiation_isrf.c's #radiation_end_force_propagation. Every field
      feeding it is physical too: the pairwise operators in
      radiation_propagation_iact.h convert their comoving-coordinate
      estimates before accumulating. */
  float u;

  /*! Snapshot of #u taken once per step (feedback_reset_part,
      cell_drift.c), before the density loop's h-iterations begin. The
      density loop's kernel mean reads it, so that value does not drift
      across h-iterations, and the end-force update rebuilds #u from it, so
      that update is idempotent. Until that update, #u still equals it. */
  float u_prev;

  /*! Hyperbolic propagation state: the tracked specific flux moment,
      mass-specific like #u. Zeroed unconditionally at
      first init (no IC field proposed for it); relaxed every step in the
      extra ghost (radiation_isrf.c's exact-relaxation update), then limited
      there against #u. Read by neighbours in the gradient loop (the old
      value, before this cell's extra ghost) and in the force loop (the new
      value: every force task runs after the extra ghosts of both cells it
      pairs, and a foreign particle is received after its own). Written to
      snapshots as
      "PESpecificFluxes"/"LWSpecificFluxes" (tracers_io.h), following
      #u's own "PESpecificEnergy(ies)" convention; no IC input
      field exists, and one added later would be the singular
      "PESpecificFlux"/"LWSpecificFlux".

      `a`-SCALING: PHYSICAL and mass-specific, like #u, with no
      scale-factor exponent of its own, which is what the output field
      declares. Its own redshift residual `-H*F` is applied in
      radiation_isrf.c's #radiation_end_gradient_propagation. */
  float specific_flux[3];

  /*! `(1/rho) div(rho F)` accumulator of this step's relaxed flux, FORCE
      loop (radiation_propagation_iact.h), whose dispatch fires both sides of
      a pair whenever either kernel reaches, so the mirrored pair is never
      split at h_i != h_j. Consumed by #radiation_end_force_propagation.
      Zeroed once per step by #radiation_end_gradient_propagation, not at the
      drift, so a snapshot holds the last step's value. PHYSICAL: the force
      loop converts its comoving-coordinate estimate before accumulating, so
      the snapshot output's declared `0.f` exponent (tracers_io.h) is
      correct. */
  float div_specific_flux;

  /*! Negativity-triggered artificial-dissipation source term, FORCE loop
      (radiation_propagation_iact.h): pairwise signal-velocity conductivity
      on the #u jump, #u still holding `u^n` there, applied by
      #radiation_end_force_propagation. The force loop's dispatch fires
      both sides of a pair whenever either kernel reaches, which is what
      keeps the mirrored credit/debit pair whole at h_i != h_j. Scratch:
      zeroed once per step by radiation_snapshot_part_propagation, like
      #grad_u, since the force loop runs exactly once per step.
      PHYSICAL, like every accumulator the pairwise operators fill. */
  float dissipation_u;

  /*! `(1/rho) grad(rho u)` accumulator, gradient loop
      (radiation_propagation_iact.h). Scratch: zeroed once per step by
      radiation_snapshot_part_propagation, since the gradient loop runs
      exactly once per step (never re-run across h-iterations). PHYSICAL:
      per physical length, not per comoving one. */
  float grad_u[3];

  /*! Mass-specific per-band emission dose still owed to this particle by
      every star that has touched it (#radiation_iact_nonsym_feedback_apply),
      not yet injected into #u. Persistent, dumped with #part like
      #specific_flux; zero at first init, no IC field (an IC has no
      notion of "dose in flight"). Drained once per step, for active
      particles only, by #radiation_snapshot_part_propagation into
      #u_source_rate; every star's touch only ever
      adds to it, so any number of stars on any time bins superpose without
      losing or double-counting emission. */
  float u_dose_reservoir;

  /*! This step's mass-specific per-band source rate, drawn down from
      #u_dose_reservoir by
      #radiation_snapshot_part_propagation and consumed by
      #radiation_end_force_propagation's exact-relaxation update. Scratch:
      recomputed every step for active particles, not restart-critical (an
      inactive particle recomputes it correctly the moment it next becomes
      active), but dumped anyway since it lives in #part alongside the
      persistent fields above. */
  float u_source_rate;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Most negative #u written by #radiation_end_force_propagation since the
      previous snapshot, 0 if none was negative. Reset on the first update
      after a snapshot (#feedback_part_data.u_min_snapshot_index). "Since the
      previous snapshot" means since the last increment of
      #engine.snapshot_output_count, which also happens when a FOF seeding
      catalogue is dumped (FOF:dump_catalogue_when_seeding, engine.c), not
      only at a real snapshot dump. Written as
      "PEMinimumSpecificEnergies"/"LWMinimumSpecificEnergies". PHYSICAL, like
      #u. */
  float u_min_since_snapshot;

  /*! Cumulative mass-specific dose this particle has been handed by the
      dose reservoir since first init, RESCALED exactly as
      #radiation_end_force_propagation rescales it into #u
      (`c_hyp/c`) but NOT relaxed by #radiation_relaxation_phi_factor: the
      raw amount attempted this step, `dt_prev*(c_hyp/c)*u_source_rate`,
      summed step over step. Never reset (unlike #u_min_since_snapshot): an
      energy-conservation check reads this as a running total at every
      snapshot, so a mid-run reset would break its own conservation
      identity. Written as
      "PECumulativeInjectedSpecificEnergies"/
      "LWCumulativeInjectedSpecificEnergies". PHYSICAL, like #u. */
  float cumulative_injected;

  /*! Cumulative mass-specific energy this particle's #u update has
      attributed to relaxation (dust absorption plus the cosmological
      redshift term folded into the same decay) and to the fraction of
      this step's source/transport terms that never reached #u because
      the step was optically thick (`(1-phi)` of each), since first init.
      Exactly `(u_prev + dt_prev*phi*dissipation_u)*(1-e) +
      ((c_hyp/c)*u_source_rate - div_specific_flux)*dt_prev*(1-phi)`,
      `e = exp(-a)`, `phi = radiation_relaxation_phi_factor(a)`, `a =
      (c_hyp*kappa + H)*dt_prev` -- the same `e`/`phi`/`a` #u's own update
      uses this step, read before #u is overwritten. This is NOT the pure
      dust-extinction loss alone: #div_specific_flux (transport) and
      #dissipation_u (the artificial-dissipation source) are folded in
      too, because the exact-relaxation update mixes all three under one
      `phi`. Summed with #cumulative_injected and the current #u at a
      snapshot, `E + Abs - Inj` isolates exactly the part of the update
      the closed-form split above does not attribute to injection or the
      surviving field: the transport and dissipation residual, which the
      SPH divergence's kernel-sum identity and the dissipation's pairwise
      antisymmetry drive to ~0 when summed over the whole particle set.
      Never reset, for the same reason as #cumulative_injected. Written as
      "PECumulativeAbsorbedSpecificEnergies"/
      "LWCumulativeAbsorbedSpecificEnergies". PHYSICAL, like #u. */
  float cumulative_absorbed;
#endif
};

/**
 * @brief Per-operator ISRF band-physics coefficients carried by each hydro
 * particle, one instance per #radiation_isrf_operator in
 * #feedback_part_data.isrf_operator. A moment reads its own operator
 * through #radiation_isrf_moment_to_operator.
 */
struct feedback_isrf_operator_data {

  /*! Band-specific local linear dust absorption rate (see
      #radiation_get_part_linear_absorption_rate), cached once per step
      (radiation_snapshot_part_propagation) so the propagation loops do not
      recompute it, and the same unit conversion, per neighbour pair.
      The raw physical rate: distinct from and NOT interchangeable with the
      injection-side extinction's own, independently-computed kappa
      (#radiation_get_part_ISRF_extinction_factors). */
  float kappa;

  /*! M1 closure tensor `D(f)` from the owning moment's own
     #feedback_isrf_moment_data.u, #feedback_isrf_moment_data.specific_flux
     (#radiation_isrf_operator_owner) and #feedback_part_data.c_hyp, cached
     by #radiation_cache_m1_closure_part (drift-time reset and, once #c_hyp
     itself is known, the density ghost) so the gradient loop reads it per
     pair without rebuilding it. */
  float m1_closure_D[3][3];

  /*! Kernel-mean of the neighbours' |rho_prev*u_prev|, density loop
      (radiation_propagation_iact.h): the local field-scale reference the
      negativity trigger (#radiation_end_gradient_propagation)
      divides an undershoot by. Scratch: zeroed every h-iteration alongside
      #feedback_isrf_moment_data.div_specific_flux. */
  float ngb_mean_abs_u_V;

  /*! Negativity-triggered artificial-dissipation coefficient, REACTIVE
      component: raised by the
      negativity trigger and decayed otherwise, updated once per step in
      #radiation_end_gradient_propagation (not the density ghost, which
      re-runs across h-iterations). Persistent, dumped with #part like
      #feedback_isrf_moment_data.specific_flux; zero at first init, no IC
      field. Read by THIS step's force loop, which dissipates `u^n`, the
      same state the trigger read, so an undershoot is corrected one step
      after it appears. Applied to a pair UNGATED, as
      `max(trigger_i, trigger_j)`: the trigger only ever fires on a
      particle that is already locally wrong, so it is local by
      construction.

      `a`-SCALING: dimensionless, exponent 0. */
  float dissipation_alpha_trigger;

  /*! Negativity-triggered artificial-dissipation coefficient, ANTICIPATORY
      component:
      the `h/lambda`-gated floor (#radiation_dissipation_alpha_floor_band),
      which supplies dissipation on a positive front the negativity trigger
      is structurally blind to. Written alongside the trigger component
      above, in the same once-per-step ghost, and persistent for the same
      reason.

      Kept SEPARATE from the trigger rather than pre-combined with max(),
      even though the force loop combines them unconditionally: the two are
      produced by different mechanisms on different conditions (reactive
      undershoot response versus anticipatory resolution gating), so a run
      that dissipates too much or too little is only diagnosable when the
      two contributions can be read apart.

      `a`-SCALING: dimensionless, exponent 0. */
  float dissipation_alpha_floor;
};

/**
 * @brief Feedback fields carried by each hydro particles
 *
 * Carries the HII ionization tag core (radiation.c's
 * radiation_tag_part_as_ionized() and friends). A struct part field, not
 * struct xpart, so it rides the particle's normal MPI exchange and restart
 * dump automatically, unlike xpart's owner-only payload
 * (feedback_xpart_data.HII_region: excess_photon_energy_HI,
 * photoionization_rate_HI), which stays local to the owning rank.
 */
struct feedback_part_data {

  /*! Tag to mark the particle as ionized. */
  char is_ionized;

  /*! Largest #part.time_bin among this particle and every neighbour
      accumulated by the ISRF density loop this h-iteration
      (radiation_propagation_iact.h's runner_iact_isrf_propagation and
      runner_iact_nonsym_isrf_propagation), reset to #part.time_bin at
      the start of each h-iteration by radiation_init_part_propagation.
      Drives #c_hyp: the propagation speed is set from the SLOWEST
      particle in the kernel, not this particle's own step, so that the
      receiver's CFL condition holds for every pair by construction (see
      radiation_end_density_propagation). Placed here, right after
      #is_ionized, to land in that field's own compiler padding rather
      than growing #part (verified by a standalone sizeof/offsetof
      probe, not by inspection: `struct part` stays 640 bytes). */
  timebin_t max_ngb_time_bin;

  /*! Id of the star that ionized this particle. */
  long long star_id;

  /*! Simulation time until which this particle stays flagged as ionized. */
  double end_time;

  /*! Neutral hydrogen mass fraction, cached by the cooling step right after
      its species update (grackle_1+: HI_frac; grackle_0: 1.0f
      unconditionally, no species tracked). Not read by anything yet: a
      future MPI scheme's F3 latch must consume the PREVIOUS pass's value,
      not the current step's.

      The cache write is skipped on cooling_new_energy()'s early-return
      paths (subgrid-ionized floor, or pinned by
      IONIZATION_FEEDBACK_DEBUG_FIXED_*_TEMPERATURE_K), so on a particle's
      first such step this field is still its zero-init default, 0.0f, the
      OPPOSITE of the 1.0f "no data" sentinel above. A reader must not treat
      0.0f as "fully ionized" without checking the cache has actually been
      written for that particle. */
  float neutral_H_frac;

  /*! Per-moment ISRF transport state, indexed by #radiation_isrf_moment. */
  struct feedback_isrf_moment_data isrf_moment[ISRF_MOMENT_COUNT];

  /*! Per-operator ISRF band-physics coefficients, indexed by
      #radiation_isrf_operator; a moment reaches its own operator through
      #radiation_isrf_moment_to_operator. */
  struct feedback_isrf_operator_data isrf_operator[ISRF_OPERATOR_COUNT];

  /*! Comoving density snapshot, cached once per step by
      radiation_snapshot_part_propagation at the same call site as
      #feedback_isrf_moment_data.u_prev (before this step's density accumulators
      are reset), so it holds the previous step's fully-converged comoving
      density. Needed because the density loop's kernel-mean accumulation
      (radiation_propagation_iact.h) runs interleaved with SPH's own density
      sum: `p->rho` is a partial accumulator there, not a density, until the
      density ghost finalizes it. The gradient loop's `grad(u)` and the force
      loop's `div(F)` use this SAME snapshot rather than the by-then-available,
      more current ghost-finalized density, because the staggered time
      integrator's stability on a disordered particle distribution depends
      on `grad` being minus the adjoint of `div` in the `m*rho` inner
      product: that identity only holds when both operators are built from
      the same `rho_i`/`rho_j`. Never legitimately 0 or negative: seeded to
      1.0f at first init (before any real density has ever been computed,
      see radiation_isrf.c) rather than 0.0f, since `grad(u)`'s own formula
      needs `rho_i` itself (not just `1/rho_i`), and a 0.0f seed would turn
      into +inf under any reciprocal taken from it. A placeholder value
      is safe there regardless, since `u`/`F` are also still 0 at that
      point, so every term the placeholder feeds into is itself 0. */
  float rho_prev;

  /*! This particle's kernel-local hyperbolic propagation speed
      (`c_hyp_i = min(C_hyp*h_i/dt_max(i), c)`, `dt_max(i)` the longest
      timestep among this particle and every neighbour in its kernel,
      #max_ngb_time_bin), cached for active particles by
      radiation_end_density_propagation (the density ghost, after the
      h-iteration converges), so every receiver's CFL condition
      `c_i*dt_j <= C_hyp*h_i` holds by construction for the pairs the
      force loop reaches. An inactive particle's value is simply last
      active step's, like #time_bin itself. Shared by both bands (unlike
      #feedback_isrf_operator_data.kappa): the propagation speed is a
      property of the particle's resolution and its kernel's slowest
      clock, not of its dust opacity. */
  float c_hyp;

  /*! This particle's own physical timestep (`dt_i`, NOT #max_ngb_time_bin's
      `dt_max(i)`: only #c_hyp uses the kernel maximum), cached by
      radiation_snapshot_part_propagation (the drift, once per step, every
      particle whether active or not) so the exact-relaxation finalizes
      (radiation_isrf.c) and the dose-reservoir drawdown/restore
      (radiation_snapshot_part_propagation/radiation_part_has_no_neighbours)
      do not need to recompute it a second and third time. */
  float dt_prev;

#ifdef SWIFT_DEBUG_CHECKS
  /*! #engine.snapshot_output_count at the last write of
      #feedback_isrf_moment_data.u_min_since_snapshot: the index of the snapshot
      those values belong to. Incremented by a FOF seeding catalogue dump as
      well as a real snapshot; see
     #feedback_isrf_moment_data.u_min_since_snapshot. */
  int u_min_snapshot_index;
#endif

  /*! With ISRF_propagation off: simulation step (#engine.ti_current)
      #feedback_isrf_moment_data.u were last written at.
      radiation_iact_nonsym_feedback_apply compares this against the current
      step: a match means some star already wrote this step, so a further touch
      (a second illuminating star) sums into the existing value; a mismatch
      means this is the first touch this step, so #feedback_isrf_moment_data.u
     are zeroed before summing. This is what makes the field an instantaneous
      strength rather than an ever-growing total, while still summing multiple
      simultaneously-illuminating stars correctly within one step. With
      ISRF_propagation on, this is only bookkeeping (the last step any star
      touched this particle): the dose-reservoir form never resets
      #feedback_isrf_moment_data.u, so no consumer relies on it there.
      feedback_first_init_part sets this to -1 (never a valid step) so the very
      first touch of a particle's life also resets rather than summing onto
      uninitialized memory. */
  integertime_t ISRF_last_touch_ti;

  /*! Has this particle been illuminated (any band's
     #feedback_isrf_moment_data.u nonzero) by any star's injection pass, and is
     that illumination episode still live? Dedicated flag, not inferred from
     #feedback_isrf_moment_data.u itself, since those now reset every step a
     star touches this particle and so cannot signal "newly illuminated" via a
     zero-crossing. Mirrors #is_ionized's claimed/not-claimed cycle, including
     the reset half: gates a first-touch-only timestep_sync_part call in
      radiation_iact_nonsym_feedback_apply, mirroring
      feedback_hii_claim_part/feedback_iact_HII_maintain_ionized_part's own
      claim-vs-maintain split, and is cleared once #ISRF_illumination_end_ti
      lapses (radiation_gas.c:radiation_reset_part_ISRF_illumination_tag,
      called from feedback_reset_part), exactly as cooling clears #is_ionized
      once its own end_time lapses (cooling_gear_subgrid.h). A particle that
      leaves every illuminating star's kernel, then re-enters one later (a
      star re-approaches, a new star's kernel reaches it, or its own h
      changes), therefore gets a fresh sync on re-illumination instead of
      being silently skipped forever. */
  char is_illuminated_ISRF;

  /*! Absolute integer time (#engine.ti_current units) until which
      #is_illuminated_ISRF stays set. Renewed to
      `ti_current + RADIATION_ISRF_TAG_LIFETIME_INTERVALS * ti_step` on
      every injection touch (ti_step = the illuminating star's own
      integer timestep), whether or not this is the particle's first touch
      this episode. Mirrors #feedback_iact_HII_maintain_ionized_part's
      per-pass renewal of the HII tag's own end_time. The
      RADIATION_ISRF_TAG_LIFETIME_INTERVALS buffer (radiation.h) keeps
      the window from lapsing between two touches by a star on a coarser
      time bin than this particle's own once-per-step expiry check
      (feedback_reset_part). feedback_first_init_part sets this to -1 so a
      never-illuminated particle's garbage/zero-initialized state can never
      read as "still illuminated". */
  integertime_t ISRF_illumination_end_ti;

  /*! Absolute integer time (#engine.ti_current units) by which every dose
      currently held in #feedback_isrf_moment_data.u_dose_reservoir must have
      been fully drained. Extended to `ti_current + ti_step_star` on every
      star touch (never reset), so it always covers the latest-finishing
      contributing star's own step. feedback_first_init_part sets this to -1,
      like #ISRF_illumination_end_ti, so a never-touched particle's
      reservoir is never mistaken for one with a live horizon. */
  integertime_t ISRF_reservoir_end_ti;
};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! specific energy received from supernovae */
  float delta_u;

  /*! Momemtum received from a supernovae */
  float delta_p[3];

  /*! Radiation struct */
  struct {

    /*! Momemtum received from a radiation_pressure */
    float delta_p[3];

    /*! Lifetime-cumulative |delta_p| from radiation pressure (scalar sum,
        not vector: isotropic kicks would else cancel). Snapshot field
        CumulativeMomentumFromRadiationPressure (feedback_io.h); moved here
        from tracers_xpart_data so it no longer needs --with-tracers=GEAR. */
    float cumulative_momentum;

    /*! Largest single-event radiation-pressure kick velocity (outflow
        diagnostic). Snapshot field MaxKickVelocityFromRadiationPressure;
        same relocation as #cumulative_momentum above. */
    float max_kick_velocity;
  } radiation;

  /*! HII ionization owner-computed payload, local to the owning rank (the
      tag core itself lives on struct part's feedback_data, see
      feedback_part_data's own doxygen). */
  struct {

    /*! Mean photon energy above the 13.6 eV HI ionization threshold for
        the tagging star, frozen at tag time (only set when
        GEARFeedback:HII_couple_ionization_rate is on; 0 otherwise). Stored
        in cgs (erg), not internal units, since the absolute per-particle
        value underflows float precision in this project's internal unit
        system. */
    float excess_photon_energy_HI;

    /*! Photoionization rate coefficient Gamma_HI from the tagging star at
        this particle's location, frozen at tag time (internal 1/time;
        only set when GEARFeedback:HII_couple_ionization_rate is on, 0
        otherwise). */
    float photoionization_rate_HI;

  } HII_region;

  /*! Indicator if the particule receive energy from SN specifically */
  char hit_by_SN;

  /*! Indicator if the particle receives energy from SW specifically */
  char hit_by_winds;

  /*! Indicator if the particle receives energy from radiation pressure
   * specifically */
  char hit_by_radiation;
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Is the star dead? */
  int is_dead;

  /*! Gas density at the star location. This is also the inverse of
   * normalisation factor used for the enrichment. */
  float enrichment_weight;

  /*! Does the particle needs to go through the feedback loops? */
  char will_do_feedback;

  /*! Does the particle needs to go through the HII ionization loop? */
  char will_do_HII_ionization;

  /*! Integer number of neighbours */
  int num_ngbs;

  /*! Gas density gradient at the star location */
  float grad_rho_star[3];

  /*! Gas metallicity at the star location */
  float Z_star;

  /*! Number of Ia supernovae */
  float number_snia;

  /*! Number of II supernovae */
  float number_snii;

  /* Supernovae data struct */
  struct {

    /*! Energy injected in the surrounding particles */
    float energy_ejected;

    /*! Total mass ejected by the supernovae */
    float mass_ejected;

  } supernovae;

  /*! Chemical composition of the mass ejected */
  double metal_mass_ejected[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Stellar winds data struct */
  struct {

    /*! Energy injected in the surrounding particles */
    float energy_ejected;

    /*! Mass injected in the surrounding particles */
    float mass_ejected;

  } winds;

  /*! Radiation data structs */
  struct {

    /*! Bolometric luminosity (physical units) from the stellar evolution */
    double L_bol;

    /*! Number of ionizing photons per unit time, split evenly across the
        n_HII_pixels active angular pixels (this is a HUGE number, so must
        be a double) (physical units) from the stellar evolution. A pure
        rate: never debited, so it stays valid all pass for the rate-coupled
        flux calculation (GEARFeedback:HII_couple_ionization_rate). */
    double dot_N_ion_pix[HII_MAX_ANGULAR_PIXELS];

    /*! dot_N_ion_pix as it stood at the previous HII rebuild pass, cached so
        radiation_open_ionizing_photon_budget() can integrate the emission
        rate over dt_back with a trapezoid rule (average of the rate at the
        two ends of the interval) instead of a rectangle rule at the rate
        "now" alone. A monotonically-declining SSP emission rate makes the
        rectangle rule systematically under-issue photons (biased, not
        noise -- it does not average out over passes). Negative is the
        sentinel for "no previous pass yet" (set at star formation), which
        falls back to the pre-fix rectangle rule for a star's first pass. */
    double dot_N_ion_pix_prev[HII_MAX_ANGULAR_PIXELS];

    /*! Photon *count* spendable per pixel this HII rebuild pass:
        dot_N_ion_pix * dt_back (elapsed time since the last pass), debited
        by radiation_consume_ionizing_photons. Budgeting counts over the real
        elapsed interval, instead of comparing bare rates, is what makes the
        ionized extent independent of the rebuild cadence. */
    double N_ion_budget_pix[HII_MAX_ANGULAR_PIXELS];

    /*! Number of active angular pixels this star is currently using
        (1 = spherical/HEALPix disabled) */
    int n_HII_pixels;

    /*! Mass in the HII region generated by this star particle */
    float mass_HII_region;

    /*! Co-moving HII region radius before the star died or was not
        eligible to form HII regions anymore. Same algorithm's bookkeeping
        caveat as the live #h_hii it is retired from. */
    float final_HII_radius;

    /*! Ionized gas mass of that same final HII region. */
    float final_HII_mass;

    /*! Star age at this HII region's last rebuild pass. Double, since
        dt_back (age now minus this) scales the whole photon budget above. */
    double HII_region_last_rebuild;

    /*! Star age the last time this star was given a chance to rebuild its
        HII region, whether or not gas was found (unlike
        HII_region_last_rebuild, which only advances on an actual rebuild).
        Anchors the photon-budget interval dt_back instead, so a pass
        skipped by a gas-free working-level cell does not make the next
        real pass look like it covers the whole gap. See
        runner_radiation_feedback.c. */
    double HII_region_last_attempt;

    /*! Mean photon energy above the 13.6 eV HI ionization threshold,
        cached once per HII rebuild pass (only computed when
        GEARFeedback:HII_couple_ionization_rate is on; 0 otherwise). Stored
        in cgs (erg), not internal units, since the absolute per-particle
        value underflows float precision in this project's internal unit
        system. */
    float mean_excess_photon_energy_HI;

    /*! Moment luminosity (physical units), indexed by #radiation_isrf_moment:
        non-ionizing PE, 6-11.2 eV, and Lyman-Werner, 11.2-13.6 eV (H2
        photodissociating photons). Read from the radiation table's own
        L_PE/L_LW (or Integrated_L_PE/Integrated_L_LW) datasets, which
        carry the band split directly. Feeds the injection term; only
        computed when GEARFeedback:with_interstellar_radiation_field is on, 0
        otherwise. #ISRF_MOMENT_LW_PHOTON is set equal to #ISRF_MOMENT_LW: it
        is energy-equivalent in erg/s like every other entry, not a photon
        rate, so no new unit handling applies to it. */
    double L_band[ISRF_MOMENT_COUNT];

    /*! Photospheric effective temperature (internal units), a
        stellar-evolution diagnostic written to the snapshot and not used
        by any feedback channel. For an IMF-population particle it is the
        value at the upper mass bound of the stars still alive, i.e. the
        hottest surviving star, not an IMF average. 0 when the radiation
        table carries no "Teff" dataset. */
    float teff;

    /*! This star's feedback timestep (proper time, internal units),
        cached once per step by feedback_prepare_radiation_feedback so
        radiation_iact_nonsym_feedback_apply does not repeat the
        cosmology table lookup for every gas neighbour. */
    float Delta_t;

#ifdef SWIFT_DEBUG_CHECKS
    /*! Integer time at the start of the step #Delta_t was cached for.
        radiation_iact_nonsym_feedback_apply recomputes it and asserts a
        match before trusting the cached value. */
    integertime_t Delta_t_cached_ti_begin;
#endif

  } radiation;
};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
