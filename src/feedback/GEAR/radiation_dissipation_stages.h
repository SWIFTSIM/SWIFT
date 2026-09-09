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
#ifndef SWIFT_RADIATION_DISSIPATION_STAGES_GEAR_H
#define SWIFT_RADIATION_DISSIPATION_STAGES_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_dissipation_stages.h
 * @brief Compile-time switches and constants for the LW/FUV
 * artificial-dissipation stages (design-lw-fuv-design-b-dissipation.md
 * Section 5.2). This is the only place a stage is turned on.
 *
 * Stage 1 (the triggered scalar conductivity) is unconditional and is
 * configured by the GEARFeedback:LW_FUV_dissipation_* runtime parameters.
 * Stages 2-4 are compile-time options, off by default: uncomment one or
 * more of the three #defines below to build with them.
 *
 * Each stage's formulas are compiled in both states, following the MAGMA2
 * precedent (Section 5.1 item 2): the macro selects a `const int` that a
 * runtime `if` tests, so a disabled stage still goes through the compiler
 * and cannot rot. The macros additionally guard the per-particle fields the
 * stages need, so a build without them pays no memory for them; that is why
 * this is a leaf header with no includes of its own, rather than a block
 * inside radiation.h. struct feedback_part_data
 * (GEAR_thermal/feedback_struct.h) is reached through part.h, which
 * radiation.h itself includes, so it cannot include radiation.h back.
 */

/*! Stage 2: slope-limited linear reconstruction of the Stage-1 jump to the
    pair midpoint, in MAGMA2's compiled van Leer form (Rosswog 2026
    Eq. 149-150, Rosswog 2020b Eq. 21-23, Chan et al. 2021 Eq. 31). Keeps the
    dissipation small on a resolved gradient, so it earns its place together
    with a trigger that fires on fronts (Stage 4) rather than on its own.
    Adds 24 bytes per gas particle. */
/* #define RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION */

/*! Stage 3: anisotropic artificial dissipation of the flux moment along
    `n n` with the `d(div F)/dt` switch (Chan et al. 2021 Eq. 32-37, coded in
    src/rt/SPHM1RT). Covers the flux-driven regime the Stage-1 scalar term
    provably cannot reach, where the `c_hyp^2/rho^2` amplification at a
    strongly rarefied particle drives the growth. Adds 40 bytes per gas
    particle. */
/* #define RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX */

/*! Stage 4: anticipatory noise trigger blended into the Stage-1 coefficient,
    `alpha_aim = max(alpha_negativity, alpha_noise)` (Rosswog 2015a Eq. 83,
    Eq. 87-89 transcribed from `div v` to `div F`). Stage 1's own trigger is
    reactive: it can only raise the coefficient once the field has already
    gone negative, which leaves a window in which the cooling module reads a
    negative value. This one fires on the estimator inconsistency between a
    particle and its neighbours, before the sign change happens. Adds 16
    bytes per gas particle. */
/* #define RADIATION_LW_FUV_DISSIPATION_ANTICIPATORY_TRIGGER */

/*! Stage 2: separation, in units of the pair's larger smoothing length,
    above which the van Leer limiter is no longer suppressed by its Gaussian
    factor. MAGMA2's own `const_viscosity_eta_crit` value; the Gaussian's
    width, `exp(-25*(eta-eta_crit)^2)`, is hard-coded to match it. */
#define RADIATION_LW_FUV_DISSIPATION_ETA_CRIT 1.0f

/*! Stage 3: amplitude `A` of the `d(div F)/dt` flux switch, Chan et al. 2021
    Eq. 36, as shipped in src/rt/SPHM1RT/rt.h. */
#define RADIATION_LW_FUV_DISSIPATION_FLUX_SWITCH_AMPLITUDE 200.f

/*! Stage 4: reference value the noise indicator `N` is compared against in
    `alpha_noise = alpha_max*N/(N + this)`. PROVISIONAL AND UNGAUGED: the
    design (Section 5.2) takes it from where a glass steady state's own
    `div F` sign noise sits, which the Tier-1 five-corner runs measure once
    the accumulators exist. The value here is a conservative placeholder
    only: with `N` in [0, 1] it caps `alpha_noise` at `alpha_max/2` and gives
    a few percent of `alpha_max` at the few-percent sign noise a glass fixed
    point is expected to carry. Do not enable Stage 4 in a production run
    before this constant has been measured. */
#define RADIATION_LW_FUV_DISSIPATION_NOISE_REFERENCE 1.0f

#endif /* SWIFT_RADIATION_DISSIPATION_STAGES_GEAR_H */
