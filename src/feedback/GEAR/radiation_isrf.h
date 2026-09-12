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
#ifndef SWIFT_RADIATION_ISRF_GEAR_H
#define SWIFT_RADIATION_ISRF_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_isrf.h
 * @brief Receiver-side LW/FUV dust extinction and hyperbolic-relaxation
 * propagation physics for GEAR: gas-side opacity, extinction, local
 * absorption rate, and the propagation mixing fraction.
 */

/*! Absolute guard on the denominator of the pair-jump contrast ratio
    `q_ij` (#radiation_dissipation_force_accumulate_band), in the internal
    units of `rho*u`. NOT a free physics parameter: the relative term
    `1e-3*u_V_scale` dominates it in every case where the pair carries any
    field at all, and this constant only decides the answer when both
    particles' `rho*u` AND both their kernel-mean references are zero -- in
    which case the numerator is zero too, so any positive value gives
    `q_ij = 0` (floor fully on), the behaviour of the scheme before the gate
    existed. That state is the whole box's initial condition, not a corner
    case, so the value satisfies a two-sided bound rather than being taken
    as the smallest representable float.

    Upper: it must stay far below the faintest genuine per-particle
    `|rho*u|`, or it suppresses `q_ij` and pins the floor fully on in the
    diffuse low-field phase the gate exists to treat correctly. The faintest
    values measured across the shipped examples are `~3e-07`
    (ISRFHyperbolicPropagation, FUV) and `~7e-14` (HomogeneousBox, LW), so
    `1e-20` is close to seven orders below the faintest real signal either
    band produces. It carries no `a^dim` weight of its own, unlike the
    comoving `rho*u` it is added to, so that margin (~7e6, i.e. ~6.85
    orders of magnitude, from the binding LW figure above) is spent by
    `a^dim` itself around `a ~ (1/7e6)^(1/3) ~ 0.005` and is gone below
    it. This is an honestly-disclosed, unresolved comoving-scaling
    assumption, not a proven bound. What happens there is the
    degradation named above, not an error: the guard starts to dominate,
    `q_ij` goes to zero and the floor applies ungated, as it did before the
    gate existed.

    Lower: the denominator it guards is divided into by `q0`
    (#LW_FUV_dissipation_pair_gate_q0, calibrated in [0.2, 0.5]), and
    `-ffast-math` reassociates that into `|d_ij| / (denominator * q0)`. With
    flush-to-zero also active, a denormal `denominator * q0` becomes exactly
    zero and the all-zero-field case evaluates `0/0`, filling every
    particle's `u_FUV`/`u_LW` with NaN from the first step. `FLT_MIN` has no
    headroom at all (`FLT_MIN * 0.3` is denormal); `1e-20` keeps the product
    normal for any `q0` at or above `1e-4`. The `min(q_ij_raw, 1.f)` clamp
    at the point of use is NOT a backstop for this: under `-ffast-math` it
    is folded away in the optimized binary (confirmed by disassembly), so
    this constant's own magnitude is what prevents the failure, not the
    clamp. */
#define RADIATION_LW_FUV_DISSIPATION_U_V_ABSOLUTE_FLOOR 1e-20f

/*! Module-scope mirror of #feedback_props.LW_FUV_dissipation_pair_gate_q0,
    written once by feedback_props_init()/feedback_struct_restore() and
    read-only afterwards.

    The pair-gate knee is needed inside
    #radiation_dissipation_force_accumulate_band, which is reached from
    runner_iact_[nonsym_]isrf_dissipation. Those two hooks are called from
    ~26 sites in src/runner_doiact_functions_hydro.h, a template shared by
    every SPH module, and their signatures carry no #feedback_props. Passing
    the parameter down that template would push a GEAR-only property through
    module-neutral core code; storing it per particle would spend memory on
    a value identical for every particle. A single write-once module-owned
    scalar does neither.

    Zero means "not yet written" and is treated as gate-disabled (`G_ij =
    1`, floor fully on) at the point of use, so a code path that somehow
    reaches the force loop before either writer ran degrades to the ungated
    floor rather than silently switching the floor off. */
extern float radiation_lw_fuv_dissipation_pair_gate_q0;

struct part;
struct xpart;
struct cosmology;
struct unit_system;
struct hydro_props;
struct engine;
struct cooling_function_data;

void radiation_first_init_part(struct part *restrict p);
void radiation_snapshot_part_propagation(struct part *p,
                                         const struct engine *e);
void radiation_init_part_propagation(struct part *p);
void radiation_end_density_propagation(struct part *p, const struct engine *e);
void radiation_part_has_no_neighbours(struct part *p, const struct engine *e);
void radiation_end_gradient_propagation(struct part *p, const struct engine *e);
void radiation_end_force_propagation(struct part *p, const struct engine *e);
float radiation_get_comoving_gas_column_density_at_part(const struct part *p);
void radiation_get_part_LW_FUV_extinction_factors(
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct part *p, float Z, const struct cooling_function_data *cooling,
    float *extinction_FUV, float *extinction_LW);
float radiation_get_part_linear_absorption_rate(const struct unit_system *us,
                                                float Z, float rho_p,
                                                float sigma_d_band_cgs,
                                                float local_dust_to_gas_ratio);
float radiation_relaxation_phi_factor(float a);

#endif /* SWIFT_RADIATION_ISRF_GEAR_H */
