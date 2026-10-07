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
#ifndef SWIFT_FEEDBACK_GEAR_RADIATION_H
#define SWIFT_FEEDBACK_GEAR_RADIATION_H

#include "../../feedback_properties.h"
#include "cooling.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

/**
 * @file src/feedback/GEAR/feedback_radiation.h
 * @brief Gas-side functions of the GEAR subgrid radiation (HII tag and ISRF
 * fields of a #part). Only the GEAR thermal feedback module carries this
 * state.
 */

struct engine;

char feedback_part_can_be_ionized(const struct part *p, const struct xpart *xp,
                                  const struct engine *e);
void feedback_iact_HII_ionization(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling,
    const struct feedback_props *feedback_props, const integertime_t ti_begin,
    const double time, const double dt_back);

double feedback_iact_HII_maintain_ionized_part(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const double time,
    const double dt_back);

char feedback_is_part_tagged_as_ionized(const struct part *p,
                                        const struct xpart *xp);
long long feedback_get_part_ionized_star_id(const struct part *p,
                                            const struct xpart *xp);

double feedback_get_part_u_PE(const struct part *p);
double feedback_get_part_u_LW(const struct part *p);
double feedback_get_part_u_LW_PHOTON(const struct part *p);
float feedback_get_part_dissipation_alpha_PE(const struct part *p);
float feedback_get_part_dissipation_alpha_LW(const struct part *p);
float feedback_get_part_dissipation_alpha_LW_PHOTON(const struct part *p);
float feedback_get_part_div_specific_flux_PE(const struct part *p);
float feedback_get_part_div_specific_flux_LW(const struct part *p);
float feedback_get_part_div_specific_flux_LW_PHOTON(const struct part *p);
void feedback_get_part_specific_flux_PE(const struct part *p, float *ret);
void feedback_get_part_specific_flux_LW(const struct part *p, float *ret);
void feedback_get_part_specific_flux_LW_PHOTON(const struct part *p,
                                               float *ret);
float feedback_get_part_u_min_since_snapshot_PE(const struct part *p,
                                                const struct engine *e);
float feedback_get_part_u_min_since_snapshot_LW(const struct part *p,
                                                const struct engine *e);
float feedback_get_part_u_min_since_snapshot_LW_PHOTON(const struct part *p,
                                                       const struct engine *e);
float feedback_get_part_cumulative_injected_PE(const struct part *p);
float feedback_get_part_cumulative_injected_LW(const struct part *p);
float feedback_get_part_cumulative_injected_LW_PHOTON(const struct part *p);
float feedback_get_part_cumulative_absorbed_PE(const struct part *p);
float feedback_get_part_cumulative_absorbed_LW(const struct part *p);
float feedback_get_part_cumulative_absorbed_LW_PHOTON(const struct part *p);
float feedback_get_part_c_hyp(const struct part *p);
float feedback_get_part_pending_specific_energy(const struct part *p, int m);

#endif /* SWIFT_FEEDBACK_GEAR_RADIATION_H */
