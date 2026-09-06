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
 * @brief Receiver-side LW/FUV dust extinction and Yukawa propagation
 * physics for GEAR: gas-side opacity, extinction, local absorption rate,
 * and the propagation mixing fraction.
 */

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
float radiation_get_comoving_gas_column_density_at_part(const struct part *p);
void radiation_get_part_LW_FUV_extinction_factors(
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct part *p, float Z, const struct cooling_function_data *cooling,
    float *extinction_FUV, float *extinction_LW);
float radiation_get_part_linear_absorption_rate(const struct unit_system *us,
                                                float Z, float rho_p,
                                                float sigma_d_band_cgs,
                                                float local_dust_to_gas_ratio);
float radiation_get_isrf_propagation_alpha(float h, float kappa_i, float w_min);
float radiation_compute_yukawa_w_min(const struct hydro_props *hydro_props);
float radiation_compute_yukawa_kernel_second_moment(
    const struct hydro_props *hydro_props);

#endif /* SWIFT_RADIATION_ISRF_GEAR_H */
