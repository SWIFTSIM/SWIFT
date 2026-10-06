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
 * @brief Receiver-side LW/PE dust extinction and hyperbolic-relaxation
 * propagation physics for GEAR: opacity, extinction and absorption.
 */

#include "feedback_struct.h"

struct part;
struct xpart;
struct cosmology;
struct unit_system;
struct hydro_props;
struct engine;
struct cooling_function_data;
struct feedback_props;
struct phys_const;
struct radiation;
struct stellar_model;

/*! Photon-number-weighted mean Lyman-Werner photon energy of the radiation
 * table in erg, or 0 while radiation is inactive. Diagnostic only: the H2
 * photodissociation rate uses #RADIATION_SIGMA_H2_OVER_E_LW_CGS. Set once at
 * start-up and on restart by #radiation_set_lw_photon_energy_cgs. */
extern double radiation_lw_photon_energy_cgs;

void radiation_set_lw_photon_energy_cgs(const struct radiation *rad,
                                        const struct stellar_model *sm);

/**
 * @brief Set the band-edge weights of @p fb_props from the radiation table.
 *
 * @param fb_props (output) The #feedback_props to set.
 * @param rad The main stellar model's #radiation.
 * @param sm The main #stellar_model, for its IMF mass range.
 */
void radiation_set_band_edge_coefficients(struct feedback_props *fb_props,
                                          const struct radiation *rad,
                                          const struct stellar_model *sm);

void radiation_first_init_part(struct part *restrict p);
void radiation_snapshot_part_propagation(struct part *p,
                                         const struct engine *e);
void radiation_init_part_propagation(struct part *p);
void radiation_end_density_propagation(struct part *p, const struct engine *e);
void radiation_part_has_no_neighbours(struct part *p, const struct engine *e);
void radiation_end_gradient_propagation(struct part *p, const struct engine *e);
void radiation_end_force_propagation(struct part *p, const struct engine *e);
float radiation_isrf_part_timestep(const struct part *restrict p,
                                   const struct engine *e);
float radiation_get_comoving_extinction_path(
    const struct feedback_props *fb_props, const struct part *p,
    const struct xpart *xp, const float r, const struct cosmology *cosmo,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling);
float radiation_get_comoving_gas_column_density_at_part(
    const struct part *p, const float extinction_path);
void radiation_get_part_ISRF_extinction_factors(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct part *p, float Z,
    const struct cooling_function_data *cooling, const float extinction_path,
    float extinction[ISRF_OPERATOR_COUNT]);
float radiation_get_part_linear_absorption_rate(
    const struct unit_system *us, const struct phys_const *phys_const, float Z,
    float rho_p, float sigma_d_band_cgs, float local_dust_to_gas_ratio);
float radiation_relaxation_phi_factor(float a);

#endif /* SWIFT_RADIATION_ISRF_GEAR_H */
