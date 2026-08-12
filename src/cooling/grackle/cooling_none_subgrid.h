/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_NONE_COOLING_SUBGRID_H
#define SWIFT_NONE_COOLING_SUBGRID_H

/**
 * @file src/cooling/grackle/cooling_none_subgrid.h
 * @brief Subgrid model for none cooling, independent from grackle.
 */

/**
 * @brief Update the gas properties with subgrid physics.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current simulation time.
 * @return Always 0 -- this subgrid model never floors a particle's energy.
 */
INLINE static int cooling_update_part_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double dt, double dt_therm, double time) {
  return 0;
}

/**
 * @brief Compute Grackle's RT_heating_rate/RT_HI_ionization_rate for a
 * particle tagged by GEAR's rate-coupled HII feedback. No-op here: this
 * subgrid model never tags particles as rate-coupled.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param time_units Internal-time-to-cgs-time conversion factor.
 * @param heating_rate_cgs (return) Heating rate, in raw cgs. Unused.
 * @param HI_ionization_rate (return) Photoionization rate coefficient, in
 *        internal 1/time. Unused.
 * @return Always 0 -- no per-particle rates from this subgrid model.
 */
INLINE static int cooling_get_rate_coupled_RT_fields_subgrid(
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp, double time_units, double *heating_rate_cgs,
    double *HI_ionization_rate) {
  return 0;
}

/**
 * @brief Debug-only: hold every non-ionized particle fixed at a fixed
 * temperature. No-op here: this subgrid model never forces a particle's
 * energy.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param u_out (return) The forced internal energy. Unused.
 * @return Always 0 -- this subgrid model never forces a particle's energy.
 */
INLINE static int cooling_debug_fix_neutral_temperature_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, float *u_out) {
  return 0;
}

/**
 * @brief Expire a GEAR rate-coupled HII tag once Grackle has consumed it
 * for this step's solve. No-op here: this subgrid model never tags
 * particles as rate-coupled.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param time The current simulation time.
 */
INLINE static void cooling_expire_rate_coupled_tag_subgrid(
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double time) {}

/**
 * @brief Cache this step's neutral hydrogen mass fraction on the
 * particle's feedback-model core. No-op here: this subgrid model has no
 * tag core to cache into.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
INLINE static void cooling_cache_neutral_H_fraction_subgrid(
    const struct cooling_function_data *cooling, struct part *p,
    const struct xpart *xp) {}

#endif /* SWIFT_NONE_COOLING_SUBGRID_H */
