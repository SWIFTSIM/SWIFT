/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TRACERS_NONE_H
#define SWIFT_TRACERS_NONE_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cooling.h"
#include "part.h"
#include "tracers_struct.h"

/**
 * @brief Update the particle tracers just after it has been initialised at the
 * start of a step.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 */
static INLINE void tracers_after_init(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after it has been drifted.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 */
static INLINE void tracers_after_drift(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after its time-step has been
 * computed.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 */
static INLINE void tracers_after_timestep(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after its time-step has been
 * computed.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 */
static INLINE void tracers_first_init_xpart(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling) {}

#endif /* SWIFT_TRACERS_NONE_H */
