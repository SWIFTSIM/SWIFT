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
#ifndef SWIFT_TRACERS_NONE_H
#define SWIFT_TRACERS_NONE_H

/* Config parameters. */
#include <config.h>

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
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_part)
 */
static INLINE void tracers_after_timestep_part(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time,
    const double time_step_length, const int *const tracers_triggers_started) {}

/**
 * @brief Update the star particle tracers just after its time-step has been
 * computed.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_spart)
 */
static INLINE void tracers_after_timestep_spart(
    struct spart *sp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const double time_step_length,
    const int *const tracers_triggers_started) {}

/**
 * @brief Update the black hole particle tracers just after its time-step has
 * been computed.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_bpart)
 */
static INLINE void tracers_after_timestep_bpart(
    struct bpart *bp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const double time_step_length,
    const int *const tracers_triggers_started) {}

/**
 * @brief Initialise the tracer data at the start of a calculation.
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

/**
 * @brief Initialise the star tracer data at the start of a calculation.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 */
static INLINE void tracers_first_init_spart(struct spart *sp,
                                            const struct unit_system *us,
                                            const struct phys_const *phys_const,
                                            const struct cosmology *cosmo) {}

/**
 * @brief Initialise the black hole tracer data at the start of a calculation.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 */
static INLINE void tracers_first_init_bpart(struct bpart *bp,
                                            const struct unit_system *us,
                                            const struct phys_const *phys_const,
                                            const struct cosmology *cosmo) {}

/**
 * @brief Update the particles' tracer data after a stellar feedback
 * event.
 *
 * Nothing to do here.
 *
 * @param xp The extended particle data.
 */
static INLINE void tracers_after_feedback(struct xpart *xp) {}

/**
 * @brief Update the particles' tracer data with values before an AGN feedback
 * event. Note: this function is called in `black_holes_iact.h` before the
 * particle data are updated.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the basic particle data.
 * @param xp The extended particle data.
 * (internal physical units)
 */
static INLINE void tracers_before_black_holes_feedback(
    const struct part *p, struct xpart *xp, const float scale_factor) {}

/**
 * @brief Update the particles' tracer data after an AGN feedback
 * event.
 *
 * Nothing to do here.
 *
 * @param p The particle data.
 * @param xp The extended particle data.
 * @param with_cosmology Are we running with cosmology?
 * @param scale_factor The current scale-factor (if running with cosmo)
 * @param time The current time (if running without cosmo)
 * @param Amount of energy injected in the feedback event (internal physical
 * units)
 */
static INLINE void tracers_after_black_holes_feedback(
    const struct part *p, struct xpart *xp, const int with_cosmology,
    const float scale_factor, const double time, const double delta_energy) {}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * Nothing to do here.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
static INLINE void tracers_after_snapshot_part(const struct part *p,
                                               struct xpart *xp) {}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * Nothing to do here.
 *
 * @param sp the #spart.
 */
static INLINE void tracers_after_snapshot_spart(struct spart *sp) {}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * Nothing to do here.
 *
 * @param sp the #spart.
 */
static INLINE void tracers_after_snapshot_bpart(struct bpart *bp) {}

/**
 * @brief Split the tracer content of a particle into n pieces
 *
 * Nothing to do here.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void tracers_split_part(
    struct part *p, struct xpart *xp, const double n) {}
#endif /* SWIFT_TRACERS_NONE_H */
