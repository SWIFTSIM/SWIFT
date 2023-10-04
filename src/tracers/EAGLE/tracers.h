/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2021 Edo Altamura (edoardo.altamura@manchester.ac.uk)
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
#ifndef SWIFT_TRACERS_EAGLE_H
#define SWIFT_TRACERS_EAGLE_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "black_holes.h"
#include "cooling.h"
#include "engine.h"
#include "part.h"
#include "star_formation.h"
#include "tracers_struct.h"

/**
 * @brief Update the particle tracers just after it has been initialised at the
 * start of a step.
 *
 * Nothing to do here in the EAGLE model.
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
 * In EAGLE we record the highest temperature reached and the average SFR.
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
    const double time_step_length, const int *const tracers_triggers_started) {

  /* Current temperature */
  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, p, xp);

  /* New record? */
  if (temperature > xp->tracers_data.maximum_temperature) {

    xp->tracers_data.maximum_temperature = temperature;

    if (with_cosmology) {
      xp->tracers_data.maximum_temperature_scale_factor = cosmo->a;
    } else {
      xp->tracers_data.maximum_temperature_time = time;
    }
  }

  /* Accumulate average SFR */
  for (int i = 0; i < num_snapshot_triggers_part; ++i) {
    if (tracers_triggers_started[i])
      xp->tracers_data.averaged_SFR[i] +=
          star_formation_get_SFR(p, xp) * time_step_length;
  }
}

/**
 * @brief Update the star particle tracers just after its time-step has been
 * computed.
 *
 * In EAGLE, nothing to do.
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
 * In EAGLE, we record the average accr. rate.
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
    const int *const tracers_triggers_started) {

  const float accr_rate = black_holes_get_accretion_rate(bp);

  /* Accumulate average accretion rate */
  for (int i = 0; i < num_snapshot_triggers_part; ++i) {
    if (tracers_triggers_started[i])
      bp->tracers_data.averaged_accretion_rate[i] +=
          accr_rate * time_step_length;
  }
}

/**
 * @brief Initialise the tracer data at the start of a calculation.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 */
static INLINE void tracers_first_init_xpart(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling) {

  xp->tracers_data.maximum_temperature = -1.f;
  xp->tracers_data.maximum_temperature_time = -1.f;
  xp->tracers_data.hit_by_SNII_feedback = 0;
  xp->tracers_data.hit_by_AGN_feedback = 0;
  xp->tracers_data.AGN_feedback_energy = 0.f;

  xp->tracers_data.density_before_last_AGN_feedback_event = -1.f;
  xp->tracers_data.entropy_before_last_AGN_feedback_event = -1.f;
  xp->tracers_data.density_at_last_AGN_feedback_event = -1.f;
  xp->tracers_data.entropy_at_last_AGN_feedback_event = -1.f;

  xp->tracers_data.last_AGN_injection_scale_factor = -1.f;
  xp->tracers_data.density_at_last_AGN_feedback_event = -1.f;

  xp->tracers_data.hit_by_jet_feedback = 0;
  xp->tracers_data.jet_feedback_energy = 0.f;
  xp->tracers_data.last_AGN_jet_feedback_scale_factor = 0.f;
  xp->tracers_data.last_AGN_jet_feedback_time = 0.f;
  xp->tracers_data.last_jet_kick_velocity = 0.f;
}

/**
 * @brief Initialise the star tracer data at the start of a calculation.
 *
 * In EAGLE, nothing to do.
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
                                            const struct cosmology *cosmo) {
  for (int i = 0; i < num_snapshot_triggers_bpart; ++i)
    bp->tracers_data.averaged_accretion_rate[i] = 0.f;
}

/**
 * @brief Update the particles' tracer data after a stellar feedback
 * event.
 *
 * @param xp The extended particle data.
 */
static INLINE void tracers_after_feedback(struct xpart *xp) {

  xp->tracers_data.hit_by_SNII_feedback++;
}

/**
 * @brief Update the particles' tracer data with values before an AGN feedback
 * event. Note: this function is called in `black_holes_iact.h` before the
 * particle data are updated.
 *
 * @param p Pointer to the basic particle data.
 * @param xp The extended particle data.
 * (internal physical units)
 */
static INLINE void tracers_before_black_holes_feedback(
    const struct part *p, struct xpart *xp, const float scale_factor) {

  xp->tracers_data.density_before_last_AGN_feedback_event =
      hydro_get_comoving_density(p) /
      (scale_factor * scale_factor * scale_factor);

  /* Physical entropy (NB entropy has no scale-factor dependence) */
  xp->tracers_data.entropy_before_last_AGN_feedback_event =
      hydro_get_comoving_entropy(p, xp);
}

/**
 * @brief Update the particles' tracer data after an AGN feedback
 * event.
 *
 * @param xp The extended particle data.
 * @param with_cosmology Are we running with cosmology?
 * @param scale_factor The current scale-factor (if running with cosmo)
 * @param time The current time (if running without cosmo)
 * @param delta_energy Amount of energy injected in the feedback event
 * (internal physical units)
 */
static INLINE void tracers_after_black_holes_feedback(
    const struct part *p, struct xpart *xp, const int with_cosmology,
    const float scale_factor, const double time, const double delta_energy) {

  if (with_cosmology)
    xp->tracers_data.last_AGN_injection_scale_factor = scale_factor;
  else
    xp->tracers_data.last_AGN_injection_time = time;

  xp->tracers_data.density_at_last_AGN_feedback_event =
      hydro_get_comoving_density(p) /
      (scale_factor * scale_factor * scale_factor);

  /* Physical entropy (NB entropy has no scale-factor dependence) */
  xp->tracers_data.entropy_at_last_AGN_feedback_event =
      hydro_get_comoving_entropy(p, xp);

  xp->tracers_data.hit_by_AGN_feedback++;
  xp->tracers_data.AGN_feedback_energy += delta_energy;
}

static INLINE void tracers_after_jet_feedback(
    const struct part *p, struct xpart *xp, const int with_cosmology,
    const float scale_factor, const double time, const double delta_energy,
    const float vel_kick) {

  if (with_cosmology)
    xp->tracers_data.last_AGN_jet_feedback_scale_factor = scale_factor;
  else
    xp->tracers_data.last_AGN_jet_feedback_time = time;
  xp->tracers_data.hit_by_jet_feedback++;
  xp->tracers_data.jet_feedback_energy += delta_energy;
  xp->tracers_data.last_jet_kick_velocity = vel_kick;
}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
static INLINE void tracers_after_snapshot_part(const struct part *p,
                                               struct xpart *xp) {

  for (int i = 0; i < num_snapshot_triggers_part; ++i)
    xp->tracers_data.averaged_SFR[i] = 0.f;
}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * @param sp the #spart.
 */
static INLINE void tracers_after_snapshot_spart(struct spart *sp) {

  for (int i = 0; i < num_snapshot_triggers_part; ++i)
    sp->tracers_data.averaged_SFR[i] = 0.f;
}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * @param bp the #bpart.
 */
static INLINE void tracers_after_snapshot_bpart(struct bpart *bp) {

  for (int i = 0; i < num_snapshot_triggers_bpart; ++i)
    bp->tracers_data.averaged_accretion_rate[i] = 0.f;
}

/**
 * @brief Split the tracer content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void tracers_split_part(
    struct part *p, struct xpart *xp, const double n) {

  xp->tracers_data.AGN_feedback_energy /= n;
}

#endif /* SWIFT_TRACERS_EAGLE_H */
