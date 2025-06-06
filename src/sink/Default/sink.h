/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_SINK_H
#define SWIFT_DEFAULT_SINK_H

#include <float.h>

/* Local includes */
#include "minmax.h"
#include "random.h"
#include "sink_part.h"
#include "sink_properties.h"

/**
 * @brief Computes the time-step of a given sink particle.
 *
 * Note: In the Default sink, no time-step limit is applied.
 *
 * @param sink Pointer to the sink-particle data.
 * @param sink_properties Properties of the sink model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param cosmo The current cosmological model (used if running with
 * cosmology).
 * @param grav_props The current gravity properties.
 * @param time The current time (used if running without cosmology).
 * @param time_base The time base.
 */
__attribute__((always_inline)) INLINE static float sink_compute_timestep(
    const struct sink* const sink, const struct sink_props* sink_properties,
    const int with_cosmology, const struct cosmology* cosmo,
    const struct gravity_props* grav_props, const double time,
    const double time_base) {

  return FLT_MAX;
}

/**
 * @brief Initialises the sink-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 * @param sink_props The properties of the sink particles scheme.
 * @param e The #engine
 */
__attribute__((always_inline)) INLINE static void sink_first_init_sink(
    struct sink* sp, const struct sink_props* sink_props,
    const struct engine* e) {}

/**
 * @brief Prepares a particle for the sink calculation.
 *
 * @param p The #part to act upon.
 * @param sink_props The properties of the sink particles scheme.
 */
__attribute__((always_inline)) INLINE static void sink_init_part(
    struct part* restrict p, const struct sink_props* sink_props) {}

/**
 * @brief Prepares a sink-particle for its interactions
 *
 * @param sp The #sink particle to act upon.
 */
__attribute__((always_inline)) INLINE static void sink_init_sink(
    struct sink* sp) {

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;

#ifdef SWIFT_SINK_DENSITY_CHECKS
  sp->N_check_density = 0;
  sp->N_check_density_exact = 0;
  sp->rho_check = 0.f;
  sp->rho_check_exact = 0.f;
  sp->n_check = 0.f;
  sp->n_check_exact = 0.f;
  sp->inhibited_check_exact = 0;
#endif
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sp The #sink.
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void sink_predict_extra(
    struct sink* restrict sp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The #sink particle.
 */
__attribute__((always_inline)) INLINE static void sink_reset_predicted_values(
    struct sink* restrict sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The #sink particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void sink_kick_extra(
    struct sink* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on sinks
 *
 * @param si The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_end_density(
    struct sink* si, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = si->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  si->density.wcount *= h_inv_dim;
  si->density.wcount_dh *= h_inv_dim_plus_one;

#ifdef SWIFT_SINK_DENSITY_CHECKS
  si->rho_check *= h_inv_dim;
  si->n_check *= h_inv_dim;
#endif
}

/**
 * @brief Sets all particle fields to sensible values when the #sink has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sinks_sink_has_no_neighbours(
    struct sink* restrict sp, const struct cosmology* cosmo) {

  warning(
      "Sink particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      sp->id, sp->h, sp->density.wcount);

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  sp->density.wcount = kernel_root * h_inv_dim;
  sp->density.wcount_dh = 0.f;
}

/**
 * @brief Compute the accretion rate of the sink and any quantities
 * required swallowing based on an accretion rate
 *
 * Adapted from black_holes_prepare_feedback
 *
 * @param si The sink particle.
 * @param props The properties of the sink scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param cooling Properties of the cooling model.
 * @param floor_props Properties of the entropy floor.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 * @param ti_begin Integer time value at the beginning of timestep
 */
__attribute__((always_inline)) INLINE static void sink_prepare_swallow(
    struct sink* restrict si, const struct sink_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {}

/**
 * @brief Calculate if the gas has the potential of becoming
 * a sink.
 *
 * Return 0 if no sink formation should occur.
 * Note: called in runner_do_sink_formation
 *
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * @param entropy_floor The entropy_floor properties.
 *
 */
INLINE static int sink_is_forming(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct sink_props* sink_props, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const struct entropy_floor_properties* restrict entropy_floor) {

  return 0;
}

/**
 * @brief Decides whether a particle should be converted into a
 * sink or not.
 *
 * No SF should occur, so return 0.
 * Note: called in runner_do_sink_formation
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param sink_props The properties of the sink model.
 * @param e The #engine (for random numbers).
 * @param dt_sink The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int sink_should_convert_to_sink(
    const struct part* p, const struct xpart* xp,
    const struct sink_props* sink_props, const struct engine* e,
    const double dt_sink) {

  return 0;
}

/**
 * @brief Copies the properties of the gas particle over to the
 * sink particle.
 *
 * Nothing to do here.
 *
 * @param p The gas particles.
 * @param xp The additional properties of the gas particles.
 * @param sink the new created #sink particle.
 * @param e The #engine.
 * @param sink_props The sink properties to use.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The hydro properties to use.
 * @param us The internal unit system.
 * @param cooling The cooling function to use.
 */
INLINE static void sink_copy_properties(
    const struct part* p, const struct xpart* xp, struct sink* sink,
    const struct engine* e, const struct sink_props* sink_props,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* Set a smoothing length */
  sink->h = p->h;
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a gas particle.
 *
 * @param sp The #sink to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_part(
    struct sink* sp, const struct part* p, const struct xpart* xp,
    const struct cosmology* cosmo) {}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a sink particle.
 *
 * @param spi The #sink to update.
 * @param spj The #sink that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_sink(
    struct sink* spi, const struct sink* spj, const struct cosmology* cosmo) {}

/**
 * @brief Should the sink spawn a star particle?
 *
 * Nothing to do here.
 *
 * @param sink the sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param cosmo The cosmological parameters and properties.
 * @param with_cosmology If we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 */
INLINE static int sink_spawn_star(struct sink* sink, const struct engine* e,
                                  const struct sink_props* sink_props,
                                  const struct cosmology* cosmo,
                                  const int with_cosmology,
                                  const struct phys_const* phys_const,
                                  const struct unit_system* restrict us) {
  return 0;
}

/**
 * @brief Copy the properties of the sink particle towards the new star.
 * This function also needs to update the sink particle.
 *
 * Nothing to do here.
 *
 * @param sink The #sink particle.
 * @param sp The star particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param cosmo The cosmological parameters and properties.
 * @param with_cosmology If we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 */
INLINE static void sink_copy_properties_to_star(
    struct sink* sink, struct spart* sp, const struct engine* e,
    const struct sink_props* sink_props, const struct cosmology* cosmo,
    const int with_cosmology, const struct phys_const* phys_const,
    const struct unit_system* restrict us) {}

/**
 * @brief Update the #sink particle properties before spawning a star.
 *
 * @param sink The #sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param phys_const The physical constants in internal units.
 */
INLINE static void sink_update_sink_properties_before_star_formation(
    struct sink* sink, const struct engine* e,
    const struct sink_props* sink_props, const struct phys_const* phys_const) {}

/**
 * @brief Update the #sink particle properties right after spawning a star.
 *
 * @param sink The #sink particle that spawed stars.
 * @param sp The #spart particle spawned.
 * @param e The #engine
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param star_counter The star loop counter.
 */
INLINE static void sink_update_sink_properties_during_star_formation(
    struct sink* sink, const struct spart* sp, const struct engine* e,
    const struct sink_props* sink_props, const struct phys_const* phys_const,
    int star_counter) {}

/**
 * @brief Update the #sink particle properties after star formation.
 *
 * @param sink The #sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param phys_const The physical constants in internal units.
 */
INLINE static void sink_update_sink_properties_after_star_formation(
    struct sink* sink, const struct engine* e,
    const struct sink_props* sink_props, const struct phys_const* phys_const) {}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The sink data of a gas particle.
 * @param gp The part's #gpart.
 */
__attribute__((always_inline)) INLINE static void sink_store_potential_in_part(
    struct sink_part_data* p_data, const struct gpart* gp) {}

/**
 * @brief Compute all quantities required for the formation of a sink such as
 * kinetic energy, potential energy, etc. This function works on the
 * neighbouring gas particles.
 *
 * Nothing to do here.
 *
 * @param e The #engine.
 * @param p The #part for which we compute the quantities.
 * @param xp The #xpart data of the particle #p.
 * @param pi A neighbouring #part of #p.
 * @param xpi The #xpart data of the particle #pi.
 * @param cosmo The cosmological parameters and properties.
 * @param sink_props The sink properties to use.
 */
INLINE static void sink_prepare_part_sink_formation_gas_criteria(
    struct engine* e, struct part* restrict p, struct xpart* restrict xp,
    struct part* restrict pi, struct xpart* restrict xpi,
    const struct cosmology* cosmo, const struct sink_props* sink_props) {}

/**
 * @brief Compute all quantities required for the formation of a sink. This
 * function works on the neighbouring sink particles.
 *
 * Nothing to do here.
 *
 * @param e The #engine.
 * @param p The #part for which we compute the quantities.
 * @param xp The #xpart data of the particle #p.
 * @param si A neighbouring #sink of #p.
 * @param cosmo The cosmological parameters and properties.
 * @param sink_props The sink properties to use.
 */
INLINE static void sink_prepare_part_sink_formation_sink_criteria(
    struct engine* e, struct part* restrict p, struct xpart* restrict xp,
    struct sink* restrict si, const int with_cosmology,
    const struct cosmology* cosmo, const struct sink_props* sink_props,
    const double time) {}

#endif /* SWIFT_DEFAULT_SINK_H */
