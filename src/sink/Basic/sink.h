/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
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
#ifndef SWIFT_BASIC_SINK_H
#define SWIFT_BASIC_SINK_H

#include <float.h>

/* Put pragma if gsl around here */
#ifdef HAVE_LIBGSL
#include <gsl/gsl_cdf.h>
#endif

/* Local includes */
#include "active.h"
#include "chemistry.h"
#include "cooling.h"
#include "feedback.h"
#include "minmax.h"
#include "random.h"
#include "sink_part.h"
#include "sink_properties.h"
#include "star_formation.h"

/**
 * @brief Computes the time-step of a given sink particle.
 *
 * @param sp Pointer to the sink-particle data.
 */
__attribute__((always_inline)) INLINE static float sink_compute_timestep(
    const struct sink *const sink, const struct sink_props *sink_properties,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props, const double time,
    const double time_base) {

  return FLT_MAX;
}

/**
 * @brief Initialises the sink-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The #sink particle to act upon.
 * @param sink_props The properties of the sink particles scheme.
 * @param e The #engine
 */
__attribute__((always_inline)) INLINE static void sink_first_init_sink(
    struct sink *sp, const struct sink_props *sink_props,
    const struct engine *e) {

  sp->time_bin = 0;

  sp->number_of_gas_swallows = 0;
  sp->number_of_direct_gas_swallows = 0;
  sp->number_of_sink_swallows = 0;
  sp->number_of_direct_sink_swallows = 0;
  sp->swallowed_angular_momentum[0] = 0.f;
  sp->swallowed_angular_momentum[1] = 0.f;
  sp->swallowed_angular_momentum[2] = 0.f;

  /* Initially set the subgrid mass equal to the dynamical mass read from the
   * ICs. */
  sp->subgrid_mass = sp->mass;

  sink_mark_sink_as_not_swallowed(&sp->merger_data);
}

/**
 * @brief Initialisation of particle data before the hydro density loop.
 * Note: during initalisation (space_init)
 *
 * @param p The #part to act upon.
 * @param sink_props The properties of the sink particles scheme.
 */
__attribute__((always_inline)) INLINE static void sink_init_part(
    struct part *restrict p, const struct sink_props *sink_props) {}

/**
 * @brief Initialisation of sink particle data before sink loops.
 * Note: during initalisation (space_init_sinks)
 *
 * @param sp The #sink particle to act upon.
 */
__attribute__((always_inline)) INLINE static void sink_init_sink(
    struct sink *sp) {

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
  sp->rho_gas = 0.f;
  sp->sound_speed_gas = 0.f;
  sp->velocity_gas[0] = 0.f;
  sp->velocity_gas[1] = 0.f;
  sp->velocity_gas[2] = 0.f;
  sp->ngb_mass = 0.f;
  sp->num_ngbs = 0;
  sp->accretion_rate = 0.f;
  sp->mass_at_start_of_step = sp->mass; /* sp->mass may grow in nibbling mode */

#ifdef DEBUG_INTERACTIONS_SINKS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_accretion[i] = -1;
  sp->num_ngb_accretion = 0;

  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_merger[i] = -1;
  sp->num_ngb_merger = 0;

  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_formation[i] = -1;
  sp->num_ngb_formation = 0;
#endif

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
    struct sink *restrict sp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The #sink particle.
 */
__attribute__((always_inline)) INLINE static void sink_reset_predicted_values(
    struct sink *restrict sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The #sink particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void sink_kick_extra(
    struct sink *sp, float dt) {}

/**
 * @brief Finishes the calculation of density on sinks
 *
 * @param si The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_end_density(
    struct sink *si, const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = si->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  si->density.wcount *= h_inv_dim;
  si->density.wcount_dh *= h_inv_dim_plus_one;

  /* Finish the density calculation */
  si->rho_gas *= h_inv_dim;

  /* For the following, we also have to undo the mass smoothing
   * (N.B.: bp->velocity_gas is in BH frame, in internal units). */
  const float rho_inv = 1.f / si->rho_gas;
  si->sound_speed_gas *= h_inv_dim * rho_inv;
  si->velocity_gas[0] *= h_inv_dim * rho_inv;
  si->velocity_gas[1] *= h_inv_dim * rho_inv;
  si->velocity_gas[2] *= h_inv_dim * rho_inv;

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
    struct sink *restrict sp, const struct cosmology *cosmo) {

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
    struct sink *restrict si, const struct sink_props *props,
    const struct phys_const *constants, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling,
    const struct entropy_floor_properties *floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {

  if (dt == 0. || si->rho_gas == 0.) return;

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;

  /* (Subgrid) mass of the sink (internal units) */
  const double sink_mass = si->subgrid_mass;

  /* We can now compute the accretion rate (internal units) */
  /* Standard approach: compute accretion rate for all gas simultaneously.
   *
   * Convert the quantities we gathered to physical frame (all internal
   * units). Note: velocities are already in black hole frame. */
  const double gas_rho_phys = si->rho_gas * cosmo->a3_inv;
  const double gas_c_phys = si->sound_speed_gas * cosmo->a_factor_sound_speed;
  const double gas_c_phys2 = gas_c_phys * gas_c_phys;
  const double gas_v_phys[3] = {si->velocity_gas[0] * cosmo->a_inv,
                                si->velocity_gas[1] * cosmo->a_inv,
                                si->velocity_gas[2] * cosmo->a_inv};
  const double gas_v_norm2 = gas_v_phys[0] * gas_v_phys[0] +
                             gas_v_phys[1] * gas_v_phys[1] +
                             gas_v_phys[2] * gas_v_phys[2];

  const double denominator2 = gas_v_norm2 + gas_c_phys2;

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure that the denominator is strictly positive */
  if (denominator2 <= 0)
    error(
        "Invalid denominator for sink particle %lld in Bondi rate "
        "calculation.",
        si->id);
#endif

  const double denominator_inv = 1. / sqrt(denominator2);

  double accr_rate = 4. * M_PI * G * G * sink_mass * sink_mass * gas_rho_phys *
                     denominator_inv * denominator_inv * denominator_inv;

  /* Integrate forward in time */
  si->subgrid_mass += accr_rate * dt;
}

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
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct sink_props *sink_props, const struct phys_const *phys_const,
    const struct cosmology *cosmo,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cooling_function_data *restrict cooling,
    const struct entropy_floor_properties *restrict entropy_floor) {

  /* Sink formation is not implemented in this model. */
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
    const struct part *p, const struct xpart *xp,
    const struct sink_props *sink_props, const struct engine *e,
    const double dt_sink) {

  return 0;
}

/**
 * @brief Copies the properties of the gas particle over to the
 * sink particle.
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
    const struct part *p, const struct xpart *xp, struct sink *sink,
    const struct engine *e, const struct sink_props *sink_props,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct phys_const *phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cooling_function_data *restrict cooling) {}

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
    struct sink *sp, const struct part *p, const struct xpart *xp,
    const struct cosmology *cosmo) {

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float sink_mass = sp->mass;

  /* Increase the dynamical mass of the sink. */
  sp->mass += gas_mass;
  sp->gpart->mass += gas_mass;

  /* Comoving and physical distance between the particles */
  const float dx[3] = {sp->x[0] - p->x[0], sp->x[1] - p->x[1],
                       sp->x[2] - p->x[2]};
  const float dx_physical[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a,
                                dx[2] * cosmo->a};

  /* Relative velocity between the sink and the part */
  const float dv[3] = {sp->v[0] - p->v[0], sp->v[1] - p->v[1],
                       sp->v[2] - p->v[2]};

  const float a = cosmo->a;
  const float H = cosmo->H;
  const float a2H = a * a * H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                  a2H * dx[2] + dv[2]};

  /* Compute the physical relative velocity between the particles */
  const float dv_physical[3] = {v_plus_H_flow[0] * cosmo->a_inv,
                                v_plus_H_flow[1] * cosmo->a_inv,
                                v_plus_H_flow[2] * cosmo->a_inv};

  /* Collect the swallowed angular momentum */
  sp->swallowed_angular_momentum[0] +=
      gas_mass *
      (dx_physical[1] * dv_physical[2] - dx_physical[2] * dv_physical[1]);
  sp->swallowed_angular_momentum[1] +=
      gas_mass *
      (dx_physical[2] * dv_physical[0] - dx_physical[0] * dv_physical[2]);
  sp->swallowed_angular_momentum[2] +=
      gas_mass *
      (dx_physical[0] * dv_physical[1] - dx_physical[1] * dv_physical[0]);

  /* Update the sink momentum */
  const float sink_mom[3] = {sink_mass * sp->v[0] + gas_mass * p->v[0],
                             sink_mass * sp->v[1] + gas_mass * p->v[1],
                             sink_mass * sp->v[2] + gas_mass * p->v[2]};

  sp->v[0] = sink_mom[0] / sp->mass;
  sp->v[1] = sink_mom[1] / sp->mass;
  sp->v[2] = sink_mom[2] / sp->mass;
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];

  /* This sink swallowed a gas particle */
  sp->number_of_gas_swallows++;
  sp->number_of_direct_gas_swallows++;

#ifdef SWIFT_DEBUG_CHECKS
  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message(
      "sink %lld swallow gas particle %lld. "
      "(Mass = %e, "
      "Delta_v = [%f, %f, %f] U_V, "
      "Delta_x = [%f, %f, %f] U_L, "
      "Delta_v_rad = %f)",
      sp->id, p->id, sp->mass, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);
#endif
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a sink particle.
 *
 * @param spi The #sink to update.
 * @param spj The #sink that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_sink(
    struct sink *spi, const struct sink *spj, const struct cosmology *cosmo) {

  /* Get the current dynamical masses */
  const float spi_dyn_mass = spi->mass;
  const float spj_dyn_mass = spj->mass;

  /* Increase the masses of the sink. */
  spi->mass += spj->mass;
  spi->gpart->mass += spj->mass;
  spi->subgrid_mass += spj->subgrid_mass;

  /* Collect the swallowed angular momentum */
  spi->swallowed_angular_momentum[0] += spj->swallowed_angular_momentum[0];
  spi->swallowed_angular_momentum[1] += spj->swallowed_angular_momentum[1];
  spi->swallowed_angular_momentum[2] += spj->swallowed_angular_momentum[2];

  /* Update the sink momentum */
  const float sink_mom[3] = {
      spi_dyn_mass * spi->v[0] + spj_dyn_mass * spj->v[0],
      spi_dyn_mass * spi->v[1] + spj_dyn_mass * spj->v[1],
      spi_dyn_mass * spi->v[2] + spj_dyn_mass * spj->v[2]};

  spi->v[0] = sink_mom[0] / spi->mass;
  spi->v[1] = sink_mom[1] / spi->mass;
  spi->v[2] = sink_mom[2] / spi->mass;
  spi->gpart->v_full[0] = spi->v[0];
  spi->gpart->v_full[1] = spi->v[1];
  spi->gpart->v_full[2] = spi->v[2];

  /* This sink swallowed a sink particle */
  spi->number_of_sink_swallows++;
  spi->number_of_direct_sink_swallows++;

  /* Add all other swallowed particles swallowed by the swallowed sink */
  spi->number_of_sink_swallows += spj->number_of_sink_swallows;
  spi->number_of_gas_swallows += spj->number_of_gas_swallows;

  message("sink %lld swallow sink particle %lld. New mass: %e.", spi->id,
          spj->id, spi->mass);
}

/**
 * @brief Should the sink spawn a star particle?
 *
 * @param sink the sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param cosmo The cosmological parameters and properties.
 * @param with_cosmology If we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 */
INLINE static int sink_spawn_star(struct sink *sink, const struct engine *e,
                                  const struct sink_props *sink_props,
                                  const struct cosmology *cosmo,
                                  const int with_cosmology,
                                  const struct phys_const *phys_const,
                                  const struct unit_system *restrict us) {

  /* Star formation from sinks is disabled in this model. */
  return 0;
}

/**
 * @brief Copy the properties of the sink particle towards the new star. Also,
 * give the stars some properties such as position and velocity.
 *
 * This function also needs to update the sink particle.
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
    struct sink *sink, struct spart *sp, const struct engine *e,
    const struct sink_props *sink_props, const struct cosmology *cosmo,
    const int with_cosmology, const struct phys_const *phys_const,
    const struct unit_system *restrict us) {}

/**
 * @brief Update the #sink particle properties before spawning a star.
 *
 * In GEAR, we check if the sink had an IMF change from pop III to pop II
 * during the last gas/sink accretion loops. If so, we draw a new target mass
 * with the correct IMF so that stars have the right metallicities.
 *
 * @param sink The #sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param phys_const The physical constants in internal units.
 */
INLINE static void sink_update_sink_properties_before_star_formation(
    struct sink *sink, const struct engine *e,
    const struct sink_props *sink_props, const struct phys_const *phys_const) {}

/**
 * @brief Update the #sink particle properties right after spawning a star.
 *
 * In GEAR: Important properties that are updated are the sink mass and the
 * sink->target_mass_Msun to draw the next star mass.
 *
 * @param sink The #sink particle that spawed stars.
 * @param sp The #spart particle spawned.
 * @param e The #engine
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param star_counter The star loop counter.
 */
INLINE static void sink_update_sink_properties_during_star_formation(
    struct sink *sink, const struct spart *sp, const struct engine *e,
    const struct sink_props *sink_props, const struct phys_const *phys_const,
    int star_counter) {}

/**
 * @brief Update the #sink particle properties after star formation.
 *
 * In GEAR, this is unused.
 *
 * @param sink The #sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param phys_const The physical constants in internal units.
 */
INLINE static void sink_update_sink_properties_after_star_formation(
    struct sink *sink, const struct engine *e,
    const struct sink_props *sink_props, const struct phys_const *phys_const) {}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The sink data of a gas particle.
 * @param gp The part's #gpart.
 */
__attribute__((always_inline)) INLINE static void sink_store_potential_in_part(
    struct sink_part_data *p_data, const struct gpart *gp) {}

/**
 * @brief Compute all quantities required for the formation of a sink such as
 * kinetic energy, potential energy, etc. This function works on the
 * neighbouring gas particles.
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
    struct engine *e, struct part *restrict p, struct xpart *restrict xp,
    struct part *restrict pi, struct xpart *restrict xpi,
    const struct cosmology *cosmo, const struct sink_props *sink_props) {}

/**
 * @brief Compute all quantities required for the formation of a sink. This
 * function works on the neighbouring sink particles.
 *
 * @param e The #engine.
 * @param p The #part for which we compute the quantities.
 * @param xp The #xpart data of the particle #p.
 * @param si A neighbouring #sink of #p.
 * @param cosmo The cosmological parameters and properties.
 * @param sink_props The sink properties to use.
 */
INLINE static void sink_prepare_part_sink_formation_sink_criteria(
    struct engine *e, struct part *restrict p, struct xpart *restrict xp,
    struct sink *restrict si, const int with_cosmology,
    const struct cosmology *cosmo, const struct sink_props *sink_props,
    const double time) {}

#endif /* SWIFT_BASIC_SINK_H */
