/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
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
#ifndef SWIFT_FORCING_BOUNDARY_PARTICLES_H
#define SWIFT_FORCING_BOUNDARY_PARTICLES_H

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <float.h>

/* Local includes. */
#include "engine.h"
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief Forcing Term Properties
 */
struct forcing_terms {

  /*! Particles with ID <= boundary_particle_max_id are treated as boundary
   * particles. */
  int boundary_particle_max_id;

  /*! Whether or not boundary particles are fixed in place (0 (false, default)
   * or 1 (true))  */
  int enable_fixed_position;

  /*! Whether or not boundary particles have hydrodynamic accelerations enabled
   * (0 (false, default) or 1 (true))  */
  int enable_hydro_acceleration;

  /*! Whether or not boundary particles have gravitational accelerations enabled
   * (0 (false, default) or 1 (true))  */
  int enable_grav_acceleration;
};

/**
 * @brief Computes the hydrodynamic forcing terms.
 *
 * Resets the hydrodynamic accelerations of boundary particles to zero.
 * If using fixed boundary particles, their velocities are also set to zero.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param s The #space we act on.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_hydro_terms_apply(
    const double time, const struct forcing_terms *terms, const struct space *s,
    const struct phys_const *phys_const, struct part *p, struct xpart *xp) {

  if (terms->enable_hydro_acceleration) {
    /* Skip if not resetting hydro accelerations. */
    return;
  }

  if (p->id <= terms->boundary_particle_max_id) {
    /* Reset the hydro accelerations of boundary particles. */
    p->a_hydro[0] = 0.0f;
    p->a_hydro[1] = 0.0f;
    p->a_hydro[2] = 0.0f;

#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)

    /* Some values need to be reset in the Gizmo case. */
    hydro_prepare_force(p, &c->hydro.xparts[k], cosmo, e->hydro_properties,
                        e->pressure_floor_props,
                        /*dt_alpha=*/0, /*dt_therm=*/0);
    rt_prepare_force(p);
#endif

    if (terms->enable_fixed_position) {
      /* Set velocity of fixed boundary particle to zero. */
      xp->v_full[0] = 0.f;
      xp->v_full[1] = 0.f;
      xp->v_full[2] = 0.f;
    }
  }
}

/**
 * @brief Computes the gravitational forcing terms.
 *
 * Resets the gravitational accelerations of boundary particles to zero.
 * If using fixed boundary particles, their velocities are also set to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_grav_terms_apply(
    const long long id, const struct forcing_terms *terms, struct gpart *gp) {

  if (terms->enable_grav_acceleration) {
    /* Skip if not resetting grav accelerations. */
    return;
  }

  if (id <= terms->boundary_particle_max_id) {
    /* Reset the grav accelerations of boundary particles. */
    gp->a_grav[0] = 0.0f;
    gp->a_grav[1] = 0.0f;
    gp->a_grav[2] = 0.0f;

    if (terms->enable_fixed_position) {
      /* Set velocity of fixed boundary particle to zero. */
      gp->v_full[0] = 0.f;
      gp->v_full[1] = 0.f;
      gp->v_full[2] = 0.f;
    }
  }
}

/**
 * @brief Sets the forcing of gparts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_gpart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct gpart *gp) {

  if (id <= terms->boundary_particle_max_id && terms->enable_fixed_position) {
    /* Set velocity of fixed boundary particle to zero. */
    gp->v_full[0] = 0.f;
    gp->v_full[1] = 0.f;
    gp->v_full[2] = 0.f;
  }
}

/**
 * @brief Sets the forcing of parts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_part_drift_apply(
    const long long id, const struct forcing_terms *terms, struct part *p,
    struct xpart *xp) {

  if (id <= terms->boundary_particle_max_id && terms->enable_fixed_position) {
    /* Set velocity of fixed boundary particle to zero. */
    xp->v_full[0] = 0.f;
    xp->v_full[1] = 0.f;
    xp->v_full[2] = 0.f;
  }
}

/**
 * @brief Sets the forcing of sparts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_spart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct spart *sp) {

  if (id <= terms->boundary_particle_max_id && terms->enable_fixed_position) {
    /* Set velocity of fixed boundary particle to zero. */
    sp->v[0] = 0.f;
    sp->v[1] = 0.f;
    sp->v[2] = 0.f;
  }
}

/**
 * @brief Sets the forcing of bparts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param bp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_bpart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct bpart *bp) {

  if (id <= terms->boundary_particle_max_id && terms->enable_fixed_position) {
    /* Set velocity of fixed boundary particle to zero. */
    bp->v[0] = 0.f;
    bp->v[1] = 0.f;
    bp->v[2] = 0.f;
  }
}

/**
 * @brief Computes the time-step condition due to the forcing terms.
 *
 * Nothing to do here. --> Return FLT_MAX.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static float forcing_terms_timestep(
    double time, const struct forcing_terms *terms,
    const struct phys_const *phys_const, const struct part *p,
    const struct xpart *xp) {

  return FLT_MAX;
}

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms *terms) {

  message(
      "Forcing terms is 'Boundary particles'. Max boundary particle ID: %i.",
      terms->boundary_particle_max_id);
  message(
      "Run with options: Enable fixed position: %i, Enable hydro acceleration: "
      "%i, Enable grav acceleration: %i",
      terms->enable_fixed_position, terms->enable_hydro_acceleration,
      terms->enable_grav_acceleration);
}

/**
 * @brief Initialises the forcing term properties
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space object.
 * @param terms The forcing term properties to initialize
 */
static INLINE void forcing_terms_init(struct swift_params *parameter_file,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct space *s,
                                      struct forcing_terms *terms) {

  terms->boundary_particle_max_id = parser_get_param_int(
      parameter_file, "BoundaryParticles:boundary_particle_max_id");
  terms->enable_fixed_position = parser_get_opt_param_int(
      parameter_file, "BoundaryParticles:enable_fixed_position", 0);
  terms->enable_hydro_acceleration = parser_get_opt_param_int(
      parameter_file, "BoundaryParticles:enable_hydro_acceleration", 0);
  terms->enable_grav_acceleration = parser_get_opt_param_int(
      parameter_file, "BoundaryParticles:enable_grav_acceleration", 0);

  if (terms->enable_fixed_position != 0 && terms->enable_fixed_position != 1) {
    error(
        "BoundaryParticles:enable_fixed_position must be either 0 (false) or 1 "
        "(true).");
  }

  if (terms->enable_fixed_position) {
    /* If using fixed boundary particles, both hydro and grav accelerations must
     * be reset to 0. */
    if (terms->enable_hydro_acceleration != 0) {
      error(
          "BoundaryParticles:enable_hydro_acceleration must be 0 (false) when "
          "using fixed boundary particles.");
    }
    if (terms->enable_grav_acceleration != 0) {
      error(
          "BoundaryParticles:enable_grav_acceleration must be 0 (false) when "
          "using fixed boundary particles.");
    }

  } else {
    /* If not using fixed boundary particles, hydro and grav accelerations can
     * be enabled independently. */
    if (terms->enable_hydro_acceleration != 0 &&
        terms->enable_hydro_acceleration != 1) {
      error(
          "BoundaryParticles:enable_hydro_acceleration must be either 0 "
          "(false) or 1 (true).");
    }
    if (terms->enable_grav_acceleration != 0 &&
        terms->enable_grav_acceleration != 1) {
      error(
          "BoundaryParticles:enable_grav_acceleration must be either 0 (false) "
          "or 1 (true).");
    }
  }
}

#endif /* SWIFT_FORCING_BOUNDARY_PARTICLES_H */
