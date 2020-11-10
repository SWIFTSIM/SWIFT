/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_QLA_STAR_FORMATION_H
#define SWIFT_QLA_STAR_FORMATION_H

/* Local includes */
#include "cosmology.h"
#include "engine.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "stars.h"
#include "units.h"

/**
 * @file src/star_formation/EAGLE/star_formation.h
 * @brief Star formation model used in the EAGLE model
 */

#define star_formation_need_update_dx_max 0

/**
 * @brief Properties of the EAGLE star formation model.
 */
struct star_formation {

  /*! Critical overdensity */
  double over_density;
};

/**
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * @param entropy_floor_props The entropy floor assumed in this run.
 */
INLINE static int star_formation_is_star_forming(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* entropy_floor_props) {

  /* Minimal density (converted from mean baryonic density)
     for star formation */
  const double rho_mean_b_times_min_over_den =
      cosmo->mean_density_Omega_b * starform->over_density;

  /* Physical density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Verify whether we are above the over-density threshold */
  return (physical_density > rho_mean_b_times_min_over_den);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart.
 *
 * Nothing to do here. Particles that pass the SF criterion get automcatically
 * converted to a star. No need to compute or store a star formation rate.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    const struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

  xp->sf_data.convert_to_star = 1;
}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.
 *
 * Equation 21 of Schaye & Dalla Vecchia 2008.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct engine* e,
    const double dt_star) {

  return xp->sf_data.convert_to_star;
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * Nothing to do here in the quick Lyman-alpha model.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct star_formation* starform, const int with_cosmology) {}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 * @param starform the star formation law properties to use.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * @param convert_part Did we convert a part(or spawned one)?
 */
INLINE static void star_formation_copy_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cooling_function_data* cooling,
    const int convert_part) {

  /* Store the current mass */
  sp->mass = hydro_get_mass(p);

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units.
 * @param hydro_props The propertis of the hydro model.
 * @param cosmo The current cosmological model.
 * @param entropy_floor The properties of the entropy floor used in this
 * simulation.
 * @param starform the star formation law properties to initialize
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    const struct cosmology* cosmo,
    const struct entropy_floor_properties* entropy_floor,
    struct star_formation* starform) {

  /* Read the critical density contrast from the parameter file*/
  starform->over_density =
      parser_get_param_double(parameter_file, "QLAStarFormation:over_density");
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  message("Star formation law is Quick Lyman-alpha");
  message("Over-density for star formation: %f", starform->over_density);
}

/**
 * @brief Finishes the density calculation.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the quick Lyman-alpha star formation model.
 *
 * @param p The particle to act upon
 * @param xp The extra particle to act upon
 * @param cd The global star_formation information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_formation_end_density(
    struct part* p, struct xpart* xp, const struct star_formation* cd,
    const struct cosmology* cosmo) {}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the quick Lyman-alpha star formation model.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(struct part* p, struct xpart* xp,
                                      const struct star_formation* cd,
                                      const struct cosmology* cosmo) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * state to start the density loop.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the quick Lyman-alpha star formation model.
 *
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    struct part* p, const struct star_formation* data) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state at the beginning of the simulation after the ICs have been read.
 *
 * Mark the particles as not needing to be converted to stars.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_part(const struct phys_const* phys_const,
                               const struct unit_system* us,
                               const struct cosmology* cosmo,
                               const struct star_formation* data,
                               const struct part* p, struct xpart* xp) {

  xp->sf_data.convert_to_star = 0;
}

/**
 * @brief Split the star formation content of a particle into n pieces
 *
 * Nothing to do here.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void star_formation_split_part(
    struct part* p, struct xpart* xp, const double n) {}

/**
 * @brief Deal with the case where no spart are available for star formation.
 *
 * @param e The #engine.
 * @param p The #part.
 * @param xp The #xpart.
 */
__attribute__((always_inline)) INLINE static void
star_formation_no_spart_available(const struct engine* e, const struct part* p,
                                  const struct xpart* xp) {
  /* Nothing to do since we just turn gas particles into DM. */
}

/**
 * @brief Decides whether a new particle should be created or if the hydro
 * particle needs to be transformed.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 *
 * @return 1 if a new spart needs to be created.
 */
INLINE static int star_formation_should_spawn_spart(
    struct part* p, struct xpart* xp, const struct star_formation* starform) {
  return 0;
}

/**
 * @brief Compute some information for the star formation model based
 * on all the particles that were read in.
 *
 * This is called once on start-up of the code.
 *
 * Nothing to do here for the quick Lyman-alpha model.
 *
 * @param star_form The #star_formation structure.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_stats(struct star_formation* star_form,
                                const struct engine* e) {}

#endif /* SWIFT_QLA_STAR_FORMATION_H */
