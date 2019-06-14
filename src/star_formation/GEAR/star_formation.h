/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
#ifndef SWIFT_GEAR_STAR_FORMATION_H
#define SWIFT_GEAR_STAR_FORMATION_H

/* Local includes */
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "random.h"
#include "star_formation_struct.h"
#include "units.h"

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
INLINE static float get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  /* Physical constants */
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  /* Gas properties */
  const double T_transition = hydro_props->hydrogen_ionization_temperature;
  const double mu_neutral = hydro_props->mu_neutral;
  const double mu_ionised = hydro_props->mu_ionised;

  /* Particle temperature */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Temperature over mean molecular weight */
  const double T_over_mu = hydro_gamma_minus_one * u * m_H / k_B;

  /* Are we above or below the HII -> HI transition? */
  if (T_over_mu > (T_transition + 1.) / mu_ionised)
    return T_over_mu * mu_ionised;
  else if (T_over_mu < (T_transition - 1.) / mu_neutral)
    return T_over_mu * mu_neutral;
  else
    return T_transition;
}
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
 *
 */
INLINE static int star_formation_is_star_forming(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const struct entropy_floor_properties* restrict entropy_floor) {
  /* Some useful constants */
  const double G = phys_const->const_newton_G;
  const double kb = phys_const->const_boltzmann_k;
  const double mH = phys_const->const_proton_mass;
  /* We first check if we are supposed to include turbulence estimation
   * otherewise we keep 0 */
  float sigma2 = 0.f;
  if (starform->with_sigma > 0) {
    sigma2 = xp->sf_data.sigma2;
  }
  /* Compute the temperature */
  double const T =
      get_temperature(phys_const, hydro_props, us, cosmo, cooling, p, xp);
  /* Other useful values */
  const int N = starform->Njeans;
  const double h = p->h;
  /* We suppose that the gas is neutral */
  const double mu = hydro_props->mu_neutral;
  /* Maximum temperature allowed (temperature criterion) */
  const double T0 = starform->Max_temperature;
  /* Compute density */
  const double physical_density = hydro_get_physical_density(p, cosmo);
  /* We compute the minimal density for star formation (see Revaz & Jablonka,
   * 2018 eq (3)) */
  const double rho_sfr = M_PI * (hydro_gamma * kb * T / mu / mH + sigma2) / h /
                         h / 4. / G / pow(N, 2. / 3.);
  /* Temperature criterion for star formation eligibility */
  if (T > T0) {
    return 0;
  }
  /* Density criterion */
  else {
    if (physical_density > rho_sfr) {
      /* The particle is eligible (satisfy both conditions) */
      return 1;
    } else {
      return 0;
    }
  }
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const double dt_star) {
  if (dt_star == 0.) {
    xp->sf_data.SFR = 0.;
  } else {
    const double G = phys_const->const_newton_G;
    const double c_star = starform->star_formation_rate;
    const double physical_density = hydro_get_physical_density(p, cosmo);
    if (physical_density != 0) {
      /* We compute the star formation rate and store it (see Revaz & Jablonka,
       * 2012, eq. (5)) */
      /* Units are mass/time/distance³ */
      xp->sf_data.SFR = c_star * physical_density *
                        sqrt(physical_density * 32.f * G) / sqrt(3 * M_PI);
    } else {
      xp->sf_data.SFR = 0.;
    }
  }
  return;
}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    struct part* p, struct xpart* xp, const struct star_formation* starform,
    const struct engine* e, const double dt_star) {
  /* Calculate the propability of forming a star, see Revaz & Jablonka (2012),
   * eq. (4) */
  const double prob = 1. - exp(xp->sf_data.SFR * dt_star * (-1.) /
                               hydro_get_physical_density(p, e->cosmology));
  /* Get a unique random number between 0 and 1 for star formation */
  const double random_number =
      random_unit_interval(p->id, e->ti_current, random_number_star_formation);
  if (random_number > prob) {
    return 0;
  } else {
    return 1;
  }
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct star_formation* starform, const int with_cosmology) {
  /* Check if it is the first time steps after star formation */
  if (xp->sf_data.SFR > 0.) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      xp->sf_data.SFR = -(e->cosmology->a);
    } else {
      xp->sf_data.SFR = -(e->time);
    }
  }
}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle.
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 * @param starform the star formation law properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
INLINE static void star_formation_copy_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* Store the current mass */
  sp->mass = hydro_get_mass(p);

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }
  /* Store the chemistry struct in the star particle */
  sp->chemistry_data = p->chemistry_data;

  /* Store the tracers data */
  sp->tracers_data = xp->tracers_data;

  /* Store the birth density in the star particle */
  sp->birth_density = hydro_get_physical_density(p, cosmo);
  /* Store the birth temperature*/
  sp->birth_temperature =
      get_temperature(starform->phys_const, starform->hydro_props, starform->us,
                      cosmo, e->cooling_func, p, xp);
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param starform the star formation law properties to initialize
 *
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    struct star_formation* starform) {
  /* Jeans Number used in starformation elligibility criterion */
  starform->Njeans =
      parser_get_param_int(parameter_file, "GEARStarFormation:jeans_number");
  /* Starformation efficiency (used in SFR calculation) */
  starform->star_formation_rate = parser_get_param_double(
      parameter_file, "GEARStarFormation:star_formation_efficiency");
  /* Wheter we take into account local turbulence */
  starform->with_sigma = parser_get_param_int(
      parameter_file, "GEARStarFormation:with_turbulence_estimation");
  /* Maximum temperature for starformation */
  starform->Max_temperature =
      parser_get_param_double(parameter_file,
                              "GEARStarFormation:Max_temperature") *
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  ;
  starform->hydro_props = hydro_props;
  starform->us = us;
  starform->phys_const = phys_const;
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {
  message("Star formation law is 'GEAR'");
}

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd The global star_formation information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_formation_end_density(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* cd, const struct cosmology* cosmo) {
  /* To finish the turbulence estimation we devide by the density */
  xp->sf_data.sigma2 /= hydro_get_physical_density(p, cosmo);
}
/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(struct part* restrict p,
                                      struct xpart* restrict xp,
                                      const struct star_formation* cd,
                                      const struct cosmology* cosmo) {
  /* If part has 0 neighbours, the estimation of turbulence is 0 */
  xp->sf_data.sigma2 = 0.f;
}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_part(const struct phys_const* restrict phys_const,
                               const struct unit_system* restrict us,
                               const struct cosmology* restrict cosmo,
                               const struct star_formation* data,
                               struct part* restrict p) {

  xp->sf_data.sigma2 = 0.f;
}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to extended particle data
 * @param data The global star_formation information.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* data) {
  xp->sf_data.sigma2 = 0.f;
}

#endif /* SWIFT_GEAR_STAR_FORMATION_H */
