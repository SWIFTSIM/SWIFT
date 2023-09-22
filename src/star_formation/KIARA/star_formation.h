/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
 *               2023 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_KIARA_STAR_FORMATION_H
#define SWIFT_KIARA_STAR_FORMATION_H

/* Local includes */
#include "adiabatic_index.h"
#include "chemistry.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "exp10.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "random.h"
#include "stars.h"
#include "units.h"

#define star_formation_need_update_dx_max 0

/**
 * @file src/star_formation/KIARA/star_formation.h
 * @brief Star formation model used in the KIARA model (currently same as KIARA)
 */

/**
 * @brief Functional form of the star formation law
 */
enum star_formation_H2_model {
  simba_star_formation_density_thresh, /*<! All eligible gas is fully star-forming (H2_frac=1)*/
  simba_star_formation_kmt_model, /*<! Use Krumholz+Gnedin 2011 subgrid model for H2 */
  simba_star_formation_grackle_model /*<! Use H2_frac computed by grackle (or 1-HI_frac) */
};

/**
 * @brief Properties of the KIARA star formation model.
 */
struct star_formation {

  /* The Schmidt law parameters */
  struct {

    /*! Star formation efficiency */
    double sfe;

    /*! star formation efficiency over free fall time constant (cm^1.5 g^-.5
     * s^-1) */
    double mdot_const;

  } schmidt_law;

  /*! Critical overdensity above which SF is allowed */
  double min_over_den;

  struct {

    /*! (Subgrid) temperature threshold for SF to use combined with the density
     * threshold */
    double T_threshold;

    /*! (Subgrid) Hydrogen number density threshold for SF */
    double nH_threshold;

  } subgrid_thresh;

  /* H2 model --------------------------------------------------------- */

  /*! Which H2 model are we using? */
  enum star_formation_H2_model H2_model;

  /*! Scaling factor for the clumping factor in the KMT H2 model */
  float clumping_factor_scaling;

  /* Convert g/cm^2 to Msun/pc^2 */
  double surface_rho_to_Msun_per_parsec2;

  /* Convert code_mass/code_area to g/cm^2 */
  double conv_factor_surface_rho_to_cgs;

  /* Total metal mass fraction in the Sun */
  float Z_solar;

};


/**
 * @brief Calculate if the satisfies the conditions for star formation.
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
INLINE static int star_formation_is_star_forming_subgrid(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* entropy_floor_props) {

  const double number_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Get the Hydrogen mass fraction */
  const double XH = chemistry_get_metal_mass_fraction_for_star_formation(
      p)[chemistry_element_H];

  /* Get the subgrid properties
   * Note these are both in physical frame already */
  const double subgrid_T_cgs = cooling_get_subgrid_temperature(p, xp);
  const double subgrid_rho = cooling_get_subgrid_density(p, xp);
  const double subgrid_n_H = subgrid_rho * XH / phys_const->const_proton_mass;
  const double subgrid_n_H_cgs = subgrid_n_H * number_density_to_cgs;

  /* Now, determine whether we are very cold or (cold and dense) enough
   *
   * This would typically be (T < 10^3 OR (T < 10^4.5 AND n_H > 10))
   * with T and n_H subgrid properties.
   *
   * Recall that particles above the EoS have T_sub = T and rho_sub = rho.
   */
  return ( (subgrid_T_cgs < starform->subgrid_thresh.T_threshold) &&
           (subgrid_n_H_cgs > starform->subgrid_thresh.nH_threshold) );
}

/**
 * @brief Calculate if the gas particle satisfies the conditions for star formation.
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

  /* Decide whether we should form stars or not */

  /* No star formation for particles in the wind */
  if (p->feedback_data.decoupling_delay_time > 0.f) return 0;

  /* No star formation for particles that can't cool */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f) return 0;

  /* Minimal density (converted from mean baryonic density)
   * for star formation */
  const double rho_mean_b_times_min_over_den =
      cosmo->mean_density_Omega_b * starform->min_over_den;

  /* Physical density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Check overdensity criterion */
  if (physical_density < rho_mean_b_times_min_over_den) return 0;

  //message("SF %d\n",star_formation_is_star_forming_subgrid( p, xp, starform, phys_const, cosmo, hydro_props, us, cooling, entropy_floor_props));

  return star_formation_is_star_forming_subgrid(
    p, xp, starform, phys_const, cosmo, hydro_props,
    us, cooling, entropy_floor_props);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart. The star formation is calculated as a simple
 * Schmidt law with an efficiency per free-fall time that can be specified,
 * the free-fall time is based on the total SPH density.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR_schmidt_law(
    struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

  /* Mass density of this particle */
  const float physical_density = cooling_get_subgrid_density(p, xp);

  /* Calculate the SFR per gas mass */
  const double SFRpergasmass =
      starform->schmidt_law.mdot_const * sqrt(physical_density);

  /* Store the SFR */
  p->sf_data.SFR = p->sf_data.H2_fraction * SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart.
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
    struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

  /* Abort early if time-step size is 0 */
  if (dt_star == 0.) {
    p->sf_data.SFR = 0.f;
    return;
  }

  /* Hydrogen number density of this particle (assuming primordial H abundance)
   */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  if (starform->H2_model == simba_star_formation_density_thresh) {
    p->sf_data.H2_fraction = 1.f;
  }
  else if (starform->H2_model == simba_star_formation_kmt_model) {
    /* gas_sigma is double because we do some cgs conversions */
    double gas_sigma = 0.f;
    float gas_Z = 0.f;
    float chi = 0.f;
    float s = 0.f;
    float clumping_factor = 30.f;
    float gas_gradrho_mag = 0.f;

    p->sf_data.H2_fraction = 0.f;

    gas_Z = p->chemistry_data.metal_mass_fraction_total;
    gas_Z /= starform->Z_solar;
    if (gas_Z < 0.01f) {
      gas_Z = 0.01f;
    }

    if (physical_density > 0.f) {
      gas_gradrho_mag = sqrtf(
        p->rho_gradient[0] * p->rho_gradient[0] +
        p->rho_gradient[1] * p->rho_gradient[1] +
        p->rho_gradient[2] * p->rho_gradient[2]
      );

      if (gas_gradrho_mag > 0) {
        gas_sigma = (p->rho * p->rho) / gas_gradrho_mag;

        /* surface density must be in Msun/pc^2 */
        gas_sigma *= starform->surface_rho_to_Msun_per_parsec2 
                      * cosmo->a2_inv;

        /* Lower clumping factor with higher resolution 
          (CF = 30 @ ~1 kpc resolution) */
        clumping_factor *= starform->clumping_factor_scaling;
        if (clumping_factor < 1.f) {
          clumping_factor = 1.f;
        }

        /* chi ~ 1/R ~ 1/clump from KG11 eq. 3 */
        chi = 0.756f * (1.f + 3.1f * powf(gas_Z, 0.365f)) 
                * (30.f / clumping_factor);
        s = logf(1.f + 0.6f * chi + 0.01f * chi * chi) 
              / (0.0396f * powf(clumping_factor, 2.f / 3.f) 
                  * gas_Z * gas_sigma);

        if (s > 0.f && s < 2.f) {
          p->sf_data.H2_fraction = 1.f - 0.75f * (s / (1.f + 0.25f * s));
        }
      }
    }
  }
  else if (starform->H2_model == simba_star_formation_grackle_model) {
#if COOLING_GRACKLE_MODE >= 2
    p->sf_data.H2_fraction = (xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac);
#else
    p->sf_data.H2_fraction = (1. - xp->cooling_data.HI_frac);
#endif
  }
  else {
      error("Invalid H2 model in star formation!!!");
  }

  star_formation_compute_SFR_schmidt_law(p, xp, starform, phys_const,
                                             hydro_props, cosmo, dt_star);
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
 * @param star_prob The probability of converting to a star particle.
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct engine* e,
    const double dt_star,
    double *star_prob) {

  /* Calculate the propability of forming a star */
  const double prob = p->sf_data.SFR * dt_star / hydro_get_mass(p);
  *star_prob = prob;

  /* Get a unique random number between 0 and 1 for star formation */
  const double random_number =
      random_unit_interval(p->id, e->ti_current, random_number_star_formation);

  /* Have we been lucky and need to form a star? */
  return (prob > random_number);
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
  if (p->sf_data.SFR > 0.f) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      p->sf_data.SFR = -e->cosmology->a;
    } else {
      p->sf_data.SFR = -e->time;
    }
  }
}

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
 * @param convert_part Did we convert a part (or spawned one)?
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

  /* Store the current mass as the initial mass */
  sp->mass_init = hydro_get_mass(p);

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Move over the splitting data */
  sp->split_data = xp->split_data;

  /* Store the chemistry struct in the star particle */
  sp->chemistry_data = p->chemistry_data;

  /* Store the tracers data */
  sp->tracers_data = xp->tracers_data;

  /* Store the birth density in the star particle */
  sp->birth_density = cooling_get_subgrid_density(p, xp);

  /* Store the birth temperature in the star particle */
  sp->birth_temperature = cooling_get_subgrid_temperature(p, xp);

  /* Flag that this particle has not done feedback yet */
  sp->feedback_data.feedback_mass_to_launch = 0.f;
  sp->feedback_data.feedback_energy_reservoir = 0.f;
  sp->last_enrichment_time = sp->birth_time;
  sp->count_since_last_enrichment = -1;
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

  /* Get the Gravitational constant */
  const double G_newton = phys_const->const_newton_G;

  /* Get the surface density unit Msun / pc^2 in internal units */
  const double Msun_per_pc2 =
      phys_const->const_solar_mass /
      (phys_const->const_parsec * phys_const->const_parsec);

  starform->surface_rho_to_Msun_per_parsec2 = 1. / Msun_per_pc2;
  starform->conv_factor_surface_rho_to_cgs = 
            units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) *
              units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);

  /* Read the H2 model we are using */
  char H2_model[32];
  parser_get_param_string(parameter_file,
                              "KIARAStarFormation:H2_model", H2_model);

  if (strstr(H2_model, "Thresh") != NULL) {
    starform->H2_model = simba_star_formation_density_thresh;
  }
  else if (strstr(H2_model, "KMT") != NULL) {
    starform->H2_model = simba_star_formation_kmt_model;
  }
  else if (strstr(H2_model, "Grackle") != NULL) {
    starform->H2_model = simba_star_formation_grackle_model;
  }
  else {
    error("Invalid H2 model in SF params %s", H2_model);
  }

  /* Read the ISM subgrid clumping factor value at the resolved scale (KMT model only) */
  starform->clumping_factor_scaling =
        parser_get_opt_param_double(parameter_file,
                           "KIARAStarFormation:clumping_factor_scaling", 30.f);

  /* Read the total metal mass fraction of the Sun */
  starform->Z_solar = 
        parser_get_opt_param_double(parameter_file,
                           "KIARAStarFormation:Z_solar", 0.0134f);

  /* Get the star formation efficiency */
  starform->schmidt_law.sfe = parser_get_param_double(
        parameter_file, "KIARAStarFormation:star_formation_efficiency");

  /* Calculate the ff constant */
  const double ff_const = sqrt(3.0 * M_PI / (32.0 * G_newton));

  /* Calculate the constant */
  starform->schmidt_law.mdot_const = starform->schmidt_law.sfe / ff_const;

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den = parser_get_param_double(
      parameter_file, "KIARAStarFormation:min_over_density");

  /* Read threshold properties */
  starform->subgrid_thresh.T_threshold = parser_get_param_double(
  parameter_file, "KIARAStarFormation:threshold_temperature_K");

  starform->subgrid_thresh.nH_threshold = parser_get_param_double(
  parameter_file, "KIARAStarFormation:threshold_number_density_H_p_cm3");

}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  message("Star formation model is KIARA");

  message("Star formation law is a Schmidt law: Star formation efficiency = %e",
          starform->schmidt_law.sfe);
  message("Particles are star-forming if their properties obey "
          "T < %e K AND n_H > %e cm^-3, and have overdensity > %e",
  starform->subgrid_thresh.T_threshold,
  starform->subgrid_thresh.nH_threshold,
  starform->min_over_den);

}

/**
 * @brief Return the star formation rate of a particle.
 * Remember that SFR can be <0 because it stores last expansion factor
 * when it was SF.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
INLINE static float star_formation_get_SFR(const struct part* p,
                                           const struct xpart* xp) {
  if (p->sf_data.SFR <= 0.)
    return 0.f;
  else
    return p->sf_data.SFR;
}

/**
 * @brief Finishes the density calculation.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the KIARA star formation model.
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
 * density loop for the KIARA star formation model.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(struct part* p, struct xpart* xp,
                                      const struct star_formation* cd,
                                      const struct cosmology* cosmo) {
}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * state to start the density loop.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the KIARA star formation model.
 *
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    struct part* p, const struct star_formation* data) {
}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state at the beginning of the simulation after the ICs have been read.
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
                               struct part* p, struct xpart* xp) {
  /* This may need to be updated elsewhere */
  p->sf_data.H2_fraction = 0.f;
  star_formation_init_part(p, data);
}

/**
 * @brief Split the star formation content of a particle into n pieces
 *
 * We only need to split the SFR if it is positive, i.e. it is not
 * storing the redshift/time of last SF event.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void star_formation_split_part(
    struct part* p, struct xpart* xp, const double n) {

  if (p->sf_data.SFR > 0.) p->sf_data.SFR /= n;
}

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
  /* Nothing to do, we just skip it and deal with it next step */
}

/**
 * @brief Compute some information for the star formation model based
 * on all the particles that were read in.
 *
 * This is called once on start-up of the code.
 *
 * Nothing to do here for KIARA.
 *
 * @param star_form The #star_formation structure.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_stats(struct star_formation* star_form,
                                const struct engine* e) {}

#endif /* SWIFT_KIARA_STAR_FORMATION_H */
