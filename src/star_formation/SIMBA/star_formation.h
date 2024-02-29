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
#ifndef SWIFT_SIMBA_STAR_FORMATION_H
#define SWIFT_SIMBA_STAR_FORMATION_H

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
 * @file src/star_formation/SIMBA/star_formation.h
 * @brief Star formation model used in the SIMBA model (currently same as SIMBA)
 */

/**
 * @brief Functional form of the star formation law
 */
enum star_formation_law {
  simba_star_formation_schmidt_law, /*<! Schmidt law */
  simba_star_formation_pressure_law /*<! Pressure law */
};

/**
 * @brief Choice of star formation threshold
 */
enum star_formation_threshold {
  simba_star_formation_threshold_Z_dep,  /*<! SF threshold based on Schaye+2004
                                            Z-dependence */
  simba_star_formation_threshold_subgrid /*<! SF threshold based on subgrid
                                            properties */
};

/**
 * @brief Functional form of the star formation law
 */
enum star_formation_H2_model {
  simba_star_formation_density_thresh, /*<! All eligible gas is fully star-forming (H2_frac=1)*/
  simba_star_formation_kmt_model, /*<! Use Krumholz+Gnedin 2011 subgrid model for H2 */
  simba_star_formation_grackle_model /*<! Use H2_frac computed by grackle (or 1-HI_frac) */
};

/**
 * @brief Properties of the SIMBA star formation model.
 */
struct star_formation {

  /* SF law ------------------------------------------------------------*/

  /*! Which SF law are we using? */
  enum star_formation_law SF_law;

  /* The Schmidt model parameters */
  struct {

    /*! Star formation efficiency */
    double sfe;

    /*! star formation efficiency over free fall time constant (cm^1.5 g^-.5
     * s^-1) */
    double mdot_const;

  } schmidt_law;

  /* The pressure law model parameters */
  struct {

    /*! Normalization of the KS star formation law (internal units) */
    double KS_normalization;

    /*! Normalization of the KS star formation law (Msun / kpc^2 / yr) */
    double KS_normalization_MSUNpYRpKPC2;

    /*! Slope of the KS law */
    double KS_power_law;

    /*! Slope of the high density KS law */
    double KS_high_den_power_law;

    /*! KS law High density (H number density) threshold (internal units) */
    double KS_high_den_thresh;

    /*! KS high density normalization (internal units) */
    double KS_high_den_normalization;

    /*! KS high density normalization (H atoms per cm^3)  */
    double KS_high_den_thresh_HpCM3;

    /*! gas fraction */
    double fgas;

    /*! Star formation law slope */
    double SF_power_law;

    /*! star formation normalization (internal units) */
    double SF_normalization;

    /*! star formation high density slope */
    double SF_high_den_power_law;

    /*! Star formation high density normalization (internal units) */
    double SF_high_den_normalization;

  } pressure_law;

  /* Density for direct conversion to star -------------------------------- */

  /*! Max physical density (H atoms per cm^3)*/
  double gas_density_direct_HpCM3;

  /*! Max physical density (internal units) */
  double gas_density_direct;

  /* SF threshold --------------------------------------------------------- */

  enum star_formation_threshold SF_threshold;

  /*! Critical overdensity above which SF is allowed */
  double min_over_den;

  struct {

    /*! Density threshold to form stars (internal units) */
    double density_threshold;

    /*! Density threshold to form stars in user units */
    double density_threshold_HpCM3;

    /*! Maximum density threshold to form stars (internal units) */
    double density_threshold_max;

    /*! Maximum density threshold to form stars (H atoms per cm^3) */
    double density_threshold_max_HpCM3;

    /*! Reference metallicity for metal-dependant threshold */
    double Z0;

    /*! Inverse of reference metallicity */
    double Z0_inv;

    /*! critical density Metallicity power law (internal units) */
    double n_Z0;

    /*! Dalla Vecchia & Schaye entropy differnce criterion */
    double entropy_margin_threshold_dex;

    /*! 10^Tdex of Dalla Vecchia & Schaye entropy difference criterion */
    double ten_to_entropy_margin_threshold_dex;

  } Z_dep_thresh;

  struct {

    /*! (Subgrid) temperature threshold for SF to use on its own */
    double T_threshold1;

    /*! (Subgrid) temperature threshold for SF to use combined with the density
     * threshold */
    double T_threshold2;

    /*! (Subgrid) Hydrogen number density threshold for SF */
    double nH_threshold;

    /*! (Subgrid) fraction of ISM gas assumed to be in cold (SF) phase */
    double cold_ISM_frac;

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
 * @param entropy_floor_props The entropy floor assumed in this run.
 */
INLINE static int star_formation_is_star_forming_Z_dep(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* entropy_floor_props) {

  /* Physical density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Get the Hydrogen number density (assuming primordial H abundance) */
  const double rho_H = physical_density * hydro_props->hydrogen_mass_fraction;

  /* Get the density threshold for star formation */
  const double Z =
      chemistry_get_total_metal_mass_fraction_for_star_formation(p);

  double density_threshold;

  /* Schaye (2004), eq. 19 and 24 */
  if (Z > 0.) {
    density_threshold =
        starform->Z_dep_thresh.density_threshold *
        powf(Z * starform->Z_dep_thresh.Z0_inv, starform->Z_dep_thresh.n_Z0);
    density_threshold =
        min(density_threshold, starform->Z_dep_thresh.density_threshold_max);
  } else {
    density_threshold = starform->Z_dep_thresh.density_threshold_max;
  }

  /* Convert to mass density */
  density_threshold *= phys_const->const_proton_mass;

  /* Check if it exceeded the minimum density */
  if (rho_H < density_threshold) return 0;

  /* Calculate the entropy of the particle */
  const double entropy = hydro_get_physical_entropy(p, xp, cosmo);

  /* Calculate the entropy that will be used to calculate
   * the off-set from the EoS */
  const double entropy_eos = entropy_floor(p, cosmo, entropy_floor_props);

  /* Check the Schaye & Dalla Vecchia 2012 EOS-based temperature criterion */
  return (entropy <
          entropy_eos *
              starform->Z_dep_thresh.ten_to_entropy_margin_threshold_dex);
}

/**
 * @brief Calculate if the satisfies the conditions for star formation.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param entropy_floor_props The entropy floor assumed in this run.
 */
INLINE static int star_formation_is_star_forming_subgrid(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* entropy_floor_props) {

  /* Get the Hydrogen mass fraction */
  const double XH = chemistry_get_metal_mass_fraction_for_star_formation(
      p)[chemistry_element_H];

  /* Get the subgrid properties
   * Note these are both in physical frame already */
  const double subgrid_T_cgs = cooling_get_subgrid_temperature(p, xp);
  const double subgrid_rho = cooling_get_subgrid_density(p, xp);
  const double subgrid_n_H = subgrid_rho * XH / phys_const->const_proton_mass;

  float T_lim = entropy_floor_temperature(p, cosmo, entropy_floor_props);

  /* Now, determine whether we are very cold or (cold and dense) enough
   *
   * This would typically be (T < 10^3 OR (T < 10^4.5 AND n_H > 10))
   * with T and n_H subgrid properties.
   *
   * Recall that particles above the EoS have T_sub = T and rho_sub = rho.
   */
  return (subgrid_T_cgs < T_lim &&
           subgrid_n_H > starform->subgrid_thresh.nH_threshold);
}

/**
 * @brief Calculate if the satisfies the conditions for star formation.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param entropy_floor_props The entropy floor assumed in this run.
 */
INLINE static int star_formation_is_star_forming(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct unit_system* us, struct cooling_function_data* cooling, 
    const struct entropy_floor_properties* entropy_floor_props) {

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

  /* Decide whether we should form stars or not */

  /* Are we above the threshold for automatic star formation? */
  if (physical_density >
      starform->gas_density_direct * phys_const->const_proton_mass) return 1;

  /* Check overdensity criterion */
  if (physical_density < rho_mean_b_times_min_over_den) return 0;

  /* Check density/temperature/entropy criteria */
  if (starform->SF_threshold == simba_star_formation_threshold_subgrid) {
    return star_formation_is_star_forming_subgrid( p, xp, starform, 
		    phys_const, cosmo, hydro_props, entropy_floor_props);
  } else if (starform->SF_threshold == simba_star_formation_threshold_Z_dep) {
    return star_formation_is_star_forming_Z_dep( p, xp, starform, 
		    phys_const, cosmo, hydro_props, entropy_floor_props);
  } else {
    error("Invalid SF threshold model, aborting.\n");
    return 0;
  }
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
  p->sf_data.SFR = starform->subgrid_thresh.cold_ISM_frac * p->sf_data.H2_fraction * SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart. The star formation is calculated using a pressure
 * law based on Schaye and Dalla Vecchia (2008), the star formation
 * rate is calculated as:
 *
 * \dot{m}_\star = A (1 Msun / pc^-2)^-n m_gas (\gamma/G * f_g *
 *                 pressure)**((n-1)/2)
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR_pressure_law(
    struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

  /* Hydrogen number density of this particle (assuming primordial H abundance)
   */
  const double physical_density = cooling_get_subgrid_density(p, xp);
  const double nH = hydro_props->hydrogen_mass_fraction * physical_density /
                    phys_const->const_proton_mass;

  /* Get the pressure used for the star formation, this is
   * the maximum the physical pressure of the particle and the
   * floor pressure. The floor pressure is used implicitly
   * when getting the physical pressure. */
  const double pressure = hydro_get_physical_pressure(p, cosmo);

  /* Calculate the specific star formation rate */
  double SFRpergasmass;
  if (nH < starform->pressure_law.KS_high_den_thresh) {

    SFRpergasmass = starform->pressure_law.SF_normalization *
                    pow(pressure, starform->pressure_law.SF_power_law);

  } else {

    SFRpergasmass = starform->pressure_law.SF_high_den_normalization *
                    pow(pressure, starform->pressure_law.SF_high_den_power_law);
  }

  /* Store the SFR */
  p->sf_data.SFR = p->sf_data.H2_fraction * SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart. Only called if particle satisfies SF criteria.
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

  /* Compute SFR based on selected method */  
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

    if (hydro_get_physical_density(p, cosmo) > 0.f) {
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
#if defined(COOLING_GRACKLE_MODE)
  else if (starform->H2_model == simba_star_formation_grackle_model) {
#if COOLING_GRACKLE_MODE >= 2
    p->sf_data.H2_fraction = (xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac);
#else
    p->sf_data.H2_fraction = (1. - xp->cooling_data.HI_frac);
#endif
  }
#endif
  else {
      error("Invalid H2 model in star formation (%d)!!!",starform->H2_model);
  }


  /* Determine which star formation model to use */
  switch (starform->SF_law) {

    case simba_star_formation_schmidt_law:
      star_formation_compute_SFR_schmidt_law(p, xp, starform, phys_const,
                                             hydro_props, cosmo, dt_star);
      break;
    case simba_star_formation_pressure_law:
      star_formation_compute_SFR_pressure_law(p, xp, starform, phys_const,
                                              hydro_props, cosmo, dt_star);
      break;
    default:
      error("Invalid star formation model!!!");
  }
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
  sp->birth_density = hydro_get_physical_density(p, cosmo);

  /* Store the birth temperature in the star particle */
  sp->birth_temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                  cosmo, cooling, p, xp);

  /* Flag that this particle has not done feedback yet */
  sp->f_E = -1.f;
  sp->number_of_SNII_events = 0;
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

  /* Get the SF surface density unit Msun / kpc^2 / yr in internal units */
  const double kpc = 1000. * phys_const->const_parsec;
  const double Msun_per_kpc2_per_year =
      phys_const->const_solar_mass / (kpc * kpc) / phys_const->const_year;

  /* Conversion of number density from cgs */
  const double number_density_from_cgs =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Check if we are using the Schmidt law for the star formation rate,
   * defaults to pressure law if is not explicitely set to a Schmidt law */
  char temp[32];
  parser_get_param_string(parameter_file, "SIMBAStarFormation:SF_model", temp);

  /* Read the H2 model we are using */
  char H2_model[32];
  parser_get_param_string(parameter_file,
                              "SIMBAStarFormation:H2_model", H2_model);

  if (strstr(H2_model, "Thresh") != NULL) {
    starform->H2_model = simba_star_formation_density_thresh;
  }
  else if (strstr(H2_model, "KMT") != NULL) {
    starform->H2_model = simba_star_formation_kmt_model;
  }
  else if (strstr(H2_model, "Grackle") != NULL) {
#if defined(COOLING_GRACKLE_MODE)
    starform->H2_model = simba_star_formation_grackle_model;
#else
    error("Cannot use %s H2 model without GRACKLE being on", H2_model);
#endif
  }
  else {
    error("Invalid H2 model in SF params %s", H2_model);
  }

  /* Read the clumping factor scaling, should be in resolution in kpc */
  starform->clumping_factor_scaling =
        parser_get_param_double(parameter_file,
                                  "SIMBAStarFormation:clumping_factor_scaling");

  /* Read the total metal mass fraction of the Sun */
  starform->Z_solar = 
        parser_get_opt_param_double(parameter_file,
                                      "SIMBAStarFormation:Z_solar", 0.0134f);

  if (strcmp(temp, "SchmidtLaw") == 0) {

    /* Schmidt model */
    starform->SF_law = simba_star_formation_schmidt_law;

    /* Get the star formation efficiency */
    starform->schmidt_law.sfe = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:star_formation_efficiency");

    /* Calculate the ff constant */
    const double ff_const = sqrt(3.0 * M_PI / (32.0 * G_newton));

    /* Calculate the constant */
    starform->schmidt_law.mdot_const = starform->schmidt_law.sfe / ff_const;

  } else if (strcmp(temp, "PressureLaw") == 0) {

    /* Pressure model */
    starform->SF_law = simba_star_formation_pressure_law;

    /* Read the gas fraction from the file */
    starform->pressure_law.fgas = parser_get_opt_param_double(
        parameter_file, "SIMBAStarFormation:gas_fraction", 1.);

    /* Read the Kennicutt-Schmidt power law exponent */
    starform->pressure_law.KS_power_law = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:KS_exponent");

    /* Calculate the power law of the corresponding star formation Schmidt law
     */
    starform->pressure_law.SF_power_law =
        (starform->pressure_law.KS_power_law - 1.) / 2.;

    /* Read the normalization of the KS law in KS law units */
    starform->pressure_law.KS_normalization_MSUNpYRpKPC2 =
        parser_get_param_double(parameter_file,
                                "SIMBAStarFormation:KS_normalisation");

    /* Convert to internal units */
    starform->pressure_law.KS_normalization =
        starform->pressure_law.KS_normalization_MSUNpYRpKPC2 *
        Msun_per_kpc2_per_year;

    /* Calculate the starformation pre-factor (eq. 12 of Schaye & Dalla Vecchia
     * 2008) */
    starform->pressure_law.SF_normalization =
        starform->pressure_law.KS_normalization *
        pow(Msun_per_pc2, -starform->pressure_law.KS_power_law) *
        pow(hydro_gamma * starform->pressure_law.fgas / G_newton,
            starform->pressure_law.SF_power_law);

    /* Read the high density Kennicutt-Schmidt power law exponent */
    starform->pressure_law.KS_high_den_power_law = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:KS_high_density_exponent");

    /* Calculate the SF high density power law */
    starform->pressure_law.SF_high_den_power_law =
        (starform->pressure_law.KS_high_den_power_law - 1.) / 2.;

    /* Read the high density criterion for the KS law in number density per cm^3
     */
    starform->pressure_law.KS_high_den_thresh_HpCM3 = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:KS_high_density_threshold_H_p_cm3");

    /* Transform the KS high density criterion to simulation units */
    starform->pressure_law.KS_high_den_thresh =
        starform->pressure_law.KS_high_den_thresh_HpCM3 *
        number_density_from_cgs;

    /* Convert to a mass density assuming primordial abundances */
    const float KS_high_den_thresh_rho =
        starform->pressure_law.KS_high_den_thresh *
        phys_const->const_proton_mass / hydro_props->hydrogen_mass_fraction;

    /* Pressure on the entropy floor at the high-density threshold
     *
     * Note that we use FLT_MAX as the comoving density to make sure
     * the floor is applied no matter what redshift we are at. This will
     * always be a density above the comoving density threashold for the floor
     * to be used.*/
    const double EOS_high_den_pressure = entropy_floor_gas_pressure(
        KS_high_den_thresh_rho, FLT_MAX, cosmo, entropy_floor);

    /* Calculate the KS high density normalization
     * We want the SF law to be continous so the normalisation of the second
     * power-law is the value of the first power-law at the high-density
     * threshold
     */
    starform->pressure_law.KS_high_den_normalization =
        starform->pressure_law.KS_normalization *
        pow(Msun_per_pc2, starform->pressure_law.KS_high_den_power_law -
                              starform->pressure_law.KS_power_law) *
        pow(hydro_gamma * EOS_high_den_pressure * starform->pressure_law.fgas /
                G_newton,
            (starform->pressure_law.KS_power_law -
             starform->pressure_law.KS_high_den_power_law) *
                0.5f);

    /* Calculate the SF high density normalization */
    starform->pressure_law.SF_high_den_normalization =
        starform->pressure_law.KS_high_den_normalization *
        pow(Msun_per_pc2, -starform->pressure_law.KS_high_den_power_law) *
        pow(hydro_gamma * starform->pressure_law.fgas / G_newton,
            starform->pressure_law.SF_high_den_power_law);
  } else {
    error("Invalid SF law model: '%s'", temp);
  }

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den = parser_get_param_double(
      parameter_file, "SIMBAStarFormation:min_over_density");

  /* Get the maximum physical density for SF */
  starform->gas_density_direct_HpCM3 = parser_get_opt_param_double(
      parameter_file, "SIMBAStarFormation:density_direct_H_p_cm3", FLT_MAX);

  /* Convert the maximum physical density to internal units */
  starform->gas_density_direct =
      starform->gas_density_direct_HpCM3 * number_density_from_cgs;

  /* Check if we are using the Schmidt law for the star formation rate,
   * defaults to pressure law if is not explicitely set to a Schmidt law */
  char temp_SF[32];
  parser_get_param_string(parameter_file, "SIMBAStarFormation:SF_threshold",
                          temp_SF);

  if (strcmp(temp_SF, "Zdep") == 0) {

    /* Z-dep (Schaye+2004) model */
    starform->SF_threshold = simba_star_formation_threshold_Z_dep;

    starform->Z_dep_thresh.entropy_margin_threshold_dex =
        parser_get_opt_param_double(parameter_file,
                                    "SIMBAStarFormation:EOS_entropy_margin_dex",
                                    0.5);

    starform->Z_dep_thresh.ten_to_entropy_margin_threshold_dex =
        exp10(starform->Z_dep_thresh.entropy_margin_threshold_dex);

    /* Read the normalization of the metallicity dependent critical
     * density*/
    starform->Z_dep_thresh.density_threshold_HpCM3 = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:threshold_norm_H_p_cm3");

    /* Convert to internal units */
    starform->Z_dep_thresh.density_threshold =
        starform->Z_dep_thresh.density_threshold_HpCM3 *
        number_density_from_cgs;

    /* Read the scale metallicity Z0 */
    starform->Z_dep_thresh.Z0 = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:threshold_Z0");
    starform->Z_dep_thresh.Z0_inv = 1. / starform->Z_dep_thresh.Z0;

    /* Read the power law of the critical density scaling */
    starform->Z_dep_thresh.n_Z0 = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:threshold_slope");

    /* Read the maximum allowed density for star formation */
    starform->Z_dep_thresh.density_threshold_max_HpCM3 =
        parser_get_param_double(
            parameter_file, "SIMBAStarFormation:threshold_max_density_H_p_cm3");

    /* Convert to internal units */
    starform->Z_dep_thresh.density_threshold_max =
        starform->Z_dep_thresh.density_threshold_max_HpCM3 *
        number_density_from_cgs;

  } else if (strcmp(temp_SF, "Subgrid") == 0) {

#ifdef COOLING_SIMBA
    error(
        "The 'Subgrid' SF threshold in the SIMBA star formation model cannot "
        "be used in combination with SIMBA cooling. A cooling model with "
        "subgrid quantities (such as 'COLIBRE' using the Ploeckinger tables) "
        "must be used. Alternatively, the 'Zdep' threshold should be used as "
        "it can be combined with any cooling model.");
#endif

    /* Subgrid quantities based model */
    starform->SF_threshold = simba_star_formation_threshold_subgrid;

    /* Read threshold properties */
    starform->subgrid_thresh.T_threshold1 = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:threshold_temperature1_K");
    starform->subgrid_thresh.T_threshold2 = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:threshold_temperature2_K");
    starform->subgrid_thresh.nH_threshold = parser_get_param_double(
        parameter_file, "SIMBAStarFormation:threshold_number_density_H_p_cm3");
    starform->subgrid_thresh.nH_threshold *= number_density_from_cgs;

    /* When using subgrid model, need the cold ISM fraction */
    starform->subgrid_thresh.cold_ISM_frac = parser_get_opt_param_double(
        parameter_file, "SIMBACooling:cold_ISM_frac", 1.);

  } else {
    error("Invalid SF threshold model: '%s'", temp_SF);
  }
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  message("Star formation model is SIMBA/KIARA");

  switch (starform->SF_threshold) {
    case simba_star_formation_threshold_Z_dep:

      message("Density threshold follows Schaye (2004)");
      message(
          "the normalization of the density threshold is given by"
          " %e #/cm^3, with metallicity slope of %e, and metallicity "
          "normalization"
          " of %e, the maximum density threshold is given by %e #/cm^3",
          starform->Z_dep_thresh.density_threshold_HpCM3,
          starform->Z_dep_thresh.n_Z0, starform->Z_dep_thresh.Z0,
          starform->Z_dep_thresh.density_threshold_max_HpCM3);
      message(
          "Temperature threshold is given by Dalla Vecchia and Schaye (2012)");
      message(
          "The temperature threshold offset from the EOS is given by: %e dex",
          starform->Z_dep_thresh.entropy_margin_threshold_dex);

      break;
    case simba_star_formation_threshold_subgrid:

      message("Density threshold uses subgrid quantities");
      message(
          "Particles are star-forming if their properties obey (T_sub < %e K "
          "OR (T_sub < %e K AND n_H,sub > %e cm^-3))",
          starform->subgrid_thresh.T_threshold1,
          starform->subgrid_thresh.T_threshold2,
          starform->subgrid_thresh.nH_threshold);
      message("Cold ISM fraction = %e", starform->subgrid_thresh.cold_ISM_frac);

      break;
    default:
      error("Invalid star formation threshold!!!");
  }

  switch (starform->SF_law) {
    case simba_star_formation_schmidt_law:
      message(
          "Star formation law is a Schmidt law: Star formation efficiency = %e",
          starform->schmidt_law.sfe);
      break;
    case simba_star_formation_pressure_law:
      message(
          "Star formation law is a pressure law (Schaye & Dalla Vecchia "
          "2008): ");
      message(
          "With properties: normalization = %e Msun/kpc^2/yr, slope of the"
          "Kennicutt-Schmidt law = %e and gas fraction = %e ",
          starform->pressure_law.KS_normalization_MSUNpYRpKPC2,
          starform->pressure_law.KS_power_law, starform->pressure_law.fgas);
      message("At densities of %e H/cm^3 the slope changes to %e.",
              starform->pressure_law.KS_high_den_thresh_HpCM3,
              starform->pressure_law.KS_high_den_power_law);
      break;
    default:
      error("Invalid star formation law!!!");
  }

  message("Running with a direct conversion density of: %e #/cm^3",
          starform->gas_density_direct_HpCM3);
}

/**
 * @brief Return the star formation rate of a particle.
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
 * density loop for the SIMBA star formation model.
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
 * density loop for the SIMBA star formation model.
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
 * density loop for the SIMBA star formation model.
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
 * Nothing to do here.
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
 * Nothing to do here for SIMBA.
 *
 * @param star_form The #star_formation structure.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_stats(struct star_formation* star_form,
                                const struct engine* e) {}

#endif /* SWIFT_SIMBA_STAR_FORMATION_H */
