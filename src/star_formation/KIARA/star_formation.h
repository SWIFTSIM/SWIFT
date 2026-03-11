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
#include "fof.h"
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
 * @brief Star formation model used in the KIARA model
 */

/**
 * @brief Functional form of the star formation law and H2 model
 */
enum star_formation_H2_model {
  /*<! All eligible gas is fully star-forming (H2_frac=1) */
  kiara_star_formation_density_thresh,
  /*<! Use Krumholz+Gnedin 2011 subgrid model for H2 */
  kiara_star_formation_kmt_model,
  /*<! Use H2_frac computed by grackle (or 1-HI_frac) */
  kiara_star_formation_grackle_model
};

enum star_formation_SF_model {
  /*<! Use Kennicutt-Schmidt law based on free-fall time with const efficiency
   */
  kiara_star_formation_SchmidtLaw,
  /*<! Use model based on Wada & Norman 2007 to compute SFR */
  kiara_star_formation_WadaNorman,
  /*<! Use Schmidt Law with efficiency based on lognormal threshold density */
  kiara_star_formation_lognormal
};

/**
 * @brief Properties of the KIARA star formation model.
 */
struct star_formation {

  /*! Which SF model are we using? */
  enum star_formation_SF_model SF_model;

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

  /*! Convert g/cm^2 to Msun/pc^2 */
  double surface_rho_to_Msun_per_parsec2;

  /*! Convert code_mass/code_area to g/cm^2 */
  double conv_factor_surface_rho_to_cgs;

  /*! Total metal mass fraction in the Sun */
  float Z_solar;

  /* The lognormal & WadaNorman SF model parameters */
  struct {
    /*! Characteristic density (approx 10^-1.5 Mo/pc^3) */
    double rho0;

    /*! Critical number density for SF (input parameter) */
    double ncrit;

    /*! Critical density in code units for SF */
    double rhocrit;

    /*! for WN07, time unit conversion factor to scale efficiency */
    double time_to_year_inverse;

    /*! for WN07, M* unit conversion factor to scale efficiency */
    double to_solar_mass;

    /*! for WN07, SFR unit conversion factor to scale efficiency */
    double to_msun_per_yr;

    /*! The epsilon value to control sub-grid factors, multiplies SFR */
    double epsilon;

    /*! For lognormal, conversion factor for free fall time */
    double ff_const_inv;

    /*! For WN07, flag to set method used to compute epsc */
    int wn07_epsc_method;

  } lognormal;
};

/**
 * @brief Determine the mass loading factor for
 *        a given star particle based on its host galaxy.
 *
 * @param group_stellar_mass The stellar mass of the host galaxy
 * @param minimum_galaxy_stellar_mass Floor for stellar mass in eta calculation
 * @param FIRE_eta_normalization Normalization of eta at FIRE_eta_break
 * @param FIRE_eta_break M* at which eta(M*) slope changes
 * @param FIRE_eta_lower_slope Slope below break
 * @param FIRE_eta_upper_slope Slope above break
 */
__attribute__((always_inline)) INLINE static float feedback_mass_loading_factor(
    const struct cosmology *cosmo, const double group_stellar_mass,
    const float minimum_galaxy_stellar_mass, const float FIRE_eta_normalization,
    const float FIRE_eta_break, const float FIRE_eta_lower_slope,
    const float FIRE_eta_upper_slope, const float FIRE_eta_lower_slope_EOR,
    const float wind_eta_suppression_redshift) {

  const float m_star = fmax(group_stellar_mass, minimum_galaxy_stellar_mass);

  float slope = FIRE_eta_lower_slope;
  if (cosmo->z > 6) slope = FIRE_eta_lower_slope_EOR;
  if (m_star > FIRE_eta_break) {
    slope = FIRE_eta_upper_slope;
  }

  float eta = FIRE_eta_normalization * powf(m_star / FIRE_eta_break, slope);

  const float a_suppress_inv = (1.f + fabs(wind_eta_suppression_redshift));
  if (wind_eta_suppression_redshift > 0 &&
      cosmo->z > wind_eta_suppression_redshift) {
    eta *= cosmo->a * cosmo->a * a_suppress_inv * a_suppress_inv;
  } else if (wind_eta_suppression_redshift < 0) {
    eta *= expf(-powf(cosmo->a * a_suppress_inv, -3.f));
  }

  return fmax(eta, 0.f);
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
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * @param entropy_floor_props The entropy floor assumed in this run.
 */
INLINE static int star_formation_is_star_forming_subgrid(
    const struct part *p, const struct xpart *xp,
    const struct star_formation *starform, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const struct entropy_floor_properties *entropy_floor_props) {

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
  return ((subgrid_T_cgs < starform->subgrid_thresh.T_threshold) &&
          (subgrid_n_H_cgs > starform->subgrid_thresh.nH_threshold));
}

/**
 * @brief Calculate if the gas particle satisfies the conditions for star
 *        formation.
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
    const struct part *p, const struct xpart *xp,
    const struct star_formation *starform, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const struct entropy_floor_properties *entropy_floor_props) {

  /* Decide whether we should form stars or not */

  /* No star formation for particles in the wind */
  if (p->decoupled) return 0;

  /* No star formation when outside subgrid model */
  if (cooling_get_subgrid_temperature(p, xp) <= 0.f) return 0;

  /* Minimal density (converted from mean baryonic density)
   * for star formation */
  const float rho_mean_b_times_min_over_den =
      cosmo->mean_density_Omega_b * starform->min_over_den;

  /* Physical density of the particle */
  const float physical_density = hydro_get_physical_density(p, cosmo);

  /* Check overdensity criterion */
  if (physical_density < rho_mean_b_times_min_over_den) return 0;

  return star_formation_is_star_forming_subgrid(p, xp, starform, phys_const,
                                                cosmo, hydro_props, us, cooling,
                                                entropy_floor_props);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #part. The star formation is calculated as a simple
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
    struct part *p, struct xpart *xp, const struct star_formation *starform,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct cosmology *cosmo, const double dt_star) {

  /* Mass density of this particle */
  // const float physical_density = cooling_get_subgrid_density(p, xp);
  const float physical_density = cooling_get_subgrid_density(p, xp);

  /* Calculate the SFR per gas mass */
  const float SFRpergasmass =
      starform->schmidt_law.mdot_const * sqrt(physical_density);

  /* Store the SFR */
  p->sf_data.SFR = p->sf_data.H2_fraction * SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #part. The star formation is calculated based on
 * the model of Wada & Norman 2007, eq. 17
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR_wn07(
    struct part *p, struct xpart *xp, const struct star_formation *starform,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct cosmology *cosmo, const double dt_star) {

  p->sf_data.SFR = 0.f;
  /* Mean density of the gas described by a lognormal PDF (i.e. dense gas). */
  const double rho_V = cooling_get_subgrid_density(p, xp);
  const double rho_0 = starform->lognormal.rho0;

  /* Collect information about galaxy that the particle belongs to */
  const float galaxy_mstar = p->galaxy_data.stellar_mass;
  const float galaxy_ssfr = p->galaxy_data.specific_sfr;
  const float galaxy_sfr = galaxy_mstar * galaxy_ssfr;

  /* Density is too low, so no SF */
  if (rho_V <= 1.001 * rho_0) return;

  /* Calculate SFR efficiency, which WN07 suggests should
   * scale from 0.1 for starbursts to 0.001 for normal galaxies
   */
  double epsc = 0.01;
  /* if it's not in a galaxy, assume it's a small proto-galaxy so starburst */
  if (galaxy_mstar == 0.f) {
    epsc = 0.1;
  } else if (starform->lognormal.wn07_epsc_method == 0) { /* constant value */
    epsc = 0.01;
  }
  /* deviation from all-galaxies main sequence taken from Koprowski+24 */
  else if (starform->lognormal.wn07_epsc_method == 1 ||
           starform->lognormal.wn07_epsc_method == 2) {
    const double mstar = galaxy_mstar * starform->lognormal.to_solar_mass;
    const double sfrmax_data = pow(10., 3.69 - 3.81 * exp(-0.47 * cosmo->z));
    const double M0_data = pow(10., 11.91 - 2.48 * exp(-0.44 * cosmo->z));
    const double sfr_data = sfrmax_data / (1. + (M0_data / mstar));
    const double sfr_msun_per_yr =
        galaxy_sfr * starform->lognormal.to_msun_per_yr;
    if (starform->lognormal.wn07_epsc_method == 1) {
      epsc = 0.01 * sfr_msun_per_yr / sfr_data;
    } else {
      epsc = 0.01 * sqrtf(sfr_msun_per_yr / sfr_data);
    }
  }
  /* deviation from SF-galaxies main sequence data taken from Koprowski+24 */
  else if (starform->lognormal.wn07_epsc_method == 3 ||
           starform->lognormal.wn07_epsc_method == 4) {
    const double mstar = galaxy_mstar * starform->lognormal.to_solar_mass;
    const double sfrmax_data = pow(10., 3.47 - 3.13 * exp(-0.56 * cosmo->z));
    const double M0_data = pow(10., 11.69 - 1.66 * exp(-0.53 * cosmo->z));
    const double sfr_data = sfrmax_data / (1. + (M0_data / mstar));
    const double sfr_msun_per_yr =
        galaxy_sfr * starform->lognormal.to_msun_per_yr;
    epsc = 0.01 * sqrt(sfr_msun_per_yr / sfr_data);
  }
  /* based on direct scaling with sSFR */
  else if (starform->lognormal.wn07_epsc_method == 5) {
    epsc =
        min(galaxy_ssfr * starform->lognormal.time_to_year_inverse * 1.e7, 0.1);
  }
  /* Scale with galaxy SFR instead of SSFR */
  else if (starform->lognormal.wn07_epsc_method == 6) {
    const float sfr_msun_per_yr =
        galaxy_sfr * starform->lognormal.to_msun_per_yr;
    epsc = min(0.001 * cbrt(sfr_msun_per_yr * 300.), 0.1);
    epsc = max(epsc, 0.001);
  } else {
    error("wn07_epsc_method value of %d not recognised",
          starform->lognormal.wn07_epsc_method);
  }

  /* Restrict range to that specified in WN07 */
  epsc = min(epsc, 0.1);
  epsc = max(epsc, 0.001);

  /* Calculate parameters in WN07 model */
  const double sigma = sqrt(2. * log(rho_V / rho_0));
  const double z_num = log(starform->lognormal.rhocrit / rho_0) - sigma * sigma;
  const double z_den = sqrt(2.) * sigma;
  const double z = z_num / z_den;

  /* fraction of dense gas */
  const double fc = 0.5 * erfc(z);

  /* This is the SFR density from eq. 17, except use actual
   * density rho_V not estimated density rho_c. 3pi/32=0.294524.
   */
  const double rhosfr = epsc * sqrt(0.294524 * phys_const->const_newton_G * rho_V) * fc;

  /* multiply by dense gas effective volume to get SFR (rho_V appears in both
   * this eqn and previous one so it is cancelled out for efficiency) */
  const double sfr =
      rhosfr * (p->cooling_data.subgrid_fcold * hydro_get_mass(p));

  /* Multiply by the H2 fraction */
  p->sf_data.SFR = starform->lognormal.epsilon * sfr * p->sf_data.H2_fraction;
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #part. The star formation is calculated by assuming
 * a lognormal density distribution with a mean density given by the
 * subgrid density, above a critical density threshold that is an
 * input parameter.  Lognormal params based on sims by Wada & Norman 2007.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR_lognormal(
    struct part *p, struct xpart *xp, const struct star_formation *starform,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct cosmology *cosmo, const double dt_star) {

  /* Limit for Eq. 12 in Wada+Norman 2007 */
  const double rho_limit = 1.001 * starform->lognormal.rho0;

  /* H2 fraction in the particle */
  const double f_H2 = p->sf_data.H2_fraction;

  /* Mean density of the gas described by a lognormal PDF
   * (i.e. subgrid ISM gas) PHYSICAL */
  double rho_V = cooling_get_subgrid_density(p, xp);

  /* SF cannot occur below characteristic
   * density~1 cm^-3 (formulae below give nans)
   */
  if (rho_V < rho_limit || f_H2 <= 0.) {
    p->sf_data.SFR = 0.f;
    return;
  }

  /* Calculate the SFR based on a lognormal density distribution at rho0 with a
   * threshold density for star formation rhocrit.  sigma comes from WN07 model.
   */
  const double rho_0 = starform->lognormal.rho0;
  /* Mass-averaged density for cold phase */
  const double sigma = sqrt(log(2. * rho_V / rho_0));
  const double z_num =
      log(starform->lognormal.rhocrit / rho_0) - (sigma * sigma);
  const double z_den = sqrt(2.) * sigma;
  const double z = z_num / z_den;

  /* Calculate lognormal fraction from the WN07 model */
  const double f_c = 0.5 * erfc(z);

  const double rho_phys = hydro_get_physical_density(p, cosmo);

  /* Calculate the SFR per gas mass, using lognormal mass fraction above
   * rhocrit as efficiency
   */
  const double sSFR = f_c * starform->lognormal.ff_const_inv * sqrt(rho_phys);

  const double mass = f_H2 * p->cooling_data.subgrid_fcold * hydro_get_mass(p);

  /* Store the SFR */
  p->sf_data.SFR = starform->lognormal.epsilon * sSFR * mass;
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #part.  This involves computing H2 fraction, then using
 * the chosen SF model to compute the SFR.
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
    struct part *p, struct xpart *xp, const struct star_formation *starform,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct cosmology *cosmo, const double dt_star) {

  /* Abort early if time-step size is 0 */
  if (dt_star == 0.) {
    p->sf_data.SFR = 0.f;
    return;
  }

  /* Physical gas density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Compute the H2 fraction of the particle */
  switch (starform->H2_model) {
    case kiara_star_formation_density_thresh:
      p->sf_data.H2_fraction = 1.f;
      break;
    case kiara_star_formation_kmt_model:
      p->sf_data.H2_fraction = 0.f;

      /* gas_sigma is double because we do some cgs conversions */
      double gas_sigma = 0.f;
      float gas_Z = 0.f;
      float chi = 0.f;
      float s = 0.f;
      float clumping_factor = 30.f;
      float gas_gradrho_mag = 0.f;

      gas_Z = p->chemistry_data.metal_mass_fraction_total;
      gas_Z /= starform->Z_solar;
      if (gas_Z < 0.01f) {
        gas_Z = 0.01f;
      }

      if (physical_density > 0.f) {
        gas_gradrho_mag = sqrtf(p->rho_gradient[0] * p->rho_gradient[0] +
                                p->rho_gradient[1] * p->rho_gradient[1] +
                                p->rho_gradient[2] * p->rho_gradient[2]);

        if (gas_gradrho_mag > 0) {
          gas_sigma = (p->rho * p->rho) / gas_gradrho_mag;

          /* surface density must be in Msun/pc^2 */
          gas_sigma *=
              starform->surface_rho_to_Msun_per_parsec2 * cosmo->a2_inv;

          /* Lower clumping factor with higher resolution
            (CF = 30 @ ~1 kpc resolution) */
          clumping_factor *= starform->clumping_factor_scaling;
          if (clumping_factor < 1.f) {
            clumping_factor = 1.f;
          }

          /* chi ~ 1/R ~ 1/clump from KG11 eq. 3 */
          chi = 0.756f * (1.f + 3.1f * powf(gas_Z, 0.365f)) *
                (30.f / clumping_factor);
          s = logf(1.f + 0.6f * chi + 0.01f * chi * chi) /
              (0.0396f * powf(clumping_factor, 2.f / 3.f) * gas_Z * gas_sigma);

          if (s > 0.f && s < 2.f) {
            p->sf_data.H2_fraction = 1.f - 0.75f * (s / (1.f + 0.25f * s));
          }
        }
      }
      break;
    case kiara_star_formation_grackle_model:
#if COOLING_GRACKLE_MODE >= 2
      p->sf_data.H2_fraction =
          p->cooling_data.subgrid_fcold *
          (xp->cooling_data.H2I_frac + xp->cooling_data.H2II_frac);
#else
      p->sf_data.H2_fraction = 1. - xp->cooling_data.HI_frac;
#endif
      break;
    default:
      error("Invalid H2 model in star formation!");
      break;
  }

  /* Now compute the star formation rate and save it to the particle */
  switch (starform->SF_model) {
    case kiara_star_formation_SchmidtLaw:
      star_formation_compute_SFR_schmidt_law(p, xp, starform, phys_const,
                                             hydro_props, cosmo, dt_star);
      break;
    case kiara_star_formation_WadaNorman:
      star_formation_compute_SFR_wn07(p, xp, starform, phys_const, hydro_props,
                                      cosmo, dt_star);
      break;
    case kiara_star_formation_lognormal:
      star_formation_compute_SFR_lognormal(p, xp, starform, phys_const,
                                           hydro_props, cosmo, dt_star);
      break;
    default:
      error("Invalid SF model in star formation!!!");
      break;
  }
}

/**
 * @brief Returns the number of new star particles to create per SF event.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 *
 * @return The number of extra star particles to generate per gas particles.
 *        (return 0 if the gas particle itself is to be converted)
 */
INLINE static int star_formation_number_spart_to_spawn(
    struct part *p, struct xpart *xp, const struct star_formation *starform) {

  return 0;
}

/**
 * @brief Returns the number of particles to convert per SF event.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 *
 * @return The number of particles to generate per gas particles.
 *        (This has to be 0 or 1)
 */
INLINE static int star_formation_number_spart_to_convert(
    const struct part *p, const struct xpart *xp,
    const struct star_formation *starform) {

  return 1;
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
    const struct part *p, const struct xpart *xp,
    const struct star_formation *starform, const struct engine *e,
    const double dt_star) {

  /* Calculate the propability of forming a star */
  const double prob = max(p->sf_data.SFR, 0.f) * dt_star / hydro_get_mass(p);

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
    struct part *p, struct xpart *xp, const struct star_formation *starform) {
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
    struct part *p, struct xpart *xp, const struct engine *e,
    const struct star_formation *starform, const int with_cosmology) {

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
 * @param chem_data The global properties of the chemistry scheme.
 * @param convert_part Did we convert a part (or spawned one)?
 * @param (return) displacement The 3D displacement vector of the star with
 * respect to the sink position.
 */
INLINE static void star_formation_copy_properties(
    const struct part *p, const struct xpart *xp, struct spart *sp,
    const struct engine *e, const struct star_formation *starform,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling,
    const struct chemistry_global_data *chem_data, const int convert_part,
    float displacement[3]) {

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

  /* Flag that this particle has not done feedback yet
  sp->feedback_data.physical_energy_reservoir = 0.;
  sp->feedback_data.N_launched = 0;
  sp->feedback_data.mass_to_launch = 0.f;
  sp->feedback_data.total_mass_kicked = 0.f;
  sp->feedback_data.wind_velocity = 0.f;
  sp->feedback_data.eta_suppression_factor = 1.f;*/
  sp->last_enrichment_time = sp->birth_time;
  sp->count_since_last_enrichment = 0;

  /* Copy FoF galaxy data from spawning particle */
  sp->galaxy_data.stellar_mass = p->galaxy_data.stellar_mass;
  sp->galaxy_data.gas_mass = p->galaxy_data.gas_mass;
  sp->galaxy_data.specific_sfr = p->galaxy_data.specific_sfr;

  /* Slightly displace spawned particle to avoid zeros */
  if (1) {

    const float max_displacement = 0.1;
    const double delta_x =
        2.f * random_unit_interval(sp->id, e->ti_current,
                                   (enum random_number_type)0) -
        1.f;
    const double delta_y =
        2.f * random_unit_interval(sp->id, e->ti_current,
                                   (enum random_number_type)1) -
        1.f;
    const double delta_z =
        2.f * random_unit_interval(sp->id, e->ti_current,
                                   (enum random_number_type)2) -
        1.f;

    /* Update the displacement */
    displacement[0] = delta_x * max_displacement * p->h;
    displacement[1] = delta_y * max_displacement * p->h;
    displacement[2] = delta_z * max_displacement * p->h;

    /* Move the spart */
    sp->x[0] += displacement[0];
    sp->x[1] += displacement[1];
    sp->x[2] += displacement[2];

    /* Copy the position to the gpart */
    sp->gpart->x[0] = sp->x[0];
    sp->gpart->x[1] = sp->x[1];
    sp->gpart->x[2] = sp->x[2];
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
    struct swift_params *parameter_file, const struct phys_const *phys_const,
    const struct unit_system *us, const struct hydro_props *hydro_props,
    const struct cosmology *cosmo,
    const struct entropy_floor_properties *entropy_floor,
    struct star_formation *starform) {

  /* Get the Gravitational constant */
  const double G_newton = phys_const->const_newton_G;

  /* CGS density conversion */
  const double rho_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  /* Get the surface density unit Msun / pc^2 in internal units */
  const double Msun_per_pc2 =
      phys_const->const_solar_mass /
      (phys_const->const_parsec * phys_const->const_parsec);

  starform->surface_rho_to_Msun_per_parsec2 = 1. / Msun_per_pc2;
  starform->conv_factor_surface_rho_to_cgs =
      rho_to_cgs * units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);

  /* Read the SF model we are using */
  char SF_model[32];
  parser_get_param_string(parameter_file, "KIARAStarFormation:SF_model",
                          SF_model);

  if (strstr(SF_model, "SchmidtLaw") != NULL) {
    starform->SF_model = kiara_star_formation_SchmidtLaw;
  } else if (strstr(SF_model, "WadaNorman") != NULL) {
    starform->SF_model = kiara_star_formation_WadaNorman;
  } else if (strstr(SF_model, "lognormal") != NULL) {
    starform->SF_model = kiara_star_formation_lognormal;
  } else {
    error("Invalid SF model in SF params %s", SF_model);
  }

  /* Read the H2 model we are using */
  char H2_model[32];
  parser_get_param_string(parameter_file, "KIARAStarFormation:H2_model",
                          H2_model);

  if (strstr(H2_model, "Thresh") != NULL) {
    starform->H2_model = kiara_star_formation_density_thresh;
  } else if (strstr(H2_model, "KMT") != NULL) {
    starform->H2_model = kiara_star_formation_kmt_model;
  } else if (strstr(H2_model, "Grackle") != NULL) {
    starform->H2_model = kiara_star_formation_grackle_model;
  } else {
    error("Invalid H2 model in SF params %s", H2_model);
  }

  /* Read the ISM subgrid clumping factor value at the resolved scale
   * (KMT model only)
   */
  starform->clumping_factor_scaling = parser_get_opt_param_double(
      parameter_file, "KIARAStarFormation:clumping_factor_scaling", 30.f);

  /* Read the total metal mass fraction of the Sun */
  starform->Z_solar = parser_get_opt_param_double(
      parameter_file, "KIARAStarFormation:Z_solar", 0.0134f);

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den = parser_get_param_double(
      parameter_file, "KIARAStarFormation:min_over_density");

  /* Read threshold properties */
  starform->subgrid_thresh.T_threshold = parser_get_param_double(
      parameter_file, "KIARAStarFormation:threshold_temperature_K");

  starform->subgrid_thresh.nH_threshold = parser_get_param_double(
      parameter_file, "KIARAStarFormation:threshold_number_density_H_p_cm3");

  /* Calculate the ff constant */
  const double ff_const = sqrt(3. * M_PI / (32. * G_newton));

  if (starform->SF_model == kiara_star_formation_SchmidtLaw) {
    /* Get the star formation efficiency */
    starform->schmidt_law.sfe = parser_get_param_double(
        parameter_file, "KIARAStarFormation:star_formation_efficiency");

    /* Calculate the constant */
    starform->schmidt_law.mdot_const = starform->schmidt_law.sfe / ff_const;
  } else if (starform->SF_model == kiara_star_formation_WadaNorman ||
             starform->SF_model == kiara_star_formation_lognormal) {

    /* Critical density for SF (in physical cm^-3) for lognormal SF model */
    starform->lognormal.ncrit = parser_get_opt_param_double(
        parameter_file, "KIARAStarFormation:lognormal_critical_density", 1.e3);

    /* code units */
    starform->lognormal.rhocrit =
        starform->lognormal.ncrit * 1.673e-24 / rho_to_cgs;

    /* Set characeristic density of ISM lognormal (Table 1 of Wada+Norman07)
     * in code units */
    starform->lognormal.rho0 =
        pow(10., -1.5) * 1.98841e33 / (rho_to_cgs * pow(3.08567758149e18, 3.));

    /* Like the star formation efficiency but to control for sub-grid factors */
    starform->lognormal.epsilon = parser_get_param_double(
        parameter_file, "KIARAStarFormation:lognormal_epsilon");

    /* used to scale epsilon_c (efficiency) in lognormal model to sSFR */
    starform->lognormal.time_to_year_inverse =
        (365.25 * 24. * 60. * 60.) /
        units_cgs_conversion_factor(us, UNIT_CONV_TIME);

    /* used to scale epsilon_c (efficiency) in lognormal model to SFR */
    starform->lognormal.to_solar_mass =
        units_cgs_conversion_factor(us, UNIT_CONV_MASS) / 1.98841e33;

    /* used to scale epsilon_c (efficiency) in lognormal model to SFR */
    starform->lognormal.to_msun_per_yr =
        units_cgs_conversion_factor(us, UNIT_CONV_SFR) / 1.98841e33 *
        (365.25 * 24. * 60. * 60.);

    /* Calculate the constant needed for the free-fall time */
    starform->lognormal.ff_const_inv = 1. / ff_const;

    if (starform->SF_model == kiara_star_formation_WadaNorman) {
      starform->lognormal.wn07_epsc_method = parser_get_opt_param_int(
          parameter_file, "KIARAStarFormation:wn07_epsc_method", 0);
    }
  }
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation *starform) {

  message("Star formation model is KIARA");

  message(
      "Particles are eligible for star formation if their have "
      "T < %e K, n_H > %e cm^-3, and overdensity > %e",
      starform->subgrid_thresh.T_threshold,
      starform->subgrid_thresh.nH_threshold, starform->min_over_den);

  if (starform->SF_model == kiara_star_formation_SchmidtLaw) {
    message(
        "Star formation law is a Schmidt law: Star formation "
        "efficiency = %e",
        starform->schmidt_law.sfe);
  } else if (starform->SF_model == kiara_star_formation_WadaNorman) {
    message(
        "Star formation based on Wada+Norman 2007: critical "
        "density (code units) = %e",
        starform->lognormal.rhocrit);
  } else if (starform->SF_model == kiara_star_formation_lognormal) {
    message(
        "Star formation based on lognormal density pdf: "
        "critical density (code units) = %e"
        "efficiency = %e",
        starform->lognormal.rhocrit, starform->lognormal.epsilon);
  } else {
    error("Invalid SF model in star formation!!!");
  }
}

/**
 * @brief Return the star formation rate of a particle.
 * Remember that SFR can be <0 because it stores last expansion factor
 * when it was SF.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
INLINE static float star_formation_get_SFR(const struct part *p,
                                           const struct xpart *xp) {
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
    struct part *p, struct xpart *xp, const struct star_formation *cd,
    const struct cosmology *cosmo) {}

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
star_formation_part_has_no_neighbours(struct part *p, struct xpart *xp,
                                      const struct star_formation *cd,
                                      const struct cosmology *cosmo) {}

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
    struct part *p, const struct star_formation *data) {}

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
star_formation_first_init_part(const struct phys_const *phys_const,
                               const struct unit_system *us,
                               const struct cosmology *cosmo,
                               const struct star_formation *data,
                               struct part *p, struct xpart *xp) {
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
    struct part *p, struct xpart *xp, const double n) {

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
star_formation_no_spart_available(const struct engine *e, const struct part *p,
                                  const struct xpart *xp) {
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
star_formation_first_init_stats(struct star_formation *star_form,
                                const struct engine *e) {}

#endif /* SWIFT_KIARA_STAR_FORMATION_H */
