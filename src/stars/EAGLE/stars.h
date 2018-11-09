/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_STARS_H
#define SWIFT_EAGLE_STARS_H

#include <float.h>
#include "minmax.h"
#include "imf.h"

static const float SNII_min_mass_msun = 100; // temporary guess.
static const float imf_max_mass_msun = 100; // temporary guess.
static const float log_min_metallicity = 0.0; // temporary guess.

/**
 * @brief Computes the gravity time-step of a given star particle.
 *
 * @param sp Pointer to the s-particle data.
 */
__attribute__((always_inline)) INLINE static float stars_compute_timestep(
    const struct spart* const sp) {

  return FLT_MAX;
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp) {

  sp->time_bin = 0;
}

/**
 * @brief Prepares a s-particle for its interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_init_spart(
    struct spart* sp) {

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    sp->ids_ngbs_density[i] = -1;
  sp->num_ngb_density = 0;
#endif

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The particle.
 */
__attribute__((always_inline)) INLINE static void stars_reset_predicted_values(
    struct spart* restrict sp) {}

/**
 * @brief Finishes the calculation of (non-gravity) forces acting on stars
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_end_force(
    struct spart* sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void stars_kick_extra(
    struct spart* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_density(
    struct spart* sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->density.wcount *= h_inv_dim;
  sp->density.wcount_dh *= h_inv_dim_plus_one;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_spart_has_no_neighbours(
    struct spart* restrict sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  sp->density.wcount = kernel_root * h_inv_dim;
  sp->density.wcount_dh = 0.f;
}

// -------------------- Work in progress ------------------------------

// This really needs to be somewhere else
inline float interpol_1d(float *table, int i, float dx) {
  float result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

inline float interpol_2d(float **table, int i, int j, float dx, float dy) {
  float result;

  result = (1 - dx) * (1 - dy) * table[i][j] + (1 - dx) * dy * table[i][j + 1] +
           dx * (1 - dy) * table[i + 1][j] + dx * dy * table[i + 1][j + 1];

  return result;
}



inline float dying_mass_msun(float age_Gyr, float metallicity,
			     const struct stars_props *restrict star_properties) {

  // check out units for all these quantities and name them accordingly
  float mass = 0, d_metal, d_time1 = 0, d_time2 = 0, log_age_yr, mass1, mass2;

  int metal_index, i, index_time1 = -1, index_time2 = -1;

  // What do we do with the constants here?
  switch(star_properties->stellar_lifetime_flag) {
    case 0:
      /* padovani & matteucci 1993 */

      if (age_Gyr > 0.039765318659064693)
        mass = exp(M_LN10 * 
            (7.764 - (1.79 -
	             (1.338 - 0.1116 * (9 + log10(age_Gyr))) * 
                     (1.338 - 0.1116 * (9 + log10(age_Gyr)))) /
                        0.2232));
      else if (age_Gyr > 0.003)
        mass = pow((age_Gyr - 0.003) / 1.2, -1 / 1.85);
      else
        mass = imf_max_mass_msun;
      break;
     
    case 1:
      /* maeder & meynet 1989 */

      if (age_Gyr >= 8.4097378)
        mass = exp(M_LN10 * (1 - log10(age_Gyr)) / 0.6545);
      else if (age_Gyr >= 0.35207776)
        mass = exp(M_LN10 * (1.35 - log10(age_Gyr)) / 3.7);
      else if (age_Gyr >= 0.050931493)
        mass = exp(M_LN10 * (0.77 - log10(age_Gyr)) / 2.51);
      else if (age_Gyr >= 0.010529099)
        mass = exp(M_LN10 * (0.17 - log10(age_Gyr)) / 1.78);
      else if (age_Gyr >= 0.0037734787)
        mass = exp(M_LN10 * (-0.94 - log10(age_Gyr)) / 0.86);
      else if (age_Gyr > 0.003)
        mass = pow((age_Gyr - 0.003) / 1.2, -0.54054053);
      else
        mass = imf_max_mass_msun;
      break;

    case 2:
      /* portinari et al. 1998 */

      if (age_Gyr <= 0) {
        mass = imf_max_mass_msun;
        break;
      }

      log_age_yr = log10(age_Gyr * 1.e9);

      // Can we simplify this whole thing somehow to be more readable?
      if (metallicity <= star_properties->lifetimes.metallicity[0]) {
        metal_index = 0;
        d_metal = 0.0;
      } else if (metallicity >= star_properties->lifetimes.metallicity[star_properties->lifetimes.n_z - 1]) {
        metal_index = star_properties->lifetimes.n_z - 2;
        d_metal = 1.0;
      } else {
        for (metal_index = 0; metal_index < star_properties->lifetimes.n_z - 1; metal_index++)
          if (star_properties->lifetimes.metallicity[metal_index + 1] > metallicity) break;

        d_metal = (metallicity - star_properties->lifetimes.metallicity[metal_index]) /
                  (star_properties->lifetimes.metallicity[metal_index + 1] -
                   star_properties->lifetimes.metallicity[metal_index]);
      }

      if (log_age_yr >= star_properties->lifetimes.dyingtime[metal_index][0]) {
        index_time1 = 0;
        d_time1 = 0.0;
      } else if (log_age_yr <= star_properties->lifetimes.dyingtime[metal_index][star_properties->lifetimes.n_mass - 1]) {
        index_time1 = star_properties->lifetimes.n_mass - 2;
        d_time1 = 1.0;
      }

      if (log_age_yr >= star_properties->lifetimes.dyingtime[metal_index + 1][0]) {
        index_time2 = 0;
        d_time2 = 0.0;
      } else if (log_age_yr <=
                 star_properties->lifetimes.dyingtime[metal_index + 1][star_properties->lifetimes.n_mass - 1]) {
        index_time2 = star_properties->lifetimes.n_mass - 2;
        d_time2 = 1.0;
      }

      i = star_properties->lifetimes.n_mass;
      while (i >= 0 && (index_time1 == -1 || index_time2 == -1)) {
        i--;
        if (star_properties->lifetimes.dyingtime[metal_index][i] >= log_age_yr && index_time1 == -1) {
          index_time1 = i;
          d_time1 = (log_age_yr - star_properties->lifetimes.dyingtime[metal_index][index_time1]) /
                    (star_properties->lifetimes.dyingtime[metal_index][index_time1 + 1] -
                     star_properties->lifetimes.dyingtime[metal_index][index_time1]);
        }
        if (star_properties->lifetimes.dyingtime[metal_index + 1][i] >= log_age_yr && index_time2 == -1) {
          index_time2 = i;
          d_time2 = (log_age_yr - star_properties->lifetimes.dyingtime[metal_index + 1][index_time2]) /
                    (star_properties->lifetimes.dyingtime[metal_index + 1][index_time2 + 1] -
                     star_properties->lifetimes.dyingtime[metal_index + 1][index_time2]);
        }
      }

      mass1 = interpol_1d(star_properties->lifetimes.mass, index_time1, d_time1);
      mass2 = interpol_1d(star_properties->lifetimes.mass, index_time2, d_time2);

      mass = (1.0 - d_metal) * mass1 + d_metal * mass2;
      break;

    default:
      error("stellar lifetimes not defined\n");

  }
  if (mass > imf_max_mass_msun) mass = imf_max_mass_msun;

  return mass;
}

inline float lifetime_in_Gyr(float mass, float metallicity,
			     const struct stars_props *restrict star_properties) {

  double time = 0, d_mass, d_metal;

  int mass_index, metal_index;

  switch (star_properties->stellar_lifetime_flag) {
    case 0:
      /* PM93 (Padovani & Matteucci 1993) */

      if (mass <= 0.6)
        time = 160.0;
      else if (mass <= 6.6)
        time =
            pow(10.0, (0.334 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) /
                          0.1116);
      else
        time = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 1:
      /* MM89 (Maeder & Meynet 1989) */

      if (mass <= 1.3)
        time = pow(10.0, -0.6545 * log10(mass) + 1.0);
      else if (mass <= 3.0)
        time = pow(10.0, -3.7 * log10(mass) + 1.35);
      else if (mass <= 7.0)
        time = pow(10.0, -2.51 * log10(mass) + 0.77);
      else if (mass <= 15.0)
        time = pow(10.0, -1.78 * log10(mass) + 0.17);
      else if (mass <= 60.0)
        time = pow(10.0, -0.86 * log10(mass) - 0.94);
      else
        time = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 2:
      /* P98 (Portinari et al. 1998) */

      if (mass <= star_properties->lifetimes.mass[0]) {
        mass_index = 0;
        d_mass = 0.0;
      } else if (mass >= star_properties->lifetimes.mass[star_properties->lifetimes.n_mass - 1]) {
        mass_index = star_properties->lifetimes.n_mass - 2;
        d_mass = 1.0;
      } else {
        for (mass_index = 0; mass_index < star_properties->lifetimes.n_mass - 1; mass_index++)
          if (star_properties->lifetimes.mass[mass_index + 1] > mass) break;

        d_mass = (mass - star_properties->lifetimes.mass[mass_index]) /
                 (star_properties->lifetimes.mass[mass_index + 1] - star_properties->lifetimes.mass[mass_index]);
      }

      if (metallicity <= star_properties->lifetimes.metallicity[0]) {
        metal_index = 0;
        d_metal = 0.0;
      } else if (metallicity >= star_properties->lifetimes.metallicity[star_properties->lifetimes.n_z - 1]) {
        metal_index = star_properties->lifetimes.n_z - 2;
        d_metal = 1.0;
      } else {
        for (metal_index = 0; metal_index < star_properties->lifetimes.n_z - 1; metal_index++)
          if (star_properties->lifetimes.metallicity[metal_index + 1] > metallicity) break;

        d_metal = (metallicity - star_properties->lifetimes.metallicity[metal_index]) /
                  (star_properties->lifetimes.metallicity[metal_index + 1] -
                   star_properties->lifetimes.metallicity[metal_index]);
      }

      /* time in Gyr */
      time = exp(M_LN10 * interpol_2d(star_properties->lifetimes.dyingtime, metal_index, mass_index,
                                       d_metal, d_mass)) / 1.0e9;

      break;

    default:
      error("stellar lifetimes not defined");
  }

  return time;
}

inline void determine_bin_yield(int *iz_low, int *iz_high, float *dz, float log_metallicity, 
				const struct stars_props *restrict star_properties){
  if (log_metallicity > log_min_metallicity) {
    int j;
    for (j = 0; j < star_properties->yield_AGB.N_Z - 1 &&
                         log_metallicity > star_properties->yield_AGB.Metallicity[j + 1];
         j++);
    *iz_low = j;
    *iz_high = *iz_low + 1;

    if (*iz_high >= star_properties->yield_AGB.N_Z) *iz_high = star_properties->yield_AGB.N_Z - 1;

    if (log_metallicity >= star_properties->yield_AGB.Metallicity[0] &&
        log_metallicity <= star_properties->yield_AGB.Metallicity[star_properties->yield_AGB.N_Z - 1])
      *dz = log_metallicity - star_properties->yield_AGB.Metallicity[*iz_low];
    else
      *dz = 0;

    float deltaz = star_properties->yield_AGB.Metallicity[*iz_high] - star_properties->yield_AGB.Metallicity[*iz_low];

    if (deltaz > 0)
      *dz = *dz / deltaz;
    else
      dz = 0;
  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0;
  }
}


inline static void evolve_SNIa(void){

}

inline static void evolve_SNII(void){

}

inline static void evolve_AGB(float log_min_mass, float log_max_mass,
                              float log_metallicity, float *initial_metals,
			      const struct stars_props *restrict stars,
			      struct spart *restrict sp){
  // come up with more descriptive index names
  int ilow, ihigh, imass, i = 0;

  //if (AGB_MassTransferOn == 0) return;

  float metallicity = exp(M_LN10 * log_metallicity);
  message("metallicity %.5e", metallicity); // TEMPORARY to help compiling

  /* determine integration range, limiting to stars that become AGB stars */
  if (log_max_mass > log10(SNII_min_mass_msun))
    log_max_mass = log10(SNII_min_mass_msun);

  if (log_min_mass >= log_max_mass) return;

  /* determine which mass bins will contribute */
  determine_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh, stars);
  float *stellar_yield;
  stellar_yield = (float *)malloc(ihigh-ilow*sizeof(float));
  if (stellar_yield == NULL) error("could not allocate stellar_yield");

  /* determine yield of these bins (not equally spaced bins) */
  int iz_low, iz_high;
  float dz;
  determine_bin_yield(&iz_low, &iz_high, &dz, log_metallicity, stars);

  /* compute stellar_yield as function of mass */
  float metals[chemistry_element_count], mass;
  for (i = 0; i < chemistry_element_count; i++) {
    for (imass = ilow; imass < ihigh + 1; imass++)
      /* yield_AGB.SPH refers to elements produced, yield_AGB.Ejecta_SPH refers
       * to elements already in star */
      stellar_yield[imass] =
          (1 - dz) * (stars->yield_AGB.SPH[iz_low][i][imass] +
                      initial_metals[i] * stars->yield_AGB.Ejecta_SPH[iz_low][imass]) +
          dz * (stars->yield_AGB.SPH[iz_high][i][imass] +
                initial_metals[i] * stars->yield_AGB.Ejecta_SPH[iz_high][imass]);

// temporarily commented ifdef until sorted out how to treat this.
// #ifdef BG_DOUBLE_IMF
//     if (phys_dens > All.IMF_PhysDensThresh)
//       metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 6, stellar_yield);
//     else
//       metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 2, stellar_yield);
// #else
    metals[i] = integrate_imf(log_min_mass, log_max_mass, 0.0, 2, stellar_yield,stars);
// #endif
  }

  for (imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] =
        (1 - dz) * (stars->yield_AGB.total_metals_SPH[iz_low][imass] +
                    metallicity * stars->yield_AGB.Ejecta_SPH[iz_low][imass]) +
        dz * (stars->yield_AGB.total_metals_SPH[iz_high][imass] +
              metallicity * stars->yield_AGB.Ejecta_SPH[iz_high][imass]);

// temporarily commented ifdef until sorted out how to treat this.
// #ifdef BG_DOUBLE_IMF
//   if (phys_dens > All.IMF_PhysDensThresh)
//     mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 6, stellar_yield);
//   else
//     mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2, stellar_yield);
// #else
  mass = integrate_imf(log_min_mass, log_max_mass, 0.0, 2, stellar_yield, stars);
// #endif

  /* yield normalization */
  float norm0, norm1;

  /* zero all negative values */
  for (i = 0; i < chemistry_element_count; i++)
    if (metals[i] < 0) metals[i] = 0;

  if (mass < 0) mass = 0;

  /* get the total mass ejected from the table */
  for (imass = ilow; imass < ihigh + 1; imass++)
    stellar_yield[imass] = (1 - dz) * stars->yield_AGB.Ejecta_SPH[iz_low][imass] +
                           dz * stars->yield_AGB.Ejecta_SPH[iz_high][imass];

  norm0 = integrate_imf(log_min_mass, log_max_mass, 0.0, 2, stellar_yield, stars);

  /* compute the total mass ejected */
  norm1 = mass + metals[chemistry_element_H] + metals[chemistry_element_He];

  /* normalize the yields */
  if (norm1 > 0) {
    for (i = 0; i < chemistry_element_count; i++) metals[i] *= (norm0 / norm1);
    mass *= (norm0 / norm1);

    // Is this really the right way around? mass += metals and metals += mass? Maybe rename variables?
    for (i = 0; i < chemistry_element_count; i++) {
      sp->metals_released[i] += metals[i];
      sp->mass_from_agb += metals[i];
    }
    sp->metal_mass_released += mass;
    sp->metals_from_agb += mass;
  } else {
    error("wrong normalization!!!! norm1 = %e\n", norm1);
  }

  free(stellar_yield);
}

// Analogue of eagle_do_stellar_evolution...
inline static void compute_stellar_evolution(const struct stars_props *restrict star_properties, 
					     struct spart *restrict sp) {
  
  // temporary declarations. need to be passed in.
  float age_of_star_Gyr = 0;
  float dt_Gyr = 0;
  float metallicity = 1.0;

  // set max and min mass of dying stars
  float log10_max_dying_mass =
      log10(dying_mass_msun(age_of_star_Gyr, metallicity, star_properties));
  float log10_min_dying_mass = log10(
      dying_mass_msun(age_of_star_Gyr + dt_Gyr, metallicity, star_properties));

  if (log10_min_dying_mass > log10_max_dying_mass) error("min dying mass is greater than max dying mass");

  /* integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_msun */
  if (log10_min_dying_mass == log10_max_dying_mass) return;

  float log_metallicity;
  if (metallicity > 0)
    log_metallicity = log10(metallicity);
  else
    log_metallicity = log_min_metallicity;


  // dummy declaration, needs to be done properly...
  float initial_metals[10];
  for(int j = 0; j < 10; j++) initial_metals[j] = 0.0;

  // Evolve SNIa, SNII, AGB

  evolve_SNIa();

  evolve_SNII(); 

  evolve_AGB(0,0,0,initial_metals,star_properties,sp);

}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 * @param stars_properties The #stars_props
 */
__attribute__((always_inline)) INLINE static void stars_evolve_spart(
    struct spart* restrict sp, const struct stars_props* stars_properties,
    const struct cosmology* cosmo) {

    // Sort out what happens with first_call static

    // get mass and initial mass of particle to pass to evolution function
    // except we're working in mass fraction so maybe not necessary...
    //float mass = hydro_get_mass(p);
    //float initial_mass = ;

    // Set elements released to zero

    // Evolve the star
    compute_stellar_evolution(stars_properties, sp);
    
    // set_particle_metal_content
}

//-------------------- Star formation work in progress --------------------

//inline static void compute_sfr(struct part* restrict p, 
//	const struct cosmology* restrict cosmo,
//	const struct unit_system* us, float dt) {
//  
//  const dt_cgs = dt*units_cgs_conversion_factor(us, UNIT_CONV_TIME);
//  
//  int N_stars_generated;
//  double sfr = 0;
//
//  // get flags (lines 196-197 in eagle_sfr.c)
//
//  // Find out how many stars generated
//  N_stars_generated = 0; // set to zero for now
//
//  // Mass of stars (CHECK what's meant to be happening here...) 
//  double particle_mass_cgs = hydro_get_mass(p) * units_cgs_conversion_factor(us,UNIT_CONV_MASS);
//  double generations = 1.0; // copied from parameter files in eagle. find out what it means...
//  double mass_of_star = particle_mass_cgs/(generations - N_stars_generated);
//
//  // get pressure
//  double pressure_cgs = hydro_get_physical_pressure(p, cosmo) * units_cgs_conversion_factor(us,UNIT_CONV_PRESSURE);
//
//  // get density
//  double density_cgs = hydro_get_physical_density(p,cosmo) * units_cgs_conversion_factor(us,UNIT_CONV_DENSITY);
//
//  // set ks parameters. RHS need to be read in (TODO).
//  //double ks_norm_cgs = KS_norm_cgs;
//  //double ks_exponent = SF_SchmidtLawExponent;
//  //if (density_cgs > KSHighDens_thres_cgs) {
//  //  ks_norm_cgs = KS_norm_HighDens_cgs;
//  //  ks_exponent = All.SF_SchmidtLawHighDensExponent;
//  //}
//
//  // inverse star formation timescale in seconds
//  double timescale_sfr_inv_cgs = ks_norm_cgs * pow(pressure_cgs * hydro_gamma / const_newton_G_cgs, 0.5*(ks_exponent - 1)); // CHECK if possible values of ks_exponent can allow for not using pow
//
//  // calculate star formation rate, and stellar mass formed
//  double sm_dt_cgs = particle_mass_cgs * inv_tsfr_cgs;
//  double sm_cgs = sm_dt_cgs * dt_cgs;
//
//  // if density is high enough convert all gas mass to stars right away (CHECK: there is condition in EAGLE to bypass this check)
//  if (particle_mass_cgs > max_physical_density) {
//    sm_cgs = particle_mass_cgs/(generations - N_stars_generated);
//    sm_dt_cgs = sm_cgs/dt_cgs;
//  }
//
//  // EAGLE has a check for if sm_cgs > particle_mass_cgs. do we need this?
//
//  // EAGLE adds the stellar mass to a running total. This would be over all the particles, so maybe need to do this higher up. Also, is it necessary?
//  // total_sm += sm;
//  
//  // convert to formation rate to stellar mass per year
//  sfr = sm_dt_cgs * sfr_conversion_factor;
//  // add sfr to appropriate time bin?
//  //TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
//
//  // calculate probability of turning part into spart
//  // does this need an extra switch as in EAGLE?
//  float part_to_spart_prob;
//  if (density_cgs > max_density) {
//    part_to_spart_prob = 1;
//  } else {
//    part_to_spart_prob = sm / mass_of_star;
//  }
//
//  // Get random number with gsl (see Eagle/system.c:155) 
//  // random_n = ...
//  
//  if (star_formation == 0) random_n = part_to_spart_prob + 1;
//
//  if (random_n < part_to_spart_prob) ; // Make star...
//
//}

#endif /* SWIFT_EAGLE_STARS_H */
