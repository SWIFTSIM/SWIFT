/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
/**
 * @file src/cooling/KIARA/cooling.c
 * @brief Cooling using the GRACKLE 3.x library.
 */

#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Include header */
#include "cooling.h"

/* Some standard headers. */
#include <assert.h>
#include <fenv.h>
#include <float.h>
#include <math.h>

/* The grackle library itself */
#include <grackle.h>
extern chemistry_data *grackle_data;

/* Local includes. */
#include "chemistry.h"
#include "cooling_io.h"
#include "entropy_floor.h"
#include "error.h"
#include "fof.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "star_formation.h"
#include "units.h"

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * Also calls the additional H reionisation energy injection if need be.
 *
 * @param cosmo The current cosmological model.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The space data, including a pointer to array of particles
 * @param time The current system time
 */
void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time) {

  /* set current time */
  if (cooling->redshift == -1) {
    cooling->units.a_value = cosmo->a;
  } else {
    cooling->units.a_value = 1. / (1. + cooling->redshift);
  }
}

/**
 * @brief Print the chemical network
 *
 * @param xp The #xpart to print
 */
void cooling_print_fractions(const struct xpart *restrict xp) {

  const struct cooling_xpart_data tmp = xp->cooling_data;
#if COOLING_GRACKLE_MODE > 0
  message("HI %g, HII %g, HeI %g, HeII %g, HeIII %g, e %g", tmp.HI_frac,
          tmp.HII_frac, tmp.HeI_frac, tmp.HeII_frac, tmp.HeIII_frac,
          tmp.e_frac);
#endif

#if COOLING_GRACKLE_MODE > 1
  message("HM %g, H2I %g, H2II %g", tmp.HM_frac, tmp.H2I_frac, tmp.H2II_frac);
#endif

#if COOLING_GRACKLE_MODE > 2
  message("DI %g, DII %g, HDI %g", tmp.DI_frac, tmp.DII_frac, tmp.HDI_frac);
#endif
  message("Metal: %g", tmp.metal_frac);
}

/**
 * @brief Initializes grackle particle quantities
 * assuming purely neutral gas
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE void cooling_grackle_init_part(
    const struct cooling_function_data *cooling, struct part *restrict p,
    struct xpart *restrict xp) {

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;
  zero = 0.f;

  /* primordial chemistry >= 1: Start with everything neutral (as in dark ages)
   */
  xp->cooling_data.HI_frac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  xp->cooling_data.HII_frac = zero;
  xp->cooling_data.HeI_frac =
      p->chemistry_data.metal_mass_fraction[chemistry_element_He];
  xp->cooling_data.HeII_frac = zero;
  xp->cooling_data.HeIII_frac = zero;
  xp->cooling_data.e_frac = xp->cooling_data.HII_frac +
                            0.25 * xp->cooling_data.HeII_frac +
                            0.5 * xp->cooling_data.HeIII_frac;
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  xp->cooling_data.HM_frac = zero;
  xp->cooling_data.H2I_frac = zero;
  xp->cooling_data.H2II_frac = zero;
#endif  // MODE >= 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  xp->cooling_data.DI_frac = grackle_data->DeuteriumToHydrogenRatio *
                             grackle_data->HydrogenFractionByMass;
  xp->cooling_data.DII_frac = zero;
  xp->cooling_data.HDI_frac = zero;
#endif  // MODE >= 3
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_first_init_part(const struct phys_const *restrict phys_const,
                             const struct unit_system *restrict us,
                             const struct hydro_props *hydro_props,
                             const struct cosmology *restrict cosmo,
                             const struct cooling_function_data *cooling,
                             struct part *restrict p,
                             struct xpart *restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.time_last_event = -cooling->thermal_time;

  /* Initialize grackle ionization fractions */
  cooling_grackle_init_part(cooling, p, xp);

  p->cooling_data.subgrid_fcold = 0.f;

  /* Initialize dust properties */
#if COOLING_GRACKLE_MODE >= 2
  p->cooling_data.dust_mass = 0.f;
  for (int i = 0; i < chemistry_element_count; i++) {
    p->cooling_data.dust_mass_fraction[i] = 0.f;
  }

  p->cooling_data.dust_temperature = 0.f;
#endif
}

/**
 * @brief Perform additional init on the cooling properties of the
 * (x-)particles that requires the density to be known.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE void cooling_post_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *cooling, const struct part *restrict p,
    struct xpart *restrict xp) {}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
float cooling_get_radiated_energy(const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Computes H and H2 self-shielding for G0 calculation.
 * Based on Schauer et al. 2015 eqs 8,9.
 *
 * @param p The particle to act upon.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static float
cooling_compute_self_shielding(const struct part *restrict p,
                               const struct cooling_function_data *cooling) {

  float fH2_shield = 1.f;
#if COOLING_GRACKLE_MODE >= 2
  float T_ism = p->cooling_data.subgrid_temp;
  if (T_ism > 0.f) {
    /* Compute self-shielding from H */
    const float a = cooling->units.a_value;
    const float a3_inv = 1.f / (a * a * a);
    const double rho_grad_norm2 = p->rho_gradient[0] * p->rho_gradient[0] +
                                  p->rho_gradient[1] * p->rho_gradient[1] +
                                  p->rho_gradient[2] * p->rho_gradient[2];
    const double rho_grad_norm_inv =
        (rho_grad_norm2 > 0.) ? 1. / sqrt(rho_grad_norm2) : 0.;
    const double rho_com = hydro_get_comoving_density(p);
    double L_eff_com = rho_com * rho_grad_norm_inv;
    const double L_eff_com_max = kernel_gamma * p->h;
    const double L_eff_com_min = MIN_SHIELD_H_FRAC * p->h;
    L_eff_com = fmin(L_eff_com, L_eff_com_max);
    L_eff_com = fmax(L_eff_com, L_eff_com_min);
    const double L_eff_in_cm = L_eff_com * a * cooling->units.length_units;
    const double rho_to_n_cgs =
        cooling->units.density_units * 5.97729e23 * 0.75;
    const double rho_cgs_phys = rho_com * a3_inv * rho_to_n_cgs;
    const double NH_cgs = rho_cgs_phys * L_eff_in_cm;
    const double xH = NH_cgs * 3.50877e-24;
    const double fH_shield = pow(1.f + xH, -1.62) * exp(-0.149 * xH);

    fH2_shield *= fH_shield;
    /* Extra self-shielding from H2 if present - DON'T DO THIS HERE SINCE IT IS IN CRACKLE 
    const float fH2 = p->sf_data.H2_fraction; 
    if (fH2 > 0.f) { 
      const double NH2_cgs = fH2 * NH_cgs; 
      const double DH2_cgs = 1.e-5 * sqrt(2.*1.38e-16 * T_ism * 2.98864e23); 
      const double xH2 = NH2_cgs * 1.18133e-14; 
      fH2_shield *= 0.9379 * pow(1.f + xH2 / DH2_cgs, -1.879) +
          0.03465 * pow(1.f + xH2, -0.473) * exp(-2.293e-4 * sqrt(1.f + xH2));
    } */
  }
#endif

  return fH2_shield;
}

/**
 * @brief Returns the value of G0 for given particle p
 *
 * @param p Pointer to the particle data.
 * @param rho Physical density in system units.
 * @param cooling The properties of the cooling function.
 * @param dt The cooling timestep.
 *
 */
__attribute__((always_inline)) INLINE static float cooling_compute_G0(
    const struct part *restrict p, const float rho,
    const struct cooling_function_data *cooling, const float mstar,
    const float ssfr, const double dt) {

  float G0 = 0.f;
  float fH2_shield = 1.f;
  /* Determine ISRF in Habing units based on chosen method */
  if (cooling->G0_computation_method == 0) {
    G0 = 0.f;
  } 
  else if (cooling->G0_computation_method == 1) {
    fH2_shield = cooling_compute_self_shielding(p, cooling);
    G0 = fH2_shield * p->chemistry_data.local_sfr_density * cooling->G0_factor1;
  } 
  else if (cooling->G0_computation_method == 2) {
    G0 = ssfr * cooling->G0_factor2;
  } 
  else if (cooling->G0_computation_method == 3) {
    if (ssfr > 0.) {
      G0 = ssfr * cooling->G0_factor2;
    } 
    else {
      fH2_shield = cooling_compute_self_shielding(p, cooling);
      G0 = fH2_shield * p->chemistry_data.local_sfr_density *
           cooling->G0_factor1;
    }
  } 
  else if (cooling->G0_computation_method == -3) {
    if (p->chemistry_data.local_sfr_density > 0.) {
      fH2_shield = cooling_compute_self_shielding(p, cooling);
      G0 = fH2_shield * p->chemistry_data.local_sfr_density *
           cooling->G0_factor1;
    } 
    else {
      G0 = ssfr * cooling->G0_factor2;
    }
  }
#if COOLING_GRACKLE_MODE >= 2
  else if (cooling->G0_computation_method == 4) {
    /* Remember SNe_ThisTimeStep stores SN **rate** */
    G0 = p->cooling_data.SNe_ThisTimeStep * cooling->G0_factorSNe * dt;
  } 
  else if (cooling->G0_computation_method == 5) {
    float pssfr = max(p->sf_data.SFR, 0.f);
    pssfr /= max(mstar, 8. * p->mass);
    G0 = max(ssfr, pssfr) * cooling->G0_factor2 +
         p->cooling_data.SNe_ThisTimeStep * cooling->G0_factorSNe * dt;
  }
#endif
  else {
    error("G0_computation_method %d not recognized\n",
          cooling->G0_computation_method);
  }

  /* Scale G0 by user-input value */
  G0 *= cooling->G0_multiplier;

  if (mstar * 1.e10 > 1.e9 && p->id % 100000 == 0 && p->chemistry_data.local_sfr_density > 0) {
    message("G0: id=%lld M*=%g SFR=%g rho_sfr=%g SNe=%g Td=%g fshield=%g G0=%g",
            p->id,
            mstar * 1.e10,
            mstar * 1.e10 *
                ssfr / (1.e6 * cooling->time_to_Myr),
            p->chemistry_data.local_sfr_density * 0.002 / 1.6,
	    p->cooling_data.SNe_ThisTimeStep,
            p->cooling_data.dust_temperature,
            fH2_shield,
            G0);
  }

  return G0;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  if (cooling->chemistry.use_grackle == 1) {
    message("Cooling function is 'Grackle', mode = %i",
            cooling->chemistry.primordial_chemistry);
  } else if (cooling->chemistry.use_grackle == 2) {
    message("Cooling function is 'Crackle', mode = %i",
            cooling->chemistry.primordial_chemistry);
  }

  message("CloudyTable = %s", cooling->cloudy_table);
  message("Redshift = %g", cooling->redshift);
  message("UV background flag = %d", cooling->with_uv_background);
  message("Metal cooling flag = %i", cooling->chemistry.metal_cooling);
  message("Self Shielding flag = %i", cooling->self_shielding_method);
  message("Max subgrid density (internal units) = %g",
          cooling->max_subgrid_density);
  message("Thermal time = %g", cooling->thermal_time);
  message("Specific Heating Rates flag = %i",
          cooling->provide_specific_heating_rates);
  message("Volumetric Heating Rates flag = %i",
          cooling->provide_volumetric_heating_rates);
  message("Units:");
  message("\tComoving = %i", cooling->units.comoving_coordinates);
  message("\tLength = %g", cooling->units.length_units);
  message("\tNumber Density = %g", cooling->units.density_units);
  message("\tTime = %g", cooling->units.time_units);
  message("\tScale Factor = %g (units: %g)", cooling->units.a_value,
          cooling->units.a_units);
}

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE >= 1
void cooling_copy_to_grackle1(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  /* HI */
  species_densities[0] = xp->cooling_data.HI_frac * rho;
  data->HI_density = &species_densities[0];

  /* HII */
  species_densities[1] = xp->cooling_data.HII_frac * rho;
  data->HII_density = &species_densities[1];

  /* HeI */
  species_densities[2] = xp->cooling_data.HeI_frac * rho;
  data->HeI_density = &species_densities[2];

  /* HeII */
  species_densities[3] = xp->cooling_data.HeII_frac * rho;
  data->HeII_density = &species_densities[3];

  /* HeIII */
  species_densities[4] = xp->cooling_data.HeIII_frac * rho;
  data->HeIII_density = &species_densities[4];

  /* electrons */
  species_densities[5] = xp->cooling_data.e_frac * rho;
  data->e_density = &species_densities[5];
}
#else
void cooling_copy_to_grackle1(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  data->HI_density = NULL;
  data->HII_density = NULL;
  data->HeI_density = NULL;
  data->HeII_density = NULL;
  data->HeIII_density = NULL;
  data->e_density = NULL;
}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE >= 2
void cooling_copy_to_grackle2(
    grackle_field_data *data, const struct part *p, const struct xpart *xp,
    const struct cooling_function_data *restrict cooling, const double dt,
    gr_float rho, gr_float species_densities[N_SPECIES]) {
  /* HM */
  species_densities[6] = xp->cooling_data.HM_frac * rho;
  data->HM_density = &species_densities[6];

  /* H2I */
  species_densities[7] = xp->cooling_data.H2I_frac * rho;
  data->H2I_density = &species_densities[7];

  /* H2II */
  species_densities[8] = xp->cooling_data.H2II_frac * rho;
  data->H2II_density = &species_densities[8];

  /* Dust model */
  if (cooling->use_grackle_dust_evol == 1) {
    /* Load dust and metal info */
    species_densities[20] =
        p->cooling_data.dust_mass / p->mass * species_densities[12];
    data->dust_density = &species_densities[20];
    species_densities[21] =
        p->cooling_data.SNe_ThisTimeStep * dt / p->mass * species_densities[12];
    /* need to pass the number of SNe per volume of particle;
     * recall SNe_ThisTimeStep is the SNe rate */

    // if (species_densities[21] > 1.e15) message("SNe_density: %g %g %g
    // %g\n",p->mass, species_densities[12], p->cooling_data.SNe_ThisTimeStep,
    // species_densities[21]);
    data->SNe_ThisTimeStep = &species_densities[21];
    // if( chemistry_get_total_metal_mass_fraction_for_cooling(p)>0.f)
    // message("Zsm= %g Zp= %g Z= %g Zd=
    // %g",chemistry_get_total_metal_mass_fraction_for_cooling(p),
    // p->chemistry_data.metal_mass_fraction_total, species_densities[19],
    // species_densities[20]);

    /* Determine ISRF in Habing units based on chosen method, -1 == non-ISM */
    if (p->cooling_data.subgrid_temp == 0.f) {
      species_densities[22] = -1.f;
    } else {
      species_densities[22] = p->cooling_data.G0;
    }

    data->isrf_habing = &species_densities[22];

    const double rho_grad_norm2 = p->rho_gradient[0] * p->rho_gradient[0] +
                                  p->rho_gradient[1] * p->rho_gradient[1] +
                                  p->rho_gradient[2] * p->rho_gradient[2];
    const double rho_grad_norm_inv =
        (rho_grad_norm2 > 0.) ? 1. / sqrt(rho_grad_norm2) : 0.;
    const double rho_com = hydro_get_comoving_density(p);
    double L_eff_com = rho_com * rho_grad_norm_inv;
    const double L_eff_com_max = kernel_gamma * p->h;
    const double L_eff_com_min = MIN_SHIELD_H_FRAC * p->h;
    L_eff_com = fmin(L_eff_com, L_eff_com_max);
    L_eff_com = fmax(L_eff_com, L_eff_com_min);

    species_densities[23] = L_eff_com * cooling->units.a_value;
    data->H2_self_shielding_length = &species_densities[23];

    /* Load gas metallicities */
    for (int i = 0; i < chemistry_element_count; i++) {
      species_densities[24 + i] = max(
          p->chemistry_data.metal_mass_fraction[i] * species_densities[12], 0.);
      species_densities[24 + chemistry_element_count + i] =
          max(p->cooling_data.dust_mass_fraction[i] * species_densities[20], 0);
      // if (i>0) printf("dust densities: %d %g %g
      // %g\n",i,species_densities[23+chemistry_element_count+i],species_densities[23+chemistry_element_count+i]/data->dust_density[0],p->cooling_data.dust_mass_fraction[i]*p->mass
      // / p->cooling_data.dust_mass);
    }

    data->He_gas_metalDensity = &species_densities[25];
    data->C_gas_metalDensity = &species_densities[26];
    data->N_gas_metalDensity = &species_densities[27];
    data->O_gas_metalDensity = &species_densities[28];
    data->Ne_gas_metalDensity = &species_densities[29];
    data->Mg_gas_metalDensity = &species_densities[30];
    data->Si_gas_metalDensity = &species_densities[31];
    data->S_gas_metalDensity = &species_densities[32];
    data->Ca_gas_metalDensity = &species_densities[32];
    data->Fe_gas_metalDensity = &species_densities[34];
    /* Load dust metallicities */
    data->He_dust_metalDensity = &species_densities[36];
    data->C_dust_metalDensity = &species_densities[37];
    data->N_dust_metalDensity = &species_densities[38];
    data->O_dust_metalDensity = &species_densities[39];
    data->Ne_dust_metalDensity = &species_densities[40];
    data->Mg_dust_metalDensity = &species_densities[41];
    data->Si_dust_metalDensity = &species_densities[42];
    data->S_dust_metalDensity = &species_densities[43];
    data->Ca_dust_metalDensity = &species_densities[44];
    data->Fe_dust_metalDensity = &species_densities[45];
  } else {
    data->dust_density = NULL;
    data->SNe_ThisTimeStep = NULL;
    data->isrf_habing = NULL;
    data->He_gas_metalDensity = NULL;
    data->C_gas_metalDensity = NULL;
    data->N_gas_metalDensity = NULL;
    data->O_gas_metalDensity = NULL;
    data->Ne_gas_metalDensity = NULL;
    data->Mg_gas_metalDensity = NULL;
    data->Si_gas_metalDensity = NULL;
    data->S_gas_metalDensity = NULL;
    data->Ca_gas_metalDensity = NULL;
    data->Fe_gas_metalDensity = NULL;
    data->He_dust_metalDensity = NULL;
    data->C_dust_metalDensity = NULL;
    data->N_dust_metalDensity = NULL;
    data->O_dust_metalDensity = NULL;
    data->Ne_dust_metalDensity = NULL;
    data->Mg_dust_metalDensity = NULL;
    data->Si_dust_metalDensity = NULL;
    data->S_dust_metalDensity = NULL;
    data->Ca_dust_metalDensity = NULL;
    data->Fe_dust_metalDensity = NULL;
  }
}
#else
void cooling_copy_to_grackle2(
    grackle_field_data *data, const struct part *p, const struct xpart *xp,
    const struct cooling_function_data *restrict cooling, const double dt,
    gr_float rho, gr_float species_densities[N_SPECIES]) {
  data->HM_density = NULL;
  data->H2I_density = NULL;
  data->H2II_density = NULL;
}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part
 * @param xp The #xpart
 * @param rho Particle density
 */
#if COOLING_GRACKLE_MODE >= 3
void cooling_copy_to_grackle3(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  /* DI */
  species_densities[9] = xp->cooling_data.DI_frac * rho;
  data->DI_density = &species_densities[9];

  /* DII */
  species_densities[10] = xp->cooling_data.DII_frac * rho;
  data->DII_density = &species_densities[10];

  /* HDI */
  species_densities[11] = xp->cooling_data.HDI_frac * rho;
  data->HDI_density = &species_densities[11];
}
#else
void cooling_copy_to_grackle3(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
  data->DI_density = NULL;
  data->DII_density = NULL;
  data->HDI_density = NULL;
}
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
#if COOLING_GRACKLE_MODE >= 1
void cooling_copy_from_grackle1(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {

  const double rhoinv = 1.f / rho;
  /* HI */
  xp->cooling_data.HI_frac = *data->HI_density * rhoinv;

  /* HII */
  xp->cooling_data.HII_frac = *data->HII_density * rhoinv;

  /* HeI */
  xp->cooling_data.HeI_frac = *data->HeI_density * rhoinv;

  /* HeII */
  xp->cooling_data.HeII_frac = *data->HeII_density * rhoinv;

  /* HeIII */
  xp->cooling_data.HeIII_frac = *data->HeIII_density * rhoinv;

  /* e */
  xp->cooling_data.e_frac = *data->e_density * rhoinv;
}
#else
void cooling_copy_from_grackle1(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {}
#endif

/**
 * @brief copy the grackle data to a #xpart.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
#if COOLING_GRACKLE_MODE >= 2
void cooling_copy_from_grackle2(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho) {
  double rhoinv = 1.f / rho;

  /* HM */
  xp->cooling_data.HM_frac = *data->HM_density * rhoinv;
  /* H2I */
  xp->cooling_data.H2I_frac = *data->H2I_density * rhoinv;
  /* H2II */
  xp->cooling_data.H2II_frac = *data->H2II_density * rhoinv;

  /* Dust model */
  if (cooling->use_grackle_dust_evol == 1) {
    /* Load gas metallicities */
    p->chemistry_data.metal_mass_fraction_total = *data->metal_density * rhoinv;
    p->chemistry_data.metal_mass_fraction[1] =
        *data->He_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[2] =
        *data->C_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[3] =
        *data->N_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[4] =
        *data->O_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[5] =
        *data->Ne_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[6] =
        *data->Mg_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[7] =
        *data->Si_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[8] =
        *data->S_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[9] =
        *data->Ca_gas_metalDensity * rhoinv;
    p->chemistry_data.metal_mass_fraction[10] =
        *data->Fe_gas_metalDensity * rhoinv;

    /* Load dust metallicities */
    p->cooling_data.dust_mass = *data->dust_density * p->mass * rhoinv;

    if (p->cooling_data.dust_mass > 0.5 * p->mass) {
      warning("DUST > METALS Mg=%g Zg=%g mdust=%g mmet=%g\n", p->mass,
              chemistry_get_total_metal_mass_fraction_for_cooling(p),
              p->cooling_data.dust_mass,
              p->mass * chemistry_get_total_metal_mass_fraction_for_cooling(p));
    }

    p->cooling_data.dust_mass_fraction[0] = 0.f;
    if (*data->dust_density > 0.f) {
      rhoinv = 1. / *data->dust_density;

      const double rhodust[chemistry_element_count - 1] = {
          *data->He_dust_metalDensity, *data->C_dust_metalDensity,
          *data->N_dust_metalDensity,  *data->O_dust_metalDensity,
          *data->Ne_dust_metalDensity, *data->Mg_dust_metalDensity,
          *data->Si_dust_metalDensity, *data->S_dust_metalDensity,
          *data->Ca_dust_metalDensity, *data->Fe_dust_metalDensity};

      /* no Helium metal density */
      assert(rhodust[0] == 0.);

      for (int i = 1; i < chemistry_element_count; i++) {
        p->cooling_data.dust_mass_fraction[i] = rhodust[i - 1] * rhoinv;
        p->cooling_data.dust_mass_fraction[0] +=
            p->cooling_data.dust_mass_fraction[i];
      }
    } else {
      for (int i = 1; i < chemistry_element_count; i++) {
        p->cooling_data.dust_mass_fraction[i] = 0.f;
      }
    }
  }
}
#else
void cooling_copy_from_grackle2(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho) {}
#endif

/**
 * @brief copy the grackle data to a #xpart
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
#if COOLING_GRACKLE_MODE > 2
void cooling_copy_from_grackle3(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {

  const double rhoinv = 1.f / rho;
  /* DI */
  xp->cooling_data.DI_frac = *data->DI_density * rhoinv;

  /* DII */
  xp->cooling_data.DII_frac = *data->DII_density * rhoinv;

  /* HDI */
  xp->cooling_data.HDI_frac = *data->HDI_density * rhoinv;
}
#else
void cooling_copy_from_grackle3(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho) {}
#endif

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
void cooling_copy_to_grackle(
    grackle_field_data *data, const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling, const struct part *p,
    const struct xpart *xp, const double dt, const double T_warm,
    gr_float species_densities[N_SPECIES], gr_float *iact_rates, int mode) {

  int i;
  /* set values */
  /* grid */
  data->grid_dx = 0.f;
  data->grid_rank = GRACKLE_RANK;

  data->grid_dimension = malloc(GRACKLE_RANK * sizeof(int));
  data->grid_start = malloc(GRACKLE_RANK * sizeof(int));
  data->grid_end = malloc(GRACKLE_RANK * sizeof(int));
  for (i = 0; i < 3; i++) {
    /* The active dimension not including ghost zones */
    data->grid_dimension[i] = 1;
    data->grid_start[i] = 0;
    data->grid_end[i] = 0;
  }

  data->grid_dimension[0] = GRACKLE_NPART;
  data->grid_end[0] = GRACKLE_NPART - 1;

  /* get particle density, internal energy in
     physical coordinates (still code units) */
  double T_subgrid = cooling_get_subgrid_temperature(p, xp);
  /* mode 1 is cooling time, here we want non-subgrid values in all cases */
  if (mode == 1) {
    species_densities[12] = hydro_get_physical_density(p, cosmo);
    species_densities[13] = hydro_get_physical_internal_energy(p, xp, cosmo);
    species_densities[14] = cooling->T_CMB_0 * (1.f + cosmo->z);
    species_densities[15] = 0.f;
    data->grid_end[0] = -1;  // this signals to crackle to turn off UVB
  }
  /* non-subgrid case, here we set the floor temperature
   * by the EoS (if applicable).  Note the cold fraction has a small
   * limit otherwise one can get underflows in crackle. */
  else if (p->cooling_data.subgrid_temp == 0. ||
           p->cooling_data.subgrid_fcold <= 1.e-6) {
    species_densities[12] = hydro_get_physical_density(p, cosmo);
    species_densities[13] = hydro_get_physical_internal_energy(p, xp, cosmo);
    species_densities[14] = T_warm;
    /* specific_heating_rate has to be in cgs units;
       no unit conversion done within grackle */
    species_densities[15] =
        hydro_get_physical_internal_energy_dt(p, cosmo) * cooling->dudt_units;
  }
  /* subgrid ISM case, use subgrid ISM values and set floor by T_CMB */
  else {
    /* Physical sub-grid density*/
    species_densities[12] =
        cooling_get_subgrid_density(p, xp) * p->cooling_data.subgrid_fcold;
    /* Physical internal energy */
    species_densities[13] = cooling_convert_temp_to_u(
        T_subgrid, xp->cooling_data.e_frac, cooling, p);
    /* CMB temp is floor*/
    species_densities[14] = cooling->T_CMB_0 * (1.f + cosmo->z);
    /* If tracking H2, turn off specific heating rate in ISM. */
    species_densities[15] = 0.f;
    if (cooling->ism_adiabatic_heating_method == 1) {
      species_densities[15] +=
          hydro_get_physical_internal_energy_dt(p, cosmo) * cooling->dudt_units;
    }
  }
  /* load into grackle structure */
  data->density = &species_densities[12];
  data->internal_energy = &species_densities[13];
  data->temperature_floor = &species_densities[14];
  data->specific_heating_rate = &species_densities[15];

  /* velocity (maybe not needed?) */
  species_densities[16] = xp->v_full[0] * cosmo->a_inv;
  species_densities[17] = xp->v_full[1] * cosmo->a_inv;
  species_densities[18] = xp->v_full[2] * cosmo->a_inv;
  data->x_velocity = &species_densities[16];
  data->y_velocity = &species_densities[17];
  data->z_velocity = &species_densities[18];

  cooling_copy_to_grackle1(data, p, xp, species_densities[12],
                           species_densities);
  cooling_copy_to_grackle2(data, p, xp, cooling, dt, species_densities[12],
                           species_densities);
  cooling_copy_to_grackle3(data, p, xp, species_densities[12],
                           species_densities);

  /* RT heating and ionisation rates */
  data->RT_heating_rate = &iact_rates[0];
  data->RT_HI_ionization_rate = &iact_rates[1];
  data->RT_HeI_ionization_rate = &iact_rates[2];
  data->RT_HeII_ionization_rate = &iact_rates[3];
  data->RT_H2_dissociation_rate = &iact_rates[4];

  species_densities[19] =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) *
      species_densities[12];
  data->metal_density = &species_densities[19];

  for (i = 0; i < N_SPECIES; i++) {
    if (fpclassify(species_densities[i]) == FP_NAN ||
        fpclassify(species_densities[i]) == FP_INFINITE) {
      error(
          "Passing a non-finite value to grackle! "
          "i=%d / %d, species_densities[i]=%g\n",
          i, N_SPECIES, species_densities[i]);
    }
  }
}

/**
 * @brief copy a #xpart to the grackle data
 *
 * Warning this function creates some variable, therefore the grackle call
 * should be in a block that still has the variables.
 *
 * @param data The grackle_field_data structure from grackle.
 * @param p The #part.
 * @param xp The #xpart.
 * @param rho The particle density.
 */
void cooling_copy_from_grackle(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho) {

  cooling_copy_from_grackle1(data, p, xp, rho);
  cooling_copy_from_grackle2(data, p, xp, cooling, rho);
  cooling_copy_from_grackle3(data, p, xp, rho);
}

/**
 * @brief free memory associated with grackle driver
 *
 * @param data The grackle_field_data structure from grackle.
 */
void cooling_grackle_free_data(grackle_field_data *data) {

  free(data->grid_dimension);
  free(data->grid_start);
  free(data->grid_end);
}

/**
 * @brief Compute the energy of a particle after dt and update the particle
 * chemistry data
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param mode 0=energy, 1=cooling time, 2=temperature, 3=pressure, 4=gamma
 *
 * @return desired quantity based on mode
 */
gr_float cooling_grackle_driver(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    double dt, double T_warm, int mode) {

  /* set current units for conversion to physical quantities */
  code_units units = cooling->units;

  /* initialize data to send to grackle */
  gr_float *species_densities;
  species_densities = (gr_float *)calloc(N_SPECIES, sizeof(gr_float));
  grackle_field_data data;
  // cooling_grackle_malloc_fields(&data, 1, cooling->chemistry.use_dust_evol);

  /* load particle information from particle to grackle data */
  cooling_copy_to_grackle(&data, us, cosmo, cooling, p, xp, dt, T_warm,
                          species_densities, iact_rates, mode);

  /* Run Grackle in desired mode */
  gr_float return_value = 0.f;
  double t_dust = 0.f;

  switch (mode) {
    case 0:
      /* solve chemistry, advance thermal energy by dt */
      if (solve_chemistry(&units, &data, dt) == 0) {
        error("Error in Grackle solve_chemistry.");
      }
      // if (solve_chemistry(&units, &data, -dt) == 0) {
      //   error("Error in Crackle solve_chemistry.");
      // }
      /* copy from grackle data to particle */
      cooling_copy_from_grackle(&data, p, xp, cooling, species_densities[12]);
      return_value = data.internal_energy[0];
#if COOLING_GRACKLE_MODE >= 2
      /* Compute dust temperature */
      t_dust = p->cooling_data.dust_temperature;
      if (calculate_dust_temperature(&units, &data, &t_dust) == 0) {
        error("Error in Grackle calculate dust temperature.");
      }

      p->cooling_data.dust_temperature = t_dust;

      /* Reset accumulated local variables to zero */
      p->cooling_data.SNe_ThisTimeStep = 0.f;
#endif
      break;

    case 1:
      /* compute cooling time */
      if (calculate_cooling_time(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_cooling_time.");
      }
      break;

    case 2:
      /* compute temperature */
      if (calculate_temperature(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_temperature.");
      }
      break;

    case 3:
      /* compute pressure */
      if (calculate_pressure(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_pressure.");
      }
      break;
    case 4:
      /* compute gamma */
      if (calculate_gamma(&units, &data, &return_value) == 0) {
        error("Error in Grackle calculate_gamma.");
      }
      break;
  }
  cooling_grackle_free_data(&data);
  free(species_densities);

  return return_value;
}

/**
 * @brief Compute the cooling time. Optionally uses input values
 * for rho and u, but leaves all other properties of p the same.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data.
 * @param rhocool Density used in tcool calculation.
 * @param ucool Density used in tcool calculation.
 *
 * @return cooling time
 */
gr_float cooling_time(const struct phys_const *restrict phys_const,
                      const struct unit_system *restrict us,
                      const struct hydro_props *hydro_properties,
                      const struct cosmology *restrict cosmo,
                      const struct cooling_function_data *restrict cooling,
                      const struct part *restrict p, struct xpart *restrict xp,
                      const float rhocool, const float ucool) {

  /* Removes const in declaration*/
  struct part p_temp = *p;

  if (rhocool > 0.f) p_temp.rho = rhocool;
  if (ucool > 0.f) p_temp.u = ucool;

  gr_float iact_rates[5] = {0., 0., 0., 0., 0.};

  gr_float cooling_time =
      cooling_grackle_driver(phys_const, us, cosmo, hydro_properties, cooling,
                             &p_temp, xp, iact_rates, 0., 0., 1);

  return cooling_time;
}

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
float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_properties,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp) {

  /* Remove const in declaration*/
  struct part p_temp = *p;
  struct xpart xp_temp = *xp;

  gr_float iact_rates[5] = {0., 0., 0., 0., 0.};

  const float temperature =
      cooling_grackle_driver(phys_const, us, cosmo, hydro_properties, cooling,
                             &p_temp, &xp_temp, iact_rates, 0., 0., 2);

  return temperature;
}

/**
 * @brief Do dust sputtering and return metals back to gas phase
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the #xpart data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE void cooling_sputter_dust(
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, const double dt) {

  /* Do dust destruction in stream particle */
  if (cooling->use_grackle_dust_evol && p->cooling_data.dust_mass > 0.f) {
    const float u_phys = hydro_get_physical_internal_energy(p, xp, cosmo);
    const float Tstream =
        cooling_convert_u_to_temp(u_phys, xp->cooling_data.e_frac, cooling, p);
    const double Tstream_K =
        Tstream * units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    /* Sputtering negligible at low-T */
    if (Tstream_K > 1.e4) {
      const double rho_cgs = hydro_get_physical_density(p, cosmo) *
                             units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

      /* sputtering timescale, Tsai & Mathews (1995) */
      const double tsp = 1.7e8 * 3.15569251e7 /
                         units_cgs_conversion_factor(us, UNIT_CONV_TIME) *
                         (cooling->dust_grainsize / 0.1) * (1.e-27 / rho_cgs) *
                         (pow(2.e6 / Tstream, 2.5) + 1.0);

      const float dust_mass_old = p->cooling_data.dust_mass;

      /* Update dust mass but limit the destruction */
      p->cooling_data.dust_mass -=
          dust_mass_old * (1.f - exp(-3.f * min(dt / tsp, 5.f)));
      if (p->cooling_data.dust_mass < 0.5f * dust_mass_old) {
        p->cooling_data.dust_mass = 0.5f * dust_mass_old;
      }

#ifdef FIREHOSE_DEBUG_CHECKS
      message(
          "FIREHOSE_SPUT: id=%lld mdust=%g mdustnew=%g T=%g rho=%g "
          "tsp_dt_ratio=%g",
          p->id, dust_mass_old, p->cooling_data.dust_mass, Tstream,
          rho_cgs / 1.673e-24, /* units of mp */
          tsp / dt);
#endif
      /* factor by which dust mass changed */
      const float dust_mass_new = p->cooling_data.dust_mass;
      const float dust_mass_ratio = dust_mass_new / dust_mass_old;
      p->chemistry_data.metal_mass_fraction_total = 0.f;

      for (int elem = chemistry_element_He; elem < chemistry_element_count;
           ++elem) {
        const float Z_dust_elem_old = p->cooling_data.dust_mass_fraction[elem];
        const float Z_dust_elem_new = Z_dust_elem_old * dust_mass_ratio;
        const float Z_elem_old = p->chemistry_data.metal_mass_fraction[elem];
        const float elem_mass_old = Z_elem_old * hydro_get_mass(p);

        /* This is the positive amount of metal mass to add since we
         * are losing dust mass when sputtering. */
        const float delta_metal_mass_elem =
            (Z_dust_elem_old * dust_mass_old - Z_dust_elem_new * dust_mass_new);
        const float elem_mass_new = elem_mass_old + delta_metal_mass_elem;
        const float Z_elem_new = elem_mass_new / hydro_get_mass(p);

        p->chemistry_data.metal_mass_fraction[elem] = Z_elem_new;
        p->cooling_data.dust_mass_fraction[elem] *= dust_mass_ratio;

        /* Sum up to get the new Z value */
        if (elem != chemistry_element_H && elem != chemistry_element_He) {
          p->chemistry_data.metal_mass_fraction_total += Z_elem_new;
        }
      }

      /* Make sure that X + Y + Z = 1 */
      const float Y_He =
          p->chemistry_data.metal_mass_fraction[chemistry_element_He];
      p->chemistry_data.metal_mass_fraction[chemistry_element_H] =
          1.f - Y_He - p->chemistry_data.metal_mass_fraction_total;

      /* Make sure H fraction does not go out of bounds */
      if (p->chemistry_data.metal_mass_fraction[chemistry_element_H] > 1.f ||
          p->chemistry_data.metal_mass_fraction[chemistry_element_H] < 0.f) {
        for (int i = chemistry_element_H; i < chemistry_element_count; i++) {
          warning("\telem[%d] is %g", i,
                  p->chemistry_data.metal_mass_fraction[i]);
        }

        error(
            "Hydrogen fraction exeeds unity or is negative for"
            " particle id=%lld due to dust sputtering",
            p->id);
      }
    }
  }
}

/**
 * @brief Compute particle quantities for the firehose model
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the #xpart data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE void firehose_cooling_and_dust(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *restrict hydro_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, const double dt) {

  /* Initialize cooling time */
  p->cooling_data.mixing_layer_cool_time = 0.f;

  const float u = hydro_get_comoving_internal_energy(p, xp);

  /* If it's not a firehose particle, just compute particle cooling time */
  if (p->chemistry_data.radius_stream <= 0.f ||
      p->chemistry_data.rho_ambient <= 0.f) {
    p->cooling_data.mixing_layer_cool_time = cooling_time(
        phys_const, us, hydro_props, cosmo, cooling, p, xp, p->rho, u);
    return;
  }

  /* It's a firehose particles, so compute the cooling rate
   * in the mixing layer */
  const float rhocool = 0.5f * (p->chemistry_data.rho_ambient + p->rho);
  const float ucool = 0.5f * (p->chemistry_data.u_ambient + u);

  /* +ive if heating -ive if cooling*/
  p->cooling_data.mixing_layer_cool_time = cooling_time(
      phys_const, us, hydro_props, cosmo, cooling, p, xp, rhocool, ucool);
#ifdef FIREHOSE_DEBUG_CHECKS
  message("FIREHOSE_COOLING: id=%lld nH=%g T=%g rhoamb=%g Tamb=%g tcool=%g",
          p->id,
          hydro_get_physical_density(p, cosmo) * cooling->units.density_units *
              0.75 / 1.673e-24,
          u_old * cosmo->a_factor_internal_energy / cooling->temp_to_u_factor,
          p->chemistry_data.rho_ambient,
          p->chemistry_data.u_ambient * cosmo->a_factor_internal_energy /
              cooling->temp_to_u_factor,
          p->cooling_data.mixing_layer_cool_time);
#endif

  cooling_sputter_dust(us, cosmo, cooling, p, xp, dt);
}

/**
 * @brief Set subgrid ISM properties and initialise chemistry
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 */
void cooling_init_chemistry(
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p) {

  /* If it's eligible for SF and metal-free, crudely self-enrich to
     very small level; needed to kick-start Grackle dust. Arbitrary
     amount to kick-start the model */
  const float init_dust_to_gas = 0.2f;
  const float total_Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);
  const float self_Z =
      (1.f - init_dust_to_gas) * cooling->self_enrichment_metallicity;
  if (p->cooling_data.subgrid_temp > 0.f && total_Z < self_Z) {
    float Z_sun = 0.f;
    for (int i = 1; i < 10; i++) {
      Z_sun += cooling->chemistry.SolarAbundances[i];
    }

    /* Distribute the self-enrichment metallicity among elements
       assuming solar abundance ratios*/
    p->chemistry_data.metal_mass_fraction_total = 0.f;
    p->cooling_data.dust_mass = 0.f;
    /* Offset index for the SolarAbundaces array (starts at He -> Fe) */
    int j = 1;
    for (int i = chemistry_element_C; i < chemistry_element_count; i++) {
      /* fraction of gas mass in each element */
      p->chemistry_data.metal_mass_fraction[i] =
          (1.f - init_dust_to_gas) * cooling->chemistry.SolarAbundances[j] *
          cooling->self_enrichment_metallicity;
      p->chemistry_data.metal_mass_fraction[i] /= Z_sun;

      /* Update to the new metal mass fraction */
      p->chemistry_data.metal_mass_fraction_total +=
          p->chemistry_data.metal_mass_fraction[i];

      /* fraction of dust mass in each element */
      p->cooling_data.dust_mass_fraction[i] =
          init_dust_to_gas * cooling->chemistry.SolarAbundances[j];
      p->cooling_data.dust_mass_fraction[i] /= Z_sun;

      /* Sum up all of the dust mass */
      p->cooling_data.dust_mass += p->cooling_data.dust_mass_fraction[i] *
                                   p->chemistry_data.metal_mass_fraction[i] *
                                   p->mass;
      j++;
    }

    p->chemistry_data.metal_mass_fraction[chemistry_element_He] =
        cooling->chemistry.SolarAbundances[0];
    /* Since He is fixed at SolarAbundances[0], make sure the hydrogen
       fraction makes sense, i.e. X_H + Y_He + Z = 1. */
    p->chemistry_data.metal_mass_fraction[chemistry_element_H] =
        1.f - p->chemistry_data.metal_mass_fraction[chemistry_element_He] -
        p->chemistry_data.metal_mass_fraction_total;

    if (p->chemistry_data.metal_mass_fraction[chemistry_element_H] > 1.f ||
        p->chemistry_data.metal_mass_fraction[chemistry_element_H] < 0.f) {
      for (int i = chemistry_element_H; i < chemistry_element_count; i++) {
        warning("\telem[%d] is %g", i,
                p->chemistry_data.metal_mass_fraction[i]);
      }

      error(
          "Hydrogen fraction exeeds unity or is negative for"
          " particle id=%lld",
          p->id);
    }
  }
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param iact_rates Interaction rates for radiative transfer (if used)
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
void cooling_do_grackle_cooling(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    const double dt, const double dt_therm) {

  /* Self-enrich gas if very low metallicity */
  cooling_init_chemistry(cooling, p);

  /* Collect information about galaxy that the particle belongs to */
  const float galaxy_mstar = p->galaxy_data.stellar_mass;
  const float galaxy_ssfr = p->galaxy_data.specific_sfr;

  /* Compute the ISRF */
  p->cooling_data.G0 =
      fmax(cooling_compute_G0(p, p->cooling_data.subgrid_dens, cooling,
                              galaxy_mstar, galaxy_ssfr, dt),
           0.);

  /* Compute the entropy floor */
  // const double T_warm = entropy_floor_temperature(p, cosmo, floor_props);
  const double T_warm = warm_ISM_temperature(p, cooling, phys_const, cosmo);
  const double u_warm =
      cooling_convert_temp_to_u(T_warm, xp->cooling_data.e_frac, cooling, p);

  /* Do grackle cooling */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
  // const float T_old = cooling_get_temperature( phys_const, hydro_props, us,
  // cosmo, cooling, p, xp); // for debugging only
  gr_float u_new = u_old;
  u_new = cooling_grackle_driver(phys_const, us, cosmo, hydro_props, cooling, p,
                                 xp, iact_rates, dt, T_warm, 0);

  /* Apply simulation-wide minimum temperature */
  u_new = max(u_new, hydro_props->minimal_internal_energy);

  /* Assign new thermal energy to particle */
  /* Calculate the cooling rate */
  float cool_du_dt = (u_new - u_old) / dt_therm;

  if (p->cooling_data.subgrid_temp == 0.) {
    /* Normal cooling; check that we are not going to go
     * below any of the limits */
    if (u_new > GRACKLE_HEATLIM * u_old) u_new = GRACKLE_HEATLIM * u_old;
    if (u_new < GRACKLE_COOLLIM * u_old) u_new = GRACKLE_COOLLIM * u_old;
    // u_new = max(u_new, u_warm);

    /* Rennehan: Recompute the actual thermal evolution after setting min/max */
    cool_du_dt = (u_new - u_old) / dt_therm;

    /* Update the internal energy time derivative,
     * which will be evolved later */
    hydro_set_physical_internal_energy_dt(p, cosmo, cool_du_dt);

    /* If there is any dust outside of the ISM, sputter it
     * back into gas phase metals */
    cooling_sputter_dust(us, cosmo, cooling, p, xp, dt);
  } else {
    /* Particle is in subgrid mode; result is stored in subgrid_temp */
    p->cooling_data.subgrid_temp =
        cooling_convert_u_to_temp(u_new, xp->cooling_data.e_frac, cooling, p);

    /* Set the subgrid cold ISM fraction for particle */
    /* Get H number density */
    const double rho = hydro_get_physical_density(p, cosmo);
    const float X_H =
        chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
    const double n_H = rho * X_H / phys_const->const_proton_mass;

    const double fcold_max = cooling_compute_cold_ISM_fraction(n_H, cooling);

    /* Compute cooling time in warm ISM component */
    float rhocool = p->rho * (1.f - p->cooling_data.subgrid_fcold);
    rhocool = fmax(rhocool, 0.001f * p->rho);
    const float u = hydro_get_comoving_internal_energy(p, xp);
    float tcool = cooling_time(phys_const, us, hydro_props, cosmo, cooling, p,
                               xp, rhocool, u);

    /* Evolve fcold upwards by cooling from warm ISM on
     * relevant cooling timescale */
    if (tcool < 0.f) {
      p->cooling_data.subgrid_fcold +=
          (fcold_max - p->cooling_data.subgrid_fcold) * (1.f - exp(dt / tcool));
    }

    /* Compare the adiabatic heating to the heating required to heat the
     * gas back up to the equation of state line, and assume this is the
     * fraction of cold cloud mass destroyed. */
    if (cooling->ism_adiabatic_heating_method == 2) {
      /* Use adiabatic du/dt to evaporate cold gas clouds, into warm phase */
      const double f_evap =
          hydro_get_physical_internal_energy_dt(p, cosmo) * dt_therm /
          (hydro_get_physical_internal_energy(p, xp, cosmo) - u_new);

      /* If it's in the ISM of a galaxy, suppress cold fraction */
      if (f_evap > 0.f && galaxy_mstar > 0.f) {
        p->cooling_data.subgrid_fcold *= max(1. - f_evap, 0.f);
      }
    }

    /* Set internal energy time derivative to 0 for overall particle */
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);

    /* No cooling in warm phase since it is fixed on the EoS */
    cool_du_dt = 0.f;

    /* Force the overall particle to lie on the equation of state
    hydro_set_physical_internal_energy(p, xp, cosmo, u_warm);*/

    /* set subgrid properties for use in SF routine */
    cooling_set_particle_subgrid_properties(phys_const, us, cosmo, hydro_props,
                                            floor_props, cooling, p, xp);

    /* Overall particle u is combination of EOS u and subgrid u */
    const float u_part = p->cooling_data.subgrid_fcold * u_new +
                         (1. - p->cooling_data.subgrid_fcold) * u_warm;
    hydro_set_physical_internal_energy(p, xp, cosmo, u_part);
  } /* subgrid mode */

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cool_du_dt * dt_therm;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle' extended data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
void cooling_cool_part(const struct phys_const *restrict phys_const,
                       const struct unit_system *restrict us,
                       const struct cosmology *restrict cosmo,
                       const struct hydro_props *hydro_props,
                       const struct entropy_floor_properties *floor_props,
                       const struct pressure_floor_props *pressure_floor_props,
                       const struct cooling_function_data *restrict cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const double dt, const double dt_therm,
                       const double time) {

  /* Compute cooling time and other quantities needed for firehose */
  firehose_cooling_and_dust(phys_const, us, cosmo, hydro_props, cooling, p, xp,
                            dt);

  /* Update the subgrid properties */
  cooling_set_particle_subgrid_properties(phys_const, us, cosmo, hydro_props,
                                          floor_props, cooling, p, xp);

  /* No cooling if particle is decoupled */
  if (p->decoupled) return;

  if (cooling->do_cooling_in_rt) return;

  /* No cooling happens over zero time */
  if (dt == 0.f || dt_therm == 0.f) return;

  /* Interaction rates for RT; not used here */
  gr_float iact_rates[5] = {0., 0., 0., 0., 0.};

  /* Do the cooling and chemistry */
  cooling_do_grackle_cooling(phys_const, us, cosmo, hydro_props, floor_props,
                             cooling, p, xp, iact_rates, dt, dt_therm);

  /* Record this cooling event */
  xp->cooling_data.time_last_event = time;
}

/**
 * @brief Set the subgrid properties (rho, T) of the gas particle for use
 *        in SF routine
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_set_particle_subgrid_properties(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp) {

  /* No subgrid ISM if particle is decoupled */
  if (p->decoupled) {
    /* Make sure these are always set for the wind particles */
    p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    p->cooling_data.subgrid_temp = 0.;
    p->cooling_data.subgrid_fcold = 0.f;

    return;
  }

  /* Get temperature of overall particle */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float temperature =
      cooling_convert_u_to_temp(u, xp->cooling_data.e_frac, cooling, p);

  /* Get density */
  const double rho = hydro_get_physical_density(p, cosmo);

  /* Subgrid model is on if particle is in the Jeans EOS regime */
  const double T_warm = warm_ISM_temperature(p, cooling, phys_const, cosmo);
  // entropy_floor_gas_temperature( rho, rho_com, cosmo, floor_props);
  const double u_warm =
      cooling_convert_temp_to_u(T_warm, xp->cooling_data.e_frac, cooling, p);

  /* Check if it is in subgrid mode: Must be in Jeans EoS regime
   * and have nonzero cold gas */
  if (T_warm > 0 && u < u_warm * cooling->entropy_floor_margin) {
    /* YES: If first time in subgrid, set temperature to particle T,
       otherwise limit to particle T */
    if (p->cooling_data.subgrid_temp == 0.) {
      p->cooling_data.subgrid_temp = temperature;
      /* Reset grackle subgrid quantities assuming neutral gas */
      cooling_grackle_init_part(cooling, p, xp);
      /* Initialize ISM cold fraction */
      p->cooling_data.subgrid_fcold = cooling->cold_ISM_frac;
    } else {
      /* Subgrid temperature should be no higher
       * than overall particle temperature */
      p->cooling_data.subgrid_temp =
          min(p->cooling_data.subgrid_temp, temperature);
    }

    /* Compute subgrid density assuming pressure equilibrium */
    const float X_H =
        chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
    const double n_H = rho * X_H / phys_const->const_proton_mass;
    p->cooling_data.subgrid_dens = cooling_compute_subgrid_density(
        rho, n_H, temperature, p->cooling_data.subgrid_temp, cooling);
  } else {
    /* NO: subgrid density is the actual particle's physical density */
    p->cooling_data.subgrid_dens = rho;

    /* set subgrid temperature to 0 indicating it's not in subgrid mode */
    p->cooling_data.subgrid_temp = 0.f;

    /* No more cold gas! */
    p->cooling_data.subgrid_fcold = 0.f;
  }
}

/*
 * @brief Returns the subgrid temperature of a particle.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @return The subgrid temperature in internal units.
 */
float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp) {
  return p->cooling_data.subgrid_temp;
}

/*
 * @brief Returns the subgrid density of a particle.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @return The subgrid density in physical internal units.
 */
float cooling_get_subgrid_density(const struct part *p,
                                  const struct xpart *xp) {
  return p->cooling_data.subgrid_dens;
}

/**
 * @brief Compute the y-Compton contribution of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
double Cooling_get_ycompton(const struct phys_const *phys_const,
                            const struct hydro_props *hydro_props,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            const struct part *p, const struct xpart *xp) {

  return 0.;
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step,
 * since Grackle sub-cycles the cooling as needed.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The #cosmology.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 */
float cooling_timestep(const struct cooling_function_data *restrict cooling,
                       const struct phys_const *restrict phys_const,
                       const struct cosmology *restrict cosmo,
                       const struct unit_system *restrict us,
                       const struct hydro_props *hydro_props,
                       const struct part *restrict p,
                       const struct xpart *restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part *p, struct xpart *xp, double n) {

  xp->cooling_data.radiated_energy /= n;
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param phys_const The #phys_const.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_units(const struct unit_system *us,
                        const struct phys_const *phys_const,
                        struct cooling_function_data *cooling) {

  /* These are conversions from code units to cgs. */

  /* first cosmo:
     units for the expansion factor,
     use cosmological redshift!
     a_value is updated in cooling_update() */
  cooling->units.a_units = 1.0;
  if (cooling->redshift == -1) {
    cooling->units.a_value = 0.01;
  } else {
    cooling->units.a_value = 1.0;
  }

  /* We assume here all quantities to be in physical coordinate */
  cooling->units.comoving_coordinates = 0;

  /* CMB temperature */
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* then units */
  /* converts physical density to cgs number density for H */
  cooling->units.density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->units.length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  cooling->units.time_units = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->units.velocity_units =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  cooling->temp_to_u_factor =
      phys_const->const_boltzmann_k /
      (hydro_gamma_minus_one * phys_const->const_proton_mass *
       units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
  cooling->dudt_units =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* converts galaxy sSFR into G0 by scaling to MW values */
  const double time_to_yr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
                            (365.25f * 24.f * 60.f * 60.f);
  const double mass_to_solar_mass = 1.f / phys_const->const_solar_mass;
  const double length_to_pc =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) / 3.08567758e18f;

  cooling->time_to_Myr = time_to_yr * 1.e-6;

  /* G0 for MW=1.6 (Parravano etal 2003).  */
  /* Scaled to SFR density in solar neighborhood =0.002 Mo/Gyr/pc^3
     (J. Isern 2019) */
  cooling->G0_factor1 = 1.6f * mass_to_solar_mass /
                        (0.002f * time_to_yr * 1.e-9) /
                        (length_to_pc * length_to_pc * length_to_pc);

  /* Calibrated to sSFR for MW=2.71e-11 (Licquia etal 2015) */
  cooling->G0_factor2 = 1.6f / (2.71e-11f * time_to_yr);

  /* Calibrated to 0.015 SNe per year for MW (BC Reed 2005) */
  cooling->G0_factorSNe = 1.6f / 0.015;
}

/**
 * @brief Initialises Grackle.
 *
 * @param cooling The cooling properties to initialize
 */
void cooling_init_grackle(struct cooling_function_data *cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

  chemistry_data *chemistry = &cooling->chemistry;

  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters(chemistry) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry & cooling

  // Flag to activate the grackle machinery:
  chemistry->use_grackle = 2;  // 1=original grackle, 2=crackle
  /* Flag to include radiative cooling and actually update the thermal energy
   * during the chemistry solver. If off, the chemistry species will still be
   * updated. The most common reason to set this to off is to iterate the
   * chemistry network to an equilibrium state. */
  chemistry->with_radiative_cooling = 1;

  /* Flag to control which primordial chemistry network is used (set by Config
   * file) */
  chemistry->primordial_chemistry = COOLING_GRACKLE_MODE;

  /** Flag to enable H2 formation on dust grains, dust cooling, and dust-gas
   * heat transfer follow Omukai (2000). This assumes that the dust to gas ratio
   * scales with the metallicity. Default: 0. */
  chemistry->h2_on_dust = 0;

  /* Flag to enable metal cooling using the Cloudy tables. If enabled, the
   * cooling table to be used must be specified with the grackle_data_file
   * parameter. Default: 0. */
  chemistry->metal_cooling = cooling->with_metal_cooling;

  /* Flag to enable an effective CMB temperature floor. This is implemented by
   * subtracting the value of the cooling rate at TCMB from the total cooling
   * rate. Default: 1. */
  chemistry->cmb_temperature_floor = 1;

  /* Flag to enable a UV background. If enabled, the cooling table to be used
   * must be specified with the grackle_data_file parameter. Default: 0. */
  chemistry->UVbackground = cooling->with_uv_background;

  /* Path to the data file containing the metal cooling and UV background
   * tables: */
  chemistry->grackle_data_file = cooling->cloudy_table;

  /* The ratio of specific heats for an ideal gas. A direct calculation for the
   * molecular component is used if primordial_chemistry > 1. Default: 5/3. */
  chemistry->Gamma = hydro_gamma;

  /* Flag to control which three-body H2 formation rate is used. */
  chemistry->three_body_rate = 0;

  /* Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel
   * (2004). Default: 0. */
  chemistry->cie_cooling = 0;

  /* Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004).
   * Default: 0 */
  chemistry->h2_optical_depth_approximation = 0;
  /* Flag to enable a spatially uniform heating term approximating
   * photo-electric heating from dust from Tasker & Bryan (2008). Default: 0. */
  chemistry->photoelectric_heating = 0;
  /* photo-electric on [but not adjusted to local background, beware!] */
  chemistry->photoelectric_heating_rate = 8.5e-26;

  /* Flag to enable Compton heating from an X-ray background following Madau &
   * Efstathiou (1999). Default: 0. */
  chemistry->Compton_xray_heating = 0;

  /* Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
   * in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0. */
  chemistry->LWbackground_intensity = 0;

  /* Flag to enable suppression of Lyman-Werner flux due to Lyman-series
   * absorption
   *    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000).
   *    Default: 0. */
  chemistry->LWbackground_sawtooth_suppression = 0;

  /* volumetric heating rates is being provided in the volumetric_heating_rate
   * field of grackle_field_data */
  chemistry->use_volumetric_heating_rate =
      cooling->provide_volumetric_heating_rates;

  /* specific heating rates is being provided in the specific_heating_rate field
   * of grackle_field_data */
  chemistry->use_specific_heating_rate =
      cooling->provide_specific_heating_rates;

  /* Set parameters of temperature floor: 0=none, 1=provide scalar,
   * 2=provide array */
  chemistry->use_temperature_floor = 2;

  /* arrays of ionization and heating rates from radiative transfer solutions
   * are being provided */
  chemistry->use_radiative_transfer = 0;

  /* must be enabled to couple the passed radiative transfer fields to the
   * chemistry solver */
  chemistry->radiative_transfer_coupled_rate_solver = 0;

  /* enable intermediate stepping in applying radiative transfer fields to
   * chemistry solver. */
  chemistry->radiative_transfer_intermediate_step = 0;

  /* only use hydrogen ionization and heating rates from the radiative transfer
   * solutions. */
  chemistry->radiative_transfer_hydrogen_only = 0;

  /* Use Rahmati+13 self-shielding; 0=none, 1=HI only, 2=HI+HeI, 3=HI+HeI but
   * set HeII rates to 0 */
  chemistry->self_shielding_method = cooling->self_shielding_method;

  /* control behaviour of Grackle sub-step integrator */
  chemistry->max_iterations = cooling->max_step;
  chemistry->exit_after_iterations_exceeded = 0;
  chemistry->accuracy = cooling->timestep_accuracy;

  /* control behaviour of Grackle sub-step integration damping */
  if (cooling->grackle_damping_interval > 0) {
    chemistry->use_subcycle_timestep_damping = 1;
    chemistry->subcycle_timestep_damping_interval =
        cooling->grackle_damping_interval;
  } else {
    chemistry->use_subcycle_timestep_damping = 0;
    chemistry->subcycle_timestep_damping_interval = 0;
  }
  /* run on a single thread since Swift sends each particle to a single thread
   *chemistry->omp_nthreads = 1; */

  /* Turn on Li+ 2019 dust evolution model */
  chemistry->use_dust_evol = cooling->use_grackle_dust_evol;
  chemistry->use_dust_density_field = cooling->use_grackle_dust_evol;

  /* Load dust evolution parameters */
  if (cooling->use_grackle_dust_evol == 1) {
    chemistry->dust_destruction_eff = cooling->dust_destruction_eff;
    chemistry->sne_coeff = cooling->dust_sne_coeff;
    chemistry->sne_shockspeed = cooling->dust_sne_shockspeed;
    chemistry->dust_grainsize = cooling->dust_grainsize;
    chemistry->dust_growth_densref = cooling->dust_growth_densref;
    chemistry->dust_growth_tauref = cooling->dust_growth_tauref;
    /* Enable dust temperature calculation using ISRF */
    chemistry->metal_cooling = 1;
    chemistry->dust_chemistry = 1;
    chemistry->h2_on_dust = 1;
    chemistry->use_isrf_field = 1;
    chemistry->H2_self_shielding = 4;
    /* 2 means we specify the H2 shielding
       length ourselves (the gas smoothing length) */
    chemistry->H2_custom_shielding = 2;
    /* Solar abundances to pass to Grackle:
       He  (10.93 in units where log[H]=12, so photospheric mass fraction
          -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
       C   (8.43 -> 2.38e-3, AG=3.18e-3)
       N   (7.83 -> 0.70e-3, AG=1.15e-3)
       O   (8.69 -> 5.79e-3, AG=9.97e-3)
       Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
       Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
       Si  (7.51 -> 6.71e-4, AG=7.30e-4)
       S   (7.12 -> 3.12e-4, AG=3.80e-4)
       Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
       Fe (7.50 -> 1.31e-3, AG=1.92e-3)
       */
    chemistry->SolarAbundances[0] = 0.2485;
    chemistry->SolarAbundances[1] = 2.38e-3;
    chemistry->SolarAbundances[2] = 0.70e-3;
    chemistry->SolarAbundances[3] = 5.79e-3;
    chemistry->SolarAbundances[4] = 1.26e-3;
    chemistry->SolarAbundances[5] = 7.14e-4;
    chemistry->SolarAbundances[6] = 6.71e-3;
    chemistry->SolarAbundances[7] = 3.12e-4;
    chemistry->SolarAbundances[8] = 0.65e-4;
    chemistry->SolarAbundances[9] = 1.31e-3;
  } else {
    chemistry->use_dust_evol = 0;
  }

  cooling->use_grackle_h2_form =
      cooling->use_grackle_dust_evol && COOLING_GRACKLE_MODE >= 2;

  /* Initialize the chemistry object. */
  if (initialize_chemistry_data(&cooling->units) == 0) {
    error("Error in initialize_chemistry_data.");
  }
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling) {

  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");

  /* read parameters */
  cooling_read_parameters(parameter_file, cooling, phys_const, us);

  /* Set up the units system. */
  cooling_init_units(us, phys_const, cooling);

  /* Set up grackle */
  cooling_init_grackle(cooling);
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data *cooling) {
  //_free_chemistry_data(&cooling->chemistry, &grackle_rates);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * Nothing to do beyond writing the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream) {
  restart_write_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                       stream, "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Nothing to do beyond reading the structure from the stream.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo) {
  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  /* Set up grackle */
  cooling_init_grackle(cooling);
}
