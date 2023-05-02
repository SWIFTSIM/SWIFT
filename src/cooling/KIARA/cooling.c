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

#include "../config.h"

/* Include header */
#include "cooling.h"

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <math.h>

/* The grackle library itself */
#include <grackle.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling_io.h"
#include "entropy_floor.h"
#include "star_formation.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* define heating and cooling limits on thermal energy, per timestep */
#define GRACKLE_HEATLIM 1000.0
#define GRACKLE_COOLLIM 0.01

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
void cooling_update(const struct cosmology* cosmo,
		    struct pressure_floor_props *pressure_floor_props,
                    struct cooling_function_data* cooling, struct space* s) {
  /* set current time */
  if (cooling->redshift == -1)
    cooling->units.a_value = cosmo->a;
  else
    cooling->units.a_value = 1. / (1. + cooling->redshift);
}

/**
 * @brief Print the chemical network
 *
 * @param xp The #xpart to print
 */
void cooling_print_fractions(const struct xpart* restrict xp) {

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
void cooling_first_init_part(const struct phys_const* restrict phys_const,
                             const struct unit_system* restrict us,
                             const struct hydro_props* hydro_props,
                             const struct cosmology* restrict cosmo,
                             const struct cooling_function_data* cooling,
                             const struct part* restrict p,
                             struct xpart* restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
  xp->cooling_data.time_last_event = -cooling->thermal_time;

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;
  zero = 0.;

  /* primordial chemistry >= 1: Start with everything neutral (as in dark ages)
   */
  xp->cooling_data.HI_frac = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  xp->cooling_data.HII_frac = zero;
  xp->cooling_data.HeI_frac = p->chemistry_data.metal_mass_fraction[chemistry_element_He];
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
  xp->cooling_data.dust_mass = zero;
  for (int i=0; i<chemistry_element_count; i++) xp->cooling_data.dust_mass_fraction[i] = zero;
  xp->cooling_data.dust_temperature = zero;
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
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
float cooling_get_radiated_energy(const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print_backend(const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle', mode = %i", cooling->chemistry.primordial_chemistry);
  message("CloudyTable = %s", cooling->cloudy_table);
  message("Redshift = %g", cooling->redshift);
  message("UV background flag = %d", cooling->with_uv_background);
  message("Metal cooling flag = %i", cooling->chemistry.metal_cooling);
  message("Self Shielding flag = %i", cooling->self_shielding_method);
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
void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, gr_float rho,
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

  /* HeII */
  species_densities[5] = xp->cooling_data.e_frac * rho;
  data->e_density = &species_densities[5];
}
#else
void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, gr_float rho,
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
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, 
                              const struct cooling_function_data* restrict cooling,
			      gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
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
      species_densities[20] = xp->cooling_data.dust_mass / p->mass * species_densities[12];
      data->dust_density = &species_densities[20];
      species_densities[21] = p->feedback_data.SNe_ThisTimeStep;
      data->SNe_ThisTimeStep = &species_densities[21];
      species_densities[22] = p->chemistry_data.smoothed_G0;
      if (p->id%10000==0 && p->chemistry_data.G0>0) message("%lld  z=%g  rho=%g  G0=%g  SNe=%g\n",p->id,1./cooling->units.a_value-1.,species_densities[12]*cooling->units.density_units/1.673e-24,p->chemistry_data.G0, p->feedback_data.SNe_ThisTimeStep);
      data->isrf_habing = &species_densities[22];
      for (int i=0; i<chemistry_element_count; i++) {
          species_densities[23+i] = p->chemistry_data.metal_mass_fraction[i] * species_densities[12];
          species_densities[23+chemistry_element_count+i] = xp->cooling_data.dust_mass_fraction[i] * species_densities[12];
	  //if (species_densities[23+i]<0 || species_densities[23+chemistry_element_count+i]<0) message("NEGATIVE MASS FRACTION! %lld %d %g %g %g %g\n",p->id, i, p->chemistry_data.metal_mass_fraction[i], species_densities[23+i], species_densities[23+chemistry_element_count+i], xp->cooling_data.dust_mass_fraction[i]);
      }

      /* Load gas metallicities */
      data->He_gas_metalDensity = &species_densities[24];
      data->C_gas_metalDensity = &species_densities[25];
      data->N_gas_metalDensity = &species_densities[26];
      data->O_gas_metalDensity = &species_densities[27];
      data->Ne_gas_metalDensity = &species_densities[28];
      data->Mg_gas_metalDensity = &species_densities[29];
      data->Si_gas_metalDensity = &species_densities[30];
      data->S_gas_metalDensity = &species_densities[31];  
      data->Ca_gas_metalDensity = &species_densities[32];
      data->Fe_gas_metalDensity = &species_densities[33];
      /* Load dust metallicities */
      data->He_dust_metalDensity = &species_densities[35];
      data->C_dust_metalDensity = &species_densities[36];
      data->N_dust_metalDensity = &species_densities[37];
      data->O_dust_metalDensity = &species_densities[38];
      data->Ne_dust_metalDensity = &species_densities[39];
      data->Mg_dust_metalDensity = &species_densities[40];
      data->Si_dust_metalDensity = &species_densities[41];
      data->S_dust_metalDensity = &species_densities[42];
      data->Ca_dust_metalDensity = &species_densities[43];
      data->Fe_dust_metalDensity = &species_densities[44];
  }
  else {
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
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, 
                              const struct cooling_function_data* restrict cooling,
			      gr_float rho,
                              gr_float species_densities[N_SPECIES]) {
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
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, gr_float rho,
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
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, gr_float rho,
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
void cooling_copy_from_grackle1(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {

  /* HI */
  xp->cooling_data.HI_frac = *data->HI_density / rho;

  /* HII */
  xp->cooling_data.HII_frac = *data->HII_density / rho;

  /* HeI */
  xp->cooling_data.HeI_frac = *data->HeI_density / rho;

  /* HeII */
  xp->cooling_data.HeII_frac = *data->HeII_density / rho;

  /* HeIII */
  xp->cooling_data.HeIII_frac = *data->HeIII_density / rho;

  /* e */
  xp->cooling_data.e_frac = *data->e_density / rho;
}
#else
void cooling_copy_from_grackle1(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {}
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
void cooling_copy_from_grackle2(grackle_field_data* data, struct part* p,
                                struct xpart* xp, 
                                const struct cooling_function_data* restrict cooling,
			        gr_float rho) {
  /* HM */
  xp->cooling_data.HM_frac = *data->HM_density / rho;
  /* H2I */
  xp->cooling_data.H2I_frac = *data->H2I_density / rho;
  /* H2II */
  xp->cooling_data.H2II_frac = *data->H2II_density / rho;

  /* Dust model */
  if (cooling->use_grackle_dust_evol == 1) {
      xp->cooling_data.dust_mass = *data->dust_density * p->mass / rho;
      /* Load gas metallicities */
      p->chemistry_data.metal_mass_fraction[1] = *data->He_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[2] = *data->C_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[3] = *data->N_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[4] = *data->O_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[5] = *data->Ne_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[6] = *data->Mg_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[7] = *data->Si_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[8] = *data->S_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[9] = *data->Ca_gas_metalDensity / rho;
      p->chemistry_data.metal_mass_fraction[10] = *data->Fe_gas_metalDensity / rho;
      /* Load dust metallicities */
      xp->cooling_data.dust_mass_fraction[1] = *data->He_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[2] = *data->C_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[3] = *data->N_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[4] = *data->O_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[5] = *data->Ne_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[6] = *data->Mg_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[7] = *data->Si_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[8] = *data->S_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[9] = *data->Ca_dust_metalDensity / rho;
      xp->cooling_data.dust_mass_fraction[10] = *data->Fe_dust_metalDensity / rho;
      //message("DUST: m %g  rho %g  m %g  C %g  Fe %g\n",p->mass, rho, xp->cooling_data.dust_mass, xp->cooling_data.dust_mass_fraction[1], xp->cooling_data.dust_mass_fraction[8]);
  }
}
#else
void cooling_copy_from_grackle2(grackle_field_data* data, struct part* p,
                                struct xpart* xp, 
                                const struct cooling_function_data* restrict cooling,
			        gr_float rho) {}
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
void cooling_copy_from_grackle3(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {

  /* DI */
  xp->cooling_data.DI_frac = *data->DI_density / rho;

  /* DII */
  xp->cooling_data.DII_frac = *data->DII_density / rho;

  /* HDI */
  xp->cooling_data.HDI_frac = *data->HDI_density / rho;
}
#else
void cooling_copy_from_grackle3(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho) {}
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
void cooling_copy_to_grackle(grackle_field_data* data,
    			     const struct unit_system* restrict us,
                             const struct cosmology* restrict cosmo,
                             const struct cooling_function_data* restrict cooling,
                             const struct part* p, const struct xpart* xp,
			     const double dt, const double u_floor,
			     gr_float species_densities[N_SPECIES]) {

  int i;
  /* set values */
  /* grid */
  data->grid_dx = 0.f;
  data->grid_rank = GRACKLE_RANK;

  data->grid_dimension = malloc(GRACKLE_RANK * sizeof(int));
  data->grid_start = malloc(GRACKLE_RANK * sizeof(int));
  data->grid_end = malloc(GRACKLE_RANK * sizeof(int));
  for (i = 0; i < 3; i++) {
      data->grid_dimension[i] = 1; // the active dimension not including ghost zones.
      data->grid_start[i] = 0;
      data->grid_end[i] = 0;
  }
  data->grid_dimension[0] = GRACKLE_NPART;
  data->grid_end[0] = GRACKLE_NPART - 1;

  /* specific_heating_rate has to be in cgs units; no unit conversion done within grackle */
  species_densities[15] = hydro_get_physical_internal_energy_dt(p, cosmo) * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS) / units_cgs_conversion_factor(us, UNIT_CONV_TIME); 
  data->specific_heating_rate = &species_densities[15];

  /* get particle density, internal energy in physical coordinates (still code units) */
  double T_subgrid = cooling_get_subgrid_temperature(p, xp); 
  if ( p->cooling_data.subgrid_temp == 0. ) {  // normal cooling mode
      species_densities[12] = hydro_get_physical_density(p, cosmo);
      species_densities[13] = hydro_get_physical_internal_energy(p, xp, cosmo);
      /* volumetric_heating_rate stores the minimum thermal energy for this particle */
      species_densities[14] = u_floor;
  }
  else {  // subgrid ISM model
      species_densities[12] = cooling_get_subgrid_density(p, xp); // physical subgrid density
      species_densities[13] = cooling_convert_temp_to_u(T_subgrid, xp->cooling_data.e_frac, cooling, p); // physical internal energy
      float T_cmb = 2.73*(1.f+cosmo->z);
      species_densities[14] = cooling_convert_temp_to_u(T_cmb, 0.f, cooling, p); // CMB temp is floor
      /* If tracking H2, turn off specific heating rate in ISM because it ruins H2 fraction */
      species_densities[15] = 0.;
  }
  /* load into grackle structure */
  data->density = &species_densities[12];
  data->internal_energy = &species_densities[13];
  data->volumetric_heating_rate = &species_densities[14]; /* This is actually the minimum thermal energy for this particle */

  /* velocity (maybe not needed?) */
  species_densities[16] = p->v_full[0] * cosmo->a_inv;
  species_densities[17] = p->v_full[1] * cosmo->a_inv;
  species_densities[18] = p->v_full[2] * cosmo->a_inv;
  data->x_velocity = &species_densities[16];
  data->y_velocity = &species_densities[17];
  data->z_velocity = &species_densities[18];

  cooling_copy_to_grackle1(data, p, xp, species_densities[12], species_densities);
  cooling_copy_to_grackle2(data, p, xp, cooling, species_densities[12], species_densities);
  cooling_copy_to_grackle3(data, p, xp, species_densities[12], species_densities);

  data->RT_heating_rate = NULL;
  data->RT_HI_ionization_rate = NULL;
  data->RT_HeI_ionization_rate = NULL;
  data->RT_HeII_ionization_rate = NULL;
  data->RT_H2_dissociation_rate = NULL;

  species_densities[19] = chemistry_get_total_metal_mass_fraction_for_cooling(p) * species_densities[12];
  data->metal_density = &species_densities[19];
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
void cooling_copy_from_grackle(grackle_field_data* data, struct part* p,
                               struct xpart* xp, 
			       const struct cooling_function_data* restrict cooling,
			       gr_float rho) {

  cooling_copy_from_grackle1(data, p, xp, rho);
  cooling_copy_from_grackle2(data, p, xp, cooling, rho);
  cooling_copy_from_grackle3(data, p, xp, rho);
}

/**
 * @brief free memory associated with grackle driver
 *
 * @param data The grackle_field_data structure from grackle.
 */
void cooling_grackle_free_data(grackle_field_data* data) {

  free(data->grid_dimension);
  free(data->grid_start);
  free(data->grid_end);
  //free(data->metal_density);
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
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, double dt,
    double u_floor, int mode) {

  /* set current units for conversion to physical quantities */
  code_units units = cooling->units;

  /* initialize data to send to grackle */
  gr_float *species_densities;
  species_densities = (gr_float *)calloc(N_SPECIES, sizeof(gr_float));
  grackle_field_data data;

  /* load particle information from particle to grackle data */
  cooling_copy_to_grackle(&data, us, cosmo, cooling, p, xp, dt, u_floor, species_densities);

  /* Run Grackle in desired mode */
  gr_float return_value = 0.f;
  switch (mode) {
    case 0:
      if( *data.dust_density> *data.metal_density ) if (cooling_get_subgrid_temperature(p, xp)>0) message("SUBGRID: %lld before nH=%g  u=%g  T=%g  fH2=%g  Mdust=%g  Tdust=%g DTM=%g %g\n",p->id, species_densities[12]*cooling->units.density_units/1.673e-24, species_densities[13], cooling_get_subgrid_temperature(p, xp), xp->cooling_data.H2I_frac+xp->cooling_data.H2II_frac, xp->cooling_data.dust_mass, xp->cooling_data.dust_temperature,  *data.dust_density, *data.metal_density);
      /* solve chemistry, advance thermal energy by dt */
      if (solve_chemistry(&units, &data, dt) == 0) {
        error("Error in Grackle solve_chemistry.");
      }
      /* copy from grackle data to particle */
      cooling_copy_from_grackle(&data, p, xp, cooling, species_densities[12]);
      return_value = data.internal_energy[0];
#if COOLING_GRACKLE_MODE >= 2
      /* Compute dust temperature */
      double t_dust = xp->cooling_data.dust_temperature;
      if (calculate_dust_temperature(&units, &data, &t_dust) == 0) {
        error("Error in Grackle calculate dust temperature.");
      }
      xp->cooling_data.dust_temperature = t_dust;
      /* Reset accumulated local variables to zero */
      p->feedback_data.SNe_ThisTimeStep = 0.f;
      if( *data.dust_density > *data.metal_density ) if (cooling_get_subgrid_temperature(p, xp)>0) message("SUBGRID: %lld after nH=%g  u=%g  T=%g  fH2=%g  Mdust=%g  Tdust=%g DTM=%g %g\n",p->id, species_densities[12]*cooling->units.density_units/1.673e-24, return_value, cooling_get_subgrid_temperature(p, xp), xp->cooling_data.H2I_frac+xp->cooling_data.H2II_frac, xp->cooling_data.dust_mass, xp->cooling_data.dust_temperature, *data.dust_density, *data.metal_density);
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
 * @brief Compute the cooling time
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 *
 * @return cooling time
 */
gr_float cooling_time(const struct phys_const* restrict phys_const,
                      const struct unit_system* restrict us,
                      const struct hydro_props* hydro_properties,
                      const struct cosmology* restrict cosmo,
                      const struct cooling_function_data* restrict cooling,
                      const struct part* restrict p,
                      struct xpart* restrict xp) {

  struct part p_temp = *p;  // gets rid of const in declaration
  gr_float cooling_time = cooling_grackle_driver(
      phys_const, us, cosmo, hydro_properties, cooling, &p_temp, xp, 0., 0., 1);
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
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_properties,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  struct part p_temp = *p;  // gets rid of const in declaration
  struct xpart xp_temp = *xp;  // gets rid of const in declaration
  float temperature = cooling_grackle_driver(
      phys_const, us, cosmo, hydro_properties, cooling, &p_temp, &xp_temp, 0., 0., 2);
  /* const float mu = 4. / (1. + 3. * hydro_properties->hydrogen_mass_fraction);  // testing, for neutral gas only
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo); // * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  const float temperature = hydro_gamma_minus_one * mu * u * phys_const->const_proton_mass/phys_const->const_boltzmann_k; */
  return temperature;
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
void cooling_cool_part(const struct phys_const* restrict phys_const,
                       const struct unit_system* restrict us,
                       const struct cosmology* restrict cosmo,
                       const struct hydro_props* hydro_props,
                       const struct entropy_floor_properties* floor_props,
		       const struct pressure_floor_props *pressure_floor_props,
                       const struct cooling_function_data* restrict cooling,
                       struct part* restrict p, struct xpart* restrict xp,
                       const double dt, const double dt_therm,
                       const double time) {

  /* No cooling if particle is decoupled */
  if (p->feedback_data.decoupling_delay_time > 0.f
        || p->feedback_data.cooling_shutoff_delay_time > 0.f) {
    return;
  }

  /* Update the subgrid properties */
  cooling_set_particle_subgrid_properties( phys_const, us, 
	  cosmo, hydro_props, floor_props, cooling, p, xp);

  /* No cooling happens over zero time */
  if (dt == 0.f || dt_therm == 0.f ) {
    return;
  }

  /* Current energy */
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
  //const float T_old = cooling_get_temperature( phys_const, hydro_props, us, cosmo, cooling, p, xp); // for debugging only

  /* Compute the entropy floor */
  const double T_floor = entropy_floor_temperature(p, cosmo, floor_props);
  const double u_floor = cooling_convert_temp_to_u(T_floor, xp->cooling_data.e_frac, cooling, p);

  /* Do grackle cooling, if it's been more than thermal_time since last cooling */
  gr_float u_new = u_old;
  u_new = cooling_grackle_driver(phys_const, us, cosmo, hydro_props, cooling,
                                   p, xp, dt_therm, u_floor, 0);

  /* Apply simulation-wide minimum temperature */
  u_new = max(u_new, hydro_props->minimal_internal_energy);

  /* Assign new thermal energy to particle */
  float cool_du_dt = 0.;

  if (p->cooling_data.subgrid_temp == 0.) {  
    /* Normal cooling; check that we are not going to go below any of the limits */
    if (u_new > GRACKLE_HEATLIM * u_old) u_new = GRACKLE_HEATLIM * u_old;
    if (u_new < GRACKLE_COOLLIM * u_old) u_new = GRACKLE_COOLLIM * u_old;
    u_new = max(u_new, u_floor);

    /* Calculate the cooling rate */
    cool_du_dt = (u_new - u_old) / dt_therm;

    /* Update the internal energy time derivative */
    hydro_set_physical_internal_energy_dt(p, cosmo, cool_du_dt);
  }
  else {
    /* Particle is in subgrid mode; result is stored in subgrid_temp */
    p->cooling_data.subgrid_temp = cooling_convert_u_to_temp(u_new, xp->cooling_data.e_frac, cooling, p);
//    if( p->id%1000==0) message("SUBGRID: copy %lld u_new=%g  u_floor=%g  T=%g\n",p->id, u_new, u_floor, cooling_get_subgrid_temperature(p, xp));

    /* Set internal energy time derivative to 0 for overall particle */
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);

    /* Force the overall particle to lie on the equation of state */
    hydro_set_physical_internal_energy(p, xp, cosmo, u_floor);
  }

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cool_du_dt * dt_therm;

  /* Record this cooling event */
  xp->cooling_data.time_last_event = time;

  /* set subgrid properties for use in SF routine */
  cooling_set_particle_subgrid_properties(
      phys_const, us, cosmo, hydro_props, floor_props, cooling, p, xp);
}

/**
 * @brief Set the subgrid properties (rho, T) of the gas particle for use in SF routine
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
    const struct cooling_function_data *cooling, 
    struct part *p, struct xpart *xp) {

  /* Get temperature of overall particle */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float temperature = cooling_convert_u_to_temp(u, xp->cooling_data.e_frac, cooling, p);

  /* Get H number density */
  const double rho_com = hydro_get_comoving_density(p);
  const double rho = hydro_get_physical_density(p, cosmo);

  /* Subgrid model is on if particle is in the Jeans EOS regime */
  const float T_floor = entropy_floor_Jeans_temperature( rho, rho_com, cosmo, floor_props);
  if (T_floor > 0) {
    /* YES: If first time in subgrid, set temperature to particle T, otherwise limit to particle T */
    if (p->cooling_data.subgrid_temp == 0. ) p->cooling_data.subgrid_temp = temperature;
    /* Subgrid temperature should be no higher than overall particle temperature */
    else p->cooling_data.subgrid_temp = min(p->cooling_data.subgrid_temp, temperature);

    /* We set the subgrid density based on pressure equilibrium with overall particle */
    p->cooling_data.subgrid_dens = rho * temperature / p->cooling_data.subgrid_temp;

//    if( p->id%1000==0) message("SUBGRID: set %lld, rho=%g  Tfloor=%g  u=%g  Told=%g  T=%g  dudt=%g\n",p->id, rho*cooling->units.density_units/1.673e-24, T_floor, u, temperature, p->cooling_data.subgrid_temp,hydro_get_physical_internal_energy_dt(p, cosmo));

  }
  else {
    /* NO: subgrid density is the actual particle's physical density */
    p->cooling_data.subgrid_dens = rho;
    
    /* set subgrid temperature to 0 indicating it's not in subgrid mode */
    p->cooling_data.subgrid_temp = 0.;
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
double Cooling_get_ycompton(const struct phys_const* phys_const,
                            const struct hydro_props* hydro_props,
                            const struct unit_system* us,
                            const struct cosmology* cosmo,
                            const struct cooling_function_data* cooling,
                            const struct part* p, const struct xpart* xp) {

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
float cooling_timestep(const struct cooling_function_data* restrict cooling,
                       const struct phys_const* restrict phys_const,
                       const struct cosmology* restrict cosmo,
                       const struct unit_system* restrict us,
                       const struct hydro_props* hydro_props,
                       const struct part* restrict p,
                       const struct xpart* restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Split the coolong content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part* p, struct xpart* xp, double n) {

  xp->cooling_data.radiated_energy /= n;
}

/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param phys_const The #phys_const.
 * @param cooling The cooling properties to initialize
 */
void cooling_init_units(const struct unit_system* us,
                        const struct phys_const* phys_const,
                        struct cooling_function_data* cooling) {

  /* These are conversions from code units to cgs. */

  /* first cosmo */
  cooling->units.a_units = 1.0;  // units for the expansion factor
  if (cooling->redshift == -1) {  // use cosmological redshift
      cooling->units.a_value = 0.01;  // arbitrary; gets reset in cooling_update()
  }
  else cooling->units.a_value = 1.0;  

  /* We assume here all quantities to be in physical coordinate */
  cooling->units.comoving_coordinates = 0;

  /* then units */
  /* converts physical density to cgs number density for H */
  cooling->units.density_units =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  cooling->units.length_units =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  cooling->units.time_units = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  cooling->units.velocity_units =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  cooling->temp_to_u_factor = phys_const->const_boltzmann_k / (hydro_gamma_minus_one * phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
}

/**
 * @brief Initialises Grackle.
 *
 * @param cooling The cooling properties to initialize
 */
void cooling_init_grackle(struct cooling_function_data* cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

  chemistry_data* chemistry = &cooling->chemistry;

  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters(chemistry) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry & cooling

  // Flag to activate the grackle machinery:
  chemistry->use_grackle = 1;  // grackle on (duh)
  // Flag to include radiative cooling and actually update the thermal energy
  // during the chemistry solver. If off, the chemistry species will still be
  // updated. The most common reason to set this to off is to iterate the
  // chemistry network to an equilibrium state. Default: 1.
  chemistry->with_radiative_cooling = 1; // cooling on
  // Flag to control which primordial chemistry network is used (set by Config
  // file)
  chemistry->primordial_chemistry = COOLING_GRACKLE_MODE;
  // Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat
  // transfer follow Omukai (2000). This assumes that the dust to gas ratio
  // scales with the metallicity. Default: 0.
  chemistry->h2_on_dust = 0;  // dust cooling/chemistry on
  // Flag to enable metal cooling using the Cloudy tables. If enabled, the
  // cooling table to be used must be specified with the grackle_data_file
  // parameter. Default: 0.
  chemistry->metal_cooling = cooling->with_metal_cooling;  // metal cooling on
  // Flag to enable an effective CMB temperature floor. This is implemented by
  // subtracting the value of the cooling rate at TCMB from the total cooling
  // rate. Default: 1.
  chemistry->cmb_temperature_floor = 1;
  // Flag to enable a UV background. If enabled, the cooling table to be used
  // must be specified with the grackle_data_file parameter. Default: 0.
  chemistry->UVbackground = cooling->with_uv_background;
  // Path to the data file containing the metal cooling and UV background
  // tables:
  chemistry->grackle_data_file = cooling->cloudy_table;  // data file
  // The ratio of specific heats for an ideal gas. A direct calculation for the
  // molecular component is used if primordial_chemistry > 1. Default: 5/3.
  chemistry->Gamma = hydro_gamma;  // our eos set in Config.sh
  // Flag to control which three-body H2 formation rate is used.
  chemistry->three_body_rate = 0;
  // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel
  // (2004). Default: 0.
  chemistry->cie_cooling = 0;
  // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004).
  // Default: 0
  chemistry->h2_optical_depth_approximation = 0;
  // Flag to enable a spatially uniform heating term approximating
  // photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
  chemistry->photoelectric_heating =
      0;  // photo-electric on [but not adjusted to local background, beware!]
  chemistry->photoelectric_heating_rate = 8.5e-26;
  // Flag to enable Compton heating from an X-ray background following Madau &
  // Efstathiou (1999). Default: 0.
  chemistry->Compton_xray_heating = 0;
  // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
  // in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
  chemistry->LWbackground_intensity = 0;
  // Flag to enable suppression of Lyman-Werner flux due to Lyman-series
  // absorption
  //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000).
  //    Default: 0.
  chemistry->LWbackground_sawtooth_suppression = 0;
  // volumetric heating rates is being provided in the volumetric_heating_rate
  // field of grackle_field_data
  chemistry->use_volumetric_heating_rate = cooling->provide_volumetric_heating_rates;
  // specific heating rates is being provided in the specific_heating_rate field
  // of grackle_field_data
  chemistry->use_specific_heating_rate = cooling->provide_specific_heating_rates;
  // arrays of ionization and heating rates from radiative transfer solutions
  // are being provided
  chemistry->use_radiative_transfer = 0;
  // must be enabled to couple the passed radiative transfer fields to the
  // chemistry solver
  chemistry->radiative_transfer_coupled_rate_solver = 0;
  // enable intermediate stepping in applying radiative transfer fields to
  // chemistry solver.
  chemistry->radiative_transfer_intermediate_step = 0;
  // only use hydrogen ionization and heating rates from the radiative transfer
  // solutions.
  chemistry->radiative_transfer_hydrogen_only = 0;
  // Use Rahmati+13 self-shielding; 0=none, 1=HI only, 2=HI+HeI, 3=HI+HeI but
  // set HeII rates to 0
  chemistry->self_shielding_method = cooling->self_shielding_method;
  // control behaviour of Grackle sub-step integrator
  //chemistry->max_iterations = cooling->max_step;
  //chemistry->exit_after_iterations_exceeded = 0;
  // run on a single thread since Swift sends each particle to a single thread
  //chemistry->omp_nthreads = 1;

  // Turn on Li+ 2019 dust evolution model
  chemistry->use_dust_evol = cooling->use_grackle_dust_evol;

  // Load dust evolution parameters
  if (cooling->use_grackle_dust_evol == 1) {
    chemistry->dust_destruction_eff = cooling->dust_destruction_eff;
    chemistry->sne_coeff = cooling->dust_sne_coeff;
    chemistry->sne_shockspeed = cooling->dust_sne_shockspeed;
    chemistry->dust_grainsize = cooling->dust_grainsize;
    chemistry->dust_growth_densref = cooling->dust_growth_densref;
    chemistry->dust_growth_tauref = cooling->dust_growth_tauref;
    // Enable dust temperature calculation using ISRF
    chemistry->metal_cooling = 1;
    chemistry->dust_chemistry = 1;
    chemistry->h2_on_dust = 1;
    chemistry->use_isrf_field = 1;
    chemistry->H2_self_shielding = 3;
  } else {
    chemistry->use_dust_evol = 0;
  }

  cooling->use_grackle_h2_form = cooling->use_grackle_dust_evol && COOLING_GRACKLE_MODE >= 2;

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
void cooling_init_backend(struct swift_params* parameter_file,
                          const struct unit_system* us,
                          const struct phys_const* phys_const,
                          const struct hydro_props* hydro_props,
                          struct cooling_function_data* cooling) {

  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");

  /* read parameters */
  cooling_read_parameters(parameter_file, cooling, phys_const);

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
void cooling_clean(struct cooling_function_data* cooling) {
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
void cooling_struct_dump(const struct cooling_function_data* cooling,
                         FILE* stream) {
  restart_write_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
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
void cooling_struct_restore(struct cooling_function_data* cooling, FILE* stream,
                            const struct cosmology* cosmo) {
  restart_read_blocks((void*)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  /* Set up grackle */
  cooling_init_grackle(cooling);
}
