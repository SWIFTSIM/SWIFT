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
#ifndef SWIFT_COOLING_KIARA_H
#define SWIFT_COOLING_KIARA_H

/**
 * @file src/cooling/KIARA/cooling.h
 * @brief Cooling using the GRACKLE 3.1.1 library.
 */

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <math.h>

/* The grackle library itself */
#include <grackle.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling_properties.h"
#include "entropy_floor.h"
#include "error.h"
#include "fof.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3
#if COOLING_GRACKLE_MODE >= 2
#define N_SPECIES 46 /* This further includes properties for dust model */
#else
#define N_SPECIES                                 \
  21 /* This includes extra values at end to hold \
        rho,u,dudt,vx,vy,vz,u_floor,mZ,dummyvar */
#endif

/* define heating and cooling limits on thermal energy, per timestep */
#define GRACKLE_HEATLIM 1000.f
#define GRACKLE_COOLLIM 0.01f
#define MAX_COLD_ISM_FRACTION 0.9f
/* Minimum particle column length as a fraction of p->h.
 * Should be <= mean interparticle spacing.
 * For 48 Ngb mean spacing ~ 0.887
 * For 57 Ngb mean spacing ~ 0.837
 * For 114 Ngb mean spacing ~ 0.664
 *
 * Basically a limiter on rho / |grad rho| for a depth
 * for shielding, and how bad you think |grad rho|
 * really is at estimating the depth. */
#define MIN_SHIELD_H_FRAC 0.332f

void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time);

void cooling_print_fractions(const struct xpart *restrict xp);
void cooling_first_init_part(const struct phys_const *restrict phys_const,
                             const struct unit_system *restrict us,
                             const struct hydro_props *hydro_properties,
                             const struct cosmology *restrict cosmo,
                             const struct cooling_function_data *cooling,
                             struct part *restrict p,
                             struct xpart *restrict xp);
void cooling_post_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp);

void cooling_print_backend(const struct cooling_function_data *cooling);

void cooling_copy_to_grackle1(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]);
void cooling_copy_to_grackle2(
    grackle_field_data *data, const struct part *p, const struct xpart *xp,
    const struct cooling_function_data *restrict cooling, const double dt,
    gr_float rho, gr_float species_densities[N_SPECIES]);
void cooling_copy_to_grackle3(grackle_field_data *data, const struct part *p,
                              const struct xpart *xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]);
void cooling_copy_from_grackle1(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho);
void cooling_copy_from_grackle2(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho);
void cooling_copy_from_grackle3(grackle_field_data *data, const struct part *p,
                                struct xpart *xp, gr_float rho);
void cooling_copy_to_grackle(
    grackle_field_data *data, const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling, const struct part *p,
    const struct xpart *xp, const double dt, const double T_floor,
    gr_float species_densities[N_SPECIES], gr_float *iact_rates, int mode);
void cooling_copy_from_grackle(
    grackle_field_data *data, struct part *p, struct xpart *xp,
    const struct cooling_function_data *restrict cooling, gr_float rho);
void cooling_grackle_free_data(grackle_field_data *data);
gr_float cooling_grackle_driver(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_properties,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    double dt, double T_floor, int mode);
void cooling_do_grackle_cooling(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, gr_float *iact_rates,
    const double dt, const double dt_therm);
gr_float cooling_time(const struct phys_const *restrict phys_const,
                      const struct unit_system *restrict us,
                      const struct hydro_props *hydro_properties,
                      const struct cosmology *restrict cosmo,
                      const struct cooling_function_data *restrict cooling,
                      const struct part *restrict p, struct xpart *restrict xp,
                      const float rhocool, const float ucool);

float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *hydro_properties,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp);

void firehose_cooling_and_dust(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct hydro_props *restrict hydro_props,
    const struct cooling_function_data *restrict cooling,
    struct part *restrict p, struct xpart *restrict xp, const double dt);

void cooling_cool_part(const struct phys_const *restrict phys_const,
                       const struct unit_system *restrict us,
                       const struct cosmology *restrict cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct pressure_floor_props *pressure_floor_props,
                       const struct cooling_function_data *restrict cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const double dt, const double dt_therm,
                       const double time);

void cooling_set_particle_subgrid_properties(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp);

float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp);

float cooling_get_subgrid_density(const struct part *p, const struct xpart *xp);

float cooling_get_radiated_energy(const struct xpart *restrict xp);

double cooling_get_ycompton(const struct phys_const *phys_const,
                            const struct hydro_props *hydro_props,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            const struct part *p, const struct xpart *xp);

float cooling_timestep(const struct cooling_function_data *restrict cooling,
                       const struct phys_const *restrict phys_const,
                       const struct cosmology *restrict cosmo,
                       const struct unit_system *restrict us,
                       const struct hydro_props *hydro_properties,
                       const struct part *restrict p,
                       const struct xpart *restrict xp);

void cooling_split_part(struct part *p, struct xpart *xp, double n);

void cooling_init_units(const struct unit_system *us,
                        const struct phys_const *phys_const,
                        struct cooling_function_data *cooling);
void cooling_init_grackle(struct cooling_function_data *cooling);

void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling);

void cooling_clean(struct cooling_function_data *cooling);
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream);
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo);
void cooling_sputter_dust(const struct unit_system *restrict us,
                          const struct cosmology *restrict cosmo,
                          const struct cooling_function_data *restrict cooling,
                          struct part *restrict p, struct xpart *restrict xp,
                          const double dt);

/**
 * @brief Compute the electron pressure of a #part based on the cooling
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
INLINE static double cooling_get_electron_pressure(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {
  return 0;
}

/**
 * @brief Compute the specific thermal energy (physical) for a given
 * temperature.
 *
 * Converts T to u (internal physical units) for a given particle.
 *
 * @param temperature Particle temperature in K
 * @param ne Electron number density relative to H atom density
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 */
INLINE static double cooling_convert_temp_to_u(
    const double temperature, const double ne,
    const struct cooling_function_data *cooling, const struct part *p) {

  const float X_H =
      chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
  const float yhelium = (1. - X_H) / (4. * X_H);
  const float mu = (1. + yhelium) / (1. + ne + 4. * yhelium);

  return temperature * mu * cooling->temp_to_u_factor;
}

/**
 * @brief Compute the temperature for a given physical specific energy
 *
 * Converts T to u (internal physical units) for a given particle.
 *
 * @param u Physical specific energy
 * @param ne Electron number density relative to H atom density
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 */
INLINE static double cooling_convert_u_to_temp(
    const double u, const double ne,
    const struct cooling_function_data *cooling, const struct part *p) {

  const float X_H =
      chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
  const float yhelium = (1. - X_H) / (4. * X_H);
  const float mu = (1. + yhelium) / (1. + ne + 4. * yhelium);

  return u / (mu * cooling->temp_to_u_factor);
}

/**
 * @brief Compute the cold ISM fraction at a given factor above subgrid
 * threshold density
 *
 * Compute the cold ISM fraction at a given factor above subgrid threshold
 * density. This uses a fit to the density vs. cold gas fraction relation from
 * Springel+Hernquist 2003.
 *
 * @param dens_fac Density factor above threshold density
 * @param cooling #cooling_function_data struct.
 */
INLINE static double cooling_compute_cold_ISM_fraction(
    const double n_H, const struct cooling_function_data *cooling) {

  float fc = cooling->cold_ISM_frac;
  float dens_fac = n_H * cooling->subgrid_threshold_n_H_inv;
  if (dens_fac > 1.) {
    fc = cooling->cold_ISM_frac +
         (1. - cooling->cold_ISM_frac) * (1. - exp(-log10(dens_fac)));
    fc = fmin(fc, MAX_COLD_ISM_FRACTION);
  }
  return fc;
}

/**
 * @brief Compute the subgrid density based on pressure equilibrium in a 2-phase
 * ISM model
 *
 * We set the subgrid density based on pressure equilibrium with overall
 * particle. The pressure is set by 1-cold_ISM_frac of the mass in the warm
 * phase.
 *
 * @param rho SPH (non-subgrid) physical particle density.
 * @param n_H SPH (non-subgrid) physical particle H number density.
 * @param temp SPH (non-subgrid) particle temperature.
 * @param subgrid_temp Subgrid particle temperature.
 * @param cooling #cooling_function_data struct.
 */
INLINE static double cooling_compute_subgrid_density(
    const double rho, const double n_H, const double temp,
    const double subgrid_temp, const struct cooling_function_data *cooling) {

  const double ism_frac = cooling_compute_cold_ISM_fraction(
      n_H * cooling->subgrid_threshold_n_H_inv, cooling);
  double subgrid_dens =
      (1.f - ism_frac) * rho * temp / (ism_frac * subgrid_temp);

  /* Cap at max value which should be something vaguely like GMC densities */
  subgrid_dens = fmin(subgrid_dens, cooling->max_subgrid_density);
  return subgrid_dens;
}

/**
 * @brief Return warm ISM temperature if above SF threshold density, otherwise
 * 0.
 *
 * @param p Pointer to the particle data.
 * @param cooling The properties of the cooling function.
 * @param us The unit system.
 * @param phys_const The physical constant in internal units.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static float warm_ISM_temperature(
    const struct part *restrict p, const struct cooling_function_data *cooling,
    const struct phys_const *phys_const, const struct cosmology *cosmo) {

  float temperature = 0.f;

  /* Mean baryon density in co-moving internal units for over-density condition
   * (Recall cosmo->critical_density_0 is 0 in a non-cosmological run,
   * making the over-density condition a no-op) */
  const float rho_crit_0 = cosmo->critical_density_0;
  const float rho_crit_baryon = cosmo->Omega_b * rho_crit_0;
  const double rho_com = hydro_get_comoving_density(p);
  const double rho_phys = hydro_get_physical_density(p, cosmo);

  if (rho_com >= rho_crit_baryon * 100.f) {
    const double n_H = rho_phys * 0.75 / phys_const->const_proton_mass;
    if (n_H * cooling->subgrid_threshold_n_H_inv > 1.f) {
      const float EOS_slope = cooling->subgrid_warm_ism_EOS;
      temperature = cooling->subgrid_threshold_T *
                    pow(n_H * cooling->subgrid_threshold_n_H_inv, EOS_slope);
    }
  }
  return temperature;
}

/**
 * @brief Copy all grackle fields into a new set.  THIS IS ONLY USED FOR
 * DEBUGGING.
 *
 * @param my_fields The target (new) set of grackle particle properties.
 * @param old_fields The original (old) set of grackle particle properties.
 * @param field_size Number of particles to copy.
 */
INLINE static void cooling_copy_grackle_fields(grackle_field_data *my_fields,
                                               grackle_field_data *old_fields,
                                               int field_size) {
  int i;

  for (i = 0; i < field_size; i++) {

    printf("loop copy_grackle_fields %g %p\n", old_fields->density[0],
           my_fields->density);

    my_fields->density[i] = old_fields->density[i];
    my_fields->HI_density[i] = old_fields->HI_density[i];
    my_fields->HII_density[i] = old_fields->HII_density[i];
    my_fields->HM_density[i] = old_fields->HM_density[i];
    my_fields->HeI_density[i] = old_fields->HeI_density[i];
    my_fields->HeII_density[i] = old_fields->HeII_density[i];
    my_fields->HeIII_density[i] = old_fields->HeIII_density[i];
    my_fields->H2I_density[i] = old_fields->H2I_density[i];
    my_fields->H2II_density[i] = old_fields->H2II_density[i];
    my_fields->DI_density[i] = old_fields->DI_density[i];
    my_fields->DII_density[i] = old_fields->DII_density[i];
    my_fields->HDI_density[i] = old_fields->HDI_density[i];
    my_fields->e_density[i] = old_fields->e_density[i];
    // solar metallicity
    my_fields->metal_density[i] = old_fields->metal_density[i];

    my_fields->x_velocity[i] = 0.0;
    my_fields->y_velocity[i] = 0.0;
    my_fields->z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    my_fields->internal_energy[i] = old_fields->internal_energy[i];

    my_fields->volumetric_heating_rate[i] =
        old_fields->volumetric_heating_rate[i];
    my_fields->specific_heating_rate[i] = old_fields->specific_heating_rate[i];
    my_fields->temperature_floor[i] = old_fields->temperature_floor[i];

    my_fields->isrf_habing[i] = old_fields->isrf_habing[i];
    my_fields->RT_HI_ionization_rate[i] = old_fields->RT_HI_ionization_rate[i];
    my_fields->RT_HeI_ionization_rate[i] =
        old_fields->RT_HeI_ionization_rate[i];
    my_fields->RT_HeII_ionization_rate[i] =
        old_fields->RT_HeII_ionization_rate[i];
    my_fields->RT_H2_dissociation_rate[i] =
        old_fields->RT_H2_dissociation_rate[i];
    my_fields->RT_heating_rate[i] = old_fields->RT_heating_rate[i];

    if (grackle_data->use_dust_evol) {
      my_fields->dust_density[i] = old_fields->dust_density[i];
      my_fields->He_gas_metalDensity[i] = old_fields->He_gas_metalDensity[i];
      my_fields->C_gas_metalDensity[i] = old_fields->C_gas_metalDensity[i];
      my_fields->N_gas_metalDensity[i] = old_fields->N_gas_metalDensity[i];
      my_fields->O_gas_metalDensity[i] = old_fields->O_gas_metalDensity[i];
      my_fields->Ne_gas_metalDensity[i] = old_fields->Ne_gas_metalDensity[i];
      my_fields->Mg_gas_metalDensity[i] = old_fields->Mg_gas_metalDensity[i];
      my_fields->Si_gas_metalDensity[i] = old_fields->Si_gas_metalDensity[i];
      my_fields->S_gas_metalDensity[i] = old_fields->S_gas_metalDensity[i];
      my_fields->Ca_gas_metalDensity[i] = old_fields->Ca_gas_metalDensity[i];
      my_fields->Fe_gas_metalDensity[i] = old_fields->Fe_gas_metalDensity[i];
      my_fields->He_dust_metalDensity[i] = old_fields->He_dust_metalDensity[i];
      my_fields->C_dust_metalDensity[i] = old_fields->C_dust_metalDensity[i];
      my_fields->N_dust_metalDensity[i] = old_fields->N_dust_metalDensity[i];
      my_fields->O_dust_metalDensity[i] = old_fields->O_dust_metalDensity[i];
      my_fields->Ne_dust_metalDensity[i] = old_fields->Ne_dust_metalDensity[i];
      my_fields->Mg_dust_metalDensity[i] = old_fields->Mg_dust_metalDensity[i];
      my_fields->Si_dust_metalDensity[i] = old_fields->Si_dust_metalDensity[i];
      my_fields->S_dust_metalDensity[i] = old_fields->S_dust_metalDensity[i];
      my_fields->Ca_dust_metalDensity[i] = old_fields->Ca_dust_metalDensity[i];
      my_fields->Fe_dust_metalDensity[i] = old_fields->Fe_dust_metalDensity[i];
      my_fields->SNe_ThisTimeStep[i] = old_fields->SNe_ThisTimeStep[i];
    }
  }
  printf("done copy_grackle_fields\n");

  return;
}

/**
 * @brief Allocate a new set of grackle fields in memory. THIS IS ONLY USED FOR
 * DEBUGGING.
 *
 * @param my_fields The target (new) set of grackle particle properties.
 * @param field_size Number of particles to copy.
 * @param dust_flag Are we using grackle's dust model (Jones, Smith, Dave 2024)?
 */
INLINE static void cooling_grackle_malloc_fields(grackle_field_data *my_fields,
                                                 int field_size,
                                                 int dust_flag) {
  my_fields->density = malloc(field_size * sizeof(gr_float));
  my_fields->internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields->x_velocity = malloc(field_size * sizeof(gr_float));
  my_fields->y_velocity = malloc(field_size * sizeof(gr_float));
  my_fields->z_velocity = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 1
  my_fields->HI_density = malloc(field_size * sizeof(gr_float));
  my_fields->HII_density = malloc(field_size * sizeof(gr_float));
  my_fields->HeI_density = malloc(field_size * sizeof(gr_float));
  my_fields->HeII_density = malloc(field_size * sizeof(gr_float));
  my_fields->HeIII_density = malloc(field_size * sizeof(gr_float));
  my_fields->e_density = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 2
  my_fields->HM_density = malloc(field_size * sizeof(gr_float));
  my_fields->H2I_density = malloc(field_size * sizeof(gr_float));
  my_fields->H2II_density = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 3
  my_fields->DI_density = malloc(field_size * sizeof(gr_float));
  my_fields->DII_density = malloc(field_size * sizeof(gr_float));
  my_fields->HDI_density = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  my_fields->metal_density = malloc(field_size * sizeof(gr_float));

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields->volumetric_heating_rate = malloc(field_size * sizeof(gr_float));
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields->specific_heating_rate = malloc(field_size * sizeof(gr_float));
  my_fields->temperature_floor = malloc(field_size * sizeof(gr_float));

  // radiative transfer ionization / dissociation rate fields (provide in units
  // [1/s])
  my_fields->RT_HI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_HeI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
  // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  my_fields->RT_heating_rate = malloc(field_size * sizeof(gr_float));

  // H2 model
  my_fields->H2_self_shielding_length = malloc(field_size * sizeof(gr_float));
  my_fields->H2_custom_shielding_factor = malloc(field_size * sizeof(gr_float));
  my_fields->isrf_habing = malloc(field_size * sizeof(gr_float));

  if (dust_flag) {
    my_fields->dust_density = malloc(field_size * sizeof(gr_float));
    my_fields->He_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->C_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->N_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->O_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Ne_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Mg_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Si_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->S_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Ca_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Fe_gas_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->He_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->C_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->N_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->O_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Ne_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Mg_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Si_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->S_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Ca_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->Fe_dust_metalDensity = malloc(field_size * sizeof(gr_float));
    my_fields->SNe_ThisTimeStep = malloc(field_size * sizeof(gr_float));
  }
  return;
}

/**
 * @brief Free a set of grackle fields from memory. THIS IS ONLY USED FOR
 * DEBUGGING.
 *
 * @param my_fields The target (new) set of grackle particle properties.
 * @param dust_flag Are we using grackle's dust model (Jones, Smith, Dave 2024)?
 */
INLINE static void cooling_grackle_free_fields(grackle_field_data *my_fields,
                                               int dust_flag) {
  free(my_fields->grid_dimension);
  free(my_fields->grid_start);
  free(my_fields->grid_end);

  free(my_fields->density);
  free(my_fields->internal_energy);
  free(my_fields->x_velocity);
  free(my_fields->y_velocity);
  free(my_fields->z_velocity);
  // for primordial_chemistry >= 1
  free(my_fields->HI_density);
  free(my_fields->HII_density);
  free(my_fields->HeI_density);
  free(my_fields->HeII_density);
  free(my_fields->HeIII_density);
  free(my_fields->e_density);
  // for primordial_chemistry >= 2
  free(my_fields->HM_density);
  free(my_fields->H2I_density);
  free(my_fields->H2II_density);
  // for primordial_chemistry >= 3
  free(my_fields->DI_density);
  free(my_fields->DII_density);
  free(my_fields->HDI_density);
  // for metal_cooling = 1
  free(my_fields->metal_density);

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  free(my_fields->volumetric_heating_rate);
  // specific heating rate (provide in units [egs s^-1 g^-1]
  free(my_fields->specific_heating_rate);
  free(my_fields->temperature_floor);

  // radiative transfer ionization / dissociation rate fields (provide in units
  // [1/s])
  free(my_fields->RT_HI_ionization_rate);
  free(my_fields->RT_HeI_ionization_rate);
  free(my_fields->RT_HeII_ionization_rate);
  free(my_fields->RT_H2_dissociation_rate);
  // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  free(my_fields->RT_heating_rate);

  // H2 model
  free(my_fields->H2_self_shielding_length);
  free(my_fields->H2_custom_shielding_factor);
  free(my_fields->isrf_habing);

  if (dust_flag) {
    free(my_fields->dust_density);
    free(my_fields->He_gas_metalDensity);
    free(my_fields->C_gas_metalDensity);
    free(my_fields->N_gas_metalDensity);
    free(my_fields->O_gas_metalDensity);
    free(my_fields->Ne_gas_metalDensity);
    free(my_fields->Mg_gas_metalDensity);
    free(my_fields->Si_gas_metalDensity);
    free(my_fields->S_gas_metalDensity);
    free(my_fields->Ca_gas_metalDensity);
    free(my_fields->Fe_gas_metalDensity);
    free(my_fields->He_dust_metalDensity);
    free(my_fields->C_dust_metalDensity);
    free(my_fields->N_dust_metalDensity);
    free(my_fields->O_dust_metalDensity);
    free(my_fields->Ne_dust_metalDensity);
    free(my_fields->Mg_dust_metalDensity);
    free(my_fields->Si_dust_metalDensity);
    free(my_fields->S_dust_metalDensity);
    free(my_fields->Ca_dust_metalDensity);
    free(my_fields->Fe_dust_metalDensity);
    free(my_fields->SNe_ThisTimeStep);
  }
  return;
}
#endif /* SWIFT_COOLING_KIARA_H */
