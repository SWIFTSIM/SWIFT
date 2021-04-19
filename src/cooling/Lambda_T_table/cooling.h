/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_COOLING_LAMBDA_T_TABLE_H
#define SWIFT_COOLING_LAMBDA_T_TABLE_H

/**
 * @file src/cooling/Lambda_T_table/cooling.h
 * @brief A cooling function that is interpolated on Lambda(T) data from an
 * ASCII table.
 */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "cooling_properties.h"
#include "cosmology.h"
#include "entropy_floor.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "part.h"

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
INLINE static void cooling_update(const struct cosmology* cosmo,
                                  struct cooling_function_data* cooling,
                                  struct space* s) {
  // Add content if required.
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current time (since the Big Bang or start of the run) in
 * internal units.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, const float dt,
    const float dt_therm, const double time) {

  const float rho = p->rho;
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float P_old = gas_pressure_from_internal_energy(rho, u_old);
  const float T_old = P_old * phys_const->const_proton_mass /
                      (rho * phys_const->const_boltzmann_k);

  const double logT = log10(T_old);
  int jl = 0;
  int ju = cooling->nT;
  while (ju - jl > 1) {
    int jm = (ju + jl) >> 1;
    if (logT > cooling->logT[jm]) {
      jl = jm;
    } else {
      ju = jm;
    }
  }
  if (jl == cooling->nT - 1) {
    --jl;
  }

  const double logTm = cooling->logT[jl];
  const double logTp = cooling->logT[jl + 1];
  const double logLambdam = cooling->logLambda[jl];
  const double logLambdap = cooling->logLambda[jl + 1];
  const double fT = (logT - logTm) / (logTp - logTm);
  float cooling_du_dt = -pow(10., logLambdam * (1. - fT) + logLambdap * fT) *
                        rho; /* todo: need to divide by mu^2 */

  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  if (u_old + hydro_du_dt * dt < 0.0f) {
    cooling_du_dt = 0.0f;
  } else if (u_old + (hydro_du_dt + cooling_du_dt) * dt < 0.0f) {
    cooling_du_dt = (u_old + dt * hydro_du_dt) / dt;
  }

  hydro_set_physical_internal_energy_dt(p, cosmo, hydro_du_dt + cooling_du_dt);
}

/**
 * @brief Computes the cooling time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended data of the particle.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props, const struct part* restrict p,
    const struct xpart* restrict xp) {

  const float rho = p->rho;
  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float P_old = gas_pressure_from_internal_energy(rho, u_old);
  const float T_old = P_old * phys_const->const_proton_mass /
                      (rho * phys_const->const_boltzmann_k);

  const double logT = log10(T_old);
  int jl = 0;
  int ju = cooling->nT;
  while (ju - jl > 1) {
    int jm = (ju + jl) >> 1;
    if (logT > cooling->logT[jm]) {
      jl = jm;
    } else {
      ju = jm;
    }
  }
  if (jl == cooling->nT - 1) {
    --jl;
  }

  const double logTm = cooling->logT[jl];
  const double logTp = cooling->logT[jl + 1];
  const double logLambdam = cooling->logLambda[jl];
  const double logLambdap = cooling->logLambda[jl + 1];
  const double logdLambdadTm = cooling->logdLambdadT[jl];
  const double logdLambdadTp = cooling->logdLambdadT[jl + 1];
  const double fT = (logT - logTm) / (logTp - logTm);
  const float cooling_du_dt =
      pow(10., logLambdam * (1. - fT) + logLambdap * fT) *
      rho; /* todo: need to divide by mu^2 */
  const float cooling_d2u_dt2 =
      pow(10., logdLambdadTm * (1. - fT) + logdLambdadTp * fT) *
      rho; /* todo: need to divide by mu^2 */

  float dt = FLT_MAX;
  if (cooling_du_dt > 0.0f) {
    const float cooling_dt = u_old / cooling_du_dt;
    dt = fminf(dt, 0.1f * cooling_dt);
  }
  if (cooling_d2u_dt2 > 0.0f) {
    const float cooling_dt = sqrtf(u_old / cooling_d2u_dt2);
    dt = fminf(dt, 10. * cooling_dt);
  }
  return dt;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param data The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* data, const struct part* restrict p,
    struct xpart* restrict xp) {}

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
INLINE static float cooling_get_temperature(
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
 * @param Returns the subgrid temperature of a particle.
 *
 * This model has no subgrid quantity. We return -1.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_temperature(const struct part* p,
                                                    const struct xpart* xp) {
  return -1.f;
}

/**
 * @param Returns the subgrid density of a particle.
 *
 * This model has no subgrid quantity. We return -1.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_density(const struct part* p,
                                                const struct xpart* xp) {
  return -1.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return 0.f;
}

/**
 * @brief Split the cooling content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
static INLINE void cooling_split_part(struct part* p, struct xpart* xp,
                                      double n) {}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The cooling properties to initialize
 */
static INLINE void cooling_init_backend(struct swift_params* parameter_file,
                                        const struct unit_system* us,
                                        const struct phys_const* phys_const,
                                        const struct hydro_props* hydro_props,
                                        struct cooling_function_data* cooling) {

  char file_name[1000];
  parser_get_param_string(parameter_file, "LambdaTTableCooling:file_name",
                          file_name);

  cooling->logT = (double*)malloc(100 * sizeof(double));
  cooling->logLambda = (double*)malloc(100 * sizeof(double));
  cooling->nT = 100;

  FILE* file = fopen(file_name, "r");
  int nT = 0;
  while (fscanf(file, "%lg\t%lg", &cooling->logT[nT],
                &cooling->logLambda[nT]) == 2) {
    ++nT;
    if (cooling->nT == nT) {
      cooling->nT <<= 1;
      cooling->logT =
          (double*)realloc(cooling->logT, cooling->nT * sizeof(double));
      cooling->logLambda =
          (double*)realloc(cooling->logLambda, cooling->nT * sizeof(double));
    }
  }
  cooling->nT = nT;
  cooling->logT = (double*)realloc(cooling->logT, cooling->nT * sizeof(double));
  cooling->logLambda =
      (double*)realloc(cooling->logLambda, cooling->nT * sizeof(double));

  float basefacs[5] = {1, 5, -3, 0, 0};
  double mp = phys_const->const_proton_mass;
  double factor = units_general_cgs_conversion_factor(us, basefacs) * (mp * mp);
  message("factor: %g", factor);
  for (int i = 0; i < nT; ++i) {
    cooling->logT[i] /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
    cooling->logLambda[i] /= factor;
  }

  cooling->logdLambdadT = (double*)malloc(cooling->nT * sizeof(double));
  cooling->logdLambdadT[0] = (cooling->logLambda[1] - cooling->logLambda[0]) /
                             (cooling->logT[1] - cooling->logT[0]);
  for (int i = 1; i < nT - 1; ++i) {
    cooling->logdLambdadT[i] =
        (cooling->logLambda[i + 1] - cooling->logLambda[i - 1]) /
        (cooling->logT[i + 1] - cooling->logT[i - 1]);
  }
  cooling->logdLambdadT[nT - 1] =
      (cooling->logLambda[nT - 1] - cooling->logLambda[nT - 2]) /
      (cooling->logT[nT - 1] - cooling->logT[nT - 2]);

  for (int i = 0; i < nT; ++i) {
    if (cooling->logT[i] > 0.) {
      cooling->logT[i] = log10(cooling->logT[i]);
    } else {
      cooling->logT[i] = -99.;
    }
    if (cooling->logLambda[i] > 0.) {
      cooling->logLambda[i] = log10(cooling->logLambda[i]);
    } else {
      cooling->logLambda[i] = -99.;
    }
    cooling->logdLambdadT[i] = abs(cooling->logdLambdadT[i]);
    if (cooling->logdLambdadT[i] > 0.) {
      cooling->logdLambdadT[i] = log10(cooling->logdLambdadT[i]);
    } else {
      cooling->logdLambdadT[i] = -99.;
    }
  }

  fclose(file);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Lambda-T table'.");
  /*  for (int i = 0; i < cooling->nT; ++i) {
      message("%g\t%g", cooling->logT[i], cooling->logLambda[i]);
    }*/
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
static INLINE void cooling_clean(struct cooling_function_data* cooling) {

  free(cooling->logT);
  free(cooling->logLambda);
  free(cooling->logdLambdadT);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
static INLINE void cooling_struct_dump(
    const struct cooling_function_data* cooling, FILE* stream) {}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Empty structure so nothing to do here.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
static INLINE void cooling_struct_restore(struct cooling_function_data* cooling,
                                          FILE* stream,
                                          const struct cosmology* cosmo) {}

#endif /* SWIFT_COOLING_LAMBDA_T_TABLE_H */
