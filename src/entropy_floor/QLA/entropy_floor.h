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
 ******************************************************************************/
#ifndef SWIFT_ENTROPY_FLOOR_QLA_H
#define SWIFT_ENTROPY_FLOOR_QLA_H

#include "adiabatic_index.h"
#include "cosmology.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "parser.h"
#include "units.h"

/**
 * @file src/entropy_floor/QLA/entropy_floor.h
 * @brief Entropy floor used in the quick Lyman-alpha model
 */

/**
 * @brief Properties of the entropy floor in the QLA model.
 */
struct entropy_floor_properties {

  /*! Density threshold for the floor in Hydrogen atoms per cubic cm */
  float density_threshold_H_p_cm3;

  /*! Density threshold for the floor in internal units */
  float density_threshold;

  /*! Inverse of the density threshold for the floor in internal units */
  float density_threshold_inv;

  /*! Over-density threshold for the floor */
  float over_density_threshold;

  /*! Temperature of the floor at the density threshold in Kelvin */
  float temperature_norm_K;

  /*! Temperature of the floor at the density thresh. in internal units */
  float temperature_norm;

  /*! Pressure of the floor at the density thresh. in internal units */
  float pressure_norm;
};

/**
 * @brief Compute the entropy floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
static INLINE float entropy_floor(
    const struct part *p, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props) {

  /* Comoving density in internal units */
  const float rho_com = hydro_get_comoving_density(p);

  /* Physical density in internal units */
  const float rho_phys = hydro_get_physical_density(p, cosmo);

  /* Mean baryon density in co-moving internal units for over-density condition
   * (Recall cosmo->critical_density_0 is 0 in a non-cosmological run,
   * making the over-density condition a no-op) */
  const float rho_crit_0 = cosmo->critical_density_0;
  const float rho_crit_baryon = cosmo->Omega_b * rho_crit_0;

  /* Physical pressure */
  float pressure = 0.f;

  /* Are we in the regime of the equation of state? */
  if ((rho_com >= rho_crit_baryon * props->over_density_threshold) &&
      (rho_phys >= props->density_threshold)) {

    const float pressure_floor =
        props->pressure_norm * rho_phys * props->density_threshold_inv;

    pressure = max(pressure, pressure_floor);
  }

  /* Convert to an entropy.
   * (Recall that the entropy is the same in co-moving and phycial frames) */
  return gas_entropy_from_pressure(rho_phys, pressure);
}

/**
 * @brief Compute the temperature from the entropy floor for a given #part
 *
 * Calculate the EoS temperature, the particle is not updated.
 * This is the temperature exactly corresponding to the imposed EoS shape.
 * It only matches the entropy returned by the entropy_floor() function
 * for a neutral gas with primoridal abundance.
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
static INLINE float entropy_floor_temperature(
    const struct part *p, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props) {

  /* Comoving density in internal units */
  const float rho_com = hydro_get_comoving_density(p);

  /* Physical density in internal units */
  const float rho_phys = hydro_get_physical_density(p, cosmo);

  /* Mean baryon density in co-moving internal units for over-density condition
   * (Recall cosmo->critical_density_0 is 0 in a non-cosmological run,
   * making the over-density condition a no-op) */
  const float rho_crit_0 = cosmo->critical_density_0;
  const float rho_crit_baryon = cosmo->Omega_b * rho_crit_0;

  /* Physical */
  float temperature = 0.f;

  /* Are we in the regime of the equation of state? */
  if ((rho_com >= rho_crit_baryon * props->over_density_threshold) &&
      (rho_phys >= props->density_threshold)) {

    const float temperature_floor = props->temperature_norm;
    temperature = max(temperature, temperature_floor);
  }

  return temperature;
}

/**
 * @brief Initialise the entropy floor by reading the parameters and converting
 * to internal units.
 *
 * The input temperatures and number densities are converted to entropy and
 * density assuming a neutral gas of primoridal abundance.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_const The physical constants.
 * @param hydro_props The propoerties of the hydro scheme.
 * @param props The entropy floor properties to fill.
 */
static INLINE void entropy_floor_init(struct entropy_floor_properties *props,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct hydro_props *hydro_props,
                                      struct swift_params *params) {

  /* Read the parameters in the units they are set */
  props->density_threshold_H_p_cm3 = parser_get_param_float(
      params, "QLAEntropyFloor:density_threshold_H_p_cm3");
  props->over_density_threshold =
      parser_get_param_float(params, "QLAEntropyFloor:over_density_threshold");
  props->temperature_norm_K =
      parser_get_param_float(params, "QLAEntropyFloor:temperature_norm_K");

  /* Initial Hydrogen abundance (mass fraction) */
  const double X_H = hydro_props->hydrogen_mass_fraction;

  /* Now convert to internal units assuming primodial Hydrogen abundance */
  props->temperature_norm =
      props->temperature_norm_K /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  props->density_threshold =
      props->density_threshold_H_p_cm3 /
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY) *
      phys_const->const_proton_mass / X_H;

  /* We assume neutral gas */
  const float mean_molecular_weight = hydro_props->mu_neutral;

  /* Get the common terms */
  props->density_threshold_inv = 1.f / props->density_threshold;

  /* P_norm = (k_B * T) / (m_p * mu) * rho_threshold */
  props->pressure_norm =
      ((phys_const->const_boltzmann_k * props->temperature_norm) /
       (phys_const->const_proton_mass * mean_molecular_weight)) *
      props->density_threshold;
}

/**
 * @brief Print the properties of the entropy floor to stdout.
 *
 * @param props The entropy floor properties.
 */
static INLINE void entropy_floor_print(
    const struct entropy_floor_properties *props) {

  message("Entropy floor is 'Quick Lyman-alpha' with:");
  message("Floor with slope n=%.3f at rho=%e (%e H/cm^3) and T=%.1f K", 1.f,
          props->density_threshold, props->density_threshold_H_p_cm3,
          props->temperature_norm);
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of entropy floor to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void entropy_floor_write_flavour(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Entropy floor", "Quick Lyman-alpha");
}
#endif

/**
 * @brief Write an entropy floor struct to the given FILE as a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void entropy_floor_struct_dump(
    const struct entropy_floor_properties *props, FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct entropy_floor_properties),
                       1, stream, "entropy floor", "entropy floor properties");
}

/**
 * @brief Restore a entropy floor struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void entropy_floor_struct_restore(
    struct entropy_floor_properties *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct entropy_floor_properties), 1,
                      stream, NULL, "entropy floor properties");
}

#endif /* SWIFT_ENTROPY_FLOOR_QLA_H */
