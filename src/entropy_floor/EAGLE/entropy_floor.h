/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ENTROPY_FLOOR_EAGLE_H
#define SWIFT_ENTROPY_FLOOR_EAGLE_H

/* Code config */
#include <config.h>

/* System include */
#include <stdio.h>

/* Pre-declarations */
struct cosmology;
struct part;
struct phys_const;
struct unit_system;
struct hydro_props;
struct swift_params;

/**
 * @file src/entropy_floor/EAGLE/entropy_floor.h
 * @brief Entropy floor used in the EAGLE model
 */

/**
 * @brief Properties of the entropy floor in the EAGLE model.
 */
struct entropy_floor_properties {

  /*! Density threshold for the Jeans floor in Hydrogen atoms per cubic cm */
  float Jeans_density_threshold_H_p_cm3;

  /*! Density threshold for the Jeans floor in internal units */
  float Jeans_density_threshold;

  /*! Inverse of the density threshold for the Jeans floor in internal units */
  float Jeans_density_threshold_inv;

  /*! Over-density threshold for the Jeans floor */
  float Jeans_over_density_threshold;

  /*! Slope of the Jeans floor power-law */
  float Jeans_gamma_effective;

  /*! Temperature of the Jeans floor at the density threshold in Kelvin */
  float Jeans_temperature_norm_K;

  /*! Temperature of the Jeans floor at the density thresh. in internal units */
  float Jeans_temperature_norm;

  /*! Pressure of the Jeans floor at the density thresh. in internal units */
  float Jeans_pressure_norm;

  /*! Density threshold for the Cool floor in Hydrogen atoms per cubic cm */
  float Cool_density_threshold_H_p_cm3;

  /*! Density threshold for the Cool floor in internal units */
  float Cool_density_threshold;

  /*! Inverse of the density threshold for the Cool floor in internal units */
  float Cool_density_threshold_inv;

  /*! Over-density threshold for the Cool floor */
  float Cool_over_density_threshold;

  /*! Slope of the Cool floor power-law */
  float Cool_gamma_effective;

  /*! Temperature of the Cool floor at the density threshold in Kelvin */
  float Cool_temperature_norm_K;

  /*! Temperature of the Cool floor at the density thresh. in internal units */
  float Cool_temperature_norm;

  /*! Pressure of the Cool floor at the density thresh. in internal units */
  float Cool_pressure_norm;
};

float entropy_floor_gas_pressure(const float rho_phys, const float rho_com,
                                 const struct cosmology *cosmo,
                                 const struct entropy_floor_properties *props);

float entropy_floor(const struct part *p, const struct cosmology *cosmo,
                    const struct entropy_floor_properties *props);

float entropy_floor_gas_temperature(
    const float rho_phys, const float rho_com, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props);

float entropy_floor_temperature(const struct part *p,
                                const struct cosmology *cosmo,
                                const struct entropy_floor_properties *props);

void entropy_floor_init(struct entropy_floor_properties *props,
                        const struct phys_const *phys_const,
                        const struct unit_system *us,
                        const struct hydro_props *hydro_props,
                        struct swift_params *params);

void entropy_floor_print(const struct entropy_floor_properties *props);

#ifdef HAVE_HDF5

void entropy_floor_write_flavour(hid_t h_grp);
#endif

void entropy_floor_struct_dump(const struct entropy_floor_properties *props,
                               FILE *stream);

void entropy_floor_struct_restore(struct entropy_floor_properties *props,
                                  FILE *stream);

#endif /* SWIFT_ENTROPY_FLOOR_EAGLE_H */
