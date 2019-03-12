/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_HYDRO_PROPERTIES
#define SWIFT_HYDRO_PROPERTIES

/**
 * @file hydro_properties.h
 * @brief Contains all the constants and parameters of the hydro scheme
 */

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local includes. */
#include "restart.h"

/* Forward declarations */
struct cosmology;
struct swift_params;
struct gravity_props;
struct phys_const;
struct unit_system;

/**
 * @brief Contains all the constants and parameters of the hydro scheme
 */
struct hydro_props {

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal smoothing length */
  float h_max;

  /*! Minimal smoothing length expressed as ratio to softening length */
  float h_min_ratio;

  /*! Minimal smoothing length */
  float h_min;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Time integration properties */
  float CFL_condition;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /*! Minimal temperature allowed */
  float minimal_temperature;

  /*! Minimal physical internal energy per unit mass */
  float minimal_internal_energy;

  /*! Initial temperature */
  float initial_temperature;

  /*! Initial physical internal energy per unit mass */
  float initial_internal_energy;

  /*! Primordial hydrogen mass fraction for initial energy conversion */
  float hydrogen_mass_fraction;

  /*! Temperature of the neutral to ionized transition of Hydrogen */
  float hydrogen_ionization_temperature;

  /*! Mean molecular weight below hydrogen ionization temperature */
  float mu_neutral;

  /*! Mean molecular weight above hydrogen ionization temperature */
  float mu_ionised;

  /*! Artificial viscosity parameters */
  struct {
    /*! For the fixed, simple case. Also used to set the initial AV
        coefficient for variable schemes. */
    float alpha;

    /*! Artificial viscosity (max) for the variable case (e.g. M&M) */
    float alpha_max;

    /*! Artificial viscosity (min) for the variable case (e.g. M&M) */
    float alpha_min;

    /*! The decay length of the artificial viscosity (used in M&M, etc.) */
    float length;
  } viscosity;

  /*! Thermal diffusion parameters */
  struct {

    /*! Initialisation value, or the case for constant thermal diffusion coeffs
     */
    float alpha;

    /*! Tuning parameter for speed of ramp up/down */
    float beta;

    /*! Maximal value for alpha_diff */
    float alpha_max;

    /*! Minimal value for alpha_diff */
    float alpha_min;

  } diffusion;
};

void hydro_props_print(const struct hydro_props *p);
void hydro_props_init(struct hydro_props *p,
                      const struct phys_const *phys_const,
                      const struct unit_system *us,
                      struct swift_params *params);

void hydro_props_update(struct hydro_props *p, const struct gravity_props *gp,
                        const struct cosmology *cosmo);

#if defined(HAVE_HDF5)
void hydro_props_print_snapshot(hid_t h_grpsph, const struct hydro_props *p);
#endif

/* Dump/restore. */
void hydro_props_struct_dump(const struct hydro_props *p, FILE *stream);
void hydro_props_struct_restore(const struct hydro_props *p, FILE *stream);

/* Setup for tests */
void hydro_props_init_no_hydro(struct hydro_props *p);

#endif /* SWIFT_HYDRO_PROPERTIES */
