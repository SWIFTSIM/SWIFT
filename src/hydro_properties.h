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
#include "hydro_parameters.h"
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

  /* ------ Smoothing lengths parameters ---------- */

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal smoothing length (internal units) */
  float h_max;

  /*! Minimal smoothing length expressed as ratio to softening length */
  float h_min_ratio;

  /*! Minimal smoothing length (internal units) */
  float h_min;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /* ------ Neighbour number definition ------------ */

  /*! Are we using the mass-weighted definition of neighbour number? */
  int use_mass_weighted_num_ngb;

  /* ------ Time integration parameters ------------ */

  /*! Time integration properties */
  float CFL_condition;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /* ------ Temperature parameters ----------------- */

  /*! Minimal temperature allowed */
  float minimal_temperature;

  /*! Minimal physical internal energy per unit mass (internal units) */
  float minimal_internal_energy;

  /*! Initial temperature */
  float initial_temperature;

  /*! Initial physical internal energy per unit mass (internal units) */
  float initial_internal_energy;

  /*! Primordial hydrogen mass fraction for initial energy conversion */
  float hydrogen_mass_fraction;

  /*! Temperature of the neutral to ionized transition of Hydrogen */
  float hydrogen_ionization_temperature;

  /*! Mean molecular weight below hydrogen ionization temperature */
  float mu_neutral;

  /*! Mean molecular weight above hydrogen ionization temperature */
  float mu_ionised;

  /* ------ Particle splitting parameters ---------- */

  /*! Is particle splitting activated? */
  int particle_splitting;

  /*! Mass above which particles get split (internal units) */
  float particle_splitting_mass_threshold;

  /*! Are we generating random IDs when splitting particles? */
  int generate_random_ids;

  /* ------ Viscosity and diffusion ---------------- */

  /*! Artificial viscosity parameters */
  struct viscosity_global_data viscosity;

  /*! Thermal diffusion parameters */
  struct diffusion_global_data diffusion;
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
