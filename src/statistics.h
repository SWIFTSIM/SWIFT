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
#ifndef SWIFT_STATISTICS_H
#define SWIFT_STATISTICS_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "lock.h"

/* Some standard headers. */
#include <stdio.h>

/* Pre-declarations */
struct phys_const;
struct space;
struct unit_system;

/**
 * @brief Quantities collected for physics statistics
 */
struct statistics {

  /*! Kinetic energy (internal units)*/
  double E_kin;

  /*! Internal energy (internal units)*/
  double E_int;

  /*! Self potential energy (internal units)*/
  double E_pot_self;

  /*! External potential energy (internal units)*/
  double E_pot_ext;

  /*! Radiated energy (internal units) */
  double E_rad;

  /*! Entropy (internal units) */
  double entropy;

  /*! Total mass (internal units)*/
  double total_mass;

  /*! Total dm mass (internal units)*/
  double dm_mass;

  /*! Total gas mass (internal units)*/
  double gas_mass;

  /*! Total sink mass (internal units)*/
  double sink_mass;

  /*! Total stellar mass (internal units)*/
  double star_mass;

  /*! Total BH mass (internal units)*/
  double bh_mass;

  /*! Total BH subgrid mass (internal units)*/
  double bh_subgrid_mass;

  /*! Total metal mass in gas (internal units)*/
  double gas_Z_mass;

  /*! Total metal mass in stars (internal units)*/
  double star_Z_mass;

  /*! Total metal mass in BH (internal units)*/
  double bh_Z_mass;

  /*! Sum of instantaneous accretion rate of all BHs (internal units)*/
  double bh_accretion_rate;

  /*! Total accreted mass of all BHs (internal units)*/
  double bh_accreted_mass;

  /* Total BH bolometric luminosity of all BHs (internal units) */
  double bh_bolometric_luminosity;

  /* Total jet power of all BHs (internal units) */
  double bh_jet_power;

  /*! Momentum (internal units)*/
  double mom[3];

  /*! Angular momentum (internal units) */
  double ang_mom[3];

  /*! Centre of mass (internal units)*/
  double centre_of_mass[3];

  /*! Total gas mass that is in Hydrogen (all species) */
  double gas_H_mass;

  /*! Total gas mass that is in Molecular Hydrogen */
  double gas_H2_mass;

  /*! Total gas mass that is in Atomic Hydrogen */
  double gas_HI_mass;

  /*! Total gas mass that is in Helium (all species) */
  double gas_He_mass;

  /*! Total Magnetic Energy */
  double E_mag;

  /*! Total divB error */
  double divB_error;

  /*! Total Cross Helicity */
  double H_cross;

  /*! Total Magnetic helicity */
  double H_mag;

  /*! Lock for threaded access */
  swift_lock_type lock;
};

void stats_collect(const struct space* s, struct statistics* stats);
void stats_add(struct statistics* a, const struct statistics* b);
void stats_write_file_header(FILE* file, const struct unit_system* us,
                             const struct phys_const* phys_const);
void stats_write_to_file(FILE* file, const struct statistics* stats,
                         const double time, const double a, const double z,
                         const int step);
void stats_init(struct statistics* s);
void stats_finalize(struct statistics* s);

#ifdef WITH_MPI
extern MPI_Datatype statistics_mpi_type;
extern MPI_Op statistics_mpi_reduce_op;

void stats_create_mpi_type(void);
void stats_free_mpi_type(void);
#endif

#endif /* SWIFT_STATISTICS_H */
