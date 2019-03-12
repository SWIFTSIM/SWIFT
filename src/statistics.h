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
#ifndef SWIFT_STATISTICS_H
#define SWIFT_STATISTICS_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "lock.h"
#include "space.h"

/**
 * @brief Quantities collected for physics statistics
 */
struct statistics {

  /*! Kinetic energy */
  double E_kin;

  /*! Internal energy */
  double E_int;

  /*! Self potential energy */
  double E_pot_self;

  /*! External potential energy */
  double E_pot_ext;

  /*! Radiative energy */
  double E_rad;

  /*! Entropy */
  double entropy;

  /*! Mass */
  double mass;

  /*! Momentum */
  double mom[3];

  /*! Angular momentum */
  double ang_mom[3];

  /*! Centre of mass */
  double centre_of_mass[3];

  /*! Lock for threaded access */
  swift_lock_type lock;
};

void stats_collect(const struct space* s, struct statistics* stats);
void stats_add(struct statistics* a, const struct statistics* b);
void stats_print_to_file(FILE* file, const struct statistics* stats,
                         double time);
void stats_init(struct statistics* s);
void stats_finalize(struct statistics* s);

#ifdef WITH_MPI
extern MPI_Datatype statistics_mpi_type;
extern MPI_Op statistics_mpi_reduce_op;

void stats_create_mpi_type(void);
#endif

#endif /* SWIFT_STATISTICS_H */
