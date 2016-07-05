/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_CELL_H
#define SWIFT_CELL_H

/* Includes. */
#include <stddef.h>

/* Local includes. */
#include "lock.h"
#include "multipole.h"
#include "part.h"
#include "task.h"

/* Avoid cyclic inclusions */
struct space;

/* Max tag size set to 2^29 to take into account some MPI implementations
 * that use 2^31 as the upper bound on MPI tags and the fact that
 * cell_next_tag is multiplied by 2 when passed to an MPI function.
 * The maximum was lowered by a further factor of 2 to be on the safe side.*/
#define cell_max_tag (1 << 29)

/* Global variables. */
extern int cell_next_tag;

/* Packed cell. */
struct pcell {

  /* Stats on this cell's particles. */
  double h_max;
  int ti_end_min, ti_end_max;

  /* Number of particles in this cell. */
  int count, gcount;

  /* tag used for MPI communication. */
  int tag;

  /* Relative indices of the cell's progeny. */
  int progeny[8];
};

/* Structure to store the data of a single cell. */
struct cell {

  /* The cell location on the grid. */
  double loc[3];

  /* The cell dimensions. */
  double h[3];

  /* Max radii in this cell. */
  double h_max;

  /* Minimum and maximum end of time step in this cell. */
  int ti_end_min, ti_end_max;

  /* Minimum dimension, i.e. smallest edge of this cell. */
  float dmin;

  /* Maximum slack allowed for particle movement. */
  float slack;

  /* Maximum particle movement in this cell. */
  float dx_max;

  /* The depth of this cell in the tree. */
  int depth, split, maxdepth;

  /* Nr of parts. */
  int count, gcount;

  /* Pointers to the particle data. */
  struct part *parts;

  /* Pointers to the extra particle data. */
  struct xpart *xparts;

  /* Pointers to the gravity particle data. */
  struct gpart *gparts;

  /* Pointers for the sorted indices. */
  struct entry *sort, *gsort;
  unsigned int sorted, gsorted;

  /* Pointers to the next level of cells. */
  struct cell *progeny[8];

  /* Parent cell. */
  struct cell *parent;

  /* Super cell, i.e. the highest-level supercell that has interactions. */
  struct cell *super;

  /* The task computing this cell's sorts. */
  struct task *sorts, *gsorts;
  int sortsize, gsortsize;

  /* The tasks computing this cell's density. */
  struct link *density, *force, *grav;
  int nr_density, nr_force, nr_grav;

  /* The hierarchical tasks. */
  struct task *ghost, *init, *drift, *kick;

  /* Task receiving data. */
  struct task *recv_xv, *recv_rho, *recv_ti;

  /* Task send data. */
  struct link *send_xv, *send_rho, *send_ti;

  /* Tasks for gravity tree. */
  struct task *grav_up, *grav_down;

  /* Task for external gravity */
  struct task *grav_external;

  /* Number of tasks that are associated with this cell. */
  int nr_tasks;

  /* Is the data of this cell being used in a sub-cell? */
  int hold, ghold;

  /* Spin lock for various uses. */
  swift_lock_type lock, glock;

  /* ID of the previous owner, e.g. runner. */
  int owner;

  /* Momentum of particles in cell. */
  double mom[3], ang_mom[3];

  /* Mass, potential, internal  and kinetic energy of particles in this cell. */
  double mass, e_pot, e_int, e_kin;

  /* Number of particles updated in this cell. */
  int updated, g_updated;

  /* Linking pointer for "memory management". */
  struct cell *next;

  /* ID of the node this cell lives on. */
  int nodeID;

  /* Bit mask of the proxies this cell is registered with. */
  unsigned long long int sendto;

  /* Pointer to this cell's packed representation. */
  struct pcell *pcell;
  int pcell_size;
  int tag;

  /* This cell's multipole. */
  struct multipole multipole;

} __attribute__((aligned(64)));

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/* Function prototypes. */
void cell_split(struct cell *c, ptrdiff_t parts_offset);
int cell_locktree(struct cell *c);
void cell_unlocktree(struct cell *c);
int cell_glocktree(struct cell *c);
void cell_gunlocktree(struct cell *c);
int cell_pack(struct cell *c, struct pcell *pc);
int cell_unpack(struct pcell *pc, struct cell *c, struct space *s);
int cell_pack_ti_ends(struct cell *c, int *ti_ends);
int cell_unpack_ti_ends(struct cell *c, int *ti_ends);
int cell_getsize(struct cell *c);
int cell_link_parts(struct cell *c, struct part *parts);
int cell_link_gparts(struct cell *c, struct gpart *gparts);
void cell_init_parts(struct cell *c, void *data);
void cell_init_gparts(struct cell *c, void *data);
void cell_convert_hydro(struct cell *c, void *data);
void cell_clean_links(struct cell *c, void *data);

#endif /* SWIFT_CELL_H */
